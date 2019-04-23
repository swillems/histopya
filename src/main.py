#!venv/bin/python


import src.io
import src.ions
import src.aggregates
import src.peptides


def main(parameters):
    with src.io.Log(parameters["LOG_FILE_NAME"]) as log:
        with log.newSection("Starting analysis"):
            ions = src.ions.importAllFromCsv(parameters, log)
            pseudo_aggregate_ions_indices = src.ions.getPseudoAggregatesIndices(ions, parameters, log)
            calibration_aggregates = ions[pseudo_aggregate_ions_indices][::2]
            ions = src.ions.calibrateAll(
                ions,
                calibration_aggregates,
                parameters,
                log
            )
            estimation_aggregates = ions[pseudo_aggregate_ions_indices][1::2]
            ion_alignment_parameters = src.ions.estimateAlignmentParameters(
                estimation_aggregates,
                parameters,
                log,
            )
            # TODO plotting
            src.io.plotEstimates(
                estimation_aggregates,
                parameters,
                log
            )
            ions, neighbors = src.ions.detectAllIonNeighbors(
                ions,
                ion_alignment_parameters,
                parameters,
                log
            )
            # importlib.reload(src.ions)
            anchors, anchor_ions, ions = src.aggregates.defineFromIons(
                ions,
                parameters,
                log
            )
            # TODO
            src.io.plotAnchorCounts(parameters, anchors, ions, log)
            ions = src.ions.calibrateIntensities(
                ions,
                anchor_ions,
                parameters,
                log
            )
            quick_isotopic_pairs = src.aggregates.findQuickIsotopes(
                ions,
                anchors,
                anchor_ions,
                ion_alignment_parameters,
                parameters,
                log
            )
            anchor_alignment_parameters = src.aggregates.estimateAlignmentParameters(
                ions,
                quick_isotopic_pairs,
                anchor_ions,
                parameters,
                log
            )
            anchors, ions = src.aggregates.calibrateChannelDriftShift(
                ions,
                anchor_ions,
                anchors,
                ion_alignment_parameters,
                anchor_alignment_parameters,
                parameters,
                log
            )
            # anchor_alignment_parameters['DT'] = (np.array(anchor_alignment_parameters['DT'])*3).tolist()
            # anchor_alignment_parameters['RT'] = (np.array(anchor_alignment_parameters['RT'])*3).tolist()
            neighbors = src.aggregates.detectAllAnchorNeighbors(
                ions,
                anchors,
                anchor_ions,
                anchor_alignment_parameters,
                parameters,
                log
            )
            neighbors += neighbors.T
            # src.aggregates.writeMgf(neighbors, anchors, anchor_ions, ions, parameters, log)
            base_mass_dict = src.peptides.loadBaseMassDict(parameters, log)
            proteins, total_protein_sequence, ptms, ptm_matrix = src.peptides.importProteinsAndPtms(parameters, log)
            peptides, peptide_index_matrix, digestion_matrix = src.peptides.digestProteins(
                proteins,
                total_protein_sequence,
                ptm_matrix,
                parameters,
                log,
            )
            peptide_masses, fragments = src.peptides.calculateMasses(
                peptides,
                peptide_index_matrix,
                base_mass_dict,
                total_protein_sequence,
                parameters,
                log,
            )
            anchor_boundaries, fragment_peptide_indices, fragment_indices = src.aggregates.matchAnchorsToFragments(
                fragments,
                anchors,
                base_mass_dict,
                parameters,
                log
            )
            anchor_peptide_scores, anchor_peptide_match_counts = src.aggregates.getAnchorPeptideMatrix(
                anchors,
                neighbors,
                peptides,
                anchor_boundaries,
                fragment_peptide_indices,
                parameters,
                log
            )
            anchor_fragment_indices = src.aggregates.getAnchorFragmentIndices(
                anchor_peptide_match_counts,
                anchor_boundaries,
                fragment_indices,
                fragment_peptide_indices,
                parameters,
                log
            )
            precursor_indices = src.aggregates.findFragmentPrecursors(
                anchor_peptide_match_counts,
                anchors,
                neighbors,
                anchor_alignment_parameters,
                peptide_masses,
                base_mass_dict,
                parameters,
                log
            )
            annotation_data = src.aggregates.writePercolatorFile(
                anchors,
                base_mass_dict,
                anchor_peptide_match_counts,
                fragments,
                anchor_fragment_indices,
                neighbors,
                peptides,
                proteins,
                peptide_masses,
                precursor_indices,
                anchor_peptide_scores,
                peptide_index_matrix,
                total_protein_sequence,
                parameters,
                log
            )
            src.io.runPercolator(parameters, log)
        log.printMessage("Analysis completed")


if __name__ == "__main__":
    import src.parameters
    parameters = src.parameters.parseFromCommandLine()
    main(parameters)


def runAllWorkFlowSteps(settings):
    with src.io.Log(settings["LOG_FILE_NAME"]) as log:
        with log.newSection("Running HistoPyA workflow"):
            variables = {}
            parameters = settings["PARAMETERS"]
            parameters.update(settings["PARAMETERS"]["IO"]["PATHS"])
            parameters.update(settings["PARAMETERS"]["IO"]["STATIC"])
            parameters.update(settings["PARAMETERS"]["IO"]["DYNAMIC"])
            for step in settings["WORKFLOW"]:
                variables = runSingleWorkflowStep(
                    step,
                    variables,
                    parameters,
                    log
                )
            log.printMessage("Finished all workflow steps")


def runSingleWorkflowStep(
    step,
    variables,
    parameters,
    log
):
    if step == "IMPORT_APEX_FILES":
        variables["ions"] = src.ions.importAllFromCsv(parameters, log)
    elif step == "DETECT_PSEUDO_AGGREGATES":
        variables["pseudo_aggregate_ions_indices"] = src.ions.getPseudoAggregatesIndices(
            variables["ions"],
            parameters,
            log
        )
    elif step == "PARTITION_PSEUDO_AGGREGATES":
        # calibration_aggregates = ions[pseudo_aggregate_ions_indices][::2]
        variables["calibration_aggregate_indices"] = variables["pseudo_aggregate_ions_indices"][::2]
        variables["estimation_aggregate_indices"] = variables["pseudo_aggregate_ions_indices"][1::2]
    elif step == "CALIBRATE_SAMPLE_RTS":
        ions = src.ions.calibrateAll(
            ions,
            calibration_aggregates,
            parameters,
            log
        )
        estimation_aggregates = ions[pseudo_aggregate_ions_indices][1::2]
        ion_alignment_parameters = src.ions.estimateAlignmentParameters(
            estimation_aggregates,
            parameters,
            log,
        )
        # TODO plotting
        src.io.plotEstimates(
            estimation_aggregates,
            parameters,
            log
        )
        ions, neighbors = src.ions.detectAllIonNeighbors(
            ions,
            ion_alignment_parameters,
            parameters,
            log
        )
        # importlib.reload(src.ions)
        anchors, anchor_ions, ions = src.aggregates.defineFromIons(
            ions,
            parameters,
            log
        )
        # TODO
        src.io.plotAnchorCounts(parameters, anchors, ions, log)
        ions = src.ions.calibrateIntensities(
            ions,
            anchor_ions,
            parameters,
            log
        )
        quick_isotopic_pairs = src.aggregates.findQuickIsotopes(
            ions,
            anchors,
            anchor_ions,
            ion_alignment_parameters,
            parameters,
            log
        )
        anchor_alignment_parameters = src.aggregates.estimateAlignmentParameters(
            ions,
            quick_isotopic_pairs,
            anchor_ions,
            parameters,
            log
        )
        anchors, ions = src.aggregates.calibrateChannelDriftShift(
            ions,
            anchor_ions,
            anchors,
            parameters,
            log
        )
        # anchor_alignment_parameters['DT'] = (np.array(anchor_alignment_parameters['DT'])*3).tolist()
        # anchor_alignment_parameters['RT'] = (np.array(anchor_alignment_parameters['RT'])*3).tolist()
        neighbors = src.aggregates.detectAllAnchorNeighbors(
            ions,
            anchors,
            anchor_ions,
            anchor_alignment_parameters,
            parameters,
            log
        )
        neighbors += neighbors.T
        # src.aggregates.writeMgf(neighbors, anchors, anchor_ions, ions, parameters, log)
        base_mass_dict = src.peptides.loadBaseMassDict(parameters, log)
        proteins, total_protein_sequence, ptms, ptm_matrix = src.peptides.importProteinsAndPtms(parameters, log)
        peptides, peptide_index_matrix, digestion_matrix = src.peptides.digestProteins(
            proteins,
            total_protein_sequence,
            ptm_matrix,
            parameters,
            log,
        )
        peptide_masses, fragments = src.peptides.calculateMasses(
            peptides,
            peptide_index_matrix,
            base_mass_dict,
            total_protein_sequence,
            parameters,
            log,
        )
        anchor_boundaries, fragment_peptide_indices, fragment_indices = src.aggregates.matchAnchorsToFragments(
            fragments,
            anchors,
            base_mass_dict,
            parameters,
            log
        )
        anchor_peptide_scores, anchor_peptide_match_counts = src.aggregates.getAnchorPeptideMatrix(
            anchors,
            neighbors,
            peptides,
            anchor_boundaries,
            fragment_peptide_indices,
            parameters,
            log
        )
        anchor_fragment_indices = src.aggregates.getAnchorFragmentIndices(
            anchor_peptide_match_counts,
            anchor_boundaries,
            fragment_indices,
            fragment_peptide_indices,
            parameters,
            log
        )
        precursor_indices = src.aggregates.findFragmentPrecursors(
            anchor_peptide_match_counts,
            anchors,
            neighbors,
            anchor_alignment_parameters,
            peptide_masses,
            base_mass_dict,
            parameters,
            log
        )
        annotation_data = src.aggregates.writePercolatorFile(
            anchors,
            base_mass_dict,
            anchor_peptide_match_counts,
            fragments,
            anchor_fragment_indices,
            neighbors,
            peptides,
            proteins,
            peptide_masses,
            precursor_indices,
            anchor_peptide_scores,
            peptide_index_matrix,
            total_protein_sequence,
            parameters,
            log
        )
        src.io.runPercolator(parameters, log)
    return variables
