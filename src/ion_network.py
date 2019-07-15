#!venv/bin/python


import src.io
import src.ions
import src.aggregates
import src.peptides
import src.parameters
import src.gui


class IonNetwork(object):

    def __init__(self, input_file_name, pre_exists):
        self.is_created = False
        self.is_annotated = False
        if not pre_exists:
            self.parameters = src.parameters.updateParameters(
                src.parameters.getDefaultParameters(),
                input_file_name
            )
        else:
            self.parameters = src.io.loadParametersFromINET(input_file_name)

    def create(self):
        with src.io.Log(self.parameters["LOG_FILE_NAME"]) as log:
            with log.newSection("Creating ion-network"):
                ions = src.ions.importAllFromCsv(self.parameters, log)
                pseudo_aggregate_ions_indices = src.ions.getPseudoAggregatesIndices(
                    ions,
                    self.parameters,
                    log
                )
                calibration_aggregates = ions[pseudo_aggregate_ions_indices][::2]
                ions = src.ions.calibrateAll(
                    ions,
                    calibration_aggregates,
                    self.parameters,
                    log
                )
                src.io.saveArray(
                    ions[pseudo_aggregate_ions_indices],
                    "PSEUDO_AGGREGATE_IONS_FILE_NAME",
                    self.parameters,
                    log
                )
                estimation_aggregates = ions[pseudo_aggregate_ions_indices][1::2]
                ion_alignment_parameters = src.ions.estimateAlignmentParameters(
                    estimation_aggregates,
                    self.parameters,
                    log,
                )
                # TODO plotting
                src.io.plotEstimates(
                    estimation_aggregates,
                    self.parameters,
                    log
                )
                ions, neighbors = src.ions.detectAllIonNeighbors(
                    ions,
                    ion_alignment_parameters,
                    self.parameters,
                    log,
                    save=self.parameters["SAVE_UNFILTERED_IONS"]
                )
                # importlib.reload(src.ions)
                anchors, anchor_ions, ions = src.aggregates.defineFromIons(
                    ions,
                    self.parameters,
                    log
                )
                # TODO
                src.io.plotAnchorCounts(self.parameters, anchors, ions, log)
                ions = src.ions.calibrateIntensities(
                    ions,
                    anchor_ions,
                    self.parameters,
                    log
                )
                quick_isotopic_pairs = src.aggregates.findQuickIsotopes(
                    ions,
                    anchors,
                    anchor_ions,
                    ion_alignment_parameters,
                    self.parameters,
                    log
                )
                anchor_alignment_parameters = src.aggregates.estimateAlignmentParameters(
                    ions,
                    quick_isotopic_pairs,
                    anchor_ions,
                    self.parameters,
                    log
                )
                anchors, ions = src.aggregates.calibrateChannelDriftShift(
                    ions,
                    anchor_ions,
                    anchors,
                    ion_alignment_parameters,
                    anchor_alignment_parameters,
                    self.parameters,
                    log
                )
                # anchor_alignment_parameters['DT'] = (np.array(anchor_alignment_parameters['DT'])*3).tolist()
                # anchor_alignment_parameters['RT'] = (np.array(anchor_alignment_parameters['RT'])*3).tolist()
                neighbors = src.aggregates.detectAllAnchorNeighbors(
                    ions,
                    anchors,
                    anchor_ions,
                    anchor_alignment_parameters,
                    self.parameters,
                    log
                )
            log.printMessage("Created ion-network")
        self.is_created = True

    def annotate(self):
        parameters = self.parameters
        with src.io.Log(parameters["LOG_FILE_NAME"]) as log:
            with log.newSection("Annotating ion-network"):
                anchors = src.io.loadArray("ANCHORS_FILE_NAME", parameters)
                anchor_ions = src.io.loadMatrix(
                    "ANCHOR_IONS_FILE_NAME",
                    parameters,
                )
                ions = src.io.loadArray("IONS_FILE_NAME", parameters)
                ion_alignment_parameters = src.io.loadJSON(
                    "ION_ALIGNMENT_PARAMETERS_FILE_NAME",
                    parameters,
                )
                anchor_alignment_parameters = src.io.loadJSON(
                    "ANCHOR_ALIGNMENT_PARAMETERS_FILE_NAME",
                    parameters,
                )
                neighbors = src.io.loadMatrix(
                    "ANCHOR_NEIGHBORS_FILE_NAME",
                    parameters,
                )
                neighbors += neighbors.T
                # src.aggregates.writeMgf(neighbors, anchors, anchor_ions, ions, parameters, log)
                base_mass_dict = src.peptides.loadBaseMassDict(parameters, log)
                proteins, total_protein_sequence, ptms, ptm_matrix = src.peptides.importProteinsAndPtms(
                    parameters,
                    log
                )
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
        self.is_annotated = True

    def browse(self):
        parameters = self.parameters
        dataset = src.gui.Dataset(parameters)
        gui = src.gui.GUI(dataset)

    def exportAggregates(self):
        print("Dummy for exportAggregates function")

    def exportIons(self):
        print("Dummy for exportIons function")

    def exportEdges(self):
        print("Dummy for exportEdges function")

    def exportAnnotations(self):
        print("Dummy for exportAnnotations function")
