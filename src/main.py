#!venv/bin/python


import src.io
import src.ions
import src.aggregates
import src.peptides
import src.gui


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
            src.io.saveArray(
                ions[pseudo_aggregate_ions_indices],
                "PSEUDO_AGGREGATE_IONS_FILE_NAME",
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
                log,
                save=parameters["SAVE_UNFILTERED_IONS"]
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



def create(parameters):
    with src.io.Log(parameters["LOG_FILE_NAME"]) as log:
        with log.newSection("Creating ion-network"):
            ions = src.ions.importAllFromCsv(parameters, log)
            pseudo_aggregate_ions_indices = src.ions.getPseudoAggregatesIndices(ions, parameters, log)
            calibration_aggregates = ions[pseudo_aggregate_ions_indices][::2]
            ions = src.ions.calibrateAll(
                ions,
                calibration_aggregates,
                parameters,
                log
            )
            src.io.saveArray(
                ions[pseudo_aggregate_ions_indices],
                "PSEUDO_AGGREGATE_IONS_FILE_NAME",
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
                log,
                save=parameters["SAVE_UNFILTERED_IONS"]
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
        log.printMessage("Analysis completed")


def annotate(parameters):
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


if __name__ == "__main__":
    import src.parameters
    parameters, action = src.parameters.parseFromCommandLine()
    if action == "F":
        main(parameters)
    if action == "C":
        create(parameters)
    if action == "A":
        annotate(parameters)
    if action == "B":
        dataset = src.gui.Dataset(parameters)
        gui = src.gui.GUI(dataset)


#
# def runAllWorkFlowSteps(settings):
#     with src.io.Log(settings["LOG_FILE_NAME"]) as log:
#         with log.newSection("Running HistoPyA workflow"):
#             variables = {}
#             parameters = settings["PARAMETERS"]
#             parameters.update(settings["PARAMETERS"]["IO"]["PATHS"])
#             parameters.update(settings["PARAMETERS"]["IO"]["STATIC"])
#             parameters.update(settings["PARAMETERS"]["IO"]["DYNAMIC"])
#             for step in settings["WORKFLOW"]:
#                 variables = runSingleWorkflowStep(
#                     step,
#                     variables,
#                     parameters,
#                     log
#                 )
#             log.printMessage("Finished all workflow steps")
#
#
# def runSingleWorkflowStep(
#     step,
#     variables,
#     parameters,
#     log
# ):
#     if step == "IMPORT_APEX_FILES":
#         variables["ions"] = src.ions.importAllFromCsv(parameters, log)
#     elif step == "DETECT_PSEUDO_AGGREGATES":
#         variables["pseudo_aggregate_ions_indices"] = src.ions.getPseudoAggregatesIndices(
#             variables["ions"],
#             parameters,
#             log
#         )
#     elif step == "PARTITION_PSEUDO_AGGREGATES":
#         # calibration_aggregates = ions[pseudo_aggregate_ions_indices][::2]
#         variables["calibration_aggregate_indices"] = variables["pseudo_aggregate_ions_indices"][::2]
#         variables["estimation_aggregate_indices"] = variables["pseudo_aggregate_ions_indices"][1::2]
#     elif step == "CALIBRATE_SAMPLE_RTS":
#         ions = src.ions.calibrateAll(
#             ions,
#             calibration_aggregates,
#             parameters,
#             log
#         )
#         estimation_aggregates = ions[pseudo_aggregate_ions_indices][1::2]
#         ion_alignment_parameters = src.ions.estimateAlignmentParameters(
#             estimation_aggregates,
#             parameters,
#             log,
#         )
#         # TODO plotting
#         src.io.plotEstimates(
#             estimation_aggregates,
#             parameters,
#             log
#         )
#         ions, neighbors = src.ions.detectAllIonNeighbors(
#             ions,
#             ion_alignment_parameters,
#             parameters,
#             log
#         )
#         # importlib.reload(src.ions)
#         anchors, anchor_ions, ions = src.aggregates.defineFromIons(
#             ions,
#             parameters,
#             log
#         )
#         # TODO
#         src.io.plotAnchorCounts(parameters, anchors, ions, log)
#         ions = src.ions.calibrateIntensities(
#             ions,
#             anchor_ions,
#             parameters,
#             log
#         )
#         quick_isotopic_pairs = src.aggregates.findQuickIsotopes(
#             ions,
#             anchors,
#             anchor_ions,
#             ion_alignment_parameters,
#             parameters,
#             log
#         )
#         anchor_alignment_parameters = src.aggregates.estimateAlignmentParameters(
#             ions,
#             quick_isotopic_pairs,
#             anchor_ions,
#             parameters,
#             log
#         )
#         anchors, ions = src.aggregates.calibrateChannelDriftShift(
#             ions,
#             anchor_ions,
#             anchors,
#             parameters,
#             log
#         )
#         # anchor_alignment_parameters['DT'] = (np.array(anchor_alignment_parameters['DT'])*3).tolist()
#         # anchor_alignment_parameters['RT'] = (np.array(anchor_alignment_parameters['RT'])*3).tolist()
#         neighbors = src.aggregates.detectAllAnchorNeighbors(
#             ions,
#             anchors,
#             anchor_ions,
#             anchor_alignment_parameters,
#             parameters,
#             log
#         )
#         neighbors += neighbors.T
#         # src.aggregates.writeMgf(neighbors, anchors, anchor_ions, ions, parameters, log)
#         base_mass_dict = src.peptides.loadBaseMassDict(parameters, log)
#         proteins, total_protein_sequence, ptms, ptm_matrix = src.peptides.importProteinsAndPtms(parameters, log)
#         peptides, peptide_index_matrix, digestion_matrix = src.peptides.digestProteins(
#             proteins,
#             total_protein_sequence,
#             ptm_matrix,
#             parameters,
#             log,
#         )
#         peptide_masses, fragments = src.peptides.calculateMasses(
#             peptides,
#             peptide_index_matrix,
#             base_mass_dict,
#             total_protein_sequence,
#             parameters,
#             log,
#         )
#         anchor_boundaries, fragment_peptide_indices, fragment_indices = src.aggregates.matchAnchorsToFragments(
#             fragments,
#             anchors,
#             base_mass_dict,
#             parameters,
#             log
#         )
#         anchor_peptide_scores, anchor_peptide_match_counts = src.aggregates.getAnchorPeptideMatrix(
#             anchors,
#             neighbors,
#             peptides,
#             anchor_boundaries,
#             fragment_peptide_indices,
#             parameters,
#             log
#         )
#         anchor_fragment_indices = src.aggregates.getAnchorFragmentIndices(
#             anchor_peptide_match_counts,
#             anchor_boundaries,
#             fragment_indices,
#             fragment_peptide_indices,
#             parameters,
#             log
#         )
#         precursor_indices = src.aggregates.findFragmentPrecursors(
#             anchor_peptide_match_counts,
#             anchors,
#             neighbors,
#             anchor_alignment_parameters,
#             peptide_masses,
#             base_mass_dict,
#             parameters,
#             log
#         )
#         annotation_data = src.aggregates.writePercolatorFile(
#             anchors,
#             base_mass_dict,
#             anchor_peptide_match_counts,
#             fragments,
#             anchor_fragment_indices,
#             neighbors,
#             peptides,
#             proteins,
#             peptide_masses,
#             precursor_indices,
#             anchor_peptide_scores,
#             peptide_index_matrix,
#             total_protein_sequence,
#             parameters,
#             log
#         )
#         src.io.runPercolator(parameters, log)
#     return variables


# taskset -c 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86 wine64 /media/proteomics/MISC/Sander/UW/UniversalWorkflow_Software/Apx2D/Apex2D_v3.1.0.9.5/Apex3D64.exe -pRawDirName /media/proteomics/RAW/External/20130423_Tenzer_UDMSE_LFQ/S130423_05.raw -outputDirName /media/proteomics/MISC/Sander/APEX/20130423_Tenzer_UDMSE_LFQ_apex2d -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 1 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -leCentroidThresholdCounts 1 -heCentroidThresholdCounts 1 -startingRTMin 95 -endingRTMin 99



# for i in {10..14}; do taskset -c 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86 wine64 /media/proteomics/MISC/Sander/UW/UniversalWorkflow_Software/Apx2D/Apex2D_v3.1.0.9.5/Apex3D64.exe -pRawDirName /media/proteomics/RAW/External/20130423_Tenzer_UDMSE_LFQ/S130423_"$i".raw -outputDirName /media/proteomics/MISC/Sander/APEX/20130423_Tenzer_UDMSE_LFQ_apex2d/"$i" -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 1 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -leCentroidThresholdCounts 1 -heCentroidThresholdCounts 1; done

# for i in {1..5}; do taskset -c 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86 wine64 /media/proteomics/MISC/Sander/UW/UniversalWorkflow_Software/Apx2D/Apex2D_v3.1.0.9.5/Apex3D64.exe -pRawDirName "/media/proteomics/Thesis 2019/RAW data/190513_SynaptData/190508_ProteomeComparison_K562_0"$i"_A.raw" -outputDirName /media/proteomics/MISC/Sander/APEX/20190508_Thesis -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -bEnableCentroids 0; done

# for i in {1..3}; do apex3d -pRawDirName "/media/proteomics/RAW/External/20181203_Lagies_HDMSE_Metabolomics/20181203_Panc1_"$i".raw" -outputDirName /media/proteomics/MISC/Sander/APEX/20181203_Lagies_HDMSE_Metabolomics -lockMassZ1 556.2772 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -bEnableCentroids 0 -leThresholdCounts 1 -heThresholdCounts 1 -function 2; done

# for i in {5..9}; do taskset -c 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86 wine64 /media/proteomics/MISC/Sander/UW/UniversalWorkflow_Software/Apx2D/Apex2D_v3.1.0.9.5/Apex3D64.exe -pRawDirName /media/proteomics/RAW/External/20130423_Tenzer_UDMSE_LFQ/S130423_0"$i".raw -outputDirName /media/proteomics/MISC/Sander/APEX/20130423_Tenzer_UDMSE_LFQ_default -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -bEnableCentroids 0 -startingRTMin 90 -endingRTMin 95 ; done

# for i in {10..14}; do apex3d -pRawDirName /media/proteomics/RAW/External/20130423_Tenzer_UDMSE_LFQ/S130423_"$i".raw -outputDirName /media/proteomics/MISC/Sander/APEX/20130423_Tenzer_UDMSE_LFQ_default -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -bEnableCentroids 0; done

# import src.gui; import importlib; importlib.reload(src.gui); d=src.gui.Dataset("data/tenzer/parameters.json"); g=src.gui.GUI(d)

# for i in {10..14}; do awk '(($4 > 90)&& ($4 < 95)) || NR==1' data/tenzer/apex/S130423_"$i"_Apex3DIons.csv > data/test2/apex/sample_a_"$i".csv; done

# percolator percolator.csv --results-psms percolator/target_pims.csv --results-peptides percolator/target_peptides.csv --decoy-results-psms percolator/decoy_pims.csv -I "concatenated" -A --results-proteins percolator/target_proteins.csv -P DECOY_ -x -D 15 2>&1 | tee percolator/percolator_log.txt

# awk '$3 < 0.01' data/480_PharmaFluidics/results/percolator/target_pims.csv | cut -f 6 | cut -d "_" -f 2 | sort | uniq -c



# for i in {0..9}; do awk '(($4 > 50)&& ($4 < 55)) || NR==1' data/ecoli_sonar/apex/28Oct2016_06"$i"_Apex3DIons.csv > data/test_sonar/apex/sample_"$i".csv; done
#
# for i in {10..14}; do awk '(($4 > 90)&& ($4 < 99)) || NR==1' data/tenzer/apex/S130423_"$i"_Apex3DIons.csv > data/test_lfq/apex/sample_A_"$i".csv; done
#
# for i in {1..2}; do awk '(($4 > 90)&& ($4 < 99)) || NR==1' data/lfq_swim_190327/apex/LFQnano_SWIMDIAvsUDMSE_SynaptG2Si_SWIMDIA_Condition_B_0"$i"_Apex3DIons.csv > data/test_lfq/apex/sample_B_"$i".csv; done
#
# for i in {1..1}; do awk '(($4 > 90)&& ($4 < 99)) || NR==1' data/lfq_swim_190327/apex/LFQnano_SWIMDIAvsUDMSE_SynaptG2Si_SWIMDIA_QC_0"$i"_Apex3DIons.csv > data/test_lfq/apex/sample_QC_"$i".csv; done
#
#
# sftp://sander@locutus.ugent.be/home/sander/proteomics/histopya/data/tenzer/apex/S130423_0"$i"_Apex3DIons.csv
