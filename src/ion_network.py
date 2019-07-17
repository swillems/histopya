#!venv/bin/python


import src.io
import src.ions
import src.aggregates
import src.peptides
import src.parameters
import src.gui
import numpy as np


class IonNetwork(object):

    def __init__(self, input_file_name):
        self.parameters = src.parameters.importParameters(input_file_name)

    def hasIons(self):
        return hasattr(self, "ions")

    def hasQuickClusters(self):
        # TODO define
        return hasattr(self, "quick_clusters")

    def isRTCalibrated(self):
        if not self.hasIons():
            return False
        if self.ions[0]["RT"] != self.ions[0]["CALIBRATED_RT"]:
            return True
        return np.any(self.ions["RT"] != self.ions["CALIBRATED_RT"])

    def isDTCalibrated(self):
        if not self.hasIons():
            return False
        if np.any(self.ions[0]["DT"] != self.ions[0]["CALIBRATED_DT"]):
            return True
        return np.any(self.ions["DT"] != self.ions["CALIBRATED_DT"])

    def isMZCalibrated(self):
        if not self.hasIons():
            return False
        if np.any(self.ions[0]["MZ"] != self.ions[0]["CALIBRATED_MZ"]):
            return True
        return np.any(self.ions["MZ"] != self.ions["CALIBRATED_MZ"])

    def isPairwiseAligned(self):
        if not self.hasIons():
            return False
        if np.any(self.ions[0]["AGGREGATE_INDEX"] > 0):
            return True
        return np.any(self.ions["AGGREGATE_INDEX"] > 0)

    def isNormalized(self):
        if not self.hasIons():
            return False
        if np.any(self.ions[0]["INTENSITY"] != self.ions[0]["CALIBRATED_INTENSITY"]):
            return True
        return np.any(self.ions["INTENSITY"] != self.ions["CALIBRATED_INTENSITY"])

    def hasNodes(self):
        return (
            hasattr(
                self,
                "aggregates"
            ) and hasattr(
                self,
                "aggregate_ions"
            )
        )

    def hasEdges(self):
        return hasattr(self, "neighbors")

    def isCreated(self):
        return (self.hasNodes() and self.hasEdges())

    def hasFragmentTargets(self):
        return hasattr(self, "fragments")

    def isAnnotated(self):
        return (
            hasattr(
                self,
                "anchor_peptide_scores"
            ) and hasattr(
                self,
                "anchor_peptide_match_counts"
            )
        )

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
                # src.io.plotEstimates(
                #     estimation_aggregates,
                #     self.parameters,
                #     log
                # )
                # src.io.plotCalibrationResults(estimation_aggregates, self.parameters, log)
                ions, neighbors = src.ions.detectAllIonNeighbors(
                    ions,
                    ion_alignment_parameters,
                    self.parameters,
                    log,
                    save=self.parameters["SAVE_UNFILTERED_IONS"]
                )
                # importlib.reload(src.ions)
                aggregates, aggregate_ions, ions = src.aggregates.defineFromIons(
                    ions,
                    self.parameters,
                    log
                )
                # TODO
                # src.io.plotAnchorCounts(self.parameters, aggregates, ions, log)
                ions = src.ions.calibrateIntensities(
                    ions,
                    aggregate_ions,
                    self.parameters,
                    log
                )
                quick_isotopic_pairs = src.aggregates.findQuickIsotopes(
                    ions,
                    aggregates,
                    aggregate_ions,
                    ion_alignment_parameters,
                    self.parameters,
                    log
                )
                anchor_alignment_parameters = src.aggregates.estimateAlignmentParameters(
                    ions,
                    quick_isotopic_pairs,
                    aggregate_ions,
                    self.parameters,
                    log
                )
                aggregates, ions = src.aggregates.calibrateChannelDriftShift(
                    ions,
                    aggregate_ions,
                    aggregates,
                    ion_alignment_parameters,
                    anchor_alignment_parameters,
                    self.parameters,
                    log
                )
                # anchor_alignment_parameters['DT'] = (np.array(anchor_alignment_parameters['DT'])*3).tolist()
                # anchor_alignment_parameters['RT'] = (np.array(anchor_alignment_parameters['RT'])*3).tolist()
                neighbors = src.aggregates.detectAllAnchorNeighbors(
                    ions,
                    aggregates,
                    aggregate_ions,
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
                aggregates = src.io.loadArray("ANCHORS_FILE_NAME", parameters)
                anchor_alignment_parameters = src.io.loadJSON(
                    "ANCHOR_ALIGNMENT_PARAMETERS_FILE_NAME",
                    parameters,
                )
                neighbors = src.io.loadMatrix(
                    "ANCHOR_NEIGHBORS_FILE_NAME",
                    parameters,
                )
                neighbors += neighbors.T
                # src.aggregates.writeMgf(neighbors, aggregates, aggregate_ions, ions, parameters, log)
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
                    aggregates,
                    base_mass_dict,
                    parameters,
                    log
                )
                anchor_peptide_scores, anchor_peptide_match_counts = src.aggregates.getAnchorPeptideMatrix(
                    aggregates,
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
                    aggregates,
                    neighbors,
                    anchor_alignment_parameters,
                    peptide_masses,
                    base_mass_dict,
                    parameters,
                    log
                )
                annotation_data = src.aggregates.writePercolatorFile(
                    aggregates,
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
        # TODO
        print("Dummy for exportAggregates function")

    def exportIons(self):
        # TODO
        print("Dummy for exportIons function")

    def exportEdges(self):
        # TODO
        print("Dummy for exportEdges function")

    def exportAnnotations(self):
        # TODO
        print("Dummy for exportAnnotations function")
