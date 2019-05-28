#!venv/bin/python


import src.parameters
import src.io
import src.ions
import src.aggregates
import src.peptides
import numpy as np
import pandas as pd
import scipy.sparse
from matplotlib import pyplot as plt

# Loading data
parameter_file_name = "data/tenzer/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "_interactive.txt")
with log.newSection("Loading data"):
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
    base_mass_dict = src.peptides.loadBaseMassDict(parameters, log)
    proteins, total_protein_sequence, ptms, ptm_matrix = src.peptides.importProteinsAndPtms(parameters, log)
    peptides, peptide_index_matrix, digestion_matrix = src.peptides.digestProteins(
        proteins,
        total_protein_sequence,
        ptm_matrix,
        parameters,
        log,
    )
    anchor_peptide_scores = src.io.loadMatrix(
        "ANCHOR_PEPTIDE_SCORES_FILE_NAME",
        parameters,
    )
    anchor_peptide_match_counts = src.io.loadMatrix(
        "ANCHOR_PEPTIDE_MATCH_COUNTS_FILE_NAME",
        parameters,
    )

# Unfiltered LFQ calculations
with log.newSection("Calculating unfiltered aggregate LFQ"):
    parameters["UNFILTERED_IONS_FILE_NAME"] = "data/tenzer/results/ions_unfiltered_tmp.npy"
    unfiltered_ions = src.io.loadArray("UNFILTERED_IONS_FILE_NAME", parameters)
    unfiltered_anchors, unfiltered_anchor_ions, unfiltered_ions = src.aggregates.defineFromIons(
        unfiltered_ions,
        parameters,
        log,
        save=False,
        remove_noise=False,
        order_anchors=False
    )
    conditions = {
        "A": slice(0, None, 2),
        "B": slice(1, None, 2),
    }
    anchors_per_condition = {}
    anchor_ions_per_condition = {}
    sample_ions = unfiltered_anchor_ions.T.tocsr()
    for condition, condition_slice in conditions.items():
        condition_anchors = np.empty(
            len(unfiltered_anchors),
            dtype=[
                ("INTENSITY", np.float),
                ("CV", np.float),
                ("ION_COUNT", np.int),
            ]
        )
        condition_anchor_ions = sample_ions[condition_slice].T.tocsr()
        anchor_sizes = np.diff(condition_anchor_ions.indptr)
        intensities = np.zeros(len(unfiltered_anchors), dtype=np.float)
        cvs = np.zeros(len(unfiltered_anchors), dtype=np.float)
        for size in range(1, np.max(anchor_sizes) + 1):
            anchor_subset = np.flatnonzero(anchor_sizes == size)
            anchor_subset_matrix = np.concatenate(
                [
                    condition_anchor_ions.data[
                        condition_anchor_ions.indptr[anchor_index]: condition_anchor_ions.indptr[anchor_index + 1]
                    ] for anchor_index in anchor_subset
                ]
            ).reshape(-1, size)
            intensity_matrix = unfiltered_ions["CALIBRATED_INTENSITY"][anchor_subset_matrix]
            cvs[anchor_subset] = scipy.stats.variation(intensity_matrix, axis=1)
            intensities[anchor_subset] = np.average(intensity_matrix, axis=1)
        condition_anchors["ION_COUNT"] = anchor_sizes
        condition_anchors["CV"] = cvs
        condition_anchors["INTENSITY"] = intensities
        anchors_per_condition[condition] = condition_anchors
        anchor_ions_per_condition[condition] = condition_anchor_ions
    anchors_per_condition["QC"] = {
        "INTENSITY": anchors_per_condition["A"]["INTENSITY"] + anchors_per_condition["B"]["INTENSITY"],
    }

# Plotting unfiltered aggregate_intensities
with log.newSection("Plotting unfiltered aggregate intensities and counts"):
    fig, ax = plt.subplots(2, 1, sharex=True)
    tmp = plt.subplots_adjust(hspace=0.1)
    for sample in range(parameters["SAMPLE_COUNT"]):
        sample_anchors = unfiltered_ions["AGGREGATE_INDEX"][unfiltered_ions["SAMPLE"] == sample]
        frequencies = np.bincount(unfiltered_anchors["ION_COUNT"][sample_anchors])
        start = np.flatnonzero(frequencies > 0)[0]
        tmp = ax[0].plot(
            np.arange(start, len(frequencies)),
            frequencies[start:]
        )
    tmp = ax[0].legend(parameters["APEX_FILE_NAMES"])
    tmp = ax[0].set_ylabel("Frequency")
    tmp = ax[1].boxplot(
        [
            np.log2(
                anchors_per_condition["QC"]["INTENSITY"][
                    unfiltered_anchors["ION_COUNT"] == i
                ]
            ) for i in range(
                1,
                parameters["SAMPLE_COUNT"] + 1
            )
        ],
        whis="range"
    )
    tmp = ax[1].set_ylabel("Average log2(intensity)")
    tmp = ax[1].set_xlabel("Aggregate ion count")
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "aggregate_intensities.pdf", bbox_inches='tight')
    tmp = plt.close()

# LFQ calculations
with log.newSection("Calculating aggregate LFQ"):
    conditions = {
        "A": slice(0, None, 2),
        "B": slice(1, None, 2),
    }
    anchors_per_condition = {}
    anchor_ions_per_condition = {}
    sample_ions = anchor_ions.T.tocsr()
    for condition, condition_slice in conditions.items():
        condition_anchors = np.empty(
            len(anchors),
            dtype=[
                ("INTENSITY", np.float),
                ("CV", np.float),
                ("ION_COUNT", np.int),
            ]
        )
        condition_anchor_ions = sample_ions[condition_slice].T.tocsr()
        anchor_sizes = np.diff(condition_anchor_ions.indptr)
        intensities = np.zeros(len(anchors), dtype=np.float)
        cvs = np.zeros(len(anchors), dtype=np.float)
        for size in range(1, np.max(anchor_sizes) + 1):
            anchor_subset = np.flatnonzero(anchor_sizes == size)
            anchor_subset_matrix = np.concatenate(
                [
                    condition_anchor_ions.data[
                        condition_anchor_ions.indptr[anchor_index]: condition_anchor_ions.indptr[anchor_index + 1]
                    ] for anchor_index in anchor_subset
                ]
            ).reshape(-1, size)
            intensity_matrix = ions["CALIBRATED_INTENSITY"][anchor_subset_matrix]
            cvs[anchor_subset] = scipy.stats.variation(intensity_matrix, axis=1)
            intensities[anchor_subset] = np.average(intensity_matrix, axis=1)
        condition_anchors["ION_COUNT"] = anchor_sizes
        condition_anchors["CV"] = cvs
        condition_anchors["INTENSITY"] = intensities
        anchors_per_condition[condition] = condition_anchors
        anchor_ions_per_condition[condition] = condition_anchor_ions
    logfcs = np.log(anchors_per_condition["A"]["INTENSITY"] / anchors_per_condition["B"]["INTENSITY"])
    anchors_per_condition["QC"] = {
        "INTENSITY": anchors_per_condition["A"]["INTENSITY"] + anchors_per_condition["B"]["INTENSITY"],
    }

# Plotting edge COUNTS
with log.newSection("Plotting aggregate consistent coelution counts"):
    fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [5, 1]})
    tmp = plt.subplots_adjust(hspace=0.1)
    a, b = np.unique(np.diff(neighbors.indptr), return_counts=True)
    tmp = ax[0].scatter(a, np.log(b), marker=".")
    tmp = ax[0].set_ylabel("Log(Aggregate frequency)")
    tmp = ax[1].boxplot(np.diff(neighbors.indptr), whis="range", vert=False, widths=0.5)
    tmp = ax[1].set_yticks([])
    tmp = ax[1].set_xlabel("Edge count")
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "aggregate_edge_counts.pdf", bbox_inches='tight')
    tmp = plt.close()

# Calculating edge accuracy
with log.newSection("Calculating consistent coeluton accuracy"):
    classes = np.round(logfcs)
    a_indices, b_indices = neighbors.nonzero()
    overlap = neighbors.data
    a_classes = classes[a_indices]
    b_classes = classes[b_indices]
    c = np.in1d(a_classes, [-2, 0, 1])
    c &= np.in1d(b_classes, [-2, 0, 1])
    a_classes = a_classes[c]
    b_classes = b_classes[c]
    a_indices = a_indices[c]
    b_indices = b_indices[c]
    overlap = neighbors.data[c]
    eco1 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    yea1 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    hum1 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    eco2 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    yea2 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    hum2 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    for i in range(2, parameters["SAMPLE_COUNT"] + 1):
        full_overlap = overlap == i
        a_full = a_classes[full_overlap]
        b_full = b_classes[full_overlap]
        ecoli, human, yeast = np.unique(np.concatenate([a_full, b_full]), return_counts=True)[1]
        total = ecoli + human + yeast
        ecoli_ratio = ecoli / total
        human_ratio = human / total
        yeast_ratio = yeast / total
        ecoli_expected = ecoli_ratio * ecoli_ratio
        human_expected = human_ratio * human_ratio
        yeast_expected = yeast_ratio * yeast_ratio
        equal = a_full == b_full
        ecoli_hits = equal[(a_full == -2) | (a_full == -2)]
        human_hits = equal[(a_full == 0) | (a_full == 0)]
        yeast_hits = equal[(a_full == 1) | (a_full == 1)]
        ecoli_hit_ratio = 2 * np.sum(ecoli_hits) / total
        yeast_hit_ratio = 2 * np.sum(yeast_hits) / total
        human_hit_ratio = 2 * np.sum(human_hits) / total
        eco1[i] = ecoli_hit_ratio / ecoli_ratio
        yea1[i] = yeast_hit_ratio / yeast_ratio
        hum1[i] = human_hit_ratio / human_ratio
        eco2[i] = ecoli_expected / ecoli_ratio
        yea2[i] = yeast_expected / yeast_ratio
        hum2[i] = human_expected / human_ratio

# Plotting edge accuracy
with log.newSection("Plotting consitent coelution accuracy"):
    e_ex, = plt.plot(eco1[2:], marker="o", linestyle="-", c="blue")
    h_ex, = plt.plot(hum1[2:], marker="o", linestyle="-", c="red")
    y_ex, = plt.plot(yea1[2:], marker="o", linestyle="-", c="green")
    e_th, = plt.plot(eco2[2:], marker=".", linestyle=":", c="blue")
    h_th, = plt.plot(hum2[2:], marker=".", linestyle=":", c="red")
    y_th, = plt.plot(yea2[2:], marker=".", linestyle=":", c="green")
    tmp = plt.legend(
        (
            e_ex,
            h_ex,
            y_ex,
            e_th,
            h_th,
            y_th
        ), (
            "logFC = -2 experimental (Ecoli)",
            "logFC = 0 experimental (Human)",
            "logFC = 1 experimental (Yeast)",
            "logFC = -2 expected (Ecoli)",
            "logFC = 0 expected (Human)",
            "logFC = 1 expected (Yeast)"
        ),
        loc="upper left"
    )
    tmp = plt.xlabel("Consistenly co-eluting samples")
    tmp = plt.xticks(
        np.arange(parameters["SAMPLE_COUNT"] - 1),
        np.arange(2, parameters["SAMPLE_COUNT"] + 1)
    )
    tmp = plt.ylabel("% Consistent co-eluting aggregates with equal logFC")
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "edge_accuracy.pdf", bbox_inches='tight')
    tmp = plt.close()


# Annotation accuracy
with log.newSection("Calculating aggregate annotation accuracy"):
    percolated_annotations = pd.read_csv(
        parameters["PERCOLATOR_TARGET_PIMS"],
        delimiter="\t"
    )
    percolated_fdrs = percolated_annotations.values[:, 2]
    percolated_anchors, percolated_peptides = anchor_peptide_match_counts.nonzero()
    significant_pims = percolated_annotations.values[
        percolated_fdrs <= 0.01,
        0
    ].astype(int)
    significant_anchors = percolated_anchors[significant_pims]
    significant_peptides = percolated_peptides[significant_pims]
    significant_proteins = src.peptides.getProteinAccessions(
        significant_peptides,
        peptides,
        proteins
    )
    significant_organisms = np.array(
        [
            prot.split("_")[-1] for prot in significant_proteins
        ]
    )

# Plot annotation accuracies
with log.newSection("Plotting aggreggate annotation accuracy"):
    # organism_coloring = np.repeat("grey", len(significant_organisms)).astype(np.object)
    # organism_coloring[significant_organisms == "HUMAN"] = "red"
    # organism_coloring[significant_organisms == "ECOLI"] = "blue"
    # organism_coloring[significant_organisms == "YEAST"] = "green"
    # significant_organisms[organism_coloring == "grey"] = "OTHER"
    # unique_organisms = ["OTHER", "ECOLI", "YEAST", "HUMAN"]
    unique_organisms = ["ECOLI", "YEAST", "HUMAN"]
    fig, ax = plt.subplots(1, len(unique_organisms), sharex=True, sharey=True)
    for i, c in enumerate(unique_organisms):
        scatter = ax[i].scatter(
            np.log(anchors_per_condition["QC"]["INTENSITY"][significant_anchors]),
            logfcs[significant_anchors],
            c="grey",
            s=1
        )
        s = np.flatnonzero(significant_organisms == c)
        scatter = ax[i].scatter(
            np.log(anchors_per_condition["QC"]["INTENSITY"][significant_anchors[s]]),
            logfcs[significant_anchors][s],
            c="red",
            s=1
        )
        tmp = ax[i].set_title(c)
    tmp = ax[1].set_xlabel("Log(QC)")
    tmp = ax[0].set_ylabel("Log(A/B)")
    # tmp = fig.text(0.5, 0.04, "Log(QC)", ha='center')
    # tmp = fig.text(0.04, 0.5, "Log(A/B)", va='center', rotation='vertical')
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "annotation_accuracy.pdf", bbox_inches='tight')
    tmp = plt.close()
