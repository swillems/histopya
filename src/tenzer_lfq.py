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
import matplolib as sns
from sklearn import linear_model

# Initializing
parameter_file_name = "data/tenzer/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "_interactive.txt")

# Loading data
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

# # Plot calibration
# with log.newSection("Plotting calibration results"):
#     estimation_anchors = src.io.loadArray(
#         "PSEUDO_AGGREGATE_IONS_FILE_NAME",
#         parameters,
#         log
#     )[1::2]
#     fig, ax = plt.subplots(3, 2, sharex=True, sharey="row")
#     tmp = plt.subplots_adjust(hspace=0.1, wspace=0.1)
#     quick_anchors = np.empty(
#         len(estimation_anchors),
#         dtype=estimation_anchors.dtype
#     )
#     for attribute in quick_anchors.dtype.names:
#         quick_anchors[attribute] = np.average(estimation_anchors[attribute], axis=1)
#     for attribute, x_loc, y_loc in [
#         ("MZ", 0, 0),
#         ("CALIBRATED_MZ", 1, 0),
#         ("DT", 0, 1),
#         ("CALIBRATED_DT", 1, 1),
#         ("RT", 0, 2),
#         ("CALIBRATED_RT", 1, 2),
#     ]:
#         average_attribute = quick_anchors[attribute].reshape((-1, 1))
#         data = estimation_anchors[attribute] - average_attribute
#         if attribute in parameters["RELATIVE_ATTRIBUTES"]:
#             data *= 1000000 / average_attribute
#             label = "{} PPM".format(attribute)
#         else:
#             label = attribute
#         tmp = sns.violinplot(
#             x='variable',
#             y='value',
#             data=pd.melt(
#                 pd.DataFrame(
#                     data
#                 ),
#             ),
#             ax=ax[y_loc, x_loc],
#             inner=None
#         )
#         tmp = ax[y_loc, x_loc].axes.get_xaxis().set_visible(False)
#         tmp = ax[y_loc, x_loc].axes.get_yaxis().set_visible(False)
#         tmp = ax[y_loc, x_loc].set_xlabel("")
#         tmp = ax[y_loc, x_loc].set_ylabel("")
#     tmp = ax[0, 0].set_ylabel(u"Δ MZ (ppm)")
#     tmp = ax[1, 0].set_ylabel(u"Δ DT")
#     tmp = ax[2, 0].set_ylabel(u"Δ RT")
#     tmp = ax[0, 0].axes.get_yaxis().set_visible(True)
#     tmp = ax[1, 0].axes.get_yaxis().set_visible(True)
#     tmp = ax[2, 0].axes.get_yaxis().set_visible(True)
#     tmp = ax[0, 0].set_title("Uncalibrated")
#     tmp = ax[0, 1].set_title("Calibrated")
#     tmp = ax[2, 0].axes.get_xaxis().set_visible(True)
#     tmp = ax[2, 1].axes.get_xaxis().set_visible(True)
#     # tmp = fig.text(0.5, 0.04, 'Sample', ha='center')
#     tmp = ax[2, 0].get_xaxis().set_ticklabels(
#         [
#             "_".join(x.split("/")[-1].split("_")[:-1]) for x in parameters["APEX_FILE_NAMES"]
#         ]
#     )
#     tmp = ax[2, 1].get_xaxis().set_ticklabels(
#         [
#             "_".join(x.split("/")[-1].split("_")[:-1]) for x in parameters["APEX_FILE_NAMES"]
#         ]
#     )
#     tmp = plt.setp(ax[2, 0].get_xticklabels(), rotation=45, ha="right")
#     tmp = plt.setp(ax[2, 1].get_xticklabels(), rotation=45, ha="right")
#     # plt.show()
#     tmp = plt.savefig(parameters["PLOTS_PATH"] + "calibration.pdf", bbox_inches='tight')
#     tmp = plt.close()

# Plot calibration
with log.newSection("Plotting calibration results"):
    estimation_anchors = src.io.loadArray(
        "PSEUDO_AGGREGATE_IONS_FILE_NAME",
        parameters,
        log
    )[1::2]
    fig, ax = plt.subplots(2, 3, sharex="col", sharey=True)
    tmp = plt.subplots_adjust(hspace=0.1, wspace=0.1)
    quick_anchors = np.empty(
        len(estimation_anchors),
        dtype=estimation_anchors.dtype
    )
    for attribute in quick_anchors.dtype.names:
        quick_anchors[attribute] = np.average(estimation_anchors[attribute], axis=1)
    for attribute, x_loc, y_loc in [
        ("MZ", 0, 0),
        ("CALIBRATED_MZ", 0, 1),
        ("DT", 1, 0),
        ("CALIBRATED_DT", 1, 1),
        ("RT", 2, 0),
        ("CALIBRATED_RT", 2, 1),
    ]:
        average_attribute = quick_anchors[attribute].reshape((-1, 1))
        data = estimation_anchors[attribute] - average_attribute
        if attribute in parameters["RELATIVE_ATTRIBUTES"]:
            data *= 1000000 / average_attribute
            label = "{} PPM".format(attribute)
        else:
            label = attribute
        tmp = sns.violinplot(
            x='value',
            y='variable',
            data=pd.melt(
                pd.DataFrame(
                    data
                ),
            ),
            ax=ax[y_loc, x_loc],
            inner=None,
            orient="h"
        )
        tmp = ax[y_loc, x_loc].axes.get_xaxis().set_visible(False)
        tmp = ax[y_loc, x_loc].axes.get_yaxis().set_visible(False)
        tmp = ax[y_loc, x_loc].set_xlabel("")
        tmp = ax[y_loc, x_loc].set_ylabel("")
    tmp = ax[0, 0].set_title(u"Δ MZ (ppm)")
    tmp = ax[0, 1].set_title(u"Δ DT")
    tmp = ax[0, 2].set_title(u"Δ RT")
    tmp = ax[1, 0].axes.get_xaxis().set_visible(True)
    tmp = ax[1, 1].axes.get_xaxis().set_visible(True)
    tmp = ax[1, 2].axes.get_xaxis().set_visible(True)
    tmp = ax[0, 2].yaxis.set_label_position("right")
    tmp = ax[1, 2].yaxis.set_label_position("right")
    tmp = ax[0, 2].set_ylabel("Uncalibrated")
    tmp = ax[1, 2].set_ylabel("Calibrated")
    tmp = ax[0, 2].axes.get_yaxis().set_visible(True)
    tmp = ax[1, 2].axes.get_yaxis().set_visible(True)
    tmp = ax[0, 2].tick_params(axis='y', which='both', left=False)
    tmp = ax[1, 2].tick_params(axis='y', which='both', left=False)
    tmp = ax[0, 0].get_yaxis().set_ticklabels(
        [
            "_".join(x.split("/")[-1].split("_")[:-1]) for x in parameters["APEX_FILE_NAMES"]
        ]
    )
    tmp = ax[1, 0].get_yaxis().set_ticklabels(
        [
            "_".join(x.split("/")[-1].split("_")[:-1]) for x in parameters["APEX_FILE_NAMES"]
        ]
    )
    tmp = ax[0, 0].axes.get_yaxis().set_visible(True)
    tmp = ax[1, 0].axes.get_yaxis().set_visible(True)
    # plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "calibration.pdf", bbox_inches='tight')
    tmp = plt.close()

# Unfiltered LFQ calculations
with log.newSection("Calculating unfiltered aggregate LFQ"):
    # parameters["UNFILTERED_IONS_FILE_NAME"] = "data/tenzer/results/ions_unfiltered_tmp.npy"
    unfiltered_ions = src.io.loadArray("IONS_UNFILTERED_FILE_NAME", parameters)
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

# Plotting unfiltered aggregate intensities and counts
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
    tmp = ax[0].legend(
        [
            "_".join(x.split("/")[-1].split("_")[:-1]) for x in parameters["APEX_FILE_NAMES"]
        ],
        # bbox_to_anchor=(0., 1.02, 1., .102),
        # loc=5,
        fontsize="small",
        ncol=2,
        # mode="expand",
        # borderaxespad=0.
    )
    tmp = ax[0].set_ylabel("Frequency")
    x = [
        np.log2(
            anchors_per_condition["QC"]["INTENSITY"][
                unfiltered_anchors["ION_COUNT"] == i
            ]
        ) for i in range(
            1,
            parameters["SAMPLE_COUNT"] + 1
        )
    ]
    # x[0][0]=0 # inf present
    # d = pd.DataFrame(
    #     np.stack(
    #         [
    #             np.concatenate(x),
    #             np.repeat(np.arange(parameters["SAMPLE_COUNT"]) + 1, [len(i) for i in x]),
    #         ]
    #     ).T,
    #     columns=["QC_INTENSITY", "ION_COUNT"]
    # )
    # tmp = sns.violinplot(
    #     x='ION_COUNT',
    #     y='QC_INTENSITY',
    #     data=d,
    #     inner="quartile",
    #     # gridsize=1000,
    # )
    tmp = ax[1].boxplot(x, whis="range")
    tmp = ax[1].set_ylabel("Average log2(intensity)")
    tmp = ax[1].set_xlabel("Aggregate ion count")
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "aggregate_counts_and_intensities.pdf", bbox_inches='tight')
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
    percolated_pims = pd.read_csv(
        parameters["PERCOLATOR_TARGET_PIMS"],
        delimiter="\t"
    )
    percolated_pim_fdrs = percolated_pims.values[:, 2]
    anchor_pims, peptide_pims = anchor_peptide_match_counts.nonzero()
    pim_fdr = 0.01
    significant_percolated_pims = percolated_pims.values[
        percolated_pim_fdrs <= pim_fdr,
        0
    ].astype(int)
    significant_anchors = anchor_pims[significant_percolated_pims]
    significant_peptides = peptide_pims[significant_percolated_pims]
    significant_peptide_sequences = src.peptides.getPeptideSequences(
        significant_peptides,
        peptides,
        peptide_index_matrix,
        total_protein_sequence
    )
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
    organism_coloring = np.repeat("grey", len(significant_organisms)).astype(np.object)
    organism_coloring[significant_organisms == "HUMAN"] = "red"
    organism_coloring[significant_organisms == "ECOLI"] = "blue"
    organism_coloring[significant_organisms == "YEAST"] = "green"
    significant_organisms[organism_coloring == "grey"] = "OTHER"
    unique_organisms = ["OTHER", "ECOLI", "YEAST", "HUMAN"]
    for i, c in enumerate(unique_organisms):
        s = np.flatnonzero(significant_organisms == c)
        scatter = plt.scatter(
            np.log(anchors_per_condition["QC"]["INTENSITY"][significant_anchors[s]]),
            logfcs[significant_anchors][s],
            c=organism_coloring[s],
            s=1
        )
    tmp = plt.xlabel("Log(A+B)")
    tmp = plt.ylabel("Log(A/B)")
    tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "annotation_accuracy.pdf", bbox_inches='tight')
    tmp = plt.close()

# # Plot annotation accuracies
# with log.newSection("Plotting aggreggate annotation accuracy"):
#     unique_organisms = ["ECOLI", "YEAST", "HUMAN"]
#     fig, ax = plt.subplots(1, len(unique_organisms), sharex=True, sharey=True)
#     for i, c in enumerate(unique_organisms):
#         scatter = ax[i].scatter(
#             np.log(anchors_per_condition["QC"]["INTENSITY"][significant_anchors]),
#             logfcs[significant_anchors],
#             c="grey",
#             s=1
#         )
#         s = np.flatnonzero(significant_organisms == c)
#         scatter = ax[i].scatter(
#             np.log(anchors_per_condition["QC"]["INTENSITY"][significant_anchors[s]]),
#             logfcs[significant_anchors][s],
#             c="red",
#             s=1
#         )
#         tmp = ax[i].set_title(c)
#     tmp = ax[1].set_xlabel("Log(A+B)")
#     tmp = ax[0].set_ylabel("Log(A/B)")
#     # tmp = plt.show()
#     tmp = plt.savefig(parameters["PLOTS_PATH"] + "annotation_accuracy.pdf", bbox_inches='tight')
#     tmp = plt.close()

# Plot significant aggreggate ion count enrichment
with log.newSection("Plotting significant aggreggate ion count enrichment"):
    counts = np.bincount(anchors["ION_COUNT"][significant_anchors])
    tmp = plt.plot(np.arange(2, len(counts)), counts[2:], marker=".")
    tmp = plt.xlabel("Significant aggregate ion counts")
    tmp = plt.ylabel("Frequency")
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "annotated_aggregate_enrichment.pdf", bbox_inches='tight')
    tmp = plt.close()

# Plot annotation example
with log.newSection("Plotting significant aggreggate ion annotation example"):
    top_score_index = np.argmax(anchor_peptide_scores.data)
    anchor_index = np.searchsorted(
        anchor_peptide_scores.indptr,
        top_score_index,
        "right"
    ) - 1
    anchor_candidates = fragment_peptide_indices[
        slice(*anchor_boundaries[anchor_index])
    ]
    anchor_neighbors = neighbors.indices[
        neighbors.indptr[anchor_index]: neighbors.indptr[anchor_index + 1]
    ]
    raw_neighbor_candidates = np.concatenate(
        [
            fragment_peptide_indices[
                slice(*anchor_boundaries[n])
            ] for n in anchor_neighbors
        ]
    )
    neighbor_candidates = raw_neighbor_candidates[
        np.isin(
            raw_neighbor_candidates,
            anchor_candidates
        )
    ]
    candidates, candidate_counts = np.unique(
        neighbor_candidates,
        return_counts=True
    )
    counts, frequency = np.unique(
        candidate_counts,
        return_counts=True
    )
    counts = np.concatenate([[0], counts])
    frequency = np.concatenate(
        [
            [len(anchor_candidates)],
            np.cumsum(frequency[::-1])[::-1]
        ]
    )
    ransac = linear_model.RANSACRegressor()
    tmp = ransac.fit(
        counts.reshape(-1, 1)[:-1],
        np.log(frequency).reshape(-1, 1)[:-1]
    )
    score = -ransac.predict(counts[-1])[0][0]
    # plt.bar(counts, frequency)
    # plt.scatter(counts, frequency, marker=".")
    # plt.xlabel("#Peptides with fragment count")
    # plt.ylabel("Frequency")
    # plt.show()
    tmp = plt.scatter(counts, np.log(frequency), marker=".")
    tmp = plt.plot(
        [
            0,
            counts[-1]
        ],
        [
            ransac.predict(counts[0])[0][0],
            ransac.predict(counts[-1])[0][0]
        ]
    )
    tmp = plt.xlabel("#Peptides with fragment count")
    tmp = plt.ylabel("Log frequency")
    tmp = plt.title("Score={}".format(score))
    # plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "annotation_example.pdf", bbox_inches='tight')
    tmp = plt.close()
