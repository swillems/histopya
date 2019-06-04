#!venv/bin/python


import src.parameters
import src.io
import src.ions
import src.aggregates
import src.peptides
import pandas as pd
import numpy as np
import src.parallelization as mp
import scipy.sparse
from collections import defaultdict
import json
import sys
import time
import traceback as tb
from scipy.special import binom as binom
from contextlib import contextmanager
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn import linear_model
import sklearn
import csv
# import importlib
# importlib.reload(src.aggregates)
# parameter_file_name = "data/lfq_single_A_B_QC/parameters_res_auto.json"
parameter_file_name = "data/test2/parameters.json"
# parameter_file_name = "data/lfq_swim_190327/parameters_manual.json"
# parameter_file_name = "data/lfq_udmse_190327/parameters_default_peakpicking_test.json"
# parameter_file_name = "data/k562_3_samples/parameters_res_25000.json"
# parameter_file_name = "data/comparison_k562_lfq/parameters_k562.json"
# parameter_file_name = "data/comparison_k562_lfq/parameters_lfq.json"
# parameters = src.parameters.importParameterDictFromJSON("data/test/parameters.json")
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
# parameters = src.parameters.importParameterDictFromJSON("data/pharmafluidics/parameters_30_min.json"); log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "_interactive.txt")
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "_interactive.txt")
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
anchor_peptide_scores = src.io.loadMatrix(
    "ANCHOR_PEPTIDE_SCORES_FILE_NAME",
    parameters,
)
anchor_peptide_match_counts = src.io.loadMatrix(
    "ANCHOR_PEPTIDE_MATCH_COUNTS_FILE_NAME",
    parameters,
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
parameters["PERCOLATOR_DATA_FILE_NAME"] = parameters["PERCOLATOR_DATA_FILE_NAME"][:-4] + "_interactive.csv"
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
    log,
)



















sample_types = "tenzer"

if sample_types == "in_house":
    conditions = {
        "A": slice(None, 9),
        "B": slice(9, 18),
        "QC": slice(18, None)
    }
elif sample_types == "test":
    conditions = {
        "A": 0,
        "B": 1,
        "QC": 2
    }
elif sample_types == "tenzer":
    conditions = {
        "A": slice(0, None, 2),
        "B": slice(1, None, 2),
    }
elif sample_types == "600ms":
    conditions = {
        "A": [0,1,4,5],
        "B": [2,3,6,7],
    }
elif sample_types == "swath":
    conditions = {
        "A": slice(None, 3),
        "B": slice(3, 6),
        "QC": slice(6, None)
    }


for condition in conditions.keys():
    condition_anchors = anchors_per_condition[condition]
    size_range = list(range(0, np.max(condition_anchors["ION_COUNT"]) + 1))
    anchor_subsets = [
        np.flatnonzero(
            condition_anchors["ION_COUNT"] == size
        ) for size in size_range
    ]
    cv_percentile = np.array(
        [
            np.percentile(
                condition_anchors["CV"][anchor_subsets[i]],
                range(101)
            ) for i in size_range
        ]
    )
    for cvs in cv_percentile:
        plt.plot(list(range(101))[::-1], cvs[::-1])
    plt.legend(
        [
            "Anchors found in {} samples of condition {} ({})".format(
                size,
                condition,
                len(anchor_subsets[size])
            ) for size in size_range
        ]
    )
    plt.xlabel("% of data".format(condition))
    plt.ylabel("Maximum CV of calibrated intensity")
    plt.show()

























#
# a_int = anchors["CALIBRATED_LOGINT"][a_indices]
# b_int = anchors["CALIBRATED_LOGINT"][b_indices]

#
#
# import statsmodels.stats
#
# x = np.flatnonzero(
#  (
#      ~np.isnan(diffs)
#  # ) & (
#  #     anchors["CALIBRATED_LOGINT"] > 13
#  ) & (
#       anchors["ION_COUNT"]==parameters["SAMPLE_COUNT"]
#   ) & (
#       anchors["LE"]
#  )
# )
# plt.show(plt.plot(*np.unique(np.round(diffs[x], 1), return_counts=True)))
# plt.show(sns.jointplot(means_a[x], diffs[x], kind="hex"))
# plt.show(sns.jointplot(means_a[x][::100], diffs[x][::100], kind="kde"))
#
# low = -.5
# upp = .5
# results = []
# for i in np.flatnonzero(anchors["ION_COUNT"] == 9):
#     z = a[i].data
#     x1 = z[condition_a]
#     x2 = z[condition_b]
#     # np.mean(x1), x1
#     # np.mean(x2), x2
#     result = statsmodels.stats.weightstats.ttost_ind(
#         x1,
#         x2,
#         low,
#         upp,
#         usevar='unequal',
#     )
#     results.append(result[0])
#     # result
#     # input()
#
#




def getIsotopicNeighbors(parameters, anchors, neighbors):
    PPM = 10
    isotopic_mass_delta = parameters["ISOTOPIC_DELTA_MASS"]
    max_charge = parameters["MAX_CHARGE"]
    # TODO include above in parameters?
    neighbors = scipy.sparse.load_npz(parameters["ANCHOR_NEIGHBORS_FILE_NAME"])
    neighbors += neighbors.T
    a_ind, b_ind = neighbors.nonzero()
    a_mass = anchors["CALIBRATED_MZ"][a_ind]
    b_mass = anchors["CALIBRATED_MZ"][b_ind]
    delta_mass = b_mass - a_mass
    positive = delta_mass > 0
    a_ind = a_ind[positive]
    b_ind = b_ind[positive]
    a_mass = a_mass[positive]
    b_mass = b_mass[positive]
    delta_mass = delta_mass[positive]
    # TODO include mismatches?
    z_indices = {
        i: abs(
            delta_mass * i - isotopic_mass_delta
        ) * 1000000 < PPM * a_mass for i in range(1, max_charge + 1)
    }
    z_indices2 = {
        i: abs(
            delta_mass * i - 2 * isotopic_mass_delta
        ) * 1000000 < PPM * a_mass for i in range(1, max_charge + 1)
    }
    # z_indices3 = {
    #     i: abs(
    #         delta_mass * i - 3 * isotopic_mass_delta
    #     ) * 1000000 < PPM * a_mass for i in range(1, max_charge + 1)
    # }
    # z_indices4 = {
    #     i: abs(
    #         delta_mass * i - 4 * isotopic_mass_delta
    #     ) * 1000000 < PPM * a_mass for i in range(1, max_charge + 1)
    # }
    isotopic_pairs = np.logical_or.reduce(
        list(z_indices.values()) + list(z_indices2.values())
    )
    isotopic_neighbors = scipy.sparse.csr_matrix(
        (
            isotopic_pairs,
            (
                a_ind,
                b_ind,
            )
        ),
        shape=neighbors.shape,
        dtype=bool
    )
    isotopic_neighbors.eliminate_zeros()
    return isotopic_neighbors


def getIsotopes(parameters, isotopic_neighbors):
    isotopic_neighbors += isotopic_neighbors.T
    isotope_count, anchor_labels = scipy.sparse.csgraph.connected_components(
        isotopic_neighbors,
        directed=False,
        return_labels=True
    )
    anchor_order = np.argsort(anchor_labels)
    anchor_label_breaks = np.concatenate(
        [
            [0],
            np.flatnonzero(np.diff(anchor_labels[anchor_order]) > 0) + 1,
            [len(anchor_labels)]
        ]
    )
    isotopes = []
    for i in np.flatnonzero(np.diff(anchor_label_breaks) > 1):
        isotope = anchor_order[
            anchor_label_breaks[i]: anchor_label_breaks[i + 1]
        ]
        isotopes.append(isotope)
    return isotopes


def setIsotopes(anchors, isotopes, parameters):
    # TODO SPECTRUM should be ISOTOPE
    log = parameters["LOG"]
    log.printMessage("Setting anchor isotopes")
    defined = len(np.concatenate(isotopes))
    undefined = len(anchors) - defined
    for isotope_index, isotope_anchors in enumerate(isotopes):
        # TODO +1 for isotopes only?
        anchors["SPECTRUM"][isotope_anchors] = isotope_index + undefined
    anchors["SPECTRUM"][anchors["SPECTRUM"] == 0] = np.arange(undefined)
    # log.printMessage("Saving anchor spectra")
    # np.save(parameters["ANCHOR_SPECTRUM_FILE_NAME"], anchors)
    return anchors


def defineSpectra(anchors, parameters):
    log = parameters["LOG"]
    log.printMessage("Defining isotopes")
    isotope_attributes = scipy.sparse.csr_matrix(
        (
            np.ones(len(anchors)),
            (anchors["SPECTRUM"], np.arange(len(anchors))),
        ),
        dtype=np.float
    )
    isotopes = np.empty(
        isotope_attributes.shape[0],
        dtype=[
            # Floats
            ("CALIBRATED_RT", np.float),
            ("CALIBRATED_DT", np.float),
            ("CALIBRATED_MASS", np.float),
            # Boolean
            ("LE", np.bool),
            # Integer
            ("ANCHOR_COUNT", np.int),
            ("CHARGE", np.int),
            # Pointers
            ("ANALYTE", np.int),
        ]
    )
    isotope_order = np.argsort(anchors["SPECTRUM"])
    log.printMessage("Calculating spectrum properties")
    isotopes["ANCHOR_COUNT"] = np.diff(isotope_attributes.indptr)
    isotopes["ANALYTE"] = 0
    isotopic_mass_delta = parameters["ISOTOPIC_DELTA_MASS"]
    max_charge = parameters["MAX_CHARGE"]
    isotopic_masses = np.zeros(len(isotopes), dtype=np.float)
    isotopic_charges = np.zeros(len(isotopes), dtype=np.int)
    for isotope_index in np.flatnonzero(isotopes["ANCHOR_COUNT"] > 1):
        isotopic_anchors = isotope_attributes.indices[
            slice(
                *isotope_attributes.indptr[isotope_index: isotope_index + 2]
            )
        ]
        mzs = np.sort(anchors["CALIBRATED_MZ"][isotopic_anchors])
        deltas = np.diff(mzs)
        charges = np.round(isotopic_mass_delta / deltas)
        isotopic_masses[isotope_index] = mzs[0]
        try:
            isotopic_charges[isotope_index] = np.max(charges[charges <= max_charge])
        except ValueError:
            isotopic_charges[isotope_index] = 0
    isotopes["CHARGE"] = isotopic_charges
    proton_mass = parameters["PROTON_MASS"]
    isotopes["CALIBRATED_MASS"] = isotopic_masses * isotopic_charges - isotopic_charges * proton_mass
    first_anchor_indices = isotope_order[np.cumsum(isotopes["ANCHOR_COUNT"]) - 1]
    isotopes["LE"] = anchors["LE"][first_anchor_indices]
    for attribute in [
        "CALIBRATED_RT",
        "CALIBRATED_DT",
    ]:
        isotope_attributes.data = anchors[attribute][isotope_order]
        attribute_row_sums = isotope_attributes.sum(axis=1).A1
        isotopes[attribute] = attribute_row_sums / isotopes["ANCHOR_COUNT"]
    isotope_anchors = scipy.sparse.csr_matrix(
        (
            np.ones(len(anchors)),
            (anchors["SPECTRUM"], np.arange(len(anchors))),
        ),
        dtype=np.bool
    )
    # log.printMessage("Saving isotopes")
    # np.save(parameters["ISOTOPES_FILE_NAME"], isotopes)
    # log.printMessage("Saving isotope anchors")
    # scipy.sparse.save_npz(parameters["ISOTOPE_ANCHORS_FILE_NAME"], isotope_anchors)
    return isotopes, isotope_anchors


def getIsotopeNeighbors(anchors, neighbors, isotopes):
    a_ind, b_ind = neighbors.nonzero()
    l = np.stack([anchors["SPECTRUM"][a_ind], anchors["SPECTRUM"][b_ind]])
    a, b = np.unique(l, return_counts=True, axis=1)
    # isotope_neighbors = scipy.sparse.csr_matrix(
    #     (
    #         b,
    #         a
    #     ),
    #     shape=(len(isotopes), len(isotopes)),
    #     dtype=np.int
    # )
    z = isotopes["ANCHOR_COUNT"][a]
    x = np.flatnonzero(
        (
            z[1,:] > 1
        )&(
            z[0,:] > 1
        )&(
            b > 1
        )
    )
    # j=17;i=x[j];b[i];isotopes["ANCHOR_COUNT"][a[:,i]]
    i = isotopes["ANCHOR_COUNT"][a[:,x]]
    # TODO change tolerance?
    j = np.multiply.reduce(i, axis=0) - np.maximum.reduce(i, axis=0)
    k = np.flatnonzero(b[x] >= j)
    isotope_neighbors = scipy.sparse.csr_matrix(
        (
            a[0, x[k]] != a[1, x[k]],
            a[:, x[k]]
        ),
        shape=(len(isotopes), len(isotopes)),
        dtype=np.bool
    )
    # isotope_neighbors.setdiag(0)
    isotope_neighbors.eliminate_zeros()
    return isotope_neighbors


isotopic_neighbors = getIsotopicNeighbors(parameters, anchors, neighbors)
isotopes = getIsotopes(parameters, isotopic_neighbors)
anchors = setIsotopes(anchors, isotopes, parameters)
isotopes, isotope_anchors = defineSpectra(anchors, parameters)
isotope_neighbors = getIsotopeNeighbors(anchors, neighbors, isotopes)







potential_mono_ind = a_ind[np.logical_or.reduce(list(z_indices.values()))]
non_mono_ind = b_ind[np.logical_or.reduce(list(z_indices.values()))]
monos = ~np.in1d(potential_mono_ind, non_mono_ind)
mono_ind = np.unique(potential_mono_ind[monos])

anchor_z = np.zeros(len(anchors), dtype=int)

for z in range(1, max_charge + 1):
    anchor_z[a_ind[z_indices[z]]] = z
    anchor_z[b_ind[z_indices[z]]] = z

anchor_mono = np.zeros(len(anchors), dtype=int)
anchor_mono[mono_ind] = 1


anchor_mass = (
    (
        anchor_z * anchors["CALIBRATED_MZ"] - anchor_z
    ) - isotopic_mass_delta * (1 - anchor_mono)
) * (anchor_mono > 0)
#
# plt.show(
#     plt.plot(
#         *np.unique(
#             np.round(
#                 delta_mass[
#                     np.in1d(
#                         a_ind,
#                         np.flatnonzero(
#                             (anchor_z * anchor_mono) == 1
#                         )
#                     ) & np.in1d(
#                         b_ind,
#                         np.flatnonzero(
#                             (anchor_z * anchor_mono) == 1
#                         )
#                     )
#                 ],
#                 4
#             ),
#             return_counts=True
#         )
#     )
# )

#
#
























































# aas = src.io.loadJSON("AMINO_ACID_FILE_NAME", parameters)
# aas["C"] += parameters["FIXED_MODIFICATIONS"]["C"]
# aas["K"] += parameters["FIXED_MODIFICATIONS"]["K"]
# aas["O"] = aas["K"]
# aas["U"] = aas["K"]
# del aas["K"]
# del aas["U"]
# del aas["$"]
#
# # mass_aa_dict = {round(mass, 1): aa for aa, mass in aas.items()}
# #
# # mass_aa_dict[round(aas["I"], 1)] = "[I/L]"
# dipeptide_aas = {}
# # a, b = list(zip(*aas.items()))
# for i, (aa1, mass1) in enumerate(aas.items()):
#     for aa2, mass2 in aas.items():
#         mass = mass1 + mass2
#         aa = "{}{}".format(aa1, aa2)
#         dipeptide_aas[aa] = mass
#
# # dipeptide_mass_aa_dict.update(mass_aa_dict)
#
# mass_aa_dict = {}
# for aa, mass in aas.items():
#     mass = round(mass, 1)
#     if mass in mass_aa_dict:
#         mass_aa_dict[mass] = "{}/{}".format(mass_aa_dict[mass], aa)
#     else:
#         mass_aa_dict[mass] = aa
#
#
# for aa, mass in dipeptide_aas.items():
#     mass = round(mass, 1)
#     if mass in mass_aa_dict:
#         mass_aa_dict[mass] = "{}/{}".format(mass_aa_dict[mass], aa)
#     else:
#         mass_aa_dict[mass] = aa
#
#
# for m, aa in mass_aa_dict.items():
#     if len(aa) > 1:
#         mass_aa_dict[m] = "[{}]".format(aa)
#
# first, second = neighbors.nonzero()
# delta_masses = anchors["MZ"][second] - anchors["MZ"][first]
# #
# # dm, dm_freq = np.unique(np.round(delta_masses, 3), return_counts=True)
# # tmp = plt.plot(dm, dm_freq, marker=".")
# # for i, j, k in [(m, dm_freq[np.flatnonzero(m == dm)[0]], aa) if m in dm else (m, 0, aa) for m, aa in mass_aa_dict.items()]:
# #     tmp = plt.text(i, j, k, rotation=45)
# #
# #
# # # locs, labels = plt.xticks(list(mass_aa_dict.keys()), list(mass_aa_dict.values()))
# # # tmp = plt.setp(labels, rotation=90)
# # tmp = plt.show()
#
# delta_neighbors = neighbors.copy()
# delta_neighbors.data = np.round(delta_masses, 1)
# explainable_neighbors = neighbors.copy()
# explainable_neighbors.data = np.isin(
#     np.round(delta_masses, 1),
#     np.round([mass for mass in mass_aa_dict], 1)
# )
# explainable_neighbors.eliminate_zeros()
# explainable_neighbor_counts = np.sum(explainable_neighbors.A, axis=1) > 0
#
#
#
# def extendSequence(
#     index,
#     delta_neighbors,
#     explainable_neighbors,
#     mass_aa_dict,
#     aa_count=0,
#     non_aa_count=0,
#     max_non_aa_count=0,
#     potentials=None,
#     min_aa_count=5,
#     require_clique=True
# ):
#     new_potentials = delta_neighbors.indices[
#         delta_neighbors.indptr[index]: delta_neighbors.indptr[index + 1]
#     ]
#     if require_clique and (potentials is not None):
#         new_potentials = new_potentials[np.isin(new_potentials, potentials)]
#     if np.sum(explainable_neighbor_counts[new_potentials]) + aa_count < min_aa_count:
#         return
#     elif len(new_potentials) == 0:
#         yield ""
#     for new_index in new_potentials:
#         mass = delta_neighbors[index, new_index]
#         if explainable_neighbors[index, new_index]:
#             aa = mass_aa_dict[mass]
#             aa_count += 1
#         elif non_aa_count < max_non_aa_count:
#             aa = "[{}]".format(mass, 1)
#             non_aa_count += 1
#         else:
#             continue
#         for extension in extendSequence(
#             index=new_index,
#             delta_neighbors=delta_neighbors,
#             explainable_neighbors=explainable_neighbors,
#             mass_aa_dict=mass_aa_dict,
#             aa_count=aa_count,
#             non_aa_count=non_aa_count,
#             max_non_aa_count=max_non_aa_count,
#             potentials=new_potentials,
#             min_aa_count=min_aa_count,
#             require_clique=require_clique
#         ):
#             yield aa + extension
#         if explainable_neighbors[index, new_index]:
#             aa_count -= 1
#         else:
#             non_aa_count -= 1
#
#
#
# sequence_tags = []
# for index in range(delta_neighbors.shape[0]):
# # for index in range(100):
#     mass = anchors["MZ"][index]
#     sequence_tag_start = "[{}]".format(np.round(mass, 1))
#     for sequence_tag in extendSequence(
#         index,
#         delta_neighbors,
#         explainable_neighbors,
#         mass_aa_dict,
#         max_non_aa_count=0,
#         min_aa_count=10
#     ):
#         print(sequence_tag)
#         input(sequence_tag)
#         sequence_tags.append(sequence_tag)
#
#
# sequence_tags = np.unique(sequence_tags)
# # sequence_tags = np.array([s[::-1] for s in sequence_tags])
# for i in sorted(sequence_tags): i















# def plotIntraRunDistances():
quick_isotopic_pairs = src.io.loadArray(
    "PSEUDO_ISOTOPIC_PAIRS_FILE_NAME",
    parameters,
    log
)
first_ion_set = anchor_ions[quick_isotopic_pairs[:, 0]].toarray()
second_ion_set = anchor_ions[quick_isotopic_pairs[:, 1]].toarray()
matched_isotopes = ions[np.stack([first_ion_set, second_ion_set], axis=1)]
alignment_parameters = {}
for attribute in [
    "DT",
    "RT",
]:
    ptps = np.ptp(matched_isotopes[attribute], axis=1)
    if attribute in parameters["RELATIVE_ATTRIBUTES"]:
        ptps *= 1000000 / np.min(matched_isotopes[attribute], axis=1)
        ylabel = "Maximum {} distance (ppm) between pseudo isotopic pairs".format(attribute)
    else:
        ylabel = "Maximum {} distance between pseudo isotopic pairs".format(attribute)
    ptp_distribution = np.percentile(ptps, range(101), axis=0)
    plt.show(
        plt.plot(ptp_distribution) + [
            plt.axvline(100 * parameters["ANCHOR_ALIGNMENT_PERCENTILE_THRESHOLD"], color="black"),
            plt.legend(
                [
                    "Sample {}".format(i) for i in range(parameters["SAMPLE_COUNT"])
                ]
            ),
            plt.xlabel("% of data"),
            plt.ylabel(ylabel)
        ]
    )

















































































# qc_intensity = anchors_per_condition["A"]["INTENSITY"] + anchors_per_condition["B"]["INTENSITY"]
# qc_intensity = qc_intensity[ans]
# organism_logfcs = logfcs[ans]
order = np.argsort(-qc_intensity)
avgs = [-2, 0, 1]
counts = [0, 0, 0]
clusters = np.empty(len(order), dtype=int)
for i in order:
    logfc = organism_logfcs[i]
    if not (-10 < logfc < 10):
        continue
    if logfc < (avgs[0] + avgs[1]) / 2:
        clusters[i] = 2
        avgs[0] = (avgs[0] * counts[0] + logfc) / (counts[0] + 1)
        counts[0] += 1
    elif logfc < (avgs[1] + avgs[2]) / 2:
        clusters[i] = 1
        avgs[1] = (avgs[1] * counts[1] + logfc) / (counts[1] + 1)
        counts[1] += 1
    else:
        clusters[i] = 0
        avgs[2] = (avgs[2] * counts[2] + logfc) / (counts[2] + 1)
        counts[2] += 1







mzd = anchors["MZ"][a] - anchors["MZ"][b]
mzd_0 = mzd[(clusters[aa] == clusters[bb]) & (clusters[aa] == 0)]
mzd_2 = mzd[(clusters[aa] == clusters[bb]) & (clusters[aa] == 2)]

margin = 0.001
ppm_margin = 5


mzd_ecoli = np.sort(mzd_0[mzd_0 > 0])
mzd_yeast = np.sort(mzd_2[mzd_2 > 0])

mzd_all = np.unique(np.concatenate([mzd_ecoli, mzd_yeast]))

# ecoli_left = np.searchsorted(mzd_ecoli - margin, mzd_all, "left")
# ecoli_right = np.searchsorted(mzd_ecoli + margin, mzd_all, "right")
# yeast_left = np.searchsorted(mzd_yeast - margin, mzd_all, "left")
# yeast_right = np.searchsorted(mzd_yeast + margin, mzd_all, "right")

ecoli_left = np.searchsorted(mzd_ecoli * (1 - ppm_margin/1000000), mzd_all, "left")
ecoli_right = np.searchsorted(mzd_ecoli * (1 + ppm_margin/1000000), mzd_all, "right")
yeast_left = np.searchsorted(mzd_yeast * (1 - ppm_margin/1000000), mzd_all, "left")
yeast_right = np.searchsorted(mzd_yeast * (1 + ppm_margin/1000000), mzd_all, "right")

ecoli_count = ecoli_left - ecoli_right
yeast_count = yeast_left - yeast_right

ecoli_relative = ecoli_count/len(mzd_ecoli)
yeast_relative = yeast_count/len(mzd_yeast)

ratio = np.log2(ecoli_relative / yeast_relative)
selection = np.flatnonzero(~np.isinf(ratio))

import matplotlib as mpl
import matplotlib.cm as cm

norm = mpl.colors.Normalize(vmin=np.min(ratio[selection]), vmax=np.max(ratio[selection]))
m = cm.ScalarMappable(norm=norm, cmap="RdYlGn")

# plt.scatter(mzd_all[selection], ecoli_relative[selection], marker=".", c=m.to_rgba(ratio[selection]))
plt.plot(mzd_all, ecoli_relative)
plt.plot(mzd_all, -yeast_relative)
# plt.scatter(
#     mzd_all[c],
#     np.zeros(len(c)),
#     c=np.log2(ecoli_count[c] + yeast_count[c]),
#     cmap="RdYlGn"
# )
plt.show()

# min_count = 1 / 100000
# c = ecoli_relative > min_count
# c |= yeast_relative > min_count
min_count = 100
c = ecoli_count > min_count
c |= yeast_count > min_count
c = np.flatnonzero(c)
# c = c[np.argsort(np.log2(ecoli_count[c] + yeast_count[c]))]
# plt.show(
#     plt.scatter(
#         mzd_all[c],
#         np.log2(
#             (
#                 ecoli_count[c] / np.sum(ecoli_count[c])
#             ) / (
#                 yeast_count[c] / np.sum(yeast_count[c])
#             )
#         ),
#         c=np.log2(ecoli_count[c] + yeast_count[c]),
#         cmap="RdYlGn"
#     )
# )
c = c[
    np.argsort(
        np.log2(
            (
                ecoli_count[c] / np.sum(ecoli_count[c])
            ) / (
                yeast_count[c] / np.sum(yeast_count[c])
            )
        )
    )
]
plt.show(
    plt.scatter(
        mzd_all[c],
        np.log2(ecoli_count[c] + yeast_count[c]),
        c=np.log2(
            (
                ecoli_count[c] / np.sum(ecoli_count[c])
            ) / (
                yeast_count[c] / np.sum(yeast_count[c])
            )
        ),
        cmap="RdYlGn"
    )
)


aas = src.io.loadJSON("AMINO_ACID_FILE_NAME", parameters)
aa_keys, aa_vals = list(zip(*aas.items()))
aa_mzs = np.array(aa_vals)
aa_indices = np.searchsorted(mzd_all, aa_mzs)
ecoli_aa = ecoli_count[aa_indices]
yeast_aa = yeast_count[aa_indices]

plt.bar(aa_vals, ecoli_aa / np.sum(ecoli_aa))
plt.bar(aa_vals, -yeast_aa / np.sum(yeast_aa))
plt.xticks(aa_vals, aa_keys)
plt.show()

plt.scatter(aa_vals, np.log2((ecoli_aa / np.sum(ecoli_aa)) / (yeast_aa / np.sum(yeast_aa))))











x = color != 3
x &= ~np.isnan(organism_logfcs)
x &= ~np.isinf(organism_logfcs)
np.sum(color[x] != clusters[x]) / np.sum(color[x] == clusters[x])

ax = plt.subplots(1, 4)[1]
ax[3].scatter(np.log2(qc_intensity), organism_logfcs, c=clusters, marker=".", s=1)
for i, c in enumerate(np.unique(color)):
    organism = organism_map[c]
    if organism == "UNKNOWN":
        continue
    tmp = ax[i].scatter(
        np.log2(qc_intensity),
        organism_logfcs,
        marker=".",
        c="grey",
        label=organism,
        s=1
        # s=scores
    )
    elements = color == c
    local_scores = scores[elements]
    X = np.stack([np.log2(qc_intensity[elements]), organism_logfcs[elements]]).T
    selection = ~np.any(np.isinf(X), axis=1)
    selection &= ~np.any(np.isnan(X), axis=1)
    local_scores = local_scores[selection]
    X = X[selection]
    # z = lof().fit(X)
    # outliers = z.predict(X=X)
    lof = sklearn.neighbors.LocalOutlierFactor()
    outliers = lof.fit_predict(X=X)
    fdr_estimate = np.sum(outliers!=1) / np.sum(outliers==1)
    print(fdr_estimate)
    tmp = ax[i].scatter(
        X[:,0],
        X[:,1],
        marker=".",
        c="red",
        # c=outliers,
        label=organism,
        cmap="RdYlGn",
        s=1
        # s=local_scores
    )
    # ax[i].xlabel("Log(QC))")
    # ax[i].ylabel("LogFC(A/B)")

plt.show()









proteins, protein_counts = np.unique(annotation_array[:,-1], return_counts=True)
protein_logfcs = defaultdict(list)
protein_intensities = defaultdict(list)
scores = annotation_array[:,-4].astype(float)
protein_scores = defaultdict(list)
for i, protein in enumerate(annotation_array[:,-1]):
    protein_logfcs[protein].append(organism_logfcs[i])
    protein_intensities[protein].append(qc_intensity[i])
    protein_scores[protein].append(scores[i])

protein_datasize = {p: len(i) for p, i in protein_logfcs.items()}
protein_datasize = np.array([protein_datasize[p] for p in annotation_array[:,-1]])
protein_scores = {p: np.sum(i) for p, i in protein_scores.items()}
protein_scores = np.array([protein_scores[p] for p in annotation_array[:,-1]])
protein_logfcs = {p: np.median(i) for p, i in protein_logfcs.items()}
protein_logfcs = np.array([protein_logfcs[p] for p in annotation_array[:,-1]])
protein_intensities = {p: np.median(i) for p, i in protein_intensities.items()}
protein_intensities = np.array([protein_intensities[p] for p in annotation_array[:,-1]])
for c in np.unique(color):
    organism = organism_map[c]
    if organism == "UNKNOWN":
        continue
    elements = color == c
    tmp = plt.scatter(
        np.log2(protein_intensities[elements]),
        protein_logfcs[elements],
        marker=".",
        c=color_palette[c],
        s=protein_scores[elements],
        label=organism
    )
    # finite = np.isfinite(organism_logfcs[elements])
    # finite &= np.isfinite(np.log2(qc_intensity[elements]))
    # plt.show(
    #     sns.jointplot(
    #         np.log2(qc_intensity[elements][finite]),
    #         organism_logfcs[elements][finite],
    #         kind="kde"
    #     )
    # )
    # plt.show(plt.boxplot(organism_logfcs[elements][finite]))

plt.legend()
plt.xlabel("Log(QC))")
plt.ylabel("LogFC(A/B)")
plt.show()











































# neighbor_counts = np.diff(neighbors.indptr)
# neighbor_count_subsets = [
#     np.unique(
#         neighbor_counts[anchors["ION_COUNT"]==i],
#         return_counts=True
#     ) for i in range(2, parameters["SAMPLE_COUNT"] + 1)
# ]


neighbor_count_subsets = [
    np.unique(
        np.diff((neighbors == i).indptr),
        return_counts=True
    ) for i in range(2, parameters["SAMPLE_COUNT"] + 1)
]

for a, b in neighbor_count_subsets:
    tmp = plt.plot(a, b, marker=".")

plt.legend(
    [
        "Overlap of at least size {}".format(i) for i in range(2, parameters["SAMPLE_COUNT"] + 1)
    ]
)
plt.show()




































from scipy.optimize import curve_fit

def generalfunc(lam):
    def func(data, size, scale):
        return size * scipy.stats.tukeylambda.pdf(data, lam, scale=scale)
    return func



attribute = "CALIBRATED_MZ"
avgs = np.average(estimation_anchors[attribute], axis=1)
d = (estimation_anchors[attribute] - avgs.reshape(-1, 1))
d *= 1000000 / estimation_anchors[attribute]
l, r = np.percentile(d, (1, 99))
d = d[(d > l) & (d < r)]
# d = estimation_anchors[attribute][:, 0] - estimation_anchors[attribute][:, 1]
# x, y = np.unique(np.round(d.flatten(), 3), return_counts=True)
d = d.flatten()
x, y = np.unique(np.round(d, 1), return_counts=True)
max_lambda = scipy.stats.ppcc_max(d)
func = generalfunc(max_lambda)
popt, pcov = curve_fit(func, x, y)
xx = np.linspace(np.min(x), np.max(x), 1000)
yy = popt[0] * scipy.stats.tukeylambda.pdf(xx, max_lambda, scale=popt[1])
plt.plot(x, y, marker=".")
plt.plot(xx, yy, color="r")
plt.show()






from scipy.optimize import curve_fit
def func(data, size, nc, loc, scale):
    return size * scipy.stats.ncx2.pdf(data, 2, nc, loc=loc, scale=scale)


attribute = "CALIBRATED_RT"
ptps = np.ptp(estimation_anchors[attribute], axis=1)
x, y = np.unique(np.round(ptps, 3), return_counts=True)
popt, pcov = curve_fit(func, x, y)
xx = np.linspace(np.min(x), np.max(x), 1000)
yy = popt[0] * scipy.stats.ncx2.pdf(xx, popt[1], popt[2], loc=popt[3], scale=popt[4])
yy = 4000 * scipy.stats.ncx2.pdf(xx, 4, 2)
plt.plot(x, y, marker=".")
plt.plot(xx, yy, color="r")
plt.show()



































#
#
# std_errors = np.stack(
#     [
#         np.abs(predicted_mzs[i] - masses) / np.std(predicted_mzs_errors[i]) for i in np.arange(
#             1,
#             parameters["MAXIMUM_PRECURSOR_CHARGE"] + 1
#         )
#     ]
# )
# ds = [
#     np.percentile(std_errors[i], range(101)) for i in range(3)
# ]
#
#
# errors = np.stack(
#     [
#         (predicted_mzs[i] - masses) / predicted_mzs[i] for i in np.arange(
#             1,
#             parameters["MAXIMUM_PRECURSOR_CHARGE"] + 1
#         )
#     ]
# )
#
# plt.show(plt.plot(*np.unique(np.round(errors[1], 2), return_counts=True)))











































ai = anchor_ions.copy()
ai.data = ions["CALIBRATED_INTENSITY"][ai.data]
d = {}
for i in range(0, 1000, 50):
    l = np.sum(ai > i, axis=1).A.flatten()
    d[i] = np.unique(l, return_counts=True)

z = []
for i, j in d.items():
    k = np.zeros(9)
    k[j[0]] = j[1]
    z.append(k)

z=np.stack(z)

plt.show([plt.plot(i) for i in z])





n = neighbors[significant_anchors].T.tocsr()[significant_anchors]
a,b = n.nonzero()

p1 = significant_peptides[a]
p2 = significant_peptides[b]

s = np.abs(peptide_masses[p1]-peptide_masses[p2]) < 1
s &= np.abs(peptide_masses[p1]-peptide_masses[p2]) > 0

aa = anchors[significant_anchors[a[s]]]







n = neighbors[significant_anchors]
a, b = n.nonzero()

mzd = anchors["MZ"][b] - anchors["MZ"][significant_anchors[a]]


isotopes = np.abs(mzd - base_mass_dict["atoms"]["H+"]) < 0.01
isotopes |= np.abs(mzd - 2 * base_mass_dict["atoms"]["H+"]) < 0.01
isotopes |= np.abs(mzd - 3 * base_mass_dict["atoms"]["H+"]) < 0.01



for i in range(1, 100):
    l=significant_peptides==a[-i]
    ai=anchor_ions[significant_anchors[l]]
    ai.data = ions["CALIBRATED_INTENSITY"][ai.data]
    z=ai[np.diff(ai.indptr)==10]
    zz=z.todense().A
    plt.show(plt.plot(zz.T))



s = np.flatnonzero(anchors["ION_COUNT"] == parameters["SAMPLE_COUNT"])
# rts = anchors["RT"][s]
# order = np.argsort(rts)
# rts = rts[order]
# s = s[order]
a = anchor_ions[s].todense().A
# b = (ions["CALIBRATED_DT"][a] - anchors["DT"][s].reshape(-1, 1))
b = 1000000 * (ions["CALIBRATED_MZ"][a] - anchors["MZ"][s].reshape(-1, 1)) / anchors["MZ"][s].reshape(-1, 1)
b = ions["CALIBRATED_RT"][a] - anchors["RT"][s].reshape(-1, 1)
plt.show([[plt.plot(np.percentile(b[:,i],np.arange(101)),np.arange(101)),plt.axhline(50),plt.axvline(0)] for i in range(27)])
# plt.show(plt.plot(b))
# [plt.show(plt.plot(rts, b[:,i])) for i in range(parameters["SAMPLE_COUNT"])]



















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
selected_anchor_indices = percolated_anchors[significant_pims]
for anchor_index in selected_anchor_indices[:1]:
    # anchor_index = anchor_indices[selected_anchor_index]
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
    if (len(frequency) <= 3):
        continue
    counts = np.concatenate([[0], counts])
    frequency = np.concatenate(
        [
            [len(anchor_candidates)],
            np.cumsum(frequency[::-1])[::-1]
        ]
    )
    ransac = linear_model.RANSACRegressor()
    try:
        tmp = ransac.fit(
            counts.reshape(-1, 1)[:-1],
            np.log(frequency).reshape(-1, 1)[:-1]
        )
    except ValueError:
        continue
    score = -ransac.predict(counts[-1])[0][0]
    plt.bar(counts, frequency)
    plt.scatter(counts, frequency, marker=".")
    plt.xlabel("#Peptides with fragment count")
    plt.ylabel("Frequency")
    plt.show()
    plt.scatter(counts, np.log(frequency), marker=".")
    plt.plot(
        [
            0,
            counts[-1]
        ],
        [
            ransac.predict(counts[0])[0][0],
            ransac.predict(counts[-1])[0][0]
        ]
    )
    plt.xlabel("#Peptides with fragment count")
    plt.ylabel("Log frequency")
    plt.title("Score={}".format(score))
    plt.show()





























a, b = neighbors.nonzero()
c = (
    (ions["SAMPLE"][a] == 0) & (ions["SAMPLE"][b] == 1)
)
a_rt = ions["DT"][a[c]]
b_rt = ions["DT"][b[c]]
diff_rt = a_rt - b_rt
order = np.argsort(a_rt)
X = np.stack([a_rt[order], diff_rt[order]]).T
from sklearn.neighbors import KDTree
kdt = KDTree(X, leaf_size=30, metric='euclidean')
# TODO Faster?: https://github.com/nmslib/hnswlib
max_check = 10
ind = kdt.query(X, k=max_check, return_distance=False)

path_length = np.zeros(ind.shape[0], dtype=np.int)
for i in np.flatnonzero(np.any(ind > 0, axis=1))[::-1]:
    path_length[i] = np.max(path_length[ind[i]] + 1)

idx = np.argmax(path_length)
largest_path = [idx]
while(path_length[idx] > 1):
    s = ind[idx]
    idx = s[np.argmax(path_length[s] == path_length[idx] - 1)]
    largest_path.append(idx)


def moving_average(a, n=1000):
    result = np.cumsum(a, dtype=float)
    return (result[n:] - result[:-n]) / n


s = np.array(largest_path)
tmp = plt.scatter(X[:, 0][::10], X[:, 1][::10], marker=".", c="grey")
# tmp = plt.scatter(X[:, 0][s], X[:, 1][s], c="blue", marker=".")
# tmp = plt.plot(
#     moving_average(X[:, 0][s]),
#     moving_average(X[:, 1][s]),
#     c="red"
# )
plt.show()


#
# import sklearn.neighbors
# n_neighbors = 1000
# knn = sklearn.neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
# y_ = knn.fit(
#     X[10**5:2*10**5, 0].reshape(-1, 1),
#     X[10**5:2*10**5, 1].reshape(-1, 1)
# ).predict(
#     X[10**5:2*10**5, 0].reshape(-1, 1)
# )
#
# tmp = plt.plot(
#     moving_average(X[10**5:2*10**5, 0]),
#     moving_average(y_.reshape(1,-1)),
#     c="green"
# )









































# Single sample co-elution
sample_rt_indices = src.aggregates.__indexIonRT(
    anchor_ions,
    ions,
    parameters,
    log
)
sample_rts = [
    ions["RT"][sample_rt_indices[i]] for i in range(
        parameters["SAMPLE_COUNT"]
    )
]

max_rt_error = np.max(ions["RT_ERROR"])
sample_ions = [
    np.flatnonzero(ions["SAMPLE"] == i) for i in range(
        parameters["SAMPLE_COUNT"]
    )
]
all_rt_errors = [
    np.sqrt(
        ions[sample]["RT_ERROR"]**2 + max_rt_error**2
    ) for sample in sample_ions
]
low_indices = [
    np.searchsorted(
        sample_rts[sample],
        ions["RT"][sample_ions[sample]] - parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"] * rt_errors,
        "left"
    ) for sample, rt_errors in enumerate(all_rt_errors)
]
high_indices = [
    np.searchsorted(
        sample_rts[sample],
        ions["RT"][sample_ions[sample]] + parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"] * rt_errors,
        "left"
    ) for sample, rt_errors in enumerate(all_rt_errors)
]

counts = np.zeros((len(ions), parameters["SAMPLE_COUNT"]), dtype=np.int)
for sample in range(parameters["SAMPLE_COUNT"]):
    for ion_index, low_index, high_index in zip(
        sample_ions[sample],
        low_indices[sample],
        high_indices[sample]
    ):
        ion = ions[ion_index]
        ion_rt = ion["RT"]
        ion_dt = ion["DT"]
        ion_sample = ion["SAMPLE"]
        ion_neighbors = sample_rt_indices[sample][low_index: high_index]
        neighbor_dt_error = np.abs(
            ions["DT"][ion_neighbors] - ion_dt
        )
        allowed_neighbor_dt_error = np.sqrt(
            ion["DT_ERROR"]**2 + ions["DT_ERROR"][ion_neighbors]**2
        ) * parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"]
        neighbor_rt_error = np.abs(
            ions["RT"][ion_neighbors] - ion_rt
        )
        allowed_neighbor_rt_error = np.sqrt(
            ion["RT_ERROR"]**2 + ions["RT_ERROR"][ion_neighbors]**2
        ) * parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"]
        ion_neighbor_count = np.sum(
            (
                neighbor_rt_error < allowed_neighbor_rt_error
            ) & (
                neighbor_dt_error < allowed_neighbor_dt_error
            )
        )
        counts[ion_index, sample] = ion_neighbor_count


































file_name = "/media/proteomics/MISC/Sander/APEX/20130423_Tenzer_UDMSE_LFQ_apex2d/5/func002.csv"
data = pd.read_csv(file_name, sep=",").values


def xic_old(data, mz, dt, mz_ppm=10, dt_err=1):
    RT = 2
    MZ = 3
    DT = 5
    INT = 9
    l = data[:, MZ] > mz * (1 - mz_ppm / 1000000)
    l &= data[:, MZ] < mz * (1 + mz_ppm / 1000000)
    l &= data[:, DT] > dt - dt_err
    l &= data[:, DT] < dt + dt_err
    s = data[l]
    plt.show(plt.plot(s[:, RT], s[:, INT], marker="."))



def xic(data, mz, dt, mz_ppm=10, dt_err=1):
    RT = 2
    MZ = 3
    DT = 5
    INT = 9
    l = data[:, MZ] > mz * (1 - mz_ppm / 1000000)
    l &= data[:, MZ] < mz * (1 + mz_ppm / 1000000)
    l &= data[:, DT] > dt - dt_err
    l &= data[:, DT] < dt + dt_err
    s = data[l]
    rts = np.concatenate([[-1], s[:, RT], [-1]])
    ints = np.cumsum(np.concatenate([[0, 0], s[:, INT], [0]]))
    x = np.flatnonzero(np.diff(rts) != 0)
    y = np.diff(ints[x])
    x = rts[x[1:]]
    plt.show(plt.plot(x, y, marker="."))


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm
import matplotlib.colors

def rtdt(data, mz, dt, rt, mz_ppm=10, dt_err=1, rt_err=1):
    RT = 2
    MZ = 3
    DT = 5
    INT = 9
    l = data[:, MZ] > mz * (1 - mz_ppm / 1000000)
    l &= data[:, MZ] < mz * (1 + mz_ppm / 1000000)
    l &= data[:, DT] > dt - dt_err
    l &= data[:, DT] < dt + dt_err
    l &= data[:, RT] > rt - rt_err
    l &= data[:, RT] < rt + rt_err
    s = data[l]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(
            np.min(s[:, INT]),
            np.max(s[:, INT])
        ),
        cmap="RdYlGn"
    ).to_rgba(s[:, INT])
    plt.show(ax.scatter(s[:, RT], s[:, DT], s[:, MZ], c=colors))






























percolated_peptides = pd.read_csv(
    parameters["PERCOLATOR_TARGET_PEPTIDES"],
    delimiter="\t"
)
percolated_peptide_fdrs = percolated_peptides.values[:, 2]
peptide_fdr = 0.01
significant_percolated_peptides = np.array(
    [
        p.split(".")[1] for p in percolated_peptides.values[
            percolated_peptide_fdrs <= peptide_fdr,
            4
        ]
    ]
)

percolated_proteins = pd.read_csv(
    parameters["PERCOLATOR_TARGET_PROTEINS"],
    delimiter="\t"
)
percolated_protein_fdrs = percolated_proteins.values[:, 2]
protein_fdr = 0.01
significant_percolated_proteins = percolated_proteins.values[
    percolated_protein_fdrs <= protein_fdr,
    0
]

unique_organisms = ["ECOLI", "YEAST", "HUMAN"]
select = np.isin(significant_peptide_sequences, significant_percolated_peptides)
select &= np.isin(significant_proteins, significant_percolated_proteins)
select &= np.isin(significant_organisms, unique_organisms)
significant_anchors = significant_anchors[select]
significant_peptides = significant_peptides[select]
significant_proteins = significant_proteins[select]
significant_organisms = significant_organisms[select]
significant_peptide_sequences = significant_peptide_sequences[select]
n = neighbors[significant_anchors].T.tocsr()[significant_anchors]
cluster_count, significant_clusters = scipy.sparse.csgraph.connected_components(
    n,
    directed=False,
    return_labels=True
)
anchor_intensities = anchor_ions[significant_anchors]
anchor_intensities.data = ions["INTENSITY"][anchor_intensities.data]
anchor_intensities = anchor_intensities.todense().A
quant_file_name = ".tmp/tenzer_quant/tenzer_quant.csv"
with open(quant_file_name, "w") as outfile:
    csv_file = csv.writer(outfile, delimiter="\t")
    tmp = csv_file.writerow(
        [
            "index",
            "anchor",
            "cluster",
            "peptide",
            "protein",
        ] + ["Intensity " + x.split("/")[-1] for x in parameters["APEX_FILE_NAMES"]]
    )
    for index, anchor_index in enumerate(significant_anchors):
        cluster = significant_clusters[index]
        peptide = significant_peptide_sequences[index]
        protein = significant_proteins[index]
        # ion_start_index = anchor_intensities.indptr[index]
        # ion_end_index = anchor_intensities.indptr[index + 1]
        # samples = anchor_intensities.indices[ion_start_index: ion_end_index]
        # intensities = anchor_intensities.indices[ion_start_index: ion_end_index]
        intensities = anchor_intensities[index]
        tmp = csv_file.writerow(
            [
                index,
                anchor_index,
                cluster,
                peptide,
                protein,
            ] + intensities.tolist()
        )

















def __indexIonRT(anchor_ions, ions, parameters, log):
    with log.newSection("Indexing ion RTs"):
        sample_ions = anchor_ions.T.tocsr()
        sample_rt_indices = []
        for sample in range(parameters["SAMPLE_COUNT"]):
            selected_ions = sample_ions.data[
                sample_ions.indptr[sample]: sample_ions.indptr[sample + 1]
            ]
            sample_rt_order = np.argsort(ions["RT"][selected_ions])
            sample_rt_indices.append(selected_ions[sample_rt_order])
        sample_rts = [
            ions["RT"][sample_rt_indices[i]] for i in range(
                parameters["SAMPLE_COUNT"]
            )
        ]
        log.printMessage("Estimating max rt errors")
        max_rt_error = np.max(ions["RT_ERROR"])
        sample_ions = [
            np.flatnonzero(ions["SAMPLE"] == i) for i in range(
                parameters["SAMPLE_COUNT"]
            )
        ]
        all_rt_errors = [
            np.sqrt(2) * ions[sample]["RT_ERROR"] for sample in sample_ions
        ]
        log.printMessage("Indexing lower rt borders")
        low_indices = [
            np.searchsorted(
                sample_rts[sample],
                ions["RT"][sample_ions[sample]] - parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"] * rt_errors,
                "left"
            ) for sample, rt_errors in enumerate(all_rt_errors)
        ]
        log.printMessage("Indexing upper rt borders")
        high_indices = [
            np.searchsorted(
                sample_rts[sample],
                ions["RT"][sample_ions[sample]] + parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"] * rt_errors,
                "left"
            ) for sample, rt_errors in enumerate(all_rt_errors)
        ]
        order = np.argsort(np.concatenate(sample_ions))
        low_indices = np.concatenate(low_indices)[order]
        high_indices = np.concatenate(high_indices)[order]
    return sample_rt_indices, low_indices, high_indices



sample_rt_indices, low_rt_indices, high_rt_indices = __indexIonRT(
    anchor_ions,
    ions,
    parameters,
    log
)

low_indices = low_rt_indices
high_indices = high_rt_indices

ion_index = 999999
anchor_index = ions[ion_index]["AGGREGATE_INDEX"]

ion = ions[ion_index]
ion_rt = ion["RT"]
ion_dt = ion["DT"]
ion_sample = ion["SAMPLE"]
low_index = low_indices[ion_index]
high_index = high_indices[ion_index]
ion_neighbors = sample_rt_indices[ion_sample][low_index: high_index]
ion_neighbors = ion_neighbors[
    np.abs(
        ions["DT"][ion_neighbors] - ion_dt
    ) < np.sqrt(
        ion["DT_ERROR"]**2 + ions["DT_ERROR"][ion_neighbors]**2
    ) * parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"]
]
ion_neighbors = ion_neighbors[
    ions["AGGREGATE_INDEX"][ion_neighbors] < anchor_index
]
ion_neighbors = ion_neighbors[
    np.abs(
        ions["RT"][ion_neighbors] - ion_rt
    ) < np.sqrt(
        ion["RT_ERROR"]**2 + ions["RT_ERROR"][ion_neighbors]**2
    ) * parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"]
]

# s0_ions = anchor_ions.T.tocsr()[0].data
# a0 = ions["AGGREGATE_INDEX"][s0_ions]
# n = neighbors[a0].T.tocsr()[a0]
# anchor_count, ion_labels = scipy.sparse.csgraph.connected_components(
#     n,
#     directed=False,
#     return_labels=True
# )
