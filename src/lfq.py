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
from collections import defaultdict
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn import linear_model
import csv
# import importlib
# importlib.reload(src.aggregates)
# parameter_file_name = "data/lfq_single_A_B_QC/parameters_res_auto.json"
parameter_file_name = "data/test/parameters.json"
# parameter_file_name = "data/k562_3_samples/parameters_res_25000.json"
# parameter_file_name = "data/comparison_k562_lfq/parameters_k562.json"
# parameter_file_name = "data/comparison_k562_lfq/parameters_lfq.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
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























conditions = {
    "A": slice(None, 9),
    "B": slice(9, 18),
    "QC": slice(18, None)
}
# conditions = {
#     "A": 0,
#     "B": 1,
#     "QC": 2
# }
# conditions = {
#     "A": slice(0, 10, 2),
#     "B": slice(1, 10, 2),
# }
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




# plt.show(
#     plt.boxplot(
#         [
#             anchors_per_condition["QC"]["CV"][
#                 anchors_per_condition["QC"]["ION_COUNT"] == i
#             ] for i in range(10)
#         ]
#     )
# )

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


logfcs = np.log(anchors_per_condition["A"]["INTENSITY"] / anchors_per_condition["B"]["INTENSITY"])
# good_indices = anchors_per_condition["A"]["ION_COUNT"] > 0
# good_indices &= anchors_per_condition["B"]["ION_COUNT"] > 0

organism_classification = [
    np.unique(
        np.round(
            logfcs[
                (
                    anchors_per_condition["A"]["ION_COUNT"] > i
                ) & (
                    anchors_per_condition["B"]["ION_COUNT"] > i
                ) #& anchors["LE"]
            ],
            1
        ),
        return_counts=True
    ) for i in range(9 + 1)
]
plt.show(
    [
        plt.plot(
            a, b, marker="."
        ) for (a, b) in organism_classification
    ] + [
        plt.legend(
            [
                "Found in at least {} samples of both condition A and B".format(
                    i
                ) for i in range(1, 9 + 1)
            ]
        ),
        plt.xlabel("LogFC(A/B)"),
        plt.ylabel("Frequency"),
    ]
)


classes = np.round(logfcs)

np.unique(classes[~np.isnan(classes)], return_counts=True)

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

groups = np.array(
    [
        np.flatnonzero((overlap == i)) for i in range(
            2, parameters["SAMPLE_COUNT"] + 1
        )
    ]
)

# groups = np.array(
#     [
#         np.flatnonzero(
#             (overlap == 10) & (np.round(
#                 anchors["CALIBRATED_LOGINT"][a_indices]
#             ) == i) & ~anchors["LE"][a_indices]
#         ) for i in range(
#             0, 25
#         )
#     ]
# )

group_sizes = np.array([len(i) for i in groups])

equal_classes = np.array(
    [
        np.sum(a_classes[i] == b_classes[i]) for i in groups
    ]
)

yeast_class = np.array(
    [
        (np.sum(a_classes[i] == 1) + np.sum(b_classes[i] == 1)) / 2 for i in groups
    ]
)

ecoli_class = np.array(
    [
        (np.sum(a_classes[i] == -2) + np.sum(b_classes[i] == -2)) / 2 for i in groups
    ]
)

human_class = np.array(
    [
        (np.sum(a_classes[i] == 0) + np.sum(b_classes[i] == 0)) / 2 for i in groups
    ]
)

expected_classes = np.array(
    [
        (a**2 + b**2 + c**2) / t**2 for a, b, c, t in zip(
            yeast_class,
            ecoli_class,
            human_class,
            group_sizes
        )
    ]
)

groupsize_plot = plt.plot(group_sizes / np.max(group_sizes))
yeast_plot = plt.plot(yeast_class / group_sizes)
ecoli_plot = plt.plot(ecoli_class / group_sizes)
human_plot = plt.plot(human_class / group_sizes)
equal_plot = plt.plot(equal_classes / group_sizes)
expected_plot = plt.plot(expected_classes)
plt.legend(
    ["group_sizes", "yeast", "ecoli", "human", "equal", "expected"],
)
# plt.xticks(range(25), range(0,25))
plt.xticks(range(parameters["SAMPLE_COUNT"] - 1), range(2, parameters["SAMPLE_COUNT"] + 1))
plt.show()


















neighbors = scipy.sparse.load_npz(parameters["ION_NEIGHBORS_FILE_NAME"])
neighbors += neighbors.T

def reStitchAnchors():
    a, b = neighbors.nonzero()
    s = np.stack(
        [
            ions["ANCHOR"][a],
            ions["ANCHOR"][b]
        ]
    )
    n = scipy.sparse.csr_matrix(
        (
            s[0] != s[1],
            s
        ),
        shape=(
            len(anchors),
            len(anchors)
        )
    )
    n.eliminate_zeros()
    isotope_count, anchor_labels = scipy.sparse.csgraph.connected_components(
        n,
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
    anchor_groups = []
    for i in np.flatnonzero(np.diff(anchor_label_breaks) > 1):
        isotope = anchor_order[
            anchor_label_breaks[i]: anchor_label_breaks[i + 1]
        ]
        anchor_groups.append(isotope)
    return np.array(anchor_groups)
    # defined = len(np.concatenate(anchor_groups))
    # undefined = len(anchors) - defined
    # for isotope_index, isotope_anchors in enumerate(anchor_groups):
    #     # TODO +1 for isotopes only?
    #     anchors["SPECTRUM"][isotope_anchors] = isotope_index + undefined
    # anchors["SPECTRUM"][anchors["SPECTRUM"] == 0] = np.arange(undefined)


ambiguous = np.array(
    [
        len(s) != len(np.unique(s)) for s in [
            np.concatenate(
                [
                    anchor_ions.indices[
                        anchor_ions.indptr[anchor]: anchor_ions.indptr[anchor + 1]
                    ] for anchor in anchor_group
                ]
            ) for anchor_group in anchor_groups
        ]
    ]
)










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

















# for i in range(spectrum_anchors.shape[0]):
#     a = spectrum_anchors[i].indices
#     if len(a) < 5:
#         continue
#     for anchor in anchors[a]:
#         tmp = plt.plot(
#             [anchor["CALIBRATED_MZ"]]*2,
#             [0, 2**anchor["CALIBRATED_LOGINT"]],
#             c="black"
#         )
#     plt.get_current_fig_manager().resize(width=1920, height=1080)
#     plt.show()
    # plt.show(plt.line(anchors["CALIBRATED_MZ"][a], 2**anchors["CALIBRATED_LOGINT"][a]))
    # input(np.sort(anchors[a]))


aas = {
    "G": 57.021464,
    "A": 71.037114,
    "S": 87.032028,
    "P": 97.052764,
    "V": 99.068414,
    "T": 101.047679,
    "C": 103.009185 + 57.021464, # + carbamido = 160.030649
    "I": 113.084064,
    "L": 113.084064,
    "N": 114.042927,
    "D": 115.026943,
    "Q": 128.058578,
    "K": 128.094963,
    "E": 129.042593,
    "M": 131.040485,
    "H": 137.058912,
    "F": 147.068414,
    "R": 156.101111,
    "Y": 163.063329,
    "W": 186.079313,
}


for s in spectra:
    if len(s) < 2:
        continue
    a = am[s]
    m = anchor_mass[a]
    diffs = m[:,np.newaxis] - m
    sums = m[:,np.newaxis] + m
    sums, sum_counts = np.unique(np.round(sums.flatten(), 2), return_counts=True)
    max_ind = np.argmax(sum_counts)
    np.diff(np.sort(m))
    np.sort(m)
    if sum_counts[max_ind] >= 4:
        sums[max_ind], sum_counts[max_ind], np.any(np.abs(m - sums[max_ind]) < 0.1)
    for anchor, mass in zip(anchors[a], m):
        tmp = plt.plot(
            [mass]*2,
            [0, 2**anchor["CALIBRATED_LOGINT"]],
            c="black"
        )
    plt.get_current_fig_manager().resize(width=1920, height=1080)
    plt.show()









































































































































# logpid() { while sleep 1; do ps -p $1 -o pcpu= -o pmem= ; done; }
# src/main.py -p data/test/parameters.json & logpid $! > log_cpu.txt


# import psutil
# pid = 4384
# psutil.pid_exists(pid)
# p = psutil.Process(pid)
# sum(i.cpu_percent(interval=0.1) for i in p.children(recursive=True))+p.cpu_percent(interval=0.1)
# sum(i.memory_percent(memtype="pss") for i in p.children(recursive=True))+p.memory_percent(memtype="pss")


#
#
#
#
#
#
# raise Exception("Not main function!")
#
#

# import importlib
# importlib.reload(src.aggregates)
import src.parameters
import src.io
import src.ions
import src.aggregates
import src.peptides
import pandas as pd
import numpy as np
import src.parellelization as mp
import scipy.sparse
from collections import defaultdict
import json
import sys
import time
import traceback as tb
from scipy.special import binom as binom
from contextlib import contextmanager
from collections import defaultdict
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn import linear_model
import csv
# parameter_file_name = "data/ypic/prop/parameters.json"
# parameter_file_name = "data/ypic/normal/parameters.json"
# parameter_file_name = "data/k562/parameters.json"
parameter_file_name = "data/test/parameters.json"
# parameter_file_name = "data/lfq/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"])
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
# base_mass_dict = src.peptides.loadBaseMassDict(parameters, log)
# proteins, total_protein_sequence, ptms, ptm_matrix = src.peptides.importProteinsAndPtms(parameters, log)
# peptides, peptide_index_matrix, digestion_matrix = src.peptides.digestProteins(
#     proteins,
#     total_protein_sequence,
#     ptm_matrix,
#     parameters,
#     log,
#     digestion_aas="KR"
# )
# peptide_masses, fragments = src.peptides.calculateMasses(
#     peptides,
#     peptide_index_matrix,
#     base_mass_dict,
#     total_protein_sequence,
#     parameters,
#     log,
# )
# parameters["FILTER_PRECURSOR_EXISTENCE_IN_LE"] = True
# parameters["PERCOLATOR_DATA_FILE_NAME"] = 'data/k562/results/percolator_le.csv'
# anchor_boundaries, fragment_peptide_indices, fragment_indices = src.aggregates.matchAnchorsToFragments(
#     fragments,
#     anchors,
#     base_mass_dict,
#     parameters,
#     log
# )
# anchor_peptide_scores, anchor_peptide_match_counts = src.aggregates.getAnchorPeptideMatrix(
#     anchors,
#     neighbors,
#     peptides,
#     anchor_boundaries,
#     fragment_peptide_indices,
#     parameters,
#     log
# )
# anchor_fragment_indices = src.aggregates.getAnchorFragmentIndices(
#     anchor_peptide_match_counts,
#     anchor_boundaries,
#     fragment_indices,
#     fragment_peptide_indices,
#     parameters,
#     log
# )
# # a=np.unique(anchor_peptide_scores.indices[(anchor_peptide_scores.data>=4.5)])
# # np.unique(peptides["DECOY"][a],return_counts=True)[1]
# # parameters["FILTER_PRECURSOR_EXISTENCE_IN_LE"] = True
# peptide_boundaries, anchor_indices = src.aggregates.matchPeptidesToAnchors(
#     anchors,
#     peptide_masses,
#     base_mass_dict,
#     parameters,
#     log,
# )
# precursor_indices = src.aggregates.findFragmentPrecursors(
#     anchor_peptide_match_counts,
#     anchors,
#     anchor_indices,
#     peptide_boundaries,
#     neighbors,
#     anchor_alignment_parameters,
#     parameters,
#     log
# )
# import importlib
# importlib.reload(src.io)
# importlib.reload(src.aggregates)
# src.aggregates.writePercolatorFile(
#     anchors,
#     base_mass_dict,
#     anchor_peptide_match_counts,
#     fragments,
#     anchor_fragment_indices,
#     neighbors,
#     peptides,
#     peptide_masses,
#     precursor_indices,
#     anchor_peptide_scores,
#     peptide_index_matrix,
#     total_protein_sequence,
#     parameters,
#     log
# )
#
#
#
#


























#
#
#
# fragment_mzs = np.concatenate(
#     [
#         fragments["Y_MR"] + base_mass_dict["atoms"]["H+"],
#         fragments["B_MR"] + base_mass_dict["atoms"]["H+"],
#     ]
# )
# fragment_peptide_indices = np.concatenate(
#     [fragments["PEPTIDE"]] * 2
# )
# order = np.flatnonzero(fragment_mzs > base_mass_dict["atoms"]["H+"])
# order = order[np.argsort(fragment_mzs[order])]
# fragment_mzs = fragment_mzs[order]
# fragment_peptide_indices = fragment_peptide_indices[order]
# # Id
# ppm = 20
# process_count = parameters["CPU_COUNT"]
# order = np.flatnonzero((np.diff(neighbors.indptr) >= 2) & (~anchors["LE"]))
# order = order[np.argsort(anchors["MZ"][order])]
# anchor_boundaries = np.zeros((len(anchors), 2), dtype=np.int)
#
#
# def getAnchorBounds(target_mzs, query_mzs, ppm):
#     hit_lower_boundaries = np.searchsorted(
#         target_mzs,
#         query_mzs * (1 - ppm / 1000000),
#         "left"
#     )
#     hit_upper_boundaries = np.searchsorted(
#         target_mzs,
#         query_mzs * (1 + ppm / 1000000),
#         "right"
#     )
#     # d = hit_upper_boundaries - hit_lower_boundaries
#     anchor_boundaries = np.stack(
#         [
#             hit_lower_boundaries,
#             hit_upper_boundaries
#         ]
#     ).T
#     return anchor_boundaries
#
#
# anchor_boundaries[order] = getAnchorBounds(
#     fragment_mzs,
#     anchors["MZ"][order],
#     ppm
# )
# anchor_indices = np.flatnonzero(
#     anchor_boundaries[:, 1] - anchor_boundaries[:, 0] >= 3 # TODO 3 are needed: 1 target and 2 for regression
# )
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# def multiprocessedAnnotatePerFragment(kwargs):
#     in_queue = kwargs['in_queue']
#     out_queue = kwargs['out_queue']
#     neighbors = kwargs['neighbors']
#     anchor_indices = kwargs['anchor_indices']
#     anchor_boundaries = kwargs['anchor_boundaries']
#     fragment_peptide_indices = kwargs['fragment_peptide_indices']
#     anchor_len = kwargs['anchor_len']
#     peptide_len = kwargs['peptide_len']
#     selected_anchor_indices = in_queue.get()
#     anchor_peptide_score = scipy.sparse.dok_matrix(
#         (anchor_len, peptide_len),
#         dtype=np.float
#     )
#     anchor_peptide_match_counts = scipy.sparse.dok_matrix(
#         (anchor_len, peptide_len),
#         dtype=np.int
#     )
#     # anchor_evidence = scipy.sparse.dok_matrix(
#     #     neighbors.shape,
#     #     dtype=np.int
#     # )
#     for selected_anchor_index in selected_anchor_indices:
#         anchor_index = anchor_indices[selected_anchor_index]
#         anchor_candidates = fragment_peptide_indices[
#             slice(*anchor_boundaries[anchor_index])
#         ]
#         anchor_neighbors = neighbors.indices[
#             neighbors.indptr[anchor_index]: neighbors.indptr[anchor_index + 1]
#         ]
#         raw_neighbor_candidates = np.concatenate(
#             [
#                 fragment_peptide_indices[
#                     slice(*anchor_boundaries[n])
#                 ] for n in anchor_neighbors
#             ]
#         )
#         neighbor_candidates = raw_neighbor_candidates[
#             np.isin(
#                 raw_neighbor_candidates,
#                 anchor_candidates
#             )
#         ]
#         candidates, candidate_counts = np.unique(
#             neighbor_candidates,
#             return_counts=True
#         )
#         counts, frequency = np.unique(
#             candidate_counts,
#             return_counts=True
#         )
#         if (len(frequency) < 2):# or frequency[-1] != 1: # TODO delete
#             continue
#         counts = np.concatenate([[0], counts])
#         frequency = np.concatenate(
#             [
#                 [len(anchor_candidates)],
#                 np.cumsum(frequency[::-1])[::-1]
#             ]
#         )
#         ransac = linear_model.RANSACRegressor()
#         try:
#             tmp = ransac.fit(
#                 counts.reshape(-1, 1)[:-1],
#                 np.log(frequency).reshape(-1, 1)[:-1]
#             )
#         except ValueError:
#             continue
#         score = -ransac.predict(counts[-1])[0][0]
#         # if score > 0:
#         if True:
#         #     plt.plot(counts, log_count_frequency, marker="o")
#         #     plt.plot(counts, [ransac.predict(i)[0][0] for i in counts], marker="o")
#         #     plt.show()
#             # peptide_index = candidates[candidate_counts == counts[-1]][0]
#             # anchor_peptide_score[anchor_index, peptide_index] = score
#             peptide_indices = candidates[candidate_counts == counts[-1]]
#             anchor_peptide_score[anchor_index, peptide_indices] = score
#             anchor_peptide_match_counts[anchor_index, peptide_indices] = counts[-1]
#             # anchor_evidence[
#             #     anchor_index,
#             #     anchor_neighbors
#             # ] = 1
#             # neighbor_indices = np.repeat(
#             #     anchor_neighbors,
#             #     np.diff(anchor_boundaries[anchor_neighbors], axis=1).flatten()
#             # )
#             # anchor_evidence[
#             #     anchor_index,
#             #     # neighbor_indices[raw_neighbor_candidates == peptide_index],
#             #     neighbor_indices[np.isin(raw_neighbor_candidates, peptide_indices)],
#             # ] = 2
#             # anchor_evidence[
#             #     anchor_index,
#             #     anchor_index
#             # ] = 3
#     out_queue.put(
#         # (anchor_peptide_score.tocsr(), anchor_evidence.tocsr())
#         (anchor_peptide_score.tocsr(), anchor_peptide_match_counts.tocsr())
#     )
#     out_queue.put(None)
#
#
#
#
# from sklearn import linear_model
# in_queue = mp.partitionedQueue(anchor_indices, process_count)
# anchor_peptide_scores = scipy.sparse.csr_matrix(
#     (len(anchors), len(peptides)),
#     dtype=np.float
# )
# anchor_peptide_match_counts = scipy.sparse.csr_matrix(
#     (len(anchors), len(peptides)),
#     dtype=np.int
# )
# for partial_anchor_peptide_scores, partial_anchor_match_counts in mp.parallelizedGenerator(
#     function=multiprocessedAnnotatePerFragment,
#     function_args={
#         'in_queue': in_queue,
#         'neighbors': neighbors,
#         'anchor_indices': anchor_indices,
#         'anchor_boundaries': anchor_boundaries,
#         'fragment_peptide_indices': fragment_peptide_indices,
#         'anchor_len': len(anchors),
#         'peptide_len': len(peptides)
#     },
#     process_count=process_count,
# ):
#     anchor_peptide_scores += partial_anchor_peptide_scores
#     anchor_peptide_match_counts += partial_anchor_match_counts
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# peptide_mrs = peptide_masses
# peptide_order = np.argsort(peptide_mrs)
#
#
# ppm = 10
# max_precursor_charge = 3
# process_count = parameters["CPU_COUNT"]
# order = np.flatnonzero((np.diff(neighbors.indptr) >= 2) & anchors["LE"])
# anchor_mzs = np.concatenate(
#     [
#         (anchors["MZ"][order] - base_mass_dict["atoms"]["H+"]) * i for i in range(
#             1,
#             max_precursor_charge + 1
#         )
#     ]
# )
# anchor_order = np.concatenate(
#     [order] * max_precursor_charge
# )
# order = np.argsort(anchor_mzs)
# anchor_mzs = anchor_mzs[order]
# # anchor_zs = anchor_zs[order]
# anchor_order = anchor_order[order]
#
#
# peptide_boundaries = np.zeros((len(peptides), 2), dtype=np.int)
# peptide_boundaries[peptide_order] = getAnchorBounds(
#     anchor_mzs,
#     peptide_mrs[peptide_order],
#     ppm
# )
#
#
# order = order[np.argsort(anchors["MZ"][order])]
# anchor_boundaries = np.zeros((len(anchors), 2), dtype=np.int)
# anchor_boundaries[order] = getAnchorBounds(
#     fragment_mzs,
#     anchors["MZ"][order],
#     ppm
# )
#
#
# anchor_peptide_scores2 = anchor_peptide_match_counts.astype(np.bool)
# anchor_peptide_scores2 += anchor_peptide_scores
# anchor_peptide_scores2.data -= 1
#
# match_scores = anchor_peptide_match_counts.data
# scores = anchor_peptide_scores2.data
# selected_anchor_indices, selected_peptide_indices = anchor_peptide_match_counts.nonzero()
# max_rt_difference = 0.05
# dt_shifts = []
# precursor_indices = []
# # precursor_charges = []
# for (anchor_index, peptide_index) in zip(selected_anchor_indices, selected_peptide_indices):
#     anchor = anchors[anchor_index]
#     candidate_precursor_anchors = anchor_order[
#         slice(*peptide_boundaries[peptide_index])
#     ]
#     # candidate_precursor_charges = anchor_order[
#     #     slice(*peptide_boundaries[peptide_index])
#     # ]
#     precursor_rt_difference = anchors["RT"][candidate_precursor_anchors] - anchor["RT"]
#     good_candidates = np.abs(precursor_rt_difference) < max_rt_difference # TODO, consistent in all runs
#     filtered_precursors = candidate_precursor_anchors[good_candidates]
#     # filtered_charges = candidate_precursor_charges[good_candidates]
#     if len(filtered_precursors) == 0:
#         dt_shifts.append(999)
#         precursor_indices.append(-1)
#         # precursor_charges.append(-1)
#         continue
#     dt_difference = anchors["DT"][filtered_precursors] - anchor["DT"]
#     min_dt_index = np.argmin(np.abs(dt_difference))
#     dt_shifts.append(
#         anchors["DT"][filtered_precursors[min_dt_index]] - anchor["DT"]
#     )
#     precursor_indices.append(filtered_precursors[min_dt_index])
#     # precursor_charges.append(filtered_charges[min_dt_index])
#
# dt_shifts = np.array(dt_shifts)
# precursor_indices = np.array(precursor_indices)
# # precursor_charges = np.array(precursor_charges)
#
# selected_anchor_indices, selected_peptide_indices = anchor_peptide_match_counts.nonzero()
# dt_filter = dt_shifts != 999
# new_anchor_indices = selected_anchor_indices[dt_filter]
# new_peptide_indices = selected_peptide_indices[dt_filter]
# new_scores = scores[dt_filter]
# new_match_scores = match_scores[dt_filter]
# dt_shifts = dt_shifts[dt_filter]
# precursor_indices = precursor_indices[dt_filter]
#
#
# # new_unique_anchor_indices, counts = np.unique(new_anchor_indices,return_counts=True)
# # unique_anchor_filter = np.isin(
# #     new_anchor_indices,
# #     new_unique_anchor_indices[counts==1]
# # )
# unique_anchor_filter = np.ones_like(new_anchor_indices, dtype=np.bool)
# selected_anchor_indices = new_anchor_indices[unique_anchor_filter]
# selected_peptide_indices = new_peptide_indices[unique_anchor_filter]
# selected_anchor_scores = new_scores[unique_anchor_filter]
# selected_anchor_match_scores = new_match_scores[unique_anchor_filter]
# dt_shifts = dt_shifts[unique_anchor_filter]
# precursor_indices = precursor_indices[unique_anchor_filter]
#
#
#
#
#
# peptide_fragment_borders = np.cumsum(
#     np.concatenate(
#         [
#             [0],
#             peptides["SIZE"]
#         ]
#     )
# )
# peptide_fragment_borders = np.stack(
#     [
#         peptide_fragment_borders[:-1],
#         peptide_fragment_borders[1:]
#     ]
# ).T
#
#
# # selected_anchor_indices, selected_peptide_indices = anchor_peptide_scores.nonzero()
# # selected_anchor_scores = anchor_peptide_scores.data
# # selected_raw_evidence_counts = np.diff((anchor_evidence[selected_anchor_indices] >= 1).indptr)
# selected_raw_evidence_counts = np.diff(neighbors[selected_anchor_indices].indptr)
# # selected_evidence_counts = np.diff((anchor_evidence[selected_anchor_indices] >= 2).indptr)
# selected_evidence_counts = selected_anchor_match_scores
#
# data = []
# import csv
# with open(parameters["OUTPUT_PATH"] + "/181010_ion_spectra.csv", "w") as raw_outfile:
#     outfile = csv.writer(raw_outfile, delimiter="\t")
#     header = [
#         "SpecId",
#         "Label",
#         "ScanNr",
#         "rt",
#         "dm",
#         "matches",
#         "score",
#         "dt",
#         "dt_drift",
#         # "le",
#         "anchor_count",
#         "match_ratio",
#         "anchor_mz",
#         "ion_type",
#         "precursor_mz",
#         "precursor_z",
#         "precursor_ppm",
#         # "logFC",
#         "Peptide",
#         "Proteins",
#     ]
#     tmp = outfile.writerow(header)
#     for index, anchor_index in enumerate(selected_anchor_indices):
#         if not ((1.5 < dt_shifts[index]) & (dt_shifts[index] < 4)):
#             continue
#         anchor = anchors[anchor_index]
#         peptide_index = selected_peptide_indices[index]
#         peptide = peptides[peptide_index]
#         precursor_index = precursor_indices[index]
#         precursor = anchors[precursor_index]
#         precursor_charge = int(round(peptide_masses[peptide_index] / precursor["MZ"]))
#         precursor_mr = (
#             precursor["MZ"] - base_mass_dict["atoms"]["H+"]
#         ) * precursor_charge
#         precursor_delta_mass = precursor_mr - peptide_masses[peptide_index]
#         precursor_ppm = precursor_delta_mass / precursor_mr * 1000000
#         selected_fragments = peptide_fragment_borders[peptide_index]
#         selected_fragments = fragments[
#             selected_fragments[0]: selected_fragments[1]
#         ]
#         dm_b = anchor["MZ"] - (selected_fragments["B_MR"] + base_mass_dict["atoms"]["H+"])
#         dm_y = anchor["MZ"] - (selected_fragments["Y_MR"] + base_mass_dict["atoms"]["H+"])
#         ion_type = np.concatenate(
#             [
#                 selected_fragments["B_INDEX"],
#                 -selected_fragments["Y_INDEX"]
#             ]
#         )
#         dm = np.concatenate([dm_b, dm_y])
#         dm /= anchor["MZ"]
#         dm *= 1000000
#         dm_min = np.argmin(np.abs(dm))
#         dm = dm[dm_min]
#         ion_type = ion_type[dm_min]
#         # protein_string = ";".join(
#         #     proteins["ID"][
#         #         peptide_protein_matrix.indices[
#         #             peptide_protein_matrix.indptr[peptide_index]: peptide_protein_matrix.indptr[peptide_index + 1]
#         #         ]
#         #     ]
#         # )
#         protein_string = "TODO"
#         peptide_start_index = peptide_index_matrix.indices[
#             peptide_index_matrix.indptr[peptide_index]
#         ]
#         peptide_sequence = total_protein_sequence[
#             peptide_start_index: peptide_start_index + peptide["SIZE"]
#         ]
#         row = [
#             anchor_index,
#             -1 if peptide["DECOY"] else 1,
#             anchor_index,
#             anchor["RT"],
#             dm,
#             selected_evidence_counts[index],
#             selected_anchor_scores[index],
#             anchor["DT"],
#             dt_shifts[index],
#             # 0 if anchor["LE"] else 1,
#             selected_raw_evidence_counts[index],
#             selected_evidence_counts[index] / selected_raw_evidence_counts[index],
#             anchor["MZ"],
#             ion_type,
#             precursor["MZ"],
#             precursor_charge,
#             precursor_ppm,
#             # diffs[anchor_index],
#             "-.{}.-".format(peptide_sequence),
#             # "PROTEIN",
#             protein_string,
#         ]
#         tmp = outfile.writerow(row)
#         data.append(row)
#
#
#
#
#
#
#
#
#
#
#
#









































































































#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # i=29999;anchors[i];ions[anchor_ions[i].data]
#
# selected_anchor_ions = [
#     np.concatenate(
#         [
#             anchor_ions.data[
#                 anchor_ions.indptr[j]: anchor_ions.indptr[j + 1]
#             ] for j in np.flatnonzero(anchors["ION_COUNT"] == i)
#         ]
#     ).reshape(-1, i) for i in range(3, parameters["SAMPLE_COUNT"] + 1)
# ]
# cvs_normal = [
#     scipy.stats.variation(ions["INTENSITY"][m], axis=1) for m in selected_anchor_ions
# ]
# cvs_calibrated = [
#     scipy.stats.variation(ions["CALIBRATED_INTENSITY"][m], axis=1) for m in selected_anchor_ions
# ]
# plt.show(plt.boxplot(cvs_normal))
# plt.show(plt.boxplot(cvs_calibrated))
# [
#     plt.show(
#         [
#             plt.plot(
#                 np.percentile(
#                     a,
#                     range(101)
#                 )
#             ),
#             plt.plot(
#                 np.percentile(
#                     b,
#                     range(101)
#                 )
#             )
#         ]
#     ) for a, b in zip(cvs_normal[::-1], cvs_calibrated[::-1])
# ]
#
# avg_calibrated_logints = [
#     np.log(np.average(ions["CALIBRATED_INTENSITY"][m], axis=1)) for m in selected_anchor_ions
# ]
# plt.show(plt.boxplot(avg_calibrated_logints))
#
# import seaborn as sns
# plt.show(sns.jointplot(avg_calibrated_logints[-1][::100], cvs_calibrated[-1][::100], kind="kde"))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
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

































full_anchors = np.flatnonzero(anchors["ION_COUNT"] == parameters["SAMPLE_COUNT"])

left_anchors, right_anchors = neighbors.nonzero()

b = np.concatenate([np.array([0]), np.flatnonzero(np.diff(left_anchors) > 0)])

good_pairs = np.isin(left_anchors, full_anchors)
good_pairs &= np.isin(right_anchors, full_anchors)

left_anchors = left_anchors[good_pairs]
right_anchors = right_anchors[good_pairs]

left_ion_set = anchor_ions[left_anchors].toarray()
right_ion_set = anchor_ions[right_anchors].toarray()
log_matched_intensities = np.log(ions["INTENSITY"])[np.stack([left_ion_set, right_ion_set], axis=1)]

diffs = np.ptp(log_matched_intensities, axis=1)
avg_diff = np.average(diffs, axis=1)


centered_diffs = diffs - avg_diff.reshape(-1, 1)

max_diffs = np.max(np.abs(centered_diffs),axis=1)

relative_diffs = centered_diffs / max_diffs.reshape(-1, 1)











































selected = np.flatnonzero(
    (
        anchors["ION_COUNT"] == parameters["SAMPLE_COUNT"]
    )
)
full_anchors = anchors[selected]
full_anchor_ions = anchor_ions[selected].toarray()
mz_error = ion_alignment_parameters["CALIBRATED_MZ"]
rt_error = anchor_alignment_parameters["RT"]
le_he_pairs = []
upper_index = 1
for lower_index, anchor in enumerate(full_anchors):
    candidate_upper_mass = anchor["MZ"] * (1 + mz_error / 1000000)
    try:
        while full_anchors["MZ"][upper_index] <= candidate_upper_mass:
            upper_index += 1
    except IndexError:
        upper_index = len(full_anchors)
    candidates = np.flatnonzero(
        full_anchors["LE"][lower_index: upper_index] != full_anchors["LE"][lower_index]
    )
    if len(candidates) == 0:
        continue
    aggregate_rts = ions["RT"][full_anchor_ions[lower_index]]
    matching_rts = ions["RT"][full_anchor_ions[lower_index + candidates]]
    matches = candidates[
        np.all(
            np.abs(matching_rts - aggregate_rts) < rt_error,
            axis=1
        )
    ]
    if len(matches) == 1:
        if full_anchors["LE"][lower_index]:
            le_anchor = lower_index
            he_anchor = lower_index + matches[0]
        else:
            le_anchor = lower_index + matches[0]
            he_anchor = lower_index
        le_he_pairs.append((le_anchor, he_anchor))

le_he_pairs = selected[np.array(le_he_pairs)]
anchor_index, multiplicity = np.unique(le_he_pairs, return_counts=True)
multiply_used_anchors = anchor_index[multiplicity > 1]
le_he_pairs = le_he_pairs[
    ~np.any(
        np.isin(le_he_pairs, multiply_used_anchors),
        axis=1
    )
]









#
# dt_diffs = np.diff(anchors["DT"][le_he_pairs], axis=1).ravel()
# # dt_diffs -= np.median(dt_diffs)
# mz_diffs = np.diff(anchors["MZ"][le_he_pairs], axis=1).ravel()
# dts = anchors["DT"][le_he_pairs[:, 0]]
# mzs = anchors["MZ"][le_he_pairs[:, 0]]
#
#
#
# # a, b = np.unique(np.round(mz_diffs / mzs * 1000000, 1), return_counts=True)
# a, b = np.unique(np.round(dt_diffs, 1), return_counts=True)
# plt.show(plt.plot(a, b, marker="."))
# percentiles = np.percentile(dt_diffs, range(101))
# relative_percentiles = np.percentile(dt_diffs / dts, range(101))
# plt.show(plt.plot(percentiles))
# plt.show(plt.plot(relative_percentiles))
#
# # x = dt_diffs > percentiles[10]
# # x &= dt_diffs < percentiles[90]
# x = dt_diffs / dts > relative_percentiles[10]
# x &= dt_diffs / dts < relative_percentiles[90]
# x &= dts < 190
# x &= dts > 50
# # plt.show(plt.scatter(mzs[x], dt_diffs[x] / dts[x], c=dt_diffs[x] / dts[x], cmap="RdYlGn", marker="."))
# # plt.show(plt.scatter(mzs[x], mzs[x] / dts[x], c=dt_diffs[x] / dts[x], cmap="RdYlGn", marker="."))
# plt.show(plt.scatter(dts[x], mzs[x], c=dt_diffs[x] / dts[x], cmap="RdYlGn", marker="."))
#
#
#
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
#
#
# xdata = np.stack([dts, mzs])
# ydata = 1000000 * dt_diffs / dts
# # plt.show(plt.scatter(*xdata[:, x], c=ydata[x], cmap="RdYlGn", marker="."))
#
# def func(data, size_param, angle_param, ratio_param, constant):
#     x = data[0]
#     y = data[1]
#     size = (np.sqrt(x**2 + y**2))
#     angle = np.arctan(y / x)
#     ratio = x / y
#     return size_param * size + angle_param * angle + ratio_param * ratio + constant
#     # return a * x + b * y + c * x * y + d
#
# popt, pcov = curve_fit(func, xdata[:, x][:,::2], ydata[x][::2])
# ypred = func(xdata, *popt)
# scipy.stats.pearsonr(ydata[x], ypred[x])
#
# # plt.show(plt.scatter(ydata[x], ypred[x], marker="."))
# # plt.show(sns.jointplot(ydata[x], ypred[x], kind="hex"))
# plt.show(plt.scatter(*xdata[:, x], c=ypred[x], cmap="RdYlGn", marker="."))
# relative_dt_errors = (ypred[x][1::2] - ydata[x][1::2])
# absolute_dt_errors = relative_dt_errors / 1000000 * dts[x][1::2]
# plt.show(plt.scatter(dts[x][1::2], relative_dt_errors, marker="."))
# plt.show(sns.jointplot(dts[x][1::2], relative_dt_errors, kind="hex"))
# plt.show(plt.boxplot(absolute_dt_errors))
# p = np.percentile(relative_dt_errors, range(101))
# second_derivative_p = np.gradient(np.gradient(p))
# plt.show(plt.plot(p, marker="."))
#
# new_dt = dts[x] + ypred[x] * dts[x]/1000000
# new_dt_diff = (new_dt - ions["DT"][z[:, 1]][x]) / new_dt * 1000000
# p = np.percentile(new_dt_diff, range(101))
# plt.show(plt.plot(p, marker="."))
#









import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def func(data, size_param, angle_param, ratio_param, constant):
    x = data[0]
    y = data[1]
    size = (np.sqrt(x**2 + y**2))
    angle = np.arctan(y / x)
    ratio = x / y
    return size_param * size + angle_param * angle + ratio_param * ratio + constant
    # return a * x + b * y + c * x * y + d


alldata = np.stack([ions["DT"], ions["CALIBRATED_MZ"]])
shift_parameters = []
new_pred = np.zeros(len(ions))
for sample in range(parameters["SAMPLE_COUNT"]):
    sample_ions = anchor_ions[le_he_pairs, sample].toarray()
    dt_diffs = np.diff(ions["DT"][sample_ions], axis=1).ravel()
    dts = ions["DT"][sample_ions[:, 0]]
    mzs = ions["CALIBRATED_MZ"][sample_ions[:, 0]]
    # # a, b = np.unique(np.round(mz_diffs / mzs * 1000000, 1), return_counts=True)
    # a, b = np.unique(np.round(dt_diffs, 1), return_counts=True)
    # plt.show(plt.plot(a, b, marker="."))
    # percentiles = np.percentile(dt_diffs, range(101))
    relative_percentiles = np.percentile(dt_diffs / dts, range(101))
    # plt.show(plt.plot(percentiles))
    # plt.show(plt.plot(relative_percentiles))
    # x = dt_diffs > percentiles[10]
    # x &= dt_diffs < percentiles[90]
    x = dt_diffs / dts > relative_percentiles[10] # TODO adjustable parameters!
    x &= dt_diffs / dts < relative_percentiles[90] # TODO adjustable parameters!
    x &= dts < 190 # TODO adjustable parameters!
    x &= dts > 50 # TODO adjustable parameters!
    # plt.show(plt.scatter(mzs[x], dt_diffs[x] / dts[x], c=dt_diffs[x] / dts[x], cmap="RdYlGn", marker="."))
    # plt.show(plt.scatter(mzs[x], mzs[x] / dts[x], c=dt_diffs[x] / dts[x], cmap="RdYlGn", marker="."))
    # plt.show(plt.scatter(dts[x], mzs[x], c=dt_diffs[x] / dts[x], cmap="RdYlGn", marker="."))
    xdata = np.stack([dts, mzs])
    ydata = 1000000 * dt_diffs / dts
    # plt.show(plt.scatter(*xdata[:, x], c=ydata[x], cmap="RdYlGn", marker="."))
    popt, pcov = curve_fit(func, xdata[:, x][:,::2], ydata[x][::2])
    shift_parameters.append(popt)
    sample_indices = ions["SAMPLE"] == sample
    sample_indices &= ions["LE"]
    pred = func(alldata[:, sample_indices], *shift_parameters[sample])
    new_pred[sample_indices] = pred
    # ypred = func(xdata, *popt)
    # scipy.stats.pearsonr(ydata[x], ypred[x])
    # # plt.show(plt.scatter(ydata[x], ypred[x], marker="."))
    # # plt.show(sns.jointplot(ydata[x], ypred[x], kind="hex"))
    # plt.show(plt.scatter(*xdata[:, x], c=ypred[x], cmap="RdYlGn", marker="."))
    # relative_dt_errors = (ypred[x][1::2] - ydata[x][1::2])
    # absolute_dt_errors = relative_dt_errors / 1000000 * dts[x][1::2]
    # plt.show(plt.scatter(dts[x][1::2], relative_dt_errors, marker="."))
    # plt.show(sns.jointplot(dts[x][1::2], relative_dt_errors, kind="hex"))
    # plt.show(plt.boxplot(absolute_dt_errors))
    # p = np.percentile(relative_dt_errors, range(101))
    # second_derivative_p = np.gradient(np.gradient(p))
    # plt.show(plt.plot(p, marker="."))
    # new_dt = dts[x] + ypred[x] * dts[x]/1000000
    # new_dt_diff = (new_dt - ions["DT"][sample_ions[:, 1]][x]) / new_dt * 1000000
    # p = np.percentile(new_dt_diff, range(101))
    # plt.show(plt.plot(p, marker="."))


new_dts = ions["DT"] * (1 + new_pred / 1000000)
anchor_dts = anchor_ions.copy()
anchor_dts.data = new_dts[anchor_dts.data]
anchor_dts = anchor_dts.sum(axis=1).squeeze().A.squeeze() / anchors["ION_COUNT"]

dt_diffs = (anchors["DT"] - anchor_dts) / anchors["DT"]
x = dt_diffs < 0.04
x &= dt_diffs > 0
x &= anchors["LE"]
x = np.flatnonzero(x)[::10]
plt.show(plt.scatter(anchor_dts[x], anchors["MZ"][x], c=dt_diffs[x], cmap="RdYlGn", marker="."))
plt.show(plt.scatter(anchor_dts[x], dt_diffs[x], c=dt_diffs[x], cmap="RdYlGn", marker="."))










































annotation_array = np.array(annotation_data)
# annotation_array = np.array(annotation_data[:120000])
fps = np.cumsum(annotation_array[:, 1].astype(int)!=1)
fdr = fps / np.arange(1, annotation_array.shape[0] + 1)

fdr_threshold = 0.01

annotation_array = annotation_array[:1 + np.flatnonzero(fdr < fdr_threshold)[-1]]
annotation_array = annotation_array[annotation_array[:, 1].astype(int)==1]

scores = annotation_array[:, -4].astype(float)
organisms = np.array([x.split("_")[-1] if "_" in x else x for x in annotation_array[:, -1]])
anchor_indices = annotation_array[:, 2].astype(int)
organism_logfcs = logfcs[anchor_indices]
qc_intensity = anchors_per_condition["QC"]["INTENSITY"][anchor_indices]
# qc_intensity = anchors_per_condition["A"]["INTENSITY"][anchor_indices]
organism_map = np.array(["YEAST", "HUMAN", "ECOLI", "UNKNOWN"])
predicted_organism = np.array(
    [
        0 if i > .5 else (2 if i < -1 else 1) for i in organism_logfcs
    ]
)
color_palette = np.array(["red", "blue", "green", "black"])
color = np.array([np.flatnonzero(organism_map == i)[0] if i in organism_map else 3 for i in organisms])
for c in np.unique(color)[2:3]:
    organism = organism_map[c]
    # if organism == "UNKNOWN":
    #     continue
    elements = color == c
    tmp = plt.scatter(
        np.log2(qc_intensity[elements]),
        organism_logfcs[elements],
        marker=".",
        c=color_palette[c],
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

plt.axhline(.5)
plt.axhline(-1)
plt.axvline(15)
plt.legend()
plt.xlabel("Log(QC))")
plt.ylabel("LogFC(A/B)")
plt.show()









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


x = color != 3
x &= ~np.isnan(organism_logfcs)
x &= ~np.isinf(organism_logfcs)
np.sum(color[x] != clusters[x]) / np.sum(color[x] == clusters[x])

ax = plt.subplots(1, 4)[1]
ax[3].scatter(np.log2(qc_intensity), organism_logfcs, c=clusters, marker=".")
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

























a, b = neighbors.nonzero()
c = anchors["LE"][a] != anchors["LE"][b]
n = scipy.sparse.csr_matrix(
    (
        neighbors.data[c],
        (
            a[c],
            b[c]
        )
    ),
    shape=neighbors.shape,
    dtype=np.bool
)

# n = n[anchors["LE"]]


spectra = {
    p: n.indices[
        n.indptr[p]: n.indptr[p + 1]
    ] for p in np.flatnonzero(
        anchors["LE"] & (np.diff(n.indptr) > 10)
    )
}

# for p, f in spectra.items():
#     anchors[p]
#     np.sort(anchors[f])
#     neighbors[np.ix_(f,f)].nnz / len(f) ** 2
#     input()


anchor_intensities = anchor_ions.copy()
anchor_intensities.data = ions["CALIBRATED_INTENSITY"][anchor_intensities.data]
anchor_intensities = anchor_intensities.sum(axis=1).A.squeeze() / anchors["ION_COUNT"]
with open(parameters["OUTPUT_PATH"] + "/test.mgf", "w") as outfile:
    for spectrum_index, spectrum in sorted(spectra.items()):
        tmp = outfile.write("BEGIN IONS\n")
        tmp = outfile.write(
            "TITLE=(index_{})(dt_{})(rt_{})\n".format(
                spectrum_index,
                anchors["DT"][spectrum_index],
                anchors["RT"][spectrum_index],
            )
        )
        tmp = outfile.write(
            "SCANS={}\n".format(
                spectrum_index
            )
        )
        tmp = outfile.write(
            "RTINSECONDS={}\n".format(
                anchors["RT"][spectrum_index] * 60
            )
        )
        tmp = outfile.write(
            "PEPMASS={}\n".format(
                # precursor_mz
                # 1000
                anchors["MZ"][spectrum_index]
            )
        )  # Defaults
        for anchor_index in spectrum:
            tmp = outfile.write(
                "{} {}\n".format(
                    anchors["MZ"][anchor_index],
                    anchor_intensities[anchor_index]
                    # raw_anchors[neighbor_anchor_index, ANCHOR_ION_COUNT]
                )
            )
        tmp = outfile.write("END IONS\n")
























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































int_file_name = "/home/sander/Proteomics/histopya2/tmp/random/PeakView_TTOF6600_64w_shift_iRT_extractionWindow10min_30ppm.csv"
fdr_file_name = "/home/sander/Proteomics/histopya2/tmp/random/PeakView_TTOF6600_64w_shift_iRT_extractionWindow10min_30ppm_FDRs.csv"
fdrs = pd.read_csv(fdr_file_name).values
intensities = pd.read_csv(int_file_name).values

query_fdrs = fdrs[:, -6:].astype(float)
peps = fdrs[np.any(query_fdrs, axis=1) < 0.01, 1].astype(str)

target_peptides = intensities[:, 1].astype(str)
targets = np.isin(target_peptides, peps)

target_prots = intensities[targets,0].astype(str)
target_species = np.array(
    [
        i.split("_")[-1] if "_" in i else "UNKNOWN" for i in target_prots
    ]
)
target_peps = target_peptides[targets]
target_intensities = intensities[targets, -6:].astype(float)

organism_logfcs = np.log2(target_intensities[:, 0] / target_intensities[:,3])
qc_intensity = np.log2((target_intensities[:, 0] + target_intensities[:,3]) / 2)

organisms = target_species
# qc_intensity = anchors_per_condition["A"]["INTENSITY"][anchor_indices]
organism_map = np.array(["YEAS8", "HUMAN", "ECOLI", "UNKNOWN"])

color_palette = np.array(["red", "blue", "green", "black"])
color = np.array([np.flatnonzero(organism_map == i)[0] if i in organism_map else 3 for i in organisms])
for c in np.unique(color)[::-1]:
    organism = organism_map[c]
    if organism == "UNKNOWN":
        continue
    elements = color == c
    tmp = plt.scatter(
        np.log2(qc_intensity[elements]),
        organism_logfcs[elements],
        marker=".",
        c=color_palette[c],
        label=organism
    )


plt.legend()
plt.xlabel("Log(QC))")
plt.ylabel("LogFC(A/B)")
plt.show()

predicted_organism = np.array(
    [
        0 if i > .5 else (2 if i < -1 else 1) for i in organism_logfcs
    ]
)
z=predicted_organism==color
np.sum(~z)/np.sum(z)












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



neighbor_threshold = parameters["NEIGHBOR_THRESHOLD"]
sample_count = parameters["SAMPLE_COUNT"]
# percentile_limit = parameters["ANCHOR_ALIGNMENT_PERCENTILE_THRESHOLD"]
percentile_limit = (
    1 - (
        1 - scipy.stats.norm.cdf(
            parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"]
        )
    ) * 2
) ** 2 # DT and RT are indepently done
minimum_overlap = parameters["MINIMUM_OVERLAP"]
minimum_hits = [minimum_overlap]
for total in range(1, sample_count + 1):
    select = total
    target = binom(total, select)
    target *= percentile_limit**select
    target *= (1 - percentile_limit)**(total - select)
    while target < neighbor_threshold:
    # while select > 0:
        select -= 1
        new_target = binom(total, select)
        new_target *= percentile_limit**select
        new_target *= (1 - percentile_limit)**(total - select)
        target += new_target
    minimum_hits.append(max(select, minimum_overlap))

for i in enumerate(minimum_hits): i
