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
import matplotlib
import seaborn as sns


# Initializing
parameter_file_name = "data/lfq_swim_udmse_combined/parameters_QC.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "_interactive.txt")


# Loading data
with log.newSection("Loading SWIM/HDMSE data"):
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


# Plotting cvs
with log.newSection("Plotting SWIM/HDMSE cvs"):
    full_anchor_ions = (
        anchor_ions[anchors["ION_COUNT"] == parameters["SAMPLE_COUNT"]]
    ).todense().A
    cints = ions["CALIBRATED_INTENSITY"][full_anchor_ions]
    cints_swim = cints[:, :9]
    cints_udmse = cints[:, 9:]
    cints_swim_cv = scipy.stats.variation(cints_swim, axis=1)
    cints_udmse_cv = scipy.stats.variation(cints_udmse, axis=1)
    scipy.stats.ttest_rel(cints_swim_cv, cints_udmse_cv)
    d = pd.melt(
        pd.DataFrame(
            np.stack(
                [
                    cints_swim_cv,
                    cints_udmse_cv,
                ]
            ).T,
            columns=["SWIM-DIA", "HDMSE"]
        ),
    )
    d["Y"] = 1
    d["Acquistion"] = d["variable"]
    tmp = sns.violinplot(
        x='value',
        y='Y',
        hue='Acquistion',
        split=True,
        data=d,
        inner="quartile",
        gridsize=1000,
        orient="h"
    )
    tmp = plt.ylabel("Relative Frequency")
    tmp = plt.xlabel("CV Of Fully Reproducible Aggregates")
    tmp = plt.yticks([])
    tmp = plt.xlim([0, 0.5])
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "cv_comparison.pdf", bbox_inches='tight')
    tmp = plt.close()


# Loading mgf spectra
parameter_file_name = "projects/mgf_dda_parameters.json"
# parameter_file_name = "data/searle_hela_dda/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "interactive.txt")
with log.newSection("DDA ion-network loading"):
    mgf_anchors = src.io.loadArray("ANCHORS_FILE_NAME", parameters)
    spectra = src.io.loadArray("IONS_FILE_NAME", parameters)
    mgf_neighbors = src.io.loadMatrix(
        "ANCHOR_NEIGHBORS_FILE_NAME",
        parameters,
    )



# Plotting mgf edge COUNTS
# with log.newSection("Plotting mgf peaks"):
with log.newSection("DDA mgf peak count plotting"):
    # spectrum_sizes = np.repeat(spectrum_sizes, spectrum_sizes)
    a, b = np.unique(spectra["PEAK_COUNT"], return_counts=True)
    fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [5, 1]})
    tmp = plt.subplots_adjust(hspace=0.1)
    tmp = ax[0].scatter(a, np.log2(b), marker=".")
    tmp = ax[0].set_ylabel("Log2(Spectrum Frequency)")
    tmp = ax[1].boxplot(spectra["PEAK_COUNT"], whis="range", vert=False, widths=0.5)
    tmp = ax[1].set_yticks([])
    ax[1].get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )
    tmp = ax[1].set_xlabel("Peak Count")
    # tmp = plt.show()
    tmp = plt.savefig(parameters["OUTPUT_PATH"] + "lfq_mgf_peak_counts.pdf", bbox_inches='tight')
    tmp = plt.close()


with log.newSection("Loading mgf DDA ion-network annotation"):
    database = src.peptides.loadDatabase(parameters["DATABASE_FILE_NAME"])
    base_mass_dict = database["base_mass_dict"]
    proteins = database["proteins"]
    total_protein_sequence = database["total_protein_sequence"]
    # ptms = database["ptms"]
    # ptm_matrix = database["ptm_matrix"]
    peptides = database["peptides"]
    peptide_index_matrix = database["peptide_index_matrix"]
    # digestion_matrix = database["digestion_matrix"]
    peptide_masses = database["peptide_masses"]
    fragments = database["fragments"]
    # mgf_anchor_boundaries, mgf_fragment_peptide_indices, mgf_fragment_indices = src.aggregates.matchAnchorsToFragments(
    #     fragments,
    #     mgf_anchors,
    #     base_mass_dict,
    #     parameters,
    #     log
    # )
    mgf_anchor_peptide_scores = src.io.loadMatrix(
        "ANCHOR_PEPTIDE_SCORES_FILE_NAME",
        parameters,
    )
    mgf_anchor_peptide_match_counts = src.io.loadMatrix(
        "ANCHOR_PEPTIDE_MATCH_COUNTS_FILE_NAME",
        parameters,
    )
    # mgf_anchor_fragment_indices = src.aggregates.getAnchorFragmentIndices(
    #     mgf_anchor_peptide_match_counts,
    #     mgf_anchor_boundaries,
    #     mgf_fragment_indices,
    #     mgf_fragment_peptide_indices,
    #     parameters,
    #     log
    # )
    # mgf_precursor_indices = src.aggregates.findFragmentPrecursors(
    #     mgf_anchor_peptide_match_counts,
    #     mgf_anchors,
    #     mgf_neighbors,
    #     # anchor_alignment_parameters,
    #     None,
    #     peptide_masses,
    #     base_mass_dict,
    #     parameters,
    #     log
    # )
    # annotation_data = src.aggregates.writePercolatorFile(
    #     mgf_anchors,
    #     base_mass_dict,
    #     mgf_anchor_peptide_match_counts,
    #     fragments,
    #     mgf_anchor_fragment_indices,
    #     mgf_neighbors,
    #     peptides,
    #     proteins,
    #     peptide_masses,
    #     mgf_precursor_indices,
    #     mgf_anchor_peptide_scores,
    #     peptide_index_matrix,
    #     total_protein_sequence,
    #     parameters,
    #     log
    # )
    # src.io.runPercolator(parameters, log)







# Annotation accuracy
with log.newSection("Calculating aggregate annotation accuracy"):
    percolated_pims = pd.read_csv(
        parameters["PERCOLATOR_TARGET_PIMS"],
        delimiter="\t"
    )
    percolated_pim_fdrs = percolated_pims.values[:, 2]
    anchor_pims, peptide_pims = mgf_anchor_peptide_match_counts.nonzero()
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



if True:
    percolated_peptides = pd.read_csv(
        parameters["PERCOLATOR_TARGET_PEPTIDES"],
        delimiter="\t"
    )
    percolated_peptide_fdrs = percolated_peptides.values[:, 2]
    peptide_fdr = 1
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
    protein_fdr = 1
    significant_percolated_proteins = percolated_proteins.values[
        percolated_protein_fdrs <= protein_fdr,
        0
    ]
    # unique_organisms = ["ECOLI", "YEAST", "HUMAN"]
    select = np.isin(significant_peptide_sequences, significant_percolated_peptides)
    select &= np.isin(significant_proteins, significant_percolated_proteins)
    # select &= np.isin(significant_organisms, unique_organisms)
    significant_anchors = significant_anchors[select]
    significant_peptides = significant_peptides[select]
    significant_proteins = significant_proteins[select]
    significant_organisms = significant_organisms[select]
    significant_peptide_sequences = significant_peptide_sequences[select]
    significant_percolated_pims = significant_percolated_pims[select]




# parameters["MGF_FILE"] = parameters["OUTPUT_PATH"] + "mgf_annotations.csv"
# with log.newSection("Writing results"):
#     data = []
#     header = [
#         "spectrum_index",
#         "spectrum_rt",
#         "spectrum_mz",
#         "spectrum_z",
#         "spectrum_peak_count",
#         "anchor_index",
#         "anchor_mz",
#         "anchor_intensity",
#         "pim_index",
#         "pim_match_count",
#         "pim_score",
#         "fragment_index",
#         "fragment_type",
#         "fragment_site",
#         "fragment_mr",
#         "peptide_sequence",
#         "peptide_mr",
#     ]
#     for index in np.argsort(mgf_anchors[significant_anchors]["ION_COUNT"]):
#         anchor_index = significant_anchors[index]
#         anchor = mgf_anchors[anchor_index]
#         spectrum_index = anchor["ION_COUNT"]
#         spectrum = spectra[spectrum_index]
#         peptide_index = significant_peptides[index]
#         peptide_sequence = significant_peptide_sequences[index]
#         peptide_mass = peptide_masses[peptide_index]
#         pim_index = significant_percolated_pims[index]
#         fragment_index = mgf_anchor_fragment_indices[pim_index][0]
#         if fragment_index > 0:
#             fragment = fragments[fragment_index]
#             fragment_type = "y"
#             fragment_site = fragment["Y_INDEX"]
#             fragment_mr = fragment["Y_MR"]
#         else:
#             fragment = fragments[-fragment_index]
#             fragment_type = "b"
#             fragment_site = fragment["B_INDEX"]
#             fragment_mr = fragment["B_MR"]
#         row = [
#             spectrum_index,
#             spectrum["RT"],
#             spectrum["MZ"],
#             spectrum["Z"],
#             spectrum["PEAK_COUNT"],
#             anchor_index,
#             anchor["MZ"],
#             anchor["SHIFTED_DT"],
#             pim_index,
#             mgf_anchor_peptide_match_counts.data[pim_index],
#             mgf_anchor_peptide_scores.data[pim_index],
#             fragment_index,
#             fragment_type,
#             fragment_site,
#             fragment_mr,
#             peptide_sequence,
#             peptide_mass,
#         ]
#         data.append(row)
#     src.io.saveListOfListsToCsv(
#         data,
#         "MGF_FILE",
#         parameters,
#         log,
#         header,
#         "\t"
#     )
















# chimericy in annotations
n = mgf_neighbors[significant_anchors].T.tocsr()[significant_anchors]
a, b = n.nonzero()
aa = significant_peptides[a]
bb = significant_peptides[b]
aa_peps = src.peptides.getPeptideSequences(
    aa,
    peptides,
    peptide_index_matrix,
    total_protein_sequence
)
bb_peps = src.peptides.getPeptideSequences(
    bb,
    peptides,
    peptide_index_matrix,
    total_protein_sequence
)
z = np.flatnonzero(aa_peps != bb_peps)
zz = np.stack(
    [
        aa_peps[z],
        bb_peps[z],
    ]
).T
#
# spec_ans = [[] for s in spectra]
# for anch, pep in zip(significant_anchors, significant_peptides):
#     spec = mgf_anchors[anch]["ION_COUNT"]
#     spec_ans[spec].append(pep)
#
# spec_ans = np.array([np.unique(s, return_counts=True)[1] for s in spec_ans])
#
# uni = []
# unk = []
# mis = []
# for i, s in enumerate(spec_ans):
#     c = spectra[i]["PEAK_COUNT"]
#     t = np.sum(s)
#     for j in s:
#         m = t - j
#         uni += [j / c] * j
#         unk += [(c - t) / c] * j
#         mis += [m / c] * j
#
#
# uni = np.array(uni)
# unk = np.array(unk)
# mis = np.array(mis)
