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
import os


# Initializing
parameter_file_name = "projects/hdmse_swim_qc/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "_interactive.txt")
parameters["PLOTS_PATH"] = parameters["OUTPUT_PATH"] + "figures/"
if not os.path.exists(parameters["PLOTS_PATH"]):
    os.makedirs(parameters["PLOTS_PATH"])

# Loading data
with log.newSection("Loading SWIM/HDMSE data"):
    anchor_ions = src.io.loadMatrix(
        "ANCHOR_IONS_FILE_NAME",
        parameters,
    )
    ions = src.io.loadArray("IONS_FILE_NAME", parameters)


# Plotting cvs
with log.newSection("Plotting SWIM/HDMSE cvs"):
    full_anchor_ions = (
        anchor_ions[np.diff(anchor_ions.indptr) == parameters["SAMPLE_COUNT"]]
    ).todense().A
    cints = ions["CALIBRATED_INTENSITY"][full_anchor_ions]
    cints_swim = cints[:, :9]
    cints_udmse = cints[:, 9:]
    cints_swim_cv = scipy.stats.variation(cints_swim, axis=1)
    cints_udmse_cv = scipy.stats.variation(cints_udmse, axis=1)
    results = scipy.stats.ttest_rel(cints_swim_cv, cints_udmse_cv)
    log.printMessage(
        "SWIM/UDMSE median cvs (ttest: {}, pval: {}): {} {}".format(
            results[0],
            results[1],
            np.median(cints_swim_cv),
            np.median(cints_udmse_cv),
        )
    )
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




# Plotting mgf edge COUNTS
parameter_file_name = "projects/dda/parameters.json"
# parameter_file_name = "data/searle_hela_dda/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "interactive.txt")
parameters["PLOTS_PATH"] = parameters["OUTPUT_PATH"] + "figures/"
if not os.path.exists(parameters["PLOTS_PATH"]):
    os.makedirs(parameters["PLOTS_PATH"])
with log.newSection("DDA mgf peak count plotting"):
    spectra = src.io.loadArray("IONS_FILE_NAME", parameters)
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
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "lfq_mgf_peak_counts.pdf", bbox_inches='tight')
    tmp = plt.close()


























# Initializing
# extension = "tenzer"
# extension = "udmse"
# extension = "swim"
# extension = "dda"
extension = "tenzer"
if extension == "tenzer":
    parameter_file_name = "projects/tenzer/parameters.json"
elif extension == "udmse":
    parameter_file_name = "projects/udmse/parameters.json"
elif extension == "swim":
    parameter_file_name = "projects/swim/parameters.json"
elif extension == "sonar":
    parameter_file_name = "projects/sonar/parameters.json"
elif extension == "dda":
    parameter_file_name = "projects/dda/parameters.json"


parameters = src.parameters.getDefaultParameters()
parameters = src.parameters.updateParameters(
    parameters,
    parameter_file_name
)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "_" + extension + "_lfq.txt")

with log.newSection("Loading ion-network annotation"):
    database = src.peptides.loadDatabase(parameters["DATABASE_FILE_NAME"])
    base_mass_dict = database["base_mass_dict"]
    proteins = database["proteins"]
    total_protein_sequence = database["total_protein_sequence"]
    peptides = database["peptides"]
    peptide_index_matrix = database["peptide_index_matrix"]
    peptide_masses = database["peptide_masses"]
    fragments = database["fragments"]
    anchors = src.io.loadArray("ANCHORS_FILE_NAME", parameters)
    anchor_peptide_scores = src.io.loadMatrix(
        "ANCHOR_PEPTIDE_SCORES_FILE_NAME",
        parameters,
    )
    anchor_peptide_match_counts = src.io.loadMatrix(
        "ANCHOR_PEPTIDE_MATCH_COUNTS_FILE_NAME",
        parameters,
    )
    neighbors = src.io.loadMatrix(
        "ANCHOR_NEIGHBORS_FILE_NAME",
        parameters,
    )
    neighbors += neighbors.T


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


with log.newSection("Calculating chimericy"):
    n = neighbors[significant_anchors].T.tocsr()[significant_anchors]
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
    aa_peps = np.array(
        [
            i.replace("I", "L") for i in aa_peps
        ]
    )
    bb_peps = np.array(
        [
            i.replace("I", "L") for i in bb_peps
        ]
    )
    chimeric = np.flatnonzero(aa_peps != bb_peps)
    log.printMessage(
        "Full neighbor 25, 50, 75 percentiles: {}, {}, {}".format(
            *map(int, np.percentile(np.diff(neighbors.indptr), [25, 50, 75]))
        )
    )
    log.printMessage(
        "Annotated neighbor 25, 50, 75 percentiles: {}, {}, {}".format(
            *map(int, np.percentile(np.diff(n.indptr), [25, 50, 75]))
        )
    )
    log.printMessage("Annotated edges: {} ({})".format(n.nnz, n.nnz / neighbors.nnz))
    log.printMessage("Chimericy count: {} ({})".format(len(chimeric), len(chimeric) / n.nnz))
    enrichment = np.bincount(anchors["ION_COUNT"][np.unique(significant_anchors)])
    total = np.bincount(anchors["ION_COUNT"])
    log.printMessage(
        "Annotation enrichment: two-fold {:.2f}%, fully {:.2f}%".format(
            100 * enrichment[2] / total[2],
            100 * enrichment[-1] / total[-1]
        )
    )
