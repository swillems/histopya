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

# Initializing
parameter_file_name = "data/lfq_swim_udmse_combined/parameters_QC.json"
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

# Plotting cvs
with log.newSection("Plotting cvs"):
    full_anchor_ions = (
        anchor_ions[anchors["ION_COUNT"] == parameters["SAMPLE_COUNT"]]
    ).todense().A
    cints = ions["CALIBRATED_INTENSITY"][full_anchor_ions]
    # rts = ions["CALIBRATED_RT"][full_anchor_ions]
    # dts = ions["CALIBRATED_DT"][full_anchor_ions]
    # mzs = ions["CALIBRATED_MZ"][full_anchor_ions]
    # rts = (rts - np.average(rts, axis=1).reshape(-1, 1))
    # rts /= np.std(rts, axis=1).reshape(-1, 1)
    # dts = (dts - np.average(dts, axis=1).reshape(-1, 1))
    # dts /= np.std(dts, axis=1).reshape(-1, 1)
    # mzs = (mzs - np.average(mzs, axis=1).reshape(-1, 1))
    # mzs /= np.std(mzs, axis=1).reshape(-1, 1)
    # errs = np.sqrt(rts**2 + dts**2 + mzs**2)
    # s = np.all(errs < 4, axis=1)
    # cints_swim = cints[s, :9]
    # cints_udmse = cints[s, 9:]
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
            columns=["SWIM-DIA", "UDMSE"]
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
    tmp = plt.ylabel("Distribution")
    tmp = plt.xlabel("CV of fully reproducible aggregates")
    tmp = plt.yticks([])
    tmp = plt.xlim([0, 0.5])
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + "cv_comparison.pdf", bbox_inches='tight')
    tmp = plt.close()
