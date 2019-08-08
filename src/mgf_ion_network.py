#!venv/bin/python


import src.parameters
import src.io
import src.ions
import src.aggregates
import src.peptides
import numpy as np
import scipy.sparse
from pyteomics import mgf


# Loading mgf spectra
parameter_file_name = "projects/mgf_dda_parameters.json"
mgf_file_name = 'data/lfq_qc_mgf/lfq_qc_dda.mgf'
# parameter_file_name = "data/searle_hela_dda/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "interactive.txt")
with log.newSection("DDA ion-network node creation"):
    raw_spectra = [s for s in mgf.read(mgf_file_name)]
    # raw_spectra = [s for s in mgf.read('data/searle_hela_dda/23aug2017_hela_serum_timecourse_pool_dda_001.mgf')]
    spectra = np.empty(
        len(raw_spectra),
        dtype=[
            ("RT", np.float),
            ("MZ", np.float),
            ("Z", np.int),
            ("PEAK_COUNT", np.int),
        ]
    )
    spectra['RT'] = [
        spectrum["params"]["rtinseconds"] / 60 for spectrum in raw_spectra
    ]
    spectra['MZ'] = [
        spectrum["params"]["pepmass"][0] for spectrum in raw_spectra
    ]
    spectra['Z'] = [
        spectrum["params"]["charge"][0] for spectrum in raw_spectra
    ]
    spectra["PEAK_COUNT"] = np.array([len(spectrum['m/z array']) for spectrum in raw_spectra])
    src.io.saveArray(spectra, "IONS_FILE_NAME", parameters, log)
    mgf_anchors = np.empty(
        np.sum(spectra["PEAK_COUNT"]),
        dtype=[
            ("MZ", np.float),
            ("RT", np.float),
            ("DT", np.float),
            ("SHIFTED_DT", np.float),
            ("LE", np.bool),
            ("ION_COUNT", np.int),
        ]
    )
    mgf_anchors["MZ"] = np.concatenate(
        [
            spectrum["m/z array"] for spectrum in raw_spectra
        ]
    )
    mgf_anchors["RT"] = np.repeat(
        [
            spectrum["params"]["rtinseconds"] / 60 for spectrum in raw_spectra
        ],
        spectra["PEAK_COUNT"]
    )
    mgf_anchors["DT"] = np.repeat(
        [
            spectrum["params"]["pepmass"][0] for spectrum in raw_spectra
        ],
        spectra["PEAK_COUNT"]
    )
    mgf_anchors["SHIFTED_DT"] = np.concatenate(
        [
            spectrum["intensity array"] for spectrum in raw_spectra
        ]
    )
    mgf_anchors["LE"] = 0
    mgf_anchors["ION_COUNT"] = np.repeat(np.arange(len(spectra)), spectra["PEAK_COUNT"])
    src.io.saveArray(
        mgf_anchors,
        "ANCHORS_FILE_NAME",
        parameters,
        log
    )


with log.newSection("DDA ion-network edge creation"):
    mgf_neighbors = scipy.sparse.csr_matrix(
        scipy.sparse.block_diag(
            [
                np.ones(shape=(s, s), dtype=np.bool) for s in spectra["PEAK_COUNT"]
            ]
        )
    )
    src.io.saveMatrix(
        mgf_neighbors,
        "ANCHOR_NEIGHBORS_FILE_NAME",
        parameters,
        log
    )


with log.newSection("DDA ion-network annotation"):
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
    mgf_anchor_boundaries, mgf_fragment_peptide_indices, mgf_fragment_indices = src.aggregates.matchAnchorsToFragments(
        fragments,
        mgf_anchors,
        base_mass_dict,
        parameters,
        log
    )
    mgf_anchor_peptide_scores, mgf_anchor_peptide_match_counts = src.aggregates.getAnchorPeptideMatrix(
        mgf_anchors,
        mgf_neighbors,
        peptides,
        mgf_anchor_boundaries,
        mgf_fragment_peptide_indices,
        parameters,
        log
    )
    mgf_anchor_fragment_indices = src.aggregates.getAnchorFragmentIndices(
        mgf_anchor_peptide_match_counts,
        mgf_anchor_boundaries,
        mgf_fragment_indices,
        mgf_fragment_peptide_indices,
        parameters,
        log
    )
    mgf_precursor_indices = src.aggregates.findFragmentPrecursors(
        mgf_anchor_peptide_match_counts,
        mgf_anchors,
        mgf_neighbors,
        # anchor_alignment_parameters,
        None,
        peptide_masses,
        base_mass_dict,
        parameters,
        log
    )
    annotation_data = src.aggregates.writePercolatorFile(
        mgf_anchors,
        base_mass_dict,
        mgf_anchor_peptide_match_counts,
        fragments,
        mgf_anchor_fragment_indices,
        mgf_neighbors,
        peptides,
        proteins,
        peptide_masses,
        mgf_precursor_indices,
        mgf_anchor_peptide_scores,
        peptide_index_matrix,
        total_protein_sequence,
        parameters,
        log
    )
    src.io.runPercolator(parameters, log)
