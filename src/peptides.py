#!venv/bin/python


import src.io
import numpy as np
import scipy.sparse
from collections import defaultdict


def importProteinsAndPtms(parameters, log, generate_decoy=True):
    with log.newSection("Reading protein databases"):
        sequence, protein_ids, protein_sizes, protein_ptms = __defineAminoAcidSequence(
            parameters, log, generate_decoy
        )
        proteins = __defineProteins(
            sequence,
            protein_ids,
            protein_sizes,
            generate_decoy,
            parameters,
            log
        )
        ptms, ptm_matrix = __definePTMs(
            protein_ptms,
            sequence,
            generate_decoy,
            parameters,
            log
        )
    return proteins, sequence, ptms, ptm_matrix


def __defineAminoAcidSequence(parameters, log, generate_decoy):
    with log.newSection("Reading amino acid sequence"):
        proteins = src.io.loadJSON("DATABASE_FILE_NAMES", parameters, log)
        protein_ids = sorted(proteins.keys())
        protein_sizes = np.array([len(proteins[protein][0]) for protein in protein_ids])
        sequence = "".join(proteins[protein][0] for protein in protein_ids)
        protein_offsets = np.concatenate(
            [[0], np.cumsum(protein_sizes)[:-1]]
        )
        protein_ptms = [
            (loc + protein_offsets[i], ptm) for i, protein in enumerate(protein_ids) for loc, ptm in proteins[protein][1] if loc >= 0
        ]
        if generate_decoy:
            log.printMessage("Generating reversed decoy")
            sequence += sequence[::-1]
            protein_ids += [
                "DECOY_{}".format(protein) for protein in protein_ids
            ][::-1]
            protein_sizes = np.concatenate(
                [protein_sizes, protein_sizes[::-1]]
            )
            protein_ptms += [
                (len(sequence) - loc - 1, ptm) for loc, ptm in protein_ptms[::-1]
            ]
        log.printMessage(
            "Imported {} target amino acids".format(
                len(sequence) if not generate_decoy else len(sequence) // 2
            )
        )
    return sequence, protein_ids, protein_sizes, protein_ptms


def __defineProteins(
    sequence,
    protein_ids,
    protein_sizes,
    generate_decoy,
    parameters,
    log
):
    with log.newSection("Defining proteins"):
        proteins = np.empty(
            len(protein_ids),
            dtype=[
                ("ID", object),
                ("START_INDEX", np.int),
                ("SIZE", np.int),
                ("DECOY", np.bool)
            ]
        )
        proteins["ID"] = np.array(protein_ids, dtype=object)
        proteins["START_INDEX"] = np.concatenate(
            [[0], np.cumsum(protein_sizes)[:-1]]
        )
        proteins["SIZE"] = protein_sizes
        proteins["DECOY"] = 0
        if generate_decoy:
            proteins["DECOY"][len(proteins) // 2:] = 1
        log.printMessage(
            "Imported {} target proteins".format(
                len(proteins) // 2 if generate_decoy else len(proteins)
            )
        )
    return proteins


def __definePTMs(
    protein_ptms,
    sequence,
    generate_decoy,
    parameters,
    log
):
    with log.newSection("Defining PTMs"):
        ptm_dict = defaultdict(list)
        for loc, ptm in protein_ptms:
            ptm_dict[ptm].append(loc)
        ptm_dict = {
            ptm: np.unique(ptm_locs) for ptm, ptm_locs in ptm_dict.items()
        }
        ptm_definitions = src.io.loadJSON("PTMS_FILE_NAME", parameters, log)
        ptms = np.empty(
            len(ptm_dict),
            dtype=[
                ("ID", object),
                ("MASS", np.float)
            ]
        )
        ptms["ID"] = np.array(sorted(ptm_dict.keys()), dtype=object)
        ptms["MASS"] = np.array(
            [
                ptm_definitions[ptm]["MM"] if ptm in ptm_definitions else 0 for ptm in ptms["ID"]
            ]
        )
        ptms["MASS"][np.isnan(ptms["MASS"])] = 0
        ptms = ptms[ptms["MASS"] != 0]
        ptm_sizes = np.array([len(ptm_dict[ptm]) for ptm in ptms["ID"]])
        ptm_matrix = scipy.sparse.csr_matrix(
            (
                np.ones(np.sum(ptm_sizes)),
                (
                    np.repeat(
                        np.arange(len(ptm_sizes)),
                        ptm_sizes
                    ),
                    np.concatenate([ptm_dict[ptm] for ptm in ptms["ID"]])
                )
            ),
            dtype=np.bool
        )
        log.printMessage(
            "Imported {} potential unique PTM localizations".format(
                ptm_matrix.nnz
            )
        )
    return ptms, ptm_matrix


def digestProteins(
    proteins,
    sequence,
    ptm_matrix,
    parameters,
    log,
):
    with log.newSection("Digesting proteins"):
        digestion_aas = parameters["DIGESTION_AMINO_ACIDS"]
        digestion_matrix = __calculateDigestionMatrix(
            proteins,
            sequence,
            digestion_aas,
            log
        )
        peptides, peptide_index_matrix = __definePeptides(
            proteins,
            sequence,
            digestion_matrix,
            ptm_matrix,
            log
        )
    return peptides, peptide_index_matrix, digestion_matrix


def __calculateDigestionMatrix(proteins, sequence, digestion_aas, log):
    log.printMessage("Finding cleavage points")
    protein_n_terms = proteins["START_INDEX"]
    protein_c_terms = proteins["START_INDEX"] + proteins["SIZE"] - 1
    cleavage_c_terms = np.flatnonzero(
        np.isin(list(sequence), list(digestion_aas)).astype(np.int)
    )
    cleavage_c_terms = cleavage_c_terms[
        ~np.isin(cleavage_c_terms, protein_c_terms)
    ]
    cleavage_n_terms = cleavage_c_terms + 1
    digestion_list = [
        protein_n_terms,
        protein_c_terms,
        cleavage_n_terms,
        cleavage_c_terms,
    ]
    digestion_matrix = scipy.sparse.csr_matrix(
        (
            np.ones(np.sum([len(l) for l in digestion_list])),
            (
                np.repeat(
                    np.arange(4),
                    [len(l) for l in digestion_list]
                ),
                np.concatenate(digestion_list)
            )
        ),
        dtype=np.bool
    )
    return digestion_matrix


def __definePeptides(
    proteins,
    sequence,
    digestion_matrix,
    ptm_matrix,
    log
):
    peptide_start_indices = (digestion_matrix[0] + digestion_matrix[2]).indices
    peptide_end_indices = (digestion_matrix[1] + digestion_matrix[3]).indices + 1
    peptide_sequences = np.array(
        [
            sequence[i: j] for i, j in zip(
                peptide_start_indices,
                peptide_end_indices
            )
        ],
        dtype=object
    )
    peptides, peptide_index_matrix, peptide_labels = __defineUniquePeptides(
        peptide_sequences,
        peptide_start_indices,
        sequence,
        proteins,
        log
    )
    peptides = __calculatePeptidoforms(
        peptides,
        ptm_matrix,
        peptide_start_indices,
        peptide_end_indices,
        peptide_labels,
        log
    )
    log.printMessage(
        "Found {} unique target peptide backbone sequences, and {} decoys".format(
            np.sum(~peptides["DECOY"]),
            np.sum(peptides["DECOY"]),
        )
    )
    return peptides, peptide_index_matrix


def __defineUniquePeptides(
    peptide_sequences,
    peptide_start_indices,
    sequence,
    proteins,
    log
):
    log.printMessage("Defining unique peptides")
    unique_peptides, peptide_labels = np.unique(
        peptide_sequences,
        return_inverse=True,
    )
    peptides = np.empty(
        len(unique_peptides),
        dtype=[
            ("SIZE", np.int),
            ("DECOY", np.bool),
            ("PEPTIDOFORMS", np.int),
            ("OCCURENCES", np.int),
            ("PROTEIN", np.int),
        ]
    )
    peptides["SIZE"] = [len(l) for l in unique_peptides]
    peptide_index_matrix = scipy.sparse.csr_matrix(
        (
            np.ones(len(peptide_start_indices)),
            (
                peptide_labels,
                peptide_start_indices
            )
        ),
        shape=(len(peptides), len(sequence)),
        dtype=np.int
    )
    peptides["OCCURENCES"] = np.diff(peptide_index_matrix.indptr)
    min_decoy_index = proteins[len(proteins) // 2]["START_INDEX"]
    peptides["DECOY"] = peptide_index_matrix.indices[
        peptide_index_matrix.indptr[:-1]
    ] >= min_decoy_index
    peptide_protein_origins = np.repeat(-1, len(peptides))
    index_origins = np.searchsorted(
        proteins["START_INDEX"],
        peptide_index_matrix.indices,
        "right"
    ) - 1
    selected_peptides = np.flatnonzero(peptides["OCCURENCES"] == 1)
    peptide_protein_origins[selected_peptides] = index_origins[
        peptide_index_matrix.indptr[selected_peptides]
    ]
    peptides["PROTEIN"] = peptide_protein_origins
    return peptides, peptide_index_matrix, peptide_labels


def __calculatePeptidoforms(
    peptides,
    ptm_matrix,
    peptide_start_indices,
    peptide_end_indices,
    peptide_labels,
    log
):
    log.printMessage("Calculating peptidoforms")
    modifiable_indices = ptm_matrix.sum(axis=0).A.squeeze() + 1
    peptidiform_counts = np.array(
        [
            np.product(modifiable_indices[i:j]) for i, j in zip(
                peptide_start_indices,
                peptide_end_indices
            )
        ]
    )
    peptides["PEPTIDOFORMS"] = 0
    for peptide_label, peptidiform_count in zip(
        peptide_labels, peptidiform_counts
    ):
        peptides["PEPTIDOFORMS"][peptide_label] += peptidiform_count
    return peptides


def loadBaseMassDict(parameters, log):
    log.printMessage("Importing mass definitions")
    aa_masses = src.io.loadJSON("AMINO_ACID_FILE_NAME", parameters)
    for aa, fixed_ptm in parameters["FIXED_MODIFICATIONS"].items():
        aa_masses[aa] += fixed_ptm
    atom_masses = src.io.loadJSON("ATOM_FILE_NAME", parameters)
    n_terminal_mass = parameters["TERMINI"]["N"]
    c_terminal_mass = parameters["TERMINI"]["C"]
    if n_terminal_mass is None:
        n_terminal_mass = atom_masses["H"]
    if c_terminal_mass is None:
        c_terminal_mass = atom_masses["H"] + atom_masses["O"]
    base_mass_dict = {
        "aas": aa_masses,
        "atoms": atom_masses,
        "n_terminus": n_terminal_mass,
        "c_terminus": c_terminal_mass
    }
    return base_mass_dict


def calculateMasses(
    peptides,
    peptide_index_matrix,
    base_mass_dict,
    total_protein_sequence,
    parameters,
    log,
    terminal_index=0
):
    # TODO peptidoform (fragment masses
    with log.newSection("Calculating peptide masses and fragments"):
        unique_peptide_sequence = "".join(
            [
                total_protein_sequence[start: start + size] for start, size in zip(
                    peptide_index_matrix.indices[
                        peptide_index_matrix.indptr[:-1]
                    ],
                    peptides["SIZE"]
                )
            ]
        )
        cumulative_aa_masses = np.cumsum(
            (
                np.array(
                    [base_mass_dict["aas"][i] if i in base_mass_dict["aas"] else 9999 for i in unique_peptide_sequence]
                ) * 1000000
            ).astype(np.int64)
        )
        peptide_masses = __calculatePeptideMasses(
            peptides,
            cumulative_aa_masses,
            base_mass_dict,
            log
        )
        fragments = __defineFragments(
            unique_peptide_sequence,
            base_mass_dict,
            cumulative_aa_masses,
            peptides,
            parameters,
            log,
            terminal_index
        )
    return peptide_masses, fragments


def __calculatePeptideMasses(
    peptides,
    cumulative_aa_masses,
    base_mass_dict,
    log
):
    log.printMessage("Calculating peptide masses")
    peptide_masses_cumul = cumulative_aa_masses[
        np.cumsum(peptides["SIZE"]) - 1
    ]
    peptide_masses = np.diff(
        np.concatenate(
            [[0], peptide_masses_cumul]
        )
    )
    peptide_masses += int(base_mass_dict["n_terminus"] * 1000000) + int(base_mass_dict["c_terminus"] * 1000000)
    peptide_masses = peptide_masses.astype(np.float) / 1000000
    return peptide_masses


def __defineFragments(
    unique_peptide_sequence,
    base_mass_dict,
    cumulative_aa_masses,
    peptides,
    parameters,
    log,
    terminal_index
):
    with log.newSection("Defining fragments"):
        fragments = np.empty(
            len(unique_peptide_sequence),
            dtype=[
                ("AMINO_ACID", "|S1"),
                ("PEPTIDE", np.int),
                ("B_MR", np.float),
                ("Y_MR", np.float),
                ("B_INDEX", np.int),
                ("Y_INDEX", np.int),
            ]
        )
        fragments["AMINO_ACID"] = np.array(list(unique_peptide_sequence))
        fragments["PEPTIDE"] = np.repeat(
            np.arange(len(peptides)),
            peptides["SIZE"]
        )
        fragments["B_INDEX"] = __calculateBIndices(
            peptides,
            log,
            terminal_index
        )
        fragments["Y_INDEX"] = __calculateYIndices(
            peptides,
            log,
            terminal_index
        )
        fragments["B_MR"] = __calculateBMRs(
            cumulative_aa_masses,
            peptides,
            base_mass_dict,
            log,
            terminal_index
        )
        fragments["Y_MR"] = __calculateYMRs(
            cumulative_aa_masses,
            peptides,
            base_mass_dict,
            log,
            terminal_index
        )
    return fragments


def __calculateBIndices(peptides, log, terminal_index=None):
    log.printMessage("Setting fragment b-indices")
    indices = np.arange(np.sum(peptides["SIZE"]), dtype=np.int)
    indices -= np.repeat(
        np.concatenate(
            [
                [-1],
                np.cumsum(peptides["SIZE"][:-1]) - 1
            ]
        ),
        peptides["SIZE"]
    )
    if terminal_index is not None:
        indices[np.cumsum(peptides["SIZE"]) - 1] = terminal_index
    return indices


def __calculateYIndices(peptides, log, terminal_index=None):
    log.printMessage("Setting fragment y-indices")
    indices = np.arange(np.sum(peptides["SIZE"]), dtype=np.int)
    indices -= np.repeat(
        np.cumsum(peptides["SIZE"]),
        peptides["SIZE"]
    )
    indices *= -1
    if terminal_index is not None:
        indices[np.cumsum(peptides["SIZE"][:-1])] = terminal_index
        indices[0] = terminal_index
    return indices


def __calculateBMRs(
    cumulative_aa_masses,
    peptides,
    base_mass_dict,
    log,
    terminal_index=None
):
    log.printMessage("Setting fragment b-masses")
    masses = cumulative_aa_masses.copy()
    masses -= np.repeat(
        np.concatenate(
            [
                [0],
                masses[np.cumsum(peptides["SIZE"][:-1]) - 1]
            ]
        ),
        peptides["SIZE"]
    )
    masses += int(base_mass_dict["n_terminus"] * 1000000) - int(base_mass_dict["atoms"]["H"] * 1000000)
    masses = masses.astype(np.float) / 1000000
    if terminal_index is not None:
        masses[np.cumsum(peptides["SIZE"]) - 1] = terminal_index
    return masses


def __calculateYMRs(
    cumulative_aa_masses,
    peptides,
    base_mass_dict,
    log,
    terminal_index=None
):
    log.printMessage("Setting fragment y-masses")
    masses = np.cumsum(
        np.diff(
            np.concatenate([[0], cumulative_aa_masses])
        )[::-1]
    )
    masses -= np.repeat(
        np.concatenate(
            [
                [0],
                masses[np.cumsum(peptides["SIZE"][::-1][:-1]) - 1]
            ]
        ),
        peptides["SIZE"][::-1]
    )
    masses += int(base_mass_dict["c_terminus"] * 1000000) + int(base_mass_dict["atoms"]["H"] * 1000000)
    masses = masses.astype(np.float) / 1000000
    masses = masses[::-1]
    if terminal_index is not None:
        masses[np.cumsum(peptides["SIZE"][:-1])] = terminal_index
        masses[0] = terminal_index
    return masses


def getPeptideSequenceFromIndex(
    peptide_index,
    peptides,
    peptide_index_matrix,
    total_protein_sequence
):
    if peptide_index == -1:
        return ""
    peptide_start_index = peptide_index_matrix.indices[
        peptide_index_matrix.indptr[peptide_index]
    ]
    peptide_sequence = total_protein_sequence[
        peptide_start_index: peptide_start_index + peptides[peptide_index]["SIZE"]
    ]
    return peptide_sequence


def getPeptideSequences(
    peptide_indices,
    peptides,
    peptide_index_matrix,
    total_protein_sequence
):
    return np.array(
        [
            getPeptideSequenceFromIndex(
                peptide_index,
                peptides,
                peptide_index_matrix,
                total_protein_sequence
            ) for peptide_index in peptide_indices
        ],
    )


def getProteinAccessionFromIndex(
    peptide_index,
    peptides,
    proteins
):
    if peptide_index == -1:
        return ""
    protein_index = peptides[peptide_index]["PROTEIN"]
    if protein_index != -1:
        protein_string = proteins[protein_index]["ID"]
    else:
        protein_string = "Ambiguous"
    return protein_string


def getProteinAccessions(
    peptide_indices,
    peptides,
    proteins
):
    return np.array(
        [
            getProteinAccessionFromIndex(
                peptide_index,
                peptides,
                proteins
            ) for peptide_index in peptide_indices
        ],
    )


def createDatabase(database_parameters, file_name, log):
    base_mass_dict = src.peptides.loadBaseMassDict(database_parameters, log)
    proteins, total_protein_sequence, ptms, ptm_matrix = src.peptides.importProteinsAndPtms(
        database_parameters,
        log
    )
    peptides, peptide_index_matrix, digestion_matrix = src.peptides.digestProteins(
        proteins,
        total_protein_sequence,
        ptm_matrix,
        database_parameters,
        log,
    )
    peptide_masses, fragments = src.peptides.calculateMasses(
        peptides,
        peptide_index_matrix,
        base_mass_dict,
        total_protein_sequence,
        database_parameters,
        log,
    )
    h5_dict = {
        "parameters": (database_parameters, "json"),
        "base_mass_dict": (base_mass_dict, "json"),
        "proteins": (proteins, "npy"),
        "total_protein_sequence": (total_protein_sequence, "str"),
        "ptms": (ptms, "npy"),
        "ptm_matrix": (ptm_matrix, "csr"),
        "peptides": (peptides, "npy"),
        "peptide_index_matrix": (peptide_index_matrix, "csr"),
        "digestion_matrix": (digestion_matrix, "csr"),
        "peptide_masses": (peptide_masses, "npy"),
        "fragments": (fragments, "npy"),
    }
    src.io.saveH5Dict(file_name, h5_dict)


def loadDatabase(file_name):
    h5_dict = {
        "parameters": "json",
        "base_mass_dict": "json",
        "proteins": "npy",
        "total_protein_sequence": "str",
        "ptms": "npy",
        "ptm_matrix": "csr",
        "peptides": "npy",
        "peptide_index_matrix": "csr",
        "digestion_matrix": "csr",
        "peptide_masses": "npy",
        "fragments": "npy",
    }
    result = src.io.loadH5Dict(file_name, h5_dict)
    # parameters = result["parameters"]
    # base_mass_dict = result["base_mass_dict"]
    # proteins = result["proteins"]
    # total_protein_sequence = result["total_protein_sequence"]
    # ptms = result["ptms"]
    # ptm_matrix = result["ptm_matrix"]
    # peptides = result["peptides"]
    # peptide_index_matrix = result["peptide_index_matrix"]
    # digestion_matrix = result["digestion_matrix"]
    # peptide_masses = result["peptide_masses"]
    # fragments = result["fragments"]
    return result


if __name__ == '__main__':
    pass
