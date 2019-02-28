#!venv/bin/python


import scipy.sparse
import src.io
import numpy as np
import src.parallelization as mp
from sklearn import linear_model
from scipy.optimize import curve_fit


def defineFromIons(ions, parameters, log, save=True, remove_noise=True):
    '''TODO COMMENT'''
    with log.newSection("Defining aggregate ions"):
        signal_count_threshold = parameters["SIGNAL_COUNT_THRESHOLD"]
        anchor_ions = __matchIonsToAnchors(ions, signal_count_threshold, log)
        anchors = __defineAnchorProperties(anchor_ions, ions, log)
        anchors, anchor_ions, ions = __reorderAnchorsAndIons(
            anchors,
            anchor_ions,
            ions,
            signal_count_threshold,
            log
        )
        if save:
            # src.io.saveArray(anchors, "ANCHORS_FILE_NAME", parameters, log)
            src.io.saveMatrix(
                anchor_ions,
                "ANCHOR_IONS_FILE_NAME",
                parameters,
                log
            )
            # src.io.saveArray(ions, "REORDERED_IONS_FILE_NAME", parameters, log)
    return anchors, anchor_ions, ions


def __matchIonsToAnchors(ions, signal_count_threshold, log):
    log.printMessage("Matching ions to aggregate ions")
    anchor_ions = scipy.sparse.csr_matrix(
        (
            np.arange(len(ions)),
            (ions["AGGREGATE_INDEX"], ions["SAMPLE"]),
        ),
        dtype=np.int
    )
    log.printMessage("Found {} aggregate ions".format(anchor_ions.shape[0]))
    if signal_count_threshold > 1:
        log.printMessage("Skipping noisy aggregate ions")
        anchor_sizes = np.diff(anchor_ions.indptr)
        anchor_ions = anchor_ions[
            anchor_sizes >= signal_count_threshold
        ]
        log.printMessage("Retained {} aggregate ions".format(anchor_ions.shape[0]))
    return anchor_ions


def __defineAnchorProperties(anchor_ions, ions, log=None):
    log.printMessage("Calculating ion aggregate properties")
    anchors = np.empty(
        anchor_ions.shape[0],
        dtype=[
            ("MZ", np.float),
            ("RT", np.float),
            ("DT", np.float),
            ("SHIFTED_DT", np.float),
            ("LE", np.bool),
            ("ION_COUNT", np.int),
        ]
    )
    anchors["ION_COUNT"] = np.diff(anchor_ions.indptr)
    tmp_attributes = anchor_ions.copy()
    tmp_attributes.data = ions["LE"][anchor_ions.data]
    attribute_row_sums = tmp_attributes.sum(axis=1).A1
    anchors["LE"] = attribute_row_sums > 0
    for attribute in [
        "MZ",
        "RT",
        "DT",
    ]:
        tmp_attributes.data = ions[
            "CALIBRATED_{}".format(attribute)
        ][anchor_ions.data]
        attribute_row_sums = tmp_attributes.sum(axis=1).A1
        anchors[attribute] = attribute_row_sums / anchors["ION_COUNT"]
    anchors["SHIFTED_DT"] = anchors["DT"]
    return anchors


def __reorderAnchorsAndIons(
    anchors,
    anchor_ions,
    ions,
    signal_count_threshold,
    log
):
    log.printMessage("Reordering anchors and ions")
    anchor_order = np.argsort(anchors["MZ"])
    anchors = anchors[anchor_order]
    anchor_ions = anchor_ions[anchor_order]
    ions["AGGREGATE_INDEX"] = -1
    ions["AGGREGATE_INDEX"][anchor_ions.data] = np.repeat(
        np.arange(len(anchors)),
        np.diff(anchor_ions.indptr)
    )
    if signal_count_threshold > 1:
        log.printMessage("Removing noisy ions")
        ions = ions[ions["AGGREGATE_INDEX"] > -1]
        anchor_ions.data = np.argsort(np.argsort(anchor_ions.data))
        log.printMessage("Retained {} ions".format(len(ions)))
    return anchors, anchor_ions, ions


def findQuickIsotopes(
    ions,
    anchors,
    anchor_ions,
    alignment_parameters,
    parameters,
    log,
    save=True
):
    '''TODO COMMENT'''
    with log.newSection("Detecting quick isotopic pairs"):
        selected = np.flatnonzero(
            (
                anchors["ION_COUNT"] == parameters["SAMPLE_COUNT"]
            ) & ~anchors["LE"]
        )
        full_anchors = anchors[selected]
        full_anchor_ions = anchor_ions[selected].toarray()
        isotopic_delta_mass = parameters["ISOTOPIC_DELTA_MASS"]
        isotopic_pairs = []
        upper_index = 1
        lower_index = 1
        mz_error = alignment_parameters["CALIBRATED_MZ"]
        rt_error = alignment_parameters["CALIBRATED_RT"]
        dt_error = alignment_parameters["CALIBRATED_DT"]
        for i, anchor in enumerate(full_anchors):
            isotope_mass = anchor["MZ"] + isotopic_delta_mass
            candidate_upper_mass = isotope_mass * (1 + mz_error / 1000000)
            candidate_lower_mass = isotope_mass * (1 - mz_error / 1000000)
            try:
                while full_anchors["MZ"][lower_index] <= candidate_lower_mass:
                    lower_index += 1
            except IndexError:
                lower_index = len(full_anchors)
            try:
                while full_anchors["MZ"][upper_index] <= candidate_upper_mass:
                    upper_index += 1
            except IndexError:
                upper_index = len(full_anchors)
            if upper_index == lower_index:
                continue
            a_ions = ions[full_anchor_ions[i]]
            matching_ions = ions[full_anchor_ions[lower_index: upper_index]]
            matches = np.abs(matching_ions["RT"] - a_ions["RT"]) < rt_error
            matches &= np.abs(matching_ions["DT"] - a_ions["DT"]) < a_ions["DT"] * (dt_error / 1000000)
            matches &= (matching_ions["LE"] == a_ions["LE"])
            matches = np.flatnonzero(np.all(matches, axis=1))
            if len(matches) == 1:
                isotopic_pairs.append(
                    (
                        i,
                        lower_index + matches[0]
                    )
                )
        log.printMessage("Found {} quick isotopes".format(len(isotopic_pairs)))
        isotopic_pairs = selected[np.array(isotopic_pairs)]
        if save:
            src.io.saveArray(
                isotopic_pairs,
                "PSEUDO_ISOTOPIC_PAIRS_FILE_NAME",
                parameters,
                log
            )
    return isotopic_pairs


def estimateAlignmentParameters(
    ions,
    quick_isotopic_pairs,
    anchor_ions,
    parameters,
    log,
    save=True
):
    '''TODO COMMENT'''
    with log.newSection("Estimating alignment parameters"):
        percentile_limit = int(
            50 + 50 * parameters["ANCHOR_ALIGNMENT_PERCENTILE_THRESHOLD"]
        )
        deviation_factor = parameters["ANCHOR_ALIGNMENT_DEVIATION_FACTOR"]
        #     ptps = np.ptp(estimation_anchors[attribute], axis=1)
        #     if attribute in parameters["RELATIVE_ATTRIBUTES"]:
        #         ptps *= 1000000 / np.min(estimation_anchors[attribute], axis=1)
        #     ptp_limit = np.percentile(ptps, percentile_limit) * deviation_factor
        # percentile_limit = 100 * parameters["ANCHOR_ALIGNMENT_PERCENTILE_THRESHOLD"]
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
            ptp_limit = list(np.percentile(ptps, percentile_limit, axis=0) * deviation_factor)
            # ptp_limit = list(np.percentile(ptps, percentile_limit, axis=0))
            alignment_parameters[attribute] = ptp_limit
            worst_sample = np.argmax(ptp_limit)
            log.printMessage(
                "Estimated max distance for {} is {} in sample {}".format(
                    attribute,
                    ptp_limit[worst_sample],
                    worst_sample
                )
            )
        if save:
            src.io.saveJSON(
                alignment_parameters,
                "ANCHOR_ALIGNMENT_PARAMETERS_FILE_NAME",
                parameters,
                log
            )
    return alignment_parameters


def findQuickFragmentPairs(
    ions,
    anchors,
    anchor_ions,
    ion_alignment_parameters,
    anchor_alignment_parameters,
    parameters,
    log,
    save=True
):
    '''TODO COMMENT'''
    if parameters["HE_ONLY"]:
        return np.array([])
    with log.newSection("Finding unfragmented precursors in he"):
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
        log.printMessage(
            "Found {} unfragmented precursors".format(
                le_he_pairs.shape[0]
            )
        )
        # TODO update anchor_alignment_parameters with SHIFTED_DT?
        if save:
            src.io.saveArray(
                le_he_pairs,
                "UNFRAGMENTED_PRECURSOR_PAIRS_FILE_NAME",
                parameters,
                log
            )
    return le_he_pairs


def calibrateChannelDriftShift(
    ions,
    le_he_pairs,
    anchor_ions,
    anchors,
    parameters,
    log,
    save=True
):
    '''TODO COMMENT'''
    def func(data, size_param, angle_param, ratio_param, constant):
        x = data[0]
        y = data[1]
        size = (np.sqrt(x**2 + y**2))
        angle = np.arctan(y / x)
        ratio = x / y
        return size_param * size + angle_param * angle + ratio_param * ratio + constant
    if not parameters["HE_ONLY"]:
        with log.newSection("Calibrating aggregate ions channel dts"):
            alldata = np.stack([ions["DT"], ions["CALIBRATED_MZ"]])
            shift_parameters = []
            new_pred = np.zeros(len(ions))
            # parameters["DT_SHIFT_PERCENTILE"] = 0.68
            # dt_shift_percentile = int(
            #     50 + 50 * parameters["DT_SHIFT_PERCENTILE"]
            # )
            # deviation_factor = parameters["DT_SHIFT_DEVIATION_FACTOR"]
            dt_shift_percentile = int(100 * parameters["DT_SHIFT_PERCENTILE"])
            dt_shift_lower_limit, dt_shift_upper_limit = parameters["DT_SHIFT_LIMITS"]
            for sample in range(parameters["SAMPLE_COUNT"]):
                sample_ions = anchor_ions[le_he_pairs, sample].toarray()
                dt_diffs = np.diff(ions["DT"][sample_ions], axis=1).ravel()
                dts = ions["DT"][sample_ions[:, 0]]
                mzs = ions["CALIBRATED_MZ"][sample_ions[:, 0]]
                relative_percentiles = np.percentile(dt_diffs / dts, range(101))
                x = dts < dt_shift_upper_limit
                x &= dts > dt_shift_lower_limit
                # x &= dt_diffs / dts > relative_percentiles[100 - dt_shift_percentile] * deviation_factor
                # x &= dt_diffs / dts < relative_percentiles[dt_shift_percentile] * deviation_factor
                x &= dt_diffs / dts > relative_percentiles[dt_shift_percentile]
                x &= dt_diffs / dts < relative_percentiles[100 - dt_shift_percentile]
                xdata = np.stack([dts, mzs])
                ydata = 1000000 * dt_diffs / dts
                popt, pcov = curve_fit(func, xdata[:, x], ydata[x])
                shift_parameters.append(popt)
                sample_indices = ions["SAMPLE"] == sample
                sample_indices &= ions["LE"]
                pred = func(alldata[:, sample_indices], *shift_parameters[sample])
                new_pred[sample_indices] = pred
            ion_dts = ions["DT"] * (1 + new_pred / 1000000)
            ions["SHIFTED_DT"] = ion_dts
            anchor_dts = anchor_ions.copy()
            anchor_dts.data = ion_dts[anchor_dts.data]
            anchor_dts = anchor_dts.sum(axis=1).squeeze().A.squeeze() / anchors["ION_COUNT"]
            anchors["SHIFTED_DT"] = anchor_dts
    if save:
        src.io.saveArray(
            ions,
            "IONS_FILE_NAME",
            parameters,
            log
        )
        src.io.saveArray(
            anchors,
            "ANCHORS_FILE_NAME",
            parameters,
            log
        )
    return anchors, ions


def detectAllAnchorNeighbors(
    ions,
    anchors,
    anchor_ions,
    alignment_parameters,
    parameters,
    log,
    save=True
):
    '''TODO COMMENT'''
    with log.newSection("Calculating anchor neighbors"):
        sample_rt_indices = __indexIonRT(
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
        process_count = parameters["CPU_COUNT"]
        in_queue = mp.partitionedQueue(anchors, process_count)
        neighbors = scipy.sparse.csr_matrix(
            (len(anchors), len(anchors)),
            dtype=np.int
        )
        log.printMessage("Calculating neighbors")
        for partial_neighbors in mp.parallelizedGenerator(
            function=__multiprocessedDetectAnchorNeighbors,
            function_args={
                "in_queue": in_queue,
                "parameters": parameters,
                "ions": ions,
                "anchors": anchors,
                "anchor_ions": anchor_ions,
                "sample_rt_indices": sample_rt_indices,
                "sample_rts": sample_rts,
                "alignment_parameters": alignment_parameters,
            },
            process_count=process_count,
        ):
            neighbors += partial_neighbors
        log.printMessage("Found {} anchor pairs".format(neighbors.nnz))
        if save:
            src.io.saveMatrix(
                neighbors,
                "ANCHOR_NEIGHBORS_FILE_NAME",
                parameters,
                log
            )
    return neighbors


def __multiprocessedDetectAnchorNeighbors(kwargs):
    in_queue = kwargs['in_queue']
    out_queue = kwargs['out_queue']
    ions = kwargs['ions']
    alignment_parameters = kwargs["alignment_parameters"]
    anchors = kwargs['anchors']
    anchor_ions = kwargs['anchor_ions']
    sample_rt_indices = kwargs['sample_rt_indices']
    sample_rts = kwargs['sample_rts']
    parameters = kwargs["parameters"]
    neighbor_all_channels = parameters['NEIGHBOR_ALL_CHANNELS']
    dt_ppm_error = alignment_parameters["DT"]
    rt_error = alignment_parameters["RT"]
    minimum_hits = np.array(parameters["MINIMUM_OVERLAP"])
    selected_anchors = in_queue.get()
    neighbors = scipy.sparse.dok_matrix(
        (len(anchors), len(anchors)),
        dtype=np.int
    )
    for anchor_index in selected_anchors:
        if anchors["ION_COUNT"][anchor_index] == 1:
            continue
        # if not anchors["LE"][anchor_index]:
        #     continue
        coeluting_anchors = []
        for ion_index in anchor_ions[anchor_index].data:
            ion = ions[ion_index]
            ion_rt = ion["RT"]
            ion_dt = ion["SHIFTED_DT"]
            ion_sample = ion["SAMPLE"]
            low_index = np.searchsorted(
                sample_rts[ion_sample],
                ion_rt - rt_error[ion_sample],
                "left"
            )
            high_index = np.searchsorted(
                sample_rts[ion_sample],
                ion_rt + rt_error[ion_sample],
                "right"
            )
            ion_neighbors = sample_rt_indices[ion_sample][low_index: high_index]
            if not neighbor_all_channels:
                ion_neighbors = ion_neighbors[ions["LE"][ion_neighbors] == ion["LE"]]
            ion_neighbor_dts = ions["SHIFTED_DT"][ion_neighbors]
            ion_neighbors = ion_neighbors[
                (np.abs(ion_neighbor_dts - ion_dt) * 1000000) < (
                    np.minimum(
                        ion_neighbor_dts, ion_dt
                    ) * dt_ppm_error[ion_sample] * (
                        1 + (ions["LE"][ion_neighbors] != ion["LE"])
                    )
                )
            ]
            coeluting_anchors.append(ions["AGGREGATE_INDEX"][ion_neighbors])
        try:
            coeluting_anchors = np.hstack(coeluting_anchors)
        except ValueError:
            continue
        coeluting_anchor_indices, coeluting_anchor_overlap_counts = np.unique(
            coeluting_anchors[coeluting_anchors > anchor_index],
            return_counts=True
        )
        candidate_indices = np.flatnonzero(
            coeluting_anchor_overlap_counts > 1
        )
        coeluting_anchor_indices = coeluting_anchor_indices[candidate_indices]
        coeluting_anchor_overlap_counts = coeluting_anchor_overlap_counts[candidate_indices]
        reproducible_indices = np.flatnonzero(
            coeluting_anchor_overlap_counts >= minimum_hits[
                np.diff(
            # coeluting_anchor_overlap_counts == np.diff(
                    anchor_ions[
                        np.ix_(
                            coeluting_anchor_indices,
                            anchor_ions[anchor_index].indices
                        )
                    ].indptr
                )
            ]
        )
        candidate_indices = coeluting_anchor_indices[reproducible_indices]
        neighbors[
            anchor_index, candidate_indices
        ] = coeluting_anchor_overlap_counts[reproducible_indices]
    out_queue.put(neighbors.tocsr())
    out_queue.put(None)


def __indexIonRT(anchor_ions, ions, parameters, log):
    log.printMessage("Indexing ion RTs")
    sample_ions = anchor_ions.T.tocsr()
    original_indices = []
    for sample in range(parameters["SAMPLE_COUNT"]):
        selected_ions = sample_ions.data[
            sample_ions.indptr[sample]: sample_ions.indptr[sample + 1]
        ]
        sample_rt_order = np.argsort(ions["RT"][selected_ions])
        original_indices.append(selected_ions[sample_rt_order])
    return original_indices


def matchAnchorsToFragments(
    fragments,
    anchors,
    base_mass_dict,
    parameters,
    log
):
    with log.newSection("Matching aggregates to singly-charged fragments"):
        fragment_mzs = np.concatenate(
            [
                fragments["Y_MR"] + base_mass_dict["atoms"]["H+"],
                fragments["B_MR"] + base_mass_dict["atoms"]["H+"],
            ]
        )
        fragment_peptide_indices = np.concatenate(
            [fragments["PEPTIDE"]] * 2
        )
        fragment_indices = np.concatenate(
            [np.arange(len(fragments))] * 2
        )
        fragment_indices[len(fragment_indices) // 2:] *= -1
        order = np.flatnonzero(fragment_mzs > base_mass_dict["atoms"]["H+"])
        order = order[np.argsort(fragment_mzs[order])]
        fragment_mzs = fragment_mzs[order]
        fragment_peptide_indices = fragment_peptide_indices[order]
        fragment_indices = fragment_indices[order]
        order = np.argsort(anchors["MZ"])
        anchor_boundaries = np.zeros((len(anchors), 2), dtype=np.int)
        anchor_boundaries[order] = __matchMasses(
            fragment_mzs,
            anchors["MZ"][order],
            # TODO id_ppm calibration
            parameters["IDENTIFICATION_PPM"],
            log
        )
        fragment_count = np.diff(anchor_boundaries, axis=1)
        log.printMessage(
            "Found {} fragment explanations for {} anchors".format(
                np.sum(fragment_count),
                np.sum(fragment_count > 0)
            )
        )
    return anchor_boundaries, fragment_peptide_indices, fragment_indices


def __matchMasses(target_mzs, query_mzs, ppm, log):
    log.printMessage("Matching masses")
    lower_boundaries = np.searchsorted(
        target_mzs,
        query_mzs * (1 - ppm / 1000000),
        "left"
    )
    upper_boundaries = np.searchsorted(
        target_mzs,
        query_mzs * (1 + ppm / 1000000),
        "right"
    )
    anchor_boundaries = np.stack(
        [
            lower_boundaries,
            upper_boundaries
        ]
    ).T
    return anchor_boundaries


def getAnchorPeptideMatrix(
    anchors,
    neighbors,
    peptides,
    anchor_boundaries,
    fragment_peptide_indices,
    parameters,
    log,
    save=True
):
    with log.newSection("Determining best peptide explanations per HE aggregate"):
        process_count = parameters["CPU_COUNT"]
        # he_neighbors = neighbors.copy()
        # he_neighbors.data = he_neighbors.data >= parameters["MINIMUM_OVERLAP"]
        # he_neighbors.data &= ~anchors["LE"][he_neighbors.indices]
        he_neighbors = neighbors.astype(np.bool)
        he_neighbors.data = ~anchors["LE"][he_neighbors.indices]
        he_neighbors.data &= np.repeat(~anchors["LE"], np.diff(he_neighbors.indptr))
        he_neighbors.eliminate_zeros()
        anchor_indices = np.flatnonzero(
            (
                anchor_boundaries[:, 1] - anchor_boundaries[:, 0] >= 3
                # 3 are needed: 1 target and 2 for regression
            ) & (
                np.diff(he_neighbors.indptr) >= 2
            ) & (
                ~anchors["LE"]
            )
        )
        in_queue = mp.partitionedQueue(anchor_indices, process_count)
        anchor_peptide_scores = scipy.sparse.csr_matrix(
            (len(anchors), len(peptides)),
            dtype=np.float
        )
        anchor_peptide_match_counts = scipy.sparse.csr_matrix(
            (len(anchors), len(peptides)),
            dtype=np.int
        )
        for partial_anchor_peptide_scores, partial_anchor_match_counts in mp.parallelizedGenerator(
            function=__multiprocessedAnnotatePerFragment,
            function_args={
                'in_queue': in_queue,
                'neighbors': he_neighbors,
                'anchor_indices': anchor_indices,
                'anchor_boundaries': anchor_boundaries,
                'fragment_peptide_indices': fragment_peptide_indices,
                'anchor_len': len(anchors),
                'peptide_len': len(peptides)
            },
            process_count=process_count,
        ):
            anchor_peptide_scores += partial_anchor_peptide_scores
            anchor_peptide_match_counts += partial_anchor_match_counts
        anchor_peptide_scores_with_explicit_zeros = anchor_peptide_match_counts.astype(np.bool)
        anchor_peptide_scores_with_explicit_zeros += anchor_peptide_scores
        anchor_peptide_scores_with_explicit_zeros.data -= 1
        log.printMessage(
            "Found {} anchors with {} peptide explanations".format(
                np.sum(np.diff(anchor_peptide_match_counts.indptr) > 0),
                anchor_peptide_match_counts.nnz
            )
        )
        if save:
            src.io.saveMatrix(
                anchor_peptide_match_counts,
                "ANCHOR_PEPTIDE_MATCH_COUNTS_FILE_NAME",
                parameters,
                log
            )
            src.io.saveMatrix(
                anchor_peptide_scores_with_explicit_zeros,
                "ANCHOR_PEPTIDE_SCORES_FILE_NAME",
                parameters,
                log
            )
    return anchor_peptide_scores_with_explicit_zeros, anchor_peptide_match_counts


def __multiprocessedAnnotatePerFragment(kwargs):
    in_queue = kwargs['in_queue']
    out_queue = kwargs['out_queue']
    neighbors = kwargs['neighbors']
    anchor_indices = kwargs['anchor_indices']
    anchor_boundaries = kwargs['anchor_boundaries']
    fragment_peptide_indices = kwargs['fragment_peptide_indices']
    anchor_len = kwargs['anchor_len']
    peptide_len = kwargs['peptide_len']
    selected_anchor_indices = in_queue.get()
    anchor_peptide_scores = scipy.sparse.dok_matrix(
        (anchor_len, peptide_len),
        dtype=np.float
    )
    anchor_peptide_match_counts = scipy.sparse.dok_matrix(
        (anchor_len, peptide_len),
        dtype=np.int
    )
    for selected_anchor_index in selected_anchor_indices:
        anchor_index = anchor_indices[selected_anchor_index]
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
            ransac.fit(
                counts.reshape(-1, 1)[:-1],
                np.log(frequency).reshape(-1, 1)[:-1]
            )
        except ValueError:
            continue
        score = -ransac.predict(counts[-1])[0][0]
        peptide_indices = candidates[candidate_counts == counts[-1]]
        anchor_peptide_scores[anchor_index, peptide_indices] = score
        anchor_peptide_match_counts[anchor_index, peptide_indices] = counts[-1]
    out_queue.put(
        (
            anchor_peptide_scores.tocsr(),
            anchor_peptide_match_counts.tocsr(),
        )
    )
    out_queue.put(None)


def getAnchorFragmentIndices(
    anchor_peptide_match_counts,
    anchor_boundaries,
    fragment_indices,
    fragment_peptide_indices,
    parameters,
    log
):
    with log.newSection("Retrieving fragment indices"):
        selected_anchor_indices, selected_peptide_indices = anchor_peptide_match_counts.nonzero()
        anchor_fragment_indices = [None] * anchor_peptide_match_counts.nnz
        for i, (anchor_index, peptide_index) in enumerate(
            zip(selected_anchor_indices, selected_peptide_indices)
        ):
            anchor_peptides = fragment_peptide_indices[
                slice(*anchor_boundaries[anchor_index])
            ]
            anchor_fragments = fragment_indices[
                slice(*anchor_boundaries[anchor_index])
            ]
            fragments = anchor_fragments[
                anchor_peptides == peptide_index
            ]
            anchor_fragment_indices[i] = fragments
    return anchor_fragment_indices


def matchPeptidesToAnchors(
    anchors,
    peptide_masses,
    base_mass_dict,
    parameters,
    log,
):
    with log.newSection("Matching precursors to aggregates"):
        queries_are_le = parameters["FILTER_PRECURSOR_EXISTENCE_IN_LE"]
        if queries_are_le is None:
            selected_precursor_candidates = np.arange(len(anchors))
        elif queries_are_le:
            selected_precursor_candidates = np.flatnonzero(anchors["LE"])
        elif not queries_are_le:
            selected_precursor_candidates = np.flatnonzero(~anchors["LE"])
        max_precursor_charge = parameters["MAXIMUM_PRECURSOR_CHARGE"]
        anchor_mzs = np.concatenate(
            [
                (anchors["MZ"][selected_precursor_candidates] - base_mass_dict["atoms"]["H+"]) * i for i in range(
                    1,
                    max_precursor_charge + 1
                )
            ]
        )
        anchor_indices = np.concatenate(
            [selected_precursor_candidates] * max_precursor_charge
        )
        anchor_order = np.argsort(anchor_mzs)
        anchor_mzs = anchor_mzs[anchor_order]
        anchor_indices = anchor_indices[anchor_order]
        peptide_order = np.argsort(peptide_masses)
        peptide_boundaries = np.zeros((len(peptide_masses), 2), dtype=np.int)
        peptide_boundaries[peptide_order] = __matchMasses(
            anchor_mzs,
            peptide_masses[peptide_order],
            # TODO id_ppm calibration
            parameters["IDENTIFICATION_PPM"],
            log
        )
        anchor_count = np.diff(peptide_boundaries, axis=1)
        log.printMessage(
            "Found {} potential anchor explanations for {} peptides".format(
                np.sum(anchor_count),
                np.sum(anchor_count > 0)
            )
        )
    return peptide_boundaries, anchor_indices


def findFragmentPrecursors(
    anchor_peptide_match_counts,
    anchors,
    neighbors,
    anchor_alignment_parameters,
    peptide_masses,
    base_mass_dict,
    parameters,
    log
):
    with log.newSection("Matching fragments to precursors"):
        peptide_boundaries, anchor_indices = src.aggregates.matchPeptidesToAnchors(
            anchors,
            peptide_masses,
            base_mass_dict,
            parameters,
            log,
        )
        queries_are_le = parameters["FILTER_PRECURSOR_EXISTENCE_IN_LE"]
        neighbor_all_channels = parameters['NEIGHBOR_ALL_CHANNELS']
        selected_anchor_indices, selected_peptide_indices = anchor_peptide_match_counts.nonzero()
        precursor_indices = [None] * anchor_peptide_match_counts.nnz
        if queries_are_le:
            # TODO keep not neighbor_all_channels?
            if not neighbor_all_channels:
                max_rt_difference = np.max(anchor_alignment_parameters["RT"])
                max_dt_difference = np.max(anchor_alignment_parameters["DT"])
                # TODO autocalc parameters["HE_DT_SHIFT"]
                dt_center = parameters["HE_DT_SHIFT"]
                for i, (anchor_index, peptide_index) in enumerate(
                    zip(selected_anchor_indices, selected_peptide_indices)
                ):
                    anchor = anchors[anchor_index]
                    candidate_precursor_anchors = anchor_indices[
                        slice(*peptide_boundaries[peptide_index])
                    ]
                    # TODO, consistent in all runs
                    precursor_rt_difference = anchors["RT"][candidate_precursor_anchors] - anchor["RT"]
                    good_candidates = np.abs(precursor_rt_difference) < max_rt_difference
                    candidate_dt_differences = anchors["SHIFTED_DT"][candidate_precursor_anchors]
                    good_candidates &= (
                        np.abs(
                            candidate_dt_differences - anchor["SHIFTED_DT"]
                            # (candidate_dt_differences - anchor["SHIFTED_DT"]) - dt_center
                        ) / np.maximum(
                            candidate_dt_differences,
                            anchor["SHIFTED_DT"]
                        )
                    ) < max_dt_difference
                    filtered_precursors = candidate_precursor_anchors[good_candidates]
                    precursor_indices[i] = filtered_precursors
            else:
                for i, (anchor_index, peptide_index) in enumerate(
                    zip(selected_anchor_indices, selected_peptide_indices)
                ):
                    anchor = anchors[anchor_index]
                    candidate_precursor_anchors = anchor_indices[
                        slice(*peptide_boundaries[peptide_index])
                    ]
                    anchor_neighbors = neighbors.indices[
                        neighbors.indptr[anchor_index]: neighbors.indptr[anchor_index + 1]
                    ]
                    precursors = anchor_neighbors[
                        np.isin(
                            anchor_neighbors,
                            candidate_precursor_anchors
                        )
                    ]
                    precursor_indices[i] = precursors[anchors["LE"][precursors] != anchor["LE"]]
        elif not queries_are_le:
            for i, (anchor_index, peptide_index) in enumerate(
                zip(selected_anchor_indices, selected_peptide_indices)
            ):
                anchor = anchors[anchor_index]
                candidate_precursor_anchors = anchor_indices[
                    slice(*peptide_boundaries[peptide_index])
                ]
                anchor_neighbors = neighbors.indices[
                    neighbors.indptr[anchor_index]: neighbors.indptr[anchor_index + 1]
                ]
                precursors = anchor_neighbors[
                    np.isin(
                        anchor_neighbors,
                        candidate_precursor_anchors
                    )
                ]
                precursor_indices[i] = precursors[anchors["LE"][precursors] == anchor["LE"]]
        log.printMessage(
            "Found {} aggregate ions with a precursor".format(
                np.sum([1 for aa in precursor_indices if len(aa) > 0])
            )
        )
    return precursor_indices


def writePercolatorFile(
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
):
    if parameters["REQUIRE_PHYSICAL_PRECURSOR"]:
        return __writePercolatorFileWithPhysicalPrecursor(
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
            log
        )
    else:
        return __writePercolatorFileWithoutPhysicalPrecursor(
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


def __writePercolatorFileWithPhysicalPrecursor(
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
    log
):
    with log.newSection("Creating percolator data"):
        data = []
        header = [
            "PIM_id",
            "Label",
            "ScanNr",
            "rt",
            "dm",
            "ppm",
            "reproducibility_count",
            "neighbor_count",
            "match_count",
            "match_ratio",
            # "b_ion_type",
            # "y_ion_type",
            "precursor_dm",
            "precursor_ppm",
            "precursor_z",
            # "logFC",
            "peptide_length",
            "score",
            "alternatives",
            "Peptide",
            "Proteins",
        ]
        selected_anchor_indices = np.repeat(
            np.arange(anchor_peptide_match_counts.shape[0]),
            np.diff(anchor_peptide_match_counts.indptr)
        )
        for index in np.argsort(anchor_peptide_scores.data)[::-1]:
            anchor_index = selected_anchor_indices[index]
            anchor = anchors[anchor_index]
            anchor_mr = anchor["MZ"] - base_mass_dict["atoms"]["H+"]
            fragment_index = anchor_fragment_indices[index][0]  # TODO multiple candidates?
            if fragment_index > 0:
                fragment = fragments[fragment_index]
                anchor_dm = anchor_mr - fragment["Y_MR"]
            else:
                fragment = fragments[-fragment_index]
                anchor_dm = anchor_mr - fragment["B_MR"]
            anchor_ppm = anchor_dm / anchor_mr * 1000000
            anchor_neighbors = neighbors.indices[
                neighbors.indptr[anchor_index]: neighbors.indptr[anchor_index + 1]
            ]
            match_count = anchor_peptide_match_counts.data[index]
            peptide_index = anchor_peptide_match_counts.indices[index]
            peptide = peptides[peptide_index]
            if len(precursor_indices[index]) == 0:
                continue # TODO no precursor option?
            precursor_index = precursor_indices[index][0]  # TODO multiple candidates?
            precursor = anchors[precursor_index]
            precursor_z = int(round(peptide_masses[peptide_index] / precursor["MZ"]))
            precursor_mr = (
                precursor["MZ"] - base_mass_dict["atoms"]["H+"]
            ) * precursor_z
            precursor_dm = precursor_mr - peptide_masses[peptide_index]
            precursor_ppm = precursor_dm / precursor_mr * 1000000
            score = anchor_peptide_scores.data[index]
            peptide_start_index = peptide_index_matrix.indices[
                peptide_index_matrix.indptr[peptide_index]
            ]
            peptide_sequence = total_protein_sequence[
                peptide_start_index: peptide_start_index + peptide["SIZE"]
            ]
            protein_index = peptide["PROTEIN"]
            if protein_index != -1:
                protein_string = proteins[protein_index]["ID"]
            else:
                # protein_string = ";".join(
                #     proteins["ID"][
                #         peptide_protein_matrix.indices[
                #             peptide_protein_matrix.indptr[peptide_index]: peptide_protein_matrix.indptr[peptide_index + 1]
                #         ]
                #     ]
                # )
                protein_string = "Ambiguous"
            alternatives = anchor_peptide_match_counts.indptr[anchor_index + 1] - anchor_peptide_match_counts.indptr[anchor_index]
            row = [
                index,
                -1 if peptide["DECOY"] else 1,
                anchor_index,
                anchor["RT"],
                anchor_dm,
                anchor_ppm,
                anchor["ION_COUNT"],
                len(anchor_neighbors),
                match_count,
                match_count / len(anchor_neighbors),
                # "b_ion_type",
                # "y_ion_type",
                precursor_dm,
                precursor_ppm,
                precursor_z,
                # "logFC",
                len(peptide_sequence),
                score,
                alternatives,
                "-.{}.-".format(peptide_sequence),  # TODO proper flanking?
                protein_string,
            ]
            data.append(row)
        src.io.saveListOfListsToCsv(
            data,
            "PERCOLATOR_DATA_FILE_NAME",
            parameters,
            log,
            header,
            "\t"
        )
        return data


def __writePercolatorFileWithoutPhysicalPrecursor(
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
):
    with log.newSection("Creating percolator data"):
        std_errors = __estimateMZFromDT(
            anchor_peptide_scores,
            anchors,
            peptide_masses,
            peptides,
            precursor_indices,
            parameters,
            log,
        )
        data = []
        header = [
            "PIM_id",
            "Label",
            "ScanNr",
            "rt",
            "dm",
            "ppm",
            "reproducibility_count",
            "neighbor_count",
            "match_count",
            "match_ratio",
            "estimated_z",
            "sigma_mass_distance",
            "mod_score",
            "peptide_length",
            "score",
            "alternatives",
            "Peptide",
            "Proteins",
        ]
        selected_anchor_indices = np.repeat(
            np.arange(anchor_peptide_match_counts.shape[0]),
            np.diff(anchor_peptide_match_counts.indptr)
        )
        for index in np.argsort(anchor_peptide_scores.data)[::-1]:
            anchor_index = selected_anchor_indices[index]
            anchor = anchors[anchor_index]
            anchor_mr = anchor["MZ"] - base_mass_dict["atoms"]["H+"]
            fragment_index = anchor_fragment_indices[index][0]  # TODO multiple candidates?
            if fragment_index > 0:
                fragment = fragments[fragment_index]
                anchor_dm = anchor_mr - fragment["Y_MR"]
            else:
                fragment = fragments[-fragment_index]
                anchor_dm = anchor_mr - fragment["B_MR"]
            anchor_ppm = anchor_dm / anchor_mr * 1000000
            anchor_neighbors = neighbors.indices[
                neighbors.indptr[anchor_index]: neighbors.indptr[anchor_index + 1]
            ]
            match_count = anchor_peptide_match_counts.data[index]
            peptide_index = anchor_peptide_match_counts.indices[index]
            peptide = peptides[peptide_index]
            score = anchor_peptide_scores.data[index]
            peptide_start_index = peptide_index_matrix.indices[
                peptide_index_matrix.indptr[peptide_index]
            ]
            peptide_sequence = total_protein_sequence[
                peptide_start_index: peptide_start_index + peptide["SIZE"]
            ]
            protein_index = peptide["PROTEIN"]
            if protein_index != -1:
                protein_string = proteins[protein_index]["ID"]
            else:
                protein_string = "Ambiguous"
            alternatives = anchor_peptide_match_counts.indptr[anchor_index + 1] - anchor_peptide_match_counts.indptr[anchor_index]
            best_z = np.argmin(std_errors[:, index]) + 1
            mass_sigma = std_errors[best_z - 1, index]
            mod_score = score / (.5 + mass_sigma)
            # if mass_sigma > 1:
            #     continue
            row = [
                index,
                -1 if peptide["DECOY"] else 1,
                anchor_index,
                anchor["RT"],
                # mass_sigma, #anchor_dm,
                anchor_dm,
                anchor_ppm,
                anchor["ION_COUNT"],
                len(anchor_neighbors),
                match_count,
                match_count / len(anchor_neighbors),
                # abs(6.5 - z1_ratio),
                best_z,
                # anchor_dm, #mass_sigma,
                mass_sigma,
                mod_score,
                len(peptide_sequence),
                score,
                alternatives,
                "-.{}.-".format(peptide_sequence),  # TODO proper flanking?
                protein_string,
            ]
            data.append(row)
        src.io.saveListOfListsToCsv(
            data,
            "PERCOLATOR_DATA_FILE_NAME",
            parameters,
            log,
            header,
            "\t"
        )
        return data



def writeMgf(neighbors, anchors, anchor_ions, ions, parameters, log):
    with log.newSection("Writing mgf"):
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
        spectra = {
            p: n.indices[
                n.indptr[p]: n.indptr[p + 1]
            ] for p in np.flatnonzero(
                anchors["LE"] & (np.diff(n.indptr) > parameters["MINIMUM_SPECTRUM_PEAKS"])
            )
        }
        # TODO
        anchor_intensities = anchor_ions.copy()
        anchor_intensities.data = ions["CALIBRATED_INTENSITY"][anchor_intensities.data]
        anchor_intensities = anchor_intensities.sum(axis=1).A.squeeze() / anchors["ION_COUNT"]
        with open(parameters["MGF_FILE_NAME"], "w") as outfile:
            for spectrum_index, spectrum in sorted(spectra.items()):
                outfile.write("BEGIN IONS\n")
                outfile.write(
                    "TITLE=(index_{})(dt_{})(rt_{})\n".format(
                        spectrum_index,
                        anchors["DT"][spectrum_index],
                        anchors["RT"][spectrum_index],
                    )
                )
                outfile.write("SCANS={}\n".format(spectrum_index))
                outfile.write(
                    "RTINSECONDS={}\n".format(anchors["RT"][spectrum_index] * 60)
                )
                outfile.write("PEPMASS={}\n".format(anchors["MZ"][spectrum_index]))
                for anchor_index in spectrum:
                    outfile.write(
                        "{} {}\n".format(
                            anchors["MZ"][anchor_index],
                            anchor_intensities[anchor_index]
                        )
                    )
                outfile.write("END IONS\n")


def __estimateMZFromDT(
    anchor_peptide_scores,
    anchors,
    peptide_masses,
    peptides,
    precursor_indices,
    parameters,
    log,
):
    with log.newSection("Estimating precursor mass deviation with drift times"):
        peptide_list = anchor_peptide_scores.indices
        anchor_list = np.repeat(np.arange(len(anchors)), np.diff(anchor_peptide_scores.indptr))
        scores = anchor_peptide_scores.data
        dts = anchors["DT"][anchor_list]
        masses = peptide_masses[peptide_list]
        decoys = peptides["DECOY"][peptide_list]
        hits = np.array([len(i) > 0 for i in precursor_indices])
        precursor_mz = np.array(
            [
                anchors["MZ"][i[0]] if len(i) > 0 else 0 for i in precursor_indices
            ]
        )
        precursor_charges = np.zeros(len(masses), dtype=int)
        precursor_charges[hits] = np.round(
            masses[hits] / precursor_mz[hits]
        ).astype(int)
        predicted_mzs = {}
        predicted_mzs_errors = {}
        for z in np.arange(1, parameters["MAXIMUM_PRECURSOR_CHARGE"] + 1):
            x = dts[precursor_charges == z]
            y = masses[precursor_charges == z]
            d = decoys[precursor_charges == z]
            s = scores[precursor_charges == z]
            ransac = linear_model.RANSACRegressor()
            m = np.maximum(0, np.round(s[~d])).astype(int)
            ransac.fit(
                np.repeat(x[~d], m).reshape(-1, 1),
                np.repeat(y[~d], m).reshape(-1, 1)
            )
            predicted_mzs[z] = ransac.predict(dts.reshape(-1, 1)).flatten()
            y_pred = ransac.predict(x[~d].reshape(-1, 1)).flatten()
            est = y[~d] - y_pred
            predicted_mzs_errors[z] = est
        std_errors = np.stack(
            [
                np.abs(predicted_mzs[i] - masses) / np.std(predicted_mzs_errors[i]) for i in np.arange(
                    1,
                    parameters["MAXIMUM_PRECURSOR_CHARGE"] + 1
                )
            ]
        )
        # errors = np.stack(
        #     [
        #         (predicted_mzs[i] - masses) / predicted_mzs[i] for i in np.arange(
        #             1,
        #             parameters["MAXIMUM_PRECURSOR_CHARGE"] + 1
        #         )
        #     ]
        # )
        # return errors
    return std_errors


if __name__ == "__main__":
    pass
