#!venv/bin/python


import numpy as np
import src.io
import src.parallelization as mp
import scipy.sparse
from collections import defaultdict


def importAllFromCsv(parameters, log, save=False):
    ''' Returns structured numpy array containing
        all ions sorted by sample and rt'''
    with log.newSection(
        "Importing all ions from {} samples".format(
            parameters["SAMPLE_COUNT"]
        )
    ):
        column_names, column_indices = list(
            zip(*parameters["APEX_COLUMNS"].items())
        )
        process_count = min(parameters["CPU_COUNT"], parameters["SAMPLE_COUNT"])
        in_queue = mp.partitionedQueue(parameters["APEX_FILE_NAMES"], process_count)
        ions = sorted(
            mp.parallelizedGenerator(
                function=__multiprocessedImportIonsFromApexCsv,
                function_args={
                    "parameters": parameters,
                    "column_indices": column_indices,
                    "rt_column": column_names.index("RT"),
                    'le_column': column_names.index("FUNCTION"),
                    "in_queue": in_queue,
                    "log": log
                },
                process_count=process_count,
            )
        )
        ions = __mergeIonArraysIntoStructuredArray(ions, column_names, log)
        if save:
            src.io.saveArray(ions, "RAW_IONS_FILE_NAME", parameters, log)
    return ions


def __multiprocessedImportIonsFromApexCsv(kwargs):
    in_queue = kwargs['in_queue']
    out_queue = kwargs['out_queue']
    parameters = kwargs['parameters']
    column_indices = kwargs['column_indices']
    rt_column = kwargs['rt_column']
    le_column = kwargs['le_column']
    log = kwargs['log']
    sample_indices = in_queue.get()
    while sample_indices is not None:
        for sample_index in sample_indices:
            apex_file_name = parameters["APEX_FILE_NAMES"][sample_index]
            ions = src.io.loadArrayFromCsv(
                apex_file_name, parameters,
                use_cols=column_indices,
                log=None
            )
            # ions = ions[:, column_indices]
            # if parameters["HE_ONLY"]:
            #     le_ions = ions[le_column] == 2
            #     if not np.all(le_ions):
            #         ions[le_ions] = False
            #     else:
            #         ions = ions[~ions["LE"]]
            if parameters["HE_ONLY"]:
                le_ions = ions[:, le_column] == 1
                if not np.all(le_ions):
                    ions = ions[~le_ions]
                ions[:, le_column] = 0
            ions = ions[np.argsort(ions[:, rt_column])]
            sample_column = np.full(len(ions), sample_index).reshape((-1, 1))
            ions = np.concatenate([ions, sample_column], axis=1)
            out_queue.put((sample_index, ions))
            log.printMessage(
                "Imported {} ions for sample {} from {}".format(
                    len(ions),
                    sample_index,
                    apex_file_name
                )
            )
        sample_indices = in_queue.get()
    in_queue.put(None)
    out_queue.put(None)


def __mergeIonArraysIntoStructuredArray(ions, column_names, log):
    log.printMessage(
        "Merging {} ions".format(
            sum([len(i) for s, i in ions]),
        )
    )
    ions = np.concatenate(list(zip(*ions))[1], axis=0)
    structured_ions = np.empty(
        len(ions),
        dtype=[
            ("MZ", np.float),
            ("RT", np.float),
            ("DT", np.float),
            ("MZ_ERROR", np.float),
            ("RT_ERROR", np.float),
            ("DT_ERROR", np.float),
            ("INTENSITY", np.float),
            ("CALIBRATED_MZ", np.float),
            ("CALIBRATED_RT", np.float),
            ("CALIBRATED_DT", np.float),
            ("CALIBRATED_INTENSITY", np.float),
            ("SHIFTED_DT", np.float),
            ("LE", np.bool),
            ("SAMPLE", np.int),
            ("AGGREGATE_INDEX", np.int),
        ]
    )
    structured_ions["MZ"] = ions[:, column_names.index("MZ")]
    structured_ions["RT"] = ions[:, column_names.index("RT")]
    structured_ions["DT"] = ions[:, column_names.index("DT")]
    structured_ions["MZ_ERROR"] = ions[:, column_names.index("MZ_ERROR")]
    structured_ions["RT_ERROR"] = ions[:, column_names.index("RT_ERROR")]
    structured_ions["DT_ERROR"] = ions[:, column_names.index("DT_ERROR")]
    structured_ions["INTENSITY"] = ions[:, column_names.index("INTENSITY")]
    structured_ions["CALIBRATED_RT"] = structured_ions["RT"]
    structured_ions["CALIBRATED_DT"] = structured_ions["DT"]
    structured_ions["CALIBRATED_MZ"] = structured_ions["MZ"]
    structured_ions["CALIBRATED_INTENSITY"] = structured_ions["INTENSITY"]
    structured_ions["SHIFTED_DT"] = structured_ions["DT"]
    structured_ions["LE"] = (ions[:, column_names.index("FUNCTION")] == 1)
    structured_ions["SAMPLE"] = ions[:, -1]
    structured_ions["AGGREGATE_INDEX"] = 0
    return structured_ions


def getPseudoAggregatesIndices(ions, parameters, log, save=True):
    ''''TODO comment'''
    if (not parameters["AUTO_ESTIMATE_ION_ALIGNMENT_PARAMETERS"]) and (not parameters["AUTO_CALIBRATE_ION_ALIGNMENT_PARAMETERS"]):
        return None
    with log.newSection("Detecting pseudo aggregate ion indices"):
        sample_count = parameters['SAMPLE_COUNT']
        intensity_threshold = parameters['QUICK_ANCHOR_INTENSITY_THRESHOLD']
        # calibration_selection = ions["LE"]
        # if parameters["HE_ONLY"]:
        #     calibration_selection = ~calibration_selection
        # le_indices = np.flatnonzero(
        #     calibration_selection & (
        #         ions["INTENSITY"] > 2**intensity_threshold
        #     )
        # )
        to_select_per_sample = parameters['QUICK_ANCHOR_TOP_ION_COUNT_PER_SAMPLE']
        start_indices = np.searchsorted(ions["SAMPLE"], np.arange(parameters["SAMPLE_COUNT"] + 1))
        le_indices = np.concatenate(
            [
                start + np.argpartition(
                    ions["INTENSITY"][start: end], - to_select_per_sample
                )[-to_select_per_sample:] for start, end in zip(
                    start_indices[:-1],
                    start_indices[1:]
                )
            ]
        )
        le_ions = ions[le_indices]
        le_order = np.argsort(le_ions["MZ"])
        le_ions = le_ions[le_order]
        original_ion_indices = le_indices[le_order]
        mz_diffs = np.diff(le_ions["MZ"])
        ppm_diffs = mz_diffs * 1000000 / le_ions["MZ"][:-1]
        large_breaks = np.ones(len(ppm_diffs) - sample_count, dtype=np.bool)
        for i in range(1, sample_count):
            large_left = ppm_diffs[:-i] > ppm_diffs[i:]
            large_right = ppm_diffs[sample_count - i:] > ppm_diffs[:-(sample_count - i)]
            large_breaks &= large_left[:-(sample_count - i)]
            large_breaks &= large_right[i:]
        large_breaks = np.flatnonzero(large_breaks) + 1
        anchor_indices = []
        for first_ion in large_breaks:
            anchor = le_ions[first_ion: first_ion + sample_count]
            samples = np.unique(anchor["SAMPLE"])
            if len(samples) != sample_count:
                continue
            sample_order = np.argsort(anchor["SAMPLE"])
            anchor_indices.append(sample_order + first_ion)
        anchor_indices = original_ion_indices[anchor_indices]
        selected = __removePseudoAggregateOutliers(
            ions[anchor_indices],
            parameters,
            log,
            filter_attributes=["RT", "DT"],
        )
        quick_indices = anchor_indices[selected]
        log.printMessage(
            "Found {} pseudo aggregate ions".format(np.sum(selected))
        )
        if save:
            src.io.saveArray(
                ions[quick_indices],
                "PSEUDO_AGGREGATE_IONS_FILE_NAME",
                parameters,
                log
            )
    return quick_indices


def __removePseudoAggregateOutliers(
    quick_anchors,
    parameters,
    log,
    filter_attributes=["RT", "DT"]
):
    # https://stackoverflow.com/questions/11882393/matplotlib-disregard-outliers-when-plotting/11886564#11886564
    selected_size = -1
    selected = np.ones(len(quick_anchors), dtype=np.bool)
    while selected_size != np.sum(selected):
        selected_size = np.sum(selected)
        for attribute in filter_attributes:
            s = np.flatnonzero(selected)
            ptps = np.ptp(quick_anchors[attribute][s], axis=1)
            points = ptps.reshape(-1, 1)
            median = np.median(points, axis=0)
            diff = np.sqrt(np.sum((points - median)**2, axis=-1))
            med_abs_deviation = np.median(diff)
            modified_z_scores = 0.6745 * diff / med_abs_deviation
            outliers = modified_z_scores > parameters['QUICK_ANCHOR_OUTLIER_THRESHOLD']
            selected[s[outliers]] = False
    log.printMessage(
        "Removed {} outliers".format(
            len(selected) - np.sum(selected)
        )
    )
    return selected


def calibrateAll(ions, pseudo_aggregate_ions, parameters, log, save=False):
    '''TODO comment'''
    if not parameters["AUTO_CALIBRATE_ION_ALIGNMENT_PARAMETERS"]:
        return ions
    with log.newSection("Calibrating all ions"):
        calibration_dict = __determineCalibrationTargets(
            pseudo_aggregate_ions,
            np.min(ions["RT"]),
            np.max(ions["RT"]),
            parameters,
            log
        )
        ions = __calibrate(
            ions,
            calibration_dict,
            parameters,
            log
        )
        if save:
            src.io.saveArray(ions, "CALIBRATED_IONS_FILE_NAME", parameters, log)
    return ions


def __determineCalibrationTargets(
    calibration_anchors,
    min_rt,
    max_rt,
    parameters,
    log
):
    calibration_dict = {}
    calibration_dict["RT"] = __determineRTTargets(
        calibration_anchors,
        min_rt,
        max_rt,
        parameters,
        log
    )
    for attribute in ["MZ", "DT"]:
        calibration_dict[attribute] = __determineAttributeTargets(
            calibration_anchors,
            parameters,
            attribute,
            log
        )
    return calibration_dict


def __determineRTTargets(calibration_anchors, min_rt, max_rt, parameters, log):
    with log.newSection("Calculating RT transformations"):
        sample_count = parameters["SAMPLE_COUNT"]
        rts = calibration_anchors["RT"]
        average_rts = np.average(rts, axis=1)
        rts = rts[np.argsort(average_rts)]
        sample_orders = np.argsort(rts, axis=0)
        rt_predictor_breaks = [-1]
        max_order = np.zeros(sample_count, dtype=np.int)
        for i, order in enumerate(sample_orders[:-1]):
            max_order = np.max([max_order, order], axis=0)
            # Rare situation in which rt[i, sample] == rt[i + 1, sample]
            if np.all(max_order <= i) and np.all(rts[i] < rts[i + 1]):
                rt_predictor_breaks.append(i)
        rt_predictor_breaks.append(len(calibration_anchors) - 1)
        rt_predictor_breaks = np.array(rt_predictor_breaks) + 1
        rt_predictor_rts = [np.full(sample_count, min_rt)]
        for i, first_break in enumerate(rt_predictor_breaks[:-1]):
            second_break = rt_predictor_breaks[i + 1]
            rt_predictor = np.average(rts[first_break: second_break], axis=0)
            rt_predictor_rts.append(rt_predictor)
        rt_predictor_rts.append(np.full(sample_count, max_rt))
        rt_predictor_rts = np.array(rt_predictor_rts)
        log.printMessage(
            "Found {} pseudo aggregate ion groups".format(len(rt_predictor_rts))
        )
    return rt_predictor_rts


def __determineAttributeTargets(
    calibration_anchors,
    parameters,
    attribute,
    log
):
    log.printMessage("Calculating {} offsets".format(attribute))
    average_attribute = np.average(
        calibration_anchors[attribute],
        axis=1
    ).reshape((-1, 1))
    distance = calibration_anchors[attribute] - average_attribute
    if attribute in parameters["RELATIVE_ATTRIBUTES"]:
        distance *= 1000000 / average_attribute
    sample_correction = np.median(distance, axis=0)
    return sample_correction


def __calibrate(
    ions,
    calibration_dict,
    parameters,
    log
):
    with log.newSection("Calibrating individual ions"):
        ions["CALIBRATED_RT"][:] = __calibrateRT(
            ions,
            calibration_dict,
            parameters,
            log
        )
        ion_samples = ions["SAMPLE"]
        for sample in range(parameters["SAMPLE_COUNT"]):
            sample_indices = np.flatnonzero(ion_samples == sample)
            for attribute in ["MZ", "DT"]:
                calibrated_attribute = "CALIBRATED_{}".format(attribute)
                attribute_correction = calibration_dict[attribute][sample]
                if attribute in parameters["RELATIVE_ATTRIBUTES"]:
                    ions[calibrated_attribute][sample_indices] *= (
                        1 - attribute_correction / 1000000
                    )
                else:
                    ions[calibrated_attribute][sample_indices] -= attribute_correction
    return ions


def __calibrateRT(
    ions,
    calibration_dict,
    parameters,
    log
):
    predictor_rts = calibration_dict["RT"]
    target_rts = np.average(predictor_rts, axis=1)
    ion_rts = ions["RT"].copy()
    ion_rt_order = np.argsort(ion_rts)
    inverse_ion_rt_order = np.argsort(ion_rt_order)
    ion_rts = ion_rts[ion_rt_order]
    ion_samples = ions["SAMPLE"][ion_rt_order]
    for sample in range(parameters['SAMPLE_COUNT']):
        sample_indices = np.flatnonzero(ion_samples == sample)
        sample_rts = ion_rts[sample_indices]
        new_sample_rts = sample_rts.copy()
        sample_predictor_rts = predictor_rts[:, sample]
        low_predictor_rt = sample_predictor_rts[0]
        low_index = 0
        low_target_rt = target_rts[0]
        for high_predictor_rt, high_target_rt in zip(
            sample_predictor_rts[1:],
            target_rts[1:]
        ):
            high_index = np.searchsorted(
                sample_rts,
                high_predictor_rt,
                "right"
            )
            slope = (
                high_target_rt - low_target_rt
            ) / (
                high_predictor_rt - low_predictor_rt
            )
            intercept = high_target_rt - slope * high_predictor_rt
            new_sample_rts[low_index: high_index] *= slope
            new_sample_rts[low_index: high_index] += intercept
            low_predictor_rt = high_predictor_rt
            low_index = high_index
            low_target_rt = high_target_rt
        ion_rts[sample_indices] = new_sample_rts
    return ion_rts[inverse_ion_rt_order]


def estimateAlignmentParameters(
    estimation_anchors,
    parameters,
    log,
    save=True
):
    '''TODO comment'''
    def summarizeAnchors(estimation_anchors, parameters, log):
        log.printMessage("Summarizing estimation anchors")
        anchors = np.empty(
            len(estimation_anchors),
            dtype=estimation_anchors.dtype
        )
        for attribute in anchors.dtype.names:
            anchors[attribute] = np.average(estimation_anchors[attribute], axis=1)
        return anchors
    if not parameters["AUTO_ESTIMATE_ION_ALIGNMENT_PARAMETERS"]:
        ion_alignment_parameters = src.io.loadJSON(
            "ION_ALIGNMENT_PARAMETERS_FILE_NAME",
            parameters,
        )
        return ion_alignment_parameters
    with log.newSection("Estimating alignment parameters"):
        summarized_anchors = summarizeAnchors(estimation_anchors, parameters, log)
        # percentile_limit = int(
        #     50 + 50 * parameters["ION_ALIGNMENT_PERCENTILE_THRESHOLD"]
        # )
        deviation_factor = parameters["ION_ALIGNMENT_DEVIATION_FACTOR"]
        ion_alignment_parameters = {}
        for attribute in [
            "CALIBRATED_MZ",
            "CALIBRATED_DT",
            "CALIBRATED_RT",
        ]:
            average_attribute = summarized_anchors[attribute].reshape((-1, 1))
            data = estimation_anchors[attribute] - average_attribute
            # ptps = np.ptp(data, axis=1)
            # if attribute in parameters["RELATIVE_ATTRIBUTES"]:
            #     ptps *= 1000000 / np.min(estimation_anchors[attribute], axis=1)
            # ptp_limit = np.percentile(ptps, percentile_limit) * deviation_factor
            if attribute in parameters["RELATIVE_ATTRIBUTES"]:
                data *= 1000000 / np.min(estimation_anchors[attribute], axis=1).reshape(-1, 1)
            # ptp_limit = np.max(np.percentile(data, percentile_limit, axis=0) * deviation_factor)
            # attribute_std = np.std(np.abs(data), axis=0)
            # ptp_limit = np.max(attribute_std * deviation_factor)
            attribute_std = np.sort(np.std(data, axis=0))
            ptp_limit = np.sqrt(attribute_std[-1]**2 + attribute_std[-2]**2) * deviation_factor
            # ptps = np.ptp(estimation_anchors[attribute], axis=1)
            # if attribute in parameters["RELATIVE_ATTRIBUTES"]:
            #     ptps *= 1000000 / np.min(estimation_anchors[attribute], axis=1)
            # ptp_limit = np.average(ptps) + np.std(ptps) * deviation_factor
            # if attribute in parameters["RELATIVE_ATTRIBUTES"]:
            #     data *= 1000000 / np.min(estimation_anchors[attribute], axis=1).reshape(-1, 1)
            # ptp_limit = np.percentile(data.flatten(), percentile_limit, axis=0) * deviation_factor
            ion_alignment_parameters[attribute] = ptp_limit
            log.printMessage(
                "Estimated maximum distance for {} is {}".format(
                    attribute,
                    ptp_limit
                )
            )
        if save:
            src.io.saveJSON(
                ion_alignment_parameters,
                "ION_ALIGNMENT_PARAMETERS_FILE_NAME",
                parameters,
                log
            )
    return ion_alignment_parameters


def sort(ions, attribute, log):
    '''Sort ions by attribute'''
    with log.newSection("Sorting ions by {}".format(attribute)):
        ions = ions[np.argsort(ions[attribute])]
    return ions


def detectAllIonNeighbors(
    ions,
    ion_alignment_parameters,
    parameters,
    log,
    save=False,
    pre_sort=True,
    trim=True,
    restitch=True
):
    with log.newSection("Calculating ion neighbors"):
        if pre_sort:
            ions = src.ions.sort(ions, "CALIBRATED_MZ", log)
        log.printMessage("Connecting ions")
        process_count = parameters["CPU_COUNT"]
        upper_mz_border = np.searchsorted(
            ions["CALIBRATED_MZ"],
            ions["CALIBRATED_MZ"] * (
                # TODO
                # 1 + ion_alignment_parameters["CALIBRATED_MZ"] / 1000000
                # 1 + parameters["ION_ALIGNMENT_DEVIATION_FACTOR"] * ion_alignment_parameters["CALIBRATED_MZ"] / 1000000
                1 + (
                    parameters["ION_ALIGNMENT_DEVIATION_FACTOR"] * (
                        np.sqrt(np.max(ions["MZ_ERROR"])**2 + ions["MZ_ERROR"]**2)
                    ) / 1000000
                )
            ),
            "right"
        )
        in_queue = mp.partitionedQueue(ions, process_count)
        neighbors = scipy.sparse.csr_matrix(
            (len(ions), len(ions)),
            dtype=np.bool
        )
        for partial_neighbors in mp.parallelizedGenerator(
            function=__multiprocessedDetectIonNeighbors,
            function_args={
                "in_queue": in_queue,
                "ions": ions,
                'upper_mz_border': upper_mz_border,
                # 'dt_error': ion_alignment_parameters["CALIBRATED_DT"],
                'rt_error': ion_alignment_parameters["CALIBRATED_RT"],
                'parameters': parameters
            },
            process_count=process_count,
        ):
            neighbors += partial_neighbors
        log.printMessage(
            "Found {} connected ion neighbors".format(neighbors.nnz)
        )
        if save:
            src.io.saveMatrix(
                neighbors,
                "ION_NEIGHBORS_FILE_NAME",
                parameters,
                log
            )
        log.printMessage("Making neighbors bi-directional")
        neighbors = neighbors + neighbors.T
        if trim:
            ions = trimNeighborsToAnchors(
                ions,
                neighbors,
                parameters,
                log
            )
            if restitch:
                ions = __reStitchAggregates(ions, neighbors, parameters, log)
        if save:
            src.io.saveArray(
                ions,
                "IONS_UNFILTERED_FILE_NAME",
                parameters,
                log
            )
    return ions, neighbors


def __multiprocessedDetectIonNeighbors(kwargs):
    in_queue = kwargs['in_queue']
    out_queue = kwargs['out_queue']
    ions = kwargs['ions']
    upper_mz_border = kwargs['upper_mz_border']
    # dt_error = kwargs['dt_error']
    rt_error = kwargs['rt_error']
    parameters = kwargs['parameters']
    selected_ions = in_queue.get()
    neighbors = scipy.sparse.dok_matrix(
        (len(ions), len(ions)),
        dtype=np.bool
    )
    for ion_index in selected_ions:
        ion = ions[ion_index]
        upper_index = upper_mz_border[ion_index]
        candidate_indices = ion_index + 1 + np.flatnonzero(
            np.abs(
                ions["CALIBRATED_RT"][ion_index + 1: upper_index] - ion["CALIBRATED_RT"]
            ) <= rt_error
        )
        # candidate_indices = candidate_indices[
        #     np.flatnonzero(
        #         np.abs(
        #             ions["CALIBRATED_DT"][candidate_indices] - ion["CALIBRATED_DT"]
        #         ) <= dt_error
        #     )
        # ]
        candidate_indices = candidate_indices[
            np.abs(
                ions["CALIBRATED_DT"][candidate_indices] - ion["CALIBRATED_DT"]
            ) <= parameters["ION_ALIGNMENT_DEVIATION_FACTOR"] * np.sqrt(
                ions["DT_ERROR"][candidate_indices]**2 + ion["DT_ERROR"]**2
            )
        ]
        candidate_indices = candidate_indices[
            (
                ions["CALIBRATED_MZ"][candidate_indices] - ion["CALIBRATED_MZ"]
            ) * 1000000 <= parameters["ION_ALIGNMENT_DEVIATION_FACTOR"] * np.sqrt(
                ions["MZ_ERROR"][candidate_indices]**2 + ion["MZ_ERROR"]**2
            ) * ion["CALIBRATED_MZ"]
        ]
        if not parameters["HE_ONLY"]:
            candidate_indices = candidate_indices[
                ions["LE"][candidate_indices] == ion["LE"]
            ]
        candidate_indices = candidate_indices[
            ions["SAMPLE"][candidate_indices] != ion["SAMPLE"]
        ]
        # candidate_ions = ions[candidate_indices]
        # rt_distances = (
        #     candidate_ions["CALIBRATED_RT"] - ion["CALIBRATED_RT"]
        # ) / (
        #     rt_error / parameters["ION_ALIGNMENT_DEVIATION_FACTOR"]
        # )
        # dt_distances = (
        #     candidate_ions["CALIBRATED_DT"] - ion["CALIBRATED_DT"]
        # ) / (
        #     dt_error / parameters["ION_ALIGNMENT_DEVIATION_FACTOR"]
        # )
        # # mz_distances = ( # TODO
        # #     (
        # #         candidate_ions["CALIBRATED_MZ"] - ion["CALIBRATED_MZ"]
        # #     ) / (
        # #         ion["CALIBRATED_MZ"] / parameters["ION_ALIGNMENT_DEVIATION_FACTOR"]
        # #     )
        # # )
        # mz_distances = 1
        # mahalanobis_distance = np.sqrt(
        #     rt_distances**2 + dt_distances**2 + mz_distances**2
        # )
        # candidate_indices = candidate_indices[
        #     mahalanobis_distance < parameters["ION_ALIGNMENT_DEVIATION_FACTOR"]
        # ]
        neighbors[ion_index, candidate_indices] = True
    out_queue.put(neighbors.tocsr())
    out_queue.put(None)


def trimNeighborsToAnchors(ions, neighbors, parameters, log, save=False):
    '''TODO comment'''
    with log.newSection("Trimming network into unambiguous ion aggregates"):
        ion_samples = ions["SAMPLE"]
        process_count = parameters["CPU_COUNT"]
        log.printMessage("Connecting components")
        anchor_count, ion_labels = scipy.sparse.csgraph.connected_components(
            neighbors,
            directed=False,
            return_labels=True
        )
        ion_order = np.argsort(ion_labels)
        ion_label_breaks = np.concatenate(
            [
                [0],
                np.flatnonzero(np.diff(ion_labels[ion_order]) > 0) + 1,
                [len(ion_labels)]
            ],
        )
        # pre_trimmed = np.diff(ion_label_breaks) == 1
        # anchor_list = ion_order[ion_label_breaks[:-1][pre_trimmed]].tolist()
        # in_queue = mp.partitionedQueue(ion_label_breaks[:-1][~pre_trimmed], process_count)
        in_queue = mp.partitionedQueue(ion_label_breaks[:-1], process_count)
        anchor_list = []
        log.printMessage("Trimming connected components")
        for partial_neighbors in mp.parallelizedGenerator(
            function=__multiprocessedTrimNeighborsToAnchors,
            function_args={
                "in_queue": in_queue,
                "samples": ion_samples,
                'ion_order': ion_order,
                'ion_label_breaks': ion_label_breaks,
                "neighbors": neighbors,
                "parameters": parameters,
            },
            process_count=process_count,
        ):
            anchor_list += partial_neighbors
        ions = __setAnchors(ions, anchor_list, parameters, log)
        if save:
            src.io.saveArray(ions, "ANCHORED_IONS_FILE_NAME", parameters, log)
    return ions


def __multiprocessedTrimNeighborsToAnchors(kwargs):
    in_queue = kwargs["in_queue"]
    out_queue = kwargs["out_queue"]
    samples = kwargs["samples"]
    ion_order = kwargs['ion_order']
    ion_label_breaks = kwargs['ion_label_breaks']
    neighbors = kwargs["neighbors"]
    parameters = kwargs["parameters"]
    ion_border_indices = in_queue.get()
    trimmed_anchors = []
    for i in ion_border_indices:
        anchor_ions = ion_order[
            ion_label_breaks[i]: ion_label_breaks[i + 1]
        ]
        anchor_samples = samples[anchor_ions]
        trimmed_anchors += __trimAnchor(
            anchor_ions,
            neighbors,
            anchor_samples,
            parameters
        )
    out_queue.put(trimmed_anchors)
    out_queue.put(None)


def __trimAnchor(
    anchor_ions,
    neighbors,
    anchor_samples,
    parameters
):
    if len(anchor_ions) <= parameters["SAMPLE_COUNT"]:
        if len(np.unique(anchor_samples)) == len(anchor_samples):
            return [anchor_ions]
    trimmed_anchors = []
    if len(anchor_ions) < 1000:
        M = neighbors[np.ix_(anchor_ions, anchor_ions)]
    else:
        M = neighbors[anchor_ions].T.tocsr()[anchor_ions]
    M = (M**2).multiply(M)
    distance_generator = (i for i in range(1, parameters["SAMPLE_COUNT"] + 1))
    u_v_dist = next(distance_generator)
    while True:
        component_count, component_labels = scipy.sparse.csgraph.connected_components(
            M,
            directed=False,
            return_labels=True
        )
        if component_count > 1:
            break
        distances = M.astype(np.int)
        for dist in range(2, u_v_dist + 1):
            old_distances = distances.astype(np.bool)
            new_distances = old_distances * M
            new_distances = new_distances - new_distances.multiply(old_distances)
            new_distances.setdiag(0)
            new_distances.eliminate_zeros()
            distances += new_distances * dist
        for u_v_dist in distance_generator:
            old_distances = distances.astype(np.bool)
            new_distances = old_distances * M
            new_distances = new_distances - new_distances.multiply(old_distances)
            new_distances.setdiag(0)
            new_distances.eliminate_zeros()
            distances += new_distances * u_v_dist
            columns, rows = new_distances.nonzero()
            problematic_indices = np.flatnonzero(
                anchor_samples[columns] == anchor_samples[rows]
            )
            if len(problematic_indices) > 0:
                for problem_index in problematic_indices:
                    u = columns[problem_index]
                    v = rows[problem_index]
                    if u < v:
                        continue
                    # TODO clean up
                    detoured_distances = distances[(u, v), :].todense().A
                    summed_distance = np.sum(detoured_distances, axis=0)
                    partial_distances = defaultdict(list)
                    for i in np.flatnonzero(summed_distance == u_v_dist):
                        partial_distances[detoured_distances[0, i]].append(i)
                    shortest_path_length = []
                    sorted_partial_distance_indices = [
                        j for i, j in sorted(partial_distances.items())
                    ]
                    prev = len(sorted_partial_distance_indices[0])
                    for i in sorted_partial_distance_indices[1:]:
                        short = prev * len(i)
                        shortest_path_length.append(short)
                        prev = len(i)
                    minimum = min(shortest_path_length)
                    for i in [
                        i for i, v in enumerate(shortest_path_length) if v == minimum
                    ]:
                        for j in sorted_partial_distance_indices[i]:
                            for k in sorted_partial_distance_indices[i + 1]:
                                if M[j, k] != 0:
                                    M[j, k] = 0
                                    M[k, j] = 0
                M.eliminate_zeros()
                break
    for i in range(component_count):
        connected_component = np.flatnonzero(component_labels == i)
        trimmed_anchors += __trimAnchor(
            anchor_ions[connected_component],
            neighbors,
            anchor_samples[connected_component],
            parameters
        )
    return trimmed_anchors


def __setAnchors(ions, anchors, parameters, log):
    log.printMessage("Setting ion anchors")
    for anchor_index, anchor_ions in enumerate(anchors):
        # ions["ANCHOR"][anchor_ions] = anchor_index + 1
        ions["AGGREGATE_INDEX"][anchor_ions] = anchor_index
    return ions


def __reStitchAggregates(ions, neighbors, parameters, log):
    with log.newSection("Re-stitching over-trimmed anchors"):
        a, b = neighbors.nonzero()
        while True:
            anchors, anchor_ions, ions = src.aggregates.defineFromIons(
                ions,
                parameters,
                log,
                save=False,
                remove_noise=False,
                order_anchors=False
            )
            s = np.stack(
                [
                    ions["AGGREGATE_INDEX"][a],
                    ions["AGGREGATE_INDEX"][b]
                ]
            )
            n = scipy.sparse.csr_matrix(
                (
                    (s[0] != s[1]),
                    s
                ),
            )
            n.eliminate_zeros()
            x, y = n.nonzero()
            anchor_ions.data += 1
            left = (anchor_ions[x] > 0).astype(int)
            right = (anchor_ions[y] > 0).astype(int)
            anchor_ions.data -= 1
            overlap = left + right
            mergeable = np.flatnonzero((np.sum(overlap > 1, axis=1) == 0).A)
            if len(mergeable) == 0:
                return ions
            x = x[mergeable]
            y = y[mergeable]
            with np.errstate(divide='ignore'):
                mzd = 1000000 * np.abs(
                    anchors["MZ"][x] - anchors["MZ"][y]
                ) / np.maximum(anchors["MZ"][x], anchors["MZ"][y])
                mzd /= np.std(mzd)
                rtd = np.abs(anchors["RT"][x] - anchors["RT"][y])
                rtd /= np.std(rtd)
                dtd = np.abs(anchors["DT"][x] - anchors["DT"][y])
                dtd /= np.std(dtd)
            dists = np.sqrt(mzd**2 + rtd**2 + dtd**2)
            order = np.argsort(dists)
            x = x[order]
            y = y[order]
            dists = dists[order]
            elems, firsts = np.unique(x, return_index=True)
            match = firsts + firsts % 2 * -2 + 1
            pairs = np.flatnonzero(np.isin(firsts, match))
            if len(pairs) == 0:
                return ions
            merge_x = x[firsts[pairs]]
            merge_y = x[match[pairs]]
            o = merge_x > merge_y
            merge_x = merge_x[o]
            merge_y = merge_y[o]
            selected_ions = anchor_ions[merge_x]
            new_aggregates = np.repeat(merge_y, np.diff(selected_ions.indptr))
            ions["AGGREGATE_INDEX"][selected_ions.data] = new_aggregates


def calibrateIntensities(
    ions,
    anchor_ions,
    parameters,
    log,
    save=False
):
    '''TODO comment'''
    with log.newSection("Calibrating ion intensities"):
        selected_anchor_ions = anchor_ions[
            np.diff(anchor_ions.indptr) == parameters["SAMPLE_COUNT"]
        ].todense()
        selected_intensities = ions["INTENSITY"][selected_anchor_ions]
        average_intensities = np.sum(
            selected_intensities, axis=1
        ) / parameters["SAMPLE_COUNT"]
        delta_intensities = np.log(selected_intensities)
        delta_intensities -= np.log(average_intensities.reshape(-1, 1))
        sample_corrections = np.median(-delta_intensities, axis=0)
        ions["CALIBRATED_INTENSITY"] = ions["INTENSITY"] * np.exp(
            sample_corrections
        )[ions["SAMPLE"]]
        if save:
            src.io.saveArray(ions, "IONS_FILE_NAME", parameters, log)
    return ions


if __name__ == "__main__":
    pass
