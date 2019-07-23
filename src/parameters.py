#!venv/bin/python


import os
import re
import src.io
import multiprocessing as mp
import scipy.stats
from scipy.special import binom as binom


DEFAULT_PARAMETERS_FILE_NAME = "lib/defaults/default_parameters_hdf5.json"
AUTO_COMPLETE_PATTERN = "(\[.*\])"


class ParameterError(Exception):

    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors


def parseFromCommandLine():
    import argparse
    parser = argparse.ArgumentParser(
        description="Run HistoPyA"
    )
    parser.add_argument(
        "-p",
        "--parameter_file_name",
        help="Parameter file (JSON format)",
        required=True,
    )
    parser.add_argument(
        "-a",
        "--action",
        help="(C)reate, (A)nnotate, (F)ull, (B)rowse or (G)UI (default=F)",
        choices=["C", "A", "F", "B", "G"],
        default="F",
    )
    parameter_file_name = parser.parse_args().parameter_file_name
    parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
    action = parser.parse_args().action
    return parameters, action


def getDefaultParameters():
    return src.io.loadJSON(
        DEFAULT_PARAMETERS_FILE_NAME,
        parameters=None
    )


def updateParameters(parameters, parameter_file_name=None):
    try:
        if parameter_file_name is not None:
            user_defined_parameters = src.io.loadJSON(
                parameter_file_name,
                parameters=None
            )
            parameters.update(user_defined_parameters)
        __autoCompleteAllParameters(parameters)
        __parseMultiprocessingParameters(parameters)
        __parseSamples(parameters)
        __updateDependancies(parameters)
        __setVersion(parameters)
    except Exception as e:
        raise ParameterError("No proper parameters found for input file", e)
    return parameters


def importParameterDictFromJSON(parameter_file_name, save=True):
    print("Reading parameters")
    parameters = getDefaultParameters()
    parameters = updateParameters(parameters, parameter_file_name)
    __createPaths(parameters)
    if save:
        src.io.saveJSON(
            parameters,
            "FULL_PARAMETERS_FILE_NAME",
            parameters
        )
    # except Exception:
    #     print(traceback.format_exc())
    #     print("ERROR: Parameters could not be imported.")
    #     print("       HistoPyA was aborted!")
    #     raise
    #     # sys.exit(0)
    return parameters


# def __userDefinedUpdate(key, value, parameters):
#     if isinstance(value, dict):
#         for sub_key, sub_value in value.items():
#             parameters = __userDefinedUpdate(key + [sub_key], sub_value, parameters)
#     else:
#         if len(key) == 1:
#             parameters[key[0]] = value
#         if len(key) == 2:
#             parameters[key[0]][key[1]] = value
#         if len(key) == 3:
#             parameters[key[0]][key[1]][key[2]] = value
#     return parameters


def __autoCompleteAllParameters(parameters):
    for parameter, value in parameters.items():
        if isinstance(value, str):
            __autoCompleteSingleParameter(parameters, parameter, value)
        elif isinstance(value, dict):
            for child_parameter, child_value in value.items():
                if isinstance(child_value, str):
                    __autoCompleteSingleParameter(
                        parameters,
                        child_parameter,
                        child_value,
                        parent=parameter
                    )
        elif isinstance(value, list):
            for child_parameter, child_value in enumerate(value):
                if isinstance(child_value, str):
                    __autoCompleteSingleParameter(
                        parameters,
                        child_parameter,
                        child_value,
                        parent=parameter
                    )


def __autoCompleteSingleParameter(parameters, parameter, value, parent=None):
    match = re.search(AUTO_COMPLETE_PATTERN, value)
    if match is None:
        return
    new_parameter = match.group()[1:-1]
    __autoCompleteSingleParameter(
        parameters,
        new_parameter,
        parameters[new_parameter]
    )
    if parent is None:
        parameters[parameter] = re.sub(
            AUTO_COMPLETE_PATTERN,
            parameters[new_parameter],
            value
        )
    else:
        parameters[parent][parameter] = re.sub(
            AUTO_COMPLETE_PATTERN,
            parameters[new_parameter],
            value
        )


def __createPaths(parameters):
    for parameter, value in parameters.items():
        if parameter.endswith("_PATH") and (value is not None):
            if not value.endswith("/"):
                value = value + "/"
                parameters[parameter] = value
            if not os.path.exists(value):
                os.makedirs(value)


def __parseMultiprocessingParameters(parameters):
    max_cpu_count = mp.cpu_count()
    if parameters["CPU_COUNT"] > max_cpu_count:
        parameters["CPU_COUNT"] = max_cpu_count
    elif parameters["CPU_COUNT"] <= 0:
        if parameters["CPU_COUNT"] + max_cpu_count <= 0:
            parameters["CPU_COUNT"] = 1
        else:
            parameters["CPU_COUNT"] += max_cpu_count


def __parseSamples(parameters):
    if len(parameters["APEX_FILE_NAMES"]) == 0:
        apex_file_names = []
        for file_name in sorted(os.listdir(parameters["APEX_PATH"])):
            if file_name.endswith(".csv"):
                full_file_name = os.path.join(
                    parameters["APEX_PATH"],
                    file_name
                )
                apex_file_names.append(full_file_name)
        parameters["APEX_FILE_NAMES"] = apex_file_names
    if parameters["SAMPLE_COUNT"] <= 0:
        parameters["SAMPLE_COUNT"] = len(parameters["APEX_FILE_NAMES"])
    if isinstance(parameters["MINIMUM_OVERLAP"], int):
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
                select -= 1
                new_target = binom(total, select)
                new_target *= percentile_limit**select
                new_target *= (1 - percentile_limit)**(total - select)
                target += new_target
            minimum_hits.append(max(select, minimum_overlap))
        parameters["MINIMUM_OVERLAP"] = minimum_hits


def __updateDependancies(parameters):
    if parameters["HE_ONLY"]:
        parameters["FILTER_PRECURSOR_EXISTENCE_IN_LE"] = False
        parameters['NEIGHBOR_ALL_CHANNELS'] = False
    min_overlap = min(parameters["MINIMUM_OVERLAP"])
    if parameters["SIGNAL_COUNT_THRESHOLD"] < min_overlap:
        parameters["SIGNAL_COUNT_THRESHOLD"] = min_overlap
    if parameters["USE_RT_IN_PERCOLATOR"]:
        parameters["USE_PERCOLATOR"] = True
    return parameters


def __setVersion(parameters):
    current_version = ""
    for file_name in sorted(os.listdir("docs")):
        if file_name.startswith("version_"):
            current_version = file_name
            break
    if parameters["VERSION"] not in [
        "docs/version_x.x.x",
        current_version
    ]:
        print("WARNING: versions are different!")
    parameters["VERSION"] = current_version


def importParameters(input_file_name):
    try:
        parameters = getDefaultParameters()
        parameters = src.parameters.updateParameters(
            parameters,
            input_file_name
        )
        __createPaths(parameters)
        src.io.saveJSON(
            parameters,
            "FULL_PARAMETERS_FILE_NAME",
            parameters
        )
    except ParameterError as e:
        raise e
    return parameters


if __name__ == '__main__':
    pass
