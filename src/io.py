#!venv/bin/python


import json
import sys
import time
import scipy.sparse
import numpy as np
import traceback as tb
import pandas as pd
import csv
import os
from contextlib import contextmanager
import matplotlib
matplotlib.use("Agg")
import seaborn as sns
from matplotlib import pyplot as plt


def loadArray(file_name, parameters, log=None):
    ''' Loads a numpy array from the corresponding file name in parameters'''
    if log is not None:
        log.printIO(file_name, "Loading")
    array = np.load(parameters[file_name])
    return array


def saveArray(array, file_name, parameters, log=None):
    ''' Saves a numpy array to the corresponding file name in parameters'''
    if log is not None:
        log.printIO(file_name, "Saving", parameters[file_name])
    np.save(parameters[file_name], array)


def loadMatrix(file_name, parameters, log=None):
    ''' Loads a scipy sparse matrix from
        the corresponding file name in parameters'''
    if log is not None:
        log.printIO(file_name, "Loading")
    matrix = scipy.sparse.load_npz(parameters[file_name])
    return matrix


def saveMatrix(matrix, file_name, parameters, log=None):
    ''' Saves a scipy sparse matrix to
        the corresponding file name in parameters'''
    if log is not None:
        log.printIO(file_name, "Saving", parameters[file_name])
    scipy.sparse.save_npz(parameters[file_name], matrix)


def loadJSON(file_name, parameters, log=None):
    ''' Loads a JSON dict from the corresponding file name in parameters'''
    if log is not None:
        log.printIO(file_name, "Loading")
    if parameters is None:
        with open(file_name, "r") as infile:
            json_dict = json.load(infile)
    else:
        if isinstance(parameters[file_name], list):
            json_dict = {}
            for file_name_in_list in parameters[file_name]:
                with open(file_name_in_list, "r") as infile:
                    json_dict.update(json.load(infile))
        else:
            with open(parameters[file_name], "r") as infile:
                json_dict = json.load(infile)
    return json_dict


def saveJSON(json_dict, file_name, parameters, log=None):
    ''' Saves a JSON dict to the corresponding file name in parameters'''
    if log is not None:
        log.printIO(file_name, "Saving", parameters[file_name])
    with open(parameters[file_name], "w") as outfile:
        json.dump(json_dict, outfile, indent=4, sort_keys=True)


def loadArrayFromCsv(file_name, parameters, log=None, use_cols=None):
    ''' Loads a numpy array from the corresponding CSV file name in parameters'''
    if log is not None:
        log.printIO(file_name, "Loading")
    if use_cols is None:
        array = pd.read_csv(file_name, sep=parameters["APEX_DELIMITER"]).values
    else:
        array = pd.read_csv(
            file_name,
            sep=parameters["APEX_DELIMITER"],
            usecols=use_cols
        ).values[:, np.argsort(np.argsort(use_cols))]
    return array


def saveListOfListsToCsv(
    list_of_lists,
    file_name,
    parameters,
    log=None,
    header=None,
    delimiter=","
):
    ''' Saves a numpy array to the corresponding CSV file name in parameters'''
    if log is not None:
        log.printIO(file_name, "Writing")
    with open(parameters[file_name], "w") as outfile:
        csv_file = csv.writer(outfile, delimiter=delimiter)
        if header is not None:
            csv_file.writerow(header)
            for row in list_of_lists:
                csv_file.writerow(row)


def runPercolator(parameters, log=None):
    ''' Running Percolator'''
    if not parameters["USE_PERCOLATOR"]:
        return
    log.printMessage("Running percolator")
    percolator_params = [
        parameters["PERCOLATOR_LOCATION"],
        parameters["PERCOLATOR_DATA_FILE_NAME"],
        "--results-psms", parameters["PERCOLATOR_TARGET_PIMS"],
        "--results-peptides", parameters["PERCOLATOR_TARGET_PEPTIDES"],
        "--decoy-results-psms", parameters["PERCOLATOR_DECOY_PIMS"],
        "--results-proteins", parameters["PERCOLATOR_TARGET_PROTEINS"],
        "-I", "concatenated",
        "-A",
        # "-V", "score", # Can become negative if mod_score is also included!
        "-P",
        "DECOY_",
        "-x"
    ]
    if parameters["USE_RT_IN_PERCOLATOR"]:
        percolator_params += [
            "-D",
            "15",
        ]
    percolator_params += [
        "2>&1",
        "|",
        "tee",
        parameters["PERCOLATOR_LOG"]
    ]
    os.system(" ".join(percolator_params))


class Log(object):
    '''TODO comment full class'''

    def __init__(self, file_name, indent_size=4):
        self.original_stdout = sys.stdout
        self.original_stderr = sys.stderr
        self.current_indent = 0
        self.indent_size = indent_size
        self.file_name = file_name
        file_pre_exists = os.path.isfile(file_name)
        self.file = open(self.file_name, 'a')
        sys.stdout = self
        sys.stderr = self
        if file_pre_exists:
            print("#############################")
            print("# APPENDING TO EXISTING LOG #")
            print("#############################")
        else:
            print("#############################")
            print("#     CREATING NEW LOG      #")
            print("#############################")

    def write(self, text):
        text = str(text)
        self.original_stdout.write(text)
        self.file.write(text)

    def close(self):
        self.file.close()

    def flush(self):
        self.file.flush()

    def printMessage(self, text=None):
        if text is not None:
            print(
                "{} > {}{}".format(
                    time.asctime(),
                    " " * (self.indent_size * self.current_indent),
                    text
                )
            )

    @contextmanager
    def newSection(self, text=None, newline=False):
        if text is not None:
            self.printMessage(text)
        self.current_indent += 1
        yield
        self.current_indent -= 1
        if newline:
            self.printMessage("")

    def printIO(self, parameter_name, pre_text, file_name=None):
        post_text = " ".join(parameter_name.split("_")[:-2]).lower()
        if pre_text == "Saving":
            post_text = "{} to {}".format(post_text, file_name)
        self.printMessage(
            "{} {}".format(
                pre_text,
                post_text
            )
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None:
            self.write(tb.format_exc())
            self.current_indent = 0
            self.printMessage("HistoPyA was aborted due to an error!")
        self.flush()
        self.close()
        sys.stdout = self.original_stdout
        sys.stderr = self.original_stderr

    def __str__(self):
        return self.file_name

    def __repr__(self):
        return self.file_name


if __name__ == '__main__':
    pass



































































# TODO moved from ions!!!!
def plotEstimates(
    pseudo_aggregate_ions,
    parameters,
    log
):
    plot_attributes = [
        "MZ",
        "RT",
        "DT",
        "INTENSITY",
        "CALIBRATED_MZ",
        "CALIBRATED_RT",
        "CALIBRATED_DT",
        "CALIBRATED_INTENSITY"
    ]
    plotAnchorDistributions(
        pseudo_aggregate_ions,
        plot_attributes,
        parameters,
        log
    )
    plotAnchors(pseudo_aggregate_ions, plot_attributes, parameters, log)


def plotAnchors(
    estimation_anchors,
    plot_attributes,
    parameters,
    log
):
    with log.newSection("Saving estimation aggregate plots"):
        figure_names = parameters["QUICK_ANCHOR_PLOTS"]
        plot_width = parameters["PLOT_WIDTH"]
        plot_height = parameters["PLOT_HEIGHT"]
        for attribute in plot_attributes:
            log.printMessage("Saving {} plot".format(attribute))
            for anchor in estimation_anchors:
                plt.plot(anchor[attribute], anchor["SAMPLE"])
            plt.xlabel(attribute)
            plt.ylabel("SAMPLE")
            # plt.get_current_fig_manager().resize(
            #     width=plot_width,
            #     height=plot_height
            # )
            plt.savefig(figure_names[attribute], bbox_inches='tight')
            plt.close()


def plotAnchorDistributions(
    estimation_anchors,
    plot_attributes,
    parameters,
    log
):
    def summarizeAnchors(estimation_anchors, parameters, log):
        log.printMessage("Summarizing estimation aggregates")
        anchors = np.empty(
            len(estimation_anchors),
            dtype=estimation_anchors.dtype
        )
        for attribute in anchors.dtype.names:
            anchors[attribute] = np.average(estimation_anchors[attribute], axis=1)
        return anchors
    with log.newSection("Saving estimation aggregate distribution plots"):
        anchors = summarizeAnchors(estimation_anchors, parameters, log)
        figure_names = parameters["QUICK_ANCHOR_DISTRIBUTION_PLOTS"]
        plot_width = parameters["PLOT_WIDTH"]
        plot_height = parameters["PLOT_HEIGHT"]
        for attribute in plot_attributes:
            log.printMessage("Saving {} distribution plot".format(attribute))
            average_attribute = anchors[attribute].reshape((-1, 1))
            data = estimation_anchors[attribute] - average_attribute
            if attribute in parameters["RELATIVE_ATTRIBUTES"]:
                data *= 1000000 / average_attribute
                label = "{} PPM".format(attribute)
            else:
                label = attribute
            sns.violinplot(
                x='variable',
                y='value',
                data=pd.melt(
                    pd.DataFrame(
                        data
                    ),
                )
            )
            plt.xlabel("SAMPLE")
            plt.ylabel("{} distance to aggregate averages".format(label))
            # plt.get_current_fig_manager().resize(
            #     width=plot_width,
            #     height=plot_height
            # )
            plt.savefig(figure_names[attribute], bbox_inches='tight')
            plt.close()























































#TODO moved from anchors
def plotAnchorCounts(parameters, anchors, ions, log):
    log.printMessage("Saving aggregate count plots")
    figure_names = parameters["ANCHOR_COUNT_PLOTS"]
    plot_width = parameters["PLOT_WIDTH"]
    plot_height = parameters["PLOT_HEIGHT"]
    le_ions = ions[np.flatnonzero(ions["LE"])]
    he_ions = ions[np.flatnonzero(np.logical_not(ions["LE"]))]
    for attribute, ion_set in {
        "LE": le_ions,
        "HE": he_ions
    }.items():
        for sample in range(parameters["SAMPLE_COUNT"]):
            sample_anchors = ion_set["AGGREGATE_INDEX"][ion_set["SAMPLE"] == sample]
            counts, frequencies = np.unique(
                anchors["ION_COUNT"][sample_anchors],
                return_counts=True
            )
            plt.plot(counts, frequencies)
        plt.xlabel("Aggregate size")
        plt.ylabel("Frequency")
        # plt.get_current_fig_manager().resize(
        #     width=plot_width,
        #     height=plot_height
        # )
        plt.savefig(figure_names[attribute], bbox_inches='tight')
        plt.close()

























#
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
# import seaborn as sns
# plt.show(sns.jointplot(avg_calibrated_logints[-1][::1000], cvs_calibrated[-1][::1000], kind="kde"))























#
#
#
#
#
# anchor_intensities = anchor_ions.copy()
# anchor_intensities.data = ions["CALIBRATED_INTENSITY"][anchor_intensities.data]
# anchor_intensities = np.sum(anchor_intensities, axis=1).squeeze().A.flatten()
# size_subsets = [
#     np.flatnonzero(anchors["ION_COUNT"] == i) for i in range(
#         1, parameters["SAMPLE_COUNT"] + 1
#     )
# ]
#
# plt.show(
#     [
#         plt.boxplot(
#             [
#                 np.log2(
#                     anchor_intensities[subset][
#                         anchors["LE"][subset]
#                     ]
#                 ) for subset in size_subsets
#             ]
#         )
#     ] + [
#         plt.ylabel("Log2(avg_intensity)"),
#         plt.xlabel("Aggregate ion count"),
#         plt.title("Low energy aggragates")
#     ]
# )
#
# plt.show(
#     [
#         plt.boxplot(
#             [
#                 np.log2(
#                     anchor_intensities[subset][
#                         ~anchors["LE"][subset]
#                     ]
#                 ) for subset in size_subsets
#             ]
#         )
#     ] + [
#         plt.ylabel("Log2(avg_intensity)"),
#         plt.xlabel("Aggregate ion count"),
#         plt.title("High energy aggragates")
#     ]
# )
