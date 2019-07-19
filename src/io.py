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
import h5py
import matplotlib
matplotlib.use("Agg", warn=False)
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


def loadJSON(file_name, parameters=None, log=None):
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


def saveJSON(json_dict, file_name, parameters=None, log=None):
    ''' Saves a JSON dict to the corresponding file name in parameters'''
    if parameters is None:
        parameter_file_name = file_name
    else:
        parameter_file_name = parameters[file_name]
    if log is not None:
        log.printIO(file_name, "Saving", parameter_file_name)
    with open(parameter_file_name, "w") as outfile:
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


def saveH5Dict(file_name, h5_dict):
    with h5py.File(file_name, "w") as f:
        dt = h5py.special_dtype(vlen=str)
        for key, (value, value_type) in h5_dict.items():
            # print(key, value_type, type(value))
            if value_type == "npy":
                try:
                    tmp = f.create_dataset(key, data=value, compression="lzf")
                except TypeError:
                    g = f.create_group(key)
                    g.attrs['type'] = "npy"
                    g.attrs['shape'] = (len(value), len(value.dtype))
                    for val_id, (val_subtype, tmp) in value.dtype.fields.items():
                        if val_subtype == np.dtype(object):
                            tmp = g.create_dataset(val_id, data=value[val_id], dtype=dt)
                        else:
                            tmp = g.create_dataset(val_id, data=value[val_id], compression="lzf")
            elif value_type == "csr":
                g = f.create_group(key)
                tmp = g.create_dataset('data', data=value.data, compression="lzf")
                tmp = g.create_dataset('indptr', data=value.indptr, compression="lzf")
                tmp = g.create_dataset('indices', data=value.indices, compression="lzf")
                g.attrs['shape'] = value.shape
                g.attrs['type'] = "csr"
            elif value_type == "json":
                tmp = f.create_dataset(key, data=json.dumps(value))
            elif value_type == "str":
                tmp = f.create_dataset(key, data=value)


def loadH5Dict(file_name, h5_dict):
    with h5py.File(file_name, "r") as f:
        r = {}
        for key, key_type in h5_dict.items():
            if key_type == "npy":
                if isinstance(f[key], h5py.Dataset):
                    r[key] = f[key][...]
                else:
                    g = f[key]
                    shape = g.attrs['shape']
                    elements = {k: g[k][...] for k in list(g)}
                    d = []
                    for i, e in elements.items():
                        d += [(i, e.dtype.descr[0][1])]
                    a = np.empty(shape[0], dtype=np.dtype(d))
                    for i, e in elements.items():
                        a[i] = e
                    r[key] = a
            elif key_type == "csr":
                g = f[key]
                r[key] = scipy.sparse.csr_matrix(
                    (
                        g['data'][:],
                        g['indices'][:],
                        g['indptr'][:]
                    ),
                    shape=g.attrs['shape']
                )
            elif key_type == "json":
                r[key] = json.loads(str(f[key][...]))
            elif key_type == "str":
                r[key] = str(f[key][...])
    return r


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
        self.original_stdout.flush()
        self.original_stderr.flush()

    def printMessage(self, text=None):
        if text is not None:
            print(
                "{} > {}{}".format(
                    time.asctime(),
                    " " * (self.indent_size * self.current_indent),
                    text
                )
            )
            self.flush()

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


def plotCalibrationResults(estimation_aggregates, parameters, log):
    # Plot calibration
    with log.newSection("Plotting calibration results"):
        # estimation_aggregates = loadArray(
        #     "PSEUDO_AGGREGATE_IONS_FILE_NAME",
        #     parameters,
        #     log
        # )[1::2]
        fig, ax = plt.subplots(2, 3, sharex="col", sharey=True)
        tmp = plt.subplots_adjust(hspace=0.1, wspace=0.1)
        quick_aggregates = np.empty(
            len(estimation_aggregates),
            dtype=estimation_aggregates.dtype
        )
        for attribute in quick_aggregates.dtype.names:
            quick_aggregates[attribute] = np.average(estimation_aggregates[attribute], axis=1)
        for attribute, x_loc, y_loc in [
            ("MZ", 0, 0),
            ("CALIBRATED_MZ", 0, 1),
            ("DT", 1, 0),
            ("CALIBRATED_DT", 1, 1),
            ("RT", 2, 0),
            ("CALIBRATED_RT", 2, 1),
        ]:
            average_attribute = quick_aggregates[attribute].reshape((-1, 1))
            data = estimation_aggregates[attribute] - average_attribute
            if attribute in parameters["RELATIVE_ATTRIBUTES"]:
                data *= 1000000 / average_attribute
                label = "{} PPM".format(attribute)
            else:
                label = attribute
            tmp = sns.violinplot(
                x='value',
                y='variable',
                data=pd.melt(
                    pd.DataFrame(
                        data
                    ),
                ),
                ax=ax[y_loc, x_loc],
                inner=None,
                orient="h"
            )
            tmp = ax[y_loc, x_loc].axes.get_xaxis().set_visible(False)
            tmp = ax[y_loc, x_loc].axes.get_yaxis().set_visible(False)
            tmp = ax[y_loc, x_loc].set_xlabel("")
            tmp = ax[y_loc, x_loc].set_ylabel("")
        tmp = ax[0, 0].set_title(u"Δ MZ (ppm)")
        tmp = ax[0, 1].set_title(u"Δ DT")
        tmp = ax[0, 2].set_title(u"Δ RT")
        tmp = ax[1, 0].axes.get_xaxis().set_visible(True)
        tmp = ax[1, 1].axes.get_xaxis().set_visible(True)
        tmp = ax[1, 2].axes.get_xaxis().set_visible(True)
        tmp = ax[0, 2].yaxis.set_label_position("right")
        tmp = ax[1, 2].yaxis.set_label_position("right")
        tmp = ax[0, 2].set_ylabel("Uncalibrated")
        tmp = ax[1, 2].set_ylabel("Calibrated")
        tmp = ax[0, 2].axes.get_yaxis().set_visible(True)
        tmp = ax[1, 2].axes.get_yaxis().set_visible(True)
        tmp = ax[0, 2].tick_params(axis='y', which='both', left=False)
        tmp = ax[1, 2].tick_params(axis='y', which='both', left=False)
        tmp = ax[0, 0].get_yaxis().set_ticklabels(
            [
                "Sample {}".format(x.split("_")[-2]) for x in parameters["APEX_FILE_NAMES"]
            ]
        )
        tmp = ax[1, 0].get_yaxis().set_ticklabels(
            [
                "Sample {}".format(x.split("_")[-2]) for x in parameters["APEX_FILE_NAMES"]
            ]
        )
        tmp = ax[0, 0].axes.get_yaxis().set_visible(True)
        tmp = ax[1, 0].axes.get_yaxis().set_visible(True)
        # plt.show()
        tmp = plt.savefig(parameters["PLOTS_PATH"] + "_calibration.pdf", bbox_inches='tight')
        tmp = plt.close()



if __name__ == '__main__':
    pass
