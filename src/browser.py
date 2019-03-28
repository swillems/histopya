#!venv/bin/python


import src.parameters
import src.io
import src.ions
import src.aggregates
import src.peptides
import pandas as pd
import numpy as np
import src.parallelization as mp
import scipy.sparse
from collections import defaultdict
import json
import sys
import time
import traceback as tb
from scipy.special import binom as binom
from contextlib import contextmanager
from collections import defaultdict
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn import linear_model
import sklearn
import csv
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import collections as mc
from matplotlib.widgets import Slider
parameter_file_name = "data/test/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
log = src.io.Log(parameters["LOG_FILE_NAME"][:-4] + "_interactive.txt")
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
anchor_peptide_scores = src.io.loadMatrix(
    "ANCHOR_PEPTIDE_SCORES_FILE_NAME",
    parameters,
)
anchor_peptide_match_counts = src.io.loadMatrix(
    "ANCHOR_PEPTIDE_MATCH_COUNTS_FILE_NAME",
    parameters,
)


percolated_annotations = pd.read_csv(parameters["PERCOLATOR_TARGET_PIMS"], delimiter="\t")
percolated_fdrs = percolated_annotations.values[:, 2]
percolated_anchors, percolated_peptides = anchor_peptide_match_counts.nonzero()


def plot_bg():
    selection = anchors["ION_COUNT"] >= BG_MIN_REP_COUNT
    selection &= anchors["ION_COUNT"] <= BG_MAX_REP_COUNT
    bg_sc = ax.scatter(
        anchors["RT"][selection],
        anchors["DT"][selection],
        c="lightgrey",
        marker="."
    )
    def plot_bg_after_init():
        selection = anchors["ION_COUNT"] >= BG_MIN_REP_COUNT
        selection &= anchors["ION_COUNT"] <= BG_MAX_REP_COUNT
        dt_low, dt_high = plt.ylim()
        rt_low, rt_high = plt.xlim()
        selection &= anchors["DT"] <= dt_high
        selection &= anchors["DT"] >= dt_low
        selection &= anchors["RT"] <= rt_high
        selection &= anchors["RT"] >= rt_low
        bg_sc.set_offsets(
            np.c_[
                anchors["RT"][selection],
                anchors["DT"][selection]
            ]
        )
    return plot_bg_after_init



def plot_fg():
    significant_pims = percolated_annotations.values[
        percolated_fdrs <= FG_FDR,
        0
    ].astype(int)
    significant_anchors = percolated_anchors[significant_pims]
    significant_peptides = percolated_peptides[significant_pims]
    selection = anchors["ION_COUNT"][significant_anchors] >= FG_MIN_REP_COUNT
    selection &= anchors["ION_COUNT"][significant_anchors] <= FG_MAX_REP_COUNT
    fg_sc = ax.scatter(
        anchors["RT"][significant_anchors[selection]],
        anchors["DT"][significant_anchors[selection]],
        # c=significant_peptides[selection],
        # cmap="tab10",
        marker="."
    )
    norm = mpl.colors.Normalize(
        vmin=np.min(significant_peptides),
        vmax=np.max(significant_peptides)
    )
    m = cm.ScalarMappable(norm=norm, cmap="tab10")
    clrs = m.to_rgba(significant_peptides[selection])
    fg_sc.set_color(clrs)
    def plot_fg_after_init():
        significant_pims = percolated_annotations.values[
            percolated_fdrs <= FG_FDR,
            0
        ].astype(int)
        significant_anchors = percolated_anchors[significant_pims]
        significant_peptides = percolated_peptides[significant_pims]
        selection = anchors["ION_COUNT"][significant_anchors] >= FG_MIN_REP_COUNT
        selection &= anchors["ION_COUNT"][significant_anchors] <= FG_MAX_REP_COUNT
        dt_low, dt_high = plt.ylim()
        rt_low, rt_high = plt.xlim()
        selection &= anchors["DT"][significant_anchors] <= dt_high
        selection &= anchors["DT"][significant_anchors] >= dt_low
        selection &= anchors["RT"][significant_anchors] <= rt_high
        selection &= anchors["RT"][significant_anchors] >= rt_low
        fg_sc.set_offsets(
            np.c_[
                anchors["RT"][significant_anchors[selection]],
                anchors["DT"][significant_anchors[selection]]
            ]
        )
        clrs = m.to_rgba(significant_peptides[selection])
        fg_sc.set_color(clrs)
    return plot_fg_after_init


def plot_bg_nx():
    norm = mpl.colors.Normalize(
        vmin=parameters["MINIMUM_OVERLAP"][0],
        vmax=parameters["SAMPLE_COUNT"]
    )
    m = cm.ScalarMappable(norm=norm, cmap="RdYlGn")
    lc = mc.LineCollection([], [])
    bg_nx_lc = ax.add_collection(lc)
    def plot_bg_nx_after_init():
        dt_low, dt_high = plt.ylim()
        rt_low, rt_high = plt.xlim()
        selection = anchors["ION_COUNT"] >= BG_MIN_REP_COUNT
        selection &= anchors["ION_COUNT"] <= BG_MAX_REP_COUNT
        selection &= anchors["DT"] <= dt_high
        selection &= anchors["DT"] >= dt_low
        selection &= anchors["RT"] <= rt_high
        selection &= anchors["RT"] >= rt_low
        selection = np.flatnonzero(selection)
        n = neighbors[selection].T.tocsr()[selection]
        a, b = n.nonzero()
        c = a < b
        c &= n.data >= BG_NX_MIN_OVERLAP
        c &= n.data <= BG_NX_MAX_OVERLAP
        overlap = n.data[c]
        order = np.argsort(overlap)
        a = selection[a[c]][order]
        b = selection[b[c]][order]
        overlap = overlap[order]
        start_edges = list(zip(anchors["RT"][a], anchors["DT"][a]))
        end_edges = list(zip(anchors["RT"][b], anchors["DT"][b]))
        bg_nx_lc.set_segments(list(zip(start_edges, end_edges)))
        clrs = m.to_rgba(overlap)
        bg_nx_lc.set_color(clrs)
        # lc = mc.LineCollection(list(zip(start_edges, end_edges)), colors=clrs)
    return plot_bg_nx_after_init



BG_MIN_REP_COUNT = parameters["SAMPLE_COUNT"]
BG_MAX_REP_COUNT = parameters["SAMPLE_COUNT"]
FG_FDR = 0.01
FG_MIN_REP_COUNT = parameters["SAMPLE_COUNT"]
FG_MAX_REP_COUNT = parameters["SAMPLE_COUNT"]
BG_NX_MIN_OVERLAP = parameters["SAMPLE_COUNT"]
BG_NX_MAX_OVERLAP = parameters["SAMPLE_COUNT"]


plt.ion()
fig, ax = plt.subplots()
ax.axis("off")

bottom_space = 0.3
button_count = 9
# plt.subplots_adjust(bottom=bottom_space)

axcolor = 'lightgoldenrodyellow'
BG_MIN_REP_COUNT_AX = plt.axes([0.25, bottom_space / button_count, 0.65, bottom_space / (button_count + 1)], facecolor=axcolor)
BG_MAX_REP_COUNT_AX = plt.axes([0.25, 2 * bottom_space / button_count, 0.65, bottom_space / (button_count + 1)], facecolor=axcolor)
FG_MIN_REP_COUNT_AX = plt.axes([0.25, 3 * bottom_space / button_count, 0.65, bottom_space / (button_count + 1)], facecolor=axcolor)
FG_MAX_REP_COUNT_AX = plt.axes([0.25, 4 * bottom_space / button_count, 0.65, bottom_space / (button_count + 1)], facecolor=axcolor)
BG_NX_MIN_OVERLAP_AX = plt.axes([0.25, 5 * bottom_space / button_count, 0.65, bottom_space / (button_count + 1)], facecolor=axcolor)
BG_NX_MAX_OVERLAP_AX = plt.axes([0.25, 6 * bottom_space / button_count, 0.65, bottom_space / (button_count + 1)], facecolor=axcolor)
FG_FDR_AX = plt.axes([0.25, 7 * bottom_space / button_count, 0.65, bottom_space / (button_count + 1)], facecolor=axcolor)
ax = plt.axes([0.25, bottom_space, 0.65, 0.95 - bottom_space])
plt_bg_nx = plot_bg_nx()
plt_bg = plot_bg()
plt_fg = plot_fg()

BG_MIN_REP_COUNT_SLIDER = Slider(BG_MIN_REP_COUNT_AX, "BG_MIN_REP_COUNT", 1, parameters["SAMPLE_COUNT"], BG_MIN_REP_COUNT, valstep=1)
BG_MAX_REP_COUNT_SLIDER = Slider(BG_MAX_REP_COUNT_AX, "BG_MAX_REP_COUNT", 1, parameters["SAMPLE_COUNT"], BG_MAX_REP_COUNT, valstep=1)
FG_MIN_REP_COUNT_SLIDER = Slider(FG_MIN_REP_COUNT_AX, "FG_MIN_REP_COUNT", 1, parameters["SAMPLE_COUNT"], FG_MIN_REP_COUNT, valstep=1)
FG_MAX_REP_COUNT_SLIDER = Slider(FG_MAX_REP_COUNT_AX, "FG_MAX_REP_COUNT", 1, parameters["SAMPLE_COUNT"], FG_MAX_REP_COUNT, valstep=1)
BG_NX_MIN_OVERLAP_SLIDER = Slider(BG_NX_MIN_OVERLAP_AX, "BG_NX_MIN_OVERLAP", 1, parameters["SAMPLE_COUNT"], BG_NX_MIN_OVERLAP, valstep=1)
BG_NX_MAX_OVERLAP_SLIDER = Slider(BG_NX_MAX_OVERLAP_AX, "BG_NX_MAX_OVERLAP", 1, parameters["SAMPLE_COUNT"], BG_NX_MAX_OVERLAP, valstep=1)
FG_FDR_SLIDER = Slider(FG_FDR_AX, "FG_FDR", 0, 1, valinit=FG_FDR, valstep=0.001)

def slider_update(val):
    global BG_MIN_REP_COUNT
    global BG_MAX_REP_COUNT
    global FG_MIN_REP_COUNT
    global FG_MAX_REP_COUNT
    global BG_NX_MIN_OVERLAP
    global BG_NX_MAX_OVERLAP
    global FG_FDR
    BG_MIN_REP_COUNT = BG_MIN_REP_COUNT_SLIDER.val
    BG_MAX_REP_COUNT = BG_MAX_REP_COUNT_SLIDER.val
    FG_MIN_REP_COUNT = FG_MIN_REP_COUNT_SLIDER.val
    FG_MAX_REP_COUNT = FG_MAX_REP_COUNT_SLIDER.val
    BG_NX_MIN_OVERLAP = BG_NX_MIN_OVERLAP_SLIDER.val
    BG_NX_MAX_OVERLAP = BG_NX_MAX_OVERLAP_SLIDER.val
    FG_FDR = FG_FDR_SLIDER.val
    plt_bg()
    plt_fg()
    # plot_bg_nx()

BG_MIN_REP_COUNT_SLIDER.on_changed(slider_update)
BG_MAX_REP_COUNT_SLIDER.on_changed(slider_update)
FG_MIN_REP_COUNT_SLIDER.on_changed(slider_update)
FG_MAX_REP_COUNT_SLIDER.on_changed(slider_update)
BG_NX_MIN_OVERLAP_SLIDER.on_changed(slider_update)
BG_NX_MAX_OVERLAP_SLIDER.on_changed(slider_update)
FG_FDR_SLIDER.on_changed(slider_update)
