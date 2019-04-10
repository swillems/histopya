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
from matplotlib.widgets import Slider, RadioButtons, Button
parameter_file_name = "data/tenzer/parameters.json"
#parameter_file_name = "data/ecoli_swath/parameters_manual.json"
#parameter_file_name = "data/lfq_swim_190327/parameters_manual.json"
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
proteins, total_protein_sequence, ptms, ptm_matrix = src.peptides.importProteinsAndPtms(parameters, log)
peptides, peptide_index_matrix, digestion_matrix = src.peptides.digestProteins(
    proteins,
    total_protein_sequence,
    ptm_matrix,
    parameters,
    log,
)


percolated_annotations = pd.read_csv(parameters["PERCOLATOR_TARGET_PIMS"], delimiter="\t")
percolated_fdrs = percolated_annotations.values[:, 2]
percolated_anchors, percolated_peptides = anchor_peptide_match_counts.nonzero()


def plot_bg_nx():
    norm = mpl.colors.Normalize(
        vmin=parameters["MINIMUM_OVERLAP"][0],
        vmax=parameters["SAMPLE_COUNT"]
    )
    m = cm.ScalarMappable(norm=norm, cmap="RdYlGn")
    lc = mc.LineCollection([], [])
    bg_nx_lc = ax.add_collection(lc)
    def plot_bg_nx_after_init():
        selection = np.array([], dtype=np.int)
        if NETWORK_VISIBLE == "Selected":
            if len(selected_anchors) > 0:
                selection = np.array(selected_anchors)
        if NETWORK_VISIBLE == "All":
            dt_low, dt_high = plt.ylim()
            rt_low, rt_high = plt.xlim()
            selection = anchors["ION_COUNT"] >= MIN_REP_COUNT
            selection &= anchors["ION_COUNT"] <= MAX_REP_COUNT
            selection &= anchors["DT"] <= dt_high
            selection &= anchors["DT"] >= dt_low
            selection &= anchors["RT"] <= rt_high
            selection &= anchors["RT"] >= rt_low
            selection = np.flatnonzero(selection)
        n = neighbors[selection].T.tocsr()[selection]
        a, b = n.nonzero()
        c = a > b
        a = a[c]
        b = b[c]
        overlap = n.data[c]
        order = np.argsort(overlap)
        a = selection[a[order]]
        b = selection[b[order]]
        overlap = overlap[order]
        start_edges = list(zip(anchors["RT"][a], anchors["DT"][a]))
        end_edges = list(zip(anchors["RT"][b], anchors["DT"][b]))
        bg_nx_lc.set_segments(list(zip(start_edges, end_edges)))
        clrs = m.to_rgba(overlap)
        bg_nx_lc.set_color(clrs)
        # lc = mc.LineCollection(list(zip(start_edges, end_edges)), colors=clrs)
    return plot_bg_nx_after_init


def plot_anchor_selection():
    anchor_selection = ax.scatter(
        anchors["RT"][selected_anchors],
        anchors["DT"][selected_anchors],
        facecolor="None",
        edgecolor='black',
        # s=anchors["ION_COUNT"][selected_anchors],
    )
    annotations = []
    def plot_anchor_selection_after_init():
        anchor_selection.set_offsets(
            np.c_[
                anchors["RT"][selected_anchors],
                anchors["DT"][selected_anchors]
            ]
        )
        # ax.annotate(anchors["MZ"][selected_anchors], (anchors["RT"][selected_anchors], anchors["DT"][selected_anchors]))
        for anchor_index in annotations:
            anchor_index.remove()
        del annotations[:]
        for anchor_index in selected_anchors:
            rt = anchors["RT"][anchor_index]
            dt = anchors["DT"][anchor_index]
            if LABEL_SELECTION == "m/z":
                label = round(anchors["MZ"][anchor_index], 3)
            if LABEL_SELECTION == "Peptide":
                if anchor_peptide_match_counts.indptr[anchor_index] == anchor_peptide_match_counts.indptr[anchor_index + 1]:
                    continue
                peptide_index = anchor_peptide_match_counts.indices[
                    anchor_peptide_match_counts.indptr[anchor_index]
                ]
                label = src.peptides.getSequenceFromIndex(
                    peptide_index,
                    peptides,
                    peptide_index_matrix,
                    total_protein_sequence
                )
            if LABEL_SELECTION != "None":
                ann = ax.annotate(label, (rt, dt))
                annotations.append(ann)
    return plot_anchor_selection_after_init



def plot_bg():
    selection = anchors["ION_COUNT"] >= MIN_REP_COUNT
    selection &= anchors["ION_COUNT"] <= MAX_REP_COUNT
    bg_sc = ax.scatter(
        anchors["RT"][selection],
        anchors["DT"][selection],
        c="lightgrey",
        marker=".",
        # s=anchors["ION_COUNT"][selection]
    )
    def plot_bg_after_init():
        selection = anchors["ION_COUNT"] >= MIN_REP_COUNT
        selection &= anchors["ION_COUNT"] <= MAX_REP_COUNT
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
        percolated_fdrs <= 10**LOG_FDR,
        0
    ].astype(int)
    significant_anchors = percolated_anchors[significant_pims]
    significant_peptides = percolated_peptides[significant_pims]
    selection = anchors["ION_COUNT"][significant_anchors] >= MIN_REP_COUNT
    selection &= anchors["ION_COUNT"][significant_anchors] <= MAX_REP_COUNT
    fg_sc = ax.scatter(
        anchors["RT"][significant_anchors[selection]],
        anchors["DT"][significant_anchors[selection]],
        marker=".",
        # s=anchors["ION_COUNT"][significant_anchors[selection]]
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
            percolated_fdrs <= 10**LOG_FDR,
            0
        ].astype(int)
        significant_anchors = percolated_anchors[significant_pims]
        significant_peptides = percolated_peptides[significant_pims]
        selection = anchors["ION_COUNT"][significant_anchors] >= MIN_REP_COUNT
        selection &= anchors["ION_COUNT"][significant_anchors] <= MAX_REP_COUNT
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




selected_anchors = []
LOG_FDR = -2
MIN_REP_COUNT = parameters["SAMPLE_COUNT"]
MAX_REP_COUNT = parameters["SAMPLE_COUNT"]
BG_NX_MIN_OVERLAP = parameters["MINIMUM_OVERLAP"][0]
BG_NX_MAX_OVERLAP = parameters["MINIMUM_OVERLAP"][0]
NETWORK_VISIBLE = "None"
LABEL_SELECTION = "None"


plt.ion()

bottom_space = 0.3
button_count = 5
axcolor = 'lightgoldenrodyellow'
fig, ax = plt.subplots()
ax.axis("off")
plot_width = parameters["PLOT_WIDTH"]
plot_height = parameters["PLOT_HEIGHT"]
plt.get_current_fig_manager().resize(
    width=plot_width,
    height=plot_height
)


MIN_REP_COUNT_AX = plt.axes([0.65, bottom_space / button_count, 0.3, bottom_space / (button_count + 1)], facecolor=axcolor)
MAX_REP_COUNT_AX = plt.axes([0.65, 2 * bottom_space / button_count, 0.3, bottom_space / (button_count + 1)], facecolor=axcolor)
# BG_NX_MIN_OVERLAP_AX = plt.axes([0.5, 3 * bottom_space / button_count, 0.45, bottom_space / (button_count + 1)], facecolor=axcolor)
# BG_NX_MAX_OVERLAP_AX = plt.axes([0.5, 4 * bottom_space / button_count, 0.45, bottom_space / (button_count + 1)], facecolor=axcolor)
LOG_FDR_AX = plt.axes([0.65, 3 * bottom_space / button_count, 0.3, bottom_space / (button_count + 1)], facecolor=axcolor)

MIN_REP_COUNT_SLIDER = Slider(MIN_REP_COUNT_AX, "MIN_REP_COUNT", parameters["SIGNAL_COUNT_THRESHOLD"], parameters["SAMPLE_COUNT"], MIN_REP_COUNT, valstep=1)
MAX_REP_COUNT_SLIDER = Slider(MAX_REP_COUNT_AX, "MAX_REP_COUNT", parameters["SIGNAL_COUNT_THRESHOLD"], parameters["SAMPLE_COUNT"], MAX_REP_COUNT, valstep=1)
# BG_NX_MIN_OVERLAP_SLIDER = Slider(BG_NX_MIN_OVERLAP_AX, "BG_NX_MIN_OVERLAP", parameters["MINIMUM_OVERLAP"][0], parameters["SAMPLE_COUNT"], BG_NX_MIN_OVERLAP, valstep=1)
# BG_NX_MAX_OVERLAP_SLIDER = Slider(BG_NX_MAX_OVERLAP_AX, "BG_NX_MAX_OVERLAP", parameters["MINIMUM_OVERLAP"][0], parameters["SAMPLE_COUNT"], BG_NX_MAX_OVERLAP, valstep=1)
LOG_FDR_SLIDER = Slider(LOG_FDR_AX, "LOG_FDR", -5, 0, valinit=LOG_FDR, valstep=0.01, valfmt="%1.3f")


ax_network = plt.axes(
    [
        0.45,
        bottom_space / button_count,
        0.1,
        (2 * bottom_space / button_count + bottom_space / (button_count + 1)) / 2
    ],
    facecolor=axcolor
)
network_radio = RadioButtons(ax_network, ('None', 'Selected', 'All'), active=0)


def update_network_plot(label):
    global NETWORK_VISIBLE
    NETWORK_VISIBLE = label
    plt_bg_nx()

network_radio.on_clicked(update_network_plot)




ax_label_selection = plt.axes(
    [
        0.45,
        bottom_space / button_count + (2 * bottom_space / button_count + bottom_space / (button_count + 1)) / 2,
        0.1,
        (2 * bottom_space / button_count + bottom_space / (button_count + 1)) / 2
    ],
    facecolor=axcolor
)
label_selection_radio = RadioButtons(ax_label_selection, ('None', 'm/z', 'Peptide'), active=0)


def update_label_selection_plot(label):
    global LABEL_SELECTION
    LABEL_SELECTION = label
    plt_anchor_selection()
    plt_intensities()

label_selection_radio.on_clicked(update_label_selection_plot)



ax_expand_neighbors = plt.axes([0.3, 3 * bottom_space / button_count, 0.1, bottom_space / (button_count + 1)], facecolor=axcolor)
neighbor_button = Button(ax_expand_neighbors, 'Expand neighbors', color=axcolor, hovercolor='0.975')


def expand_neighbors(val):
    selection = anchors["ION_COUNT"] >= MIN_REP_COUNT
    selection &= anchors["ION_COUNT"] <= MAX_REP_COUNT
    dt_low, dt_high = plt.ylim()
    rt_low, rt_high = plt.xlim()
    selection &= anchors["DT"] <= dt_high
    selection &= anchors["DT"] >= dt_low
    selection &= anchors["RT"] <= rt_high
    selection &= anchors["RT"] >= rt_low
    selection = np.nonzero(selection)
    n = np.unique(neighbors[selected_anchors].indices)
    n = n[np.isin(n, selection)]
    for anchor_ind in n:
        if anchor_ind not in selected_anchors:
            selected_anchors.append(anchor_ind)
    # print("Peptides: ", np.unique(anchor_peptide_match_counts[selected_anchors].indices, return_counts=True))
    plt_bg_nx()
    plt_anchor_selection()
    plt_intensities()

neighbor_button.on_clicked(expand_neighbors)


ax_refresh = plt.axes([0.3, bottom_space / button_count, 0.1, bottom_space / (button_count + 1)], facecolor=axcolor)
refresh_button = Button(ax_refresh, 'Refresh', color=axcolor, hovercolor='0.975')


def refresh(val):
    global selected_anchors
    selection = anchors["ION_COUNT"] >= MIN_REP_COUNT
    selection &= anchors["ION_COUNT"] <= MAX_REP_COUNT
    dt_low, dt_high = plt.ylim()
    rt_low, rt_high = plt.xlim()
    selection &= anchors["DT"] <= dt_high
    selection &= anchors["DT"] >= dt_low
    selection &= anchors["RT"] <= rt_high
    selection &= anchors["RT"] >= rt_low
    selection = np.nonzero(selection)
    selected_anchors = [selected_anchors[i] for i in np.flatnonzero(np.isin(selected_anchors, selection))]
    plt_bg_nx()
    plt_bg()
    plt_fg()
    plt_anchor_selection()
    plt_intensities()


refresh_button.on_clicked(refresh)


ax_clear = plt.axes([0.3, 2 * bottom_space / button_count, 0.1, bottom_space / (button_count + 1)], facecolor=axcolor)
clear_button = Button(ax_clear, 'Clear selection', color=axcolor, hovercolor='0.975')


def clear_selection(val):
    del selected_anchors[:]
    plt_bg_nx()
    plt_anchor_selection()
    plt_intensities()


clear_button.on_clicked(clear_selection)


ax_intensities = plt.axes([0.1, bottom_space, 0.15, 0.95 - bottom_space])
ax = plt.axes([0.3, bottom_space, 0.65, 0.95 - bottom_space])
plt_bg_nx = plot_bg_nx()
plt_bg = plot_bg()
plt_fg = plot_fg()
plt_anchor_selection = plot_anchor_selection()




def anchors_update(val):
    global MIN_REP_COUNT
    global MAX_REP_COUNT
    global LOG_FDR
    MIN_REP_COUNT = MIN_REP_COUNT_SLIDER.val
    MAX_REP_COUNT = MAX_REP_COUNT_SLIDER.val
    LOG_FDR = LOG_FDR_SLIDER.val
    refresh(None)

MIN_REP_COUNT_SLIDER.on_changed(anchors_update)
MAX_REP_COUNT_SLIDER.on_changed(anchors_update)
LOG_FDR_SLIDER.on_changed(anchors_update)


def onclick(event):
    if not event.dblclick:
        return
    rt = event.xdata
    dt = event.ydata
    selection = anchors["ION_COUNT"] >= MIN_REP_COUNT
    selection &= anchors["ION_COUNT"] <= MAX_REP_COUNT
    dt_low, dt_high = plt.ylim()
    rt_low, rt_high = plt.xlim()
    selection &= anchors["DT"] <= dt_high
    selection &= anchors["DT"] >= dt_low
    selection &= anchors["RT"] <= rt_high
    selection &= anchors["RT"] >= rt_low
    rt_distance = (anchors["RT"][selection] - rt) / (rt_high - rt_low)
    dt_distance = (anchors["DT"][selection] - dt) / (dt_high - dt_low)
    distance = np.sqrt(rt_distance**2 + dt_distance**2)
    min_distance_ind = np.argmin(distance)
    anchor_ind = np.flatnonzero(selection)[min_distance_ind]
    if event.key != "control":
        del selected_anchors[:]
    if anchor_ind in selected_anchors:
        sel_i = selected_anchors.index(anchor_ind)
        del selected_anchors[sel_i]
    else:
        selected_anchors.append(anchor_ind)
    # print("Peptides: ", np.unique(anchor_peptide_match_counts[selected_anchors].indices, return_counts=True))
    plt_bg_nx()
    plt_anchor_selection()
    plt_intensities()


def plot_intensities():
    norm = mpl.colors.Normalize(
        vmin=0,
        vmax=len(anchors)
    )
    m = cm.ScalarMappable(norm=norm, cmap="coolwarm")
    plots = []
    ax_intensities.set_ylim(
        np.log2(np.min(ions["CALIBRATED_INTENSITY"])),
        np.log2(np.max(ions["CALIBRATED_INTENSITY"])),
    )
    ax_intensities.set_xticks(
        list(range(parameters["SAMPLE_COUNT"]))
    )
    ax_intensities.set_xticklabels(
        parameters["APEX_FILE_NAMES"],
        rotation=45,
        ha="right"
    )
    def plot_intensities_after_init():
        for i in plots:
            i.remove()
        del plots[:]
        for anchor_index in selected_anchors:
            selected_ions = anchor_ions[anchor_index].data
            new_plot, = ax_intensities.plot(
                ions["SAMPLE"][selected_ions],
                np.log2(ions["CALIBRATED_INTENSITY"][selected_ions]),
                c=m.to_rgba(anchor_index),
                marker="."
            )
            plots.append(new_plot)
            annotation_location = (
                ions["SAMPLE"][selected_ions][-1],
                np.log2(ions["CALIBRATED_INTENSITY"][selected_ions])[-1]
            )
            if LABEL_SELECTION == "m/z":
                label = round(anchors["MZ"][anchor_index], 3)
            if LABEL_SELECTION == "Peptide":
                if anchor_peptide_match_counts.indptr[anchor_index] == anchor_peptide_match_counts.indptr[anchor_index + 1]:
                    continue
                peptide_index = anchor_peptide_match_counts.indices[
                    anchor_peptide_match_counts.indptr[anchor_index]
                ]
                label = src.peptides.getSequenceFromIndex(
                    peptide_index,
                    peptides,
                    peptide_index_matrix,
                    total_protein_sequence
                )
            if LABEL_SELECTION != "None":
                ann = ax_intensities.annotate(
                    label,
                    annotation_location
                )
                plots.append(ann)
    return plot_intensities_after_init


plt_intensities = plot_intensities()


ind = fig.canvas.mpl_connect('button_press_event', onclick)



# import importlib
# importlib.reload(src.browser)
