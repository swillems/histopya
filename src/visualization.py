#!venv/bin/python


import src.parameters
import src.io
import src.ions
import src.aggregates
import src.peptides
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import collections as mc
from matplotlib.widgets import Slider, RadioButtons, Button


class Browser(object):

    def __init__(self, parameters, **kwargs):
        self.parameters = parameters
        self.log = src.io.Log(self.parameters["LOG_FILE_NAME"][:-4] + "_browser.txt")
        self.importData(**kwargs)
        self.initWindows(**kwargs)
        # self.setInitGUIParameters(**kwargs)

    def initWindows(self, **kwargs):
        plt.ion()
        fig, tmp = plt.subplots()
        tmp.axis("off")
        self.axes = {
            "scatter": plt.axes([0.3, 0.35, 0.65, 0.6]),
            "line": plt.axes([0.05, 0.35, 0.2, 0.6]),
            "fdr_slider": plt.axes([0.3, 0.05, 0.1, 0.05]),
            "min_rep_slider": plt.axes([0.3, 0.15, 0.1, 0.05]),
            "max_rep_slider": plt.axes([0.3, 0.25, 0.1, 0.05]),
            "label_radio": plt.axes([0.45, 0.05, 0.1, 0.1]),
            "edge_radio": plt.axes([0.45, 0.05, 0.1, 0.1]),
        }
        # self.initScatter(**kwargs)
        # self.initIons(**kwargs)
        # self.initButtons(**kwargs)

    def importData(self, **kwargs):
        if "aggregates" in kwargs:
            self.anchors = kwargs["aggregates"]
        else:
            self.anchors = src.io.loadArray("ANCHORS_FILE_NAME", self.parameters)
        if "anchor_ions" in kwargs:
            self.anchor_ions = kwargs["anchor_ions"]
        else:
            self.anchor_ions = src.io.loadMatrix(
                "ANCHOR_IONS_FILE_NAME",
                self.parameters,
            )
        if "ions" in kwargs:
            self.ions = kwargs["ions"]
        else:
            self.ions = src.io.loadArray("IONS_FILE_NAME", self.parameters)
        if "neighbors" in kwargs:
            self.neighbors = kwargs["neighbors"]
        else:
            self.neighbors = src.io.loadMatrix(
                "ANCHOR_NEIGHBORS_FILE_NAME",
                self.parameters,
            )
            self.neighbors += self.neighbors.T
        if "anchor_peptide_score" in kwargs:
            self.anchor_peptide_score = kwargs["anchor_peptide_score"]
        else:
            self.anchor_peptide_scores = src.io.loadMatrix(
                "ANCHOR_PEPTIDE_SCORES_FILE_NAME",
                self.parameters,
            )
        if "anchor_peptide_match_counts" in kwargs:
            self.anchor_peptide_match_counts = kwargs["anchor_peptide_match_counts"]
        else:
            self.anchor_peptide_match_counts = src.io.loadMatrix(
                "ANCHOR_PEPTIDE_MATCH_COUNTS_FILE_NAME",
                self.parameters,
            )
        # self.proteins, self.total_protein_sequence, self.ptms, self.ptm_matrix = src.peptides.importProteinsAndPtms(
        #     self.parameters,
        #     self.log
        # )
        # self.peptides, self.peptide_index_matrix, self.digestion_matrix = src.peptides.digestProteins(
        #     self.proteins,
        #     self.total_protein_sequence,
        #     self.ptm_matrix,
        #     self.parameters,
        #     self.log
        # )
        # self.percolated_annotations = pd.read_csv(
        #     self.parameters["PERCOLATOR_TARGET_PIMS"],
        #     delimiter="\t"
        # )
        # self.percolated_fdrs = self.percolated_annotations.values[:, 2]
        # self.percolated_anchors, self.percolated_peptides = self.anchor_peptide_match_counts.nonzero()

    def setInitGUIParameters(self, **kwargs):
        self.SELECTED_ANCHORS = []
        self.LOG_FDR = -2
        self.MIN_REP_COUNT = self.parameters["SAMPLE_COUNT"]
        self.MAX_REP_COUNT = self.parameters["SAMPLE_COUNT"]
        self.NETWORK_VISIBLE = "None"
        self.LABEL_SELECTION = "None"
        self.ION_DIMENSION = "Logint"
        self.SAMPLE = "Aggregate"

    def updateVisibleNodes(self):
        dt_low, dt_high = plt.ylim()
        rt_low, rt_high = plt.xlim()
        selection = self.anchors["ION_COUNT"] >= self.MIN_REP_COUNT
        selection &= self.anchors["ION_COUNT"] <= self.MAX_REP_COUNT
        selection &= self.anchors["DT"] <= dt_high
        selection &= self.anchors["DT"] >= dt_low
        selection &= self.anchors["RT"] <= rt_high
        selection &= self.anchors["RT"] >= rt_low
        self.visible_nodes = np.flatnonzero(selection)

    def initScatter(**kwargs):
        self.initNodes()
        self.initSelectedNodes()
        self.initAnnotatedNodes()
        self.initEdges()
        plt_fg = browser.plot_fg()

    def initNodes(self):
        self.background_nodes = scatter_ax.scatter(
            self.anchors["RT"][selection],
            self.anchors["DT"][selection],
            c="lightgrey",
            marker=".",
            # s=anchors["ION_COUNT"][selection]
        )

    def initNetworkLineCollection(self):
        network_color_mapper = mpl.colors.Normalize(
            vmin=self.parameters["MINIMUM_OVERLAP"][0],
            vmax=self.parameters["SAMPLE_COUNT"]
        )
        self.network_color_mapper = cm.ScalarMappable(norm=network_color_mapper, cmap="RdYlGn")
        network_line_collection = mc.LineCollection([], [])
        self.network_line_collection = scatter_ax.add_collection(network_line_collection)

    def initButtons(**kwargs):
        pass

    def plot_bg_nx(self):
        def plot_bg_nx_after_init():
            selection = np.array([], dtype=np.int)
            if NETWORK_VISIBLE == "Selected":
                if len(self.SELECTED_ANCHORS) > 0:
                    selection = np.array(self.SELECTED_ANCHORS)
            if NETWORK_VISIBLE == "All":
                dt_low, dt_high = plt.ylim()
                rt_low, rt_high = plt.xlim()
                selection = self.anchors["ION_COUNT"] >= MIN_REP_COUNT
                selection &= self.anchors["ION_COUNT"] <= MAX_REP_COUNT
                selection &= self.anchors["DT"] <= dt_high
                selection &= self.anchors["DT"] >= dt_low
                selection &= self.anchors["RT"] <= rt_high
                selection &= self.anchors["RT"] >= rt_low
                selection = np.flatnonzero(selection)
            n = self.neighbors[selection].T.tocsr()[selection]
            a, b = n.nonzero()
            c = a > b
            a = a[c]
            b = b[c]
            overlap = n.data[c]
            order = np.argsort(overlap)
            a = selection[a[order]]
            b = selection[b[order]]
            overlap = overlap[order]
            start_edges = list(zip(self.anchors["RT"][a], self.anchors["DT"][a]))
            end_edges = list(zip(self.anchors["RT"][b], self.anchors["DT"][b]))
            bg_nx_lc.set_segments(list(zip(start_edges, end_edges)))
            clrs = m.to_rgba(overlap)
            bg_nx_lc.set_color(clrs)
            # lc = mc.LineCollection(list(zip(start_edges, end_edges)), colors=clrs)
        return plot_bg_nx_after_init

    def plot_anchor_selection(self):
        anchor_selection = ax.scatter(
            self.anchors["RT"][self.SELECTED_ANCHORS],
            self.anchors["DT"][self.SELECTED_ANCHORS],
            facecolor="None",
            edgecolor='black',
            # s=anchors["ION_COUNT"][self.SELECTED_ANCHORS],
        )
        annotations = []
        def plot_anchor_selection_after_init():
            anchor_selection.set_offsets(
                np.c_[
                    self.anchors["RT"][self.SELECTED_ANCHORS],
                    self.anchors["DT"][self.SELECTED_ANCHORS]
                ]
            )
            # ax.annotate(anchors["MZ"][self.SELECTED_ANCHORS], (anchors["RT"][self.SELECTED_ANCHORS], anchors["DT"][self.SELECTED_ANCHORS]))
            for anchor_index in annotations:
                anchor_index.remove()
            del annotations[:]
            for anchor_index in self.SELECTED_ANCHORS:
                rt = self.anchors["RT"][anchor_index]
                dt = self.anchors["DT"][anchor_index]
                if LABEL_SELECTION == "m/z":
                    label = round(self.anchors["MZ"][anchor_index], 3)
                if LABEL_SELECTION == "Peptide":
                    if self.anchor_peptide_match_counts.indptr[anchor_index] == self.anchor_peptide_match_counts.indptr[anchor_index + 1]:
                        continue
                    peptide_index = self.anchor_peptide_match_counts.indices[
                        self.anchor_peptide_match_counts.indptr[anchor_index]
                    ]
                    label = src.peptides.getSequenceFromIndex(
                        peptide_index,
                        self.peptides,
                        self.peptide_index_matrix,
                        self.total_protein_sequence
                    )
                if LABEL_SELECTION != "None":
                    ann = ax.annotate(label, (rt, dt))
                    annotations.append(ann)
        return plot_anchor_selection_after_init

    def plot_bg(self):
        def plot_bg_after_init():
            selection = self.anchors["ION_COUNT"] >= MIN_REP_COUNT
            selection &= self.anchors["ION_COUNT"] <= MAX_REP_COUNT
            dt_low, dt_high = plt.ylim()
            rt_low, rt_high = plt.xlim()
            selection &= self.anchors["DT"] <= dt_high
            selection &= self.anchors["DT"] >= dt_low
            selection &= self.anchors["RT"] <= rt_high
            selection &= self.anchors["RT"] >= rt_low
            bg_sc.set_offsets(
                np.c_[
                    self.anchors["RT"][selection],
                    self.anchors["DT"][selection]
                ]
            )
        return plot_bg_after_init

    def plot_fg(self):
        significant_pims = self.percolated_annotations.values[
            self.percolated_fdrs <= 10**LOG_FDR,
            0
        ].astype(int)
        significant_anchors = self.percolated_anchors[significant_pims]
        significant_peptides = self.percolated_peptides[significant_pims]
        selection = self.anchors["ION_COUNT"][significant_anchors] >= MIN_REP_COUNT
        selection &= self.anchors["ION_COUNT"][significant_anchors] <= MAX_REP_COUNT
        fg_sc = ax.scatter(
            self.anchors["RT"][significant_anchors[selection]],
            self.anchors["DT"][significant_anchors[selection]],
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
            significant_pims = self.percolated_annotations.values[
                self.percolated_fdrs <= 10**LOG_FDR,
                0
            ].astype(int)
            significant_anchors = self.percolated_anchors[significant_pims]
            significant_peptides = self.percolated_peptides[significant_pims]
            selection = self.anchors["ION_COUNT"][significant_anchors] >= MIN_REP_COUNT
            selection &= self.anchors["ION_COUNT"][significant_anchors] <= MAX_REP_COUNT
            dt_low, dt_high = plt.ylim()
            rt_low, rt_high = plt.xlim()
            selection &= self.anchors["DT"][significant_anchors] <= dt_high
            selection &= self.anchors["DT"][significant_anchors] >= dt_low
            selection &= self.anchors["RT"][significant_anchors] <= rt_high
            selection &= self.anchors["RT"][significant_anchors] >= rt_low
            fg_sc.set_offsets(
                np.c_[
                    self.anchors["RT"][significant_anchors[selection]],
                    self.anchors["DT"][significant_anchors[selection]]
                ]
            )
            clrs = m.to_rgba(significant_peptides[selection])
            fg_sc.set_color(clrs)
        return plot_fg_after_init

    def plot_intensities(self):
        norm = mpl.colors.Normalize(
            vmin=0,
            vmax=len(self.anchors)
        )
        m = cm.ScalarMappable(norm=norm, cmap="coolwarm")
        plots = []
        ax_intensities.set_ylim(
            np.log2(np.min(self.ions["CALIBRATED_INTENSITY"])),
            np.log2(np.max(self.ions["CALIBRATED_INTENSITY"])),
        )
        ax_intensities.set_xticks(
            list(range(self.parameters["SAMPLE_COUNT"]))
        )
        ax_intensities.set_xticklabels(
            self.parameters["APEX_FILE_NAMES"],
            rotation=45,
            ha="right"
        )
        def plot_intensities_after_init():
            for i in plots:
                i.remove()
            del plots[:]
            for anchor_index in self.SELECTED_ANCHORS:
                selected_ions = self.anchor_ions[anchor_index].data
                new_plot, = ax_intensities.plot(
                    self.ions["SAMPLE"][selected_ions],
                    np.log2(self.ions["CALIBRATED_INTENSITY"][selected_ions]),
                    c=m.to_rgba(anchor_index),
                    marker="."
                )
                plots.append(new_plot)
                annotation_location = (
                    self.ions["SAMPLE"][selected_ions][-1],
                    np.log2(self.ions["CALIBRATED_INTENSITY"][selected_ions])[-1]
                )
                if LABEL_SELECTION == "m/z":
                    label = round(self.anchors["MZ"][anchor_index], 3)
                if LABEL_SELECTION == "Peptide":
                    if self.anchor_peptide_match_counts.indptr[anchor_index] == self.anchor_peptide_match_counts.indptr[anchor_index + 1]:
                        continue
                    peptide_index = self.anchor_peptide_match_counts.indices[
                        self.anchor_peptide_match_counts.indptr[anchor_index]
                    ]
                    label = src.peptides.getSequenceFromIndex(
                        peptide_index,
                        self.peptides,
                        self.peptide_index_matrix,
                        self.total_protein_sequence
                    )
                if LABEL_SELECTION != "None":
                    ann = ax_intensities.annotate(
                        label,
                        annotation_location
                    )
                    plots.append(ann)
        return plot_intensities_after_init




parameter_file_name = "data/test/parameters.json"
parameters = src.parameters.importParameterDictFromJSON(parameter_file_name)
browser = Browser(parameters)
