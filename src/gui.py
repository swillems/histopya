#!venv/bin/python


import matplotlib
matplotlib.use("tkAgg", warn=False)
import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg as tkagg
import src.parameters
import src.io
import src.ions
import src.aggregates
import src.peptides
import pandas as pd
import numpy as np
import scipy


class GUI(object):

    def __init__(self, dataset):
        self.dataset = dataset
        self.root = tk.Tk()
        self.root.geometry("1500x800")
        self.root.winfo_toplevel().title("Ion-network Browser")
        self.root_frame = tk.Frame(self.root, bd=1, relief="raised")
        self.root_frame.pack(fill="both", expand=1)
        self.__initOptionsFrame()
        self.__initAggregatesFrame()
        self.__initIonsFrame()
        self.refresh()
        self.start()

    def start(self):
        self.root.mainloop()

    def __initOptionsFrame(self):
        self.options_frame = tk.Frame(
            self.root_frame,
            bd=1,
            relief="raised"
        )
        self.options_frame.pack(
            side="left",
            fill="x"
        )
        self.options = {}
        self.addLabelOption()
        self.addAxisOption()
        # self.addViewTypeOption()
        self.addShowEdgesOption()
        self.addMinimumSignalOption()
        self.addMaximumSignalOption()
        self.addFDROption()
        self.addSelectAllVisibleOption()
        self.addUnselectAllVisibleOption()
        self.addExpandNeighborOption()
        self.addRefreshOption()
        for option in self.options.values():
            option.pack(anchor="w", fill="x", expand=1)

    def addViewTypeOption(self):
        self.options["View type"] = tk.Frame(
            self.options_frame,
            bd=1,
            relief="raised"
        )
        tk.Label(
            self.options["View type"],
            text="View type",
        ).pack(anchor="w")
        view_types = [
            "Aggregates"
        ] + [
            name.split("/")[-1][:-4] for name in self.dataset.parameters["APEX_FILE_NAMES"]
        ]
        self.view_type = tk.IntVar()
        self.view_type.set(-1)
        for view_index, view_type in enumerate(view_types):
            tk.Radiobutton(
                self.options["View type"],
                text=view_type,
                value=view_index - 1,
                variable=self.view_type,
                command=self.refresh
            ).pack(anchor="w")

    def addRefreshOption(self):
        self.options["Refresh"] = tk.Button(
            self.options_frame,
            text="Refresh",
            command=self.refresh,
            anchor="w"
        )

    def addSelectAllVisibleOption(self):

        def onclick():
            self.dataset.selectAllVisible()
            self.dataset.updateLabelSelection(self)
            self.dataset.plotSelectedNodes(self)
            self.dataset.plotIons(self)
            self.refreshAggregateCanvas()
            self.refreshIonCanvas()

        self.options["Select all visible nodes"] = tk.Button(
            self.options_frame,
            text="Select all visible nodes",
            command=onclick,
            anchor="w"
        )

    def addUnselectAllVisibleOption(self):

        def onclick():
            self.dataset.unselectAllVisible()
            self.dataset.updateLabelSelection(self)
            self.dataset.plotSelectedNodes(self)
            self.dataset.plotIons(self)
            self.refreshAggregateCanvas()
            self.refreshIonCanvas()

        self.options["Unselect all visible nodes"] = tk.Button(
            self.options_frame,
            text="Unselect all visible nodes",
            command=onclick,
            anchor="w"
        )

    def addExpandNeighborOption(self):

        def onclick():
            self.dataset.selectVisibleNeighbors()
            self.dataset.updateLabelSelection(self)
            self.dataset.plotSelectedNodes(self)
            self.dataset.plotIons(self)
            self.refreshAggregateCanvas()
            self.refreshIonCanvas()

        self.options["Expand neighbors"] = tk.Button(
            self.options_frame,
            text="Expand neighbors",
            anchor="w",
            command=onclick
        )

    def addShowEdgesOption(self):

        def onclick():
            self.dataset.plotEdges(self)
            self.refreshAggregateCanvas()

        self.show_edges = tk.IntVar()
        self.options["Show edges"] = tk.Checkbutton(
            self.options_frame,
            text="Show edges",
            variable=self.show_edges,
            anchor="w",
            command=onclick
        )

    def addLabelOption(self):
        self.options["Label type"] = tk.Frame(
            self.options_frame,
            bd=1,
            relief="raised"
        )
        tk.Label(
            self.options["Label type"],
            text="Label type",
        ).pack(anchor="w")
        label_types = [
            "None",
            "Peptide",
            "Protein",
            "m/z",
            "rt",
            "dt",
            "Index",
        ]
        self.label_type = tk.StringVar()
        self.label_type.set(label_types[0])

        def onclick():
            self.dataset.updateLabelSelection(self)
            # self.dataset.plotSelectedNodes(self)
            # self.dataset.plotIons(self)
            self.refreshAggregateCanvas()
            self.refreshIonCanvas()

        for label_type in label_types:
            tk.Radiobutton(
                self.options["Label type"],
                text=label_type,
                variable=self.label_type,
                value=label_type,
                command=onclick
            ).pack(anchor="w")

    def addAxisOption(self):
        self.options["Ion axis type"] = tk.Frame(
            self.options_frame,
            bd=1,
            relief="raised"
        )
        tk.Label(
            self.options["Ion axis type"],
            text="Ion axis type",
        ).pack(anchor="w")
        label_types = [
            "Log intensity",
            "rt",
            "dt",
        ]
        self.axis_type = tk.StringVar()
        self.axis_type.set(label_types[0])

        def onclick():
            axis_types = {
                "Log intensity": "CALIBRATED_INTENSITY",
                "dt": "DT",
                "rt": "RT",
            }
            axis_type = axis_types[self.getAxisType()]
            axis_min_lim = np.min(self.dataset.ions[axis_type])
            axis_max_lim = np.max(self.dataset.ions[axis_type])
            if axis_type == "CALIBRATED_INTENSITY":
                axis_min_lim = np.log2(axis_min_lim)
                axis_max_lim = np.log2(axis_max_lim)
            self.ion_ax.set_ylim(
                axis_min_lim,
                axis_max_lim,
            )
            self.ion_ax.set_ylabel(axis_type)
            self.dataset.labelIons(self)
            self.dataset.plotIons(self)
            self.refreshIonCanvas()

        for label_type in label_types:
            tk.Radiobutton(
                self.options["Ion axis type"],
                text=label_type,
                variable=self.axis_type,
                value=label_type,
                command=onclick
            ).pack(anchor="w")

    def addFDROption(self):
        self.options["Log FDR threshold"] = tk.Frame(
            self.options_frame,
            bd=1,
            relief="raised"
        )
        tk.Label(
            self.options["Log FDR threshold"],
            text="Log FDR threshold",
        ).pack(anchor="w")

        def onclick(*args):
            label_type = self.getLabelType()
            if (label_type == "Peptide") or (label_type == "Protein"):
                self.dataset.updateLabelSelection(self)
                self.refreshIonCanvas()
            self.dataset.plotAnnotatedNodes(self)
            self.refreshAggregateCanvas()

        self.fdr_threshold_slider = tk.Scale(
            self.options["Log FDR threshold"],
            from_=-5,
            resolution=0.1,
            to=0,
            orient="horizontal",
            command=onclick
        )
        self.fdr_threshold_slider.set(-2)
        self.fdr_threshold_slider.pack(anchor="w")

    def addMinimumSignalOption(self):
        self.options["Minimum signal"] = tk.Frame(
            self.options_frame,
            bd=1,
            relief="raised"
        )
        tk.Label(
            self.options["Minimum signal"],
            text="Minimum signal",
        ).pack(anchor="w")
        self.minimum_replicate_count_slider = tk.Scale(
            self.options["Minimum signal"],
            from_=self.dataset.parameters["SIGNAL_COUNT_THRESHOLD"],
            to=self.dataset.parameters["SAMPLE_COUNT"],
            orient="horizontal",
            command=lambda arg: self.refresh()
        )
        self.minimum_replicate_count_slider.set(
            self.dataset.parameters["SAMPLE_COUNT"]
        )
        self.minimum_replicate_count_slider.pack(anchor="w")

    def addMaximumSignalOption(self):
        self.options["Maximum signal"] = tk.Frame(
            self.options_frame,
            bd=1,
            relief="raised"
        )
        tk.Label(
            self.options["Maximum signal"],
            text="Maximum signal",
        ).pack(anchor="w")
        self.maximum_replicate_count_slider = tk.Scale(
            self.options["Maximum signal"],
            from_=self.dataset.parameters["SIGNAL_COUNT_THRESHOLD"],
            to=self.dataset.parameters["SAMPLE_COUNT"],
            orient="horizontal",
            command=lambda arg: self.refresh()
        )
        self.maximum_replicate_count_slider.set(
            self.dataset.parameters["SAMPLE_COUNT"]
        )
        self.maximum_replicate_count_slider.pack(anchor="w")

    def __initAggregatesFrame(self):
        self.aggregate_frame = tk.Frame(
            self.root_frame,
            bd=1,
            relief="raised"
        )
        self.aggregate_frame.pack(
            side="left",
            fill="both",
            expand=True
        )
        self.aggregate_fig = plt.Figure()
        self.aggregate_ax = self.aggregate_fig.add_subplot(111)
        self.aggregate_canvas = tkagg.FigureCanvasTkAgg(self.aggregate_fig, self.aggregate_frame)
        self.aggregate_canvas.get_tk_widget().pack(
            fill='both',
            expand=True
        )
        self.aggregate_toolbar = tkagg.NavigationToolbar2Tk(
            self.aggregate_canvas,
            self.aggregate_frame
        )
        default_release = self.aggregate_toolbar.release

        def release(event):
            default_release(event)
            self.refresh()

        self.aggregate_toolbar.release = release
        self.aggregate_ax.set_xlim(
            [
                np.min(self.dataset.anchors["RT"]),
                np.max(self.dataset.anchors["RT"]),
            ]
        )
        self.aggregate_ax.set_ylim(
            [
                np.min(self.dataset.anchors["DT"]),
                np.max(self.dataset.anchors["DT"]),
            ]
        )
        self.aggregate_ax.set_xlabel("RT")
        self.aggregate_ax.set_ylabel("DT")
        self.aggregate_ax.set_title("Aggregates")

        def onclick(event):
            if (event is None) or (not event.dblclick):
                return
            self.dataset.updateSelectedNodes(self, event)
            self.dataset.plotSelectedNodes(self)
            self.dataset.plotIons(self)
            self.refreshAggregateCanvas()
            self.refreshIonCanvas()
            self.aggregate_canvas.mpl_disconnect(self.click_connection)
            self.click_connection = self.aggregate_canvas.mpl_connect(
                'button_press_event',
                onclick
            )

        self.click_connection = self.aggregate_canvas.mpl_connect(
            'button_press_event',
            onclick
        )

    def __initIonsFrame(self):
        self.ion_frame = tk.Frame(
            self.root_frame,
            bd=1,
            relief="raised"
        )
        self.ion_frame.pack(
            side="left",
            fill="y",
        )
        self.ion_fig = plt.Figure()
        self.ion_ax = self.ion_fig.add_subplot(111)
        self.ion_canvas = tkagg.FigureCanvasTkAgg(self.ion_fig, self.ion_frame)
        self.ion_toolbar = tkagg.NavigationToolbar2Tk(
            self.ion_canvas,
            self.ion_frame
        )
        self.ion_canvas.get_tk_widget().pack(
            fill='both',
            expand=True
        )
        self.getAxisType()
        self.ion_ax.set_title("Ions")

    def refresh(self):
        self.dataset.updateVisibleNodes(self)
        self.dataset.plotVisibleNodes(self)
        self.dataset.plotAnnotatedNodes(self)
        self.dataset.plotSelectedNodes(self)
        self.dataset.plotEdges(self)
        self.refreshAggregateCanvas()
        self.dataset.plotIons(self)
        self.refreshIonCanvas()

    def refreshAggregateCanvas(self):
        self.aggregate_canvas.draw()
        self.aggregate_canvas.flush_events()

    def refreshIonCanvas(self):
        self.ion_canvas.draw()
        self.ion_canvas.flush_events()

    def getMinimumReplicateCount(self):
        return self.minimum_replicate_count_slider.get()

    def getMaximumReplicateCount(self):
        return self.maximum_replicate_count_slider.get()

    def getLabelType(self):
        return self.label_type.get()

    def getViewType(self):
        try:
            return self.view_type.get()
        except AttributeError:
            return -1

    def getFDRThreshold(self):
        return 10 ** self.fdr_threshold_slider.get()

    def getAxisType(self):
        return self.axis_type.get()

    def getShowEdges(self):
        return self.show_edges.get()

    def getVisibleBoundaries(self):
        dt_low, dt_high = self.aggregate_ax.get_ylim()
        rt_low, rt_high = self.aggregate_ax.get_xlim()
        return rt_low, rt_high, dt_low, dt_high

    def dummyCommand(self, *args, **kwargs):
        print(self.getMinimumReplicateCount())
        print(self.getMaximumReplicateCount())
        print(self.getLabelType())
        print(self.getViewType())
        print(self.getFDRThreshold())
        print(self.getAxisType())
        print(self.getShowEdges())


class Dataset(object):

    def __init__(self, parameters, **kwargs):
        if isinstance(parameters, str):
            parameters = src.parameters.importParameterDictFromJSON(parameters)
        self.parameters = parameters
        self.log = src.io.Log(self.parameters["LOG_FILE_NAME"][:-4] + "_browser.txt")
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
        try:
            if "anchor_peptide_score" in kwargs:
                self.anchor_peptide_scores = kwargs["anchor_peptide_score"]
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
            self.percolated_annotations = pd.read_csv(
                self.parameters["PERCOLATOR_TARGET_PIMS"],
                delimiter="\t"
            )
            self.percolated_fdrs = self.percolated_annotations.values[:, 2]
            self.percolated_anchors, self.percolated_peptides = self.anchor_peptide_match_counts.nonzero()
        except FileNotFoundError:
            self.anchor_peptide_scores = scipy.sparse.csr_matrix([])
            self.anchor_peptide_match_counts = scipy.sparse.csr_matrix([])
            self.percolated_annotations = pd.DataFrame()
            self.percolated_fdrs = np.array([], dtype=np.int)
            self.percolated_anchors = np.array([], dtype=np.int)
            self.percolated_peptides = np.array([], dtype=np.int)
        self.proteins, self.total_protein_sequence, self.ptms, self.ptm_matrix = src.peptides.importProteinsAndPtms(
            self.parameters,
            self.log
        )
        self.peptides, self.peptide_index_matrix, self.digestion_matrix = src.peptides.digestProteins(
            self.proteins,
            self.total_protein_sequence,
            self.ptm_matrix,
            self.parameters,
            self.log
        )
        self.visible_nodes = np.array([], dtype=np.int)
        self.selected_nodes = np.array([], dtype=np.bool)
        self.node_labels = np.array([], dtype=np.bool)

    def updateVisibleNodes(self, gui):
        self.log.printMessage("Updating visible nodes")
        previous_selection = self.visible_nodes[self.selected_nodes]
        rt_low, rt_high, dt_low, dt_high = gui.getVisibleBoundaries()
        selection = self.anchors["ION_COUNT"] >= gui.getMinimumReplicateCount()
        selection &= self.anchors["ION_COUNT"] <= gui.getMaximumReplicateCount()
        selection &= self.anchors["DT"] <= dt_high
        selection &= self.anchors["DT"] >= dt_low
        selection &= self.anchors["RT"] <= rt_high
        selection &= self.anchors["RT"] >= rt_low
        self.visible_nodes = np.flatnonzero(selection)
        sample = gui.getViewType()
        if sample != -1:
            sample_ions = self.anchor_ions[
                self.visible_nodes, sample
            ].todense().A.flatten()
            self.visible_nodes = self.visible_nodes[
                np.flatnonzero(sample_ions)
            ]
        self.selected_nodes = np.isin(
            self.visible_nodes,
            previous_selection,
            assume_unique=True
        )
        self.updateLabelSelection(gui)

    def updateSelectedNodes(self, gui, event):
        self.log.printMessage("Updating selected nodes")
        rt = event.xdata
        dt = event.ydata
        rt_low, rt_high, dt_low, dt_high = gui.getVisibleBoundaries()
        rt_distance = (
            self.anchors["RT"][self.visible_nodes] - rt
        ) / (rt_high - rt_low)
        dt_distance = (
            self.anchors["DT"][self.visible_nodes] - dt
        ) / (dt_high - dt_low)
        distance = np.sqrt(rt_distance**2 + dt_distance**2)
        if event.key != "control":
            self.unselectAllVisible()
        selected_index = np.argmin(distance)
        self.selected_nodes[selected_index] = not self.selected_nodes[selected_index]
        self.updateLabelSelection(gui)

    def selectAllVisible(self):
        self.log.printMessage("Selecting all visible nodes")
        self.selected_nodes[:] = True

    def unselectAllVisible(self):
        self.log.printMessage("Unselecting all visible nodes")
        self.selected_nodes[:] = False

    def selectVisibleNeighbors(self):
        selected_neighbors = np.unique(
            self.neighbors[
                self.visible_nodes[
                    self.selected_nodes
                ]
            ].indices
        )
        self.selected_nodes |= np.isin(self.visible_nodes, selected_neighbors)

    def plotVisibleNodes(self, gui):
        self.log.printMessage("Initializing visible nodes")
        gui.background_nodes = gui.aggregate_ax.scatter(
            self.anchors["RT"][self.visible_nodes],
            self.anchors["DT"][self.visible_nodes],
            c="lightgrey",
            marker=".",
        )
        self.plotVisibleNodes = self.plotVisibleNodesPostInit

    def plotVisibleNodesPostInit(self, gui):
        self.log.printMessage("Plotting visible nodes")
        gui.background_nodes.set_offsets(
            np.c_[
                self.anchors["RT"][self.visible_nodes],
                self.anchors["DT"][self.visible_nodes]
            ]
        )

    def plotSelectedNodes(self, gui):
        self.log.printMessage("Initializing selected nodes")
        gui.selected_nodes = gui.aggregate_ax.scatter(
            self.anchors["RT"][self.visible_nodes[self.selected_nodes]],
            self.anchors["DT"][self.visible_nodes[self.selected_nodes]],
            facecolor="None",
            edgecolor='black',
        )
        self.plotSelectedNodes = self.plotSelectedNodesPostInit

    def plotSelectedNodesPostInit(self, gui):
        self.log.printMessage("Plotting selected nodes")
        gui.selected_nodes.set_offsets(
            np.c_[
                self.anchors["RT"][self.visible_nodes[self.selected_nodes]],
                self.anchors["DT"][self.visible_nodes[self.selected_nodes]]
            ]
        )

    def plotEdges(self, gui):
        self.log.printMessage("Initializing visible edges")
        gui.edge_color_mapper = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(
                vmin=self.parameters["MINIMUM_OVERLAP"][0],
                vmax=self.parameters["SAMPLE_COUNT"]
            ),
            cmap="RdYlGn"
        )
        gui.edge_collection = gui.aggregate_ax.add_collection(
            matplotlib.collections.LineCollection([], [])
        )
        self.plotEdges = self.plotEdgesPostInit

    def plotEdgesPostInit(self, gui):
        self.log.printMessage("Plotting visible edges")
        selection = np.array([], dtype=np.int)
        if gui.getShowEdges():
            selection = self.visible_nodes
        selected_neighbors = self.neighbors[selection].T.tocsr()[selection]
        a, b = selected_neighbors.nonzero()
        c = a > b
        a = a[c]
        b = b[c]
        overlap = selected_neighbors.data[c]
        order = np.argsort(overlap)
        a = selection[a[order]]
        b = selection[b[order]]
        overlap = overlap[order]
        start_edges = list(zip(self.anchors["RT"][a], self.anchors["DT"][a]))
        end_edges = list(zip(self.anchors["RT"][b], self.anchors["DT"][b]))
        gui.edge_collection.set_segments(list(zip(start_edges, end_edges)))
        gui.edge_collection.set_color(
            gui.edge_color_mapper.to_rgba(overlap)
        )

    def plotIons(self, gui):
        self.log.printMessage("Initializing ions")
        gui.ion_color_mapper = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(
                vmin=0,
                vmax=len(self.anchors)
            ),
            cmap="tab10"
        )
        gui.ion_plots = []
        gui.ion_ax.set_xticks(
            list(range(self.parameters["SAMPLE_COUNT"]))
        )
        gui.ion_ax.set_ylim(
            np.log2(np.min(self.ions["CALIBRATED_INTENSITY"])),
            np.log2(np.max(self.ions["CALIBRATED_INTENSITY"])),
        )
        gui.ion_ax.set_xticklabels(
            [name.split("/")[-1][:-4] for name in self.parameters["APEX_FILE_NAMES"]],
            rotation=45,
            ha="right"
        )
        self.plotIons = self.plotIonsPostInit

    def plotIonsPostInit(self, gui):
        self.log.printMessage("Plotting ions")
        for i in gui.ion_plots:
            i.remove()
        del gui.ion_plots[:]
        axis_types = {
            "Log intensity": "CALIBRATED_INTENSITY",
            "dt": "DT",
            "rt": "RT",
        }
        axis_type = axis_types[gui.getAxisType()]
        for anchor_index in self.visible_nodes[self.selected_nodes]:
            selected_ions = self.anchor_ions[anchor_index].data
            values = self.ions[axis_type][selected_ions]
            if axis_type == "CALIBRATED_INTENSITY":
                values = np.log2(values)
            new_plot, = gui.ion_ax.plot(
                self.ions["SAMPLE"][selected_ions],
                values,
                c=gui.ion_color_mapper.to_rgba(anchor_index),
                marker="."
            )
            gui.ion_plots.append(new_plot)
        pass

    def getSignificantPimsAnchorsAndPeptides(self, gui):
        try:
            significant_pims = self.percolated_annotations.values[
                self.percolated_fdrs <= gui.getFDRThreshold(),
                0
            ].astype(int)
        except IndexError:
            significant_pims = np.array([], dtype=np.int)
        significant_anchors = self.percolated_anchors[significant_pims]
        significant_peptides = self.percolated_peptides[significant_pims]
        return significant_pims, significant_anchors, significant_peptides

    def plotAnnotatedNodes(self, gui):
        (
            significant_pims,
            significant_anchors,
            significant_peptides
        ) = self.getSignificantPimsAnchorsAndPeptides(gui)
        gui.annotation_color_mapper = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(
                vmin=0,
                vmax=self.anchor_peptide_match_counts.shape[1]
            ),
            cmap="tab10"
        )
        selection = np.isin(significant_anchors, self.visible_nodes)
        selected_nodes = significant_anchors[selection]
        gui.annotated_nodes = gui.aggregate_ax.scatter(
            self.anchors["RT"][selected_nodes],
            self.anchors["DT"][selected_nodes],
            marker=".",
        )
        gui.annotated_nodes.set_color(
            gui.annotation_color_mapper.to_rgba(significant_peptides[selection])
        )
        self.plotAnnotatedNodes = self.plotAnnotatedNodesPostInit

    def plotAnnotatedNodesPostInit(self, gui):
        (
            significant_pims,
            significant_anchors,
            significant_peptides
        ) = self.getSignificantPimsAnchorsAndPeptides(gui)
        selection = np.isin(significant_anchors, self.visible_nodes)
        selected_nodes = significant_anchors[selection]
        gui.annotated_nodes.set_offsets(
            np.c_[
                self.anchors["RT"][selected_nodes],
                self.anchors["DT"][selected_nodes]
            ]
        )
        gui.annotated_nodes.set_color(
            gui.annotation_color_mapper.to_rgba(significant_peptides[selection])
        )

    def updateLabelSelection(self, gui):
        self.log.printMessage("Updating labels")
        label_type = gui.getLabelType()
        if label_type == "None":
            self.node_labels = np.array([""] * len(self.visible_nodes))
        if label_type == "m/z":
            self.node_labels = self.anchors["MZ"][self.visible_nodes]
        if label_type == "rt":
            self.node_labels = self.anchors["RT"][self.visible_nodes]
        if label_type == "dt":
            self.node_labels = self.anchors["DT"][self.visible_nodes]
        if label_type == "Index":
            self.node_labels = self.visible_nodes
        if (label_type == "Peptide") or (label_type == "Protein"):
            peptide_indices = self.anchor_peptide_match_counts.indices[
                self.anchor_peptide_match_counts.indptr[self.visible_nodes]
            ]
            (
                significant_pims,
                significant_anchors,
                significant_peptides
            ) = self.getSignificantPimsAnchorsAndPeptides(gui)
            selection = np.isin(self.visible_nodes, significant_anchors)
            peptide_indices[~selection] = -1
            if (label_type == "Peptide"):
                self.node_labels = src.peptides.getPeptideSequences(
                    peptide_indices,
                    self.peptides,
                    self.peptide_index_matrix,
                    self.total_protein_sequence
                )
            if (label_type == "Protein"):
                self.node_labels = src.peptides.getProteinAccessions(
                    peptide_indices,
                    self.peptides,
                    self.proteins
                )
        self.labelIons(gui)
        self.labelAggregates(gui)

    def labelIons(self, gui):
        self.log.printMessage("Labeling ions")
        for child in gui.ion_ax.get_children():
            if isinstance(
                child, matplotlib.text.Annotation
            ):
                child.remove()
        axis_types = {
            "Log intensity": "CALIBRATED_INTENSITY",
            "dt": "DT",
            "rt": "RT",
        }
        axis_type = axis_types[gui.getAxisType()]
        for label, anchor_index in zip(
            self.node_labels[self.selected_nodes],
            self.visible_nodes[self.selected_nodes]
        ):
            selected_ions = self.anchor_ions[anchor_index].data
            sample = self.ions["SAMPLE"][selected_ions][-1]
            y_coordinate = self.ions[axis_type][selected_ions][-1]
            if axis_type == "CALIBRATED_INTENSITY":
                y_coordinate = np.log2(y_coordinate)
            gui.ion_ax.annotate(
                label,
                (
                    sample,
                    y_coordinate
                )
            )

    def labelAggregates(self, gui):
        self.log.printMessage("Labeling nodes")
        for child in gui.aggregate_ax.get_children():
            if isinstance(
                child, matplotlib.text.Annotation
            ):
                child.remove()
        rts = self.anchors["RT"][self.visible_nodes[self.selected_nodes]]
        dts = self.anchors["DT"][self.visible_nodes[self.selected_nodes]]
        labels = self.node_labels[self.selected_nodes]
        for rt, dt, label in zip(rts, dts, labels):
            gui.aggregate_ax.annotate(label, (rt, dt))


if __name__ == "__main__":
    import src.parameters
    parameters, action = src.parameters.parseFromCommandLine()
    dataset = Dataset(parameters)
    gui = GUI(dataset)

# import src.gui; import importlib; importlib.reload(src.gui); d=src.gui.Dataset("data/tenzer/parameters.json"); g=src.gui.GUI(d)
