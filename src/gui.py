#!venv/bin/python

import tkinter as tk
from matplotlib import pyplot as plt
import matplotlib.backends.backend_tkagg as tkagg
import src.io
import numpy as np


class GUI(object):

    def __init__(self, dataset):
        self.dataset = dataset
        self.root = tk.Tk()
        self.root.geometry("1500x800")
        self.root.winfo_toplevel().title("HistoPyA Browser")
        self.root_frame = tk.Frame(self.root, bd=1, relief="raised")
        self.root_frame.pack(fill="both", expand=1)
        self.initOptionsFrame()
        self.initAggregatesFrame()
        self.initIonsFrame()
        self.root.mainloop()

    def initOptionsFrame(self):
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
        self.addExpandNeighborOption()
        self.addToggleNetworkOption()
        self.addAggregateLabelOption()
        self.addIonAxisOption()
        self.addFDROption()
        self.addMinimumSignalOption()
        self.addMaximumSignalOption()
        for option in self.options.values():
            option.pack(anchor="w", fill="x", expand=1)

    def addExpandNeighborOption(self):
        self.options["Expand neighbors"] = tk.Button(
            self.options_frame,
            text="Expand neighbors",
            command=lambda test=1: self.dummyCommand(test),
            anchor="w"
        )

    def addToggleNetworkOption(self):
        self.options["Toggle network"] = tk.Button(
            self.options_frame,
            text="Toggle network",
            command=lambda test=1: self.dummyCommand(test),
            anchor="w"
        )

    def addAggregateLabelOption(self):
        self.aggregate_label_type = tk.StringVar()
        self.options["Aggregate label type"] = tk.Frame(
            self.options_frame,
            bd=1,
            relief="raised"
        )
        tk.Label(
            self.options["Aggregate label type"],
            text="Aggregate label type",
        ).pack(anchor="w")
        label_types = [
            "None",
            "Peptide",
            "m/z",
            "rt",
            "dt"
        ]
        self.aggregate_label_type = tk.StringVar()
        self.aggregate_label_type.set(label_types[0])
        for label_type in label_types:
            tk.Radiobutton(
                self.options["Aggregate label type"],
                text=label_type,
                variable=self.aggregate_label_type,
                value=label_type
            ).pack(anchor="w")

    def addIonAxisOption(self):
        self.ion_axis_type = tk.StringVar()
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
            "m/z",
            "rt",
            "dt"
        ]
        self.ion_axis_type = tk.StringVar()
        self.ion_axis_type.set(label_types[0])
        for label_type in label_types:
            tk.Radiobutton(
                self.options["Ion axis type"],
                text=label_type,
                variable=self.ion_axis_type,
                value=label_type
            ).pack(anchor="w")

    def addFDROption(self):
        pass

    def addMinimumSignalOption(self):
        pass

    def addMaximumSignalOption(self):
        pass

    def initAggregatesFrame(self):
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
        fig = plt.Figure()
        self.aggregate_ax = fig.add_subplot(111)
        canvas = tkagg.FigureCanvasTkAgg(fig, self.aggregate_frame)
        canvas.get_tk_widget().pack(
            fill='both',
            expand=True
        )
        tkagg.NavigationToolbar2TkAgg(canvas, self.aggregate_frame)
        # df.plot(kind='Chart Type such as bar', legend=True, ax=ax)
        self.aggregate_ax.set_xlabel("RT")
        self.aggregate_ax.set_ylabel("DT")
        self.aggregate_ax.set_title("Aggregates")

    def initIonsFrame(self):
        self.ion_frame = tk.Frame(
            self.root_frame,
            bd=1,
            relief="raised"
        )
        self.ion_frame.pack(
            side="left",
            fill="y",
        )
        fig = plt.Figure()
        self.ion_ax = fig.add_subplot(111)
        canvas = tkagg.FigureCanvasTkAgg(fig, self.ion_frame)
        canvas.get_tk_widget().pack(
            fill='both',
            expand=True
        )
        tkagg.NavigationToolbar2TkAgg(canvas, self.ion_frame)
        self.ion_ax.set_xlabel("SAMPLE")
        # self.ion_ax.set_ylabel("DT")
        self.ion_ax.set_title("Ions")

    # def getVisibleNodes(self):
    #     dt_low, dt_high = plt.ylim()
    #     rt_low, rt_high = plt.xlim()
    #     selection = self.anchors["ION_COUNT"] >= self.MIN_REP_COUNT
    #     selection &= self.anchors["ION_COUNT"] <= self.MAX_REP_COUNT
    #     selection &= self.anchors["DT"] <= dt_high
    #     selection &= self.anchors["DT"] >= dt_low
    #     selection &= self.anchors["RT"] <= rt_high
    #     selection &= self.anchors["RT"] >= rt_low
    #     return np.flatnonzero(selection)

    def dummyCommand(self, *args, **kwargs):
        print(self.aggregate_label_type.get())


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
