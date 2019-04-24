#!src/venv/bin/python


from __future__ import print_function
import Tkinter as tk
import tkFileDialog
import tkMessageBox
from src.conversion import to_json as json_converter
from src.prioritization.prioritizer import Prioritization
from src.gui.base_gui import GenericGUI
from src.gui.tabs.curation import CurationFrame
from src.gui.tabs.prioritization import PrioritizationFrame
from src.gui.tabs.annotation import AnnotationFrame
from src.gui.tabs.sample_preparation import SamplePreparationFrame
from src.gui.tabs.ptms import PTMFrame
from src.gui.tabs.amino_acids import AminoAcidFrame
from src.gui.tabs.proteins import ProteinFrame
from collections import OrderedDict
# from time import asctime


DEFAULT_PARAMETERS = "data/sample/parameters.json"
# DEFAULT_PARAMETERS = "lib/default/parameters.json"


class PTMGUI(GenericGUI):

    def __init__(self, title):
        GenericGUI.__init__(self, title)
        self.parameters = {}
        self.loadParameters(
            use_default=True,
            ask_confirmation=False
        )
        self.setEasterEggs()

    def initSubMenus(self):
        sub_menus = OrderedDict()
        sub_menus["File"] = [
            ("New", "button", lambda: self.loadParameters(use_default=True)),
            ("Load", "button", self.loadParameters),
            ("Save", "button", self.saveParameters),
            ("Quit", "button", self.quitMenuCommand),
        ]
        sub_menus["Log"] = [
            ("Save Display To File", "checkbox", self.selectLogMenuCommand),
            ("Archive Full Analysis", "button", self.archiveMenuCommand),
        ]
        sub_menus["View"] = [
            ("Show Display", "checkbox", self.showDisplayMenuCommand),
            ("Show Advanced Parameter Tabs", "checkbox", self.showAdvancedParameterTabsMenuCommand),
        ]
        sub_menus["Help"] = [
            ("About", "button", lambda: print("TODO")),
            ("Citation", "button", lambda: print("TODO")),
        ]
        for menu_name, items in sub_menus.iteritems():
            sub_menu = self.newSubMenu(label=menu_name)
            for command_name, command_type, command in items:
                if command_type == "checkbox":
                    active = tk.BooleanVar()
                    sub_menu.add_checkbutton(
                        label=command_name,
                        variable=active,
                        command=lambda active=active, command=command: command(active)
                    )
                elif command_type == "button":
                    sub_menu.add_command(
                        label=command_name,
                        command=command
                    )

    def selectLogMenuCommand(self, active):
        if not active.get():
            self.write = self.writeTo(log_file=None)
        else:
            log_file_name = tkFileDialog.asksaveasfilename(
                title='Save log file as'
            )
            if log_file_name:
                self.write = self.writeTo(log_file=log_file_name)
            else:
                self.write = self.writeTo(log_file=None)
                active.set(False)

    def archiveMenuCommand(self):
        archive_file_name = tkFileDialog.asksaveasfilename(
            title='Save archive file as'
        )
        if archive_file_name:
            self.actOnFunction(
                self.prioritization.save,
                save_as_file_name=archive_file_name,
            )

    def showDisplayMenuCommand(self, active):
        self.toggleActiveDisplay(active.get())

    def quitMenuCommand(self):
        if tkMessageBox.askokcancel(
            "Quit?",
            "All unsaved data will be lost!"
        ):
            self.quit()

    def showAdvancedParameterTabsMenuCommand(self, active=None):
        if active is not None and active.get():
            self.parameter_types["advanced"].pack(
                side="top",
                fill="x",
                expand=1
            )
        else:
            self.parameter_types["advanced"].pack_forget()

    def loadParameters(
        self,
        use_default=False,
        ask_confirmation=True
    ):
        if ask_confirmation and not tkMessageBox.askokcancel(
            "Load parameters?",
            "All unsaved data will be lost!"
        ):
            return
        if use_default:
            parameter_file_name = DEFAULT_PARAMETERS
        else:
            parameter_file_name = tkFileDialog.askopenfilename(
                title='Load parameter file'
            )
        if parameter_file_name:
            self.prioritization = Prioritization()
            self.prioritization.loadParameters(parameter_file_name)
            self.importParameters(parameter_file_name)
            try:
                self.initDicts()
            except AttributeError:
                self.initSubMenus()
                self.initMainFrame()
                self.initDicts()
                # TODO
                # log_file = "logs/{}.txt".format(
                #     "_".join(
                #             asctime().split()
                #         )
                #     )

    def saveParameters(self):
        parameter_file_name = tkFileDialog.asksaveasfilename(
            title='Save parameter file as'
        )
        if parameter_file_name:
            self.actOnFunction(
                self.prioritization.save,
                save_as_file_name=parameter_file_name,
                item_to_save=self.prioritization.PARAMETERS
            )

    def importParameters(self, parameter_file_name):
        # TODO
        def update(
            parameters,
            parameter,
            variable,
        ):
            prev_var = parameters[parameter]
            try:
                parameters[parameter] = variable.get()
            except ValueError, e:
                if str(e) == "invalid literal for int() with base 10: ''":
                    variable.set(0)
                else:
                    variable.set(prev_var)
        for section, section_parameters in self.prioritization.PARAMETERS.iteritems():
            if section not in self.parameters:
                self.parameters[section] = {}
            for parameter, parameter_value in section_parameters.iteritems():
                if parameter not in self.parameters[section]:
                    if isinstance(parameter_value, unicode):
                        self.parameters[section][parameter] = tk.StringVar()
                    elif isinstance(parameter_value, int):
                        self.parameters[section][parameter] = tk.IntVar()
                    elif isinstance(parameter_value, bool):
                        self.parameters[section][parameter] = tk.BooleanVar()
                    elif isinstance(parameter_value, float):
                        self.parameters[section][parameter] = tk.DoubleVar()
                    else:
                        # TODO
                        self.parameters[section][parameter] = None
                        continue
                    variable = self.parameters[section][parameter]
                    self.parameters[section][parameter].trace(
                        "w",
                        lambda a, b, c, variable=variable, parameters=self.prioritization.PARAMETERS[section], parameter=parameter: update(
                            parameters=parameters,
                            parameter=parameter,
                            variable=variable
                        )
                    )
                # TODO
                if self.parameters[section][parameter] is None:
                    continue
                self.parameters[section][parameter].set(parameter_value)

    def initDicts(self):
        ptms = dict(self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"])
        self.reloadVarFixPTMs()
        self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"] = ptms
        for ptm in ptms["fixed"].itervalues():
            ptm_index = self.varfix_ptm_frame.children["all_ptms"].children['listbox'].get(0, "end").index(ptm)
            self.varfix_ptm_frame.children["all_ptms"].children['listbox'].delete(ptm_index)
            self.varfix_ptm_frame.children["fix_ptms"].children['listbox'].insert("end", ptm)
        for aa, ptms in ptms["variable"].iteritems():
            for ptm in ptms:
                ptm_index = self.varfix_ptm_frame.children["all_ptms"].children['listbox'].get(0, "end").index(ptm)
                self.varfix_ptm_frame.children["all_ptms"].children['listbox'].delete(ptm_index)
                self.varfix_ptm_frame.children["var_ptms"].children['listbox'].insert("end", ptm)
        # self.updateAADict()
        # self.updatePTMDict()
        # self.updateProteinDict()

    def getVariable(self, variable_name_tuple):
        current = self.parameters
        for name in variable_name_tuple[:-1]:
            current = current[name]
        return current[variable_name_tuple[-1]]

    def initMainFrame(self):
        self.initSelectionFrame()
        self.initParameterFrames()
        self.current_parameter_frame = tk.Frame(
            self.main_frame,
            bd=1,
            relief="raised"
        )
        self.current_parameter_frame.pack(fill="both", expand=1)
        self.showAdvancedParameterTabsMenuCommand()

    def initSelectionFrame(self):
        self.parameter_selection_frame = tk.Frame(
            self.main_frame,
            bd=1,
            relief="raised"
        )
        self.parameter_selection_frame.pack(
            side="top",
            fill="x"
        )
        self.parameter_frames = {}
        self.parameter_types = {
            "advanced": tk.Frame(self.parameter_selection_frame),
            "basic": tk.Frame(self.parameter_selection_frame),
        }
        self.parameter_types["basic"].pack(
            side="top",
            fill="x"
        )

    def initParameterFrames(self):
        for parameter_group, frame, parameter_type in [
            ("Amino Acids", AminoAcidFrame, "advanced"),
            ("PTMs", PTMFrame, "advanced"),
            ("Proteins", ProteinFrame, "advanced"),
            ("Sample Preparation", SamplePreparationFrame, "basic"),
            ("Mass Spectrometry", AnnotationFrame, "basic"),
            ("Prioritization", PrioritizationFrame, "basic"),
            ("Curation", CurationFrame, "basic"),
        ]:
            self.parameter_frames[parameter_group] = frame(
                self,
                master=self.main_frame,
                bd=2,
                relief="raised"
            )
            tk.Button(
                self.parameter_types[parameter_type],
                text=parameter_group,
                command=lambda parameter_group=parameter_group: self.switchParameterFrame(parameter_group)
            ).pack(side="left", fill="x", expand=1)

    def switchParameterFrame(self, parameter_group):
        self.current_parameter_frame.pack_forget()
        self.current_parameter_frame = self.parameter_frames[parameter_group]
        self.current_parameter_frame.pack(fill="both", expand=1)

    def createOrUpdatePTM(self):
        name = self.current_manual_ptm_fields["Name"].get()
        if name in self.prioritization.PTM_DICT and not tkMessageBox.askokcancel(
            "Existing PTM",
            "Are you sure you want to update this PTM definition?"
        ):
            return
        self.prioritization.PTM_DICT[name] = {
            "CF": self.current_manual_ptm_fields["Chemical Formula"].get(),
            "DR": self.current_manual_ptm_fields["References"].get(),
            "KW": self.current_manual_ptm_fields["Keywords"].get(),
            "MM": self.current_manual_ptm_fields["Monoisotopic Mass"].get(),
            "TG": self.current_manual_ptm_fields["Target"].get(),
        }
        self.updatePTMDict(from_file=False)

    # def deletePTM(self):
    #     name = self.current_manual_ptm_fields["Name"].get()
    #     if name not in self.prioritization.PTM_DICT:
    #         return
    #     if not tkMessageBox.askokcancel(
    #         "Existing PTM",
    #         "Are you sure you want to delete this PTM definition?"
    #     ):
    #         return
    #     del self.prioritization.PTM_DICT[name]
    #     self.current_manual_ptm_fields["Chemical Formula"].set(None)
    #     self.current_manual_ptm_fields["References"].set(None)
    #     self.current_manual_ptm_fields["Keywords"].set(None)
    #     self.current_manual_ptm_fields["Monoisotopic Mass"].set(None)
    #     self.current_manual_ptm_fields["Target"].set(None)
    #     self.updatePTMDict(from_file=False)

    def manualPTMFrame(self, parent_frame):
        frame = tk.Frame(parent_frame)
        self.current_manual_ptm_fields = OrderedDict()
        self.current_manual_ptm_fields["Name"] = tk.StringVar()
        self.current_manual_ptm_fields["Chemical Formula"] = tk.StringVar()
        self.current_manual_ptm_fields["References"] = tk.StringVar()
        self.current_manual_ptm_fields["Keywords"] = tk.StringVar()
        self.current_manual_ptm_fields["Monoisotopic Mass"] = tk.DoubleVar()
        self.current_manual_ptm_fields["Target"] = tk.StringVar()
        for row, (field, variable) in enumerate(self.current_manual_ptm_fields.items()):
            tk.Label(
                frame,
                text=field
            ).grid(row=row, column=0, sticky='w')
            tk.Entry(
                frame,
                textvariable=variable,
            ).grid(row=row, column=1, sticky='w')
        tk.Button(
            frame,
            text="Update",
            command=self.createOrUpdatePTM
        ).grid(row=len(self.current_manual_ptm_fields), column=0)
        # tk.Button(
        #     frame,
        #     text="Delete",
        #     command=self.deletePTM
        # ).grid(row=len(self.current_manual_ptm_fields), column=1)
        return frame

    def updatePTMDict(self, from_file=True):
        self.toggleActiveStatus(True)
        if self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTM_FILE_NAME"]:
            self.ptm_scroller.children["listbox"].delete(0, "end")
            if from_file:
                self.prioritization.loadPTMDict()
            for ptm in sorted(self.prioritization.PTM_DICT.iterkeys()):
                self.ptm_scroller.children["listbox"].insert("end", ptm)
        # self.ptm_scroller.children["listbox"].itemconfig(3, {'fg': 'blue'})
        self.toggleActiveStatus(False)

    def updateProteinDict(self):
        self.toggleActiveStatus(True)
        if self.prioritization.PARAMETERS["SAMPLE_PREP"]["PROTEIN_FILE_NAME"]:
            self.protein_scroller.children["listbox"].delete(0, "end")
            self.loadProteinDict()
            for protein in sorted(self.prioritization.PROTEIN_DICT.iterkeys()):
                self.protein_scroller.children["listbox"].insert("end", protein)
        self.toggleActiveStatus(False)

    def ptmVarFixFrame(self, parent_frame):
        def createSelectionFrame(parent_frame, name):
            selection_frame = tk.Frame(parent_frame, name=name)
            listbox = tk.Listbox(selection_frame, name="listbox")
            listbox.pack(side="left", fill="both", expand=True)
            scrollbar = tk.Scrollbar(selection_frame, orient="vertical")
            scrollbar.config(command=listbox.yview)
            scrollbar.pack(side="right", fill="y")
            listbox.config(yscrollcommand=scrollbar.set)
            return selection_frame
        def moveVarFixPTM(from_field, to_field, parent_frame):
            from_list = parent_frame.children[from_field].children['listbox']
            to_list = parent_frame.children[to_field].children['listbox']
            selected_index = from_list.curselection()
            selected = from_list.get(selected_index)
            from_list.delete(selected_index)
            to_list.insert("end", selected)
            aa = self.prioritization.PTM_DICT[selected]["TG"]
            if from_field == "all_ptms":
                if to_field == "fix_ptms":
                    self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"]["fixed"][aa] = selected
                elif to_field == "var_ptms":
                    if aa not in self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"]["variable"]:
                        self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"]["variable"][aa] = []
                    self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"]["variable"][aa].append(selected)
            elif to_field == "all_ptms":
                if from_field == "fix_ptms":
                    del self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"]["fixed"][aa]
                elif from_field == "var_ptms":
                    self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"]["variable"][aa].remove(selected)
                    if len(self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"]["variable"][aa]) == 0:
                        del self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"]["variable"][aa]
        frame = tk.Frame(parent_frame)
        createSelectionFrame(frame, name="all_ptms").grid(
            row=0, column=0, rowspan=6, sticky="nsew"
        )
        tk.Button(
            frame,
            text="Reload",
            name="reload_ptms",
            command=self.reloadVarFixPTMs
        ).grid(row=6, column=0, columnspan=3, sticky="ew")
        for full_name, short_name, row_offset in [
            ("Variable PTMs", "var_ptms", 0),
            ("Fixed PTMs", "fix_ptms", 3)
        ]:
            createSelectionFrame(frame, name=short_name).grid(
                row=row_offset, column=2, rowspan=3, sticky="nsew"
            )
            tk.Label(
                frame,
                text=full_name
            ).grid(row=row_offset + 0, column=1, sticky="new")
            tk.Button(
                frame,
                text="Add",
                command=lambda from_field="all_ptms", to_field=short_name: moveVarFixPTM(from_field, to_field, frame),
            ).grid(row=row_offset + 1, column=1, sticky="new")
            tk.Button(
                frame,
                text="Remove",
                command=lambda from_field=short_name, to_field="all_ptms": moveVarFixPTM(from_field, to_field, frame),
            ).grid(row=row_offset + 2, column=1, sticky="new")
        return frame

    def reloadVarFixPTMs(self):
        self.actOnFunction(self.loadPTMDict)
        selected_items = self.prioritization.PTM_DICT
        self.varfix_ptm_frame.children['all_ptms'].children['listbox'].delete(0, "end")
        self.varfix_ptm_frame.children['var_ptms'].children['listbox'].delete(0, "end")
        self.varfix_ptm_frame.children['fix_ptms'].children['listbox'].delete(0, "end")
        for ptm in sorted(selected_items):
            self.varfix_ptm_frame.children['all_ptms'].children['listbox'].insert("end", ptm)
        self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTMS"] = {"variable": {}, "fixed": {}}

    def setParameter(
        self,
        parameter_frame,
        variable_name_tuple,
        row,
    ):
        tk.Label(
            parameter_frame,
            text=variable_name_tuple[-1]
        ).grid(row=row, column=0, sticky='w')
        variable = self.getVariable(variable_name_tuple)
        tk.Entry(
            parameter_frame,
            textvariable=variable,
        ).grid(row=row, column=1, sticky='w')

    def setBrowsableParameter(
        self,
        parameter_frame,
        variable_name_tuple,
        row,
        dialog_title,
        browseOnly=True,
        saveButton=False,
        loadFileFunction=None,
        saveFileFunction=None,
    ):
        def browseFile():
            self.toggleActiveStatus(True)
            if saveButton:
                browseFileName = tkFileDialog.asksaveasfilename
            else:
                browseFileName = tkFileDialog.askopenfilename
            file_name = browseFileName(
                parent=parameter_frame,
                title=dialog_title
            )
            if file_name:
                variable.set(file_name)
                # if tkMessageBox.askyesno(
                #     "Load file",
                #     "Do you want to immediately load this file?"
                # ):
                #     self.actOnFunction(loadFileFunction)
            self.toggleActiveStatus(False)
        # def loadFile():
        #     self.actOnFunction(loadFileFunction)
        # def saveFile():
        #     self.actOnFunction(saveFileFunction)
        tk.Label(
            parameter_frame,
            text=variable_name_tuple[-1]
        ).grid(row=row, column=0, sticky='w')
        variable = self.getVariable(variable_name_tuple)
        tk.Entry(
            parameter_frame,
            textvariable=variable,
        ).grid(row=row, column=1, sticky='w')
        options = [("Browse", browseFile)]
        if not browseOnly:
            options.append(("Load", loadFile))
            options.append(("Save", saveFile))
        for column, (text, command) in enumerate(options):
            tk.Button(
                parameter_frame,
                text=text,
                command=command
            ).grid(row=row, column=column + 2, sticky='w')

    def loadAADict(self):
        file_name = self.prioritization.PARAMETERS["SAMPLE_PREP"]["AA_FILE_NAME"]
        if not file_name.endswith(".json"):
            parameter_file_name = None
            if tkMessageBox.askyesno(
                "Save file",
                "The file is not in JSON format."
                "Do you want to save it as json for quick access in the future?"
            ):
                parameter_file_name = tkFileDialog.asksaveasfilename(
                    title='Save file in .json format'
                )
            self.prioritization.AA_DICT = json_converter.loadAADict(file_name)
            if parameter_file_name:
                self.actOnFunction(
                    self.prioritization.save,
                    save_as_file_name=parameter_file_name,
                    item_to_save=self.prioritization.AA_DICT
                )
        else:
            self.prioritization.loadAADict()

    def loadPTMDict(self):
        file_name = self.prioritization.PARAMETERS["SAMPLE_PREP"]["PTM_FILE_NAME"]
        if not file_name.endswith(".json"):
            parameter_file_name = None
            if tkMessageBox.askyesno(
                "Save file",
                "The file is not in JSON format."
                "Do you want to save it as json for quick access in the future?"
            ):
                parameter_file_name = tkFileDialog.asksaveasfilename(
                    title='Save file in .json format'
                )
            self.prioritization.PTM_DICT = json_converter.loadPTMDict(file_name)
            if parameter_file_name:
                self.actOnFunction(
                    self.prioritization.save,
                    save_as_file_name=parameter_file_name,
                    item_to_save=self.prioritization.PTM_DICT
                )
        else:
            self.prioritization.loadPTMDict()

    def loadProteinDict(self):
        file_name = self.prioritization.PARAMETERS["SAMPLE_PREP"]["PROTEIN_FILE_NAME"]
        if not file_name.endswith(".json"):
            parameter_file_name = None
            if tkMessageBox.askyesno(
                "Save file",
                "The file is not in JSON format."
                "Do you want to save it as json for quick access in the future?"
            ):
                parameter_file_name = tkFileDialog.asksaveasfilename(
                    title='Save file in .json format'
                )
            self.prioritization.PROTEIN_DICT = json_converter.loadProteinDict(file_name)
            if parameter_file_name:
                self.actOnFunction(
                    self.prioritization.save,
                    save_as_file_name=parameter_file_name,
                    item_to_save=self.prioritization.PROTEIN_DICT
                )
        else:
            self.prioritization.loadProteinDict()

    def loadPrecursorList(self):
        file_name = self.prioritization.PARAMETERS["MS"]["PRECURSOR_FILE_NAME"]
        if not file_name.endswith(".json"):
            parameter_file_name = None
            if tkMessageBox.askyesno(
                "Save file",
                "The file is not in JSON format."
                "Do you want to save it as json for quick access in the future?"
            ):
                parameter_file_name = tkFileDialog.asksaveasfilename(
                    title='Save file in .json format'
                )
            self.prioritization.PRECURSOR_LIST = json_converter.loadPrecursorList(file_name)
            if parameter_file_name:
                self.actOnFunction(
                    self.prioritization.save,
                    save_as_file_name=parameter_file_name,
                    item_to_save=self.prioritization.PRECURSOR_LIST
                )
        else:
            self.prioritization.loadPrecursorList()

    def splitScrollFrame(
        self,
        parent_frame
    ):
        scroll_frame = tk.Frame(parent_frame)
        scroll_frame.pack(side="left", fill="both")
        listbox = tk.Listbox(scroll_frame, name="listbox")
        listbox.pack(side="left", fill="both")
        scrollbar = tk.Scrollbar(scroll_frame, orient="vertical")
        scrollbar.config(command=listbox.yview)
        scrollbar.pack(side="right", fill="y")
        listbox.config(yscrollcommand=scrollbar.set)
        modification_frame = tk.Frame(parent_frame)
        modification_frame.pack(side="left", fill="both", expand=1)
        return scroll_frame, modification_frame


if __name__ == "__main__":
    gui = PTMGUI("A Priori(tization) Tool")
    gui.start()





class BaseGUI(object):

    def __init__(self, title="HistoPyA"):
        self.__root = tk.Tk()
        self.__root.winfo_toplevel().title(title)
        self.__initMenuBar()
        self.__initStatusBar()
        frame = tk.Frame(self.__root, bd=1, relief="raised")
        frame.pack(fill="both", expand=1)
        self.__initMainFrame(frame)
        self.__initDisplay(frame)

    def __initMenuBar(self):
        self.__menu_bar = tk.Menu(self.__root)
        self.__root.config(menu=self.__menu_bar)

    def __initStatusBar(self):
        self.__active_status = False
        self.__display_status = tk.StringVar()
        self.toggleActiveStatus(active=False)
        tk.Label(
            self.__root,
            textvariable=self.__display_status,
            relief="sunken",
            anchor="w",
        ).pack(
            side="bottom",
            fill="x"
        )

    def __initMainFrame(self, frame):
        self.main_frame = tk.Frame(
            frame,
        )
        self.main_frame.pack(
            side="left",
            fill="both",
            expand=1
        )

    def newSubMenu(self, label):
        sub_menu = tk.Menu(self.__menu_bar, tearoff=0, name=label.lower())
        self.__menu_bar.add_cascade(
            label=label,
            menu=sub_menu
        )
        return sub_menu

    def toggleActiveStatus(self, active=None):
        if active is not None:
            self.__active_status = active
        else:
            self.__active_status = not self.__active_status
        if self.__active_status:
            self.__display_status.set("Busy")
        else:
            self.__display_status.set("Idle")

    def actOnFunction(self, function, *args, **kwargs):
        self.toggleActiveStatus(active=True)
        self.__performFunction(function, *args, **kwargs)
        self.toggleActiveStatus(active=False)

    def actOnFunctions(self, *functions):
        self.toggleActiveStatus(active=True)
        for function in functions:
            self.__performFunction(*function)
        self.toggleActiveStatus(active=False)

    def __performFunction(self, function, *args, **kwargs):
        function(*args, **kwargs)

    def start(self):
        self.write = self.writeTo(log_file=None)
        sys.stderr = self
        sys.stdout = self
        self.__root.mainloop()

    def quit(self):
        self.__root.quit()
