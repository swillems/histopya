#!venv/bin/python


import tkinter as tk
import src.parameters
import src.io
import PySimpleGUI as sg
import subprocess
import numbers


def GUIWrapper():
    size = 25
    parameter_file_name = ""
    option_list = [
        ["Select parameters", setParameters],
        ["Create ion-network", createIonNetwork],
        ["Annotate ion-network", annotateIonNetwork],
        ["Browse ion-network", browseIonNetwork],
    ]
    option_dict = dict(option_list)
    layout = [
        [sg.Button(option, size=(size, 1))] for option, action in option_list
    ]
    window = sg.Window('Ion-network GUI').Layout(layout)
    while True:
        event, values = window.Read()
        if event is None:
            break
        if event in option_dict:
            if event == "Select parameters":
                parameter_file_name = setParameters(parameter_file_name)
            else:
                if parameter_file_name == "":
                    sg.PopupError("No parameter file selected")
                else:
                    option_dict[event](parameter_file_name)


def setParameters(parameter_file_name):
    layout = [
        [
            sg.Input(parameter_file_name, key="parameter_file_name"),
            sg.FileBrowse(),
            sg.Button("New"),
            sg.Button("Modify")
        ],
        [sg.Button("OK"), sg.Button("Cancel")],
    ]
    window = sg.Window('Select parameter file').Layout(layout)
    while True:
        event, values = window.Read()
        if (event is None) or (event == "Cancel"):
            break
        if event == "New":
            parameter_file_name = readParameters()
            window.Element("parameter_file_name").Update(parameter_file_name)
        if event == "Modify":
            parameter_file_name = values["parameter_file_name"]
            parameter_file_name = readParameters(parameter_file_name)
            window.Element("parameter_file_name").Update(parameter_file_name)
        if event == "OK":
            parameter_file_name = values["parameter_file_name"]
            try:
                parameters = src.parameters.importParameterDictFromJSON(
                    parameter_file_name,
                    save=True
                )
                break
            except src.parameters.ParameterError:
                sg.PopupError("Not a valid parameter file")
    window.Close()
    return parameter_file_name


def createIonNetwork(parameter_file_name):
    command = " ".join(
        ["bash run_cmd.sh", "-i", parameter_file_name, "-a", "C"]
    )
    executeScript(command, "Create ion-network")


def annotateIonNetwork(parameter_file_name):
    command = " ".join(
        ["bash run_cmd.sh", "-i", parameter_file_name, "-a", "A"]
    )
    executeScript(command, "Annotate ion-network")


def browseIonNetwork(parameter_file_name):
    command = " ".join(
        ["bash run_cmd.sh", "-i", parameter_file_name, "-a", "B"]
    )
    executeScript(command, "Browse ion-network")


def executeScript(command, title):
    layout = [
        [sg.Output(size=(100, 30))],
    ]
    window = sg.Window(title, layout)
    window.ReadNonBlocking()
    print("Loading external script...")
    window.ReadNonBlocking()
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, bufsize=0)
    while True:
        line = p.stdout.readline()
        if not line:
            break
        print(line.decode().rstrip())
        window.ReadNonBlocking()
    print("Finished external script.")
    print("Close this window to continue...")
    while True:
        (event, value) = window.Read()
        if (event is None) or (event == "Close"):
            window.Close()
            return


def readParameters(parameter_file_name=None):
    if parameter_file_name is None:
        parameters = src.parameters.getDefaultParameters()
    else:
        try:
            parameters = src.parameters.importParameterDictFromJSON(
                parameter_file_name,
                save=False
            )
        except src.parameters.ParameterError:
            sg.PopupError("Not a valid parameter file")
            return
    must_set = [
        "APEX_PATH",
        "OUTPUT_PATH",
        "DATABASE_FILE_NAME",
    ]
    apex_path = parameters["APEX_PATH"]
    output_path = parameters["OUTPUT_PATH"]
    database_file_name = parameters["DATABASE_FILE_NAME"]
    default_dict_layout = [
        [
            sg.Text("Apex folder"),
            sg.Input(apex_path, key="APEX_PATH"),
            sg.FolderBrowse(initial_folder="data"),
        ],
        [
            sg.Text("Output folder"),
            sg.Input(output_path, key="OUTPUT_PATH"),
            sg.FolderBrowse(initial_folder="projects"),
        ],
        [
            sg.Text("Database file"),
            sg.Input(database_file_name, key="DATABASE_FILE_NAME"),
            sg.FileBrowse(initial_folder="lib/databases"),
        ],
    ]
    default_dict_layout += [
        [
            sg.Text(key, size=(50, 1)),
            sg.InputText(value, key=key, enable_events=True)
        ] for key, value in parameters.items() if (
            value is None
        ) and not (
            key in must_set
        )
    ]
    extended_dict_layout = [
        [
            sg.Text(key, size=(50, 1)),
            sg.InputText(value, key=key, enable_events=True)
        ] for key, value in parameters.items() if (
            isinstance(value, str)
        ) and not (
            key in must_set
        )
    ]
    extended_dict_layout += [
        [
            sg.Text(key, size=(50, 1)),
            sg.Input(value, key=key, enable_events=True)
        ] for key, value in parameters.items() if (
            isinstance(value, numbers.Number)
        ) and not (
            key in must_set
        )
    ]
    buttons = [
        sg.Button('Save as'),
        sg.Button('Cancel'),
        sg.Button('Show more options'),
        sg.Button('Hide more options')
    ]
    if parameter_file_name:
        buttons = [sg.Button('Save')] + buttons
    # button_dict_layout = [buttons]
    default_window = sg.Column(
        default_dict_layout,
        scrollable=True,
        vertical_scroll_only=True,
        # size=(700, 100)
    )
    extended_window = sg.Column(
        extended_dict_layout,
        visible=parameter_file_name is not None,
        scrollable=True,
        vertical_scroll_only=True,
        key="EXTENDED_OPTIONS",
        size=(700, 300)
    )
    # window = sg.Window('Window Title', resizable=True).Layout(layout).Finalize()
    window = sg.Window(
        'Read parameters',
        # size=(800, 500),
        resizable=True
    ).Layout(
        [
            *default_dict_layout,
            # [default_window],
            buttons,
            [extended_window],
        ],
    ).Finalize()
    # c = sg.Column(
    #     [
    #         [default_window],
    #         [extended_window],
    #         [button_dict_layout]
    #     ],
    #     scrollable=True,
    #     vertical_scroll_only=True,
    #     # size=(700,300)
    # )
    # window = sg.Window(
    #     'Read parameters',
    #     size=(800, 800)
    # ).Layout(
    #     [
    #         [c],
    #     ]
    # )
    while True:
        # window.Size = window.Size
        # window.ReadNonBlocking()
        event, values = window.Read()
        if (event is None) or (event == "Cancel"):
            break
        if event == "Show more options":
            window.Element('EXTENDED_OPTIONS').Update(visible=True)
            # extended_window.visible = not extended_window.visible
        if event == "Hide more options":
            window.Element('EXTENDED_OPTIONS').Update(visible=False)
        if (event == "Save as") or (event == "Save"):
            for key, value in values.items():
                if key in parameters:
                    if parameters[key] is None:
                        parameters[key] = value
                    else:
                        try:
                            if isinstance(parameters[key], bool):
                                cast_value = type(parameters[key])(int(value))
                            else:
                                cast_value = type(parameters[key])(value)
                        except ValueError:
                            sg.PopupError("Wrong format for {}".format(value))
                            break
                        else:
                            parameters[key] = cast_value
            else:
                if (event == "Save as"):
                    parameter_file_name = tk.filedialog.asksaveasfilename(
                        filetypes=(("JSON files", "*.json"),),
                        defaultextension="*.json"
                    )
                if parameter_file_name:
                    src.io.saveJSON(
                        parameters,
                        parameter_file_name,
                    )
                    break
        window.Refresh()
        window.Refresh()
        window.Size = window.Size
    window.Close()
    return parameter_file_name


if __name__ == "__main__":
    GUIWrapper()
