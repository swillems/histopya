#!venv/bin/python


# import tkinter as tk
import src.parameters
import src.io
import src.ions
import src.aggregates
import src.peptides
import PySimpleGUI as sg
import subprocess


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
            except:
                sg.PopupError("Not a valid parameter file")
    window.Close()
    return parameter_file_name


def createIonNetwork(parameter_file_name):
    command = " ".join(
        ["./run_cmd.sh", "-i", parameter_file_name, "-a", "C"]
    )
    executeScript(command, "Create ion-network")


def annotateIonNetwork(parameter_file_name):
    command = " ".join(
        ["./run_cmd.sh", "-i", parameter_file_name, "-a", "A"]
    )
    executeScript(command, "Annotate ion-network")


def browseIonNetwork(parameter_file_name):
    command = " ".join(
        ["./run_cmd.sh", "-i", parameter_file_name, "-a", "B"]
    )
    executeScript(command, "Browse ion-network")


def executeScript(command, title):
    layout = [
        [sg.Output(size=(100, 30))],
        # [sg.Button("Close")]
    ]
    # print = sg.Print
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
                save=True
            )
        except:
            sg.PopupError("Not a valid parameter file")
            return
    dict_layout = [
        [
            sg.Text(key, size=(50, 1)),
            sg.InputText(value, key=key, enable_events=True)
        ] for key, value in parameters.items() if (
            isinstance(value, str)
        )
    ]
    dict_layout += [
        [
            sg.Text(key, size=(50, 1)),
            sg.Input(value, key=key, enable_events=True)
        ] for key, value in parameters.items() if (
            isinstance(value, int)
        )
    ]
    buttons = [
        # sg.SaveAs(
        #     'Save',
        #     file_types=(("JSON files", "*.json"),),
        #     enable_events=True
        # ),
        sg.Button('Save as'),
        sg.Button('Cancel')
    ]
    if parameter_file_name:
        buttons = [sg.Button('Save')] + buttons
    dict_layout += [buttons]
    c = sg.Column(
        dict_layout,
        scrollable=True,
        vertical_scroll_only=True,
        # size=(700,300)
    )
    window = sg.Window(
        'Read parameters',
        size=(700, 800)
    ).Layout(
        [
            [c],
        ]
    )
    while True:
        event, values = window.Read()
        if (event is None) or (event == "Cancel"):
            break
        # if event in parameters:
        #     try:
        #         cast_event = type(parameters[event])(values[event])
        #         window.Element(event).Update(cast_event)
        #     except ValueError:
        #         window.Element(event).Update(parameters[event])
        if (event == "Save as") or (event == "Save"):
            for value in values:
                if value not in parameters:
                    continue
                try:
                    cast_value = type(parameters[value])(values[value])
                    parameters[value] = cast_value
                except ValueError:
                    sg.PopupError("Wrong format for {}".format(value))
                    break
            else:
                if (event == "Save as"):
                    parameter_file_name = tk.filedialog.asksaveasfilename(
                        filetypes=(("JSON files", "*.json"),),
                        defaultextension="*.json"
                    )
                    # parameter_file_name = sg.PopupGetFile(
                    #     "Save parameters as",
                    #     save_as=True,
                    #     file_types=(("JSON files", "*.json"),)
                    # )
                if parameter_file_name:
                    src.io.saveJSON(
                        parameters,
                        parameter_file_name,
                    )
                    break
    window.Close()
    return parameter_file_name


if __name__ == "__main__":
    GUIWrapper()
