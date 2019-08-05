#!venv/bin/python


import argparse
import src.ion_network


def argumentParser():
    parser = argparse.ArgumentParser(
        description="Analyze ion-networks through the commandline"
    )
    parser.add_argument(
        "-i",
        "--input_file_name",
        help="Parameter file (*.json)",
        required=True
    )
    parser.add_argument(
        "-a",
        "--actions_to_take",
        help="(A)nnotate, (B)rowse, and/or (C)reate ion-network",
        # default="CA",
    )
    # parser.add_argument(
    #     "-e",
    #     "--export_results",
    #     help="(A)ggregates, (E)dges, (I)ons, and/or (P)roteomic identifications",
    #     default="",
    # )
    return parser


if __name__ == "__main__":
    parser = argumentParser()
    input_file_name = parser.parse_args().input_file_name
    actions = parser.parse_args().actions_to_take.upper()
    # exports = parser.parse_args().export_results
    ion_network = src.ion_network.IonNetwork(input_file_name)
    if "C" in actions:
        ion_network.create()
    if "A" in actions:
        ion_network.annotate()
    if "B" in actions:
        ion_network.browse()
    # if "A" in exports:
    #     ion_network.exportAggregates()
    # if "I" in exports:
    #     ion_network.exportIons()
    # if "E" in exports:
    #     ion_network.exportEdges()
    # if "P" in exports:
    #     ion_network.exportAnnotations()
