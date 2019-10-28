# Ion-networks

## Overview

Ion-networks are a complete yet sparse representation of an MS experiment containing multiple samples. It was originally developed for LC-IMS-MS DIA data, but the IMS dimension can easily be replaced by e.g. a scanning quadrupole coordinate, since precursor m/z and drift time are correlated. The original paper describing ion-networks is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/726273v2) and all original data is available at [PRIDE PXD015318](www.ebi.ac.uk/pride/archive/projects/PXD015318).

The current implementation allows to

1. Create ion-networks
2. Annotate ion-networks (proteomics only)
3. Browse ion-networks

## Installation

The complete software suite was developed for Linux. As Windows 10 nowadays comes with a Linux subsystem (WSL), this platform is compatible as well through a slightly modified installation. All source code was written in Python 3.6 and tested on Ubuntu 18.04 and CentOS 7.6.1810. The installation includes the following:

* The full software suite to Create/Annotate/Browse ion-networks.
* Some standard indexed proteomic databases (1GB).
* [Percolator and Elude](https://github.com/percolator/percolator).
* A full demo proteomics benchmark dataset (3GB), containing an excerpt from [PXD001240](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001240) with 10 minutes from 5 samples in condition A (65% human, 15% yeast, 20% ecoli), and 5 samples in condition B (65% human, 30% yeast, 5% ecoli). This test dataset allows to recreate the main figure of our [paper](https://www.biorxiv.org/content/10.1101/726273v2), although some minor differences are present as the automatic estimation, automatic calibration and Percolator can only use the included data instead of the complete experiment.

### Windows 10

Windows users need to install a WSL for Ubuntu 18.04 by following [these steps](https://docs.microsoft.com/en-us/windows/wsl/install-win10). After the WSL has been installed and a user account has been created, copy-paste and execute the commands below line-per-line. Note that the second line takes about 30 minutes and several confirmations requiring you to manually type *y*. Once the *auto-reboot* screen prompts, you can safely click yes as this concerns only a WSL reboot without affecting Windows.

```
sudo -i
sudo apt update && sudo apt upgrade && sudo apt install python3.6 python3.6-venv python3.6-tk python3 python3-venv python3-tk libgomp1
mkdir percolator
wget --output-document percolator/percolator_3_02_01_ubuntu.tar.gz https://github.com/percolator/percolator/releases/download/rel-3-02-01/ubuntu64.tar.gz
tar xzvf percolator/percolator_3_02_01_ubuntu.tar.gz -C percolator
sudo dpkg -i percolator/elude-v3-02-linux-amd64.deb
sudo dpkg -i percolator/percolator-v3-02-linux-amd64.deb
rm -r percolator
wget https://github.com/swillems/histopya/archive/version_0.1.191028.tar.gz
tar xzvf version_0.1.191028.tar.gz
mv histopya-version_0.1.191028 histopya
cd histopya
bash install/install.sh
chmod -R 777 .
```

Next, close the WSL and download and install [MobaXterm v11.1](https://mobaxterm.mobatek.net/download-home-edition.html). Open it and press the *session* button on the top left. Select the rightmost tab *WSL* and set the Linux distribution to Ubuntu in the *Basic WSL settings* tab. Click the *Advanced WSL settings* tab and copy ```cd histopya; bash run_gui.sh``` to the *Execute the following commands at startup* window. Finally, click the *Bookmark settings* tab and change the *Session name* to e.g. *ion_network_gui*. Click the *Create a desktop shortcut to this session* button and select both options *Hide terminal on startup* and *Close MobaXterm on exit* before pressing *OK* in this popup. Confirm the session settings with *OK*. A pop-up with the GUI running should have appeared in your taskbar, allowing you to test the installation. To this end, click the *select parameters* button and subsequently *browse* in the *Select parameter file* pop-up. Navigate to the *install* folder and select the *test_lfq_parameters.json* file. Finally, perform a full test of all features with the buttons *Create/Annotate/Browse ion-network*. Each of them will show a pop-up that should complete with the final message *Finished external script*, indicating you can proceed with the next step. If the final *Browse ion-network* step correctly shows a second pop-up in which you can browse the test ion-network, you have successfully installed the complete software suite and can close MobaXTerm. For subsequent use, there is desktop-item that allows you to run the *ion_network_gui* with just a simple double-click.

NOTE: The WSL is a hidden subfolder located at e.g. ```C:\Users\sanwill\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc\LocalState\rootfs\home\histopya```. Vice-versa, your windows drives are accessible from within your WSL through the location */mnt/*. While you can copy-paste between your windows WSL and Windows 10, there might be some issues with permissions.

### Linux

First, run the following commands to install all dependencies with e.g. *apt* and *dpkg* on Ubuntu 18.04 (skip this if already satisfied):

```
sudo apt update && sudo apt upgrade && sudo apt install python3.6 python3.6-venv python3.6-tk python3 python3-venv python3-tk libgomp1
mkdir percolator
wget --output-document percolator/percolator_3_02_01_ubuntu.tar.gz https://github.com/percolator/percolator/releases/download/rel-3-02-01/ubuntu64.tar.gz
tar xzvf percolator/percolator_3_02_01_ubuntu.tar.gz -C percolator
sudo dpkg -i percolator/elude-v3-02-linux-amd64.deb
sudo dpkg -i percolator/percolator-v3-02-linux-amd64.deb
rm -r percolator
```

Next, install the source for for ion-networks itself:

```
wget https://github.com/swillems/histopya/archive/version_0.1.191028.tar.gz
tar xzvf version_0.1.191028.tar.gz
mv histopya-version_0.1.191028 histopya
cd histopya
bash install/install.sh
chmod -R 777 .
```

Finally, test if it is correctly installed by Creating/Annotating/Browsing the demo dataset with:

```
bash run_cmd.sh -i install/test_lfq_parameters.json -a CAB
```

## Input data

To create an ion-network, nothing more is needed than a single folder containing a csv per experiment with all its peak picked fragment ions. Herein, each ion needs to have the following  attributes:

1. m/z
2. DT
3. RT
4. intensity
5. m/z error (ppm)
6. DT error
7. RT error

In case of Waters data (e.g. HDMSE from a Synapt G2-Si or SONAR from a Vion), peak picked data is most easily obtainable by Waters' Apex3D software as provided by the commercial software packages MassLynx or Progenesis QI (Nonlinear Dynamics) with a command similar to:

```
Apex3D64.exe -pRawDirName sample.raw -outputDirName peak_picked_sample_folder -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0
```

Alternatively, the freely available [IMTBX and Grppr](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5826643) can be used for Waters' peak picking but this has not been tested for compatibility yet.

In case of proteomics data, a database is needed for fragment annotation. Human, Yeast, Ecoli and benchmark databases are included in the ```lib/databases/``` folder. Custum databases can be created with some simple coding (```src/peptides.py``` script), but no GUI has been provided yet.

## Usage

Analysis of ion-networks can be done both by the command-line and with a graphical user interface (GUI). To use the command-line, a *.json* file needs to be created (either manual or with the GUI) which should at least contain the following parameters:

```
{
    "APEX_PATH": "data/test_lfq/",
    "OUTPUT_PATH": "projects/test_lfq/",
    "DATABASE_FILE_NAME": "lib/databases/lfq_benchmark.hdf5"
}
```

With this parameter file, three different actions can (simultaneously) be taken with the command-line flag ```-a``` or by simple selection in the GUI.

### (C)reate ion-network

An ion-network can be created with the following command:

```
bash run_cmd.sh -i install/test_lfq_parameters.json -a C
```

While the default parameters generally suffice, the following parameters can have a significant effect:

* ```"NEIGHBOR_THRESHOLD": 0.997``` Proxy for robustness. Higher values have less deconvoluting power, but higher robustness. Set closer to 1 for experiments with many samples.
* ```"SIGNAL_COUNT_THRESHOLD": 2``` Minimum reproducibility for an aggregate, set higher for more stringent denoising.

### (A)nnotate (proteomics) ion-network

Currently, only proteomics data can be annotated with a database search approach. Human, Yeast, Ecoli and benchmark databases are available, but more can be created by the user with some custom coding.

To annotate an ion-network, use the following command:

```
bash run_cmd.sh -i install/test_lfq_parameters.json -a A
```

Significant parameters that can change annotation results include:

* ```"USE_RT_IN_PERCOLATOR": false``` Use Elude to improve Percolator. Improves annotation results, but can take a lot of time.
* ```"REQUIRE_PHYSICAL_PRECURSOR": false``` If set to true, an unfragmented precursor signal is required to consistently co-elute before an annotation can be considered.
* ```"IDENTIFICATION_PPM": 20``` Only fragments within this error range can be annotated.

### (B)rowse ion-network

An (un)annotated ion-network can be visualized with:

```
bash run_cmd.sh -i install/test_lfq_parameters.json -a B
```

Note that most ion-networks are quite vast and that it is not always wise to try and visualize everything. This is particularly true for the button *Show edges* and the slider *Minimum signal*, which generally should only be changed when zoomed in sufficiently.

Aggregates can be selected by double clicking and multiple aggregates can be selected by ctrl-double-click. Note that there is a minor bug where double clicking can force zooming when this option is still active, so we advise to first zoom to a region of interest, then turn off the zoom function and only then to start selecting aggregates of interest.
