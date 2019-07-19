# Ion-networks

## Overview

Ion-networks are a complete yet sparse representation of an MS experiment containing multiple samples. It was originally developed for LC-IMS-MS DIA data, but the IMS dimension can easily be replaced by e.g. a scanning quadrupole coordinate, since precursor m/z and drift time are correlated. The original paper describing ion-networks is available in the docs folder.

The current implementation allows to

1. Create ion-networks
2. Annotate ion-networks (proteomics only)
3. Browse ion-networks

## Installation

The complete software suite was developed for Linux. As Windows 10 nowadays comes with a Linux subsystem (WSL), this platform is compatible as well through a slightly modified installation. All source code was written in Python 3.6 and tested on Ubuntu 18.04 and CentOS 7.6.1810.

NOTE: The installation includes a full demo [proteomics benchmark dataset](URL_TODO) of 3GB, some standard proteomic databases, and [percolator](https://github.com/percolator/percolator) from external sources.

### Linux

First, run the following commands to install all dependancies on e.g. Ubuntu 18.04 (skip this if already satisfied):

```
sudo apt update && sudo apt upgrade && sudo apt install python3.6 python3.6-venv python3.6-tk python3 python3-venv python3-tk libgomp1
wget --output-document percolator/percolator_3_02_01_ubuntu.tar.gz https://github.com/percolator/percolator/releases/download/rel-3-02-01/ubuntu64.tar.gz
tar xzvf percolator/percolator_3_02_01_ubuntu.tar.gz -C install
sudo dpkg -i percolator/elude-v3-02-linux-amd64.deb
sudo dpkg -i percolator/percolator-v3-02-linux-amd64.deb
```

Finally, install the source for for ion-networks itself:

```
git clone https://github.com/swillems/histopya
cd histopya/install
bash install/install.sh
```

### Windows 10

Windows users need to install a WSL for Ubuntu 18.04 by following [these steps](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

After the WSL has been installed and a user account has been created, copy-paste the following commands (takes about an hour and some password confirmations. When prompted, click yes on the *auto-reboot* screen, as this is only a WSL reboot and not a Windows reboot):

```
sudo apt update && sudo apt upgrade && sudo apt install python3.6 python3.6-venv python3.6-tk python3 python3-venv python3-tk libgomp1
wget --output-document percolator/percolator_3_02_01_ubuntu.tar.gz https://github.com/percolator/percolator/releases/download/rel-3-02-01/ubuntu64.tar.gz
tar xzvf percolator/percolator_3_02_01_ubuntu.tar.gz -C install
sudo dpkg -i percolator/elude-v3-02-linux-amd64.deb
sudo dpkg -i percolator/percolator-v3-02-linux-amd64.deb
git clone https://github.com/swillems/histopya
cd histopya/install
bash install/install.sh
chmod -R 777 .
```

Next, download and install [MobaXterm v11.1](https://mobaxterm.mobatek.net/download-home-edition.html). Open it and press the *session* button on the top left. Select the rightmost tab *WSL* and set the Linux distribution to Ubuntu in the *Basic WSL settings* tab. Click the *Advanced WSL settings* tab and copy ```cd histopya; bash run_gui.sh``` to the *Execute the following commands at startup* window. Finally, click the *Bookmark settings* tab and change the *Session name* to e.g. *ion_network_gui*. Click the *Create a desktop shortcut to this session* button and select both options *Hide terminal on startup* and *Close MobaXterm on exit* before pressing *OK* in this popup. Confirm the session settings with *OK* and close MobaXterm.

NOTE: The WSL is hidden subfolder located at e.g. ```C:\Users\sanwill\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc\LocalState\rootfs\home\histopya```. Vice-versa, your windows drives are accessible from within your WSL through the location */mnt/*. While you can copy-paste data to and from here, there might be some issues with permissions.

## Input data

To create an ion-network, nothing more is needed than a single folder containing a csv per experiment with all its peak picked fragment ions. Herein, each ion need to have the following  attributes:

1. m/z
2. DT
3. RT
4. intensity
5. m/z error (ppm)
6. DT error
7. RT error

In case of Waters data (e.g. Synapt G2-Si or Vion), peak picked data is most easily obtainable by their Apex3D software as provided in MassLynx or in Progenesis QI (Nonlinear Dynamics) with a command similar to:

```
Apex3D64.exe -pRawDirName sample.raw -outputDirName peak_picked_sample_folder -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0
```

Alternatively, [IMTBX and Grppr](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5826643/#SD1) can be used for peak picking but this has not been tested for compatibility yet.

In case of proteomics data, a database is needed for annotation. Human, Yeast, Ecoli and benchmark databases are included in the ```lib/databases/``` folder. Custum databases can be created with some simple coding, but currently no GUI has been provided yet.

### Test installation

To test if the installation was indeed successful, run the following command to fully analyze the test dataset, containing 10 minutes from sample A (65% human, 15% yeast, 20% ecoli), sample B (65% human, 30% yeast, 5% ecoli) and a QC in which both samples are pooled.

```
bash run_histopya.sh -p data/test/parameters.json
```

## Usage

### Creation of ion-networks

All input data is required to be in *.csv format, with columns containing centroided 4D peak-picked ion sticks: MZ, RT, DT and Intensity. A fifth column containing an integer stating which *function* this is (MS1 or MS2), is currently required as well, as original development was done on Waters' HDMSE data.

In theory, the DT column can be replaced by other dimensions, such as e.g. the precursor window in a SWATH dataset, but this has not been properly tested. Note that HistoPyA's parameters allow to define column indices for each individual dimension, so *.csv files can be trimmed to only contain these 5 columns as well as containing many other columns with metadata in any order.

When using a Waters' Synapt G2-Si in SWIM-DIA or HDMSE mode for data acquisition, raw data can simply be converted with Waters' Apex3D algorithm to obtain properly formatted *.csv files when a similar command as the following is used:

```
*G:\Sander\UW\UniversalWorkflow_Software\Apx2D\Apex2D_v3.1.0.9.5\Apex3D64.exe* -pRawDirName *W:\MassLynxProjects\160923.PRO\Data\LFQ_SynaptG2Si_UDMSE_Condition_A_Sample_Alpha_01.raw* -outputDirName *G:\Sander\Datasets\LFQ\APEX\LFQ_SynaptG2Si_UDMSE_Condition_A_Sample_Alpha_01* -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0 -maxCPUs 12
```

It is suggested to create a single folder per project within HistoPyA's data folder. We advise to create an *apex* folder in this project folder, containing all *.csv files with 4D ion sticks. For convenience, all *.csv files within this folder are read by default by HistoPyA, to keep all in and output folders clean.

A *.json parameter file containing just the path to this project folder and the *.json databases to use, suffices to perform a default analysis and place all results in a separate folder within the project folder with the following command:

```
bash run_histopya.sh -p data/test/parameters.json
```

All parameters available in *libs/defaults/default_parameters.json* can be overridden by a user-defined parameter file.

Currently 4 *.json databases containing sequences and SwissProt PTMs are available: H. Sapiens, S. Cerevisiae, E. Coli and cRAP. Custom databases can easily be incorporated, but should be generated on by the user.

### Interactive graphical browser

Results can be visualized with

```
bash browse_results.sh -p data/test/parameters.json
```
