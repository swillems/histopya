# HistoPyA


## Overview

HistoPyA was developed to simultaneously analyse data from multiple proteomic samples of similar origin acquired on an LC-IMS-Q-TOF with equal acquisition parameters, preferably in SWIM-DIA mode.

It takes fully peackpicked data, i.e. ion sticks with 4 coordinates (MZ, RT, DT, Intensity), in *.csv format as input and outputs Peptide-to-Ion-Matches in *.csv format (optionally processed by Percolator).

The original paper providing more detail is available in the docs folder.


## Installation

HistoPyA was developed on a Linux Centos 6 distribution, but should work on any Linux flavor with python3.6 installed. A complete installation (including a 300 mb test data set) is done with the following terminal commands, executed from wherever you want to install histopya:

```
git clone https://github.com/swillems/histopya
sh histopya/install/install.sh
# In case python3.6 is not run by calling "python3.6",
# pass the full path as argument, e.g.:
# sh histopya/install/install.sh /usr/bin/python3.6
```

Percolator can improve results, but is not essential. If percolator is installed, but not run by calling "percolator", update the file "lib/defaults/default_parameters.json" so line 65 containing "PERCOLATOR_LOCATION" states the full path, e.g.

```
"PERCOLATOR_LOCATION": "/usr/local/bin/percolator",
```

To test if the installation was indeed successful, run the following command to fully analyze the test dataset, containing 10 minutes from sample A (65% human, 15% yeast, 20% ecoli), sample B (65% human, 30% yeast, 5% ecoli) and a QC in which both samples are pooled.

```
sh histopya/histopya.sh -p data/test/parameters.json
```

## Usage

All input data is required to be in *.csv format, with columns containing centroided 4D peak-picked ion sticks: MZ, RT, DT and Intensity. A fifth column containing an integer stating which "function" this is (MS1 or MS2), is currently required as well, as original development was done on Waters' HDMSE data.

In theory, the DT column can be replaced by other dimensions, such as e.g. the precursor window in a SWATH dataset, but this has not been properly tested. Note that HistoPyA's parameters allow to define column indices for each individual dimension, so *.csv files can be trimmed to only contain these 5 columns as well as containing many other columns with metadata in any order.

When using a Waters' Synapt G2-Si in SWIM-DIA or HDMSE mode for data acquisition, raw data can simply be converted with Waters' Apex3D algorithm to obtain properly formatted *.csv files when a similar command as the following is used:

```
"G:\Sander\UW\UniversalWorkflow_Software\Apx2D\Apex2D_v3.1.0.9.5\Apex3D64.exe" -pRawDirName "W:\MassLynxProjects\160923.PRO\Data\LFQ_SynaptG2Si_UDMSE_Condition_A_Sample_Alpha_01.raw" -outputDirName "G:\Sander\Datasets\LFQ\APEX\LFQ_SynaptG2Si_UDMSE_Condition_A_Sample_Alpha_01" -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0 -maxCPUs 12
```

It is suggested to create a single folder per project within HistoPyA's data folder. We advise to create an "apex" folder in this project folder, containing all *.csv files with 4D ion sticks. For convenience, all *.csv files within this folder are read by default by HistoPyA, to keep all in and output folders clean.

A *.json parameter file containing just the path to this project folder and the *.json databases to use, suffices to perform a default analysis and place all results in a separate folder within the project folder with the following command:

```
sh histopya.sh -p data/test/parameters.json
```

All parameters available in "libs/defaults/default_paramaters.json" can be over riden by a user-defined parameter file.

Currently 4 *.json databases containing sequences and SwissProt PTMs are available: H. Sapiens, S. Cerevisiae, E. Coli and cRAP.
