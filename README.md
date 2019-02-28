# HistoPyA


## Overview

HistoPyA was developed to simultaneously analyse data from multiple proteomic samples of similar origin acquired on an LC-IMS-Q-TOF with equal acquisition parameters, preferably in SWIM-DIA mode.

It takes fully peackpicked data, i.e. ion sticks with 4 coordinates (MZ, RT, DT, Intensity), in *.csv format as input and outputs Peptide-to-Ion-Matches in *.csv format (optionally processed by Percolator).

The original paper providing more detail is available in the docs folder.


## Installation

HistoPyA was developed on a Linux Centos 6 distribution, but should work on any Linux flavor with python3.6 installed. A complete installation (including a 100 mb zipped test data set) is done with the following terminal commands, executed from wherever you want to install histopya:

```
mkdir histopya
cd histopya
git clone TODO
sh install/install.sh
# In case python3.6 is not run by calling "python3.6", pass the full path as respectively first argument, e.g.:
# sh install/install.sh /usr/bin/python3.6
```

If percolator is not run by calling "percolator", update the file "lib/defaults/default_parameters.json" so line 65 containing "PERCOLATOR_LOCATION" states the full path, e.g.

```
"PERCOLATOR_LOCATION": "/usr/local/bin/percolator",
```

To test if the installation was indeed successful, run the following command to fully analyze the test dataset, containing 10 minutes from sample A (65% human, 15% yeast, 20% ecoli), sample B (65% human, 30% yeast, 5% ecoli) and a QC in which both samples are pooled.

```
sh histopya.sh -p data/test/parameters.json
```

## Usage

TODO Recommended usage is to create a single working folder per project, with parameters.json file and apex folder containing csvs created with Apex3D peackpicking tool (all *.csvs will be auto parsed, so no other csvs should be included).

Recommended apex command line:

```
"G:\Sander\UW\UniversalWorkflow_Software\Apx2D\Apex2D_v3.1.0.9.5\Apex3D64.exe" -pRawDirName "W:\MassLynxProjects\160923.PRO\Data\LFQ_SynaptG2Si_UDMSE_Condition_A_Sample_Alpha_01.raw" -outputDirName "G:\Sander\Datasets\LFQ\APEX\LFQ_SynaptG2Si_UDMSE_Condition_A_Sample_Alpha_01" -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0 -maxCPUs 12
```

Run HistoPyA through command line by providing the "histopya.sh -p data/test.json"
