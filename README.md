# NeuronEV
This repository contains the source code and example data for the computational analyses described in:

Single-vesicle topology reveals spatially organized neuroglial pathology in Alzheimer’s disease

The code implements the computational framework for single EV digitalization, machine learning-based classification, parameter optimization, and diagnostic performance evaluation.

Overviewhttps://github.com/yifeijiangHIM/NeuronEV/blob/main/README.md
The analysis pipeline consists of the following major steps:

EV digitization and establishment of a multi-dimensional grid
Classification of EV subpopulations, and assign them to the grid based on marker expression profiles
Generating a volcano plot and identify EV subgroups with disease relevance
For EV subgroups with AUCs above threshold, optimize grid boundaries (marker expression ranges) to further improve the diagnostic performances.
Generation of computational results, including the expression profile, P values and fold of change of the EV subgroups.
A simplified workflow is:

1.Input EV data 2.EV digitalization 3.EV subpopulation identification 4.Volcano plot analysis
5.Marker-range optimization based on ROC/AUC 6.Results and figures

Repository Structure
main/ │ ├── README.md ├── LICENSE ├── CITATION.cff

├── code/ │ ├── MAP_Neuron.m 

├── demo/ │ ├── AD Data.zip │ ├── HC Data.zip │ ├── NAD Data.zip │ ├── Biomarker List.doc 

└── results/ │ ├── Single EV List.zip │ ├── Volcano Plot.zip
