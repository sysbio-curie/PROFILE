# Instantiation of patient-specific logical models

This is a repository of data, code and analyses related to the paper "Instantiation of Patient-Specific Logical Models With Multi-Omics Data Allows Clinical Stratification of Patients". 
The paper can accessed here: TBA.

This repository can be cited with its own DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1186270.svg)](https://doi.org/10.5281/zenodo.1186270)

## Getting Started

This set of files and scripts is supposed to be self-sufficient. Please download data and scripts and follow instructions

### Requirements
- Python version 3.0 or greater
- Python's package pandas
- Perl
- R
- MaBoSS requires: flex, bison, gcc and g++

## Patient-specific instantiation pipeline

In the present pipeline, two different datasets may be used (METABRIC or TCGA) and processed for further simulations with two different logical models, either a generic one (*Fumia* model, see [paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069008)) or a breast-specific one (*Zanudo* model, see [paper](https://cancerconvergence.springeropen.com/articles/10.1186/s41236-017-0007-6)).

The following instructions will use METABRIC data and *Fumia* model.


### Generation of patient-specific profiles with R

Patient-profiles generation and all its underlying computations (funcional effect inference, normalization and binarization) have been packed in a *.Rmd* file. Please run **Fumia_META_profiles.Rmd** (in *Scripts/Profiles* folder) with R to generate profiles (in *Results/Profiles* folder) and a corresponding *.html* report (in *Scripts/Profiles* folder)

```
rmarkdown::render("Fumia_META_profiles.Rmd", "html_document")
```

### Simulation of patient-specific models with MaBoSS

Patient-specific profiles are used to instantiate patient-specific models through different methods based on node activity status, initial conditions or transition rates. Simulations are performed by *MaBoSS* software.

However, model instantiation requires to modify *MaBoSS* model files. The whole simulation process is packed in *python* script. Please note that depending on your operating system, you should choose different versions of that script (either **MaBoSS_specific_Mac** or **MaBoSS_specific_Linux**)

A minimal example of simulations can be found in the following *shell* script. This is an example of simulations using mutations and CNA information of METABRIC cohort as node activity status. Only values for model output nodes (dead-end nodes regulated by inhibitors/activators but not regulating any other node) are saved.

```
model=Fumia2013 
sim_case = META_mutations_CNA_asMutants
python3 Scripts/Simulations/MaBoSS_specific_Mac.py $model 2 "Results/Simulations/results_"$sim_case".txt" -s $sim_case -m "Results/Profiles/Fumia_META_mutCNA.csv"
```

This example is available in **example1.sh**, either in MacOS or Linux distribution versions.

Another example, computing the scores for all nodes of the model, is available in **example2.sh**, either in MacOS or Linux distribution versions. The loop structure is required to avoid an exponential increase of computation time.

## Authors

Scripts were mostly designed by Jonas Béal (jonas dot beal at curie dot fr).
See a complete list of authors in the corresponding paper



