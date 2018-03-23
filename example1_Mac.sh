#!/bin/sh

ulimit -n 3000

model=Fumia2013

#Example 1 of simulations using mutations and CNA information of METABRIC cohort as node activity status
#Only values for model output nodes (dead-end nodes regulated by inhibitors/activators but not regulating any other node) are saved 
sim_case=META_mutations_CNA_asMutants
python3 Scripts/Simulations/MaBoSS_specific.py $model -sy Mac -p 3 "Results/Simulations/results_"$sim_case".txt" -s $sim_case -m "Results/Profiles/Fumia_META_mutCNA.csv" 
