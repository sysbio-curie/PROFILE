#!/bin/sh

ulimit -n 3000

model=Cohen

#Example of simulations using mutations+CNA as node activity status and normalized RNA as transitions rates with an amplification factor of 100
#Loop on model nodes in order to save values for all nodes (using all nodes as outputs in a unique simulation is not recommended and can result in an exponential computation time)

sim_case=TCGA_mutCNA_asMutants_RNA_asTransition

list1="Apoptosis,Migration,Invasion,EMT"
#list2="Snai1,Snai2,Vim,Zeb1,Zeb2,Twist1,Cdh1,Cdh2"
#list3="p73,p63,p53"
#list4="p21,CTNNB1,DKK1,TGFbeta,NICD"


i=1

for listnodes in  $list1; do
#for listnodes in $list1 $list2 $list3 $list4 ; do

    #Simulate first WT condition:
    python3 Scripts/Simulations/MaBoSS_specific.py $model -sy Mac -p 2 "resultsN_WT_"$i".txt" -o $listnodes -s "list"$i"_WT"  

    #Simulate personalized models with mut&CNA as strict node variants (constraining nodes) and RNA as soft node variants (playing with transition rates)
    python3 Scripts/Simulations/MaBoSS_specific.py $model -sy Mac -p 2 "resultsN_"$sim_case"_"$i".txt" -o $listnodes -s "list"$i"_"$sim_case -m "Results/Profiles/Cohen_TCGA_colon_mutCNA.csv" -rb "Results/Profiles/Cohen_TCGA_colon_RNA_norm.csv" -rf 100
    i=$(echo $i+1 | bc)

done

echo "Done"

