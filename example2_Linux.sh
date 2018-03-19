#!/bin/sh

ulimit -n 3000

model=Fumia2013

#Example 2 pf simulations using mutations as node activity status and normalized RNA as transitions rates with an amplification factor of 100
#Loop on model nodes in order to save values for all nodes (using all nodes as outputs in a unique simulation is not recommended)
sim_case=META_mutations_asMutants_RNA_asTransition

list1="GLI,MAX,PTEN,eEF2K,DNA_Damage,Acidosis,Rb,GSH,CyclinB,CyclinA"
list2="ERK1_2,CyclinE,CyclinD,Apoptosis,ATM_ATR,PHDs,mTOR,FADD,MYC,Slug"
list3="WNT,p15,p14ARF,Bak,TSC1_2,Raf,TERT,FOXO,Smad2_3_E2F_p107,beta_catenin"
list4="NF_kB,HIF1,Ras,JNK,APC,UbcH10,p27,p90RSK,p21,Dsh"
list5="RAGS,Rheb,eEF2,Proliferation,RTK,PKC,NF1,TCF,Miz_1,ROS"
list6="TGFbeta,E2F_CycE,AMP_ATP,p53_Mdm2,AMPK,AKT,CHK1_2,Mdm2,GLUT1,FosJun"
list7="E2F,Hypoxia,p53_PTEN,IKK,VEGF,GSK_3_APC_AXIN,Bcl_XL,PIP3,PI3K,Nutrients"
list8="GSK_3,Caspase8,Caspase9,Lactic_acid,p53,E_cadherin,Cdc20,Smad2_3_Smad4_Miz_1,Carcinogen,Smad2_3_Smad4"
list9="TAK1,MYC_MAX,BAX,VHL,BAD,COX412,Snail,Cdh1,DNA_Repair,LDHA"
list10="GFs,CytoC_APAF1,TNFalpha,MXI1,PDK1,Cdh1_UbcH10,p70S6kab,BCL2"

i=1
for listnodes in $list1 $list2 $list3 $list4 $list5 $list6 $list7 $list8 $list9 $list10; do
    python3 Scripts/Simulations/MaBoSS_specific_Linux.py $model -n 98 -p 3 "Results/Simulations/resultsN_"$sim_case"_"$i".txt" -o $listnodes -s "list"$i"_"$sim_case -m "Results/Profiles/Fumia_META_mutCNA.csv" -rb "Results/Profiles/Fumia_META_RNA_norm.csv" -rf 100
    i=$(echo $i+1 | bc)
done

echo "Done"

