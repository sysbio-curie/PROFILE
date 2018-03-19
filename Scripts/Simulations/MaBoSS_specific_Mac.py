#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 24 2017
@author: Jonas BÉAL
jonas.beal@curie.fr
"""
#%% Imports

# Imports and tests
import sys
import os
import re
import argparse
import pandas as pd
import time
import multiprocessing as mp
import numpy

current_wd= os.getcwd()

arg = sys.argv

parser = argparse.ArgumentParser()

#Required arguments
parser.add_argument("model", help="name of MaBoSS files (without .cfg or .bnd extension)")
parser.add_argument("save_file", help="save_file is the name of the text file containing final probabilities of outputs (one simulation by line))")

#Optional arguments for computation parameters
parser.add_argument("-n","--num_nodes", type=int, help="nb of nodes in the MaBoSS exec file (ex: 100 to use MaBoSS_100n)")
parser.add_argument("-p","--num_processes", type=int, help="nb of parallel processes during simulations")

#Optional arguments for instantation parameters
parser.add_argument("-i","--inputs", help="initial probabilities of inputs (alternatively called source nodes, i.e not regulated nodes) are set to the specified value, othewise it will be 0.5")
parser.add_argument("-o","--outputs", help="outputs are marked as external nodes, whose final probabilities are saved in the result file")
parser.add_argument("-s","--suffix", help="suffix is added to all intermediate and result files")
parser.add_argument("-m","--mutants", help="name of the csv file containing perturbation profiles to define node mutants (also called node activity status): one profile/patient by line, with multiple genes separated by a comma. Binary 0/1 information (NA tolerated)")
parser.add_argument("-c","--init_cond", help="name of the csv file containing perturbation profiles to define node initial conditions: one profile/patient by line, with multiple genes separated by a comma. Binary 0/1 (NA tolerated) or continuous [0,1] information")
parser.add_argument("-rb","--rates_basic", help="name of the csv file containing perturbation profiles to define reaction rates based on node normalized state: one profile/patient by line, with multiple genes separated by a comma. Continuous [0,1] information")
parser.add_argument("-ra","--rates_advanced", help="name of the csv file containing perturbation profiles to define reaction rates based on activators/inhibitors nodes states: one profile/patient by line, with multiple genes separated by a comma. Continuous [0,1] information")
parser.add_argument("-rf","--rates_factor", help="multiplication factor for rates when using one of the 'rates' method")

args = parser.parse_args()

#%% Process Arguments

print("Arguments:\n")
print(args)
print("\n")


# path to script is used to call other scripts in the same folder
pathtoscript = os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0])))
base_path = os.path.dirname(os.path.dirname(pathtoscript))+"/"

#MaBoSS exec
if args.num_nodes is not None:
    if args.num_nodes<=64:
        maboss_exec = base_path+"MaBoSS/Mac/MaBoSS"
    elif args.num_nodes<=150:
		maboss_exec = base_path+"MaBoSS/Mac/MaBoSS_150n"
        else:
            print("Your model has more than 150 nodes, please recompile MaBoSS with the number of nodes of your model. See MaBoSS compiling help: http://maboss.curie.fr/")

if not os.path.isfile(maboss_exec):
    print("Relevant MaBoSS executable is not available")
    sys.exit(1)

print("MaBoSS executable: "+maboss_exec)

# number of parallel processes
if args.num_processes is not None:
    nbprocesses=args.num_processes
else: nbprocesses=1
print("Number of parallel processes: "+str(nbprocesses))

#Check that the model exist or translate it to maboss format if it is a .bnet file
model=args.model
path_model = base_path+"Models/"+model+"/"+model
print("Model: "+model)

if not os.path.isfile(path_model+".bnd"):
    if not os.path.isfile(path_model+".bnet"):
        print(model+".bnd and .bnet files not found")
        sys.exit(1)
    else:
        print("Convert .bnet file in .bnd file")
        sys.exit(1)
        
if not os.path.isfile(path_model+".cfg"):
    print(model+".cfg file not found")
    sys.exit(1) 

#Inputs/outputs processing
#Define constant nodes and input nodes based on .bnd file
lines = open(path_model+".bnd").readlines()
constant_nodes = dict()
input_nodes = dict()
for i in range(len(lines)):
    if re.search("Node ", lines[i]):
        node_name = lines[i].split(" ")[1]        
        if re.search("  logic = [01];", lines[i+1]):
            constant_nodes[node_name] = int(lines[i+1][-3])
        if re.search("logic = [(]?"+re.escape(node_name)+"[)]?;",lines[i+1]):
             input_nodes[node_name] = 0.5 #Unless otherwise specified, we define a 0.5 default value for input 

#Use inputs argument to modify previous dict
if args.inputs is not None:
    inputs = args.inputs.split(",")
    for input_item in inputs:
        input_nodes[input_item.split(":")[0]] = float(input_item.split(":")[1])

#Define outputs
if os.path.isfile(path_model + ".reggraph"):
    lines = open(path_model + ".reggraph").read().splitlines()
    lines=[x.split(" ") for x in lines]
    targets = [item[0] for item in lines]
    actions = [item[1] for item in lines]
    players = [item[2] for item in lines]
    model_outputs = list(set(targets) - set(players))

if args.outputs is not None:
    outputs = args.outputs.split(",")
elif os.path.isfile(path_model + ".reggraph") and len(model_outputs)>0:
    outputs = model_outputs
else:
    outputs = ["Proliferation","Apoptosis"]

if os.path.isfile(path_model + ".reggraph"):
    print("Inputs: "+str(input_nodes)+" - Constant nodes: "+str(constant_nodes)+" - Outputs: "+str(outputs)+" - Model outputs: "+str(model_outputs))
else:
    print("Inputs: "+str(input_nodes)+" - Constant nodes: "+str(constant_nodes)+" - Outputs: "+str(outputs))
    
#Define patient lists
cases_common = list()

#Define mutants    
if args.mutants is not None:
    mutants = pd.read_csv(args.mutants)
    mutants.rename(columns={'Unnamed: 0':'Name'}, inplace=True)
    mutants_dict = mutants.set_index('Name').to_dict(orient="index")
    #mutants_dict = {k: v for k, v in mutants_dict.items() if not math.isnan(v)}
    cases_mut = list(mutants_dict.keys())
    cases_common.append(cases_mut)

#Define rates
if args.rates_basic is not None:
    rates_f = float(args.rates_factor)
    rates = pd.read_csv(args.rates_basic)
    rates.rename(columns={'Unnamed: 0':'Name'}, inplace=True)
    rates_dict = rates.set_index('Name').to_dict(orient="index")
    cases_rates = list(rates_dict.keys())
    cases_common.append(cases_rates)
    
#Define rates_advanced and factor
if args.rates_advanced is not None:
    
    #Import data to define rates
    rates_a = pd.read_csv(args.rates_advanced)
    rates_a.rename(columns={'Unnamed: 0':'Name'}, inplace=True)
    rates_a_dict = rates_a.set_index('Name').to_dict(orient="index")
    cases_rates_a = list(rates_a_dict.keys())
    cases_common.append(cases_rates_a)
    rates_f = float(args.rates_factor)
    
    rates_nodes=list(rates_a)
    rates_nodes.remove('Name')
    
    #Activator/Inhibitor dictionnaries
    model_activators = {}
    for node in set(targets):
        model_activators[node] = [x for x,y,z in zip(players,actions,targets) if (x in rates_nodes) and y=='->' and z==node]
    model_inhibitors = {}
    for node in set(targets):
        model_inhibitors[node] = [x for x,y,z in zip(players,actions,targets) if (x in rates_nodes) and y=='-|' and z==node]
    
    model_activators={k: v for k, v in model_activators.items() if v}
    model_inhibitors={k: v for k, v in model_inhibitors.items() if v}
    
#Define init_cond profile
if args.init_cond is not None:
    init_cond = pd.read_csv(args.init_cond)
    init_cond.rename(columns={'Unnamed: 0':'Name'}, inplace=True)
    init_cond_dict = init_cond.set_index('Name').to_dict(orient="index")
    cases_exp = list(init_cond_dict.keys())
    cases_common.append(cases_exp)
    
#Define patient lists
if not all(v is None for v in [args.mutants, args.rates_basic, args.rates_advanced, args.init_cond]):
    cases_common = list(set.intersection(*map(set, cases_common)))

#Define save_file  and number of nodes
save_file=args.save_file
nbnodes=args.num_nodes

#Define suffix
if args.suffix is not None:
    suffix = args.suffix
else:
    suffix = "classic"
print("Suffix: "+suffix)

os.system("cp "+path_model+".bnd "+path_model+"_"+suffix+".bnd")
os.system("cp "+path_model+".cfg "+path_model+"_"+suffix+".cfg")
model=model+"_"+suffix
path_model=path_model+"_"+suffix
    
#Prepare simulation
#Define outputs as external nodes
for output in outputs:
    os.system("sed -i '' 's/^"+output+".is_internal *= *TRUE;/"+output+".is_internal=FALSE;/g' "+path_model+"'.cfg'")
    
#Define proper initial conditions for inputs and constant nodes (implicit inputs)
for input_item, input_value in dict(input_nodes, **constant_nodes).items():
    os.system("sed -i '' '/^"+input_item+".istate *=/d' "+path_model+".cfg")
    os.system("echo '["+input_item+"].istate = "+str(input_value)+"[1], "+str(1-input_value)+"[0];' >> "+path_model+".cfg")   

#Define function used to perform the simulation itself and process the output
def perform_MaBoSS_simulation(n_profile, fname, profile_name):
    path_fname=path_model+"_"+str(n_profile)
    
    string_opt=path_fname+".bnd " + path_fname + ".cfg -mb " + maboss_exec
    
    os.chdir(os.path.dirname(path_model))
    os.system(base_path+"Scripts/Simulations/MBSS_FormatTable.pl " + string_opt)
    os.chdir(current_wd)
    #subprocess.Popen([base_path+"Scripts/Simulations/MBSS_FormatTable.pl " + string_opt, cwd=os.path.dirname(path_model))
    
    # store column names (states) and last line with final state distribution in a temp file
    os.system("head -n 1 "+path_fname+"/"+fname+"_probtraj_table.csv > "+path_fname+"/"+fname+"_lastprob_table.csv")
    os.system("tail -n 1 "+path_fname+"/"+fname+"_probtraj_table.csv >> "+path_fname+"/"+fname+"_lastprob_table.csv")
    # extract the output probas from final state distribution in fname+"/"+fname+"_lastprob.csv
    os.system("Rscript "+base_path+"Scripts/Simulations/extract_output_probas.R -i "+path_fname+"/"+fname+"_lastprob_table.csv -o "+path_fname+"/"+fname+"_lastprob -n "+",".join(outputs))
	# add the number of the simulation to the line finally saved
    os.system("echo "+str(n_profile)+"\t"+profile_name+" $(tail -n 1 "+path_fname+"/"+fname+"_lastprob.csv) >> "+save_file)
    #Remove temporary files
    os.system("rm "+path_fname+"/*")
    os.system("rmdir "+path_fname)
    os.system("rm "+path_fname+".bnd")
    os.system("rm "+path_fname+".cfg")
    if profile_name is not "WT":
        previous.append(n_profile)

##Launch simulations depending on the case
os.system("echo n_profile\tPatient_ID\tTime\t"+'\t'.join(outputs)+"\tTH > "+save_file)
if all(v is None for v in [args.mutants, args.rates_basic, args.rates_advanced, args.init_cond]):
    n_profile=0
    fname=model
    perform_MaBoSS_simulation(0,fname, "WT")
    
else:
    manager=mp.Manager()
    previous=manager.list()
    processes = list()
    for i in range(len(cases_common)):
        patient_id=cases_common[i]
        fname=model+"_"+str(i)
        path_fname=path_model+"_"+str(i)
        os.system("cp "+path_model+".bnd "+path_fname+".bnd")
        os.system("cp "+path_model+".cfg "+path_fname+".cfg")

        # set init_cond profiles        
        if args.init_cond is not None:
            patient_dict = init_cond_dict[patient_id]
            patient_dict_red = { k:v for k, v in patient_dict.items() if not numpy.isnan(v) }
            for node, value in patient_dict_red.items():
                value_red = round(value,5)
                os.system("sed -i '' '/^\["+node+"\].istate/d' "+path_fname+"'.cfg'")
                os.system("echo '["+node+"].istate = "+str(value_red)+"[1], "+str(1-value_red)+"[0];' >> "+path_fname+".cfg")
        
        # set rates profiles        
        if args.rates_basic is not None:
            rates_list=rates_dict[patient_id]
            rates_list_red={ k:v for k, v in rates_list.items() if not numpy.isnan(v) }
            for node, value in rates_list_red.items():
                if not numpy.isnan(value):
                    up_value = round(10**(2*numpy.log10(rates_f)*(value-0.5)),5)
                    up_value = round(rates_f**(2*value-1),5)
                    down_value = 1/up_value
                    original_up = float(os.popen("grep -E '^\$u_"+node+" *= *' "+path_fname+".cfg | cut -d'=' -f2 | cut -d';' -f1| tr -d '\n'").read())
                    original_down = float(os.popen("grep -E '^\$d_"+node+" *= *' "+path_fname+".cfg | cut -d'=' -f2 | cut -d';' -f1| tr -d '\n'").read())
                    os.system("sed -i '' 's/u_"+node+" *= *[0-9]*\.*[0-9]*;/u_"+node+"="+str(up_value*original_up)+";/g' "+path_fname+"'.cfg'")
                    os.system("sed -i '' 's/d_"+node+" *= *[0-9]*\.*[0-9]*;/d_"+node+"="+str(down_value*original_down)+";/g' "+path_fname+"'.cfg'")
                    value_red = round(value,5)
                    os.system("sed -i '' '/^\["+node+"\].istate/d' "+path_fname+"'.cfg'")
                    os.system("echo '["+node+"].istate = "+str(value_red)+"[1], "+str(1-value_red)+"[0];' >> "+path_fname+".cfg")
                    
        # set rates_advanced profiles         
        if args.rates_advanced is not None:
            rates_dict_patient=rates_a_dict[patient_id]
            
            for node, value in model_activators.items():
                acti_value = numpy.nanmean([rates_dict_patient[node_name] for node_name in value])
                if not numpy.isnan(acti_value):
                    original_up = float(os.popen("grep -E '^\$u_"+node+" *= *' "+path_fname+".cfg | cut -d'=' -f2 | cut -d';' -f1| tr -d '\n'").read())
                    rate_value = round(rates_f**(2*acti_value-1),5)
                    os.system("sed -i '' 's/u_"+node+" *= *[0-9]*\.*[0-9]*;/u_"+node+"="+str(rate_value*original_up)+";/g' "+path_fname+"'.cfg'")
            
            for node, value in model_inhibitors.items():
                inhi_value = numpy.nanmean([rates_dict_patient[node_name] for node_name in value])
                if not numpy.isnan(inhi_value):
                    original_down = float(os.popen("grep -E '^\$d_"+node+" *= *' "+path_fname+".cfg | cut -d'=' -f2 | cut -d';' -f1| tr -d '\n'").read())
                    rate_value = round(rates_f**(2*inhi_value-1),5)
                    os.system("sed -i '' 's/d_"+node+" *= *[0-9]*\.*[0-9]*;/d_"+node+"="+str(rate_value*original_down)+";/g' "+path_fname+"'.cfg'")
            
            inputs_to_specify = list(set(input_nodes.keys()) & set(rates_nodes))
            for node in inputs_to_specify:
                value = round(rates_dict_patient[node],5)
                os.system("sed -i '' '/^\["+node+"\].istate/d' "+path_fname+"'.cfg'")
                os.system("echo '["+node+"].istate = "+str(value)+"[1], "+str(1-value)+"[0];' >> "+path_fname+".cfg")

        # set mutants profiles 
        if args.mutants is not None:
            mutants_list=mutants_dict[patient_id]
            mutants_list_red={ k:v for k, v in mutants_list.items() if not numpy.isnan(v) }
            for node, value in mutants_list_red.items():
                os.system("sed -i '' '/^\["+node+"\].istate/d' "+path_fname+"'.cfg'")
                if value==0:
                    os.system("sed -i '' 's/u_"+node+" *= *[0-9]*\.*[0-9]*;/u_"+node+"=0;/g' "+path_fname+"'.cfg'")
                    os.system("echo '["+node+"].istate=0[1], 1[0];' >> "+path_fname+".cfg")
                elif value==1:
                    os.system("sed -i '' 's/d_"+node+"= *= *[0-9]*\.*[0-9]*;/d_"+node+"=0;/g' "+path_fname+"'.cfg'")
                    os.system("echo '["+node+"].istate=1[1], 0[1];' >> "+path_fname+".cfg")
                
        while len(previous)<i-(nbprocesses-1):
            time.sleep(1)
        print(str(i)+": "+patient_id)
        p = mp.Process(target = perform_MaBoSS_simulation, args=(i,fname,patient_id))
        p.start()
        processes.append(p)
    for process in processes:
        process.join()
        
    os.system("rm "+path_model+".bnd")
    os.system("rm "+path_model+".cfg")
