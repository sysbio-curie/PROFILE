#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 24, 2017
@author: Jonas BÉAL
jonas.beal@curie.fr
Last modification on November 4, 2025 by Arnau Montagud
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

current_wd = os.getcwd()

# path to script is used to call other scripts in the same folder
pathtoscript = os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0])))
base_path = os.path.dirname(os.path.dirname(pathtoscript)) + "/"

# parser definition (do not parse args at import time)
parser = argparse.ArgumentParser()
parser.add_argument("model", help="name of MaBoSS files, without .cfg or .bnd extension (ex: 'CancerModel')")
parser.add_argument("save_file", help="save_file is the name of the text file containing final probabilities of outputs, one simulation by line (ex: 'results.txt')")
parser.add_argument("-sy", "--system", help="Computer OS, eiter Mac or Linux (ex: 'Mac' is you are using MacOS, 'Linux' is the default)")
parser.add_argument("-p", "--num_processes", type=int, help="nb of parallel processes during simulations (ex: '3' if you want to simulate 3 profiles at the same time)")
parser.add_argument("-i", "--inputs", help="initial probabilities of inputs, alternatively called source nodes (i.e not regulated nodes), are set to the specified value, otherwise it will be 0.5 (ex: 'Nutrients:0.3,Androgen:0.6' results in 'Nutrients = 0.3[1], 0.7[0]' and same for Androgen)")
parser.add_argument("-o", "--outputs", help="outputs are marked as external nodes, whose final probabilities are saved in the result file (ex: 'Proliferation,Apoptosis')")
parser.add_argument("-s", "--suffix", help="suffix is added to all intermediate and result files (ex: 'my_simulation')")
parser.add_argument("-cfg", "--CFGbypass", type=bool, help="True if you want to ignore inputs and outputs from general arguments and extract inputs, outputs, internal nodes and initial states' information directly from the provided CFG file (ex: '-cfg True' or '-cfg 1').")
parser.add_argument("-m", "--mutants", help="name of the csv file containing perturbation profiles to define node mutants (also called node activity status): one profile/patient by line, with multiple genes separated by a comma. Binary 0/1 information (NA tolerated)")
parser.add_argument("-c", "--init_cond", help="name of the csv file containing perturbation profiles to define node initial conditions: one profile/patient by line, with multiple genes separated by a comma. Binary 0/1 (NA tolerated) or continuous [0,1] information")
parser.add_argument("-rb", "--rates_basic", help="name of the csv file containing perturbation profiles to define reaction rates based on node normalized state: one profile/patient by line, with multiple genes separated by a comma. Continuous [0,1] information")
parser.add_argument("-ra", "--rates_advanced", help="name of the csv file containing perturbation profiles to define reaction rates based on activators/inhibitors nodes states: one profile/patient by line, with multiple genes separated by a comma. Continuous [0,1] information")
parser.add_argument("-rf", "--rates_factor", help="multiplication factor for rates when using one of the 'rates' methods (ex: '100' in order to have rates between 1/100 and 100)")


#%% Utility / simulation function (module level so multiprocessing can find it)
def perform_MaBoSS_simulation(n_profile, fname, profile_name, path_model_arg, maboss_exec_arg, outputs_arg, save_file_arg, previous_arg=None):
    # profile_name "WT" refers to the unmodified model
    if profile_name != "WT":
        path_fname = path_model_arg + "_" + str(n_profile)
    else:
        path_fname = path_model_arg

    string_opt = path_fname + ".bnd " + path_fname + ".cfg -mb " + maboss_exec_arg
    os.system("chmod +x " + path_fname + ".bnd")
    os.system("chmod +x " + path_fname + ".cfg")
    os.chdir(os.path.dirname(path_model_arg))
    os.system("perl " + base_path + "Scripts/Simulations/MBSS_FormatTable.pl " + string_opt)
    os.chdir(current_wd)

    # store column names (states) and last line with final state distribution in a temp file
    os.system("head -n 1 " + path_fname + "/" + fname + "_probtraj_table.csv > " + path_fname + "/" + fname + "_lastprob_table.csv")
    os.system("tail -n 1 " + path_fname + "/" + fname + "_probtraj_table.csv >> " + path_fname + "/" + fname + "_lastprob_table.csv")
    # extract the output probas from final state distribution
    os.system("Rscript " + base_path + "Scripts/Simulations/extract_output_probas.R -i " + path_fname + "/" + fname + "_lastprob_table.csv -o " + path_fname + "/" + fname + "_lastprob -n " + ",".join(outputs_arg))
    # add the number of the simulation to the line finally saved
    os.system("echo " + str(n_profile) + "\t" + profile_name + " $(tail -n 1 " + path_fname + "/" + fname + "_lastprob.csv) >> " + save_file_arg)
    if profile_name != "WT" and previous_arg is not None:
        previous_arg.append(n_profile)


def main():
    # parse args inside main so module import (by multiprocessing spawn) doesn't run the workflow
    args = parser.parse_args()

    print("Arguments:\n")
    print(args)
    print("\n")

    # Check OS
    if args.system is not None:
        system = args.system
    else:
        system = 'Linux'

    if system == 'Linux':
        sed_string = "sed -i "
    elif system == 'Mac':
        sed_string = "sed -i '' "
    else:
        print("Please use either Linux or Mac OS")
        sys.exit(1)

    # Check that the model exists
    model = args.model
    path_model_local = base_path + "Models/" + model + "/" + model
    print("Model: " + model)

    if not os.path.isfile(path_model_local + ".bnd"):
        if not os.path.isfile(path_model_local + ".bnet"):
            print(model + ".bnd and .bnet files not found")
            sys.exit(1)
        else:
            print("Convert .bnet file in .bnd file")
            sys.exit(1)

    if not os.path.isfile(path_model_local + ".cfg"):
        print(model + ".cfg file not found")
        sys.exit(1)

    # Define the CFGbypass status of the simulation
    if args.CFGbypass is not None:
        CFGbypass = args.CFGbypass
    else:
        CFGbypass = False

    # Define all nodes, constant nodes and input nodes based on .bnd file
    lines = open(path_model_local + ".bnd").readlines()
    constant_nodes = dict()
    input_nodes = dict()
    nodes = list()
    for i in range(len(lines)):
        if re.search(r"^Node .*{$", lines[i]):
            node_name = re.split(r" *{", lines[i])[0].split(" ")[1]
            nodes.append(node_name)
            if re.search(r"logic *= *\(?[01]\)?;", lines[i + 1]):
                constant_nodes[node_name] = int(lines[i + 1][-3])
            if re.search(r"logic *= *\(?" + re.escape(node_name) + r"\)?;", lines[i + 1]):
                input_nodes[node_name] = 0.5

    # Define save_file
    global save_file
    save_file = base_path + "Results/Simulations/" + args.save_file

    # Define number of nodes
    nbnodes = len(nodes)

    # MaBoSS exec selection
    global maboss_exec
    if nbnodes <= 64:
        maboss_exec = base_path + "MaBoSS/" + system + "/MaBoSS"
    elif nbnodes <= 150:
        maboss_exec = base_path + "MaBoSS/" + system + "/MaBoSS_150n"
    else:
        print("Your model has more than 150 nodes, please recompile MaBoSS with the number of nodes of your model. See MaBoSS compiling help: http://maboss.curie.fr/")

    if not os.path.isfile(maboss_exec):
        print("Relevant MaBoSS executable is not available")
        sys.exit(1)

    print("MaBoSS executable: " + maboss_exec)
    os.system("chmod +x " + maboss_exec)

    # number of parallel processes
    if args.num_processes is not None:
        nbprocesses = args.num_processes
    else:
        nbprocesses = 1
    print("Number of parallel processes: " + str(nbprocesses))

    # Use inputs argument to modify previous dict
    if args.inputs is not None:
        inputs = args.inputs.split(",")
        for input_item in inputs:
            input_nodes[input_item.split(":")[0]] = float(input_item.split(":")[1])

    # Define outputs
    if os.path.isfile(path_model_local + ".reggraph"):
        lines = open(path_model_local + ".reggraph").read().splitlines()
        lines = [x.split(" ") for x in lines]
        targets = [item[0] for item in lines]
        actions = [item[1] for item in lines]
        players = [item[2] for item in lines]
        model_outputs = list(set(targets) - set(players))

    if args.outputs is not None:
        outputs_local = args.outputs.split(",")
    elif os.path.isfile(path_model_local + ".reggraph") and len(model_outputs) > 0:
        outputs_local = model_outputs
    else:
        outputs_local = ["Proliferation", "Apoptosis"]

    global outputs
    outputs = outputs_local

    if os.path.isfile(path_model_local + ".reggraph"):
        print("Inputs (model-based and user-defined): " + str(input_nodes) + " - Constant nodes (model-based): " + str(constant_nodes) + " - Simulation outputs (user-defined): " + str(outputs) + " - Model outputs (model-based): " + str(model_outputs))
    else:
        print("Inputs (model-based and user-defined): " + str(input_nodes) + " - Constant nodes (model-based): " + str(constant_nodes) + " - Simulation outputs (user-defined): " + str(outputs))

    # Define patient lists
    cases_common = list()

    # Define mutants
    if args.mutants is not None:
        mutants = pd.read_csv(args.mutants)
        mutants.rename(columns={'Unnamed: 0': 'Name'}, inplace=True)
        mutants_dict = mutants.set_index('Name').to_dict(orient="index")
        cases_mut = list(mutants_dict.keys())
        cases_common.append(cases_mut)

    # Define rates
    if args.rates_basic is not None:
        rates_f = float(args.rates_factor)
        rates = pd.read_csv(args.rates_basic)
        rates.rename(columns={'Unnamed: 0': 'Name'}, inplace=True)
        rates_dict = rates.set_index('Name').to_dict(orient="index")
        cases_rates = list(rates_dict.keys())
        cases_common.append(cases_rates)

    # Define rates_advanced and factor
    if args.rates_advanced is not None:
        rates_a = pd.read_csv(args.rates_advanced)
        rates_a.rename(columns={'Unnamed: 0': 'Name'}, inplace=True)
        rates_a_dict = rates_a.set_index('Name').to_dict(orient="index")
        cases_rates_a = list(rates_a_dict.keys())
        cases_common.append(cases_rates_a)
        rates_f = float(args.rates_factor)

        rates_nodes = list(rates_a)
        rates_nodes.remove('Name')

        # Activator/Inhibitor dictionaries
        model_activators = {}
        for node in set(targets):
            model_activators[node] = [x for x, y, z in zip(players, actions, targets) if (x in rates_nodes) and y == '->' and z == node]
        model_inhibitors = {}
        for node in set(targets):
            model_inhibitors[node] = [x for x, y, z in zip(players, actions, targets) if (x in rates_nodes) and y == '-|' and z == node]

        model_activators = {k: v for k, v in model_activators.items() if v}
        model_inhibitors = {k: v for k, v in model_inhibitors.items() if v}

    # Define init_cond profile
    if args.init_cond is not None:
        init_cond = pd.read_csv(args.init_cond)
        init_cond.rename(columns={'Unnamed: 0': 'Name'}, inplace=True)
        init_cond_dict = init_cond.set_index('Name').to_dict(orient="index")
        cases_exp = list(init_cond_dict.keys())
        cases_common.append(cases_exp)

    # Define patient lists intersection
    if not all(v is None for v in [args.mutants, args.rates_basic, args.rates_advanced, args.init_cond]):
        cases_common = list(set.intersection(*map(set, cases_common)))

    # Define suffix
    if args.suffix is not None:
        suffix = args.suffix
    else:
        suffix = "classic"
    print("Suffix: " + suffix)

    # Prepare simulation copying model files for subsequent modifications
    os.system("cp " + path_model_local + ".bnd " + path_model_local + "_" + suffix + ".bnd")
    os.system("cp " + path_model_local + ".cfg " + path_model_local + "_" + suffix + ".cfg")
    model_local = model + "_" + suffix
    global path_model
    path_model = path_model_local + "_" + suffix

    # Define outputs as external nodes (unless CFGbypass)
    if not CFGbypass:
        for node in nodes:
            # use escaped backslashes so Python doesn't warn about invalid escape sequences
            os.system(sed_string + "'/^\\[*" + node + "\\]*\\.is_internal *= */d' " + path_model + ".cfg")
            if node in outputs:
                os.system("echo '" + node + ".is_internal = FALSE;' >> " + path_model + ".cfg")
            else:
                os.system("echo '" + node + ".is_internal = TRUE;' >> " + path_model + ".cfg")

    # Define proper initial conditions for inputs and constant nodes (implicit inputs)
    if not CFGbypass:
        for input_item, input_value in dict(input_nodes, **constant_nodes).items():
            os.system(sed_string + "'/^\\[*" + input_item + "\\]*\\.istate *=/d' " + path_model + ".cfg")
            os.system("echo '[" + input_item + "].istate = " + str(input_value) + "[1], " + str(1 - input_value) + "[0];' >> " + path_model + ".cfg")

    # Launch simulations depending on the case
    os.system("echo n_profile\tPatient_ID\tTime\t" + '\t'.join(outputs) + "\tTH > " + save_file)
    if all(v is None for v in [args.mutants, args.rates_basic, args.rates_advanced, args.init_cond]):
        n_profile = 0
        fname = model_local
        # no multiprocessing here; pass None for previous
        perform_MaBoSS_simulation(0, fname, "WT", path_model, maboss_exec, outputs, save_file, None)
    else:
        manager = mp.Manager()
        global previous
        previous = manager.list()
        processes = list()
        for i in range(len(cases_common)):
            patient_id = cases_common[i]
            fname = model_local + "_" + str(i)
            path_fname = path_model + "_" + str(i)
            os.system("cp " + path_model + ".bnd " + path_fname + ".bnd")
            os.system("cp " + path_model + ".cfg " + path_fname + ".cfg")

            # set init_cond profiles
            if args.init_cond is not None:
                patient_dict = init_cond_dict[patient_id]
                patient_dict_red = {k: v for k, v in patient_dict.items() if not numpy.isnan(v)}
                for node, value in patient_dict_red.items():
                    value_red = round(value, 5)
                    os.system(sed_string + "'/^\\[*" + node + "\\]*\\.istate/d' " + path_fname + ".cfg")
                    os.system("echo '[" + node + "].istate = " + str(value_red) + "[1], " + str(1 - value_red) + "[0];' >> " + path_fname + ".cfg")

            # set rates profiles
            if args.rates_basic is not None:
                rates_list = rates_dict[patient_id]
                rates_list_red = {k: v for k, v in rates_list.items() if not numpy.isnan(v)}
                for node, value in rates_list_red.items():
                    if not numpy.isnan(value):
                        up_value = round(rates_f ** (2 * value - 1), 5)
                        down_value = round(1 / up_value, 5)
                        original_up = float(os.popen("grep -E '^\\$u_" + node + " *= *' " + path_fname + ".cfg | cut -d'=' -f2 | cut -d';' -f1 | tr -d '\\n'").read())
                        original_down = float(os.popen("grep -E '^\\$d_" + node + " *= *' " + path_fname + ".cfg | cut -d'=' -f2 | cut -d';' -f1 | tr -d '\\n'").read())
                        os.system(sed_string + "'s/u_" + node + " *= *[0-9]*\\.*[0-9]*;/u_" + node + "=" + str(up_value * original_up) + ";/g' " + path_fname + ".cfg")
                        os.system(sed_string + "'s/d_" + node + " *= *[0-9]*\\.*[0-9]*;/d_" + node + "=" + str(down_value * original_down) + ";/g' " + path_fname + ".cfg")
                        value_red = round(value, 5)
                        os.system(sed_string + "'/^\\[*" + node + "\\]*\\.istate/d' " + path_fname + ".cfg")
                        os.system("echo '[" + node + "].istate = " + str(value_red) + "[1], " + str(1 - value_red) + "[0];' >> " + path_fname + ".cfg")

            # set rates_advanced profiles
            if args.rates_advanced is not None:
                rates_dict_patient = rates_a_dict[patient_id]

                for node, value in model_activators.items():
                    acti_value = numpy.nanmean([rates_dict_patient[node_name] for node_name in value])
                    if not numpy.isnan(acti_value):
                        original_up = float(os.popen("grep -E '^\\$u_" + node + " *= *' " + path_fname + ".cfg | cut -d'=' -f2 | cut -d';' -f1 | tr -d '\\n'").read())
                        rate_value = round(rates_f ** (2 * acti_value - 1), 5)
                        os.system(sed_string + "'s/u_" + node + " *= *[0-9]*\\.*[0-9]*;/u_" + node + "=" + str(rate_value * original_up) + ";/g' " + path_fname + ".cfg")

                for node, value in model_inhibitors.items():
                    inhi_value = numpy.nanmean([rates_dict_patient[node_name] for node_name in value])
                    if not numpy.isnan(inhi_value):
                        original_down = float(os.popen("grep -E '^\\$d_" + node + " *= *' " + path_fname + ".cfg | cut -d'=' -f2 | cut -d';' -f1 | tr -d '\\n'").read())
                        rate_value = round(rates_f ** (2 * inhi_value - 1), 5)
                        os.system(sed_string + "'s/d_" + node + " *= *[0-9]*\\.*[0-9]*;/d_" + node + "=" + str(rate_value * original_down) + ";/g' " + path_fname + ".cfg")

                inputs_to_specify = list(set(input_nodes.keys()) & set(rates_nodes))
                for node in inputs_to_specify:
                    value = round(rates_dict_patient[node], 5)
                    os.system(sed_string + "'/^\\[*" + node + "\\]*\\.istate/d' " + path_fname + ".cfg")
                    os.system("echo '[" + node + "].istate = " + str(value) + "[1], " + str(1 - value) + "[0];' >> " + path_fname + ".cfg")

            # set mutants profiles
            if args.mutants is not None:
                mutants_list = mutants_dict[patient_id]
                mutants_list_red = {k: v for k, v in mutants_list.items() if not numpy.isnan(v)}
                for node, value in mutants_list_red.items():
                    os.system(sed_string + "'/^\\[*" + node + "\\]*\\.istate/d' " + path_fname + ".cfg")
                    if value == 0:
                        os.system(sed_string + "'s/u_" + node + " *= *[0-9]*\\.*[0-9]*;/u_" + node + "=0;/g' " + path_fname + ".cfg")
                        os.system("echo '[" + node + "].istate=0[1], 1[0];' >> " + path_fname + ".cfg")
                    elif value == 1:
                        os.system(sed_string + "'s/d_" + node + " *= *[0-9]*\\.*[0-9]*;/d_" + node + "=0;/g' " + path_fname + ".cfg")
                        os.system("echo '[" + node + "].istate=1[1], 0[0];' >> " + path_fname + ".cfg")

            while len(previous) < i - (nbprocesses - 1):
                time.sleep(1)
            print(str(i) + ": " + patient_id)
            p = mp.Process(target=perform_MaBoSS_simulation, args=(i, fname, patient_id, path_model, maboss_exec, outputs, save_file, previous))
            p.start()
            processes.append(p)
        for process in processes:
            process.join()


if __name__ == '__main__':
    # On Windows multiprocessing, freeze_support is helpful; it's safe to call on Unix
    mp.freeze_support()
    main()
