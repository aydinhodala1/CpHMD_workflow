import numpy as np
import re

def lambda_group(lambda_name: str, group_name: str, atom_charges: list[list[float]], pka: float, dvdl_coeffs: list[float], init_lam: float, index: int) -> str:
    '''
    Create the md.mdp CpHMD string for a lambda group

    Inputs:
    lambda_name (str): Name for lambda group in CpHMD
    group_name (str): Name for index group in index.ndx
    atom_charges (2D list of floats): 1st list corresponds to charges at lambda = 0, 2nd corresponds to charges at lambda = 1. Note 1st number in each list must correspond to first atom in index group, second for second atom etc.
    pka (float): pka of calibration state
    dvdl_coeffs (list of floats): Coefficients for polynomial fit of V_corr
    index (int): Index number for lambda group

    Outputs:
    lambda_str (str): Lines for appending to md.mdp file for CpHMD runs. Note assumes all groups are coupled to the same group of buffer particles (charge-restraint-group-index = 1)
    '''
    #Create each line of lambda-dynamics settings
    name_str = f'lambda-dynamics-group-type{index}-name                          = {lambda_name}'
    states_str = f'lambda-dynamics-group-type{index}-n-states                      = 1'
    lam0_str = f'lambda-dynamics-group-type{index}-state-0-charges               = {str(atom_charges[0])[1:-1].replace(",","")}'
    lam1_str = f'lambda-dynamics-group-type{index}-state-1-charges               = {str(atom_charges[1])[1:-1].replace(",","")}'
    pka_str = f'lambda-dynamics-group-type{index}-state-1-reference-pka         = {pka}'
    dvdl_str = f'lambda-dynamics-group-type{index}-state-1-dvdl-coefficients     =  {dvdl_coeffs}'
    atom_str = f'lambda-dynamics-atom-set{index}-name                         = {lambda_name}'
    group_str = f'lambda-dynamics-atom-set{index}-index-group-name             = {group_name}'
    barrier_str = f'lambda-dynamics-atom-set{index}-barrier                      = 5.0'
    init_str = f'lambda-dynamics-atom-set{index}-initial-lambda               = {init_lam}'
    charge_str = f'lambda-dynamics-atom-set{index}-charge-restraint-group-index = 1'
    
    #Join all lines together, with a gap between group-type settings and atom-set settings
    lambda_str = '\n'.join(
            [name_str, states_str, lam0_str, lam1_str, pka_str, dvdl_str, 
                '',atom_str, group_str, barrier_str, init_str, charge_str]
            )

    return lambda_str

def buf_group(n_buf: int, atom_charges: list[list[float]], dvdl_coeffs: list[float], init_lam: float) -> str:
    '''
    Create the md.mdp CpHMD string for a lambda group

    Inputs:
    atom_charges (2D list of floats): 1st list corresponds to charges at lambda = 0, 2nd corresponds to charges at lambda = 1. Note 1st number in each list must correspond to first atom in index group, second for second atom etc.
    dvdl_coeffs (list of floats): Coefficients for polynomial fit of V_corr

    Outputs:
    buf_str (str): Lines for appending to md.mdp file for CpHMD runs. Note assumes all groups are coupled to the same group of buffer particles (charge-restraint-group-index = 1)
    '''
    #Create each line of lambda-dynamics settings
    name_str = f'lambda-dynamics-group-type1-name                          = BUF'
    states_str = f'lambda-dynamics-group-type1-n-states                      = 1'
    lam0_str = f'lambda-dynamics-group-type1-state-0-charges               = {str(atom_charges[0])[1:-1].replace(",","")}'
    lam1_str = f'lambda-dynamics-group-type1-state-1-charges               = {str(atom_charges[1])[1:-1].replace(",","")}'
    pka_str = f'lambda-dynamics-group-type1-state-1-reference-pka         = 7'
    dvdl_str = f'lambda-dynamics-group-type1-state-1-dvdl-coefficients     =  {dvdl_coeffs}'
    atom_str = f'lambda-dynamics-atom-set1-name                         = BUF'
    group_str = f'lambda-dynamics-atom-set1-index-group-name             = BUF'
    barrier_str = f'lambda-dynamics-atom-set1-barrier                      = 5.0'
    init_str = f'lambda-dynamics-atom-set1-initial-lambda               = {init_lam}'
    buf_str = f'lambda-dynamics-atom-set1-buffer-residue               = yes'
    nbuf_str = f'lambda-dynamics-atom-set1-buffer-residue-multiplier    = {n_buf}'
    charge_str = f'lambda-dynamics-atom-set1-charge-restraint-group-index = 1'

    #Join all lines together, with a gap between group-type settings and atom-set settings
    lambda_str = '\n'.join(
            [name_str, states_str, lam0_str, lam1_str, pka_str, dvdl_str,
                '',atom_str, group_str, barrier_str, init_str, buf_str, nbuf_str, charge_str]
            )

    return lambda_str

#Reads input from input file ("input.in"). Requires each section, both number of particles


titration_dictionary = {}
param_dictionary = {}

main_sec = False
buf_sec = False
titration_num = 1

with open("input.in", "r") as input_file:
    for line in input_file:
        line = line[:-1] 
        #Check if at the beginning of any sections
        if re.match("^MAIN", line):
            main_sec = True
        elif re.match("^BUF", line):
            buf_sec = True

        #Add parameters for a titration group to the dictionary
        elif main_sec == False and buf_sec == False:
            if re.match('^Name', line, re.IGNORECASE):
                titration_dictionary[f'name_{titration_num}'] = re.split('=|#', line)[1]
            elif re.match('^State 0', line, re.IGNORECASE):
                titration_dictionary[f'state0_{titration_num}'] = np.array(re.split(' ', re.split('=|#', line)[1][1:]), dtype = float)
            elif re.match('^State 1', line, re.IGNORECASE):
                titration_dictionary[f'state1_{titration_num}'] = np.array(re.split(' ', re.split('=|#', line)[1][1:]), dtype = float)
            elif re.match('^pka', line, re.IGNORECASE):
                titration_dictionary[f'pka_{titration_num}'] = re.split('=|#', line)[1]
            elif re.match('^correction', line, re.IGNORECASE):
                titration_dictionary[f'correction_{titration_num}'] = re.split('=|#', line)[1]
            elif re.match('^Index', line, re.IGNORECASE):
                titration_dictionary[f'indexgrp_{titration_num}'] = re.split('=|#', line)[1]
            elif re.match('^Initial lambda', line, re.IGNORECASE):
                titration_dictionary[f'initlam_{titration_num}'] = float(re.split('=|#', line)[1])
            elif re.match('^Number', line, re.IGNORECASE):
                titration_dictionary[f'num_{titration_num}'] = int(re.split('=|#', line)[1])
            elif re.match('^END', line):
                titration_num += 1

        #Gather parameters from main section for number of CpHMD active molecules
        elif main_sec == True:
            if re.match('^Buffer/surfactant', line, re.IGNORECASE):
                buf_ratio = int(re.split('=|#', line)[1])
            elif re.match('^END', line):
                main_sec = False

        #Add parameters for buffer particles
        elif buf_sec == True:
            if re.match('^State 0', line, re.IGNORECASE):
                state0_buf = float(re.split(' ', re.split('=|#', line)[1])[1])
            elif re.match('^State 1', line, re.IGNORECASE):
                state1_buf = float(re.split(' ', re.split('=|#', line)[1])[1])
            elif re.match('^correction', line, re.IGNORECASE):
                correction_buf = re.split('=|#', line)[1]
            elif re.match('^END', line):
                buf_sec = False

titration_num -=1

#Calculate required initial lambda and number of buffer particles to ensure the box remains at constant charge CpHMD
#This is assumed to be neutral at lambda = 0 for all titration sites
total_charge = 0
num_buf = 0

for site in range(1,titration_num +1):
    num_buf += buf_ratio * titration_dictionary[f'num_{site}']
    site_charge = (np.sum(titration_dictionary[f'state1_{site}']) - np.sum(titration_dictionary[f'state0_{site}'])) * titration_dictionary[f'initlam_{site}'] * titration_dictionary[f'num_{site}']
    total_charge += site_charge

charge_buf =  - total_charge / num_buf
print(charge_buf)

init_lam_buf = (charge_buf - state0_buf) / (state1_buf - state0_buf)
index = 1

param_dictionary[f'index'] = buf_group(num_buf, charge_buf, correction_buf, init_lam_buf)
index += 1

for titration_site in range(1,titration_num+1):
    try:
        for molecule in range(1,titration_dictionary[f"num_{titration_site}"]+1):
            molecule = str(molecule)
            param_dictionary[f'index'] = lambda_group(titration_dictionary[f"name_{titration_site}"] + molecule,
                         titration_dictionary[f"indexgrp_{titration_site}"] + molecule,
                         list(titration_dictionary[f"state0_{titration_site}"],titration_dictionary[f"state1_{titration_site}"]),
                         titration_dictionary[f"pka_{titration_site}"],
                         titration_dictionary[f"correction_{titration_site}"],
                         titration_dictionary[f'initlam_{titration_site}'],
                         index
                         )
            
    except:
        param_dictionary[f'index'] = lambda_group(titration_dictionary[f"name_{titration_site}"],
                titration_dictionary[f"indexgrp_{titration_site}"],
                list(titration_dictionary[f"state0_{titration_site}"],titration_dictionary[f"state1_{titration_site}"]),
                titration_dictionary[f"pka_{titration_site}"],
                titration_dictionary[f"correction_{titration_site}"],
                titration_dictionary[f'initlam_{titration_site}'],
                index
                )