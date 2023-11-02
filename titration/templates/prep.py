########### Import modules ###########

import numpy as np
import os
import re
import shutil

########### Define functions ###########

def index_group(name: str, atoms: list) -> str:
    '''
    Generate index group string for index.ndx for each titration group

    Inputs:
    name (str): Name for lambda group in index.ndx
    atoms (list): Atoms numbers in topology

    Outputs:
    index_str (str): Formatted string for appending to index.ndx file
    '''
    #Create each index group string
    title_str = f'\n[ {name} ]'
    atom_str = f'{str(atoms)[1:-1].replace(",","")}\n\n'

    index_str = '\n'.join([title_str, atom_str])
    return index_str

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
    charge_str = f'lambda-dynamics-atom-set{index}-charge-restraint-group-index = 1\n\n'
    
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
    n_buf :
    atom_charges (2D list of floats): 1st list corresponds to charges at lambda = 0, 2nd corresponds to charges at lambda = 1. Note 1st number in each list must correspond to first atom in index group, second for second atom etc.
    dvdl_coeffs (list of floats): Coefficients for polynomial fit of V_corr
    init_lam:

    Outputs:
    buf_str (str): Lines for appending to md.mdp file for CpHMD runs. Note assumes all groups are coupled to the same group of buffer particles (charge-restraint-group-index = 1)
    '''
    #Create each line of lambda-dynamics settings
    name_str = f'lambda-dynamics-group-type1-name                          = BUF'
    states_str = f'lambda-dynamics-group-type1-n-states                      = 1'
    lam0_str = f'lambda-dynamics-group-type1-state-0-charges               = {str(atom_charges[0])}'
    lam1_str = f'lambda-dynamics-group-type1-state-1-charges               = {str(atom_charges[1])}'
    pka_str = f'lambda-dynamics-group-type1-state-1-reference-pka         = 7'
    dvdl_str = f'lambda-dynamics-group-type1-state-1-dvdl-coefficients     =  {dvdl_coeffs}'
    atom_str = f'lambda-dynamics-atom-set1-name                         = BUF'
    group_str = f'lambda-dynamics-atom-set1-index-group-name             = BUF'
    barrier_str = f'lambda-dynamics-atom-set1-barrier                      = 0.0'
    init_str = f'lambda-dynamics-atom-set1-initial-lambda               = {init_lam}'
    buf_str = f'lambda-dynamics-atom-set1-buffer-residue               = yes'
    nbuf_str = f'lambda-dynamics-atom-set1-buffer-residue-multiplier    = {n_buf}'
    charge_str = f'lambda-dynamics-atom-set1-charge-restraint-group-index = 1\n\n'

    #Join all lines together, with a gap between group-type settings and atom-set settings
    lambda_str = '\n'.join(
            [name_str, states_str, lam0_str, lam1_str, pka_str, dvdl_str, '',
                atom_str, group_str, barrier_str, init_str, buf_str, nbuf_str, charge_str]
            )

    return lambda_str

########### Determine initial settings ###########

#Reads input from input file ("input.in"). Requires each section, both number of particles

titration_dictionary = {}
index_dictionary = {}
param_dictionary = {}

main_sec = False
buf_sec = False
titration_num = 1
num_ions = 0

with open("../../input.in", "r") as input_file:
    for line in input_file:
        line = line.strip("\n")
        #Check if at the beginning of any sections
        if re.match("^MAIN", line):
            main_sec = True
        elif re.match("^BUF", line):
            buf_sec = True

        #Add parameters for a titration group to the dictionary
        elif not main_sec or not buf_sec:
            if re.match('^Name', line, re.IGNORECASE):
                titration_dictionary[f'name_{titration_num}'] = re.split('=|#', line)[1].strip()
            elif re.match('^State 0', line, re.IGNORECASE):
                titration_dictionary[f'state0_{titration_num}'] = np.array(re.split(' ', re.split('=|#', line)[1][1:]), dtype = float)
            elif re.match('^State 1', line, re.IGNORECASE):
                titration_dictionary[f'state1_{titration_num}'] = np.array(re.split(' ', re.split('=|#', line)[1][1:]), dtype = float)
            elif re.match('^pka', line, re.IGNORECASE):
                titration_dictionary[f'pka_{titration_num}'] = re.split('=|#', line)[1].strip()
            elif re.match('^correction', line, re.IGNORECASE):
                titration_dictionary[f'correction_{titration_num}'] = re.split('=|#', line)[1].strip()
            elif re.match('^Index', line, re.IGNORECASE):
                titration_dictionary[f'indexgrp_{titration_num}'] = re.split('=|#', line)[1].strip()
            elif re.match('^Initial lambda', line, re.IGNORECASE):
                titration_dictionary[f'initlam_{titration_num}'] = float(re.split('=|#', line)[1])
            elif re.match('^Atom indicies', line, re.IGNORECASE):
                titration_dictionary[f'atoms_{titration_num}'] = np.array(re.split(' ', re.split('=|#', line)[1][1:]), dtype = int)
            elif re.match('^Number', line, re.IGNORECASE):
                titration_dictionary[f'num_{titration_num}'] = int(re.split('=|#', line)[1])
            elif re.match('^Atoms in molecule', line, re.IGNORECASE):
                titration_dictionary[f'mollen_{titration_num}'] = int(re.split('=|#', line)[1])
            elif re.match('^Structure', line, re.IGNORECASE):
                titration_dictionary[f'structure_{titration_num}'] = "../../"+re.split('=|#', line)[1].strip()
            elif re.match('^itpfile', line, re.IGNORECASE):
                titration_dictionary[f'itp_{titration_num}'] = re.split('=|#', line)[1].strip()
            elif re.match('^END*', line):
                titration_num += 1

        #Gather parameters from main section for number of CpHMD active molecules
        elif main_sec:
            if re.match('^Buffer/surfactant', line, re.IGNORECASE):
                buf_ratio = int(re.split('=|#', line)[1])
            elif re.match('^Number of ions', line, re.IGNORECASE):
                num_ions = int(re.split('=|#', line)[1])
            elif re.match('^END', line):
                main_sec = False

        #Add parameters for buffer particles
        elif buf_sec:
            if re.match('^State 0', line, re.IGNORECASE):
                state0_buf = float(re.split(' ', re.split('=|#', line)[1])[1])
            elif re.match('^State 1', line, re.IGNORECASE):
                state1_buf = float(re.split(' ', re.split('=|#', line)[1])[1])
            elif re.match('^correction', line, re.IGNORECASE):
                correction_buf = re.split('=|#', line)[1]
            elif re.match('^END', line):
                buf_sec = False

titration_num -=1

#Check that both or neither number of atoms in molecule and number of molecules are provided.
#Causes problems generating index.ndx downstream otherwise
for site in range(1,titration_num +1):
    if ((f'mollen_{site}' in titration_dictionary.keys()) ^ (f'num_{site}' in titration_dictionary.keys())):
        raise Exception("You must provide either both number of molecules and number of atoms in molecules or neither and declare each individually.")
    elif not (f'num_{site}' in titration_dictionary.keys()):
        titration_dictionary[f'num_{site}'] = 1

#Calculate required initial lambda and number of buffer particles to ensure the box remains at constant charge CpHMD
#This is assumed to be neutral at lambda = 0 for all titration sites
total_charge = 0
num_buf = 0

for site in range(1,titration_num+1):
    num_buf += buf_ratio * titration_dictionary[f'num_{site}']
    site_charge = (np.sum(titration_dictionary[f'state1_{site}']) - np.sum(titration_dictionary[f'state0_{site}'])) * titration_dictionary[f'initlam_{site}'] * titration_dictionary[f'num_{site}']
    total_charge += site_charge

charge_buf =  - total_charge / num_buf

init_lam_buf = (charge_buf - state0_buf) / (state1_buf - state0_buf)

########### Generate intial structure and topology ###########

#Generate initial structure in a 20 A radius sphere using packmol
with open("in.inp","w") as packmolfile:
    packmolfile.write(f"tolerance 2.0\n output out.pdb\n filetype pdb\n\n")
    for site in range(1,titration_num+1):
        packmolfile.write(f'structure {titration_dictionary[f"structure_{site}"]}\n')
        packmolfile.write(f'    number {titration_dictionary[f"num_{site}"]}\n')
        packmolfile.write(f'	inside sphere 0 0 0 20\n')
        packmolfile.write(f'end structure')
os.system('packmol < in.inp')

#Write initial topology
with open("topol.top","w") as topologyfile:
    index = 1
    topologyfile.write(f'#include "../../charmm36-mar2019-cphmd.ff/forcefield.itp"\n#include "../../charmm36-mar2019-cphmd.ff/ions.itp"\n#include "../../charmm36-mar2019-cphmd.ff/tip3p.itp"\n')
    for site in range(1,titration_num+1):
        shutil.copy2(f"../../{titration_dictionary[f'itp_{site}']}", f"../../charmm36-mar2019-cphmd.ff")
        topologyfile.write(f'#include "../../charmm36-mar2019-cphmd.ff/{titration_dictionary[f"itp_{site}"]}"\n')
    topologyfile.write(f'\n\n[ system ]\nSUF\n\n[ molecules ]\n')
    for site in range(1,titration_num+1):
        topologyfile.write(f'{titration_dictionary[f"name_{site}"]}               {titration_dictionary[f"num_{site}"]}\n')

########### Edit structure to be ready for equilibration ###########

#Enlargen box
os.system("gmx_mpi editconf -f out.pdb -o lrg.gro -c -d 2.0 -bt cubic")

#Add water
os.system("gmx_mpi solvate -cp lrg.gro -cs spc216.gro -o wet.gro -p topol.top")

#Insert buffer particles
os.system("gmx_mpi grompp -f buf.mdp -c wet.gro -o buf.tpr -p topol.top")
os.system(f"echo SOL |gmx_mpi genion -s buf.tpr -o buf.gro -pname BUF -np {num_buf} -p topol.top")

#Insert ions
os.system("gmx_mpi grompp -f buf.mdp -c buf.gro -o ion.tpr -p topol.top")
os.system(f"echo SOL |gmx_mpi genion -s ion.tpr -o ion.gro -pname NA -np {num_ions} -nname CL -nn {num_ions} -p topol.top")

########### Generate index.ndx ###########

#Generate initial index.ndx using GROMACS

os.system("echo q | gmx_mpi make_ndx -f ion.gro")

#Generate parameter strings for index.ndx file
index = 0

for site in range(1,titration_num+1):
    try:
        for molecule in range(1,titration_dictionary[f"num_{site}"]+1):
            index += 1
            index_dictionary[f'{index}'] = index_group(
                        titration_dictionary[f"indexgrp_{site}"] + "_" + str(molecule),
                        titration_dictionary[f'atoms_{site}'] + titration_dictionary[f'mollen_{site}'] * (molecule - 1)
                        )
            
    except:
        index_dictionary[f'{index}'] = index_group(
                    titration_dictionary[f"indexgrp_{site}"],
                    titration_dictionary[f'atoms_{titration_num}']
                    )

#Append index strings to index
with open("index.ndx","a") as indexfile:
    for i in range(1, index+1):
        indexfile.write(index_dictionary[f'{i}'])

########### Equilibrate structure ###########

#MPI ranks and OpenMP threads. Defaults to serial if not provided as environment variables
try:
    os.environ["OMP_NUM_THREADS"]
except:
    os.system("export OMP_NUM_THREADS=1")

try:
    mpi_ranks = os.environ["NSLOTS"]
except:
    mpi_ranks = 1


#Energy minimization
os.system("gmx_mpi grompp -f min.mdp -c ion.gro -o min.tpr -maxwarn 2")
os.system(f"mpirun -n {mpi_ranks} gmx_mpi mdrun -deffnm min -npme 0")

#NVT equilibration
os.system("gmx_mpi grompp -f nvt.mdp -c min.gro -r min.gro -o nvt.tpr -maxwarn 2")
os.system(f"mpirun -n {mpi_ranks} gmx_mpi mdrun -deffnm nvt -npme 0")

#NPT equilibration
os.system("gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -o npt.tpr -maxwarn 2")
os.system(f"mpirun -n {mpi_ranks} gmx_mpi mdrun -deffnm npt -npme 0")

########### Edit md.mdp ###########

#Generate parameter strings for md.mdp file
index = 1
param_dictionary[f'{index}'] = buf_group(num_buf, [state0_buf, state1_buf], correction_buf, init_lam_buf)

for site in range(1,titration_num+1):
    try:
        for molecule in range(1,titration_dictionary[f"num_{site}"]+1):
            index += 1
            #Check for both length and number provided
            molecule = str(molecule)
            length = str(titration_dictionary[f'mollen_{site}'])
            param_dictionary[f'{index}'] = lambda_group(
                        titration_dictionary[f"name_{site}"] + "_" + molecule,
                        titration_dictionary[f"indexgrp_{site}"] + "_" + molecule,
                        [titration_dictionary[f"state0_{site}"],titration_dictionary[f"state1_{site}"]],
                        titration_dictionary[f"pka_{site}"],
                        titration_dictionary[f"correction_{site}"],
                        titration_dictionary[f'initlam_{site}'],
                        index
                        )
            
    except:
        param_dictionary[f'{index}'] = lambda_group(
                titration_dictionary[f"name_{site}"],
                titration_dictionary[f"indexgrp_{site}"],
                [titration_dictionary[f"state0_{site}"],titration_dictionary[f"state1_{site}"]],
                titration_dictionary[f"pka_{site}"],
                titration_dictionary[f"correction_{site}"],
                titration_dictionary[f'initlam_{site}'],
                index
                )

#Write parameter strings to md.mdp
with open("md.mdp","a") as mdfile:
    mdfile.write(f'\nlambda-dynamics-number-lambda-group-types              = {index}\n')
    mdfile.write(f'lambda-dynamics-number-atom-collections                = {index}\n\n')
    for i in range(1, index+1):
        mdfile.write(param_dictionary[f'{i}'])

########### Prepare initial files for production ###########

#Create production directory and add analysis script
if os.path.exists("../production"):
    shutil.rmtree("../production")
os.mkdir("../production")
shutil.copy2("titration.py", "../production")

#Generate pH directories and populate with input files
pH_range = np.linspace(0,14,29)

for pH in pH_range:
    dirname = f"../production/{pH}"
    os.mkdir(dirname)
    shutil.copy2("npt.gro", dirname)
    shutil.copy2("index.ndx", dirname)
    shutil.copy2("production.py", dirname)

    #Point to correct force field directory
    shutil.copy2("topol.top", dirname)
    with open(f"{dirname}/topol.top","r") as topolfile:
        topol = topolfile.read()
        topol = topol.replace("../..", "../../..")
    with open(f"{dirname}/topol.top","w") as topolfile:
        topolfile.write(topol)

    #Add pH parameter to md.mdp
    shutil.copy2("md.mdp", dirname)
    with open(f"{dirname}/md.mdp","a") as mdfile:
        mdfile.write(f'\nlambda-dynamics-simulation-ph                          = {pH}\n')