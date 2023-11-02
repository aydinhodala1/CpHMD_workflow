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

def lambda_group(lambda_name: str, group_name: str, atom_charges: list[list[float]]) -> str:
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
    name_str = f'lambda-dynamics-group-type2-name                          = {lambda_name}'
    states_str = f'lambda-dynamics-group-type2-n-states                      = 1'
    lam0_str = f'lambda-dynamics-group-type2-state-0-charges               = {str(atom_charges[0])[1:-1].replace(",","")}'
    lam1_str = f'lambda-dynamics-group-type2-state-1-charges               = {str(atom_charges[1])[1:-1].replace(",","")}'
    pka_str = f'lambda-dynamics-group-type2-state-1-reference-pka         = 7'
    dvdl_str = f'lambda-dynamics-group-type2-state-1-dvdl-coefficients     =  0'
    atom_str = f'lambda-dynamics-atom-set2-name                         = {lambda_name}'
    group_str = f'lambda-dynamics-atom-set2-index-group-name             = {group_name}'
    barrier_str = f'lambda-dynamics-atom-set2-barrier                      = 5.0'
    
    #Join all lines together, with a gap between group-type settings and atom-set settings
    lambda_str = '\n'.join(
            [name_str, states_str, lam0_str, lam1_str, pka_str, dvdl_str, 
                '',atom_str, group_str, barrier_str]
            )

    return lambda_str

def buf_group() -> str:
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
    lam0_str = f'lambda-dynamics-group-type1-state-0-charges               = 0'
    lam1_str = f'lambda-dynamics-group-type1-state-1-charges               = 1'
    pka_str = f'lambda-dynamics-group-type1-state-1-reference-pka         = 7'
    dvdl_str = f'lambda-dynamics-group-type1-state-1-dvdl-coefficients     =  0'
    atom_str = f'lambda-dynamics-atom-set1-name                         = BUF'
    group_str = f'lambda-dynamics-atom-set1-index-group-name             = BUF'
    barrier_str = f'lambda-dynamics-atom-set1-barrier                      = 0.0'
    buf_str = f'lambda-dynamics-atom-set1-buffer-residue               = yes'
    nbuf_str = f'lambda-dynamics-atom-set1-buffer-residue-multiplier    = 1'

    #Join all lines together, with a gap between group-type settings and atom-set settings
    lambda_str = '\n'.join(
            [name_str, states_str, lam0_str, lam1_str, pka_str, dvdl_str, '',
                atom_str, group_str, barrier_str, buf_str, nbuf_str]
            )

    return lambda_str

########### Determine initial settings ###########

#Reads input from input file ("input.in"). Requires each section, both number of particles
main_sec = False

with open("../../input.in", "r") as input_file:
    for line in input_file:
        line = line.strip("\n")
        #Check if at the beginning of any sections
        if re.match("^MAIN", line):
            main_sec = True

        #Gather parameters from main section for number of CpHMD active molecules
        elif main_sec:
            if re.match('^END', line):
                main_sec = False

        #Add parameters for a titration group to the dictionary
        else:
            if re.match('^Name', line, re.IGNORECASE):
                name = re.split('=|#', line)[1].strip()
            elif re.match('^State 0', line, re.IGNORECASE):
                state0 = np.array(re.split(' ', re.split('=|#', line)[1][1:]), dtype = float)
            elif re.match('^State 1', line, re.IGNORECASE):
                state1 = np.array(re.split(' ', re.split('=|#', line)[1][1:]), dtype = float)
            elif re.match('^Index', line, re.IGNORECASE):
                indexgrp = re.split('=|#', line)[1].strip()
            elif re.match('^Atom indicies', line, re.IGNORECASE):
                atoms = np.array(re.split(' ', re.split('=|#', line)[1][1:]), dtype = int)
            elif re.match('^Structure', line, re.IGNORECASE):
                structure = "../../"+re.split('=|#', line)[1].strip()
            elif re.match('^itpfile', line, re.IGNORECASE):
                itpfile = re.split('=|#', line)[1].strip()


#Calculate required initial lambda and number of buffer particles to ensure the box remains at constant charge CpHMD
#This is assumed to be neutral at lambda = 0 for all titration sites


########### Generate intial structure and topology ###########

#Generate initial structure in a 20 A radius sphere using packmol
with open("in.inp","w") as packmolfile:
    packmolfile.write(f"tolerance 2.0\n output out.pdb\n filetype pdb\n\n")
    packmolfile.write(f'structure {structure}\n')
    packmolfile.write(f'    number 1\n')
    packmolfile.write(f'	inside sphere 0 0 0 20\n')
    packmolfile.write(f'end structure')
os.system('packmol < in.inp')

#Write initial topology
with open("topol.top","w") as topologyfile:
    index = 1
    topologyfile.write(f'#include "../../charmm36-mar2019-cphmd.ff/forcefield.itp"\n#include "../../charmm36-mar2019-cphmd.ff/ions.itp"\n#include "../../charmm36-mar2019-cphmd.ff/tip3p.itp"\n')
    shutil.copy2(f"../../{itpfile}", f"../../charmm36-mar2019-cphmd.ff")
    topologyfile.write(f'#include "../../charmm36-mar2019-cphmd.ff/{itpfile}"\n')
    topologyfile.write(f'\n\n[ system ]\nSUF\n\n[ molecules ]\n')
    topologyfile.write(f'{name}               1\n')

########### Edit structure to be ready for equilibration ###########

#Enlargen box
os.system("gmx_mpi editconf -f out.pdb -o lrg.gro -c -d 2.0 -bt cubic")

#Add water
os.system("gmx_mpi solvate -cp lrg.gro -cs spc216.gro -o wet.gro -p topol.top")

#Insert buffer particles
os.system("gmx_mpi grompp -f buf.mdp -c wet.gro -o buf.tpr -p topol.top")
os.system(f"echo SOL |gmx_mpi genion -s buf.tpr -o buf.gro -pname BUF -np 1 -p topol.top")

########### Generate index.ndx ###########

#Generate initial index.ndx using GROMACS
os.system("echo q | gmx_mpi make_ndx -f buf.gro")

#Append index strings to index
with open("index.ndx","a") as indexfile:
    indexfile.write(index_group(indexgrp,atoms))

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
os.system("gmx_mpi grompp -f min.mdp -c buf.gro -o min.tpr -maxwarn 2")
os.system(f"mpirun -n {mpi_ranks} gmx_mpi mdrun -deffnm min -npme 0")

#NVT equilibration
os.system("gmx_mpi grompp -f nvt.mdp -c min.gro -r min.gro -o nvt.tpr -maxwarn 2")
os.system(f"mpirun -n {mpi_ranks} gmx_mpi mdrun -deffnm nvt -npme 0")

#NPT equilibration
os.system("gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -o npt.tpr -maxwarn 2")
os.system(f"mpirun -n {mpi_ranks} gmx_mpi mdrun -deffnm npt -npme 0")

########### Edit md.mdp ###########

#Generate parameter strings for md.mdp file

#Write parameter strings to md.mdp
with open("md.mdp","a") as mdfile:
    mdfile.write(f'\nlambda-dynamics-number-lambda-group-types              = 2\n')
    mdfile.write(f'lambda-dynamics-number-atom-collections                = 2\n\n')
    mdfile.write(f'{buf_group()}\n')
    mdfile.write(lambda_group(name,indexgrp,[state0,state1]))

########### Prepare initial files for production ###########

#Create production directory and add analysis script
if os.path.exists("../production"):
    shutil.rmtree("../production")
os.mkdir("../production")
shutil.copy2("extract.py", "../production")

#Generate pH directories and populate with input files
lam_range = np.linspace(-.1,1.1,13)

for lam in lam_range:
    dirname = f"../production/{lam}"
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
        mdfile.write(f'\nlambda-dynamics-atom-set1-initial-lambda               = {lam}\nlambda-dynamics-atom-set2-initial-lambda               = {lam}\n')