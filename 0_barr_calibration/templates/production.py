########### Import modules ###########

import os
import re

########### Prepare and run  ###########

#Find number of groups
num_groups = 0

with open("md.mdp","r") as mdpfile:
    for line in mdpfile:
        if re.match("^lambda-dynamics-number-lambda-group-types", line):
            num_groups = int(line.split("=")[1].strip())
            break

#Import bash environment variables
try:
    os.environ["OMP_NUM_THREADS"]
except:
    os.system("export OMP_NUM_THREADS=1")

try:
    mpi_ranks = os.environ["NSLOTS"]
except:
    mpi_ranks = 1

#Run simulation
os.system("gmx_mpi grompp -f md.mdp -c npt.gro -p topol.top -o pro.tpr -n index.ndx -maxwarn 3")
os.system(f"mpirun -n {mpi_ranks} gmx_mpi mdrun -deffnm pro -npme 0")

#Extract CpHMD data
os.system(f"gmx_mpi cphmd -s pro.tpr -e pro.edr")