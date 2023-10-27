########### Import modules ###########

import os
import re
import matplotlib.pyplot as plt
import numpy as np

########### Prepare and run  ###########

#Find number of groups
num_groups = 0

with open("topol.top","r") as topolfile:
    while num_groups == 0:
        for line in topolfile:
            if re.match("^lambda-dynamics-number-lambda-group-types", line):
                num_groups = int(line.split("=")[1].strip())

#Import bash environment variables
omp_slots = os.environ["OMP_NUM_THREADS"]
mpi_ranks = os.environ["NSLOTS"]

#Run simulation
os.system("gmx_mpi grompp -f md.mdp -c npt.gro -p topol.top -o pro.tpr -n index.ndx -maxwarn 3")
os.system(f"mpirun -n {mpi_ranks} gmx_mpi mdrun -deffnm pro -npme 0")

#Extract CpHMD data
os.system(f"gmx_mpi cphmd -s pro.tpr -e pro.edr -numplot {num_groups}")

########### Raw data analysis ###########

eq_time = 5_000
data_dictionary = {}

time,buff,*acids = np.loadtxt("pro.edr", comments=["#","@"], unpack=True)

rounded = []
eq_steps = np.where(time == eq_time)[0][0]
for acid in acids:
    plt.plot(time,acid)
plt.xlabel("Time /ps")
plt.ylabel("$\\lambda$")
plt.savefig(f"lambdas.png")