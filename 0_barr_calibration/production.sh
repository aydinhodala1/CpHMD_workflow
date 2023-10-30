#!/bin/bash --login

#$ -N ola_NAME_VAR_REP_VAR
#$ -hold_jid ola_NAME_VAR_prep_REP_VAR
#$ -t 1-14
#$ -cwd
#$ -pe smp.pe 16

export OMP_NUM_THREADS=1

pH="4.75"
init_lam_buf=0.5
init_lam_ola=0.0
dirname=$SGE_TASK_ID

#Load GROMACS singularity
module load apps/singularity/constantph/9bb0cc2f_mpi

#Generate files and directory for simulation
mkdir ./$dirname
cp -r setup/charmm36-mar2019-cphmd.ff/ ./$dirname
cp setup/ola_sphere.ndx ./$dirname
cp setup/ola_sphere.top ./$dirname
cp setup/ola_sphere_npt.gro ./$dirname
sed s/PH_VAR/$pH/g setup/run.mdp > ./$dirname/run.mdp
sed -i s/LAM_BUF_VAR/$init_lam_buf/g ./$dirname/run.mdp
sed -i s/LAM_OLA_VAR/$init_lam_ola/g ./$dirname/run.mdp

#Run simulation and gather CpHMD data
cd ./$dirname
gmx_mpi "grompp -f run.mdp -c ola_sphere_npt.gro -p ola_sphere.top -o pro$dirname\.tpr -n ola_sphere.ndx -maxwarn 3"
mpirun -n $NSLOTS gmx_mpi "mdrun -deffnm pro$dirname -npme 0"
gmx_mpi "cphmd -s pro$dirname\.tpr -e pro$dirname\.edr -o d$dirname\.xvg -dvdl yes"
