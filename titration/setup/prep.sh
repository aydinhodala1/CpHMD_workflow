#!/bin/bash --login

#$ -N titr_prep_REP_VAR
#$ -cwd
#$ -pe smp.pe 16

export OMP_NUM_THREADS=1
#Reset envirnoment and load Gromacs
module load apps/singularity/constantph/9bb0cc2f_mpi

#Expand box for solvation
gmx_mpi "editconf -f ola_sphere.pdb -o ola_sphere_lrg.gro -c -d 2.0 -bt cubic"

#Solvate oleic acid
gmx_mpi "solvate -cp ola_sphere_lrg.gro -cs spc216.gro -o ola_sphere_wet.gro -p topol.top"
sed -i '/#include ".\/charmm36-mar2019-cphmd.ff\/forcefield.itp"/ a\#include ".\/charmm36-mar2019-cphmd.ff\/tip3p.itp"' topol.top

#Replace water with buffer particles and ions
gmx_mpi "grompp -f buf.mdp -c ola_sphere_wet.gro -p topol.top -o ola_sphere_buf.tpr"
echo SOL |gmx_mpi "genion -s ola_sphere_buf.tpr -o ola_sphere_buf.gro -p topol.top -pname BUF -np $((10*$num_ola))"

gmx_mpi "grompp -f buf.mdp -c ola_sphere_buf.gro -p topol.top -o ola_sphere_ion.tpr"
echo SOL |gmx_mpi "genion -s ola_sphere_ion.tpr -o ola_sphere_ion.gro -p topol.top -pname NA -np IONS_VAR -nnmame CL -nn IONS_VAR"

sed -i '/#include ".\/charmm36-mar2019-cphmd.ff\/forcefield.itp"/ a\#include ".\/charmm36-mar2019-cphmd.ff\/ions.itp"' topol.top

#Minimize and pre-equilbirate system without any CpHMD
gmx_mpi "grompp -f minim.mdp -c ola_sphere_buf.gro -p topol.top -o ola_sphere_min.tpr -maxwarn 2"
mpirun -n $NSLOTS gmx_mpi "mdrun -deffnm ola_sphere_min -npme 0"

gmx_mpi "grompp -f nvt.mdp -c ola_sphere_min.gro -r ola_sphere_min.gro -p topol.top -o ola_sphere_nvt.tpr -maxwarn 2"
mpirun -n $NSLOTS gmx_mpi "mdrun -deffnm ola_sphere_nvt -npme 0"

gmx_mpi "grompp -f npt.mdp -c ola_sphere_nvt.gro -r ola_sphere_nvt.gro -t ola_sphere_nvt.cpt -p topol.top -o ola_sphere_npt.tpr -maxwarn 2"
mpirun -n $NSLOTS gmx_mpi "mdrun -deffnm ola_sphere_npt -npme 0"