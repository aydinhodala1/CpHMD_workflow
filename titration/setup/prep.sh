#!/bin/bash --login

#$ -N titr_prep_REP_VAR
#$ -cwd
#$ -pe smp.pe 16

export OMP_NUM_THREADS=1

num_ola=20

#Load modules for generating initial system
module load apps/intel-17.0/obabel/2.4.1
module load apps/intel-17.0/packmol/18.169

#Generate dry oleic acid structure and begin topology file
obabel ola.sdf -O ola.pdb
sed s/NUM_VAR/$num_ola/g template.top > topol.top
sed s/NUM_VAR/$num_ola/g template.inp > sphere.inp
packmol < sphere.inp

#Reset envirnoment and load Gromacs
module purge
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

#Generate index and add titration group
echo q | gmx_mpi "make_ndx -f ola_sphere_npt.gro -o ola_sphere.ndx"
rep=0
while [ $rep -le $(($num_ola-1)) ]
do
        a=$(($rep*54+1))
        b=$(($rep*54+2))
        c=$(($rep*54+3))
        d=$(($rep*54+4))
        e=$(($rep*54+21))
        sed -i '$a\[lam'"$(($rep+1))"']' ola_sphere.ndx
        sed -i '$a\'"$a"' '"$b"' '"$c"' '"$d"' '"$e" ola_sphere.ndx
        rep=$(($rep+1))
done

#Add required CpHMD parameters to md.mdp for production run
sed s/GRP_VAR/$(($num_ola+1))/g md.mdp > run.mdp
sed -i s/BUF_NUM_VAR/$((10*$num_ola))/g run.mdp
rep=2
while [ $rep -le $(($num_ola+1)) ]
do
        sed -i '$a\ ' run.mdp
        sed -i '$a\lambda-dynamics-group-type'"$rep"'-name                          = OLA'"$(($rep-1))"'' run.mdp
        sed -i '$a\lambda-dynamics-group-type'"$rep"'-n-states                      = 1' run.mdp
        sed -i '$a\lambda-dynamics-group-type'"$rep"'-state-0-charges               =  0.75 -0.21 -0.55 -0.61 0.44' run.mdp
        sed -i '$a\lambda-dynamics-group-type'"$rep"'-state-1-charges               =  0.62 -0.28 -0.76 -0.76 0.00' run.mdp
        sed -i '$a\lambda-dynamics-group-type'"$rep"'-state-1-reference-pka         = 4.75' run.mdp
        sed -i '$a\lambda-dynamics-group-type'"$rep"'-state-1-dvdl-coefficients     =  -226.21773756   952.8559738  -1529.0831676    967.60375464  -251.28863742  -571.62201758    87.69457099' run.mdp
        sed -i '$a\ ' run.mdp
        sed -i '$a\lambda-dynamics-atom-set'"$rep"'-name                         = OLA'"$(($rep-1))"'' run.mdp
        sed -i '$a\lambda-dynamics-atom-set'"$rep"'-index-group-name             = lam'"$(($rep-1))"'' run.mdp
        sed -i '$a\lambda-dynamics-atom-set'"$rep"'-barrier                      = 5.0' run.mdp
        sed -i '$a\lambda-dynamics-atom-set'"$rep"'-initial-lambda               = LAM_OLA_VAR' run.mdp
        sed -i '$a\lambda-dynamics-atom-set'"$rep"'-charge-restraint-group-index = 1' run.mdp        
        rep=$(($rep+1))
done
