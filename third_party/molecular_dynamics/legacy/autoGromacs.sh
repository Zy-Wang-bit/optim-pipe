#!/bin/bash

# 设置基础参数和文件名
STRUCTURE=$1

#FORCEFIELD="CHAMM36"
#WATER_MODEL="tip3p"
GPU_ID=$2
NTMPI="5"
BOX_TYPE="cubic"
MaxWarn="2"
mkdir "${STRUCTURE%.*}"
cd "${STRUCTURE%.*}"

gmx pdb2gmx -f ../$STRUCTURE -o ${STRUCTURE%.*}2gmx.gro -p topol.top -ignh  <<EOF
8
1
EOF
gmx editconf -f ${STRUCTURE%.*}2gmx.gro -o ${STRUCTURE%.*}_box.gro -d 1.0 -bt cubic
gmx solvate -cp ${STRUCTURE%.*}_box.gro -o ${STRUCTURE%.*}_sol.gro -p topol.top
gmx grompp -f ../../1.minimization.mdp -c ${STRUCTURE%.*}_sol.gro -r ${STRUCTURE%.*}_sol.gro -p topol.top -o em.tpr -maxwarn $MaxWarn
gmx genion -s em.tpr -p topol.top -o system.gro -neutral -conc 0.15  <<EOF
13
EOF
gmx grompp -f ../../1.minimization.mdp -c system.gro -r system.gro -p topol.top -o em.tpr -maxwarn $MaxWarn
gmx mdrun -v -deffnm em -ntomp $NTMPI -gpu_id $GPU_ID -ntmpi 1 
gmx grompp -f ../2.1equilibration_NVT.mdp -c em.gro -r em.gro -p topol.top -o NVT.tpr -maxwarn $MaxWarn
gmx mdrun -v -deffnm NVT -ntomp $NTMPI -gpu_id $GPU_ID -ntmpi 1 
gmx grompp -f ../2.2equilibration_NPT.mdp -c NVT.gro -r NVT.gro -p topol.top -o NPT.tpr -maxwarn $MaxWarn
gmx mdrun -v -deffnm NPT -ntomp $NTMPI -gpu_id $GPU_ID -ntmpi 1 
gmx grompp -f ../3.production.mdp -c NPT.gro -r NVT.gro -p topol.top -o md.tpr
nohup gmx mdrun -s md.tpr -v -deffnm ${STRUCTURE%.*} -ntomp $NTMPI -gpu_id $GPU_ID -ntmpi 1 &