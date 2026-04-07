#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# 设置基础参数和文件名
INPUT_PDB_PATH=$1 # 例如 pdb/rank_3.pdb 或 rank_3.pdb
GPU_ID=$2
NTMPI="5"
# BOX_TYPE="cubic" # 未在脚本中使用，可以注释或删除
MaxWarn="2"

# 从输入路径中提取不含扩展名的基本名称，用于文件命名
# 例如，如果 INPUT_PDB_PATH 是 "pdb/rank_3.pdb"，PDB_BASENAME 是 "rank_3"
# 如果 INPUT_PDB_PATH 是 "rank_3.pdb"，PDB_BASENAME 也是 "rank_3"
PDB_BASENAME=$(basename "${INPUT_PDB_PATH%.pdb}")

# 创建工作目录的名称，可能包含路径
# 例如，如果 INPUT_PDB_PATH 是 "pdb/rank_3.pdb"，WORK_DIR 是 "pdb/rank_3"
# 如果 INPUT_PDB_PATH 是 "rank_3.pdb"，WORK_DIR 是 "rank_3"
WORK_DIR="${INPUT_PDB_PATH%.pdb}"

# 创建工作目录并进入
mkdir -p "$WORK_DIR" # 使用 -p 确保如果路径不存在能被创建
cd "$WORK_DIR"

# 原始PDB文件的名称 (不含路径)
ORIGINAL_PDB_FILENAME=$(basename "$INPUT_PDB_PATH")

# pdb2gmx 需要的 PDB 文件路径 (相对于当前工作目录 WORK_DIR)
# 因为我们 cd 进入了 WORK_DIR，所以原始 PDB 文件在其上一级目录
PDB_FOR_PDB2GMX="../$ORIGINAL_PDB_FILENAME"

echo "Current working directory: $(pwd)"
echo "Using PDB file: $PDB_FOR_PDB2GMX (resolved: $(realpath $PDB_FOR_PDB2GMX 2>/dev/null || echo "not found"))"
echo "Output file prefix: $PDB_BASENAME"

# 检查MDP文件是否存在 (相对于当前工作目录 WORK_DIR)
# 这些路径需要根据你实际存放MDP文件的位置来调整
MDP_MINIM="../../1.minimization.mdp"
MDP_NVT="../../2.1equilibration_NVT.mdp"
MDP_NPT="../../2.2equilibration_NPT.mdp"
MDP_PROD="../../3.production.mdp"

# 确保MDP文件存在
for mdp_file in "$MDP_MINIM" "$MDP_NVT" "$MDP_NPT" "$MDP_PROD"; do
    if [ ! -f "$mdp_file" ]; then
        echo "Error: MDP file $mdp_file (resolved: $(realpath $mdp_file 2>/dev/null || echo "not found in expected location")) not found!"
        exit 1
    fi
done

echo -e "8\n1\n1\n0\n0\n0\n0\n0\n0\n0\n0\n0" | gmx pdb2gmx -f "$PDB_FOR_PDB2GMX" -o "${PDB_BASENAME}_processed.gro" -p topol.top -ignh -ter # Fro GroFlow PDB
# echo -e "8\n1\n0\n0\n0\n0\n1\n0" | gmx pdb2gmx -f "$PDB_FOR_PDB2GMX" -o "${PDB_BASENAME}_processed.gro" -p topol.top -ignh -ter # Fro AF2 PDB

# gmx pdb2gmx -f "$PDB_FOR_PDB2GMX" -o "${PDB_BASENAME}_processed.gro" -p topol.top -ignh -ter <<EOF
# 8
# 1
# # for AF2 PDB: 000010
# 0
# 0
# 0
# 0
# 0
# 0
# 1
# 0
# 1
# 0
# EOF

gmx editconf -f "${PDB_BASENAME}_processed.gro" -o "${PDB_BASENAME}_box.gro" -d 1.0 -bt cubic
gmx solvate -cp "${PDB_BASENAME}_box.gro" -o "${PDB_BASENAME}_sol.gro" -p topol.top
gmx grompp -f "$MDP_MINIM" -c "${PDB_BASENAME}_sol.gro" -r "${PDB_BASENAME}_sol.gro" -p topol.top -o em.tpr -maxwarn $MaxWarn
gmx genion -s em.tpr -p topol.top -o "${PDB_BASENAME}_system_ions.gro" -neutral -conc 0.15 <<EOF
13
EOF
# 注意：genion 的输出文件名改为了 PDB_BASENAME 相关，以保持一致性
# 下一步grompp使用genion的输出
gmx grompp -f "$MDP_MINIM" -c "${PDB_BASENAME}_system_ions.gro" -r "${PDB_BASENAME}_system_ions.gro" -p topol.top -o em.tpr -maxwarn $MaxWarn
gmx mdrun -v -deffnm em -ntomp $NTMPI -gpu_id $GPU_ID -ntmpi 1
# NVT平衡
gmx grompp -f "$MDP_NVT" -c em.gro -r em.gro -p topol.top -o NVT.tpr -maxwarn $MaxWarn
gmx mdrun -v -deffnm NVT -ntomp $NTMPI -gpu_id $GPU_ID -ntmpi 1
# NPT平衡
gmx grompp -f "$MDP_NPT" -c NVT.gro -r NVT.gro -p topol.top -o NPT.tpr -maxwarn $MaxWarn
gmx mdrun -v -deffnm NPT -ntomp $NTMPI -gpu_id $GPU_ID -ntmpi 1
# 生产阶段
# 注意: -r NVT.gro 在这里。通常，对于成品MD，如果需要位置限制，会用 -r NPT.gro。
# 但如果你确实希望参考NVT的结构进行限制，则保持不变。
gmx grompp -f "$MDP_PROD" -c NPT.gro -r NVT.gro -p topol.top -o md.tpr -maxwarn $MaxWarn # 添加了 -maxwarn
nohup gmx mdrun -s md.tpr -v -deffnm "$PDB_BASENAME" -ntomp $NTMPI -gpu_id $GPU_ID -ntmpi 1 &

echo "MD simulation for $PDB_BASENAME submitted to background."
echo "Check nohup.out and ${PDB_BASENAME}.log for progress."
