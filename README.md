# optim-pipe

pH 依赖性抗体亲和力优化的自动化计算流水线。
集成 ProteinMPNN（序列设计）→ ESM-1b（快速评分）→ FoldX（能量计算）→ 自适应筛选，
并通过多轮湿实验闭环（R1→R2→R3…）持续优化抗体变体。

当前聚焦于抗 HBsAg 抗体 **1E62** 的 pH 敏感性改造（酸性解离），覆盖基因型 Ae、B、D1。

---

## 项目进展

| 轮次 | 状态 | 概要 |
|------|------|------|
| **1E62 R1** | 已完成 | 首轮探索，数据有限 |
| **1E62 R2** | 已完成（湿实验验证） | 20 个变体，8 个活性（40% 命中率）；com18 为最优（广谱结合 + 强酸性解离） |
| **1E62 R3** | 进行中 | 基于 R2 SAR 分析的理性组合设计，10 个变体，计算评估阶段 |
| **sdab** | 规划中 | 单域抗体 CDR His scanning 已完成，待 pH 工具就绪后启动组合设计 |

详细路线图见 [ROADMAP.md](docs/ROADMAP.md)。

---

## 环境配置

### Conda 环境

| 组件 | conda 环境名 | 用途 |
|------|-------------|------|
| **主 pipeline** | `optim-pipe` | 核心流程：接口扫描、ESM评分、FoldX调用、合并筛选 |
| ProteinMPNN | `proteinmpnn` | 序列设计（Step 2） |
| SimpleFold | `simplefold` | 结构预测（R3 评估） |
| PyRosetta | `pyrosetta` | Rosetta 分析（pH-score、dddG_elec） |

**激活主环境：**
```bash
conda activate optim-pipe
```

**主要依赖（optim-pipe 环境）：**
- Python 3.11
- PyTorch 2.8
- fair-esm 2.0
- BioPython 1.85
- pandas, numpy, scipy, PyYAML

> **注意：** 项目根目录下的 `.venv/` 是一个空的 venv 虚拟环境（未安装包），实际应使用上述 conda 环境。

### FoldX

FoldX 二进制位于 `third_party/foldx/foldx`，由 `configs/config.yaml` 中的 `paths.foldx_bin` 配置。

---

## 快速开始

```bash
conda activate optim-pipe

# 运行完整 pipeline（使用默认配置 configs/config.yaml）
bash run_pipeline.sh

# 使用指定配置
bash run_pipeline.sh configs/config_foldx_n1-0.yaml

# 清理旧产物后重新运行
CLEAN=1 bash run_pipeline.sh
```

---

## 目录结构

```
optim-pipe/
│
├── README.md                # 项目总览
├── CLAUDE.md                # AI 工作协议 + 跨对话任务管理
├── run_pipeline.sh          # 主入口：8 步 pipeline 一键运行
│
├── configs/                 # 配置文件
│   ├── config.yaml          #   主配置（路径、设计参数、筛选条件）
│   ├── config_*_n*.yaml     #   各 PDB 模板的专用配置
│   └── mpnn/                #   MPNN 链 ID 和固定位点配置
│
├── scripts/                 # 核心 pipeline 脚本（12 个模块）
│   ├── scan_interface.py    #   Step 1: 接口热点扫描
│   ├── run_mpnn_design.py   #   Step 2: MPNN 序列设计
│   ├── build_his_seeds.py   #   Step 2: His 种子库生成
│   ├── run_esm_chunk.py     #   Step 3: ESM-1b 评分
│   ├── pick_for_foldx.py    #   Step 3: 筛选候选
│   ├── repair_pdbs.py       #   Step 4: FoldX RepairPDB
│   ├── precompute_wt_ac.py  #   Step 5: 预计算 WT 相互作用能
│   ├── make_mutlist_chunk.py#   Step 6: 生成 FoldX 批次
│   ├── run_foldx_batch.py   #   Step 7: FoldX 评估
│   ├── merge_and_select.py  #   Step 8: 合并与自适应筛选
│   ├── mpnn_filter_his.py   #   辅助：过滤 His 突变
│   └── build_mutants_from_csv.py  # 辅助：从 CSV 构建突变体序列
│
├── analysis/                # 后分析模块
│   ├── pka/                 #   pKa 预测（PROPKA3 + pKAI+ 双验证）
│   ├── rosetta/             #   Rosetta pH-score 与 dddG_elec 计算
│   ├── structure_compare/   #   RMSD 结构比较
│   └── elisa_vs_computation.py  # ELISA vs 计算指标关联分析
│
├── experiments/             # 实验轮次数据与分析
│   ├── 1E62/                #   1E62 抗体系统（主线）
│   │   ├── data/            #     共享输入数据（PDB 结构、WT 序列）
│   │   ├── R1/              #     R1 实验数据
│   │   ├── R2/              #     R2 实验数据（含 wet_lab/ ELISA 结果）
│   │   └── R3/              #     R3 理性设计（含 design/, scripts/, results/）
│   ├── sdab/                #   单域抗体数据
│   ├── 3E72/                #   3E72 抗体系统（搁置中）
│   ├── 16F9/                #   16F9 抗体（搁置中）
│   └── ext_gpt/             #   GPT 辅助设计
│
├── docs/                    # 文档
│   ├── ROADMAP.md           #   后续工作路线图
│   ├── design/              #   设计文档
│   └── reports/             #   分析报告
│
├── third_party/             # 外部工具与依赖
│   ├── ProteinMPNN/         #   序列设计
│   ├── ml-simplefold/       #   结构预测（含 weights/ 模型权重）
│   ├── foldx/               #   FoldX 二进制
│   └── molecular_dynamics/  #   GROMACS 分子动力学脚本
│
├── archive/                 # 归档（旧脚本与备份）
└── workspace/               # FoldX 工作区（待确认）
```

---

## Pipeline 流程

```
Step 1: scan_interface.py     → 接口热点扫描（His 偏向位点识别）
Step 2: run_mpnn_design.py    → MPNN 序列设计 + His 种子库
Step 3: run_esm_chunk.py      → ESM-1b 评分与候选合并
Step 4: repair_pdbs.py        → FoldX RepairPDB
Step 5: precompute_wt_ac.py   → 预计算 WT 相互作用能
Step 6: make_mutlist_chunk.py → 生成 FoldX 突变批次
Step 7: run_foldx_batch.py    → FoldX BuildModel + AnalyseComplex（双 pH）
Step 8: merge_and_select.py   → 合并结果 + 自适应筛选 → final_top10k.csv
```

### R3 评估流程（理性设计专用）

R3 不经过 MPNN 重新设计，而是基于 R2 SAR 规则进行组合设计，采用独立评估流水线：

```
1. 结构生成    → PyRosetta 突变建模 + SimpleFold 集成预测
2. pKa 预测    → PROPKA3 + pKAI+ 双工具共识
3. 能量评估    → FoldX 双 pH（7.4 / 6.0）
4. pH-score    → Rosetta His 氢键网络分析
5. RMSD 计算   → CDR 区结构稳定性验证
6. 综合排名    → 多指标加权打分 → r3_evaluation_report.csv
```

---

## 配置说明

主配置文件 `configs/config.yaml` 包含以下部分：

| 配置块 | 说明 |
|--------|------|
| **paths** | 各工具和数据的路径 |
| **design** | 设计区域（链、残基窗口） |
| **interface** | 抗体-抗原接口定义（链定义、距离截断） |
| **his_bias** | His 突变偏好（最少 His 数、热点位点、扫描限制） |
| **screening** | ESM 初筛参数（保留比例、序列长度约束） |
| **foldx** | FoldX 计算参数（pH 点 [7.4, 6.0]、运行次数） |
| **resources** | 并行计算资源（批量大小、最大进程数、checkpoint 间隔） |
| **selection** | 自适应筛选（阈值、松弛策略、配额约束） |
