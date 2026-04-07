# optim-pipe

pH 依赖性抗体亲和力优化的自动化计算流水线。
采用 3 Tier 漏斗架构：Tier 1 高通量筛选（MPNN + ESM + FoldX）→ Tier 2 结构评估（PyRosetta + SimpleFold 3x）→ Tier 3 精排（dddG_elec），
通过多轮湿实验闭环（R1→R2→R3…）持续优化抗体变体。

当前聚焦于抗 HBsAg 抗体 **1E62** 的 pH 敏感性改造（酸性解离），覆盖基因型 Ae、B、D1。

---

## 项目进展

| 轮次 | 状态 | 概要 |
|------|------|------|
| **1E62 R1** | 已完成 | 筛选出约13,000序列进行酵母展示筛选，结果未出 |
| **1E62 R2** | 已完成（湿实验验证） | 20 个变体，8 个活性（40% 命中率）；com18 为最优 |
| **1E62 R3** | 进行中 | 基于 R2 SAR 分析的理性组合设计，10 个变体，计算评估阶段 |
| **sdab** | 规划中 | 单域抗体 CDR His scanning 已完成，待 pH 工具就绪后启动组合设计 |

详细路线图见 [ROADMAP.md](docs/ROADMAP.md)。

---

## 环境配置

### Conda 环境

| 组件 | conda 环境名 | 用途 |
|------|-------------|------|
| **主 pipeline** | `optim-pipe` | 核心流程：接口扫描、ESM评分、FoldX调用、筛选合并 |
| ProteinMPNN | `proteinmpnn` | 序列设计（Step 2） |
| SimpleFold | `simplefold` | 结构预测 3x 采样（Tier 2 Step 10） |
| PyRosetta | `pyrosetta` | 突变体建模 + dddG_elec + pH-score（Tier 2 Step 9a/9c） |

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
├── run_pipeline.sh          # 主入口：3 Tier / 13 Step pipeline
│
├── configs/                 # 配置文件
│   ├── config.yaml          #   主配置（路径、设计参数、筛选条件）
│   ├── config_*_n*.yaml     #   各 PDB 模板的专用配置
│   └── mpnn/                #   MPNN 链 ID 和固定位点配置
│
├── scripts/                 # 核心 pipeline 脚本（16 个模块）
│   ├── scan_interface.py    #   Step 1: 接口热点扫描
│   ├── run_mpnn_design.py   #   Step 2: MPNN 序列设计
│   ├── build_his_seeds.py   #   Step 2: His 种子库生成
│   ├── run_esm_chunk.py     #   Step 3: ESM-1b 评分
│   ├── pick_for_foldx.py    #   Step 3: 筛选候选
│   ├── repair_pdbs.py       #   Step 4: FoldX RepairPDB
│   ├── precompute_wt_ac.py  #   Step 5: 预计算 WT 相互作用能
│   ├── make_mutlist_chunk.py#   Step 6: 生成 FoldX 批次
│   ├── run_foldx_batch.py   #   Step 7: FoldX 评估
│   ├── tier1_filter.py      #   Step 8: Tier 1 自适应门槛筛选
│   ├── build_structures.py  #   Step 9a: PyRosetta 突变体建模
│   ├── run_pka.py           #   Step 9b: pKa 预测
│   ├── run_rosetta_eval.py  #   Step 9c: dddG_elec + pH-score
│   ├── run_simplefold_3x.py #   Step 10: SimpleFold 3x 采样
│   ├── run_rmsd.py          #   Step 11: CDR RMSD
│   ├── tier2_filter.py      #   Step 12: Tier 2 pKa + RMSD 筛选
│   ├── merge_and_rank.py    #   Step 13: Tier 3 精排
│   ├── merge_and_select.py  #   旧模式合并筛选（向后兼容）
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
│   └── molecular_dynamics/  #   MD 分析模块（独立）
│       ├── run_md.py            # GROMACS 全自动执行（PDB → 轨迹）
│       ├── analyze_trajectory.py# 轨迹分析入口
│       ├── compare_ph.py        # 双 pH 差异比较
│       ├── configs/             # MD 模块配置
│       ├── mdp_templates/       # Jinja2 参数模板
│       └── lib/                 # 分析器（RMSD/RMSF/氢键/SASA/接触/MMPBSA）
│
├── archive/                 # 归档（旧脚本与备份）
└── workspace/               # FoldX 工作区（待确认）
```

---

## Pipeline 流程（3 Tier / 13 Step）

```
Tier 1: 生成 + 高通量筛选                        ~10万 → ~500
  Step 1:  scan_interface.py      接口热点扫描
  Step 2:  run_mpnn_design.py     MPNN 序列设计 + His 种子库
  Step 3:  run_esm_chunk.py       ESM-1b 评分
  Step 4:  repair_pdbs.py         FoldX RepairPDB
  Step 5:  precompute_wt_ac.py    预计算 WT 相互作用能
  Step 6:  make_mutlist_chunk.py  生成 FoldX 突变批次
  Step 7:  run_foldx_batch.py     FoldX BuildModel + AnalyseComplex（双 pH）
  Step 8:  tier1_filter.py        dG + delta 自适应门槛 → tier1_candidates.csv

Tier 2: 结构评估 + 筛选（tier2.enabled 控制）       ~500 → ~100
  PyRosetta 线（复合物模板）：
    Step 9a: build_structures.py  PyRosetta 突变体建模
    Step 9b: run_pka.py           PROPKA3 + pKAI+ pKa 预测
    Step 9c: run_rosetta_eval.py  dddG_elec + pH-score
  SimpleFold 3x 线（单抗体序列）：
    Step 10: run_simplefold_3x.py SimpleFold 3B × 3 采样
    Step 11: run_rmsd.py          CDR RMSD（3x 中位数 + 异常值剔除）
  Step 12: tier2_filter.py        pKa 相对排序 + CDR RMSD 门槛

Tier 3: 精排                                     ~100 → top-N
  Step 13: merge_and_rank.py      dddG_elec 排序 + 软标记 → final_candidates.csv
```

`tier1.enabled: false` 时回退到旧 8 步流程（merge_and_select.py）。

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
| **tier1** | Tier 1 自适应筛选（目标范围、阈值、ESM 异常标记） |
| **tier2** | Tier 2 结构评估（PyRosetta/SimpleFold 配置、pKa/RMSD 筛选阈值） |
| **tier3** | Tier 3 精排（排序指标、软标记） |
| **selection** | 旧模式自适应筛选（向后兼容） |

---

## MD 分析模块

独立于 3-Tier pipeline 的分子动力学模块，提供 GROMACS 自动化执行和轨迹分析能力。

### 功能

- **MD 执行自动化**：PDB → 轨迹全自动流程，解决了 pdb2gmx 手动选择力场/水模型的问题（使用 `-ff`/`-water` 名称指定）
- **CpHMD 支持**：恒定 pH 分子动力学（λ-dynamics），His 质子化态在模拟中动态切换
- **6 项分析**：RMSD（含收敛检测）、RMSF（残基柔性）、界面氢键占有率、Buried SASA、界面接触数、MM-GBSA（可选）
- **双 pH 比较**：同一变体在 pH 7.4 vs 6.0 下的行为差异分析
- **汇总报告**：所有变体分析结果汇总为 `md_report.csv`

### 快速开始

```bash
module load gromacs/2024.2 openmpi/5.0.3
cd third_party/molecular_dynamics

# 运行 MD
python run_md.py --pdb /path/to/complex.pdb --ph 7.4 6.0 \
    --output-dir experiments/1E62/R3/md/

# 分析轨迹
python analyze_trajectory.py --traj experiments/1E62/R3/md/variant/pH_7.4/

# pH 比较
python compare_ph.py --variant-dir experiments/1E62/R3/md/variant/ \
    --base-ph 7.4 --target-ph 6.0
```

### 配置

MD 模块使用独立配置文件 `third_party/molecular_dynamics/configs/md_config.yaml`，包含 GROMACS 参数、模拟阶段设置、CpHMD 配置和分析参数。
