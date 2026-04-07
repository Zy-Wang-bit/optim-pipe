# MD 分析模块设计

## 背景

当前 pipeline 中分子动力学（MD）模拟仅用于事后验证少量最终候选和 AF3 预测结构，分析只提取 RMSD。80ns 轨迹中蕴含的结构稳定性、界面完整性、pH 敏感性信息未被充分利用。

同时，现有 MD 执行流程依赖手动操作（pdb2gmx 力场/水模型选择因 PDB 不同而索引不一致），难以批量运行。

## 目标

设计一个**独立 MD 模块**（不集成到 3-Tier pipeline），提供：
1. **MD 执行自动化**：PDB → 轨迹的全自动流程，支持 CpHMD（恒定 pH 分子动力学）
2. **轨迹分析工具集**：提取结构稳定性、界面完整性、pH 敏感性、结合亲和力指标
3. **双 pH 比较**：同一变体在 pH 7.4 vs 6.0 下的行为差异分析

## 模块结构

```
third_party/molecular_dynamics/
├── configs/
│   └── md_config.yaml          # MD 模块配置
├── mdp_templates/              # GROMACS 参数模板（Jinja2，时长可配）
│   ├── minimization.mdp.j2
│   ├── nvt.mdp.j2
│   ├── npt.mdp.j2
│   ├── production.mdp.j2
│   └── cphmd.mdp.j2           # CpHMD 专用参数
├── run_md.py                   # 执行入口：PDB → 轨迹
├── analyze_trajectory.py       # 分析入口：轨迹 → 指标报告
├── compare_ph.py               # 双 pH 差异比较
├── lib/
│   ├── gromacs_wrapper.py      # GROMACS 命令封装
│   ├── protonation.py          # His 质子化态管理（CpHMD 配置）
│   ├── analyzers/
│   │   ├── rmsd.py             # RMSD 时间序列 + 收敛判断
│   │   ├── rmsf.py             # 残基柔性
│   │   ├── hbond.py            # 界面氢键占有率
│   │   ├── sasa.py             # 溶剂可及表面积
│   │   ├── contacts.py         # 界面接触数
│   │   └── mmpbsa.py           # MM-PBSA/GBSA（可选）
│   └── report.py               # 汇总报告生成
└── legacy/                     # 现有脚本归档
    ├── autoGromacs.sh
    ├── true-autoGromacs.sh
    └── show_xvg.py
```

---

## 一、MD 执行自动化（run_md.py）

### 解决 pdb2gmx 手动选择问题

不再用 `echo` 管道选索引，改用命令行参数按名称指定：

```bash
gmx pdb2gmx -f input.pdb -ff oplsaam -water tip3p -ignh
```

`-ff` 直接指定力场目录名，`-water` 直接指定水模型名称，不依赖交互式索引。

### CpHMD（恒定 pH 分子动力学）

使用 GROMACS λ-dynamics CpHMD：

1. **可滴定残基配置**：从 `md_config.yaml` 读取参与 CpHMD 的 His 列表，或自动检测 PDB 中所有 His
2. **λ 坐标初始化**：`gmx cphmd` 生成 λ 粒子初始配置
3. **MDP 参数**：`cphmd.mdp.j2` 包含 CpHMD 参数（`cphmd-ph`、λ 更新频率等）
4. **输出**：标准轨迹 + λ 轨迹文件（记录质子化态随时间变化）

CpHMD 的优势：His 的质子化态在模拟中动态切换，无需预先假设 pH 6.0 时 His 是否质子化——这由模拟本身基于 His 在蛋白环境中的实际 pKa 决定。

### 执行流程

```
PDB 输入
  │
  ├─ pdb2gmx（-ff/-water 名称指定，CpHMD 质子化配置）
  ├─ editconf（box 创建）
  ├─ solvate（加水）
  ├─ grompp + genion（加离子，0.15M）
  ├─ minimization
  ├─ NVT equilibration
  ├─ NPT equilibration
  └─ production（CpHMD 模式或标准 MD）
  │
  ▼
轨迹文件 (.xtc + .tpr) + λ 轨迹（CpHMD）
```

### 命令行接口

```bash
# 单 PDB，单 pH（输出到默认目录）
python run_md.py --pdb input.pdb --ph 6.0

# 单 PDB，双 pH，指定输出目录
python run_md.py --pdb input.pdb --ph 7.4 6.0 --output-dir experiments/1E62/R3/md/

# 批量执行
python run_md.py --pdb-dir variants/ --ph 7.4 6.0 --output-dir experiments/1E62/R3/md/
```

### 输出组织

输出目录**可配置**，通过 `--output-dir` 命令行参数或 `md_config.yaml` 的 `output_dir` 字段指定。建议放到对应的实验目录下（如 `experiments/<system>/<round>/md/`），但不硬编码。

```
experiments/1E62/R3/md/           # 或用户指定的任意路径
├── HE1H/
│   ├── pH_7.4/
│   │   ├── topol.tpr
│   │   ├── production.xtc
│   │   ├── production.edr
│   │   └── lambda.xvg          # CpHMD λ 轨迹
│   └── pH_6.0/
│       └── ...
└── WT/
    └── ...
```

### md_config.yaml

```yaml
md:
  force_field: "oplsaam"
  water_model: "tip3p"
  box_type: "cubic"
  box_distance: 1.0              # nm
  ion_concentration: 0.15        # mol/L
  temperature: 310               # K

  stages:
    minimization:
      max_steps: 5000
    nvt:
      duration_ps: 200
    npt:
      duration_ps: 1000
    production:
      duration_ns: 50             # 可按需调整
      dt_ps: 0.002
      save_interval_ps: 100

  cphmd:
    enabled: true
    ph_values: [7.4, 6.0]
    titratable_residues: []       # 空 = 自动检测所有 His
    lambda_update_interval: 100   # steps

  gpu:
    device_id: 0
    ntmpi: 1
    ntomp: 8

  output_dir: "experiments/1E62/R3/md"  # 可配置，建议放到对应实验目录
```

---

## 二、轨迹分析工具集（analyze_trajectory.py）

### 统一接口

```bash
# 运行所有分析项
python analyze_trajectory.py --traj experiments/1E62/R3/md/HE1H/pH_7.4/

# 指定分析项
python analyze_trajectory.py --traj experiments/1E62/R3/md/HE1H/pH_7.4/ --analyses rmsd rmsf hbond

# 双 pH 比较
python compare_ph.py --variant-dir experiments/1E62/R3/md/HE1H/ --base-ph 7.4 --target-ph 6.0
```

### 分析项详细设计

#### 1. RMSD 时间序列（rmsd.py）

- **计算**：全局 Cα RMSD + CDR 区域 RMSD（CDR 定义复用 pipeline config）
- **额外功能**：自动判断收敛时间（sliding window 方差低于阈值的起始帧）
- **输出**：`rmsd_timeseries.csv`（time, global, H1, H2, H3, L1, L2, L3）+ 收敛起始时间
- **工具**：MDAnalysis

#### 2. RMSF 残基柔性（rmsf.py）

- **计算**：Per-residue Cα RMSF（取收敛后区间）
- **关注区域**：CDR 残基 + 突变位点 + 界面残基
- **输出**：`rmsf_per_residue.csv`（chain, resid, resname, rmsf）+ CDR/突变位点平均 RMSF
- **工具**：MDAnalysis

#### 3. 界面氢键占有率（hbond.py）

- **计算**：抗体-抗原链间氢键，统计每对供体-受体的占有率（出现帧数/总帧数）
- **链定义**：从 config 读取抗体链和抗原链 ID
- **输出**：`hbond_occupancy.csv`（donor_chain, donor_resid, acceptor_chain, acceptor_resid, occupancy）+ 总界面氢键数时间序列
- **工具**：MDAnalysis HydrogenBondAnalysis

#### 4. 界面 SASA（sasa.py）

- **计算**：Buried SASA = SASA(抗体) + SASA(抗原) - SASA(复合物)
- **意义**：Buried SASA 随时间减小 → 界面松散/解离趋势
- **输出**：`sasa_timeseries.csv`（time, sasa_ab, sasa_ag, sasa_complex, buried_sasa）
- **工具**：MDAnalysis 或 `gmx sasa`

#### 5. 界面接触数（contacts.py）

- **计算**：抗体-抗原原子间距离 < 阈值（默认 4.5Å）的接触对数
- **输出**：`contacts_timeseries.csv`（time, n_contacts）+ 残基对接触频率矩阵
- **工具**：MDAnalysis Contacts

#### 6. MM-GBSA 结合自由能（mmpbsa.py，可选）

- **计算**：ΔG_bind = E_complex - E_antibody - E_antigen（VdW + 电静力 + 极性溶剂化 + 非极性溶剂化）
- **帧选择**：取收敛后的最后 N 帧（默认从收敛点到轨迹末尾，每 100ps 取 1 帧）
- **输出**：`mmpbsa_summary.csv`（ΔG_bind 均值 ± 标准差）+ per-residue 能量分解
- **工具**：`gmx_MMPBSA`（需额外安装）
- **标记**：可选分析项，config 中 `mmpbsa.enabled: false` 默认关闭

---

## 三、双 pH 比较（compare_ph.py）

输入同一变体的两条轨迹分析结果，计算差异：

| 比较指标 | 计算 | 物理意义 |
|---------|------|---------|
| ΔRMSD | mean_RMSD(pH6.0) - mean_RMSD(pH7.4) | 酸性下结构更不稳定？ |
| ΔRMSF_CDR | mean_RMSF_CDR(pH6.0) - mean_RMSF_CDR(pH7.4) | 酸性下 CDR 更柔软？ |
| ΔHbond | mean_N_hbond(pH6.0) - mean_N_hbond(pH7.4) | 酸性下界面氢键减少？ |
| ΔBuried_SASA | mean_BSASA(pH6.0) - mean_BSASA(pH7.4) | 酸性下界面松散？ |
| ΔContacts | mean_N_contacts(pH6.0) - mean_N_contacts(pH7.4) | 酸性下接触减少？ |
| ΔΔG_bind | ΔG(pH6.0) - ΔG(pH7.4) | 酸性下结合变弱多少？ |

CpHMD 额外指标：**His 质子化态占有率**（从 λ 轨迹提取），直接反映每个 His 在目标 pH 下的实际质子化程度。

---

## 四、汇总报告（report.py）

将所有分析结果汇总为 `md_report.csv`，每个变体一行：

```
variant_id, ph, rmsd_mean, rmsd_converge_ns, rmsf_cdr_mean, rmsf_mutation_site,
n_hbond_mean, buried_sasa_mean, n_contacts_mean, dG_bind_mean, dG_bind_std,
his_protonation_occupancy
```

如有双 pH 数据，附加 Δ 列（delta_rmsd, delta_hbond, ...）。

---

## 五、依赖

| 工具 | 用途 | 安装方式 |
|------|------|---------|
| GROMACS (≥2023) | MD 执行 + CpHMD | 已有 |
| MDAnalysis | 轨迹分析（RMSD/RMSF/氢键/SASA/接触） | `pip install MDAnalysis` |
| Jinja2 | MDP 模板渲染 | `pip install Jinja2` |
| gmx_MMPBSA | MM-PBSA/GBSA（可选） | `pip install gmx_MMPBSA` |

---

## 六、分析配置

在 `md_config.yaml` 中新增分析段：

```yaml
analysis:
  # 链定义（用于界面分析）
  antibody_chains: ["A", "B"]    # 重链 + 轻链
  antigen_chains: ["C"]

  # CDR 定义（复用 pipeline config 或单独定义）
  cdr_regions:
    H1: { chain: "A", start: 26, end: 35 }
    H2: { chain: "A", start: 50, end: 65 }
    H3: { chain: "A", start: 95, end: 102 }
    L1: { chain: "B", start: 24, end: 34 }
    L2: { chain: "B", start: 50, end: 56 }
    L3: { chain: "B", start: 89, end: 97 }

  # 收敛判断
  convergence:
    window_ns: 5                  # sliding window 大小
    variance_threshold: 0.05      # Å²

  # 接触阈值
  contact_distance: 4.5           # Å

  # MM-GBSA（可选）
  mmpbsa:
    enabled: false
    method: "GBSA"                # PBSA 或 GBSA
    frame_interval_ps: 100
```
