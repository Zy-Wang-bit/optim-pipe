# Phase 1 · 突变位点级组合设计 + 计算评估

> 基于 R2 湿实验结果，以 com18 为骨架，从 com16/17 H 链中筛选安全位点进行交叉组合。
> 在提交湿实验前，通过多方法计算评估提供置信度排序和风险预警。
> 对应 ROADMAP § 一

---

## 1.1 问题定义

R2 的 20 个变体采用单侧设计（每个变体仅 H 链或 L 链一侧突变），实验结果表明：

- **广谱结合 + pH 敏感性**：仅 com18 同时满足（B 型 ratio 0.55，三基因型广谱，表达量 11.8 µg/µl）
- **com8/9**：广谱结合但 pH 敏感性不如 com18 显著，与 com18 差异位点有限，价值不大
- **H 链突变组**（com4/12/16/17/20）：
  - com12/20 结合活性过弱，排除
  - com4/16/17 部分基因型有结合但无目标方向的 pH 敏感性
  - 所有 H 链变体都携带 H35-37 HHH 簇（CDR1 核心），与 D1 丧失强相关

**目标**：以 com18 为唯一 L 链骨架，从 com16/17 的 H 链中拆出安全位点，通过交叉组合探索基因型选择性的驱动因素。

## 1.2 SAR 分析

### 1.2.1 排除 H35-37 HHH 簇

| 证据 | 结论 |
|------|------|
| com16/17 均携带 H35H+H36H+H37H，且 D1 结合丧失或极弱 | HHH 簇与 D1 丧失强相关 |
| H36 位为高度保守的色氨酸（W），W→H 替换破坏性大 | 直接影响 CDR1 结构 |
| com16 出现基因型依赖的 pH 反转（Ae 酸性增强 vs B 酸性减弱） | pH 效应方向不可控 |

**决定：H35/H36/H37 全部排除，不纳入任何 R3 变体。**

### 1.2.2 非 HHH 位点分层

| 位点 | WT | com16 | com17 | 分层 | SAR 依据 |
|------|-----|-------|-------|------|---------|
| H18 | L | I | I | **Tier 1** | 双变体共享同一残基 → 高置信安全 |
| H19 | R | T | T | **Tier 1** | 同上 |
| H46 | E | Q | Q | **Tier 1** | 同上 |
| H48 | V | I | I | **Tier 1** | 同上 |
| H42 | G | E | A | **差异位点** | 性质差异最大（负电荷 vs 疏水），最可能驱动基因型选择性 |
| H23 | A | Q | K | **差异位点** | 电荷差异（中性极性 vs 正电荷），可能影响界面静电环境 |
| H43 | K | S | N | **差异位点** | 侧链大小差异，辅助效应 |
| H21 | S | — | T | **com17 独有** | 保守突变（S→T），证据单一 |
| H40 | S | — | N | **com17 独有** | 证据单一 |
| H45 | L | — | M | **com17 独有** | 保守突变（L→M），证据单一 |

### 1.2.3 差异位点与基因型选择性

com16 和 com17 在基因型上的 pH 表现截然不同：

| 基因型 | com16 ratio | com17 ratio | 差异 |
|--------|------------|------------|------|
| Ae | **0.54**（酸性解离） | 0.84（不敏感） | com16 敏感，com17 不敏感 |
| B | **1.31**（酸性增强） | **0.71**（酸性解离） | 方向相反 |
| D1 | 无结合 | 极弱 | 均丧失（HHH 驱动） |

差异位点（H42/H23/H43）是两者在共享 Tier 1 和 HHH 的基础上唯一不同的位点，因此是表型差异的最可能驱动因素：

- **H42E**（com16，负电荷）可能驱动 Ae 型敏感性
- **H42A**（com17，疏水）+ **H23K**（正电荷）可能驱动 B 型敏感性
- **核心假设**：交叉组合（如 H42E + H23K）有可能同时获得 Ae + B 的敏感性

> 注意：上述 pH 表型是在 HHH 簇存在的背景下观测到的。排除 HHH 后差异位点单独能否影响 pH 行为不确定，但仍能回答框架突变对广谱结合的影响。

## 1.3 组合设计方案

**固定 L 链骨架**：com18 的 L 链突变（L24Q, L42H, L43H, L44H, L51Q, L62P）

| # | 变体名 | H 链突变（Tier 1 以外的） | 总突变数 | 设计逻辑 |
|---|--------|--------------------------|---------|---------|
| 1 | R3-base | 无（WT H 链） | 0+6 | 基线对照：com18 L 链 + WT H 链 |
| 2 | R3-T1 | — | 4+6 | Tier 1 基线：4 个高置信共享位点，检验对 com18 骨架是否中性 |
| | | **—— 单位点拆分 ——** | | |
| 3 | R3-T1-42E | + H42E | 5+6 | 隔离 com16 的 H42（负电荷）→ 检验是否驱动 Ae 敏感性 |
| 4 | R3-T1-42A | + H42A | 5+6 | 隔离 com17 的 H42（疏水）→ 检验是否驱动 B 敏感性 |
| | | **—— 差异位点全集 ——** | | |
| 5 | R3-T1-com16diff | + H23Q+H42E+H43S | 7+6 | com16 全套差异位点，对比 #3 看 H23Q 和 H43S 的增量效应 |
| 6 | R3-T1-com17diff | + H23K+H42A+H43N | 7+6 | com17 全套差异位点，对比 #4 看 H23K 和 H43N 的增量效应 |
| | | **—— 交叉组合 ——** | | |
| 7 | R3-T1-cross1 | + H42E+H23K | 6+6 | 核心交叉：com16 的 H42E + com17 的 H23K，尝试融合 Ae+B 敏感性 |
| 8 | R3-T1-cross2 | + H42E+H23K+H43N | 7+6 | 在 #7 基础上加 com17 的 H43N，测试 H43 额外贡献 |
| | | **—— 镜像交叉 + 完整模块对照 ——** | | |
| 9 | R3-T1-cross3 | + H42A+H23Q | 6+6 | 镜像交叉：com17 的 H42A + com16 的 H23Q，cross1 的对照组 |
| 10 | R3-com17clean | T1 + H21T+H23K+H40N+H42A+H43N+H45M | 10+6 | com17 全套非 HHH（10 位点），完整 com17 框架对照 |

**实验解读逻辑：**
```
#1 vs #2  → Tier 1 共享突变对 com18 骨架的影响（中性？增强？）
#3 vs #4  → H42 单点效应：负电荷(E) vs 疏水(A) 对基因型选择性的影响
#5 vs #3  → H23Q+H43S 在 H42E 背景下的增量效应
#6 vs #4  → H23K+H43N 在 H42A 背景下的增量效应
#7        → 关键交叉：H42E+H23K 能否融合 Ae+B 的 pH 敏感性？
#8 vs #7  → H43N 的额外贡献
#9 vs #7  → 镜像交叉对照：H42A+H23Q vs H42E+H23K，验证 H42/H23 的方向性
#10 vs #6 → com17 完整框架 vs 仅差异位点（含 H21T/H40N/H45M 独有位点）
```

## 1.4 计算评估方案

### 1.4.1 目标

当前 10 个 R3 变体的设计纯基于 SAR 规则推理（消去法），缺乏结构/计算层面的正面证据。在提交湿实验前，通过多方法计算评估：

1. 为每个变体附上计算置信度，辅助优先级排序
2. 排除能量严重恶化或结构不稳定的组合
3. 评估 H 链框架突变对 L 链 His（pH 开关核心）pKa 的影响

### 1.4.2 评估 Pipeline

6 步串行流程，一键运行：`bash experiments/1E62_R3/eval_r3.sh`

| Step | 方法 | 工具 | 输出 | 环境 |
|------|------|------|------|------|
| 1 | 突变体结构生成 | PyRosetta (MutateResidue + Repack + Minimize) | 10 个突变体 PDB | pyrosetta |
| 1b | 辅助结构验证 | SimpleFold (序列→结构从头预测) | 10×3 ensemble PDB | simplefold |
| 2 | pKa 预测 | PROPKA3 + pKAI+ 双工具交叉验证 | His pKa shift + consensus | optim-pipe |
| 3 | 能量评估 | FoldX AnalyseComplex (pH 7.4 + 6.0) | ddG, dddG | optim-pipe |
| 4 | pH 敏感性评分 | Rosetta pH-score + dddG_elec | His 氢键网络 + 静电能差 | pyrosetta |
| 5 | 结构稳定性 | CDR RMSD (全局 + 6 个 CDR 区域) | RMSD per CDR | optim-pipe |
| 6 | 综合排名 | 加权评分 + 筛选 + 报告 | CSV + Markdown 报告 | optim-pipe |

### 1.4.3 结构生成策略

**主要方法：PyRosetta**（MutateResidue + 8Å shell Repack + sidechain Minimize）
- 保留完整复合物（抗体+抗原），可直接用于后续 FoldX/Rosetta 评分
- 模板：`data/pdbs/n1-0_Repair.pdb`（FoldX 修复后的 AF3 结构）

**辅助验证：SimpleFold**（3B 参数，flow-matching 模型）
- 从序列从头预测，无模板偏差
- 每个变体生成 3 个 ensemble 构象
- 仅用于 RMSD 对比，不含抗原

### 1.4.4 pKa 预测方案

**核心问题**：R3 的 H 链框架突变是否改变了 L 链 His（B42/B43/B44）的 pKa？

**工具选择**：PROPKA3（经验公式）+ pKAI+（深度学习）双工具交叉验证
- 两者精度相当（MAE ~0.6-0.7），方法论互补
- 共识策略：双工具 shift 方向一致 → 高置信；矛盾 → 标记为"不确定"
- His model pKa = 6.5；理想 pH 开关窗口 6.0-7.4

**可复用 wrapper**：`analysis/pka/run_pka.py`（跨项目共享）

### 1.4.5 综合评分

**筛选规则**（基于 R2 经验）：

| 规则 | 阈值 | 依据 |
|------|------|------|
| ddG_pH7.4 | < 0.3 kcal/mol | R2: ddG<0.3 → ~78% hit rate |
| 全局 RMSD | < 1.5 Å | ROADMAP § 二假设（待 R2 验证） |
| pKa shift | \|shift\| > 1.0 触发预警 | 探索性阈值 |

**加权评分**：

| 指标 | 权重 | 理由 |
|------|------|------|
| ddG_pH7.4 | 0.50 | R2 最强预测因子 (Spearman r = -0.59 ~ -0.78) |
| dddG | 0.20 | pH 敏感性直接指标 |
| ph_score | 0.15 | His 氢键网络质量 |
| pka_avg_shift | 0.10 | 新指标，探索性 |
| dddG_elec | 0.05 | R2 相关性弱 |

**输出**：
- `r3_evaluation_report.csv` — 全量数据 + 排名 + RED/YELLOW/GREEN 推荐
- `r3_evaluation_summary.md` — 可读报告（含 caveats）

### 1.4.6 Caveats

- R2 计算预测 false positive rate 为 55%（11/20），综合评分仅供参考，不替代湿实验
- pKa 工具首次整合，为探索性指标；R3 实验结果回收后将评估其预测精度
- 结构来源为 PyRosetta 局部优化（非从头预测），远程效应可能捕捉不足

## 1.5 文件说明

### 设计文件（`experiments/1E62_R3/design/`）

| 文件 | 说明 |
|------|------|
| `r2_labels.csv` | 20 行标签文件（variant_id, h_chain_seq, l_chain_seq, active） |
| `combine_variants.py` | 输入 WT 序列 + 位点定义 → 输出 `r3_combinations.csv` |
| `r3_combinations.csv` | 10 个变体的完整定义（突变、序列、设计逻辑） |

### 评估文件（`experiments/1E62_R3/`）

| 文件 | 说明 |
|------|------|
| `config_r3_eval.yaml` | 评估配置（路径、链、CDR 定义、阈值、权重） |
| `eval_r3.sh` | 一键运行入口（自动切换 conda 环境） |
| `scripts/build_r3_structures.py` | Step 1: 多方法结构生成 (PyRosetta / SimpleFold) |
| `scripts/run_pka_batch.py` | Step 2: pKa 批量预测 |
| `scripts/run_foldx_r3.py` | Step 3: FoldX 能量评估 |
| `scripts/run_rosetta_r3.py` | Step 4: Rosetta pH-score + dddG_elec |
| `scripts/run_rmsd_r3.py` | Step 5: CDR RMSD |
| `scripts/rank_r3_variants.py` | Step 6: 综合排名 + 报告 |

### 可复用工具（`analysis/`）

| 文件 | 说明 |
|------|------|
| `analysis/pka/run_pka.py` | pKa 预测 wrapper (PROPKA3 + pKAI+)，跨项目共享 |

### 结果目录（`experiments/1E62_R3/results/`）

```
results/
├── structures/{rosetta,simplefold}/   # 突变体 PDB
├── pka/                               # pKa 预测结果
├── foldx/                             # FoldX 能量结果
├── rosetta/                           # Rosetta 评分结果
├── rmsd/                              # RMSD 结果
├── r3_evaluation_report.csv           # 综合评估表
└── r3_evaluation_summary.md           # 可读报告
```

## 1.6 验证

### 设计验证

```bash
python experiments/1E62_R3/design/combine_variants.py
```

检查项：
- 输出 10 行
- 所有 L 链与 com18 L 链一致
- H 链突变不含 H35H/H36H/H37H
- 序列长度与 WT 一致（H: 112, L: 110）

### 计算评估

```bash
bash experiments/1E62_R3/eval_r3.sh
```

检查项：
- `results/structures/rosetta/` 下 10 个 PDB，突变位点残基正确
- `results/pka/pka_summary.csv` 中 WT His pKa 在 5.5-7.5 合理范围
- `results/foldx/foldx_summary.csv` 中 R3-base 的 ddG 趋势与 R2 com18 一致
- `results/rosetta/ph_scores.csv` 中 His 数量 = 3（B42/B43/B44）
- `results/rmsd/rmsd_summary.csv` 中 R3-base RMSD 接近 0
- `results/r3_evaluation_report.csv` 中 10 行无 NaN，recommendation 列有 GREEN/YELLOW/RED
