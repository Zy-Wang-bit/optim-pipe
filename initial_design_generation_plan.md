# pH 敏感抗体下一阶段初始设计方案

## 0. 结论摘要

本轮设计应把 **1E62** 和 **sdAb** 都作为结果导向的 pH-sensitive 主库推进，而不是把 sdAb 作为机制探索库。两者最终目标一致：获得在 **pH 7.4 保留结合、pH 6.0 增强释放** 的候选。区别在于 sdAb 先验不确定性更高，因此需要更大的计算候选池、更严格的中性结合 gate，以及更强的 rescue 约束。

推荐设计规模：

| 阶段 | 1E62 | sdAb | 说明 |
|---|---:|---:|---|
| raw generated candidates | 400k–500k | 500k–600k | 生成到 post-hard-filter quota 满足为止，不固定 raw 数 |
| effective initial pool after hard filters / dedup | 250k | 300k | 本轮真正的“初始设计池” |
| Tier 1 retained | ~50k | ~60k | 用轻量工具和 surrogate 过滤 |
| Tier 2-light reviewed | ~20k | ~24k | 结构/pH 机制复核的主集合 |
| Tier 2-heavy representatives | 1k–2k | 1.5k–2.5k | 只做代表性重型复核 |
| final wet-lab library | 10k | 10k | 每个 target 的最终建库规模 |

核心比例：

```text
1E62: final 10k = performance-first, high P_hit 主导
sdAb: final 10k = performance-first, stricter neutral-retention gate
```

如果实际合成通量是“每个 target 10k 总量且 controls 也计入”，则建议：

```text
9,800 main design variants + 200 anchors/controls
```

如果 controls 可以单独做，则建议保留完整：

```text
10,000 non-control variants + independent control panel
```

---

## 1. 初始设计应该如何产生

初始设计不应从 40 aa 窗口内随机组合产生，也不应简单枚举所有可能突变。推荐采用 **seeded, constrained, result-oriented generation**：先定义可突变位点和允许氨基酸集合，再按路线生成候选，生成过程中实时执行 hard filter、去冗余和 route quota 控制。

总体流程：

```text
window + mask freeze
→ route-specific candidate generation
→ hard filter at generation time
→ dedup / cluster
→ cheap feature annotation
→ surrogate P_hit scoring
→ effective initial pool
```

### 1.1 先冻结设计空间

两个窗口固定为：

```text
1E62: VL 1–40
sdAb: VHH 72–111
```

每个位点必须先分配状态：

```text
hard_protect          绝不突变
restricted            只允许少数保守或 rescue 替换
designable_main       主设计位点
designable_rescue     只作为 rescue 位点开放
diagnostic_only       原则上不进入本轮主库，除非同时具备 P_hit 潜力
forbidden             禁止突变或禁止特定组合
```

本轮不做大规模机制网格，因此 `diagnostic_only` 位点不应消耗大量容量。

### 1.2 突变阶数和 His 数量

推荐约束：

| target | 主体突变阶数 | 上限 | His 上限 | 说明 |
|---|---:|---:|---:|---|
| 1E62 | 2–4 | 4 | 通常 ≤2，少数 ≤3 | 成熟窗口，可允许更多 His + rescue 组合 |
| sdAb | 2–3 | 3，少数强支持可到 4 | 通常 ≤2 | CDR3/FR3 风险更高，必须严控 global weakening |

sdAb 的设计原则不是保守，而是 **强约束地追求命中**。因此不应铺开 CDR3 His cluster，也不应系统性扫描 FR3；FR3 只在能支持 neutral retention 或 rescue 时进入。

---

## 2. 初始候选生成路线

### Route A：His-trigger rule design

目的：直接产生 pH-sensitive trigger。

生成内容：

```text
single-His
selected double-His
limited local His pair
His + conservative mutation
His + charge/polar rescue
known His seed expansion
```

不生成：

```text
unbounded His cluster
known forbidden pair
high-order His-only combinations
neutral-retention unsupported CDR3 His expansion
```

His seed 位点来源：

```text
wet-lab prior positive or neutral-retained signal
interface-near or conformation-regulating position
pH electrostatic environment plausible
not hard-protected
not glycan-forbidden
not known global-weakening driver
```

预期 seed 数：

```text
1E62: 8–12 个 His seed 位点
sdAb: 6–10 个 His seed 位点
```

如果 sdAb 可用 seed 少于 6 个，不建议强行扩展到更多位点；应增加 rescue 和 MPNN-seeded 路线比例。

---

### Route B：ProteinMPNN-seeded rescue

目的：不是单纯生成多样序列，而是在固定 His trigger 的情况下寻找能保留 pH 7.4 结合的周边 rescue mutation。

推荐运行方式：

```text
1. 固定一个或两个 His seed。
2. 锁定 hard_protect 位点。
3. 只开放窗口内 designable_rescue / designable_main 位点。
4. 禁止 Cys、N-X-S/T、明显 liability residues。
5. 使用 restricted amino-acid alphabet，而不是全 20 aa。
6. 多 temperature / 多 seed 采样后去冗余。
```

1E62 中 ProteinMPNN-seeded rescue 是主力路线之一。sdAb 中也应使用，但要更强地限制 CDR3 中高阶组合。

---

### Route C：wetlab-informed module expansion

目的：利用已有实验数据，把已有正/负信号转成候选生成规则。

生成内容：

```text
effective His seed + low-risk rescue
effective module + additional pH trigger
neutral-retained mutation + His trigger
failed combination minus risky mutation + rescue
known weak-binding motif exclusion
```

这个路线应优先于纯随机组合，因为它最贴近当前 assay 和展示体系。

---

### Route D：structure/interface/glycan-constrained design

目的：提高候选与真实结合界面、抗原覆盖和糖链可及性的兼容性。

生成规则：

```text
interface-near positions can be prioritized
hotspot-disrupting substitutions are removed
glycan-near bulky/aromatic/strong positive substitutions are restricted
framework positions are used mainly as rescue, not as broad exploration
single structural model cannot be the sole evidence for inclusion
```

该路线主要用于候选生成和风险标记，不单独决定最终入库。

---

### Route E：constrained high-upside candidates

目的：保留少量计算信号强但有单项风险的候选。

约束：

```text
must pass hard filter
must pass minimal neutral-retention gate
must not violate forbidden pair
must have clear pH-trigger rationale
must be sequence-diverse from main candidates
```

该路线数量应小，不能演变成机制探索库。

---

## 3. 推荐初始设计数量

这里的“初始设计数量”定义为：**通过 hard filters 和去冗余后的 effective initial pool**。raw 生成数不固定，应该生成到每条路线的 effective quota 满足为止。

### 3.1 1E62 初始设计池：250k effective candidates

| route | effective candidates | 占比 | 目的 |
|---|---:|---:|---|
| His-trigger rule design | 55k | 22% | 直接构建 pH trigger |
| ProteinMPNN-seeded His rescue | 90k | 36% | 提高 pH 7.4 保留概率 |
| wetlab-informed module expansion | 50k | 20% | 利用已有正/负信号 |
| structure/glycan-constrained design | 30k | 12% | 降低界面和糖链风险 |
| low-His / non-His rescue diversity | 20k | 8% | 增加高分候选多样性 |
| constrained high-upside | 5k | 2% | 保留少数强信号风险候选 |
| **total** | **250k** | **100%** |  |

1E62 的初始池以 **His + rescue** 为主，因为该窗口更适合直接追求性能命中。

### 3.2 sdAb 初始设计池：300k effective candidates

| route | effective candidates | 占比 | 目的 |
|---|---:|---:|---|
| constrained His-trigger rule design | 65k | 21.7% | 直接构建 pH trigger，但限制 His cluster |
| His + conservative / MPNN rescue | 125k | 41.7% | 抵消 CDR3/FR3 中性结合风险 |
| wetlab-informed module expansion | 45k | 15% | 优先使用已有实验信号 |
| neutral-retention / structure-stable variants | 40k | 13.3% | 把 sdAb 的 global weakening 风险压低 |
| coverage / glycan-safe high-score variants | 15k | 5% | 降低覆盖和糖链遮挡风险 |
| constrained high-upside | 10k | 3.3% | 保留少量高潜力但有单项风险的候选 |
| **total** | **300k** | **100%** |  |

sdAb 初始池比 1E62 多 20%，但最终湿实验容量仍然是 10k。这样做的原因是 sdAb 先验不确定性更高，应该用更多计算候选换取更严格的最终筛选，而不是消耗湿实验容量做机制探索。

---

## 4. 初始候选到最终 10k 的压缩策略

### 4.1 Tier 1：轻量全量筛选

对 effective initial pool 全量运行：

```text
hard filter validation
ESM / antibody language-model acceptability
FoldX or equivalent fast pH7.4 / pH6.0 estimate
interface / hotspot scan
liability / developability scan
data-driven P_hit surrogate
basic diversity clustering
```

Tier 1 只判断是否值得继续，不作为最终活性判断。

目标保留：

```text
1E62: 250k → ~50k
sdAb: 300k → ~60k
```

### 4.2 Tier 2-light：选择性结构和 pH 机制复核

对 Tier 1 通过候选运行轻中量复核：

```text
PyRosetta local repack / minimize
PROPKA3 + pKAI+ for His candidates
Rosetta pH-score / dddG_elec for selected top candidates
coarse CDR RMSD / geometry risk
coarse glycan proximity risk
```

目标复核规模：

```text
1E62: ~20k
sdAb: ~24k
```

### 4.3 Tier 2-heavy：只做代表性重型复核

不对所有候选运行重型工具。

对象：

```text
top P_hit clusters
sdAb CDR3-risk high-score candidates
glycan-risk but otherwise strong candidates
final library cluster representatives
```

工具可选：

```text
AF3 / multi-model complex check
explicit glycan modeling
APBS-like electrostatics
SimpleFold 3x
short MD only for final top representatives or post-hit validation
```

目标规模：

```text
1E62: 1k–2k representatives
sdAb: 1.5k–2.5k representatives
```

### 4.4 Tier 3：结果导向最终选择

最终不是选 top score，而是约束优化：

```text
maximize:
  P_hit
  P_neutral_retained
  P_acidic_release

minimize:
  P_global_weakening
  expression/display risk
  glycan/coverage risk
  sequence redundancy

constraints:
  hard filters pass
  forbidden pairs excluded
  His count within limit
  mutation count within limit
  cluster overrepresentation avoided
```

---

## 5. 最终 10k 库组成

### 5.1 1E62 final 10k

| class | count | 说明 |
|---|---:|---|
| high-confidence pH-switch candidates | 7,000 | P_hit 最高，neutral retention 和 acidic release 均强 |
| His + rescue candidates | 1,300 | His trigger 与非 His rescue 组合 |
| coverage / glycan-safe candidates | 600 | 降低抗原覆盖和糖链遮挡风险 |
| high-score diverse candidates | 500 | 在高分候选内保证序列/模块多样性 |
| constrained high-upside candidates | 300 | 单项风险存在但 pH-sensitive 信号强 |
| anchors / controls | 200 | 亲本、已知好/坏模块、assay 校准 |
| **total** | **10,000** |  |

如果 controls 可单独做，则把 `anchors / controls` 的 200 个名额并入 high-confidence pH-switch candidates。

### 5.2 sdAb final 10k

| class | count | 说明 |
|---|---:|---|
| high-confidence pH-switch candidates | 6,200 | 以结果为目标，不做大规模机制网格 |
| His + conservative / polar rescue | 1,800 | 防止 His 造成 pH 7.4 结合崩塌 |
| neutral-retention-prioritized candidates | 900 | sdAb 中性结合风险较高，必须保留该类 |
| coverage / glycan-safe candidates | 400 | 降低 coverage loss 或糖链遮挡误判 |
| high-score diverse candidates | 400 | 只在高 P_hit 候选中增加多样性 |
| constrained high-upside candidates | 200 | 有潜在 pH switch，但必须通过 hard gate |
| anchors / controls | 100 | 最小必要校准对照 |
| **total** | **10,000** |  |

sdAb 的 final 10k 中不再保留大比例 diagnostic grid。所有入库候选都必须有成为 pH-sensitive hit 的合理可能性。

---

## 6. 推荐新增计算模块

不是所有工具都必须用。本轮应采用 ROI-driven tool selection。

### 6.1 必做或强烈推荐

| 模块 | 用途 | 运行范围 |
|---|---|---|
| hard filter + liability scan | 去掉 Cys、N-X-S/T、禁配、明显 developability 风险 | all candidates |
| ESM / antibody LM | 过滤序列异常性 | all candidates |
| FoldX or equivalent fast energy | 快速估计中性结合风险和 pH delta | all candidates |
| interface / hotspot scan | 避免破坏关键界面 | all candidates |
| data-driven P_hit surrogate | 用既有湿实验校准排序 | all candidates |
| constrained optimizer | 从高分候选中构建最终 10k | Tier 3 |

### 6.2 选择性使用

| 模块 | 用途 | 运行范围 |
|---|---|---|
| PyRosetta local repack/minimize | 局部结构可行性 | Tier 1 后子集 |
| PROPKA3 + pKAI+ | His pKa 支持 | His candidates / top subset |
| Rosetta pH-score / dddG_elec | pH 依赖静电方向 | top subset |
| SimpleFold 3x | CDR 偏移和构象一致性 | sdAb 高风险候选 / top representatives |
| AF3 multi-model complex | 检查界面稳定性和假界面 | top representatives |
| explicit glycan modeling | N146 等糖链遮挡风险 | glycan-risk candidates / top representatives |
| MD / GROMACS | 最终少数候选验证 | post-hit 或最终 top candidates |

### 6.3 不建议全量使用

```text
full AF3 for all candidates
full MD or constant-pH MD for all candidates
large explicit glycan ensemble for all candidates
large diagnostic mechanism grid
single-score top 10k selection
```

---

## 7. 数据驱动 P_hit surrogate

建议新增一个轻量数据驱动模型，用既有实验结果校准计算指标。

目标输出：

```text
P_hit
P_neutral_retained
P_acidic_release
P_global_weakening
P_display_ok
P_coverage_ok
```

输入特征：

```text
mutation identity
mutation count
His count
CDR/FR location
ESM features
FoldX features
interface distance / hotspot features
pKa features when available
glycan proximity features
prior wetlab module tags
```

推荐模型：

```text
calibrated logistic regression
ridge / elastic net
gradient boosted trees
small multitask model
```

模型目的不是替代湿实验，而是把已有 assay 的经验转成候选排序信号。

---

## 8. 最终入库判定规则

每条 final 10k 候选必须满足：

```text
hard filters pass
no forbidden pair
neutral retention risk acceptable
acidic release rationale plausible
global weakening risk not high
expression/display liability acceptable
glycan/coverage risk acceptable or explicitly justified
cluster not overrepresented
```

每条入库记录必须包含：

```text
target
variant_id
sequence
mutation_list
mutation_count
His_count
source_route
candidate_class
P_hit
P_neutral_retained
P_acidic_release
P_global_weakening
P_display_ok
P_coverage_ok
structure_risk
glycan_risk
selection_reason
```

`selection_reason` 应该是结果导向的，例如：

```text
high P_hit and low neutral-risk
His-trigger with rescue support
strong pH-electrostatic support
glycan-safe high-score candidate
diverse representative from high-score cluster
```

不应使用以下理由：

```text
included for broad mechanism exploration
included to complete diagnostic grid
included because structurally interesting only
```

---

## 9. 推荐交付物

本阶段建议输出：

```text
design_config.yaml
evidence_ledger.csv
initial_pool_1E62_250k.csv
initial_pool_sdAb_300k.csv
tier1_results_1E62.csv
tier1_results_sdAb.csv
tier2_light_results_1E62.csv
tier2_light_results_sdAb.csv
tier2_heavy_representatives.csv
final_library_1E62_10k.csv
final_library_sdAb_10k.csv
control_panel.csv
selection_audit.md
wetlab_interpretation_plan.md
```

最重要的是 `selection_audit.md`：它需要说明每个 route 的生成逻辑、候选数量、过滤损失、最终入库比例，以及为什么没有把容量用于大规模机制探索。

---

## 10. 最终建议

本轮初始设计应采用：

```text
1E62: 250k effective initial candidates → 10k final
sdAb: 300k effective initial candidates → 10k final
```

raw 生成数不作为目标，应该生成到 post-hard-filter quota 满足为止。sdAb 计算候选池更大，但最终湿实验容量不扩大；它不是机制探索库，而是用更严格 gate 从更大的计算池里挑出最有机会实现 pH 敏感性的候选。

一句话总结：

> 初始设计应由 His-trigger、ProteinMPNN-seeded rescue、wetlab-informed expansion 和 structure/glycan-constrained design 四条结果导向路线产生；1E62 生成 250k effective candidates，sdAb 生成 300k effective candidates，经过轻量全量筛选、选择性结构/pH 复核和约束优化后，各自收敛到 10k final library。
