# 实验室笔记手写简版

用途：把下面内容直接抄到实验本。图建议贴 `figures/` 里的 PNG。

## 1. HBsAg 型别泛化机制

**目的**：解释 1E62 R2 变体在 Ae/B/Ce/D1 上结合差异大的原因。
**方法**：整理 R2 pH7.4 ELISA，比较 HBsAg 差异位点，并用 AF3 模板做距离扫描。
**结果**：L-only 的 com8/com9/com18 三型都强；H-only 的 com16/com20 在 D1 低。
**结论**：关键差异集中在 117-147，尤其 122/126 附近；D1 不能代表 Ce。

| 变体 | 类型 | Ae | B | D1 | 结论 |
| --- | --- | ---: | ---: | ---: | --- |
| com8 | L-only | 3.62 | 3.86 | 3.66 | 广谱强 |
| com9 | L-only | 3.94 | 3.00 | 2.94 | 广谱强 |
| com18 | L-only | 1.38 | 1.60 | 4.01 | 广谱强 |
| com16 | H-only | 1.95 | 0.78 | 0.01 | D1 失败 |
| com20 | H-only | 2.13 | 1.44 | 0.06 | D1 失败 |

配图：`fig1_r2_elisa_key_variants.png`、`fig2_hbsag_type_sensitive_map.png`

## 2. FoldX 多型别验证

**目的**：检查 FoldX 能否预测不同 HBsAg 型别上的结合强弱。
**方法**：用 Ae/B/D1 模板做 BuildModel Dif energy，并和 ELISA 表型比较。
**结果**：Ae/B/D1 AUC = 0.552/0.521/0.594。
**结论**：FoldX 计算完成，但方法没有验证成功，不能作为主筛选器。

| 型别 | AUC |
| --- | ---: |
| Ae | 0.552 |
| B | 0.521 |
| D1 | 0.594 |

配图：`fig3_foldx_auc_boundary.png`

## 3. CeS AF3 模板

**目的**：补齐 CeS 结构模板。旧 CeS 结果 light chain 只有 111 aa，不能用。
**方法**：用正确 WT heavy/light/CeS FASTA 重跑 AF3，随机 5 个 seeds。
**结果**：5 个 rank 全部通过序列核对，H=115 aa，L=113 aa，CeS=226 aa。
**结论**：corrected-light CeS WT 模板可用于后续接触扫描。

| Rank | Seed | Score |
| ---: | ---: | ---: |
| 1 | 495296019 | 0.696 |
| 2 | 635954951 | 0.512 |
| 3 | 1129861070 | 0.509 |
| 4 | 75801088 | 0.508 |
| 5 | 2026110668 | 0.469 |

配图：`fig4_af3_ces_rank_scores.png`

## 4. 1E62 R3 设计

**目的**：在 com18 L 链 His 簇背景下测试 H 链突变组合。
**方法**：设计 10 个 R3 变体，pH7.4 和 pH6.0 各做 50 ns MD。
**结果**：42E、cross1、com16diff 的酸性接触减少最明显。
**结论**：R3-T1-42E 和 R3-T1-cross1 最值得关注；com17 全套组合不稳。

| 变体 | RMSD@7.4 | ΔContacts | 结论 |
| --- | ---: | ---: | --- |
| 42E | 1.9 | -165 | 重点 |
| cross1 | 1.6 | -131 | 重点 |
| com16diff | 3.9 | -154 | 响应强但偏弱 |
| com17diff | 23.3 | +45 | 不稳 |
| com17clean | 8.3 | +27 | 不稳 |

配图：`fig5_1e62_r3_md_delta_contacts.png`

## 5. sdab 32 变体

**目的**：交付一批 pH 响应 sdab 候选。
**方法**：三条路线并行：MPNN 全设计、ELISA 模型筛选、FR-His 新机制。
**结果**：共 32 个变体；ELISA 模型和 FR-His 方向较好。
**结论**：r4_04/r4_01、A4/A1、v04/v14 是重点观察候选。

| 路线 | 数量 | 方向较好 |
| --- | ---: | ---: |
| MPNN 全设计/精简/Y37 | 12 | 3 |
| ELISA 模型 | 10 | 9 |
| FR-His | 10 | 7 |

配图：`fig6_sdab_32_variant_overview.png`

## 6. DeepScientist blind-start

**目的**：准备只含湿实验数据的启动包。
**边界**：给湿实验表、序列、目标；暂不给 FoldX/MD/模型分数和内部结论。
**状态**：任务记录显示已整理过包，但还没有启动 DeepScientist quest。正式启动前需要复核包路径。
