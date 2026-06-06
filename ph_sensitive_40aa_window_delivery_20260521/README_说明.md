# pH-sensitive antibody 40 aa window recommendation 交付说明

生成时间：2026-05-21
项目路径：`/data/ziyang/code/optim-pipe`
结果来源：`results/ph_sensitive_40aa_window/`

## 交付范围

本目录把两个抗体的 40 aa 突变窗口推荐结果、评分数据和说明文件平铺在同一目录下，便于直接查看或转交。

本包包含结果层面的全部 CSV 表、最终推荐、HTML 线性叙事和严格校验报告；不包含 150 个原始 AF3 `model.cif` 结构文件。结构证据已整理为 `af3_model_manifest.csv` 和 `residue_structure_features.csv`。如需回溯原始结构，项目内路径为：

- `results/ph_sensitive_40aa_window/structures/af3/1E62/output/`
- `results/ph_sensitive_40aa_window/structures/af3/sdAb/output/`

## 快速入口

- `result_narrative.html`：面向阅读的线性叙事报告，按 1E62 和 sdAb 分别说明 AF3、结构证据、湿实验证据、容量、糖基化护栏和最终推荐。
- `final_recommendation.md`：审阅后的最终推荐摘要。底层自动评分表仍保留在 `window_scores.csv`，最终结论不再用超过 10k 后的容量差异作为主要排序依据。
- `validation_report.md`：严格校验结果，本轮为 42 PASS / 0 WARN / 0 FAIL。
- `window_scores.csv`：最核心的数据表，包含每个 40 aa 窗口的最终评分、等级、推荐状态和证据摘要。

## 当前推荐

### 1E62

- Primary：`1E62_VL_001_040 (VL:1-40)`
- Key backup：`1E62_VL_042_081 (VL:42-81)`
- Secondary backups：`1E62_VL_045_084 (VL:45-84)`、`1E62_VL_046_085 (VL:46-85)`
- 等级：Grade A
- 证据路线：`antigen_contact_plus_CDR_His`
- 解释：AF3 父复合物模型给出稳定抗原 loop 接触证据，且容量和糖基化风险满足要求。审阅后将 `VL:42-81` 提升为关键备选，因为它覆盖 L42-L44 中性结合保留模块。实际建库需硬保护 L23C。

### sdAb

- Reviewed primary：`sdAb_VHH_072_111 (VHH:72-111)`
- Alternative primary：`sdAb_VHH_075_114 (VHH:75-114)`
- Secondary exploratory：`sdAb_VHH_057_096 (VHH:57-96)`
- 等级：Grade B
- 证据路线：`wet_module_extension`
- 解释：sdAb 已补到 100 个 AF3 结构，AF3 数量不是问题；但没有形成稳定直接抗原 loop / A determinant 接触证据，因此按湿实验模块延展路线推荐。审阅后优先选择覆盖 Q100H、G102H、V105H、E108H、D110H、Y111H 等 CDR3 His 种子的窗口。实际建库需硬保护 A96C，并禁止 V105H 与 D110H 同时出现。

## 关键数据规模

- `candidate_windows.csv`：233 个候选连续 40 aa 窗口；1E62 = 150，sdAb = 83。
- `window_scores.csv`：233 个窗口全部完成评分；无缺失评分窗口。
- `window_mutation_mask.csv`：9,320 行，等于 233 个窗口 x 40 个位置。
- `window_design_capacity_table.csv`：233 行，记录每个窗口的理论设计容量。
- `af3_model_manifest.csv`：150 个 AF3 父复合物模型；1E62 = 50，sdAb = 100。
- `residue_structure_features.csv`：23,600 行逐残基结构特征。
- `glycan_guardrail_table.csv`：9,320 行 N146 糖基化 heuristic guardrail。

## 文件说明

### 推荐与说明

- `result_narrative.html`：推荐结果的 HTML 说明。
- `final_recommendation.md`：最终 primary / backup 推荐。
- `validation_report.md`：严格校验报告。
- `README_说明.md`：本说明文件。

### 窗口与评分

- `candidate_windows.csv`：所有枚举出的连续 40 aa 窗口。
- `window_scores.csv`：每个窗口的最终评分结果。重点字段包括：
  - `grade`：A/B/C/D 等级。
  - `recommendation_status`：`recommended`、`exploratory_only` 或 `blocked`。
  - `evidence_route`：推荐证据路线。
  - `structure_summary`：AF3 结构证据摘要。
  - `wet_prior_summary`：湿实验邻近证据摘要。
  - `glycan_risk_summary`：糖基化护栏摘要。
  - `usable_noncontrol_variant_count`：扣除 control 预算后的理论可用非 control 变体数。
- `window_mutation_mask.csv`：每个窗口每个位点的设计状态、允许突变类别和保护原因。
- `window_design_capacity_table.csv`：每个窗口的理论 library 容量。

### 结构证据

- `af3_model_manifest.csv`：AF3 模型清单、seed/model index、ranking score、模型质量状态。
- `residue_structure_features.csv`：逐残基结构特征，包括抗原距离、抗原 loop 接触、A determinant 接触、N146 距离、结构簇支持和模型质量。

### 湿实验与先验

- `wet_observation_table_1e62.csv`：1E62 纳入评分的湿实验观察整理表。
- `wet_observation_table_sdab.csv`：sdAb 纳入评分的湿实验观察整理表。
- `prior_constraints_table.csv`：人工先验约束，包括 hard protect 和 soft risk。
- `legacy_md_inventory_1e62.csv`、`legacy_md_inventory_sdab.csv`：历史 MD inventory，仅作 context-only，不参与评分。

### 参考映射与抗原区域

- `reference_sequence_map.csv`：抗体和抗原参考序列编号映射。
- `hbsag_aes_regions.csv`：AeS HBsAg 区域标注，用于抗原 loop、A determinant、TM 段等判断。

## 关键字段解释

`usable_noncontrol_variant_count` 不是实际要做的变体数量，而是在当前设计规则下该窗口可供抽样/设计的理论非 control 候选空间：

```text
usable_noncontrol_variant_count = raw_feasible_variant_count - control_variant_budget
```

本轮配置为：

- 最大突变阶数：3，即只计 1、2、3 点突变组合。
- control 预算：256。
- 主库容量最低要求：10,000。

`stable_contact` 表示窗口通过结构主证据判据：多 seed 支持、结构簇支持、接触落在 HBsAg 100-164 / A determinant 等关键区域、主导簇质量可解释，且不是跨膜疏水段或异常构象驱动。

## 注意事项

- 本轮只把湿实验结果和新 AF3 结构结果纳入评分。
- FoldX、pKa、历史 MD 不参与评分；历史 MD 只作为 context-only inventory。
- 糖基化护栏是 N146 几何距离 heuristic，不等同于显式糖链构象采样。
- sdAb 的 B 级不是因为 AF3 数量不足，而是因为 100 个结构下推荐窗口仍没有形成直接稳定抗原界面证据。
- `window_scores.csv` 是自动评分结果；`final_recommendation.md` 是审阅后结论。审阅后将容量作为 10k pass/fail 门槛，而不是超过门槛后的主排序因素。
