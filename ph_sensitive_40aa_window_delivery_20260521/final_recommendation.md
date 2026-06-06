# pH 敏感抗体 40 aa 窗口审阅后最终推荐

本文件是审阅后的最终结论。底层自动逐窗口评分仍保留在 `tables/window_scores.csv`；本结论使用自动评分作为证据基础，但在容量超过 10k 后，不再把容量差异作为主要排序依据，而优先考虑湿实验信息覆盖、机制可解释性、结构证据和风险控制。

糖链说明：当前糖基化风险来自 N146 几何距离 heuristic guardrail，不是显式糖链构象采样。

## 数据基础

- 已评分候选窗口：233 个；1E62 = 150，sdAb = 83。
- AF3 亲本复合物模型：150 个；1E62 = 50，sdAb = 100。
- 严格校验：42 PASS / 0 WARN / 0 FAIL。
- 容量规则：`usable_noncontrol_variant_count >= 10,000` 作为是否满足 1w 通量的门槛；超过门槛后，容量不作为主排序因素。

## 1E62

### 首选窗口

`1E62_VL_001_040 (VL:1-40)`

- 审阅后角色：主设计窗口
- 自动评分：`A` / `recommended`
- 证据路线：`antigen_contact_plus_CDR_His`
- pH 机制假设：strong
- 容量：raw `23,380,434`，usable non-control `23,380,178`
- 结构证据：`stable_contact; seed_support=10; cluster_support=25; ag_loop_fraction=0.999; quality_ok=True`
- 糖链风险：low heuristic risk

推荐理由：

- 覆盖轻链 CDR1 及邻近框架区。
- 包含 L26H 及周边多个弱 pH 敏候选位点。
- AF3 模型稳定支持与 HBsAg 抗原性 loop 的接触。
- 现有湿实验信息支持中性结合保留和酸性释放机会的平衡。
- 容量远高于 10k 门槛。

实际建库约束：

- 必须硬保护 L23C 及其他必需框架/核心位点。
- 控制 His 总数，优先单 His 或双 His 加中性结合补偿组合。
- 避免已知或疑似负协同 His cluster。

### 备选窗口

1. `1E62_VL_042_081 (VL:42-81)`
   - 审阅后角色：关键备选 / 中性结合保留模块
   - 自动评分：`A` / `recommended`
   - 证据路线：`antigen_contact_plus_CDR_His`
   - 容量：usable non-control `23,248,941`
   - 理由：覆盖 L42-L44 相关中性结合保留模块。自动排序中它因容量略低排在 VL:45-84 和 VL:46-85 后面；审阅后认为容量已充分超过 10k，因此将其提升为更有机制意义的备选。

2. `1E62_VL_045_084 (VL:45-84)`
   - 审阅后角色：次级备选
   - 自动评分：`A` / `recommended`
   - 容量：usable non-control `23,380,178`
   - 理由：结构和容量都强，可作为 VL:42-81 之后的备选。

3. `1E62_VL_046_085 (VL:46-85)`
   - 审阅后角色：次级备选
   - 自动评分：`A` / `recommended`
   - 容量：usable non-control `23,380,178`
   - 理由：与 VL:45-84 接近；因更少覆盖 L42-L44 中性保留模块，优先级低于 VL:42-81。

### 1E62 结论

1E62 可以直接推进 10k 库设计。首选 `VL:1-40`；最有意义的备选是 `VL:42-81`；`VL:45-84` 和 `VL:46-85` 作为次级备选。正式设计前需要把 L23C 等必需位点纳入 hard protect。

## sdAb

### 审阅后首选窗口

`sdAb_VHH_072_111 (VHH:72-111)`

- 审阅后角色：主探索库窗口
- 自动评分：`B` / `recommended`
- 证据路线：`wet_module_extension`
- pH 机制假设：moderate / exploratory
- 容量：raw `21,196,183`，usable non-control `21,195,927`
- 结构证据：`insufficient_or_unstable_contact; seed_support=20; cluster_support=25; ag_loop_fraction=0; a_determinant_fraction=0; quality_ok=False`
- 糖链风险：low heuristic risk
- 软风险：`sdab_V105H_D110H_negative_pair`

推荐理由：

- 覆盖 FR3 支撑区和 CDR3 主体。
- 覆盖主要湿实验 His 种子：Q100H、G102H、V105H、E108H、D110H、Y111H。
- 比 `VHH:57-96` 更直接利用已有 sdAb 湿实验信息。
- 容量远高于 10k，即使扣除 V105H-D110H 禁止组合后仍然充足。
- AF3 不提供直接稳定抗原 loop / A determinant 接触证据，因此该窗口应定位为湿实验驱动的探索性设计，而不是结构直接界面设计。

实际建库约束：

- 必须硬保护 A96C。
- 禁止 V105H 与 D110H 同时出现。
- 控制总 His 数量，避免 CDR3 密集 His 导致中性结合崩塌。
- 尽量加入非 His 救援或中性结合保留突变。

### 可接受替代首选

`sdAb_VHH_075_114 (VHH:75-114)`

- 审阅后角色：替代主探索库窗口
- 自动评分：`B` / `recommended`
- 证据路线：`wet_module_extension`
- 容量：raw `21,196,183`，usable non-control `21,195,927`
- 风险：同样必须禁止 V105H-D110H 组合。

推荐理由：

- 覆盖与 VHH:72-111 相同的主要 CDR3 His 种子。
- 相比 VHH:72-111 略向 C 端移动，包含更多 CDR3 尾部上下文。
- 同样没有 AF3 直接稳定界面证据，因此仍是探索性湿实验驱动窗口。

### 次级探索窗口

`sdAb_VHH_057_096 (VHH:57-96)`

- 审阅后角色：次级探索窗口，不作为主库首选
- 自动评分：`B` / `recommended`
- 证据路线：`wet_module_extension`
- 容量：raw `22,858,512`，usable non-control `22,858,256`
- 结构证据：100 个 AF3 模型下仍无直接稳定接触支持。

保留理由：

- 可用于测试 FR3/CDR2 或骨架-CDR 调控假设。
- 避开高风险 CDR3 His 种子。
- 但不覆盖主要湿实验 pH 敏感种子，因此不作为下一批主 10k 库的首选。

### sdAb 结论

sdAb 可以继续设计，但应明确定位为探索性库。审阅后首选 `VHH:72-111`，`VHH:75-114` 为接近的替代首选；`VHH:57-96` 降为次级探索窗口。

## 总体结论

- 1E62 是更强的直接推进对象。主窗口使用 `1E62_VL_001_040`。
- sdAb 仍可推进，但应以湿实验信息和 FR3/CDR3 调控假设为核心。主窗口采用 `sdAb_VHH_072_111`，或根据建库便利性采用相近的 `sdAb_VHH_075_114`。
- 正式进入序列库设计前，需要更新 mutation mask / library constraints：hard-protect `1E62 L23C` 和 `sdAb A96C`，并保持 `sdAb V105H-D110H` 禁止组合。
