# 512hs121-152 sdab 初筛回溯

数据源：
- 交付清单：`docs/wet_lab_delivery/sdab_variants.csv`
- 检测结果：`docs/wet_lab_delivery/5125-hs121-152-初筛elisa-GAH-HRP检测.xlsx`
- R4 计算指标：`experiments/sdab/R4/structures/all_metrics.csv`、`experiments/sdab/R4/final/R4_delivery.csv`

## 一、检测结果对齐

- XLSX 中有检测行的编号：`512hs121`、`512hs123-152`。`512hs122` 在表内未找到检测行。
- 本报告用四个平均 OD 的最大值作为粗略 activity：`max(Ce pH6, Ce pH7.4, D1 pH6, D1 pH7.4)`。
- 粗略分级：strong `>=1.0`，weak `0.1-1.0`，background `<0.1`。

| 分组 | 数量 | strong | weak | background | missing |
|---|---:|---:|---:|---:|---:|
| S1_MPNN_R2 | 12 | 0 | 0 | 11 | 1 |
| S2_R4_CDR_only | 10 | 1 | 3 | 6 | 0 |
| S3_R3_FR_His | 10 | 0 | 0 | 10 | 0 |

## 二、Top activity 样本

| ID | variant | section | max OD | pH7.4 max | pH6 max | Ce ratio | D1 ratio | His mutations |
|---|---|---|---:|---:|---:|---:|---:|---|
| 512hs135 | r4_09 | S2_R4_CDR_only | 1.764 | 1.764 | 1.697 | 1.31 | 0.84 | TA28H; QA100H; GA102H; DA110H |
| 512hs142 | r4_14 | S2_R4_CDR_only | 0.962 | 0.962 | 0.851 | 1.44 | 1.13 | RA27H; NA54H; GA102H; DA110H |
| 512hs140 | r4_15 | S2_R4_CDR_only | 0.391 | 0.391 | 0.367 | 1.30 | 0.78 | SA31H; QA100H; GA102H; DA110H |
| 512hs141 | r4_13 | S2_R4_CDR_only | 0.367 | 0.193 | 0.367 | 0.69 | 0.53 | TA28H; NA54H; GA102H; DA110H |
| 512hs133 | r4_04 | S2_R4_CDR_only | 0.062 | 0.062 | 0.009 | 7.19 | 1.52 | GA53H; GA55H; QA100H; VA105H |
| 512hs134 | r4_01 | S2_R4_CDR_only | 0.051 | 0.051 | 0.008 | 6.61 | 1.81 | GA55H; GA56H; QA100H; VA105H |
| 512hs130 | v09_minimal | S1_MPNN_R2 | 0.025 | 0.025 | 0.011 | 2.27 | 1.26 | NA52H; DA99H; DA110H |
| 512hs147 | sdab_r3_A5 | S3_R3_FR_His | 0.023 | 0.023 | 0.013 | 1.52 | 1.70 | A61H; V93H; Y94H; D99H; V105H; D110H |

## 三、主要回溯结论

1. **R4 Track B 是唯一保留明显结合的设计族。** `hs135/r4_09` 是唯一 strong；`hs140-142` 都是 Track B 的 `G102H + D110H` 族，表现为弱到中等结合。
2. **R4 Track A 系统性失活。** `hs133/134/136/137/138/139` 都含 `Q100H + V105H + G55H` 主干，最终 MD composite 排名靠前，但实测 OD 全在背景区或接近背景。
3. **第一节 MPNN 全序列设计和第三节 FR-His 机制探索基本全灭。** 前者包含 21-51 个突变，后者把 His 放入 FR3；两者的 MD pH-contact 信号未能保证 pH7.4 结合保留。
4. **现有 active 样本没有达到理想 pH6 解离。** `hs135` 结合强，但 Ce ratio 只有 1.31，D1 ratio 为 0.84；`hs140-142` ratio 也不稳定，更多是“保留结合”而非“强 pH-switch”。

## 四、R4 失配来源

### 4.1 ElasticNet 对 Track A 的外推越界

`Q100H + V105H` 在训练集中只出现过以下记录，且都同时含 `D110H`：

| training ID | mutations | log_pH74 | log_pH6 | log_ratio |
|---|---|---:|---:|---:|
| hs63 | HQ100H;HV105H;HD110H | -0.314 | 0.718 | -1.032 |
| hs86 | HQ100H;HG102H;HV105H;HD110H | 0.013 | 0.464 | -0.451 |

也就是说，模型把 `Q100H-V105H` 学成强正 pair，但没有 D110H-free 的直接训练样本支撑。候选枚举又硬过滤掉 `V105H + D110H`，于是把这个 pair 外推到 `Q100H+V105H+G55H` 这条新主干。湿实验显示这条外推不成立。

### 4.2 Track B 有更直接的训练支撑

`Q100H + G102H + D110H` 在训练集中有直接或近邻组合：

| training ID | mutations | log_pH74 | log_pH6 | log_ratio |
|---|---|---:|---:|---:|
| hs62 | HQ100H;HG102H;HD110H | 1.022 | 0.495 | 0.527 |
| hs86 | HQ100H;HG102H;HV105H;HD110H | 0.013 | 0.464 | -0.451 |
| hs87 | HQ100H;HG102H;HE108H;HD110H | 0.520 | 0.247 | 0.273 |
| hs88 | HQ100H;HG102H;HD110H;HY111H | 0.516 | 0.095 | 0.421 |

当前 strong 样本 `r4_09 = T28H + Q100H + G102H + D110H` 正好是在 `hs62 = Q100H + G102H + D110H` 上加了 `T28H`。这比 Track A 的主干更接近已有实验支持。

### 4.3 最终 MD composite 过度看重“pH6 变松”，缺少绝对结合约束

最终 `analysis/r4/final_select.py` 只用 MD delta 指标打 composite，`model_score / delta_delta / dddG_elec / ph_score` 没有进入最终复合分。对这次实测的 R4 10 个样本，Spearman 方向大致为：

| 指标 | vs activity_max | 解读 |
|---|---:|---|
| model_score | -0.612 | 模型最高的 Track A 实测失活 |
| delta_delta | +0.624 | FoldX pH 差分更偏向实际保留结合的 Track B |
| dddG_elec | +0.663 | Rosetta 静电信号也偏向 r4_09/r4_15 |
| delta_rmsf_h3 | -0.442 | MD CDR3 柔性越强，实测反而越差 |
| delta_buried_sasa | +0.394 | Track B 的界面收紧/不打开反而对应更好结合 |

这说明当前 MD 指标更像是在回答“已结合复合物在 pH6 是否松动”，没有可靠回答“pH7.4 是否还能形成足够结合”。

## 五、建议修正

1. 之后的模型和枚举必须加 **support gate**：pair 或高阶 motif 只能在与训练集相同背景或近邻背景中使用；`Q100H+V105H` 这类只在 D110H 背景出现过的 pair，不能直接推到 D110H-free 主干。
2. 下一轮以 `r4_09/hs135` 为唯一 strong parent；优先做 `Q100H+G102H+D110H`、`T28H+Q100H+G102H+D110H`、以及 T28/S31/R27/N54 的单替换/删减对照。
3. `Q100H+V105H+G55H` Track A 主干暂时降级为高风险失活 motif，除非补做 D110H-free 训练样本证明它有效。
4. 最终筛选要把“绝对 pH7.4 结合保留”作为硬门槛或主目标；MD delta 只能作为 pH-switch 辅助信号。
5. 第一节和第三节的 MD contact/RMSF 指标在这批数据中不可作为湿实验优先级依据；需要加表达/折叠/绝对结合的早期过滤，尤其限制大规模框架重设计和 FR3 疏水核心 His 化。

合并明细表：`docs/wet_lab_delivery/sdab_512hs121_152_retrospective.csv`
