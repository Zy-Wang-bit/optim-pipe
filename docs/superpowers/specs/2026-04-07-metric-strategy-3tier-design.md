# 指标利用策略 + 3 Tier Pipeline 重构

> **日期**: 2026-04-07
> **状态**: 设计完成，待实现
> **范围**: Pipeline 结构从 4 Phase / 12 Step 重构为 3 Tier / 13 Step

---

## 1. 背景与动机

Pipeline 4 Phase / 12 Step 已全部跑通，但存在两个核心问题：

1. **指标角色不清**：14 个指标（含衍生指标）分散在 Phase B/C 中，评分公式仅用 Phase B 指标，Phase C 指标未参与筛选。当前加权评分公式有本质问题——指标与活性不显著相关时硬加权等于引入噪声。
2. **Pipeline 结构是过渡方案**：Phase C 职责过杂（两套结构生成 + pKa + dddG_elec + RMSD），筛选逻辑不明确。

本次重构将 pipeline 改为 **3 Tier 漏斗结构**，每个指标有明确的角色定位（门槛 / 排序 / 软标记），取代旧的加权评分公式。

---

## 2. 核心问题与指标分类

所有指标服务于两个核心问题：

### 问题 1：pH 7.4 结合活性 — "突变后还能不能结合？"

| 指标 | 角色 | 使用位置 | 验证依据 |
|------|------|---------|---------|
| dG_pH7.4 / ddG_pH7.4 | Tier 1 硬门槛 | Step 8 | R2 AUC=0.854 (p=0.007)，最佳单一预测因子 |
| CDR RMSD (SimpleFold 3x) | Tier 2 硬门槛 | Step 12 | R2 H1 AUC=0.818 (p=0.027)，需 3x 采样 |
| ESM | 软标记 | Step 8, 13 | R2 未见显著相关，作为序列可信度参考 |
| pH-score | 软标记 | Step 13 | 结构稳定性参考 |

### 问题 2：pH 敏感性 — "酸性下能不能解离？"

| 指标 | 角色 | 使用位置 | 验证依据 |
|------|------|---------|---------|
| delta | Tier 1 硬门槛 (>0) | Step 8 | R2 com18（最优）delta 最大 |
| pKa (propka/pkai) | Tier 2 相对排序筛选 | Step 12 | 绝对值有系统偏差，但趋势可靠 |
| dddG_elec | Tier 3 精排 | Step 13 | 物理机制最明确（His 质子化静电效应） |
| consensus | 软标记 | Step 13 | pKa 两工具一致性，可信度参考 |

### 淘汰/降级的指标

| 指标 | 处理 | 原因 |
|------|------|------|
| delta_delta | 淘汰 | 与 delta 差常数（delta_WT），单模板冗余 |
| dG_pH6.0 / ddG_pH6.0 | 中间量 | delta 的计算组件，不单独使用 |
| delta_propka/pkai | 保留计算，不做筛选 | 参考用 |
| shift_propka/pkai | 保留计算，不做筛选 | 新引入 His 时 WT 无对应残基 |
| global_rmsd | 保留计算，不做筛选 | CDR RMSD 更有区分力 |

---

## 3. Pipeline 结构：3 Tier / 13 Step

### Tier 1：生成 + 高通量筛选（~10万 → ~500）

| Step | 脚本 | 功能 | 环境 |
|------|------|------|------|
| 1 | scan_interface.py | 接口热点扫描 | optim-pipe |
| 2 | run_mpnn_design.py | MPNN 设计 + His 种子 | proteinmpnn |
| 3 | run_esm_chunk.py | ESM 评分 | optim-pipe |
| 4 | repair_pdbs.py | FoldX RepairPDB | optim-pipe |
| 5 | precompute_wt_ac.py | WT 基线能量 | optim-pipe |
| 6 | make_mutlist_chunk.py | FoldX 批次生成 | optim-pipe |
| 7 | run_foldx_batch.py | FoldX BuildModel + AnalyseComplex (pH 7.4/6.0) | optim-pipe |
| **8** | **tier1_filter.py** (NEW) | **自适应门槛筛选** | optim-pipe |

**Step 8 tier1_filter.py 筛选逻辑：**

- 硬门槛：`dG_pH7.4 < 阈值`（排除结合太差的候选）
- 硬门槛：`delta > 0`（排除无 pH 切换方向的候选）
- 软标记：ESM 异常标记（不淘汰，标记 `esm_flag`）
- 自适应调节：通过数 > 目标则收紧阈值，< 目标则放宽，直到候选量落入目标范围（默认 300-500）
- FoldX 中间文件清理：BuildModel 产生的突变体 PDB 和 fxout 在 Step 7 完成后清理，只保留 summary CSV
- 输出：`tier1_candidates.csv`

### Tier 2：结构评估 + 门槛筛选（~500 → ~100）

两条独立结构评估线，可并行执行：

**PyRosetta 线**（基于 AF3 复合物模板）：

| Step | 脚本 | 功能 | 环境 |
|------|------|------|------|
| 9a | build_structures.py | AF3 复合物 → PyRosetta 引入突变 | pyrosetta |
| 9b | run_pka.py | PROPKA3 + pKAI+ pKa 预测 | optim-pipe |
| 9c | run_rosetta_eval.py | dddG_elec + pH-score | pyrosetta |

**SimpleFold 3x 线**（基于单抗体序列，从头预测）：

| Step | 脚本 | 功能 | 环境 |
|------|------|------|------|
| 10 | run_simplefold_3x.py (NEW) | SimpleFold 3B 每变体 3 次采样 | simplefold |
| 11 | run_rmsd.py | 全局 + CDR RMSD (H1-3, L1-3) | optim-pipe |

**汇合筛选：**

| Step | 脚本 | 功能 | 环境 |
|------|------|------|------|
| **12** | **tier2_filter.py** (NEW) | **门槛筛选** | optim-pipe |

**Step 12 tier2_filter.py 筛选逻辑：**

- pKa 相对排序：按 pKa 相对于 WT/基底抗体的提升排序，取 top-N。pKa 绝对值有系统偏差，但趋势可靠——筛出比基线好的变体即可
- 硬门槛：CDR RMSD < 阈值（结合活性评估第二关，剔除结构变形过大的候选）
- 输出：`tier2_candidates.csv`（~50-100 候选）

### Tier 3：精排 + 最终输出（~100 → top-N）

| Step | 脚本 | 功能 | 环境 |
|------|------|------|------|
| **13** | **merge_and_rank.py** (重写) | **最终合并排序** | optim-pipe |

**Step 13 merge_and_rank.py 逻辑：**

- 排序：按 dddG_elec 降序（值越大 = pH 切换越强）
- 软标记：ESM 异常、pH-score 异常、pKa consensus 不一致
- 输出：`final_candidates.csv`（全部候选 + 排序 + 标记）

---

## 4. 结构策略

- **两套模板结构**：用户提供 AF3 复合物和 AF3 单抗体两个预测结构
- **PyRosetta 线用复合物模板**：pKa / dddG_elec / pH-score 需要 ab-ag 界面信息
- **SimpleFold 3x 线用单抗体序列**：从头预测，仅评估抗体 CDR 构象变化
- **Tier 1 FoldX 保持现有流程**：R2 验证表明 FoldX 在自身构建的结构上效果最佳（AUC=0.854），跨工具结构分析效果下降（R2 报告 Section 2.4.3）

---

## 5. 配置结构变化

`configs/config.yaml` 调整：

- 新增 `tier1` 段：门槛阈值、自适应目标范围
- 新增 `tier2` 段：pKa 基线、RMSD 门槛阈值
- `phase_c` 段重命名为 `tier2`，拆分 PyRosetta 线和 SimpleFold 3x 线的配置
- `selection` 段简化为 Tier 3 排序配置
- 移除 `delta_delta_weight`、`his_bonus` 等旧评分参数

---

## 6. 新增/修改文件

| 文件 | 变化 | 说明 |
|------|------|------|
| `scripts/tier1_filter.py` | NEW | Tier 1 自适应门槛筛选 |
| `scripts/run_simplefold_3x.py` | NEW | SimpleFold 3x 采样封装 |
| `scripts/tier2_filter.py` | NEW | Tier 2 门槛筛选（pKa 相对排序 + RMSD） |
| `scripts/merge_and_rank.py` | NEW (重写) | Tier 3 精排 + 最终输出 |
| `run_pipeline.sh` | 修改 | 适配 3 Tier 流程 |
| `configs/config.yaml` | 修改 | 适配新配置结构 |
| `scripts/merge_and_select.py` | 保留 | 向后兼容（旧模式时沿用） |
| `CLAUDE.md` | 修改 | 更新 pipeline 结构描述 |

---

## 7. 向后兼容

- `merge_and_select.py` 保留，旧配置（`phase_c.enabled: false`）仍可运行原 8 步 pipeline
- 新 3 Tier 流程通过新的 `run_pipeline.sh` 入口启动

---

## 8. 验证方式

1. 旧配置模式行为不变（向后兼容）
2. 使用 R3 10 变体数据端到端测试 3 Tier 流程
3. 检查 `tier1_candidates.csv`、`tier2_candidates.csv`、`final_candidates.csv` 的候选数量和指标列完整性
4. 验证 pKa 相对排序输出的合理性（比基线差的变体应被排除）

---

## 9. 未来方向（不在本次范围）

- ML 模型：结合实验数据训练预测模型，替代复杂筛选流程
- 多模板支持：delta_delta 在多 PDB 模板场景下的归一化价值
- Pareto 前沿：多目标优化替代单指标排序
