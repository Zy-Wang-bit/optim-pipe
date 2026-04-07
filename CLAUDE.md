# optim-pipe

pH 依赖性抗体亲和力优化的自动化计算流水线。

## 环境

- **主环境**: `conda activate optim-pipe` (Python 3.11, PyTorch 2.8, fair-esm 2.0, BioPython 1.85)
- **ProteinMPNN**: `conda activate proteinmpnn`
- **SimpleFold**: `conda activate simplefold`
- **PyRosetta**: `conda activate pyrosetta`

## 目录结构

| 目录 | 职责 |
|------|------|
| `scripts/` | 核心 pipeline 脚本（4 Phase / 12 Step） |
| `configs/` | YAML 配置文件 |
| `analysis/` | 后分析代码（pKa、Rosetta、RMSD、ELISA 关联） |
| `experiments/` | 按抗体系统/轮次组织的实验数据 |
| `tier2/` | Tier 2 运行时输出（structures/、pka/、rosetta/、rmsd/） |
| `phase_c/` | 旧模式 Phase C 运行时输出 |
| `docs/` | 路线图、设计文档、分析报告 |
| `third_party/` | 外部工具（FoldX、ProteinMPNN、SimpleFold、MD） |
| `archive/` | 历史归档 |
| `.tasks/` | 跨对话任务管理 |

## Pipeline 结构 (3 Tier / 13 Step)

```
Tier 1: 生成 + 高通量筛选（Step 1-8）
  1. scan_interface.py         接口热点扫描
  2. run_mpnn_design.py        MPNN 设计 + His 种子
  3. run_esm_chunk.py          ESM 评分 + 候选筛选
  4. repair_pdbs.py            RepairPDB
  5. precompute_wt_ac.py       WT 基线
  6. make_mutlist_chunk.py     批次生成
  7. run_foldx_batch.py        BuildModel + AC
  8. tier1_filter.py           dG + delta 自适应门槛 → tier1_candidates.csv

Tier 2: 结构评估 + 筛选（Step 9-12，tier2.enabled 控制）
  PyRosetta 线（复合物模板）：
    9a. build_structures.py    突变体建模
    9b. run_pka.py             pKa 预测
    9c. run_rosetta_eval.py    dddG_elec + pH-score
  SimpleFold 3x 线（单抗体序列）：
    10. run_simplefold_3x.py   3x 采样预测
    11. run_rmsd.py            CDR RMSD（3x 中位数 + 异常值剔除）
  12. tier2_filter.py          pKa 相对排序 + RMSD 门槛 → tier2_candidates.csv

Tier 3: 精排（Step 13）
  13. merge_and_rank.py        dddG_elec 排序 + 软标记 → final_candidates.csv
```

`tier1.enabled: false` 时回退到旧 4 Phase 流程（merge_and_select.py），行为与原 8 步 pipeline 完全一致。

## 任务管理

### 何时使用 .tasks/ 系统

**使用完整任务流程**（创建 `.tasks/active/<name>/` 目录）的条件——满足任意一条即可：
- 预计跨越多个对话才能完成
- 涉及多步骤实现（如新建 pipeline 模块、设计新一轮实验方案）
- 需要记录研究发现供后续对话参考

**不需要任务流程**的情况：
- 单次对话内可完成的小任务（改个 bug、加个参数、跑个脚本）
- 问答和讨论（解释代码、分析数据、比较方案）
- 一次性的文件编辑或配置修改

遇到不确定的情况，直接做。如果做到一半发现比预期复杂，再补建任务文件。

### 新对话启动
1. 检查 `.tasks/active/` 是否已有对应的活跃任务
2. 如果有，读取其 `plan.md` 和 `progress.md`，用"五问重启检查"恢复上下文
3. 如果没有，正常响应用户请求

### 完整任务流程（仅限大型任务）

**创建**：在 `.tasks/active/<task-name>/` 创建 `plan.md`、`findings.md`、`progress.md`，与用户确认后开始。

**执行**：每完成一个阶段更新 `progress.md`；重要发现记入 `findings.md`；做重大决策前重读 `plan.md`。

**完成**：将任务目录从 `active/` 移到 `done/`。

### 五问重启检查（恢复活跃任务时使用）

| 问题 | 说明 |
|------|------|
| 我在哪里？ | 当前阶段 |
| 我要去哪里？ | 下一阶段 |
| 目标是什么？ | 任务目标 |
| 我学到了什么？ | 见 findings.md |
| 我做了什么？ | 见 progress.md |

## 编码规范

- 配置路径统一使用 `configs/config.yaml`，不硬编码绝对路径
- 实验数据放 `experiments/<系统>/<轮次>/`，分析代码放 `analysis/`
- 新增外部工具放 `third_party/`
- Pipeline 运行时输出（results/, foldx/, phase_c/, logs/ 等）由 .gitignore 管理，不提交

## 注意事项

- 所有脚本必须从仓库根目录运行（使用相对路径如 `mpnn_outputs/`、`foldx/`）
- FoldX 要求 PDB 文件复制到工作目录，脚本已自动处理
- 突变命名有多种格式（R1 `A_E1H`、R3 `H1E>H`、FoldX `SA40A`）。最终数据文件统一为 `HE1H` 格式；可按“突变后氨基酸 + 链/位点 + 突变后氨基酸”的方式理解。具体示例：同一个 E1→H 突变在 R1 中可写作 `A_E1H`，在 R3 中可写作 `H1E>H`，统一后记为 `HE1H`；具体转换规则见 `analysis/naming/convert.py`
- 多 PDB 模板设计：per-template 配置在 `configs/config_foldx_n*.yaml`
- ProteinMPNN 和 SimpleFold 通过外部子进程调用，分别需要 `proteinmpnn` 和 `simplefold` conda 环境
- Tier 2 Step 9a/9c 需要 `pyrosetta` 环境，Step 10 需要 `simplefold` 环境，其余在主环境 `optim-pipe` 中运行
- Tier 配置在 `configs/config.yaml` 的 `tier1`/`tier2`/`tier3` 段；旧模式配置在 `phase_c` 段
