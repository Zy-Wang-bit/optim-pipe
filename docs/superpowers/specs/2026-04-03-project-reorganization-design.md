# optim-pipe 项目重组设计

> **日期**: 2026-04-03
> **状态**: 待实施
> **范围**: 目录结构全面重组 + 跨对话任务管理机制 + CLAUDE.md 工作协议

---

## 背景

项目经过 R1→R2→R3 三轮迭代，根目录积累了 28 个条目（18 目录 + 10 文件），存在以下问题：

1. **冗余文件**：`results/`（153M 历史中间数据）、`nano01/`（废弃模板）、`molecules/`（空目录）
2. **职责不清**：`tools/` 和 `third_party/` 分散存放外部依赖；`data/` 和 `experiments/` 边界模糊；`artifacts/` 独立于 SimpleFold 存放其权重
3. **根目录杂乱**：设计文档、路线图、临时任务文件散落在根目录
4. **缺乏跨对话连续性**：`task_plan.md` / `progress.md` / `findings.md` 在根目录且无组织规范，新对话难以恢复上下文

**目标**：根目录可见项从 28 项减至 ~11 项；建立结构化的任务管理机制；将工作协议写入 CLAUDE.md。

---

## 一、目录结构重组

### 1.1 目标根目录结构

```
optim-pipe/
├── README.md                # 项目总览（需更新路径引用）
├── CLAUDE.md                # AI 工作协议 + 跨对话机制
├── run_pipeline.sh          # 主 pipeline 入口
│
├── configs/                 # 配置文件
├── scripts/                 # 核心 pipeline 脚本（12 个模块）
├── analysis/                # 后分析代码（移走报告文件，只留代码）
│
├── experiments/             # 所有实验数据（重组，见 1.2）
├── docs/                    # 文档（新建，见 1.3）
├── third_party/             # 所有外部工具和依赖（合并，见 1.4）
│
├── archive/                 # 历史归档（保留）
├── workspace/               # FoldX 工作区（暂保留，待确认）
│
├── .tasks/                  # 跨对话任务管理（新建，见第二节）
├── .claude/                 # Claude 配置
├── .agents/                 # Agent 规则
└── .vscode/                 # VS Code 配置
```

### 1.2 experiments/ 重组

**核心变更**：为 1E62 建立上级目录，三个轮次改名 R1/R2/R3；`data/` 中的共享输入数据迁移到 `experiments/1E62/data/`。

```
experiments/
├── 1E62/                       # 1E62 抗体系统（主线）
│   ├── data/                   # 共享输入数据
│   │   ├── AF3-abag-n1.pdb     # 复合物结构（仅保留 *_Repair.pdb，重命名）
│   │   ├── AF3-abag-n2.pdb
│   │   ├── AF3-abag-n3.pdb
│   │   ├── AF3-abag-n4.pdb
│   │   ├── ab_wt.pdb           # 野生型抗体结构（原 data/wt/ab_wt.pdb）
│   │   ├── heavy.fasta         # 野生型重链序列
│   │   ├── light.fasta         # 野生型轻链序列
│   │   └── antigen.fasta       # 抗原序列
│   ├── R1/                     # 原 1E62_R1/
│   ├── R2/                     # 原 1E62_R2/（含 wet_lab/, results/, scripts/）
│   └── R3/                     # 原 1E62_R3/（含 design/, scripts/, results/）
│
├── sdab/                       # 单域抗体（保持不变）
├── 3E72/                       # 搁置中，数据有价值，保留
├── 16F9/                       # 搁置中，数据有价值，保留
└── ext_gpt/                    # GPT 辅助设计（保留）
```

**数据迁移明细**：

| 源 | 目标 | 说明 |
|----|------|------|
| `data/pdbs/n1-0_Repair.pdb` | `experiments/1E62/data/AF3-abag-n1.pdb` | 重命名 |
| `data/pdbs/n2-0_Repair.pdb` | `experiments/1E62/data/AF3-abag-n2.pdb` | 重命名 |
| `data/pdbs/n3-0_Repair.pdb` | `experiments/1E62/data/AF3-abag-n3.pdb` | 重命名 |
| `data/pdbs/n4-0_Repair.pdb` | `experiments/1E62/data/AF3-abag-n4.pdb` | 重命名 |
| `data/pdbs/n*-0.pdb`（非 Repair） | **删除** | 仅保留 Repair 版本 |
| `data/wt/ab_wt.pdb` | `experiments/1E62/data/ab_wt.pdb` | 野生型抗体结构 |
| `data/wt/heavy.fasta` | `experiments/1E62/data/heavy.fasta` | 野生型重链序列 |
| `data/wt/light.fasta` | `experiments/1E62/data/light.fasta` | 野生型轻链序列 |
| `data/wt/antigen.fasta` | `experiments/1E62/data/antigen.fasta` | 抗原序列 |
| `data/wt/original.fasta` | `experiments/1E62/data/original.fasta` | 原始序列 |
| `data/wt/wt_binding_energy.csv` | `experiments/1E62/data/wt_binding_energy.csv` | WT 结合能数据 |
| `data/pdb_AF3/` | **删除** | 不再需要 |
| `data/molecular_dynamics/` | `third_party/molecular_dynamics/` | 分子模拟脚本 |

### 1.3 docs/ 结构

```
docs/
├── ROADMAP.md                          # 路线图（原根目录）
├── design/                             # 设计文档
│   ├── Phase1-combination_design.md    # R3 组合设计方案
│   └── Phase1-binary_classification.md # ESM-2 分类模型方案
└── reports/                            # 分析报告
    ├── correlation_report.md           # 原 analysis/
    └── deep_metric_analysis.md         # 原 analysis/
```

### 1.4 third_party/ 合并

**核心变更**：将 `tools/foldx/` 和 `artifacts/` 合并进 `third_party/`，所有外部依赖统一管理。

```
third_party/
├── ProteinMPNN/                # 序列设计（不变）
├── ml-simplefold/              # 结构预测
│   ├── ...                     # 原有代码
│   └── weights/                # 新建子目录
│       ├── simplefold_3B.ckpt      # 原 artifacts/
│       ├── simplefold_1.6B.ckpt    # 原 artifacts/
│       ├── plddt_module_1.6B.ckpt  # 原 artifacts/
│       └── plddt.ckpt              # 原 artifacts/
├── foldx/                      # 原 tools/foldx/
│   └── foldx                   # FoldX 二进制
└── molecular_dynamics/         # 原 data/molecular_dynamics/
    ├── *.mdp                   # GROMACS 参数文件
    └── autoGromacs.sh          # 自动化脚本
```

### 1.5 analysis/ 简化

移走报告文件后，只保留代码：

```
analysis/
├── pka/                        # pKa 预测
│   └── run_pka.py
├── rosetta/                    # pH-score 与 dddG_elec 计算
│   ├── batch_calc_ph_score.py
│   ├── batch_calc_dddg_elec.py
│   ├── calc_ph_score_antibody.py
│   ├── calc_dddg_elec.py
│   └── rank_mutations.py
├── structure_compare/          # RMSD 比较
│   └── rmsd_ca_global.py
└── elisa_vs_computation.py     # ELISA 关联分析
```

### 1.6 删除清单

| 目标 | 大小 | 原因 |
|------|------|------|
| `results/` | 153M | 历史中间数据，无保留价值 |
| `nano01/` | 2.8M | 废弃的纳米抗体结构模板 |
| `molecules/` | 0 | 空目录 |
| `data/pdb_AF3/` | ~50M | 不再需要 |
| `data/` 目录本身 | — | 内容迁移后删除空壳 |
| `tools/` 目录本身 | — | 合并到 third_party 后删除 |
| `artifacts/` 目录本身 | — | 移入 ml-simplefold/weights/ 后删除 |

### 1.7 需要更新的代码和配置引用

| 文件 | 变更 |
|------|------|
| `configs/config.yaml` | `paths.foldx_bin`: `tools/foldx/foldx` → `third_party/foldx/foldx` |
| `configs/config.yaml` | `paths.pdb_dir` 等: `data/pdbs/` → `experiments/1E62/data/` |
| `configs/config.yaml` | `paths.wt_*`: `data/wt/` → `experiments/1E62/data/` |
| `configs/config_foldx_*.yaml` | 同上，各 PDB 专用配置中的路径 |
| `configs/config_esm_*.yaml` | ESM 配置中的路径引用 |
| SimpleFold 加载代码 | 权重路径: `artifacts/` → `third_party/ml-simplefold/weights/` |
| `run_pipeline.sh` | 若有硬编码路径需更新 |
| `README.md` | 目录结构说明全面更新 |

---

## 二、跨对话任务管理机制

### 2.1 目录结构

```
.tasks/
├── active/                 # 当前活跃任务
│   └── <task-name>/        # 每个任务一个目录
│       ├── plan.md         # 任务计划
│       ├── findings.md     # 发现与决策
│       └── progress.md     # 执行进度
└── done/                   # 已完成任务归档
    └── <task-name>/        # 完成后从 active/ 移入
```

### 2.2 文件格式规范

**plan.md**（任务计划）：
```markdown
# 任务: <任务名称>

## 目标
<一句话描述目标>

## 当前阶段
阶段 X

## 阶段
### 阶段 1: <名称>
- [ ] 具体检查项
- **状态:** pending | in_progress | complete

### 阶段 2: <名称>
...

## 关键问题
1. <问题> → <解答或待定>

## 已做决策
| 决策 | 理由 |
|------|------|
```

**progress.md**（执行进度）：
```markdown
# 进度日志

## 会话: YYYY-MM-DD

### 阶段 X: <名称>
- **状态:** in_progress | complete
- 已执行操作:
  - <操作描述>
- 创建/修改的文件:
  - <文件路径> (<说明>)

## 五问重启检查
| 问题 | 答案 |
|------|------|
| 我在哪里？ | 当前阶段 |
| 我要去哪里？ | 下一阶段 |
| 目标是什么？ | 任务目标 |
| 我学到了什么？ | 见 findings.md |
| 我做了什么？ | 见上方记录 |
```

**findings.md**（发现与决策）：
```markdown
# 发现与决策

## 需求
<需求描述>

## 研究发现
### <发现主题>
<详细内容>

## 技术决策
| 决策 | 理由 |
|------|------|

## 资源
- <关键文件或工具路径>
```

### 2.3 任务生命周期

```
新任务 → .tasks/active/<name>/ 创建 plan.md
       → 与用户确认计划
       → 执行，更新 progress.md / findings.md
       → 完成 → 移动到 .tasks/done/<name>/
```

### 2.4 迁移现有任务文件

| 源 | 目标 |
|----|------|
| 根目录 `task_plan.md` | `.tasks/done/sdab-pka/plan.md` |
| 根目录 `progress.md` | `.tasks/done/sdab-pka/progress.md` |
| 根目录 `findings.md` | `.tasks/done/sdab-pka/findings.md` |

---

## 三、CLAUDE.md 工作协议

写入以下内容到项目根目录的 `CLAUDE.md`：

### 3.1 项目概述段

- 项目名称、用途、当前状态概要
- conda 环境说明
- 关键目录职责说明

### 3.2 跨对话任务协议

1. **新对话启动** — 读取 `.tasks/active/` 下所有活跃任务的 `plan.md` 和 `progress.md`，理解当前上下文
2. **接受新任务** — 在 `.tasks/active/<task-name>/` 创建三个文件，与用户确认 plan.md 后开始执行
3. **执行中** — 每完成一个阶段更新 progress.md；重要发现记入 findings.md；做重大决策前重新读取 plan.md
4. **任务完成** — 将任务目录从 `active/` 移到 `done/`
5. **跨对话恢复** — 新对话通过 plan.md 的"当前阶段"+ progress.md 的最新会话记录恢复上下文；用"五问重启检查"验证理解

### 3.3 编码规范

- 配置路径统一使用 `configs/config.yaml`，不硬编码
- 实验数据放 `experiments/<系统>/<轮次>/`，分析代码放 `analysis/`
- 新增外部工具放 `third_party/`

---

## 四、实施顺序

建议分步实施以降低风险：

1. **创建目标目录** — `docs/`, `.tasks/`, `third_party/foldx/`, `third_party/molecular_dynamics/`, `experiments/1E62/`, `experiments/1E62/data/`, `third_party/ml-simplefold/weights/`
2. **迁移文件**（无破坏性，先复制后验证）
   - 文档文件 → `docs/`
   - 任务文件 → `.tasks/done/sdab-pka/`
   - `data/wt/*`, `data/pdbs/*_Repair.pdb` → `experiments/1E62/data/`
   - `tools/foldx/` → `third_party/foldx/`
   - `artifacts/*.ckpt` → `third_party/ml-simplefold/weights/`
   - `data/molecular_dynamics/` → `third_party/molecular_dynamics/`
   - `1E62_R1/R2/R3` → `experiments/1E62/R1/R2/R3`
   - `analysis/*.md` → `docs/reports/`
3. **更新代码引用** — 修改 config.yaml, SimpleFold 权重加载路径, run_pipeline.sh, README.md
4. **验证** — 运行 pipeline 核心步骤验证路径正确
5. **清理** — 确认迁移无误后删除源文件和空目录（`results/`, `nano01/`, `molecules/`, `data/`, `tools/`, `artifacts/`）
6. **写入 CLAUDE.md** — 工作协议最后写入，确保引用的路径都是最终版本

---

## 五、验证方案

| 验证项 | 方法 |
|--------|------|
| 路径引用完整性 | `grep -r "data/pdbs\|data/wt\|tools/foldx\|artifacts/" scripts/ configs/ analysis/ run_pipeline.sh` 应无匹配 |
| FoldX 可用 | `third_party/foldx/foldx --help` 正常输出 |
| SimpleFold 权重 | 确认 `third_party/ml-simplefold/weights/` 下 4 个文件存在且大小正确 |
| 配置文件 | `python -c "import yaml; yaml.safe_load(open('configs/config.yaml'))"` 无报错 |
| 目录结构 | 根目录 `ls` 输出不超过 15 项 |
| 任务系统 | `.tasks/done/sdab-pka/` 包含完整的 plan.md, progress.md, findings.md |
| CLAUDE.md | 非空，包含项目概述和任务协议 |
