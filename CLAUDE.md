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
| `scripts/` | 核心 pipeline 脚本（8 步） |
| `configs/` | YAML 配置文件 |
| `analysis/` | 后分析代码（pKa、Rosetta、RMSD、ELISA 关联） |
| `experiments/` | 按抗体系统/轮次组织的实验数据 |
| `docs/` | 路线图、设计文档、分析报告 |
| `third_party/` | 外部工具（FoldX、ProteinMPNN、SimpleFold、MD） |
| `archive/` | 历史归档 |
| `.tasks/` | 跨对话任务管理 |

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
- Pipeline 运行时输出（results/, foldx/, logs/ 等）由 .gitignore 管理，不提交
