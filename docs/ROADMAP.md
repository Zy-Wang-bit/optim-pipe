# 后续工作路线图

> 基于湿实验第二轮（R2）结果的规划，详见 [analysis_report.md](experiments/1E62_R2/wet_lab/analysis_report.md)

---

## 一、突变组合设计 ✅

> 目标：基于 R2 已验证的活性突变模块进行组合，快速产出下一批实验变体。

- ✅ R3 10 变体已设计并完成计算评估
- 评估报告：`experiments/1E62/R3/R3_evaluation_report.md`

## 二、结构评估能力提升 ✅

> 核心问题：如何通过结构预测判断突变体的结合活性变化？

- ✅ R2 验证完成：SimpleFold 3x 采样 + CDR RMSD，H1 RMSD AUC=0.818 (p=0.027)
- ✅ 已集成到 pipeline：Tier 2 SimpleFold 3x 线（Step 10-11），CDR RMSD 作为 Tier 2 硬门槛
- ✅ 3x 采样中位数 + VH-VL 朝向异常值剔除（global RMSD > 2.0Å）

## 三、pH 敏感性评估方法改进 ✅

- ✅ PROPKA3 + pKAI+ 双工具部署完成
- ✅ dddG_elec（Rosetta fa_elec pH 依赖静电）已实现并优化（127x 加速）
- ✅ 指标角色明确：delta（Tier 1 粗筛）→ pKa（Tier 2 相对排序）→ dddG_elec（Tier 3 精排）
- ✅ 3 Tier 漏斗 pipeline 已集成全部 pH 敏感性指标

## 四、单域抗体（sdab）突变设计

> 同一抗原的另一条改造线，目标同为酸性解离。在第三章 pH 评估方法改进完成后启动。

- **现有数据基础**（`experiments/sdab/`）：
  - CDR 区全位点 His scanning 已完成（122 个单点突变）
  - 实验数据：两个 pH 下的 KD 值（`20251203 SNB512-HS1-31 HBS-D.xlsx`）
  - 计算结果：composite rank、ddG、pH-score、ESM score（`sdab_ranked_all.csv`）
  - 结构文件 `sdab.pdb` + FoldX repaired 结构
- **数据特点**：CDR 区突变对结合活性影响较小，主要关注 KD(pH6)/KD(pH7) 比值（pH 解离效果）
- **设计策略**：
  - 从 scanning 数据中筛选 pH 解离效果最佳的位点
  - 复用改进后的 pH 评估工具（pKa 预测、改进的 pH-score）对已有 rank 重新评估
  - 将候选位点进行组合设计；sdab 为单链结构，组合自由度更高（无 H/L chain 配对约束）
- **与 1E62 的协同**：
  - pipeline 工具链完全共享
  - 注意：sdab 标签为连续 KD ratio，与 1E62 的二元活性标签不同，模型训练需做标签对齐（如 KD ratio → pH 敏感/不敏感伪标签）或在排序模型阶段统一处理

## 五、数据驱动与模型训练

> 双重目标：建立实验-算法反馈闭环 + 在 pipeline 中嵌入可训练的 AI 模块。

- **目标**：
  1. **反馈闭环**：实验结果自动化反馈到算法层面，消除人力分析瓶颈
  2. **学术定位**：需要可训练模块体现 AI 能力，当前 pipeline 缺乏 AI 闭环特征
- **模型演进路径**（共享 ESM-2 backbone，渐进式升级）：
  - **阶段一：二分类模型**（R2 数据即可启动）
    - 架构：ESM-2 embedding → 线性分类 head → sigmoid
    - 损失函数：binary cross-entropy
    - 目标：验证 ESM backbone 特征对活性预测的有效性
  - **阶段二：排序模型**（R3 数据回收后升级）
    - 架构：pairwise learning-to-rank，输入变体对 (A, B)，学习偏序关系
    - 损失函数：margin ranking loss / contrastive loss
    - 输出：连续打分值（活性变体得分 > 无活性变体）
    - 优势：不受 OD 绝对值影响，只需相对排序信息
  - **迁移成本低**：两阶段共享 ESM-2 backbone，仅替换 head + loss
- **在 pipeline 中的嵌入位置**：Step 3（ESM 评分）和 Step 8（最终筛选）之间，作为实验数据校准的中间层

---

## 六、未探索区域扫描

> 目标：系统性探索当前设计未覆盖的突变空间，发现新的 pH 敏感性增强位点。

- **背景**：R1-R3 设计集中在已知热点区域，框架区和部分 CDR 尚未系统扫描
- **方法**：
  - 对 VH/VL 全残基进行 FoldX alanine scanning，识别能量贡献大的位点
  - 对界面 5Å 内的非 His 残基进行 His 替换扫描，预测新 pH 开关位点
  - 使用 pKa 工具（PROPKA3 + pKAI+）评估候选 His 位点的 pKa 是否落在 6.0-7.4 窗口
  - 结合 ESM-2 注意力图识别序列空间中被忽略的关键位点
- **预期产出**：新候选 pH 开关位点列表，为 R4+ 设计提供计算指导
- **依赖**：pKa 工具就绪后启动（已部署）

---

## 七、时间线

> 以实验轮次为锚点，计算工作与实验制备并行推进。

```
Pipeline 3 Tier 重构完成（当前）
│
├─ 已完成 ───────────────────────────────────────────────
│   ├─ [一] ✅ R3 10 变体设计 + 计算评估
│   ├─ [二] ✅ SimpleFold 3x CDR RMSD 验证 + pipeline 集成
│   ├─ [三] ✅ PROPKA3/pKAI+/dddG_elec 部署 + pipeline 集成
│   └─ Pipeline 3 Tier / 13 Step 重构 + 全流程验证
│
├─ 近期 ─────────────────────────────────────────────────
│   ├─ [四] sdab 组合设计（pH 工具已就绪）
│   ├─ [五-阶段一] 用 R2 数据训练二分类原型
│   └─ R3 湿实验结果回收 → 验证 3 Tier 筛选效果
│
├─ R3 结果回收后 ────────────────────────────────────────
│   ├─ [五-阶段二] 合并 R2+R3 数据，升级为排序模型
│   ├─ 用实验数据校准 Tier 筛选阈值
│   └─ 设计 R4 变体
│
└─ R4+ 持续迭代 ─────────────────────────────────────────
    ├─ 1E62 + sdab 实验数据持续回流，模型持续校准
    └─ [六] 未探索区域扫描
```
