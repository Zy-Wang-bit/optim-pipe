# 3 Tier Pipeline 重构 — 实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将 pipeline 从 4 Phase / 12 Step 重构为 3 Tier / 13 Step 漏斗结构，每个指标有明确的角色定位（门槛/排序/软标记）。

**Architecture:** Tier 1 对 FoldX 结果做自适应门槛筛选（dG + delta），输出 ~500 候选。Tier 2 通过 PyRosetta 线（pKa/dddG_elec）和 SimpleFold 3x 线（CDR RMSD）并行评估后筛选。Tier 3 按 dddG_elec 精排 + 软标记输出最终候选。

**Tech Stack:** Python 3.11, PyRosetta, SimpleFold 3B, PROPKA3, pKAI+, FoldX, BioPython, pandas

**Spec:** `docs/superpowers/specs/2026-04-07-metric-strategy-3tier-design.md`

---

## File Structure

| 文件 | 操作 | 职责 |
|------|------|------|
| `configs/config.yaml` | 修改 | 新增 `tier1` / `tier2` 段，重组 `phase_c` → `tier2` |
| `scripts/tier1_filter.py` | 新建 | Tier 1 自适应门槛筛选 (dG_pH7.4 + delta + ESM 软标记) |
| `scripts/run_simplefold_3x.py` | 新建 | SimpleFold 3x 采样封装 |
| `scripts/tier2_filter.py` | 新建 | Tier 2 门槛筛选 (pKa 相对排序 + CDR RMSD) |
| `scripts/merge_and_rank.py` | 新建 | Tier 3 精排 (dddG_elec 排序 + 软标记) |
| `run_pipeline.sh` | 修改 | 适配 3 Tier 流程 |
| `scripts/build_structures.py` | 修改 | 读取 `tier1_candidates.csv` 替代 `final_top10k.csv` |
| `scripts/run_rmsd.py` | 修改 | 支持 3x 采样中位数 + 异常值剔除 |
| `CLAUDE.md` | 修改 | 更新 pipeline 结构描述 |

保留不动的文件：`scripts/merge_and_select.py`（旧模式向后兼容）、Steps 1-7 所有脚本。

---

## Task 1: 配置文件重组

**Files:**
- Modify: `configs/config.yaml`

- [ ] **Step 1: 新增 tier1 段**

在 `selection` 段之后（第 72 行后）新增：

```yaml
# ── Tier 1: 高通量筛选（Step 8）─────────────────────────────────────────────
tier1:
  enabled: true                          # false 时沿用旧 merge_and_select.py
  output: "results/tier1_candidates.csv"
  adaptive:
    target_range: [300, 500]             # 自适应目标候选数范围
    start:      { dG74_max: -3.0, delta_min: 0.5 }
    relax_step: { dG74_max: 0.5,  delta_min: -0.25 }
    floors:     { dG74_max: 1.0,  delta_min: 0.0 }
  esm:
    flag_below_percentile: 5             # ESM 低于第 5 百分位标记为异常
```

- [ ] **Step 2: 重组 phase_c 为 tier2**

将现有 `phase_c` 段重命名为 `tier2`，拆分两条结构评估线。保留 `phase_c` 作为别名兼容：

```yaml
# ── Tier 2: 结构评估 + 筛选（Steps 9-12）───────────────────────────────────
tier2:
  enabled: false                         # 设为 true 启用 Tier 2-3

  paths:
    wt_pdb: "experiments/1E62/data/pdb/ab_wt.pdb"
    wt_ab_pdb: "experiments/1E62/data/pdb/ab_wt_antibody.pdb"  # 单抗体结构（SimpleFold RMSD 用）
    template_pdb: "experiments/1E62/data/pdb/AF3-1E62-AeS-1.pdb"  # 复合物模板（PyRosetta 用）
    tier2_dir: "tier2"                   # Tier 2 所有输出的根目录

  input:
    csv: "results/tier1_candidates.csv"  # Tier 2 输入（Tier 1 输出）
    mutations_column: "mutations"
    mutations_format: "foldx"

  chains:
    binder: ["A", "B"]
    target: ["C"]
    heavy: "A"
    light: "B"

  # PyRosetta 线
  rosetta:
    repack_shell: 8.0
    minimize: true
    scorefxn: "ref2015"

  # SimpleFold 3x 线
  simplefold:
    model: "simplefold_3B"
    num_steps: 500
    nsample_per_protein: 3
    tau: 0.01
    outlier_global_rmsd: 2.0             # 全局 RMSD > 此值视为 VH-VL 朝向异常

  pka:
    tools: ["propka", "pkai+"]
    his_positions:
      - { chain: "B", resid: 42 }
      - { chain: "B", resid: 43 }
      - { chain: "B", resid: 44 }

  cdr_regions:
    H1: { chain: "A", start: 26, end: 35 }
    H2: { chain: "A", start: 50, end: 65 }
    H3: { chain: "A", start: 95, end: 102 }
    L1: { chain: "B", start: 24, end: 34 }
    L2: { chain: "B", start: 50, end: 56 }
    L3: { chain: "B", start: 89, end: 97 }

  filter:
    pka_relative: true                   # 按 pKa 相对于 WT 排序（非绝对阈值）
    pka_top_n: 200                       # pKa 排序后取 top-N
    cdr_rmsd_max: 0.5                    # CDR RMSD 硬门槛 (Å)
    cdr_rmsd_metric: "h1_rmsd"           # 使用哪个 CDR RMSD（默认 H1）
    output: "results/tier2_candidates.csv"

# ── Tier 3: 精排（Step 13）─────────────────────────────────────────────────
tier3:
  rank_by: "dddG_elec"                  # 排序指标
  rank_ascending: false                  # false = 降序（dddG_elec 越大越好）
  soft_flags:
    esm_flag: true                       # 标记 ESM 异常
    phscore_flag: true                   # 标记 pH-score 异常
    consensus_flag: true                 # 标记 pKa consensus != "agree"
  output: "results/final_candidates.csv"
```

- [ ] **Step 3: 保留旧 selection 段用于向后兼容**

不删除 `selection` 段。当 `tier1.enabled: false` 时，`merge_and_select.py` 沿用旧 `selection` 配置。

- [ ] **Step 4: Commit**

```bash
git add configs/config.yaml
git commit -m "refactor: add tier1/tier2/tier3 config sections for 3-tier pipeline"
```

---

## Task 2: tier1_filter.py — Tier 1 自适应筛选

**Files:**
- Create: `scripts/tier1_filter.py`

- [ ] **Step 1: 实现 tier1_filter.py**

```python
#!/usr/bin/env python3
"""
Step 8 (Tier 1): 自适应门槛筛选

从 FoldX + ESM 结果中按 dG_pH7.4 和 delta 做自适应门槛筛选。
输出 tier1_candidates.csv 供 Tier 2 使用。
"""
import os, sys, glob, yaml
import pandas as pd
import numpy as np


def _load_foldx_all(foldx_dir):
    """加载所有 FoldX 批次的 summary CSV。复用 merge_and_select 的逻辑。"""
    rows = []
    for csvp in glob.glob(os.path.join(foldx_dir, "batches", "*", "batch_*", "foldx_summary.csv")):
        pid = csvp.split(os.sep)[-3]
        bid = csvp.split(os.sep)[-2]
        df = pd.read_csv(csvp)
        if "mpdb" not in df.columns:
            print(f"[WARN] {csvp} 缺少 mpdb 列，跳过")
            continue
        num_cols = ["dG_pH7_4", "dG_pH6_0", "WT_dG_pH7_4", "WT_dG_pH6_0",
                    "ddG_pH7_4", "ddG_pH6_0", "delta", "delta_wt", "delta_delta"]
        for c in num_cols:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
        df["pdb_id"] = pid
        df["batch"] = bid
        rows.append(df)
    if not rows:
        raise RuntimeError("未找到 foldx_summary.csv")
    fx = pd.concat(rows, ignore_index=True)
    # 补算缺失的 delta
    if "delta" not in fx or fx["delta"].isna().all():
        if {"dG_pH6_0", "dG_pH7_4"}.issubset(fx.columns):
            fx["delta"] = fx["dG_pH6_0"] - fx["dG_pH7_4"]
    return fx


def _load_batch_meta(foldx_dir):
    """加载批次元数据（序列、突变、ESM 分数等）。"""
    rows = []
    for bcsv in glob.glob(os.path.join(foldx_dir, "batches", "*", "batch_*", "batch_seqs.csv")):
        pid = bcsv.split(os.sep)[-3]
        bid = bcsv.split(os.sep)[-2]
        df = pd.read_csv(bcsv)
        if "mpdb" not in df.columns:
            if "mpdb_base" in df.columns:
                df["mpdb"] = df["mpdb_base"].astype(str) + ".pdb"
            else:
                continue
        df["pdb_id"] = pid
        df["batch"] = bid
        rows.append(df)
    if not rows:
        return pd.DataFrame(columns=["pdb_id", "batch", "mpdb"])
    return pd.concat(rows, ignore_index=True)


def _merge_foldx_meta(fx, meta):
    """合并 FoldX 结果和元数据。"""
    keep_cols = [c for c in ["sequence", "mutations", "source", "esm_avg_logprob", "hits_hotspots"]
                 if c in meta.columns]
    merged = fx.merge(meta[["pdb_id", "batch", "mpdb"] + keep_cols],
                      on=["pdb_id", "batch", "mpdb"], how="left")
    return merged


def _flag_esm(df, percentile):
    """标记 ESM 异常值（低于指定百分位）。"""
    if "esm_avg_logprob" not in df.columns or df["esm_avg_logprob"].isna().all():
        df["esm_flag"] = False
        return df
    threshold = np.nanpercentile(df["esm_avg_logprob"].astype(float), percentile)
    df["esm_flag"] = df["esm_avg_logprob"].astype(float) < threshold
    return df


def _adaptive_filter(df, t1_cfg):
    """自适应门槛筛选：调整 dG_pH7.4 和 delta 阈值直到候选量落入目标范围。"""
    target_lo, target_hi = t1_cfg["adaptive"]["target_range"]
    start = t1_cfg["adaptive"]["start"]
    step = t1_cfg["adaptive"]["relax_step"]
    floor = t1_cfg["adaptive"]["floors"]

    dG = float(start["dG74_max"])
    dm = float(start["delta_min"])
    dG_step = float(step["dG74_max"])
    dm_step = float(step["delta_min"])
    dG_floor = float(floor["dG74_max"])
    dm_floor = float(floor["delta_min"])

    best = pd.DataFrame()
    while True:
        q = (df["dG_pH7_4"] < dG) & (df["delta"] >= dm)
        pool = df[q]

        # 去重
        dedup_key = "mutations" if "mutations" in pool.columns else "sequence" if "sequence" in pool.columns else None
        if dedup_key:
            pool = pool.drop_duplicates(subset=[dedup_key], keep="first")

        n = len(pool)
        print(f"  dG74 < {dG:.1f}, delta >= {dm:.2f} → {n} 候选")

        if target_lo <= n <= target_hi:
            best = pool
            break
        if n > target_hi:
            # 通过数太多，收紧
            if len(pool) > len(best):
                best = pool.copy()
            dG -= abs(dG_step)
            dm += abs(dm_step)
            if dG < float(start["dG74_max"]) - 5 * abs(dG_step):
                best = pool.head(target_hi)
                break
        else:
            # 通过数太少，放宽
            if len(pool) > len(best):
                best = pool.copy()
            stop_dG = dG >= dG_floor - 1e-9
            stop_dm = dm <= dm_floor + 1e-9
            if stop_dG and stop_dm:
                break
            dG = min(dG_floor, dG + abs(dG_step))
            dm = max(dm_floor, dm + dm_step)

    return best


def main(cfg_path):
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t1_cfg = cfg["tier1"]
    paths = cfg.get("paths", {})
    foldx_dir = paths.get("foldx_dir", cfg.get("foldx_dir", "foldx"))

    print("[Tier 1] 加载 FoldX 结果...")
    fx = _load_foldx_all(foldx_dir)
    print(f"  FoldX 总行数: {len(fx)}")

    meta = _load_batch_meta(foldx_dir)
    merged = _merge_foldx_meta(fx, meta)

    # ESM 异常标记
    esm_pct = t1_cfg.get("esm", {}).get("flag_below_percentile", 5)
    merged = _flag_esm(merged, esm_pct)

    print("[Tier 1] 自适应筛选...")
    filtered = _adaptive_filter(merged, t1_cfg)

    # 按 delta 降序排序
    filtered = filtered.sort_values("delta", ascending=False).reset_index(drop=True)

    out_path = t1_cfg.get("output", "results/tier1_candidates.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    out_cols = [c for c in [
        "sequence", "mutations", "pdb_id", "batch", "mpdb", "source",
        "esm_avg_logprob", "esm_flag", "hits_hotspots",
        "dG_pH7_4", "dG_pH6_0", "delta",
        "WT_dG_pH7_4", "WT_dG_pH6_0",
        "ddG_pH7_4", "ddG_pH6_0",
    ] if c in filtered.columns]
    filtered[out_cols].to_csv(out_path, index=False)
    print(f"[Tier 1] 输出: {out_path} (N={len(filtered)})")


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    main(cfg)
```

- [ ] **Step 2: 验证 — 用 R3 数据手动测试**

```bash
# 确保 FoldX 结果存在（使用已有的 R3 批次数据）
python scripts/tier1_filter.py configs/config.yaml
# 检查输出
head -5 results/tier1_candidates.csv
wc -l results/tier1_candidates.csv
```

预期：输出 CSV 包含 dG_pH7_4 < 阈值且 delta > 0 的候选，有 esm_flag 列。

- [ ] **Step 3: Commit**

```bash
git add scripts/tier1_filter.py
git commit -m "feat: add tier1_filter.py — adaptive threshold screening (Step 8)"
```

---

## Task 3: run_simplefold_3x.py — SimpleFold 3x 采样

**Files:**
- Create: `scripts/run_simplefold_3x.py`

- [ ] **Step 1: 实现 run_simplefold_3x.py**

```python
#!/usr/bin/env python3
"""
Step 10 (Tier 2): SimpleFold 3x 采样

从 tier1_candidates.csv 读取序列，用 SimpleFold 3B 每变体 3 次采样。
输出到 {tier2_dir}/structures/simplefold/

与 build_structures.py 的 SimpleFold 方法类似，但独立运行：
- 不依赖 PyRosetta
- 使用单抗体序列（不含抗原）
- 3x 采样用于后续 RMSD 中位数计算
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


def _extract_antibody_sequence(row, cfg):
    """从候选行提取单抗体序列。

    候选 CSV 的 sequence 列可能是多链格式 (A/B/C)，
    需要去掉抗原链只保留抗体链。
    """
    seq = row.get("sequence", "")
    if not seq or (isinstance(seq, float) and pd.isna(seq)):
        return None

    chains_cfg = cfg["tier2"]["chains"]
    ab_chains = [chains_cfg["heavy"], chains_cfg["light"]]

    # 多链格式: "HEAVY_SEQ/LIGHT_SEQ/ANTIGEN_SEQ"
    parts = str(seq).split("/")
    if len(parts) >= 2:
        # 按配置的链顺序取抗体链
        # 假设 A=Heavy, B=Light, C=Antigen → parts[0]/parts[1]
        return "/".join(parts[:len(ab_chains)])
    return seq


def main():
    parser = argparse.ArgumentParser(description="Tier 2 Step 10: SimpleFold 3x 采样")
    parser.add_argument("config", nargs="?", default="configs/config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    t2 = cfg["tier2"]
    sf_cfg = t2["simplefold"]
    t2_dir = t2["paths"]["tier2_dir"]
    input_csv = t2["input"]["csv"]

    out_dir = os.path.join(t2_dir, "structures", "simplefold")
    fasta_dir = os.path.join(out_dir, "fasta_input")
    os.makedirs(fasta_dir, exist_ok=True)

    df = pd.read_csv(input_csv)
    print(f"[SimpleFold 3x] 读取 {len(df)} 个候选")

    # 生成 FASTA 文件
    n_fasta = 0
    for _, row in df.iterrows():
        if "mpdb" in row and pd.notna(row["mpdb"]):
            vid = Path(str(row["mpdb"])).stem
        else:
            vid = f"var_{_:06d}"

        seq = _extract_antibody_sequence(row, cfg)
        if not seq:
            print(f"  跳过 {vid}: 无序列")
            continue

        fasta_path = os.path.join(fasta_dir, f"{vid}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{vid}\n{seq}\n")
        n_fasta += 1

    print(f"[SimpleFold 3x] 生成 {n_fasta} 个 FASTA 文件")

    cmd = [
        "simplefold",
        "--simplefold_model", sf_cfg["model"],
        "--num_steps", str(sf_cfg["num_steps"]),
        "--tau", str(sf_cfg["tau"]),
        "--nsample_per_protein", str(sf_cfg["nsample_per_protein"]),
        "--plddt",
        "--fasta_path", fasta_dir,
        "--output_dir", out_dir,
        "--output_format", "pdb",
    ]

    print(f"[SimpleFold 3x] 命令: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"[SimpleFold 3x] 完成 → {out_dir}")
    except FileNotFoundError:
        print("[SimpleFold 3x] 错误: simplefold 命令未找到，请确认 simplefold conda 环境已激活")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"[SimpleFold 3x] 错误:\n{e.stderr}")
        sys.exit(1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add scripts/run_simplefold_3x.py
git commit -m "feat: add run_simplefold_3x.py — SimpleFold 3x sampling (Step 10)"
```

---

## Task 4: 修改 run_rmsd.py — 支持 3x 采样中位数

**Files:**
- Modify: `scripts/run_rmsd.py`

- [ ] **Step 1: 阅读现有 run_rmsd.py**

```bash
# 在实现前先阅读现有代码
cat scripts/run_rmsd.py
```

- [ ] **Step 2: 添加 3x 采样中位数逻辑**

现有 `run_rmsd.py` 的 `evaluate_structure_set()` 函数处理单次采样。需要新增 `evaluate_simplefold_3x()` 函数，处理：
1. 每个变体有 3 个 PDB（sample_0, sample_1, sample_2）
2. 剔除全局 RMSD > 阈值的异常样本
3. 对每个 CDR RMSD 取中位数

在 `run_rmsd.py` 文件末尾（`main()` 之前）添加：

```python
def evaluate_simplefold_3x(wt_pdb, struct_dir, cdr_regions, outlier_threshold=2.0):
    """
    评估 SimpleFold 3x 采样的 CDR RMSD。

    每个变体有 3 个样本，剔除 global RMSD > outlier_threshold 的异常样本，
    取各 CDR RMSD 的中位数。

    Returns: [{"variant_id": str, "source": "simplefold_3x",
               "global_rmsd": float, "{cdr}_rmsd": float, "n_valid_samples": int}]
    """
    import re
    from collections import defaultdict

    # 收集所有 PDB 文件，按变体 ID 分组
    all_pdbs = sorted(glob.glob(os.path.join(struct_dir, "*.pdb")))
    variant_pdbs = defaultdict(list)
    for p in all_pdbs:
        stem = Path(p).stem
        # SimpleFold 输出格式: {variant_id}_sample_{N} 或 {variant_id}_{N}
        match = re.match(r"(.+?)(?:_sample)?_(\d+)$", stem)
        if match:
            vid = match.group(1)
            variant_pdbs[vid].append(p)
        else:
            variant_pdbs[stem].append(p)

    print(f"[RMSD 3x] {len(variant_pdbs)} 个变体，共 {len(all_pdbs)} 个 PDB")

    results = []
    for vid, pdbs in variant_pdbs.items():
        # 对每个样本计算 RMSD
        sample_results = []
        for pdb_path in pdbs:
            try:
                r = evaluate_structure_set(wt_pdb, [pdb_path], cdr_regions, "simplefold_3x")
                if r:
                    sample_results.append(r[0])
            except Exception as e:
                print(f"  {vid} 样本 {pdb_path} 失败: {e}")

        if not sample_results:
            continue

        # 剔除 global RMSD 异常值
        valid = [r for r in sample_results if r["global_rmsd"] <= outlier_threshold]
        if not valid:
            print(f"  {vid}: 所有 {len(sample_results)} 样本均为异常值 (global RMSD > {outlier_threshold}Å)")
            continue

        # 取中位数
        median_result = {"variant_id": vid, "source": "simplefold_3x",
                         "n_valid_samples": len(valid)}
        rmsd_keys = [k for k in valid[0] if k.endswith("_rmsd")]
        for k in rmsd_keys:
            vals = [r[k] for r in valid if k in r]
            median_result[k] = float(np.median(vals))

        results.append(median_result)

    return results
```

- [ ] **Step 3: 修改 main() 支持 3x 模式**

在 `main()` 中新增判断：当 SimpleFold 结构目录存在且配置了 `tier2.simplefold` 时，使用 `evaluate_simplefold_3x()`：

```python
# 在 main() 的 source 循环之后添加
# SimpleFold 3x 采样模式
sf_dir = os.path.join(t2_dir, "structures", "simplefold")
if os.path.isdir(sf_dir):
    outlier_th = cfg.get("tier2", cfg.get("phase_c", {})).get(
        "simplefold", {}).get("outlier_global_rmsd", 2.0)
    sf3x_results = evaluate_simplefold_3x(wt_pdb, sf_dir, cdr_regions, outlier_th)
    all_results.extend(sf3x_results)
    print(f"[RMSD] SimpleFold 3x: {len(sf3x_results)} 变体 (中位数)")
```

- [ ] **Step 4: Commit**

```bash
git add scripts/run_rmsd.py
git commit -m "feat: add SimpleFold 3x median RMSD with outlier removal"
```

---

## Task 5: tier2_filter.py — Tier 2 pKa 相对排序 + RMSD 筛选

**Files:**
- Create: `scripts/tier2_filter.py`

- [ ] **Step 1: 实现 tier2_filter.py**

```python
#!/usr/bin/env python3
"""
Step 12 (Tier 2): 门槛筛选

合并 PyRosetta 线 (pKa, dddG_elec, pH-score) 和 SimpleFold 3x 线 (CDR RMSD) 的结果。
按 pKa 相对于 WT 排序，取 top-N；按 CDR RMSD 硬门槛过滤。
输出 tier2_candidates.csv 供 Tier 3 使用。
"""
import os, sys, yaml
import pandas as pd
import numpy as np
from pathlib import Path


def _load_tier1(cfg):
    """加载 Tier 1 候选列表。"""
    t1_out = cfg["tier1"].get("output", "results/tier1_candidates.csv")
    df = pd.read_csv(t1_out)
    # 添加 variant_id
    if "mpdb" in df.columns:
        df["variant_id"] = df["mpdb"].apply(
            lambda x: Path(str(x)).stem if pd.notna(x) else None)
    return df


def _attach_pka(df, t2_dir):
    """左连接 pKa 结果。"""
    pka_path = os.path.join(t2_dir, "pka", "pka_summary.csv")
    if not os.path.exists(pka_path):
        print(f"[Tier 2] pKa 文件不存在: {pka_path}")
        return df
    pka = pd.read_csv(pka_path)
    if "variant_name" in pka.columns and "variant_id" not in pka.columns:
        pka = pka.rename(columns={"variant_name": "variant_id"})
    pka_cols = [c for c in ["variant_id", "avg_shift_propka", "avg_shift_pkai",
                            "overall_consensus", "pKa_propka", "pKa_pkai"] if c in pka.columns]
    df = df.merge(pka[pka_cols], on="variant_id", how="left")
    print(f"[Tier 2] pKa: {pka_path} ({len(pka)} 行)")
    return df


def _attach_rosetta(df, t2_dir):
    """左连接 Rosetta 结果 (dddG_elec + pH-score)。"""
    for filename, label in [("dddg_elec.csv", "dddG_elec"), ("ph_scores.csv", "pH-score")]:
        path = os.path.join(t2_dir, "rosetta", filename)
        if not os.path.exists(path):
            continue
        sub = pd.read_csv(path)
        if "variant_name" in sub.columns and "variant_id" not in sub.columns:
            sub = sub.rename(columns={"variant_name": "variant_id"})
        cols = [c for c in sub.columns if c == "variant_id" or c in ["dddG_elec", "ph_score",
                "ddG_elec_pH7", "ddG_elec_pH5"]]
        df = df.merge(sub[cols], on="variant_id", how="left")
        print(f"[Tier 2] {label}: {path} ({len(sub)} 行)")
    return df


def _attach_rmsd(df, t2_dir):
    """左连接 RMSD 结果（优先 SimpleFold 3x）。"""
    rmsd_path = os.path.join(t2_dir, "rmsd", "rmsd_summary.csv")
    if not os.path.exists(rmsd_path):
        print(f"[Tier 2] RMSD 文件不存在: {rmsd_path}")
        return df
    rmsd = pd.read_csv(rmsd_path)
    if "variant_name" in rmsd.columns and "variant_id" not in rmsd.columns:
        rmsd = rmsd.rename(columns={"variant_name": "variant_id"})
    # 优先 simplefold_3x 来源
    if "source" in rmsd.columns:
        sf3x = rmsd[rmsd["source"] == "simplefold_3x"]
        if len(sf3x) > 0:
            rmsd = sf3x
        else:
            # 回退到 primary
            primary = rmsd[rmsd["source"] == rmsd["source"].iloc[0]]
            rmsd = primary
    rmsd_cols = ["variant_id"] + [c for c in rmsd.columns if c.endswith("_rmsd")]
    df = df.merge(rmsd[rmsd_cols], on="variant_id", how="left")
    print(f"[Tier 2] RMSD: {rmsd_path} ({len(rmsd)} 行)")
    return df


def _filter_pka_relative(df, wt_pdb, t2_dir, filter_cfg):
    """按 pKa 相对于 WT 排序，取 top-N。"""
    # 计算 WT pKa 基线
    wt_pka_path = os.path.join(t2_dir, "pka", "wt_pka.csv")
    wt_pka_avg = None
    if os.path.exists(wt_pka_path):
        wt_pka = pd.read_csv(wt_pka_path)
        if "pKa_propka" in wt_pka.columns:
            wt_pka_avg = wt_pka["pKa_propka"].mean()

    # 排序指标: pKa_propka 的绝对值或相对于 WT 的提升
    sort_col = None
    if "pKa_propka" in df.columns and df["pKa_propka"].notna().any():
        if wt_pka_avg is not None:
            df["pka_uplift"] = df["pKa_propka"] - wt_pka_avg
            sort_col = "pka_uplift"
        else:
            sort_col = "pKa_propka"
    elif "avg_shift_propka" in df.columns and df["avg_shift_propka"].notna().any():
        sort_col = "avg_shift_propka"

    if sort_col is None:
        print("[Tier 2] 无可用 pKa 数据，跳过 pKa 排序")
        return df

    top_n = filter_cfg.get("pka_top_n", 200)
    # 排除比 WT 差的（uplift <= 0）
    if "pka_uplift" in df.columns:
        df = df[df["pka_uplift"] > 0]
        print(f"[Tier 2] pKa 相对排序: 排除 uplift <= 0 后剩 {len(df)} 候选")

    df = df.sort_values(sort_col, ascending=False).head(top_n)
    print(f"[Tier 2] pKa top-{top_n}: {len(df)} 候选 (排序列: {sort_col})")
    return df


def _filter_rmsd(df, filter_cfg):
    """CDR RMSD 硬门槛。"""
    metric = filter_cfg.get("cdr_rmsd_metric", "h1_rmsd")
    threshold = float(filter_cfg.get("cdr_rmsd_max", 0.5))

    if metric not in df.columns:
        print(f"[Tier 2] RMSD 列 {metric} 不存在，跳过 RMSD 筛选")
        return df

    before = len(df)
    df = df[df[metric] <= threshold]
    print(f"[Tier 2] RMSD 筛选 ({metric} <= {threshold}Å): {before} → {len(df)}")
    return df


def main(cfg_path):
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t2 = cfg["tier2"]
    t2_dir = t2["paths"]["tier2_dir"]
    filter_cfg = t2["filter"]

    print("[Tier 2] 加载 Tier 1 候选...")
    df = _load_tier1(cfg)
    print(f"  Tier 1 候选: {len(df)}")

    # 合并 Phase C 评估结果
    df = _attach_pka(df, t2_dir)
    df = _attach_rosetta(df, t2_dir)
    df = _attach_rmsd(df, t2_dir)

    # 筛选
    wt_pdb = t2["paths"].get("wt_pdb", "")
    df = _filter_pka_relative(df, wt_pdb, t2_dir, filter_cfg)
    df = _filter_rmsd(df, filter_cfg)

    # 输出
    out_path = filter_cfg.get("output", "results/tier2_candidates.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"[Tier 2] 输出: {out_path} (N={len(df)})")


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    main(cfg)
```

- [ ] **Step 2: Commit**

```bash
git add scripts/tier2_filter.py
git commit -m "feat: add tier2_filter.py — pKa relative ranking + CDR RMSD gate (Step 12)"
```

---

## Task 6: merge_and_rank.py — Tier 3 精排

**Files:**
- Create: `scripts/merge_and_rank.py`

- [ ] **Step 1: 实现 merge_and_rank.py**

```python
#!/usr/bin/env python3
"""
Step 13 (Tier 3): 精排 + 最终输出

从 tier2_candidates.csv 读取候选，按 dddG_elec 排序，添加软标记。
输出 final_candidates.csv。
"""
import os, sys, json, yaml
import pandas as pd
import numpy as np


def _add_soft_flags(df, t3_cfg):
    """添加软标记列。"""
    flags = t3_cfg.get("soft_flags", {})

    # ESM 异常标记（从 Tier 1 继承）
    if flags.get("esm_flag") and "esm_flag" in df.columns:
        pass  # 已有
    elif flags.get("esm_flag"):
        df["esm_flag"] = False

    # pH-score 异常标记（z-score > 2）
    if flags.get("phscore_flag") and "ph_score" in df.columns:
        ph = df["ph_score"].astype(float)
        mean, std = ph.mean(), ph.std()
        if std > 0:
            df["phscore_flag"] = (ph - mean).abs() > 2 * std
        else:
            df["phscore_flag"] = False
    else:
        df["phscore_flag"] = False

    # pKa consensus 标记
    if flags.get("consensus_flag") and "overall_consensus" in df.columns:
        df["consensus_flag"] = df["overall_consensus"].apply(
            lambda x: x not in ("agree", "neutral") if pd.notna(x) else True)
    else:
        df["consensus_flag"] = False

    return df


def main(cfg_path):
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t2_filter = cfg["tier2"]["filter"]
    t3_cfg = cfg["tier3"]

    input_csv = t2_filter.get("output", "results/tier2_candidates.csv")
    df = pd.read_csv(input_csv)
    print(f"[Tier 3] 输入: {input_csv} ({len(df)} 候选)")

    # 排序
    rank_col = t3_cfg.get("rank_by", "dddG_elec")
    rank_asc = t3_cfg.get("rank_ascending", False)

    if rank_col in df.columns and df[rank_col].notna().any():
        df = df.sort_values(rank_col, ascending=rank_asc).reset_index(drop=True)
        print(f"[Tier 3] 按 {rank_col} {'升序' if rank_asc else '降序'}排序")
    else:
        print(f"[Tier 3] 排序列 {rank_col} 不存在或全为空，保持原顺序")

    # 软标记
    df = _add_soft_flags(df, t3_cfg)

    # 输出
    out_path = t3_cfg.get("output", "results/final_candidates.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # 选择输出列
    out_cols = [c for c in [
        "variant_id", "sequence", "mutations", "pdb_id", "mpdb",
        # Tier 1 指标
        "dG_pH7_4", "ddG_pH7_4", "delta", "esm_avg_logprob",
        # Tier 2 指标
        "pKa_propka", "pKa_pkai", "avg_shift_propka", "avg_shift_pkai",
        "overall_consensus", "pka_uplift",
        "dddG_elec", "ddG_elec_pH7", "ddG_elec_pH5", "ph_score",
        "global_rmsd", "h1_rmsd", "h2_rmsd", "h3_rmsd",
        "l1_rmsd", "l2_rmsd", "l3_rmsd",
        # 标记
        "esm_flag", "phscore_flag", "consensus_flag",
    ] if c in df.columns]
    df[out_cols].to_csv(out_path, index=False)

    # 审计日志
    audit = {
        "tier3_input": len(pd.read_csv(input_csv)),
        "tier3_output": len(df),
        "rank_by": rank_col,
        "flags": {
            "esm_flagged": int(df["esm_flag"].sum()) if "esm_flag" in df.columns else 0,
            "phscore_flagged": int(df["phscore_flag"].sum()) if "phscore_flag" in df.columns else 0,
            "consensus_flagged": int(df["consensus_flag"].sum()) if "consensus_flag" in df.columns else 0,
        }
    }
    audit_dir = os.path.join(os.path.dirname(out_path), "audit")
    os.makedirs(audit_dir, exist_ok=True)
    with open(os.path.join(audit_dir, "tier3_manifest.json"), "w") as f:
        json.dump(audit, f, indent=2)

    print(f"[Tier 3] 输出: {out_path} (N={len(df)})")
    n_flagged = sum(1 for _, r in df.iterrows()
                    if r.get("esm_flag") or r.get("phscore_flag") or r.get("consensus_flag"))
    print(f"[Tier 3] 软标记: {n_flagged} 候选至少有 1 个 flag")


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    main(cfg)
```

- [ ] **Step 2: Commit**

```bash
git add scripts/merge_and_rank.py
git commit -m "feat: add merge_and_rank.py — dddG_elec ranking + soft flags (Step 13)"
```

---

## Task 7: Phase C 脚本适配 tier2 配置

**Files:**
- Modify: `scripts/build_structures.py`
- Modify: `scripts/run_pka.py`
- Modify: `scripts/run_rosetta_eval.py`
- Modify: `scripts/run_rmsd.py`

这四个脚本当前从 `cfg["phase_c"]` 读取配置。需要统一添加 `tier2` 优先、`phase_c` 回退的兼容逻辑。

- [ ] **Step 1: 添加配置兼容辅助函数**

在每个脚本的 `main()` 顶部，将配置读取改为：

```python
pc = cfg.get("tier2", cfg.get("phase_c", {}))
```

需要修改的具体位置：

- `build_structures.py` — `load_candidates()` 第 87 行和 `main()` 第 294 行
- `run_pka.py` — `main()` 中所有 `cfg["phase_c"]` 引用
- `run_rosetta_eval.py` — `main()` 中所有 `cfg["phase_c"]` 引用
- `run_rmsd.py` — `main()` 中所有 `cfg["phase_c"]` 引用

同时注意子键映射：
- 旧: `phase_c.paths.phase_c_dir` → 新: `tier2.paths.tier2_dir`
- 旧: `phase_c.structure_generation.primary` → 新: 直接由 `tier2.rosetta` 存在与否判断
- 旧: `phase_c.structure_generation.rosetta.repack_shell` → 新: `tier2.rosetta.repack_shell`

对于路径，添加回退：
```python
t2_dir = pc.get("paths", {}).get("tier2_dir", pc.get("paths", {}).get("phase_c_dir", "tier2"))
```

- [ ] **Step 2: Commit**

```bash
git add scripts/build_structures.py scripts/run_pka.py scripts/run_rosetta_eval.py scripts/run_rmsd.py
git commit -m "refactor: phase_c scripts support tier2 config with phase_c fallback"
```

---

## Task 8: 更新 run_pipeline.sh

**Files:**
- Modify: `run_pipeline.sh`

- [ ] **Step 1: 重写 pipeline 脚本为 3 Tier 结构**

```bash
#!/usr/bin/env bash
# 一键顺序跑：3 Tier / 13 Step 漏斗 pipeline
# 用法：
#   bash run_pipeline.sh [CONFIG]             # 正常跑
#   CLEAN=1 bash run_pipeline.sh [CONFIG]     # 先清理旧产物再跑
#
# Tier 模式判断：
#   tier1.enabled: true  → 使用 3 Tier 流程
#   tier1.enabled: false → 沿用旧 4 Phase 流程 (merge_and_select.py)

set -euo pipefail

CFG="${1:-configs/config.yaml}"
TS="$(date +%Y%m%d_%H%M%S)"
LOGDIR="logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/pipeline_${TS}.log"

MPNN_ARGS=(
  --num-per-pdb 20000
  --shards 5
  --temps 0.10,0.15,0.20,0.25,0.30
  --seed 4242
  --build-his-seeds
)

run(){ echo -e "\n[RUN] $*" | tee -a "$LOG"; "$@" 2>&1 | tee -a "$LOG"; }
note(){ echo -e "\n[NOTE] $*" | tee -a "$LOG"; }
ok(){ echo -e "[OK] $*" | tee -a "$LOG"; }

# 0) 可选清理
if [[ "${CLEAN:-0}" == "1" ]]; then
  note "清理旧产物..."
  rm -rf mpnn_outputs/* his_seeds/* esm_scores/* \
         foldx/repaired/* foldx/batches/* \
         results/screening/* results/final_top10k.csv results/merged_all.csv \
         results/tier1_candidates.csv results/tier2_candidates.csv results/final_candidates.csv \
         tier2/structures/* tier2/pka/* tier2/rosetta/* tier2/rmsd/* \
         phase_c/structures/* phase_c/pka/* phase_c/rosetta/* phase_c/rmsd/* \
         logs/processed_batches.json \
         configs/mpnn/chain_id.jsonl configs/mpnn/fixed_positions.jsonl || true
  ok "清理完成。"
fi

# 基础检查
[[ -f "$CFG" ]] || { echo "[ERR] 找不到配置：$CFG" | tee -a "$LOG"; exit 1; }
mkdir -p results/screening results/audit foldx/repaired foldx/batches \
         esm_scores mpnn_outputs his_seeds configs/mpnn \
         tier2/structures tier2/pka tier2/rosetta tier2/rmsd

# 读取模式
TIER1_ENABLED=$(python3 -c "
import yaml; c=yaml.safe_load(open('${CFG}'))
print(str(c.get('tier1',{}).get('enabled',False)).lower())
")
TIER2_ENABLED=$(python3 -c "
import yaml; c=yaml.safe_load(open('${CFG}'))
print(str(c.get('tier2',{}).get('enabled',False)).lower())
")

# ═══════════════════════════════════════════════════════════════════════════
# Tier 1: 生成 + 高通量筛选
# ═══════════════════════════════════════════════════════════════════════════

note "Step 1: 接口热点扫描"
run python scripts/scan_interface.py "$CFG"
ok "his_hotspots.csv 就绪。"

note "Step 2: MPNN 设计"
run python scripts/run_mpnn_design.py "$CFG" "${MPNN_ARGS[@]}"
ok "MPNN 输出就绪。"

note "Step 3: ESM 评分"
run python scripts/run_esm_chunk.py "$CFG"
ok "for_foldx.csv 就绪。"

note "Step 4: FoldX RepairPDB"
run python scripts/repair_pdbs.py "$CFG"
ok "Repaired PDB 就绪。"

note "Step 5: 预计算 WT 相互作用能"
run python scripts/precompute_wt_ac.py "$CFG"
ok "WT_ac.csv 就绪。"

note "Step 6: 生成 FoldX 批次"
run python scripts/make_mutlist_chunk.py "$CFG"
ok "FoldX 批次就绪。"

note "Step 7: FoldX 评估"
run python scripts/run_foldx_batch.py "$CFG"
ok "foldx_summary.csv 就绪。"

# FoldX 中间文件清理（保留 summary CSV）
note "清理 FoldX 中间文件..."
find foldx/batches -name "*.pdb" -not -name "*_Repair.pdb" -delete 2>/dev/null || true
find foldx/batches -name "*.fxout" -delete 2>/dev/null || true
ok "FoldX 中间文件已清理。"

if [[ "$TIER1_ENABLED" == "true" ]]; then
    # ── 3 Tier 模式 ─────────────────────────────────────────────────────
    note "Step 8: Tier 1 自适应筛选"
    run python scripts/tier1_filter.py "$CFG"
    ok "tier1_candidates.csv 就绪。"

    if [[ "$TIER2_ENABLED" == "true" ]]; then
        # ═══════════════════════════════════════════════════════════════
        # Tier 2: 结构评估 + 筛选
        # ═══════════════════════════════════════════════════════════════

        # PyRosetta 线
        note "Step 9a: PyRosetta 突变体建模"
        run conda run -n pyrosetta python scripts/build_structures.py "$CFG"
        ok "突变体结构就绪。"

        note "Step 9b: pKa 预测"
        run python scripts/run_pka.py "$CFG"
        ok "pKa 预测就绪。"

        note "Step 9c: Rosetta 评分 (dddG_elec + pH-score)"
        run conda run -n pyrosetta python scripts/run_rosetta_eval.py "$CFG"
        ok "Rosetta 评分就绪。"

        # SimpleFold 3x 线
        note "Step 10: SimpleFold 3x 采样"
        run conda run -n simplefold python scripts/run_simplefold_3x.py "$CFG"
        ok "SimpleFold 3x 就绪。"

        note "Step 11: CDR RMSD"
        run python scripts/run_rmsd.py "$CFG"
        ok "CDR RMSD 就绪。"

        # Tier 2 筛选
        note "Step 12: Tier 2 筛选 (pKa + RMSD)"
        run python scripts/tier2_filter.py "$CFG"
        ok "tier2_candidates.csv 就绪。"

        # ═══════════════════════════════════════════════════════════════
        # Tier 3: 精排
        # ═══════════════════════════════════════════════════════════════

        note "Step 13: Tier 3 精排"
        run python scripts/merge_and_rank.py "$CFG"
        ok "final_candidates.csv 就绪。"
    else
        note "Tier 2 已禁用，跳过 Steps 9-13。"
        note "Tier 1 结果: results/tier1_candidates.csv"
    fi
else
    # ── 旧模式 (4 Phase) ─────────────────────────────────────────────
    PHASE_C=$(python3 -c "
    import yaml; c=yaml.safe_load(open('${CFG}'))
    print(str(c.get('phase_c',{}).get('enabled',False)).lower())
    ")

    if [[ "$PHASE_C" == "true" ]]; then
        note "Step 7.5: 预合并（生成 Phase C 候选列表）"
        run python scripts/merge_and_select.py "$CFG"
        ok "Phase C 候选列表就绪。"

        note "Step 8: Phase C - 结构生成"
        run conda run -n pyrosetta python scripts/build_structures.py "$CFG"
        ok "突变体结构就绪。"

        note "Step 9: Phase C - pKa 预测"
        run python scripts/run_pka.py "$CFG"
        ok "pKa 预测就绪。"

        note "Step 10: Phase C - Rosetta 评分"
        run conda run -n pyrosetta python scripts/run_rosetta_eval.py "$CFG"
        ok "Rosetta 评分就绪。"

        note "Step 11: Phase C - CDR RMSD"
        run python scripts/run_rmsd.py "$CFG"
        ok "CDR RMSD 就绪。"
    else
        note "Phase C 已禁用，跳过 Steps 8-11。"
    fi

    note "Step 12: 合并与筛选"
    run python scripts/merge_and_select.py "$CFG"
    ok "最终结果 -> results/final_top10k.csv"
fi

note "全部完成 详细日志：$LOG"
```

- [ ] **Step 2: Commit**

```bash
git add run_pipeline.sh
git commit -m "refactor: run_pipeline.sh supports 3-tier mode with tier1/tier2 fallback"
```

---

## Task 9: 更新 CLAUDE.md + 文档

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: 更新 Pipeline 结构描述**

将 CLAUDE.md 中的 "Pipeline 结构 (4 Phase / 12 Step)" 部分替换为：

```markdown
## Pipeline 结构

### 3 Tier / 13 Step 模式（`tier1.enabled: true`）

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
    11. run_rmsd.py            CDR RMSD（3x 中位数）
  12. tier2_filter.py          pKa 相对排序 + RMSD 门槛 → tier2_candidates.csv

Tier 3: 精排（Step 13）
  13. merge_and_rank.py        dddG_elec 排序 + 软标记 → final_candidates.csv
```

### 旧模式 (4 Phase / 12 Step，`tier1.enabled: false`)

行为与原 8 步 pipeline 一致，Phase C 通过 `phase_c.enabled` 控制。
```

- [ ] **Step 2: 更新目录结构表**

将 `phase_c/` 更新为同时包含 `tier2/`：

```markdown
| `tier2/` | Tier 2 运行时输出（structures/、pka/、rosetta/、rmsd/） |
| `phase_c/` | 旧模式 Phase C 运行时输出 |
```

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md for 3-tier pipeline structure"
```

---

## Task 10: 更新 .gitignore + .tasks

**Files:**
- Modify: `.gitignore`
- Modify: `.tasks/active/pipeline-update/progress.md`

- [ ] **Step 1: 添加 tier2/ 到 .gitignore**

```
# Tier 2 运行时输出
tier2/
```

- [ ] **Step 2: 更新 progress.md**

在 progress.md 末尾添加 3 Tier 重构的进度记录。

- [ ] **Step 3: Commit**

```bash
git add .gitignore .tasks/active/pipeline-update/progress.md
git commit -m "chore: add tier2/ to .gitignore, update task progress"
```

---

## Task 11: 端到端验证

- [ ] **Step 1: 验证旧模式向后兼容**

确保 `tier1.enabled: false` 时行为不变：

```bash
python scripts/merge_and_select.py configs/config.yaml
# 应输出 results/final_top10k.csv，行为与重构前一致
```

- [ ] **Step 2: 验证 Tier 1 筛选**

```bash
# 临时将 tier1.enabled 改为 true
python scripts/tier1_filter.py configs/config.yaml
wc -l results/tier1_candidates.csv
head -3 results/tier1_candidates.csv
# 检查: 有 esm_flag 列，delta > 0，dG_pH7.4 在阈值内
```

- [ ] **Step 3: 验证完整 3 Tier 流程（如有 R3 数据）**

```bash
# 需要 FoldX 结果 + PyRosetta + SimpleFold 环境
# 完整测试需在有计算资源的环境中执行
bash run_pipeline.sh configs/config.yaml
# 检查输出文件:
ls -la results/tier1_candidates.csv results/tier2_candidates.csv results/final_candidates.csv
```

- [ ] **Step 4: 最终 commit**

```bash
git add -A
git commit -m "feat: complete 3-tier pipeline restructure (13 steps)"
```
