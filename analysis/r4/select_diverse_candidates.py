#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Track B: diverse exploration selector for sdab R4 (50 candidates).

Picks 50 candidates from the 9650-row filtered pool such that:
- All pairwise Hamming distances (inside Track B + between Track B and Track A) ≥ 2
- Both model score and Hamming diversity are traded off via MMR
- 5 hypothesis-driven probes are forced in to cover specific design questions
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

PRED_DIR = REPO_ROOT / "experiments/sdab/R4/predictions"
MODEL_DIR = REPO_ROOT / "experiments/sdab/R4/model"
ALL_CSV = PRED_DIR / "candidates_all.csv"
TOP50_CSV = PRED_DIR / "candidates_top50.csv"
FEATURES_JSON = MODEL_DIR / "features.json"

OUT_B = PRED_DIR / "candidates_exploration_top50.csv"
OUT_COMBINED = PRED_DIR / "candidates_combined_100.csv"
OUT_REPORT = PRED_DIR / "exploration_report.md"

SCORE_FLOOR = 0.18         # log(1.2) ≈ 0.18; fallback to 0 if pool too small
MMR_K = 45                 # MMR picks 45; hypothesis picks up to 5 → total 50
MMR_LAMBDA = 0.4           # λ·score + (1-λ)·dist, slightly favor diversity
TRACK_B_TARGET = 50
MIN_HAMMING_TO_A = 3       # stricter: rules out "Track A + 1 mutation" super/sub-set variants
MIN_HAMMING_TO_B = 2       # intra-B diversity


def load_features():
    with open(FEATURES_JSON) as f:
        feats = json.load(f)
    singles = feats["singles"]                    # 31 names in fixed order
    pos_to_idx = {name: i for i, name in enumerate(singles)}
    return singles, pos_to_idx


def muts_to_vec(mut_str, pos_to_idx, n=31):
    """'HG26H;HD110H' -> length-31 binary vector."""
    v = np.zeros(n, dtype=np.int8)
    for m in mut_str.split(";"):
        if m in pos_to_idx:
            v[pos_to_idx[m]] = 1
    return v


def build_vecs(df, pos_to_idx):
    return np.stack([muts_to_vec(s, pos_to_idx) for s in df["mutations_unified"]])


def hamming_matrix(A, B):
    """Pairwise Hamming between rows of A (nA×31) and B (nB×31) → nA×nB."""
    return (A[:, None, :] != B[None, :, :]).sum(axis=2)


def mmr_select(pool_df, pool_vecs, A_vecs, k, lam):
    """Greedy MMR: pick k candidates from pool that maximize λ·score + (1-λ)·Hamming.

    Enforces: Hamming ≥ MIN_HAMMING to Track A and to already-selected.
    """
    scores = pool_df["score"].to_numpy()
    if len(scores) > 1:
        s_min, s_max = scores.min(), scores.max()
        score_norm = (scores - s_min) / (s_max - s_min + 1e-9)
    else:
        score_norm = np.ones_like(scores)

    # Precompute min Hamming to Track A for each pool member
    hdist_to_A = hamming_matrix(pool_vecs, A_vecs).min(axis=1)  # shape (n_pool,)

    # Only candidates far enough from Track A are eligible
    eligible = hdist_to_A >= MIN_HAMMING_TO_A
    eligible_idx = np.where(eligible)[0]
    if len(eligible_idx) == 0:
        return [], pool_vecs[:0], "no eligible candidates (all within Hamming<2 of Track A)"

    # Track min Hamming to S (the selected set); initialize as large
    min_hd_to_S = np.full(len(pool_df), 99, dtype=np.int32)

    selected = []
    # Seed with the highest-scoring eligible candidate
    seed = eligible_idx[np.argmax(scores[eligible_idx])]
    selected.append(seed)
    # update min_hd_to_S
    seed_dists = (pool_vecs != pool_vecs[seed]).sum(axis=1)
    min_hd_to_S = np.minimum(min_hd_to_S, seed_dists)
    min_hd_to_S[seed] = 0  # don't pick self

    while len(selected) < k:
        # Eligible: not picked, Hamming ≥ MIN to both Track A and S
        mask = eligible & (min_hd_to_S >= MIN_HAMMING_TO_B)
        mask[selected] = False
        if not mask.any():
            break
        # MMR objective: λ·score_norm + (1-λ)·min(hdist_to_A, min_hd_to_S)/31
        combined_dist = np.minimum(hdist_to_A, min_hd_to_S) / 31.0
        obj = lam * score_norm + (1 - lam) * combined_dist
        obj[~mask] = -np.inf
        pick = int(np.argmax(obj))
        selected.append(pick)
        new_dists = (pool_vecs != pool_vecs[pick]).sum(axis=1)
        min_hd_to_S = np.minimum(min_hd_to_S, new_dists)
        min_hd_to_S[pick] = 0

    return selected, pool_vecs[selected], None


def find_hypothesis_probes(all_df, all_vecs, A_vecs, pos_to_idx,
                            already_picked_idx, already_picked_vecs):
    """Return list of (all_df_idx, hypothesis_label) for H1-H5.

    Each probe must satisfy Hamming ≥ MIN_HAMMING to Track A AND to already-picked MMR set.
    """
    already = set(already_picked_idx)

    # Helper: min Hamming to Track A and to already-picked (intra-B)
    hdist_A = hamming_matrix(all_vecs, A_vecs).min(axis=1)
    if len(already_picked_vecs) > 0:
        hdist_B = hamming_matrix(all_vecs, already_picked_vecs).min(axis=1)
    else:
        hdist_B = np.full(len(all_df), 99)

    # Indices
    Q100 = pos_to_idx["HQ100H"]
    V105 = pos_to_idx["HV105H"]
    E108 = pos_to_idx["HE108H"]
    Y111 = pos_to_idx["HY111H"]
    D110 = pos_to_idx["HD110H"]
    G55 = pos_to_idx["HG55H"]
    I51 = pos_to_idx["HI51H"]

    cdr1_idx = set(range(0, 8))
    cdr2_idx = set(range(8, 16))
    cdr3_idx = set(range(16, 31))

    def active_positions(vec):
        return {i for i, v in enumerate(vec) if v}

    def qualifies(idx_into_all, predicate):
        if idx_into_all in already:
            return False
        if hdist_A[idx_into_all] < MIN_HAMMING_TO_A:
            return False
        if hdist_B[idx_into_all] < MIN_HAMMING_TO_B:
            return False
        pos = active_positions(all_vecs[idx_into_all])
        return predicate(pos)

    def pick_top(predicate):
        order = np.argsort(-all_df["score"].to_numpy())  # high to low
        for idx in order:
            idx = int(idx)
            if qualifies(idx, predicate):
                return idx
        return None

    probes = []

    # H1: D110H present; V105H, E108H, Y111H all absent
    h1 = pick_top(lambda p: D110 in p and not (p & {V105, E108, Y111}))
    if h1 is not None:
        probes.append((h1, "H1_D110H_no_collapse_partners"))

    # H2: Exactly {Q100H, V105H} 2-mer (is it in pool? already in Track A)
    h2 = pick_top(lambda p: p == {Q100, V105})
    if h2 is not None:
        probes.append((h2, "H2_Q100H_V105H_2mer_alone"))

    # H3: ≥3 CDR1 positions, 0 CDR2+CDR3
    h3 = pick_top(lambda p: len(p & cdr1_idx) >= 3 and not (p & (cdr2_idx | cdr3_idx)))
    if h3 is not None:
        probes.append((h3, "H3_CDR1_triple"))

    # H4: exactly one CDR1 + one CDR2 + one CDR3 (3-mer)
    h4 = pick_top(lambda p: len(p) == 3
                             and len(p & cdr1_idx) == 1
                             and len(p & cdr2_idx) == 1
                             and len(p & cdr3_idx) == 1)
    if h4 is not None:
        probes.append((h4, "H4_cross_CDR_3mer"))

    # H5: contains V105H or I51H as pH-selective single, but does NOT contain G55H
    h5 = pick_top(lambda p: (V105 in p or I51 in p) and G55 not in p
                             and Q100 not in p)  # also not the Q100H+V105H core
    if h5 is not None:
        probes.append((h5, "H5_pH_selective_not_G55H"))

    return probes


def main():
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    singles, pos_to_idx = load_features()

    all_df = pd.read_csv(ALL_CSV).reset_index(drop=True)
    top_df = pd.read_csv(TOP50_CSV).reset_index(drop=True)
    print(f"Loaded {len(all_df)} all candidates, {len(top_df)} Track A candidates.")

    all_vecs = build_vecs(all_df, pos_to_idx)
    A_vecs = build_vecs(top_df, pos_to_idx)

    # Track A mutations as set-of-strings for "is in Track A" check
    trackA_muts = set(top_df["mutations_unified"])

    # Score floor filter
    def build_pool(floor):
        mask = all_df["score"] > floor
        return all_df[mask].reset_index(drop=True), np.where(mask.to_numpy())[0]

    pool_df, pool_all_idx = build_pool(SCORE_FLOOR)
    if len(pool_df) < MMR_K * 2:
        print(f"Pool too small at floor {SCORE_FLOOR:.2f} ({len(pool_df)} rows); relaxing to 0.")
        pool_df, pool_all_idx = build_pool(0.0)
    print(f"MMR pool: {len(pool_df)} candidates (score > {pool_df['score'].min():.3f})")

    pool_vecs = all_vecs[pool_all_idx]

    # MMR selection
    mmr_picks_local, mmr_vecs, warn = mmr_select(pool_df, pool_vecs, A_vecs, MMR_K, MMR_LAMBDA)
    if warn:
        print(f"[WARN] {warn}")
    print(f"MMR picked {len(mmr_picks_local)} / {MMR_K} candidates.")

    mmr_df = pool_df.iloc[mmr_picks_local].copy()
    mmr_df["selection_method"] = "mmr"
    mmr_df["hypothesis"] = ""
    # Map back to all_df index for hypothesis dedup
    mmr_all_idx = [int(pool_all_idx[i]) for i in mmr_picks_local]

    # Hypothesis probes (search against full all_df, not just pool)
    mmr_vecs_stacked = all_vecs[np.array(mmr_all_idx, dtype=int)] if mmr_all_idx else np.zeros((0, 31), dtype=np.int8)
    probes = find_hypothesis_probes(all_df, all_vecs, A_vecs, pos_to_idx,
                                     already_picked_idx=mmr_all_idx,
                                     already_picked_vecs=mmr_vecs_stacked)
    probe_idx_set = set()
    probe_rows = []
    for idx_all, label in probes:
        if idx_all in probe_idx_set:
            continue
        probe_idx_set.add(idx_all)
        row = all_df.iloc[idx_all].copy()
        row["selection_method"] = "hypothesis"
        row["hypothesis"] = label
        probe_rows.append(row)
    probe_df = pd.DataFrame(probe_rows) if probe_rows else pd.DataFrame()
    print(f"Hypothesis probes: {len(probe_df)} / 5 hypotheses covered.")

    # Combine Track B
    frames = [mmr_df]
    if len(probe_df) > 0:
        frames.append(probe_df)
    B = pd.concat(frames, ignore_index=True)
    # Drop intra-B duplicates by mutations_unified
    B = B.drop_duplicates(subset=["mutations_unified"]).reset_index(drop=True)

    # Backfill with more MMR picks if B < target (walk down pool by MMR objective)
    if len(B) < TRACK_B_TARGET:
        picked_muts = set(B["mutations_unified"])
        # Re-run MMR with current B as starting point, continue picking
        # Simple approach: pick by objective from remaining pool until target reached
        B_vecs_list = [muts_to_vec(s, pos_to_idx) for s in B["mutations_unified"]]
        B_vecs_arr = np.array(B_vecs_list) if B_vecs_list else np.zeros((0, 31), dtype=np.int8)
        scores = pool_df["score"].to_numpy()
        s_min, s_max = scores.min(), scores.max()
        score_norm = (scores - s_min) / (s_max - s_min + 1e-9)
        hdist_A = hamming_matrix(pool_vecs, A_vecs).min(axis=1)
        while len(B) < TRACK_B_TARGET:
            # Distance to current B
            hd_B = hamming_matrix(pool_vecs, B_vecs_arr).min(axis=1) if len(B_vecs_arr) > 0 else \
                np.full(len(pool_df), 99)
            eligible = (hdist_A >= MIN_HAMMING_TO_A) & (hd_B >= MIN_HAMMING_TO_B)
            # Exclude already-picked by mut string
            excl = pool_df["mutations_unified"].isin(picked_muts).to_numpy()
            eligible &= ~excl
            if not eligible.any():
                print(f"[WARN] Cannot reach {TRACK_B_TARGET}: pool exhausted at {len(B)}")
                break
            combined = np.minimum(hdist_A, hd_B) / 31.0
            obj = MMR_LAMBDA * score_norm + (1 - MMR_LAMBDA) * combined
            obj[~eligible] = -np.inf
            pick = int(np.argmax(obj))
            row = pool_df.iloc[pick].copy()
            row["selection_method"] = "mmr_backfill"
            row["hypothesis"] = ""
            B = pd.concat([B, pd.DataFrame([row])], ignore_index=True)
            picked_muts.add(row["mutations_unified"])
            new_vec = muts_to_vec(row["mutations_unified"], pos_to_idx)[None, :]
            B_vecs_arr = np.vstack([B_vecs_arr, new_vec])

    if len(B) > TRACK_B_TARGET:
        B = B.head(TRACK_B_TARGET).reset_index(drop=True)

    # Add min_hamming_to_trackA column
    B_vecs = build_vecs(B, pos_to_idx)
    hd_B_to_A = hamming_matrix(B_vecs, A_vecs).min(axis=1)
    B["min_hamming_to_trackA"] = hd_B_to_A

    # Output column order
    cols = ["mutations_unified", "order",
            "log_pH74_pred", "log_pH6_pred", "log_ratio_pred", "score",
            "confidence", "positions_cdr1", "positions_cdr2", "positions_cdr3",
            "contains_D110H", "selection_method", "hypothesis", "min_hamming_to_trackA"]
    B = B[cols]
    B.to_csv(OUT_B, index=False)

    # Combined
    top_df["track"] = "A"
    top_df["selection_method"] = "trackA_greedy"
    top_df["hypothesis"] = ""
    top_df["min_hamming_to_trackA"] = 0
    B["track"] = "B"
    combined = pd.concat([top_df, B], ignore_index=True)
    combined_cols = ["track"] + cols
    combined = combined[combined_cols]
    combined.to_csv(OUT_COMBINED, index=False)

    # Report
    lines = [f"# sdab R4 Track B (Exploration) Report", ""]
    lines.append(f"- Track B size: **{len(B)}** / {TRACK_B_TARGET}")
    lines.append(f"- Combined (A+B): **{len(combined)}**")
    lines.append(f"- MMR picks: {(B['selection_method'].isin(['mmr', 'mmr_backfill'])).sum()}")
    lines.append(f"- Hypothesis picks: {(B['selection_method'] == 'hypothesis').sum()}")
    lines.append("")
    lines.append("## Diversity")
    lines.append(f"- Min Hamming to Track A: mean={hd_B_to_A.mean():.2f}, min={hd_B_to_A.min()}, max={hd_B_to_A.max()}")
    pairwise_B = hamming_matrix(B_vecs, B_vecs)
    np.fill_diagonal(pairwise_B, 99)
    lines.append(f"- Track B internal pairwise min Hamming: {pairwise_B.min()}")
    all100_vecs = np.vstack([A_vecs, B_vecs])
    pairwise_all = hamming_matrix(all100_vecs, all100_vecs)
    np.fill_diagonal(pairwise_all, 99)
    lines.append(f"- Combined 100 pairwise Hamming: min={pairwise_all.min()}, median={np.median(pairwise_all[pairwise_all < 99]):.1f}, max={pairwise_all.max()}")
    lines.append("")
    lines.append("## CDR Coverage (Track B)")
    for tag, col in [("CDR1", "positions_cdr1"), ("CDR2", "positions_cdr2"), ("CDR3", "positions_cdr3")]:
        lines.append(f"- Track B with {tag} mutation: {(B[col] > 0).sum()} / {len(B)}")
    lines.append("")
    lines.append("## Composition")
    lines.append(f"- Track B containing D110H: {int(B['contains_D110H'].sum())} / {len(B)}")
    contains_q100_v105 = B["mutations_unified"].apply(
        lambda s: "HQ100H" in s.split(";") and "HV105H" in s.split(";")).sum()
    lines.append(f"- Track B containing Q100H+V105H core: {int(contains_q100_v105)} / {len(B)}")
    lines.append(f"- Track B order distribution: {dict(B['order'].value_counts().sort_index())}")
    lines.append("")
    lines.append("## Hypothesis coverage")
    hypo_rows = B[B["selection_method"] == "hypothesis"]
    for _, r in hypo_rows.iterrows():
        lines.append(f"- **{r['hypothesis']}**: `{r['mutations_unified']}` (score={r['score']:.3f}, Hamming to A={r['min_hamming_to_trackA']})")
    lines.append("")
    lines.append("## Sample of Track B candidates (first 10, last 10)")
    lines.append("### First 10")
    lines.append(B.head(10)[["mutations_unified", "order", "score", "selection_method", "hypothesis", "min_hamming_to_trackA"]].to_string(index=False))
    lines.append("\n### Last 10")
    lines.append(B.tail(10)[["mutations_unified", "order", "score", "selection_method", "hypothesis", "min_hamming_to_trackA"]].to_string(index=False))

    OUT_REPORT.write_text("\n".join(lines) + "\n")

    print(f"\nWrote {OUT_B} ({len(B)} rows)")
    print(f"Wrote {OUT_COMBINED} ({len(combined)} rows)")
    print(f"Wrote {OUT_REPORT}")

    # Completion criteria
    print("\n=== Completion criteria ===")
    checks = []
    checks.append((len(B) == TRACK_B_TARGET, f"Track B = {TRACK_B_TARGET}: {len(B)}"))
    checks.append((len(combined) == 100, f"Combined = 100: {len(combined)}"))
    checks.append((pairwise_B.min() >= MIN_HAMMING_TO_B,
                   f"Track B internal min Hamming ≥ {MIN_HAMMING_TO_B}: {pairwise_B.min()}"))
    checks.append((hd_B_to_A.min() >= MIN_HAMMING_TO_A,
                   f"Track B → Track A min Hamming ≥ {MIN_HAMMING_TO_A}: {hd_B_to_A.min()}"))
    cdr_counts = [(B["positions_cdr1"] > 0).sum(),
                  (B["positions_cdr2"] > 0).sum(),
                  (B["positions_cdr3"] > 0).sum()]
    checks.append((all(c >= 5 for c in cdr_counts), f"Each CDR has ≥ 5 Track B candidates: {cdr_counts}"))
    checks.append((int(B["contains_D110H"].sum()) >= 10,
                   f"Track B with D110H ≥ 10: {int(B['contains_D110H'].sum())}"))
    non_core = len(B) - int(contains_q100_v105)
    checks.append((non_core >= 20,
                   f"Track B not containing Q100H+V105H core ≥ 20: {non_core}"))
    checks.append((hd_B_to_A.mean() >= 3.0,
                   f"Track B avg min Hamming to A ≥ 3: {hd_B_to_A.mean():.2f}"))

    all_pass = True
    for ok, msg in checks:
        print(f"{'[PASS]' if ok else '[FAIL]'} {msg}")
        if not ok:
            all_pass = False

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
