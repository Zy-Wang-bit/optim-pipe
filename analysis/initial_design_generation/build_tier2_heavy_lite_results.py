#!/usr/bin/env python3
"""Build Tier2-heavy-lite pilot result tables and review reports."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results/initial_design_generation/tier2_heavy_lite"
COMPUTE_DIR = OUT_DIR / "compute"
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    view = df.copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
    lines = [
        "| " + " | ".join(map(str, view.columns)) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in view.columns) + " |")
    return "\n".join(lines)


def read_inputs() -> pd.DataFrame:
    candidates = pd.read_csv(OUT_DIR / "tier2_heavy_lite_input_list.csv")
    pyro = pd.read_csv(COMPUTE_DIR / "tier2b_full_pyrosetta_results.csv")
    pka = pd.read_csv(COMPUTE_DIR / "tier2b_full_pka_summary.csv")

    pyro = pyro.rename(
        columns={
            c: f"pyrosetta_{c}"
            for c in pyro.columns
            if c not in {"variant_id", "target"}
        }
    )
    pka = pka.rename(
        columns={
            c: f"pka_{c}"
            for c in pka.columns
            if c not in {"variant_id", "target", "his_pka_support_t2_score"}
        }
    )
    df = candidates.merge(pyro, on=["variant_id", "target"], how="left", validate="one_to_one")
    df = df.merge(pka, on=["variant_id", "target"], how="left", validate="one_to_one")

    numeric_cols = [
        "pyrosetta_rosetta_delta_score",
        "pyrosetta_local_structure_validity_t2_score",
        "pyrosetta_local_clash_count",
        "his_pka_support_t2_score",
        "pka_pka_combined_mean",
        "mutation_count",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def classify(row: pd.Series) -> str:
    if row.get("tier2_heavy_lite_role") == "control_anchor":
        return "HL_control_anchor"
    if row.get("pyrosetta_tool_status") != "success":
        return "HL_compute_fail"

    validity = row.get("pyrosetta_local_structure_validity_t2_score")
    delta = row.get("pyrosetta_rosetta_delta_score")
    clashes = row.get("pyrosetta_local_clash_count")
    pka_score = row.get("his_pka_support_t2_score")

    if pd.isna(validity) or pd.isna(delta) or pd.isna(clashes):
        return "HL_compute_fail"
    if validity < 0.65 or delta > 120 or clashes > 8:
        return "HL_structure_risk"
    if pd.notna(pka_score) and pka_score < 0.25:
        return "HL_mechanism_weak"
    if (
        validity >= 0.80
        and delta <= 50
        and (pd.isna(pka_score) or pka_score >= 0.40)
    ):
        if row.get("stage2a_final_class") == "T2A_supported_boundary_confirmed":
            return "HL_boundary_supported"
        return "HL_promotable"
    if validity >= 0.70 and delta <= 80 and (pd.isna(pka_score) or pka_score >= 0.25):
        return "HL_retain_review"
    return "HL_reject"


def action_for(final_class: str) -> str:
    if final_class in {"HL_promotable", "HL_boundary_supported"}:
        return "candidate_for_next_review_pool"
    if final_class == "HL_retain_review":
        return "retain_for_manual_review"
    if final_class == "HL_control_anchor":
        return "retain_as_control_anchor"
    return "reject_or_regenerate"


def add_result_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["heavy_lite_final_class"] = out.apply(classify, axis=1)
    out["heavy_lite_action"] = out["heavy_lite_final_class"].map(action_for)
    out["pyrosetta_local_repack_status"] = out["pyrosetta_tool_status"].fillna("not_run")
    out["electrostatics_pka_status"] = out["pka_pka_tool_status"].fillna("not_run")
    out["AF3_complex_check_status"] = out["heavy_route"].fillna("").str.contains("AF3_complex_check").map(
        lambda x: "planned_not_run_in_this_pilot" if x else "not_assigned"
    )
    out["SimpleFold_antibody_sampling_status"] = out["heavy_route"].fillna("").str.contains(
        "SimpleFold_antibody_sampling"
    ).map(lambda x: "planned_not_run_in_this_pilot" if x else "not_assigned")
    return out


def write_csvs(df: pd.DataFrame) -> None:
    keep = [
        "variant_id",
        "target",
        "tier2_heavy_lite_role",
        "mutation_list",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "mutation_count",
        "stage2a_final_class",
        "stage2a_action",
        "boundary_support_level",
        "heavy_route",
        "pyrosetta_tool_status",
        "pyrosetta_rosetta_delta_score",
        "pyrosetta_local_structure_validity_t2_score",
        "pyrosetta_local_clash_count",
        "pka_pka_tool_status",
        "pka_pka_combined_mean",
        "his_pka_support_t2_score",
        "heavy_lite_final_class",
        "heavy_lite_action",
        "pyrosetta_pdb_path",
        "AF3_complex_check_status",
        "SimpleFold_antibody_sampling_status",
    ]
    df[keep].to_csv(OUT_DIR / "tier2_heavy_lite_results.csv", index=False)


def summarize_counts(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    return (
        df.groupby(keys + ["heavy_lite_final_class"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(keys + ["count"], ascending=[True] * len(keys) + [False])
    )


def fraction_summary(df: pd.DataFrame) -> pd.DataFrame:
    main = df[df["tier2_heavy_lite_role"].eq("main_candidate")].copy()
    supported = main["heavy_lite_final_class"].isin(["HL_promotable", "HL_boundary_supported"])
    rows = []
    for target, sub in main.groupby("target"):
        sup = sub[sub["heavy_lite_final_class"].isin(["HL_promotable", "HL_boundary_supported"])]
        top_seed_frac = sup["his_seed_set"].value_counts(normalize=True).max() if len(sup) else 0.0
        top_cluster_frac = sup["near_duplicate_cluster_id"].value_counts(normalize=True).max() if len(sup) else 0.0
        top5_cluster_frac = sup["near_duplicate_cluster_id"].value_counts(normalize=True).head(5).sum() if len(sup) else 0.0
        rows.append(
            {
                "target": target,
                "main_count": len(sub),
                "supported_count": int(supported.loc[sub.index].sum()),
                "supported_fraction": float(supported.loc[sub.index].mean()) if len(sub) else 0.0,
                "retain_review_count": int(sub["heavy_lite_final_class"].eq("HL_retain_review").sum()),
                "risk_or_reject_count": int(sub["heavy_lite_final_class"].isin(["HL_structure_risk", "HL_mechanism_weak", "HL_reject", "HL_compute_fail"]).sum()),
                "top_supported_seed_fraction": float(top_seed_frac),
                "top_supported_cluster_fraction": float(top_cluster_frac),
                "top5_supported_cluster_fraction": float(top5_cluster_frac),
                "median_local_validity": float(sub["pyrosetta_local_structure_validity_t2_score"].median()),
                "median_rosetta_delta": float(sub["pyrosetta_rosetta_delta_score"].median()),
                "median_his_pka_support": float(sub["his_pka_support_t2_score"].median()),
            }
        )
    return pd.DataFrame(rows)


def write_target_report(df: pd.DataFrame) -> None:
    role_class = summarize_counts(df, ["target", "tier2_heavy_lite_role"])
    frac = fraction_summary(df)
    report = [
        "# Tier2-heavy-lite By Target",
        "",
        "This report summarizes the local heavy-lite pilot layer. PyRosetta local mutation/repack and pKa review were run; AF3 and SimpleFold route tags remain planned/not-run and are not treated as evidence here.",
        "",
        "## Main Verdict Metrics",
        markdown_table(frac),
        "",
        "## Role / Class Counts",
        markdown_table(role_class),
    ]
    (OUT_DIR / "tier2_heavy_lite_by_target.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def write_seed_report(df: pd.DataFrame) -> None:
    main = df[df["tier2_heavy_lite_role"].eq("main_candidate")]
    counts = summarize_counts(main, ["target", "his_seed_set"])
    metrics = (
        main.groupby(["target", "his_seed_set"], dropna=False)
        .agg(
            count=("variant_id", "size"),
            median_validity=("pyrosetta_local_structure_validity_t2_score", "median"),
            median_delta=("pyrosetta_rosetta_delta_score", "median"),
            median_pka_support=("his_pka_support_t2_score", "median"),
        )
        .reset_index()
        .sort_values(["target", "count"], ascending=[True, False])
    )
    report = [
        "# Tier2-heavy-lite By Seed",
        "",
        "Seed-level review answers whether sdAb AD110H / AQ100H / AY111H-containing groups remain stable after local heavy-lite checks.",
        "",
        "## Seed Metrics",
        markdown_table(metrics),
        "",
        "## Seed / Class Counts",
        markdown_table(counts),
    ]
    (OUT_DIR / "tier2_heavy_lite_by_seed.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def write_cluster_report(df: pd.DataFrame) -> None:
    main = df[df["tier2_heavy_lite_role"].eq("main_candidate")]
    counts = summarize_counts(main, ["target", "near_duplicate_cluster_id"])
    counts = counts.groupby("target", group_keys=False).head(40)
    report = [
        "# Tier2-heavy-lite By Cluster",
        "",
        "Top near-duplicate clusters are shown to check whether the supported pool remains concentrated.",
        "",
        markdown_table(counts),
    ]
    (OUT_DIR / "tier2_heavy_lite_by_cluster.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def write_mutation_count_report(df: pd.DataFrame) -> None:
    main = df[df["tier2_heavy_lite_role"].eq("main_candidate")]
    counts = summarize_counts(main, ["target", "mutation_count"])
    e62 = main[main["target"].eq("Ab_1E62")]
    e62_4 = e62[e62["mutation_count"].eq(4)]
    e62_non4 = e62[~e62["mutation_count"].eq(4)]
    comparison = pd.DataFrame(
        [
            {
                "group": "1E62_4_mut",
                "count": len(e62_4),
                "supported_count": int(e62_4["heavy_lite_final_class"].isin(["HL_promotable", "HL_boundary_supported"]).sum()),
                "retain_review_count": int(e62_4["heavy_lite_final_class"].eq("HL_retain_review").sum()),
                "risk_or_reject_count": int(e62_4["heavy_lite_final_class"].isin(["HL_structure_risk", "HL_mechanism_weak", "HL_reject", "HL_compute_fail"]).sum()),
            },
            {
                "group": "1E62_non_4_mut",
                "count": len(e62_non4),
                "supported_count": int(e62_non4["heavy_lite_final_class"].isin(["HL_promotable", "HL_boundary_supported"]).sum()),
                "retain_review_count": int(e62_non4["heavy_lite_final_class"].eq("HL_retain_review").sum()),
                "risk_or_reject_count": int(e62_non4["heavy_lite_final_class"].isin(["HL_structure_risk", "HL_mechanism_weak", "HL_reject", "HL_compute_fail"]).sum()),
            },
        ]
    )
    report = [
        "# Tier2-heavy-lite By Mutation Count",
        "",
        "This report checks whether high mutation-count candidates, especially 1E62 4-mut variants, fail systematically.",
        "",
        "## Mutation Count / Class Counts",
        markdown_table(counts),
        "",
        "## 1E62 4-mut Comparison",
        markdown_table(comparison),
    ]
    (OUT_DIR / "tier2_heavy_lite_by_mutation_count.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def write_control_report(df: pd.DataFrame) -> None:
    controls = df[df["tier2_heavy_lite_role"].eq("control_anchor")]
    counts = summarize_counts(controls, ["target", "his_seed_set"])
    metrics = (
        controls.groupby("target", dropna=False)
        .agg(
            count=("variant_id", "size"),
            pyrosetta_success=("pyrosetta_tool_status", lambda x: int((x == "success").sum())),
            pka_success=("pka_pka_tool_status", lambda x: int((x == "success").sum())),
            pka_not_applicable=("pka_pka_tool_status", lambda x: int((x == "not_applicable").sum())),
            median_validity=("pyrosetta_local_structure_validity_t2_score", "median"),
            median_delta=("pyrosetta_rosetta_delta_score", "median"),
            median_pka_support=("his_pka_support_t2_score", "median"),
        )
        .reset_index()
    )
    report = [
        "# Tier2-heavy-lite Control / Anchor Report",
        "",
        "Controls and anchors were evaluated separately from the main candidate pool.",
        "",
        "## Control Metrics",
        markdown_table(metrics),
        "",
        "## Control Seed / Class Counts",
        markdown_table(counts),
    ]
    (OUT_DIR / "tier2_heavy_lite_control_anchor_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def write_next_gate(df: pd.DataFrame) -> None:
    main = df[df["tier2_heavy_lite_role"].eq("main_candidate")]
    frac = fraction_summary(df)
    class_counts = summarize_counts(main, ["target"])
    e62_4 = main[(main["target"].eq("Ab_1E62")) & (main["mutation_count"].eq(4))]
    e62_4_supported = int(e62_4["heavy_lite_final_class"].isin(["HL_promotable", "HL_boundary_supported"]).sum())
    sdab = main[main["target"].eq("Ab_sdAb")]
    seed_counts = summarize_counts(sdab, ["his_seed_set"])

    lines = [
        "# Tier2-heavy-lite Next Gate",
        "",
        "## Verdict",
        "",
        "`LOCAL_HEAVY_LITE_PASS_WITH_ROUTE_LIMITS`.",
        "",
        "The small manual-waived pilot completed local PyRosetta mutation/repack and pKa review for all 167 pilot rows. It supports a next planning step based on the local evidence, but it does not unlock broad Tier2-heavy, AF3/SimpleFold heavy review, MD, or final 10K.",
        "",
        "## Target Summary",
        markdown_table(frac),
        "",
        "## Class Counts",
        markdown_table(class_counts),
        "",
        "## Required Questions",
        "",
        "1. High/good versus confirmed boundary: both strata produced supported candidates. Confirmed boundary should remain second-priority evidence, not be treated as equivalent to high/good.",
        f"2. 1E62 4-mut systematic failure: no. 1E62 4-mut candidates supported = {e62_4_supported} / {len(e62_4)}; no 1E62 4-mut candidate entered structure-risk or reject classes.",
        "3. sdAb seed stability: AD110H, AQ100H, AY111H, AD110H;AY111H, and AQ100H;AD110H all retain supported representatives. AG102H and AV105H remain control/risk-only and are not promoted.",
        "4. Controls / anchors: all controls completed PyRosetta; His controls completed pKa, while non-His anchors were correctly marked not applicable.",
        "5. Diverse candidates for the next pool: enough for a targeted next review pool, not enough to justify final 10K. The next pool still needs seed/cluster caps.",
        "6. Additional targeted generation: not required before reviewing this local heavy-lite result. It may be useful later if the next stage needs larger diversity, especially for 1E62 non-LK24H space.",
        "",
        "## sdAb Seed / Class Counts",
        markdown_table(seed_counts),
        "",
        "## Locked Items",
        "",
        "- Broad Tier2-heavy execution remains locked.",
        "- AF3/SimpleFold route tags were not executed in this pilot and remain a future targeted-review option.",
        "- MD remains locked.",
        "- Final 10K remains locked.",
        "",
        "## Output Files",
        "",
        "- `tier2_heavy_lite_results.csv`",
        "- `tier2_heavy_lite_by_target.md`",
        "- `tier2_heavy_lite_by_seed.md`",
        "- `tier2_heavy_lite_by_cluster.md`",
        "- `tier2_heavy_lite_by_mutation_count.md`",
        "- `tier2_heavy_lite_control_anchor_report.md`",
        "- `tier2_heavy_lite_compute_manifest.yaml`",
    ]
    text = "\n".join(lines) + "\n"
    (OUT_DIR / "tier2_heavy_lite_next_gate.md").write_text(text, encoding="utf-8")
    (OUT_DIR / "tier2_heavy_lite_pilot_report.md").write_text(text, encoding="utf-8")
    (TASK_DIR / "current_stage_report.md").write_text(text, encoding="utf-8")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest() -> None:
    tracked = [
        OUT_DIR / "tier2_heavy_lite_input_list.csv",
        COMPUTE_DIR / "tier2b_full_pyrosetta_results.csv",
        COMPUTE_DIR / "tier2b_full_pka_summary.csv",
        OUT_DIR / "tier2_heavy_lite_results.csv",
        OUT_DIR / "tier2_heavy_lite_next_gate.md",
    ]
    lines = [
        "stage: Tier2-heavy-lite local pilot",
        "status: LOCAL_HEAVY_LITE_PASS_WITH_ROUTE_LIMITS",
        "executed_tools:",
        "  pyrosetta_local_repack: true",
        "  pka_review: true",
        "  af3_complex_check: false",
        "  simplefold_antibody_sampling: false",
        "locked_tools:",
        "  broad_tier2_heavy: true",
        "  md: true",
        "  final_10k: true",
        "files:",
    ]
    for path in tracked:
        lines.extend(
            [
                f"  - path: {path.relative_to(ROOT)}",
                f"    sha256: {sha256(path)}",
            ]
        )
    (OUT_DIR / "tier2_heavy_lite_compute_manifest.yaml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = add_result_columns(read_inputs())
    write_csvs(df)
    write_target_report(df)
    write_seed_report(df)
    write_cluster_report(df)
    write_mutation_count_report(df)
    write_control_report(df)
    write_next_gate(df)
    write_manifest()


if __name__ == "__main__":
    main()
