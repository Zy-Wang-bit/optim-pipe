from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from active_learning import model
from active_learning import tools


def _build_candidate_pool(args: argparse.Namespace) -> None:
    out_csv = Path(args.out_csv)
    tools.reject_output_overwrite(out_csv, [Path(args.config_yaml), Path(args.wt_fasta)])
    config = tools.load_round_config(Path(args.config_yaml))
    pool = tools.build_candidate_pool_from_config(
        config,
        Path(args.wt_fasta),
        system_id=args.system_id,
        round_id=args.round_id,
        candidate_mode=args.candidate_mode,
    )
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pool.to_csv(out_csv, index=False)


def _build_sdab_warm_start(args: argparse.Namespace) -> None:
    out_csv = Path(args.out_csv)
    tools.reject_output_overwrite(
        out_csv,
        [
            Path(args.config_yaml),
            Path(args.wt_fasta),
            Path(args.legacy_r4_training),
            Path(args.direct_wet_labels) if args.direct_wet_labels else None,
        ],
    )
    config = tools.load_round_config(Path(args.config_yaml))
    rows = tools.build_sdab_warm_start_rows(
        config,
        Path(args.wt_fasta),
        Path(args.legacy_r4_training),
        Path(args.direct_wet_labels) if args.direct_wet_labels else None,
    )
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(out_csv, index=False)


def _run_round(args: argparse.Namespace) -> None:
    config = tools.load_round_config(Path(args.config_yaml))
    config["round_id"] = args.round_id or config.get("round_id")
    if args.top_k is not None:
        config["acquisition"]["top_k"] = args.top_k
    round_id = config["round_id"]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    output_paths = model.round_output_paths(out_dir, round_id)
    tools.reject_output_paths_overwrite_inputs(
        output_paths,
        [
            path
            for path in [
                Path(args.config_yaml),
                Path(args.training_csv),
                Path(args.candidate_csv),
                Path(args.embedding_csv),
                Path(args.failed_synthesis_csv) if args.failed_synthesis_csv else None,
                Path(args.invalid_candidates_csv) if args.invalid_candidates_csv else None,
                Path(args.previously_selected_csv) if args.previously_selected_csv else None,
            ]
            if path is not None
        ],
    )

    artifacts = model.run_active_learning_round(
        config,
        pd.read_csv(args.training_csv),
        pd.read_csv(args.candidate_csv),
        tools.load_embedding_table(Path(args.embedding_csv)),
        failed_synthesis=tools.read_optional_csv(args.failed_synthesis_csv),
        invalid_candidates=tools.read_optional_csv(args.invalid_candidates_csv),
        previously_selected=tools.read_optional_csv(args.previously_selected_csv),
    )

    artifacts["training"].to_csv(output_paths["training"], index=False)
    artifacts["candidate"].to_csv(output_paths["candidate"], index=False)
    artifacts["prediction"].to_csv(output_paths["prediction"], index=False)
    artifacts["selected"].to_csv(output_paths["selected"], index=False)
    tools.write_json(output_paths["model_config"], artifacts["model_config"])
    tools.write_json(
        output_paths["provenance"],
        tools.build_round_provenance(
            config,
            {
                "wt_fasta": Path(config["wt_fasta_path"]),
                "training_csv": Path(args.training_csv),
                "candidate_pool": Path(args.candidate_csv),
                "embedding_csv": Path(args.embedding_csv),
            },
        ),
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Active-learning utilities for pH-switch sdAb design.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    pool = subparsers.add_parser("build-candidate-pool", help="Build a full single-mutant candidate pool.")
    pool.add_argument("--config-yaml", required=True)
    pool.add_argument("--wt-fasta", required=True)
    pool.add_argument("--system-id")
    pool.add_argument("--round-id")
    pool.add_argument("--candidate-mode")
    pool.add_argument("--out-csv", required=True)
    pool.set_defaults(func=_build_candidate_pool)

    warm = subparsers.add_parser("build-sdab-warm-start", help="Build sdAb warm-start training rows.")
    warm.add_argument("--config-yaml", required=True)
    warm.add_argument("--wt-fasta", required=True)
    warm.add_argument("--legacy-r4-training", required=True)
    warm.add_argument("--direct-wet-labels")
    warm.add_argument("--out-csv", required=True)
    warm.set_defaults(func=_build_sdab_warm_start)

    run = subparsers.add_parser("run-round", help="Run one active-learning selection round.")
    run.add_argument("--config-yaml", required=True)
    run.add_argument("--training-csv", required=True)
    run.add_argument("--candidate-csv", required=True)
    run.add_argument("--embedding-csv", required=True)
    run.add_argument("--failed-synthesis-csv")
    run.add_argument("--invalid-candidates-csv")
    run.add_argument("--previously-selected-csv")
    run.add_argument("--round-id")
    run.add_argument("--top-k", type=int)
    run.add_argument("--out-dir", required=True)
    run.set_defaults(func=_run_round)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
