#!/usr/bin/env python3
"""Build evidence ledger as part of the dry-run pipeline."""

from analysis.initial_design_generation.run_dry_run import build_evidence_ledger, load_inputs, read_config


def main() -> None:
    build_evidence_ledger(read_config(), load_inputs())


if __name__ == "__main__":
    main()
