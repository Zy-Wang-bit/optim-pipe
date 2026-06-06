#!/usr/bin/env python3
"""Validate minimum evidence-ledger constraints."""

import sys

from analysis.initial_design_generation.run_dry_run import validate_evidence_ledger


def main() -> None:
    failures = validate_evidence_ledger()
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        sys.exit(1)
    print("PASS: minimum evidence-ledger constraints present")


if __name__ == "__main__":
    main()
