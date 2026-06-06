#!/usr/bin/env python3
"""Report dry-run validation status from the generated report."""

from pathlib import Path


def main() -> None:
    path = Path("results/initial_design_generation/reports/dry_run_validation_report.md")
    print(path.read_text())


if __name__ == "__main__":
    main()
