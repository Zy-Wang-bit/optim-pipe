import subprocess
import sys


def test_package_cli_lists_subcommands():
    result = subprocess.run(
        [sys.executable, "-m", "active_learning", "--help"],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert "build-candidate-pool" in result.stdout
    assert "build-sdab-warm-start" in result.stdout
    assert "run-round" in result.stdout
