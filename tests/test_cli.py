"""CLI smoke tests: exercise the user-facing entry points and check that the
canonical numbers surface correctly. These cover the path a reviewer runs first
(`rlds-paired verify`)."""

import os
import subprocess
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run(args, cwd=None):
    env = {
        **os.environ,
        "PYTHONPATH": _REPO + os.pathsep + os.environ.get("PYTHONPATH", ""),
    }
    return subprocess.run(
        [sys.executable, "-m", "rlds_paired.cli", *args],
        capture_output=True,
        text=True,
        encoding="utf-8",
        cwd=cwd,
        env=env,
    )


def test_cli_verify_passes_and_shows_canonical_c():
    r = _run(["verify"])
    assert r.returncode == 0, r.stderr
    assert "All checks passed" in r.stdout
    assert "-0.194" in r.stdout
    assert "0.97" in r.stdout


def test_cli_spectral_table():
    r = _run(["spectral", "--table"])
    assert r.returncode == 0, r.stderr
    # Fenichel bound violated across the table
    assert "Fenichel" in r.stdout
    assert "NO" in r.stdout


def test_cli_scan_runs(tmp_path):
    # scan saves a figure in the cwd; run inside tmp_path so it is cleaned up
    r = _run(["scan", "--n", "8"], cwd=str(tmp_path))
    assert r.returncode == 0, r.stderr
    assert "c  =" in r.stdout
    assert "coarse grid" in r.stdout  # low-n hint shown for n < 40
