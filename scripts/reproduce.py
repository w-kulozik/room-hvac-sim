"""Regenerate every figure and result table, then check the docs against them.

Usage::

    python scripts/reproduce.py            # run all studies (about 5 minutes)
    python scripts/reproduce.py --check    # also fail if any regenerated file differs
                                           # from the committed version (needs git)

The studies are deterministic: on a clean checkout, ``--check`` reports no difference.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STUDIES = [
    "run_baseline.py",
    "compare_rated.py",
    "compare_controllers.py",
    "compare_humidity.py",
    "compare_mpc.py",
    "compare_ventilation.py",
    "compare_season.py",
]


def run(cmd: list[str]) -> None:
    """Run a command from the repository root and stop on the first failure."""
    result = subprocess.run(cmd, cwd=ROOT, check=False)
    if result.returncode != 0:
        sys.exit(f"failed: {' '.join(cmd)}")


def main() -> None:
    """Run the studies, the docs-consistency test and, optionally, the git check."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="fail if regenerated figures differ from the committed ones")
    args = parser.parse_args()

    for name in STUDIES:
        start = time.perf_counter()
        print(f"== {name}", flush=True)
        run([sys.executable, str(ROOT / "scripts" / name)])
        print(f"   done in {time.perf_counter() - start:.0f} s", flush=True)

    run([sys.executable, "-m", "pytest", "-q", "tests/test_docs_consistency.py"])

    if args.check:
        diff = subprocess.run(["git", "status", "--porcelain", "--", "figures"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout
        if diff:
            sys.exit(f"regenerated figures differ from the committed ones:\n{diff}")
        print("All figures and tables match the committed versions.")


if __name__ == "__main__":
    main()
