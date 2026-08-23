"""Run the entire study end to end.

    python3 scripts/run_all.py            # full pipeline
    python3 scripts/run_all.py --identify # also re-run the parameter identification

The identification stage is skipped by default because its result is already stored in
config/vehicle.yaml; re-running it takes several minutes and returns the same values.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The manuscript pipeline. Scripts 02-08 belong to the SEPARATE 3^3 factorial study
# (a future experimental paper) and are deliberately NOT part of it.
STAGES = [
    ("11_archive_studies.py", "manuscript studies A1-A7 (archive-locked model)"),
    ("12_archive_figures.py", "manuscript figures"),
]


def run(script: str, label: str) -> None:
    print(f"\n=== {label} ({script}) ===", flush=True)
    t0 = time.time()
    subprocess.run([sys.executable, str(ROOT / "scripts" / script)], check=True)
    print(f"--- done in {time.time() - t0:.1f} s", flush=True)


def main() -> None:
    # the archive-extraction step must run first: everything downstream reads it
    subprocess.run([sys.executable, str(ROOT / "tools" / "parse_mathcad_archive.py")],
                   check=True)
    for script, label in STAGES:
        run(script, label)
    print("\n=== running the test suite ===", flush=True)
    subprocess.run([sys.executable, "-m", "pytest", str(ROOT / "tests"), "-q"], check=True)
    print("\nAll stages completed. Results in results/, figures in results/figures/.")


if __name__ == "__main__":
    main()
