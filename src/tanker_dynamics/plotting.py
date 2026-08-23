"""Publication figure style. All figures are generated from stored result files."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "axes.linewidth": 0.7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.frameon": False,
    "lines.linewidth": 1.2,
    "figure.dpi": 160,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "axes.grid": True,
    "grid.linewidth": 0.4,
    "grid.alpha": 0.35,
    "axes.spines.top": False,
    "axes.spines.right": False,
}

# colour-blind safe
C_RIGID = "#B2182B"
C_COMPLIANT = "#2166AC"
C_ACCENT = "#0F7B6C"
C_GREY = "#4D4D4D"

# Elsevier single / 1.5 / double column widths (inches)
W1, W15, W2 = 3.54, 5.51, 7.48


def use_style():
    plt.rcParams.update(STYLE)


def save(fig, name: str, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"{name}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    print("figure:", name)
