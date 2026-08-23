"""Figures for the manuscript. Every panel reads a file produced by 09_paper_studies.py."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.plotting import (use_style, save, C_RIGID, C_COMPLIANT,  # noqa: E402
                                      C_ACCENT, C_GREY, W1, W15, W2)

P = ROOT / "results" / "paper"
FIG = ROOT / "results" / "paper_figures"
TAB = ROOT / "results" / "tables"


def fig_workflow():
    fig, ax = plt.subplots(figsize=(W1 * 1.1, 3.6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 16); ax.axis("off")
    steps = ["Source planar formulation\n(Lagrange, eq. (1)-(2))",
             "Reduced-order longitudinal model",
             "Identification of 4 unreported\nparameters + 1 damping term",
             "Reproduction of the published results",
             "Damping-term ablation\n(literal vs augmented)",
             "Forward comparison\n(matched engine)",
             "Inverse comparison\n(matched manoeuvre)",
             "Mechanism: effective accelerated mass"]
    y = 15.2
    for i, s in enumerate(steps):
        col = C_ACCENT if i in (5, 6, 7) else C_GREY
        ax.add_patch(plt.Rectangle((0.4, y - 1.25), 9.2, 1.25, fc=col, alpha=0.10,
                                   ec=col, lw=0.9))
        ax.text(5.0, y - 0.63, s, ha="center", va="center", fontsize=7.0)
        if i < len(steps) - 1:
            ax.annotate("", xy=(5.0, y - 1.72), xytext=(5.0, y - 1.28),
                        arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C_GREY))
        y -= 1.9
    save(fig, "pf03_workflow", FIG)


def fig_ablation():
    d = pd.read_csv(P / "s2_damping_ablation.csv")
    fig, axes = plt.subplots(1, 2, figsize=(W15, 2.3), constrained_layout=True)
    for ax, cfgn, lab in zip(axes, ("rigid", "compliant"),
                             ("rigid mounting", "segmented mounting")):
        for (mod, g), col, mk in zip(d.groupby("model"), (C_RIGID, C_COMPLIANT), ("o", "s")):
            ax.plot(g.m_v_kg, g[f"err_{cfgn}_pct"], marker=mk, ms=4, color=col,
                    label="literal ($D_{kch}=0$)" if mod.startswith("A") else "augmented")
        ax.axhspan(-5, 5, color=C_ACCENT, alpha=0.12)
        ax.axhline(0, color=C_GREY, lw=0.7)
        ax.set_xlabel(r"cargo per tank $m_v$, kg")
        ax.set_ylabel("error vs published, %")
        ax.set_title(lab, fontsize=8)
    axes[0].legend(fontsize=6.5)
    save(fig, "pf04_damping_ablation", FIG)


def fig_forward():
    z = np.load(P / "s3_timeseries.npz")
    s = pd.read_csv(P / "s3_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.2), constrained_layout=True)
    for cfgn, col in (("rigid", C_RIGID), ("compliant", C_COMPLIANT)):
        t, v, F, a, u1, u2 = z[cfgn]
        axes[0].plot(t, F / 1e3, color=col, label=cfgn.replace("compliant", "segmented"))
        axes[1].plot(t, v, color=col)
    axes[0].set_ylabel(r"$F_{kch}$, kN"); axes[0].set_xlabel("time, s"); axes[0].legend()
    axes[1].set_ylabel(r"$v$, m s$^{-1}$"); axes[1].set_xlabel("time, s")
    axes[0].set_xlim(0, 4); axes[1].set_xlim(0, 8)
    w = 260
    axes[2].bar(s.cargo_per_tank_kg - w, s.rigid / 1e3, width=2 * w, color=C_RIGID,
                label="rigid")
    axes[2].bar(s.cargo_per_tank_kg + w, s.compliant / 1e3, width=2 * w,
                color=C_COMPLIANT, label="segmented")
    axes[2].set_xlabel(r"cargo per tank, kg"); axes[2].set_ylabel(r"$F_{kch}^{peak}$, kN")
    axes[2].set_xticks(s.cargo_per_tank_kg); axes[2].legend(fontsize=6.5)
    for a_, t_ in zip(axes, ("(a)", "(b)", "(c)")):
        a_.set_title(t_, loc="left", fontweight="bold")
    save(fig, "pf05_forward", FIG)


def fig_inverse():
    z = np.load(P / "s4_inverse_timeseries.npz")
    t = z["t"]
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.2), constrained_layout=True)
    axes[0].plot(t, z["a"], color=C_GREY)
    axes[0].set_ylabel(r"prescribed $a_{ref}$, m s$^{-2}$")
    axes[1].plot(t, z["F_rigid"] / 1e3, color=C_RIGID, label="rigid")
    axes[1].plot(t, z["F_seg"] / 1e3, color=C_COMPLIANT, label="segmented")
    axes[1].set_ylabel(r"required $F$, kN"); axes[1].legend()
    axes[2].plot(t, z["P_rigid"] / 1e3, color=C_RIGID)
    axes[2].plot(t, z["P_seg"] / 1e3, color=C_COMPLIANT)
    axes[2].set_ylabel(r"required $P$, kW")
    for a_, t_ in zip(axes, ("(a)", "(b)", "(c)")):
        a_.set_xlabel("time, s"); a_.set_title(t_, loc="left", fontweight="bold")
    save(fig, "pf06_inverse_demand", FIG)


def fig_effective_mass():
    z = np.load(P / "s4_inverse_timeseries.npz")
    t, me = z["t"], z["m_eff_seg"]
    m_tot = float(np.nanmax(z["m_eff_rigid"]))
    fig, axes = plt.subplots(1, 2, figsize=(W15, 2.3), constrained_layout=True)
    ax = axes[0]
    ax.axhline(m_tot / 1e3, color=C_RIGID, lw=1.2, label="rigid: $m_{tot}$")
    ax.plot(t, me / 1e3, color=C_COMPLIANT, label="segmented: $m_{eff}(t)$")
    ax.set_xlabel("time, s"); ax.set_ylabel(r"effective accelerated mass, t")
    ax.set_xlim(0, 9); ax.legend(fontsize=6.5)
    ax = axes[1]
    ax.plot(t, z["u1"] * 1e3, color=C_COMPLIANT, label=r"$u_1$")
    ax.plot(t, z["u2"] * 1e3, color=C_ACCENT, ls="--", label=r"$u_2$")
    ax.set_xlabel("time, s"); ax.set_ylabel("relative travel, mm"); ax.legend(fontsize=6.5)
    for a_, t_ in zip(axes, ("(a)", "(b)")):
        a_.set_title(t_, loc="left", fontweight="bold")
    save(fig, "pf07_effective_mass", FIG)


def fig_no_benefit():
    d = pd.read_csv(P / "s6_duration_sweep.csv")
    p7 = pd.read_csv(P / "s7_parametric.csv")
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.2), constrained_layout=True)
    ax = axes[0]
    ax.axhline(0, color=C_GREY, lw=0.8)
    ax.plot(d.T_over_Tn, d.dF_peak_pct, "o-", ms=3, color=C_COMPLIANT, label="force")
    ax.plot(d.T_over_Tn, d.dP_peak_pct, "s--", ms=3, color=C_RIGID, label="power")
    ax.set_xscale("log"); ax.set_xlabel(r"$T/T_n$")
    ax.set_ylabel("demand reduction, %"); ax.legend(fontsize=6.5)
    ax = axes[1]
    ax.axhline(0, color=C_GREY, lw=0.8)
    for (sw, g), col, mk in zip(p7.groupby("sweep"),
                                (C_RIGID, C_COMPLIANT, C_ACCENT, C_GREY),
                                ("o", "s", "^", "d")):
        ax.plot(np.arange(len(g)), g.dF_peak_pct, marker=mk, ms=3.5, lw=0.9,
                color=col, label=sw)
    ax.set_xlabel("sweep index"); ax.set_ylabel("force reduction, %")
    ax.legend(fontsize=5.8)
    ax = axes[2]
    a = pd.read_csv(P / "s8_asymmetric.csv")
    ax.axhline(0, color=C_GREY, lw=0.8)
    ax.bar(np.arange(len(a)), a.dF_peak_pct, color=C_COMPLIANT)
    ax.set_xticks(np.arange(len(a)))
    ax.set_xticklabels([f"{int(100*r.front_fill)}/{int(100*r.rear_fill)}"
                        for _, r in a.iterrows()], fontsize=6, rotation=45)
    ax.set_xlabel("front/rear fill, %"); ax.set_ylabel("force reduction, %")
    for a_, t_ in zip(axes, ("(a)", "(b)", "(c)")):
        a_.set_title(t_, loc="left", fontweight="bold")
    save(fig, "pf08_matched_no_benefit", FIG)


def main():
    use_style()
    FIG.mkdir(parents=True, exist_ok=True)
    fig_workflow(); fig_ablation(); fig_forward()
    fig_inverse(); fig_effective_mass(); fig_no_benefit()
    # reuse the concept and model schematics
    import shutil
    for n in ("fig01_system_concept", "fig02_lumped_model"):
        for ext in ("pdf", "png"):
            shutil.copy(ROOT / "results" / "figures" / f"{n}.{ext}",
                        FIG / f"pf{'01' if n.startswith('fig01') else '02'}_"
                              f"{n.split('_', 1)[1]}.{ext}")
    print("figures ->", FIG)


if __name__ == "__main__":
    main()
