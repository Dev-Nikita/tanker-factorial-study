"""Manuscript figures for the archive-locked study. All data from results/archive_study/."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.plotting import (use_style, save, C_RIGID, C_COMPLIANT,  # noqa: E402
                                      C_ACCENT, C_GREY, W1, W15, W2)
from tanker_dynamics.archive_model import V1, PT1, ArchiveParams            # noqa: E402

P = ROOT / "results" / "archive_study"
FIG = ROOT / "results" / "archive_figures"


def _spring(ax, x0, x1, y, n=6, amp=0.05, lw=1.0, color="k"):
    x = np.linspace(x0, x1, 240)
    s = np.zeros_like(x)
    m = (x > x0 + 0.18 * (x1 - x0)) & (x < x1 - 0.18 * (x1 - x0))
    if m.sum() > 2:
        tt = (x[m] - x[m][0]) / (x[m][-1] - x[m][0])
        s[m] = amp * np.sin(2 * np.pi * n * tt)
    ax.plot(x, y + s, color=color, lw=lw)


def fig01_concept():
    fig, axes = plt.subplots(2, 1, figsize=(W15, 3.3))
    for ax, mode in zip(axes, ("rigid", "segmented")):
        ax.set_xlim(-0.2, 7.0); ax.set_ylim(-0.4, 1.35); ax.axis("off")
        ax.plot([-0.2, 7.0], [0, 0], color=C_GREY, lw=1.2)
        for xx in np.arange(-0.15, 7.0, 0.3):
            ax.plot([xx, xx - 0.12], [0, -0.15], color=C_GREY, lw=0.5)
        for cx in (0.9, 4.6, 5.5):
            ax.add_patch(plt.Circle((cx, 0.30), 0.30, fc="white", ec="k", lw=1.2))
        ax.add_patch(Rectangle((0.35, 0.60), 6.0, 0.17, fc="#D0D0D0", ec="k", lw=1.1))
        ax.text(0.15, 0.68, r"$M_{ch}$", ha="right", va="center", fontsize=8.5)
        col = C_RIGID if mode == "rigid" else C_COMPLIANT
        ax.add_patch(Rectangle((1.5, 0.79), 4.2, 0.20, fc="#B8B8B8", ec="k", lw=1.1))
        ax.text(6.0, 0.89, r"$M_{kuz}$", ha="left", va="center", fontsize=8.5)
        ax.add_patch(Rectangle((1.9, 1.01), 3.4, 0.26, fc=col, alpha=0.22,
                               ec=col, lw=1.3))
        ax.text(3.6, 1.14, r"$M_{v}$", ha="center", va="center", fontsize=8.5)
        if mode == "rigid":
            for xx in (1.7, 5.5):
                ax.plot([xx, xx], [0.77, 0.79], color="k", lw=3)
            for xx in (2.1, 5.1):
                ax.plot([xx, xx], [0.99, 1.01], color="k", lw=3)
            ax.text(0.0, 1.30, "(a) rigid: all three bodies move as one",
                    fontsize=8, color=C_RIGID)
        else:
            _spring(ax, 0.55, 1.50, 0.885, amp=0.045)
            _spring(ax, 1.00, 1.90, 1.135, amp=0.045)
            for cx in (2.3, 3.3, 4.3, 5.3):
                ax.add_patch(plt.Circle((cx, 0.735), 0.05, fc="white", ec="k", lw=0.7))
            for cx in (2.6, 3.6, 4.6):
                ax.add_patch(plt.Circle((cx, 0.995), 0.05, fc="white", ec="k", lw=0.7))
            ax.text(0.0, 1.30, r"(b) segmented: $C_{CHK}$, $C_{KV}$ on longitudinal guides",
                    fontsize=8, color=C_COMPLIANT)
        ax.add_patch(FancyArrowPatch((0.3, -0.30), (1.6, -0.30), arrowstyle="-|>",
                                     mutation_scale=8, color=C_ACCENT, lw=1.0))
        ax.text(1.7, -0.30, "direction of travel", va="center", fontsize=7,
                color=C_ACCENT)
    save(fig, "pa01_concept", FIG)


def fig02_model():
    fig, ax = plt.subplots(figsize=(W15, 2.3))
    ax.set_xlim(-0.4, 10.0); ax.set_ylim(-1.15, 1.9); ax.axis("off")
    ax.plot([-0.4, 10.0], [-0.45, -0.45], color=C_GREY, lw=1.0)
    ax.add_patch(plt.Circle((0.85, 0.15), 0.5, fc="white", ec="k", lw=1.3))
    ax.plot([0.85, 1.2], [0.15, 0.5], color="k", lw=1.0)
    ax.text(0.85, 0.15, r"$I_{tr}$", ha="center", va="center", fontsize=8.5)
    ax.annotate(r"$\varphi_k$", (1.4, 0.60), fontsize=9)
    ax.annotate("", xy=(2.0, 0.15), xytext=(1.35, 0.15),
                arrowprops=dict(arrowstyle="-|>", lw=1.0))
    ax.text(1.7, 0.32, r"$C_{KCH}$", ha="center", fontsize=7.5)
    for x0, w, lab, col in ((2.1, 2.2, r"$M_{ch}$", C_GREY),
                            (5.2, 1.5, r"$M_{kuz}$", C_COMPLIANT),
                            (7.9, 1.7, r"$M_{v}$", C_COMPLIANT)):
        ax.add_patch(Rectangle((x0, -0.18), w, 0.9, fc=col, alpha=0.16, ec=col, lw=1.2))
        ax.text(x0 + w / 2, 0.27, lab, ha="center", va="center", fontsize=9)
    _spring(ax, 4.30, 5.20, 0.40)
    _spring(ax, 6.70, 7.90, 0.40)
    ax.text(4.75, 0.68, r"$C_{CHK}$", ha="center", fontsize=7.5)
    ax.text(7.30, 0.68, r"$C_{KV}$", ha="center", fontsize=7.5)
    for x, lab in ((3.2, r"$X_{ch}$"), (5.95, r"$u_1$"), (8.75, r"$u_2$")):
        ax.annotate("", xy=(x + 0.5, 1.25), xytext=(x - 0.5, 1.25),
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C_ACCENT))
        ax.text(x, 1.42, lab, ha="center", color=C_ACCENT, fontsize=8.5)
    ax.text(0.0, -0.72, r"$A_1=\alpha_d M_{\Sigma}g\,$sgn$\,\dot X_{ch}$   "
                        r"$A_2=\alpha_{ch}(M_{kuz}{+}M_v)g\,$sgn$\,\dot u_1$   "
                        r"$A_3=\alpha_k M_v g\,$sgn$\,\dot u_2$",
            fontsize=7.2, color=C_GREY)
    ax.text(0.0, -1.02, r"$F_{kch}=C_{KCH}(\varphi_k-X_{ch}/R_k)/R_k$   "
                        r"(no viscous damping in the archived worksheet)", fontsize=7.2)
    save(fig, "pa02_model", FIG)


def fig03_workflow():
    fig, ax = plt.subplots(figsize=(W1 * 1.1, 3.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 14.5); ax.axis("off")
    steps = ["Archived ETAL Mathcad worksheet\n(XMCD)",
             "Programmatic extraction of scalars,\nfunctions and traction table",
             "Reduced-order longitudinal model\n(zero calibrated parameters)",
             "Reproduction of published\nstage-1 and braking results",
             "Forward comparison\n(archived traction characteristic)",
             "Matched-manoeuvre comparison\n(inverse dynamics)",
             "Parameter and trajectory sweeps"]
    y = 13.8
    for i, s in enumerate(steps):
        col = C_ACCENT if i in (4, 5) else C_GREY
        ax.add_patch(Rectangle((0.4, y - 1.3), 9.2, 1.3, fc=col, alpha=0.10,
                               ec=col, lw=0.9))
        ax.text(5.0, y - 0.65, s, ha="center", va="center", fontsize=7.0)
        if i < len(steps) - 1:
            ax.annotate("", xy=(5.0, y - 1.78), xytext=(5.0, y - 1.33),
                        arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C_GREY))
        y -= 1.95
    save(fig, "pa03_workflow", FIG)


def fig04_reproduction():
    d = pd.read_csv(P / "a1_reproduction.csv")
    z = np.load(P / "a1_timeseries.npz")
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.2), constrained_layout=True)
    ax = axes[0]
    p = ArchiveParams()
    ax.plot(V1 / p.Rk, PT1 * p.Rk / 1e3, "o-", ms=3, color=C_ACCENT)
    ax.set_xlabel(r"$\omega$, rad s$^{-1}$"); ax.set_ylabel(r"$M_{tag}$, kN m")
    ax.set_title("(a) archived traction table", loc="left", fontweight="bold", fontsize=7.5)
    ax = axes[1]
    ax.plot(z["t"], z["v"], color=C_ACCENT, label="stage 1")
    ax.plot(z["tb"], z["vb"], color=C_RIGID, label="stage 2")
    ax.set_xlabel("time, s"); ax.set_ylabel(r"$V_{ch}$, m s$^{-1}$"); ax.legend(fontsize=6.5)
    ax.set_title("(b) reproduced manoeuvre", loc="left", fontweight="bold", fontsize=7.5)
    ax = axes[2]
    lbl = [q.split(",")[0] for q in d.quantity]
    col = [C_ACCENT if s.startswith("stage 1") else C_RIGID for s in d.stage]
    ax.barh(np.arange(len(d)), d.rel_error_pct, color=col)
    ax.axvspan(-2, 2, color=C_COMPLIANT, alpha=0.12)
    ax.axvline(0, color=C_GREY, lw=0.8)
    ax.set_yticks(np.arange(len(d))); ax.set_yticklabels(lbl, fontsize=6)
    ax.set_xlabel("error vs published, %")
    ax.set_title("(c) reproduction error", loc="left", fontweight="bold", fontsize=7.5)
    save(fig, "pa04_reproduction", FIG)


def fig05_forward():
    z = np.load(P / "a3_timeseries.npz")
    s = pd.read_csv(P / "a3_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.2), constrained_layout=True)
    for k, col in (("rigid", C_RIGID), ("segmented", C_COMPLIANT)):
        t, v, F, u1, u2 = z[k]
        axes[0].plot(t, F / 1e3, color=col, label=k)
        axes[1].plot(t, v, color=col)
    axes[0].set_ylabel(r"$F_{kch}$, kN"); axes[0].legend(fontsize=6.5)
    axes[0].set_xlim(0, 3)
    axes[1].set_ylabel(r"$V_{ch}$, m s$^{-1}$")
    w = 260
    axes[2].bar(s.m_v_kg - w, s.rigid / 1e3, width=2 * w, color=C_RIGID, label="rigid")
    axes[2].bar(s.m_v_kg + w, s.segmented / 1e3, width=2 * w, color=C_COMPLIANT,
                label="segmented")
    axes[2].set_xticks(s.m_v_kg); axes[2].set_ylabel(r"$F_{kch}^{peak}$, kN")
    axes[2].set_xlabel(r"$M_v$, kg"); axes[2].legend(fontsize=6.5)
    for a, t_ in zip(axes[:2], ("(a)", "(b)")):
        a.set_xlabel("time, s")
    for a, t_ in zip(axes, ("(a)", "(b)", "(c)")):
        a.set_title(t_, loc="left", fontweight="bold")
    save(fig, "pa05_forward", FIG)


def fig06_matched():
    z = np.load(P / "a4_timeseries.npz")
    t = z["t"]
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.2), constrained_layout=True)
    axes[0].plot(t, z["a"], color=C_GREY)
    axes[0].set_ylabel(r"prescribed $a_{ref}$, m s$^{-2}$")
    axes[1].plot(t, z["F_rigid"] / 1e3, color=C_RIGID, label="rigid")
    axes[1].plot(t, z["F_seg"] / 1e3, color=C_COMPLIANT, label="segmented")
    axes[1].set_ylabel(r"required $F$, kN"); axes[1].legend(fontsize=6.5)
    axes[2].plot(t, z["P_rigid"] / 1e3, color=C_RIGID)
    axes[2].plot(t, z["P_seg"] / 1e3, color=C_COMPLIANT)
    axes[2].set_ylabel(r"required $P$, kW")
    for a, t_ in zip(axes, ("(a)", "(b)", "(c)")):
        a.set_xlabel("time, s"); a.set_title(t_, loc="left", fontweight="bold")
    save(fig, "pa06_matched", FIG)


def fig07_effective_mass():
    z = np.load(P / "a4_timeseries.npz")
    fig, axes = plt.subplots(1, 2, figsize=(W15, 2.4), constrained_layout=True)
    ax = axes[0]
    ax.axhline(11.6, color=C_RIGID, lw=1.3, label=r"rigid: $M_\Sigma$ = 11.6 t")
    for tt, me, col, lab in ((z["tf"], z["m_eff_fast"], C_ACCENT, r"$T$ = 3 s"),
                             (z["t"], z["m_eff"], C_COMPLIANT, r"$T$ = 8 s")):
        ok = np.isfinite(me)
        ax.plot(tt[ok], me[ok] / 1e3, color=col, lw=1.2, label=lab)
    ax.set_xlim(0, 10); ax.set_ylim(6, 14)
    ax.set_xlabel("time, s"); ax.set_ylabel("effective accelerated mass, t")
    ax.legend(fontsize=6.5)
    ax.text(0.98, 0.04, r"undefined where $|a_{ref}|<0.05$ m s$^{-2}$",
            transform=ax.transAxes, fontsize=5.6, color=C_GREY, ha="right")
    ax = axes[1]
    ax.plot(z["tf"], z["u1_fast"] * 1e3, color=C_ACCENT, label=r"$u_1$, $T$=3 s")
    ax.plot(z["tf"], z["u2_fast"] * 1e3, color=C_ACCENT, ls="--",
            label=r"$u_2$, $T$=3 s")
    ax.plot(z["t"], z["u1"] * 1e3, color=C_COMPLIANT, label=r"$u_1$, $T$=8 s")
    ax.plot(z["t"], z["u2"] * 1e3, color=C_COMPLIANT, ls="--", label=r"$u_2$, $T$=8 s")
    ax.set_xlim(0, 20)
    ax.set_xlabel("time, s"); ax.set_ylabel("relative travel, mm")
    ax.legend(fontsize=6, ncol=2)
    for a, t_ in zip(axes, ("(a)", "(b)")):
        a.set_title(t_, loc="left", fontweight="bold")
    save(fig, "pa07_effective_mass", FIG)


def fig08_no_benefit():
    d = pd.read_csv(P / "a5_duration.csv")
    p6 = pd.read_csv(P / "a6_parametric.csv")
    p7 = pd.read_csv(P / "a7_trajectory.csv")
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.3), constrained_layout=True)
    ax = axes[0]
    ax.axhline(0, color=C_GREY, lw=0.9)
    ax.plot(d.T_over_Tn, d.dF_pct, "o-", ms=3, color=C_COMPLIANT, label="peak force")
    ax.plot(d.T_over_Tn, d.dP_pct, "s--", ms=3, color=C_RIGID, label="peak power")
    ax.set_xscale("log"); ax.set_xlabel(r"$T/T_n$")
    ax.set_ylabel("demand reduction, %"); ax.legend(fontsize=6.5)
    ax = axes[1]
    ax.axhline(0, color=C_GREY, lw=0.9)
    marks = {"cargo mass": "o", "body mass": "s", "coupling stiffness": "^",
             "guide friction": "d", "viscous damping ratio": "v"}
    cols = {"cargo mass": C_RIGID, "body mass": C_COMPLIANT,
            "coupling stiffness": C_ACCENT, "guide friction": C_GREY,
            "viscous damping ratio": "#8B5FBF"}
    for sw, g in p6.groupby("sweep"):
        x = np.linspace(0, 1, len(g))
        ax.plot(x, g.dF_pct, marker=marks[sw], ms=3.5, lw=0.9, color=cols[sw], label=sw)
    ax.set_xlabel("normalised sweep position"); ax.set_ylabel("peak-force reduction, %")
    ax.legend(fontsize=5.5)
    ax = axes[2]
    ax.axhline(0, color=C_GREY, lw=0.9)
    for sh, g in p7.groupby("shape"):
        ax.plot(g.T_s, g.dF_pct, "o-", ms=3.5, label=sh)
    ax.set_xlabel(r"$T$, s"); ax.set_ylabel("peak-force reduction, %")
    ax.legend(fontsize=6)
    for a, t_ in zip(axes, ("(a)", "(b)", "(c)")):
        a.set_title(t_, loc="left", fontweight="bold")
    save(fig, "pa08_no_benefit", FIG)


def main():
    use_style()
    FIG.mkdir(parents=True, exist_ok=True)
    fig01_concept(); fig02_model(); fig03_workflow(); fig04_reproduction()
    fig05_forward(); fig06_matched(); fig07_effective_mass(); fig08_no_benefit()
    print("->", FIG)


if __name__ == "__main__":
    main()
