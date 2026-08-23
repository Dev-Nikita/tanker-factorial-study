"""Generate every manuscript figure from stored result files.

No numerical result is hard-coded here; each panel reads a CSV/NPZ produced by
the preceding pipeline stages.
"""
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
from tanker_dynamics import regression as rg                                  # noqa: E402

TDIR = ROOT / "results" / "tables"
FDIR = ROOT / "results" / "figures"
GDIR = ROOT / "data" / "generated"
SCALE = "full_scale"


# --------------------------------------------------------------------------- #
def _spring(ax, x0, x1, y, n=7, amp=0.045, lw=1.0, color="k"):
    x = np.linspace(x0, x1, 200)
    s = np.zeros_like(x)
    m = (x > x0 + 0.15 * (x1 - x0)) & (x < x1 - 0.15 * (x1 - x0))
    tt = (x[m] - x[m][0]) / (x[m][-1] - x[m][0])
    s[m] = amp * np.sin(2 * np.pi * n * tt)
    ax.plot(x, y + s, color=color, lw=lw)


def _damper(ax, x0, x1, y, h=0.05, lw=1.0, color="k"):
    xm = 0.5 * (x0 + x1)
    ax.plot([x0, xm - 0.02 * (x1 - x0)], [y, y], color=color, lw=lw)
    ax.add_patch(Rectangle((xm - 0.10 * (x1 - x0), y - h), 0.16 * (x1 - x0), 2 * h,
                           fill=False, ec=color, lw=lw))
    ax.plot([xm + 0.02 * (x1 - x0), xm + 0.02 * (x1 - x0)], [y - h, y + h], color=color, lw=lw)
    ax.plot([xm + 0.02 * (x1 - x0), x1], [y, y], color=color, lw=lw)


def fig01_concept():
    """Fig. 1 - rigid vs compliant-damped cargo mounting."""
    fig, axes = plt.subplots(2, 1, figsize=(W15, 3.5))
    for ax, mode in zip(axes, ("rigid", "compliant")):
        ax.set_xlim(-0.2, 6.6); ax.set_ylim(-0.35, 1.25); ax.axis("off")
        ax.plot([-0.2, 6.6], [0, 0], color=C_GREY, lw=1.2)
        for xx in np.arange(-0.15, 6.6, 0.28):
            ax.plot([xx, xx - 0.12], [0, -0.14], color=C_GREY, lw=0.5)
        # tractor
        ax.add_patch(Rectangle((0.15, 0.30), 1.05, 0.55, fc="#E8E8E8", ec="k", lw=1.0))
        # chassis / frame
        ax.add_patch(Rectangle((1.25, 0.30), 5.1, 0.16, fc="#D9D9D9", ec="k", lw=1.0))
        for cx in (0.55, 1.05, 2.1, 5.6, 6.05):
            ax.add_patch(plt.Circle((cx, 0.16), 0.16, fc="white", ec="k", lw=1.0))
        # tanks
        col = C_RIGID if mode == "rigid" else C_COMPLIANT
        for i, x0 in enumerate((2.55, 4.45)):
            ax.add_patch(Rectangle((x0, 0.48), 1.55, 0.52, fc=col, alpha=0.18, ec=col, lw=1.2))
            ax.text(x0 + 0.775, 0.74, f"$m_{{v{2-i}}}$" if i == 0 else f"$m_{{v{2-i}}}$",
                    ha="center", va="center", fontsize=8)
        if mode == "rigid":
            for x0 in (2.55, 4.45):
                for xx in (x0 + 0.15, x0 + 1.40):
                    ax.plot([xx, xx], [0.46, 0.48], color="k", lw=2.5)
            ax.text(3.3, 1.10, "rigid longitudinal attachment", fontsize=8, color=C_RIGID)
        else:
            _spring(ax, 6.00, 6.35, 0.86)
            _spring(ax, 4.10, 4.45, 0.86)
            _damper(ax, 6.00, 6.35, 0.62)
            _damper(ax, 4.10, 4.45, 0.62)
            ax.text(3.0, 1.10, "series compliant-damped coupling "
                               r"($C_{chv1},\alpha_{v1}$; $C_{v1v2},\alpha_{v2}$)",
                    fontsize=8, color=C_COMPLIANT)
            for x0 in (2.55, 4.45):
                ax.add_patch(plt.Circle((x0 + 0.25, 0.42), 0.055, fc="white", ec="k", lw=0.8))
                ax.add_patch(plt.Circle((x0 + 1.30, 0.42), 0.055, fc="white", ec="k", lw=0.8))
        ax.add_patch(FancyArrowPatch((0.2, -0.26), (1.5, -0.26), arrowstyle="-|>",
                                     mutation_scale=8, color=C_ACCENT, lw=1.0))
        ax.text(1.6, -0.26, "direction of travel $X_0$", va="center", fontsize=7, color=C_ACCENT)
        ax.text(-0.15, 1.12, "(a)" if mode == "rigid" else "(b)", fontweight="bold", fontsize=9)
    save(fig, "fig01_system_concept", FDIR)


def fig02_model():
    """Fig. 2 - lumped-parameter model with coordinates."""
    fig, ax = plt.subplots(figsize=(W15, 2.5))
    ax.set_xlim(-0.4, 9.4); ax.set_ylim(-1.0, 2.0); ax.axis("off")
    ax.plot([-0.4, 9.4], [-0.35, -0.35], color=C_GREY, lw=1.0)
    ax.add_patch(plt.Circle((0.9, 0.15), 0.5, fc="white", ec="k", lw=1.3))
    ax.plot([0.9, 1.25], [0.15, 0.5], color="k", lw=1.0)
    ax.text(0.9, 0.15, r"$I_k$", ha="center", va="center")
    ax.annotate(r"$\varphi_k$", (1.45, 0.62), fontsize=9)
    ax.annotate("", xy=(2.1, 0.15), xytext=(1.4, 0.15),
                arrowprops=dict(arrowstyle="-|>", lw=1.0))
    ax.text(1.75, 0.30, r"$C_{kch},\,D_{kch}$", ha="center", fontsize=7.5)
    boxes = [(2.2, r"$m_{ch}+m_{tr}$", C_GREY), (5.0, r"$m_{v1}$", C_COMPLIANT),
             (7.5, r"$m_{v2}$", C_COMPLIANT)]
    for x0, lab, col in boxes:
        w = 2.1 if x0 == 2.2 else 1.5
        ax.add_patch(Rectangle((x0, -0.15), w, 0.85, fc=col, alpha=0.15, ec=col, lw=1.2))
        ax.text(x0 + w / 2, 0.28, lab, ha="center", va="center", fontsize=9)
    _spring(ax, 4.30, 5.00, 0.55); _damper(ax, 4.30, 5.00, 0.05)
    _spring(ax, 6.50, 7.50, 0.55); _damper(ax, 6.50, 7.50, 0.05)
    ax.text(4.65, 0.80, r"$C_{chv1},\ \alpha_{v1}$", ha="center", fontsize=7.5)
    ax.text(7.00, 0.80, r"$C_{v1v2},\ \alpha_{v2}$", ha="center", fontsize=7.5)
    for x, lab in ((3.25, r"$X_{ch}$"), (5.75, r"$X_{v1}$"), (8.25, r"$X_{v2}$")):
        ax.annotate("", xy=(x + 0.55, 1.35), xytext=(x - 0.55, 1.35),
                    arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C_ACCENT))
        ax.text(x, 1.48, lab, ha="center", color=C_ACCENT, fontsize=8.5)
    ax.annotate("", xy=(2.0, -0.65), xytext=(0.4, -0.65),
                arrowprops=dict(arrowstyle="-|>", lw=1.0, color=C_ACCENT))
    ax.text(2.1, -0.62, r"$F_{kch}=\left(C_{kch}\theta+D_{kch}\dot\theta\right)/r_k$,"
                        r"  $\theta=\varphi_k-X_{ch}/r_k$", va="center", fontsize=7.2)
    ax.text(2.1, -0.92, r"road resistance $f_d m g$;  guide resistance $f_{ch},\,f_k$",
            va="center", fontsize=7.2, color=C_GREY)
    save(fig, "fig02_lumped_model", FDIR)


def fig03_workflow():
    """Fig. 3 - computational workflow."""
    fig, ax = plt.subplots(figsize=(W1 * 1.15, 4.6))
    ax.set_xlim(0, 10); ax.set_ylim(0, 22); ax.axis("off")
    steps = ["Source model\n(Lagrange eq. (1)-(2))",
             "Parameter identification\n(4 unreported parameters)",
             "Python / Radau implementation",
             "Reproduction of SRC-A\nTable 1 and stage-1 integrals",
             r"$3^3$ full-factorial design (27 runs)",
             "Rigid and compliant runs\nstart-off + braking",
             "Second-order response surfaces\n(10 coefficients, 17 d.o.f.)",
             "ANOVA, curvature and\ninteraction assessment",
             "Uncertainty propagation\n(Monte Carlo)",
             "Physical-test protocol\n(81-run scaled rig)"]
    y = 21.0
    for i, s in enumerate(steps):
        col = C_ACCENT if i in (4, 6) else C_GREY
        ax.add_patch(Rectangle((0.4, y - 1.35), 9.2, 1.35, fc=col, alpha=0.10,
                               ec=col, lw=0.9))
        ax.text(5.0, y - 0.68, s, ha="center", va="center", fontsize=7.2)
        if i < len(steps) - 1:
            ax.annotate("", xy=(5.0, y - 1.85), xytext=(5.0, y - 1.38),
                        arrowprops=dict(arrowstyle="-|>", lw=0.9, color=C_GREY))
        y -= 2.1
    save(fig, "fig03_workflow", FDIR)


def fig04_design():
    """Fig. 4 - the 3^3 design cube."""
    d = pd.read_csv(TDIR / "factorial_design.csv")
    fig = plt.figure(figsize=(W1, 3.0))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(d.x1_coded, d.x2_coded, d.x3_coded, s=26, c=d.total_movable_mass_kg,
               cmap="viridis", depthshade=False, edgecolor="k", linewidth=0.3)
    ax.set_xlabel(r"$x_1$ front fill"); ax.set_ylabel(r"$x_2$ rear fill")
    ax.set_zlabel(r"$x_3$ speed")
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_ticks([-1, 0, 1])
    ax.view_init(22, 38)
    ax.set_box_aspect((1, 1, 0.9))
    save(fig, "fig04_design_cube", FDIR)


def fig05_reproduction():
    """Fig. 5 - reproduction of the source Table 1."""
    r = pd.read_csv(TDIR / "original_reproduction.csv")
    t1 = r[r.block.str.startswith("Table 1 - peak")]
    m = [1250, 2500, 3750, 5000]
    fig, axes = plt.subplots(1, 2, figsize=(W15, 2.4))
    ax = axes[0]
    w = 260
    for k, (cfgname, col) in enumerate((("rigid", C_RIGID), ("compliant", C_COMPLIANT))):
        sub = t1[t1.quantity.str.contains(cfgname)]
        ax.bar(np.array(m) + (k - 0.5) * 2 * w - w / 2, sub.source_value / 1e4, width=w,
               color=col, alpha=0.45, label=f"source, {cfgname}")
        ax.bar(np.array(m) + (k - 0.5) * 2 * w + w / 2, sub.reproduced / 1e4, width=w,
               color=col, label=f"reproduced, {cfgname}")
    ax.set_xlabel(r"cargo mass per unit $m_v$, kg")
    ax.set_ylabel(r"$F_{kch}^{peak}\times10^{-4}$, N")
    ax.set_xticks(m); ax.legend(ncol=2, fontsize=6)
    ax = axes[1]
    ax.axhline(0, color=C_GREY, lw=0.8)
    ax.bar(range(len(t1)), t1.rel_error_pct, color=[C_RIGID if "rigid" in q else C_COMPLIANT
                                                    for q in t1.quantity])
    ax.axhspan(-5, 5, color=C_ACCENT, alpha=0.12)
    ax.set_xticks(range(len(t1)))
    ax.set_xticklabels([q.replace("m_v = ", "").replace(" kg, ", "\n") for q in t1.quantity],
                       fontsize=5.5, rotation=90)
    ax.set_ylabel("relative error, %")
    save(fig, "fig05_reproduction", FDIR)


def _main_effects(ax, df, y, scale=1.0, ylabel=""):
    labels = {"x1_coded": r"$x_1$ front fill", "x2_coded": r"$x_2$ rear fill",
              "x3_coded": r"$x_3$ speed"}
    for c, col, mk in zip(labels, (C_RIGID, C_COMPLIANT, C_ACCENT), ("o", "s", "^")):
        g = df.groupby(c)[y].mean() * scale
        ax.plot(g.index, g.values, marker=mk, ms=4, color=col, label=labels[c])
    ax.set_xticks([-1, 0, 1]); ax.set_xlabel("coded factor level")
    ax.set_ylabel(ylabel)


def fig06_main_effects(df):
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.1), constrained_layout=True)
    _main_effects(axes[0], df, "P_peak_W", 1e-3, r"$P_e^{peak}$, kW")
    _main_effects(axes[1], df, "P_mean_W", 1e-3, r"$\bar P_e$, kW")
    _main_effects(axes[2], df, "E_engine_J", 1e-6, r"$E_e$, MJ")
    axes[0].legend(loc="best")
    for a, t in zip(axes, ("(a)", "(b)", "(c)")):
        a.set_title(t, loc="left", fontweight="bold")
    save(fig, "fig06_main_effects_power", FDIR)


def fig07_interactions(df):
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.1), constrained_layout=True)
    pairs = [("x1_coded", "x2_coded"), ("x1_coded", "x3_coded"), ("x2_coded", "x3_coded")]
    names = {"x1_coded": r"$x_1$", "x2_coded": r"$x_2$", "x3_coded": r"$x_3$"}
    for ax, (a, b) in zip(axes, pairs):
        for lv, col in zip((-1, 0, 1), (C_RIGID, C_GREY, C_COMPLIANT)):
            g = df[df[b] == lv].groupby(a)["P_peak_W"].mean() / 1e3
            ax.plot(g.index, g.values, marker="o", ms=3.5, color=col,
                    label=f"{names[b]} = {lv:+d}")
        ax.set_xticks([-1, 0, 1]); ax.set_xlabel(names[a])
        ax.set_ylabel(r"$P_e^{peak}$, kW"); ax.legend(fontsize=6)
    save(fig, "fig07_interactions_power", FDIR)


def fig08_surfaces(df):
    """Response surface x1 x x2 at each speed level, for peak wheel force."""
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.2), sharey=True, constrained_layout=True)
    res = rg.fit(df, "F_kch_peak_N")
    g = np.linspace(-1, 1, 60)
    X1, X2 = np.meshgrid(g, g)
    for ax, lv in zip(axes, (-1, 0, 1)):
        Z = (res.params["const"] + res.params["x1"] * X1 + res.params["x2"] * X2
             + res.params["x3"] * lv + res.params["x1x2"] * X1 * X2
             + res.params["x1x3"] * X1 * lv + res.params["x2x3"] * X2 * lv
             + res.params["x1^2"] * X1**2 + res.params["x2^2"] * X2**2
             + res.params["x3^2"] * lv**2) / 1e4
        cs = ax.contourf(X1, X2, Z, levels=14, cmap="viridis")
        ax.contour(X1, X2, Z, levels=14, colors="k", linewidths=0.3)
        sub = df[df.x3_coded == lv]
        ax.scatter(sub.x1_coded, sub.x2_coded, s=14, c="white", ec="k", lw=0.5, zorder=5)
        ax.set_xlabel(r"$x_1$"); ax.set_title(rf"$x_3={lv:+d}$", fontsize=8)
        fig.colorbar(cs, ax=ax, fraction=0.046, pad=0.03)
    axes[0].set_ylabel(r"$x_2$")
    fig.suptitle(r"$F_{kch}^{peak}\times10^{-4}$, N", fontsize=8, y=1.02)
    save(fig, "fig08_surface_force", FDIR)


def fig09_loads(paired):
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.3), constrained_layout=True)
    ax = axes[0]
    ax.plot(paired.run_id, paired["F_kch_peak_N__rigid"] / 1e3, "o-", ms=3,
            color=C_RIGID, label="rigid")
    ax.plot(paired.run_id, paired["F_kch_peak_N__compliant"] / 1e3, "s-", ms=3,
            color=C_COMPLIANT, label="compliant")
    ax.set_ylabel(r"$F_{kch}^{peak}$, kN"); ax.legend()
    axes[1].plot(paired.run_id, paired["F_kch_peak_N__reduction_pct"], "^-", ms=3,
                 color=C_ACCENT)
    axes[1].set_ylabel("peak-load reduction, %")
    axes[2].plot(paired.run_id, paired["F_chv1_peak_N__compliant"] / 1e3, "s-", ms=3,
                 color=C_COMPLIANT, label=r"$F_{chv1}$")
    axes[2].plot(paired.run_id, paired["F_v1v2_peak_N__compliant"] / 1e3, "d-", ms=3,
                 color=C_ACCENT, label=r"$F_{v1v2}$")
    axes[2].set_ylabel("coupling force, kN"); axes[2].legend()
    for a, t in zip(axes, ("(a)", "(b)", "(c)")):
        a.set_xlabel("design point"); a.set_title(t, loc="left", fontweight="bold")
    save(fig, "fig09_loads_rigid_vs_compliant", FDIR)


def fig10_accel_jerk(df):
    fig, axes = plt.subplots(1, 2, figsize=(W15, 2.1), constrained_layout=True)
    _main_effects(axes[0], df, "a_peak_m_s2", 1.0, r"$a_x^{peak}$, m s$^{-2}$")
    _main_effects(axes[1], df, "j_peak_m_s3", 1.0, r"$j_x^{peak}$, m s$^{-3}$")
    axes[0].legend()
    save(fig, "fig10_acceleration_jerk", FDIR)


def fig11_relative(df):
    fig, axes = plt.subplots(1, 2, figsize=(W15, 2.1), constrained_layout=True)
    _main_effects(axes[0], df, "u1_peak_m", 1e3, r"$|u_1|^{peak}$, mm")
    _main_effects(axes[1], df, "u2_peak_m", 1e3, r"$|u_2|^{peak}$, mm")
    axes[0].legend()
    save(fig, "fig11_relative_displacement", FDIR)


def fig12_timeseries():
    p = GDIR / f"timeseries_{SCALE}.npz"
    if not p.exists():
        return
    z = np.load(p)
    keys = sorted({k.split("_", 1)[1] for k in z.files}, key=int)
    k = keys[len(keys) // 2]
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.1), constrained_layout=True)
    for cfgname, col in (("rigid", C_RIGID), ("compliant", C_COMPLIANT)):
        a = z[f"{cfgname}_{k}"]
        t, v, F, acc, jerk, u1, u2 = a
        axes[0].plot(t, F / 1e3, color=col, label=cfgname)
        axes[1].plot(t, acc, color=col)
        axes[2].plot(t, u1 * 1e3, color=col)
        if cfgname == "compliant":
            axes[2].plot(t, u2 * 1e3, color=C_ACCENT, ls="--", label=r"$u_2$")
    axes[0].set_ylabel(r"$F_{kch}$, kN"); axes[1].set_ylabel(r"$a_x$, m s$^{-2}$")
    axes[2].set_ylabel(r"$u_1$, mm")
    for a in axes:
        a.set_xlabel("time, s")
    axes[0].legend(); axes[2].legend()
    fig.suptitle(f"design point {k}", fontsize=8, y=1.03)
    save(fig, "fig12_time_histories", FDIR)


def fig13_parity(df):
    fig, ax = plt.subplots(figsize=(W1, 2.6))
    for y, col, mk in (("P_peak_W", C_RIGID, "o"), ("F_kch_peak_N", C_COMPLIANT, "s"),
                       ("a_peak_m_s2", C_ACCENT, "^")):
        res = rg.fit(df, y)
        obs, pred = df[y].to_numpy(float), res.fittedvalues
        n = (obs - obs.mean()) / obs.std(ddof=1)
        pn = (pred - obs.mean()) / obs.std(ddof=1)
        ax.scatter(n, pn, s=16, color=col, marker=mk, alpha=0.8,
                   label=f"{y} ($R^2$={res.rsquared:.3f})")
    lim = [-2.6, 2.6]
    ax.plot(lim, lim, color=C_GREY, lw=0.8, ls="--")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("standardised simulated"); ax.set_ylabel("standardised predicted")
    ax.legend(fontsize=6)
    save(fig, "fig13_parity", FDIR)


def fig14_residuals(df):
    from scipy import stats
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.1), constrained_layout=True)
    res = rg.fit(df, "F_kch_peak_N")
    r = res.resid
    axes[0].scatter(res.fittedvalues / 1e3, r, s=14, color=C_COMPLIANT)
    axes[0].axhline(0, color=C_GREY, lw=0.8)
    axes[0].set_xlabel(r"fitted $F_{kch}^{peak}$, kN"); axes[0].set_ylabel("residual, N")
    stats.probplot(r, plot=axes[1])
    axes[1].get_lines()[0].set_color(C_COMPLIANT); axes[1].get_lines()[0].set_markersize(3)
    axes[1].get_lines()[1].set_color(C_GREY)
    axes[1].set_title(""); axes[1].set_ylabel("ordered residual")
    infl = res.get_influence().cooks_distance[0]
    axes[2].stem(np.arange(1, len(infl) + 1), infl, basefmt=" ")
    axes[2].set_xlabel("design point"); axes[2].set_ylabel("Cook's distance")
    save(fig, "fig14_residual_diagnostics", FDIR)


def fig15_effects():
    frames = []
    for cfgname in ("rigid", "compliant"):
        f = TDIR / f"effect_ranking_{SCALE}_{cfgname}.csv"
        if f.exists():
            d = pd.read_csv(f); d["configuration"] = cfgname
            frames.append(d)
    if not frames:
        return
    d = pd.concat(frames)
    show = ["P_peak_W", "F_kch_peak_N", "a_peak_m_s2", "j_peak_m_s3", "u1_peak_m"]
    d = d[(d.configuration == "compliant") & (d.response.isin(show))]
    piv = d.pivot(index="term", columns="response", values="std_coefficient")
    piv = piv.reindex(["x1", "x2", "x3", "x1x2", "x1x3", "x2x3", "x1^2", "x2^2", "x3^2"])
    fig, ax = plt.subplots(figsize=(W15, 2.6))
    im = ax.imshow(piv.to_numpy(float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns, rotation=25,
                                                               ha="right", fontsize=6.5)
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index, fontsize=7)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.iloc[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5,
                        color="white" if abs(v) > 0.55 else "black")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="standardised coefficient")
    ax.grid(False)
    save(fig, "fig15_standardised_effects", FDIR)


def fig16_uncertainty():
    f = TDIR / "uncertainty_summary.csv"
    if not f.exists():
        return
    mc = pd.read_csv(ROOT / "results" / "tables" / "uncertainty_samples.csv")
    fig, axes = plt.subplots(1, 2, figsize=(W15, 2.2), constrained_layout=True)
    axes[0].violinplot([mc[mc.configuration == c].F_kch_peak_N / 1e3
                        for c in ("rigid", "compliant")], showmedians=True)
    axes[0].set_xticks([1, 2]); axes[0].set_xticklabels(["rigid", "compliant"])
    axes[0].set_ylabel(r"$F_{kch}^{peak}$, kN")
    red = mc.pivot_table(index="sample", columns="configuration", values="F_kch_peak_N")
    r = 100 * (red["rigid"] - red["compliant"]) / red["rigid"]
    axes[1].hist(r, bins=28, color=C_COMPLIANT, alpha=0.85)
    axes[1].axvline(r.mean(), color=C_RIGID, lw=1.2,
                    label=f"mean {r.mean():.2f} %")
    axes[1].set_xlabel("peak-load reduction, %"); axes[1].set_ylabel("count")
    axes[1].legend()
    save(fig, "fig16_uncertainty", FDIR)


def main():
    use_style()
    fig01_concept(); fig02_model(); fig03_workflow(); fig04_design()
    if (TDIR / "original_reproduction.csv").exists():
        fig05_reproduction()
    design = pd.read_csv(TDIR / f"factorial_design_{SCALE}.csv")
    resp = pd.read_csv(TDIR / f"factorial_responses_{SCALE}.csv")
    resp = resp[resp.scenario == "start_off"]
    df = design.merge(resp[resp.configuration == "compliant"], on="run_id")
    paired = pd.read_csv(TDIR / f"factorial_paired_{SCALE}.csv")
    fig06_main_effects(df); fig07_interactions(df); fig08_surfaces(df)
    fig09_loads(paired); fig10_accel_jerk(df); fig11_relative(df)
    fig12_timeseries(); fig13_parity(df); fig14_residuals(df)
    fig15_effects(); fig16_uncertainty()


if __name__ == "__main__":
    main()
