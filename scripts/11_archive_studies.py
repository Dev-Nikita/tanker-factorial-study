"""Manuscript studies on the archive-locked ETAL model (A1-A7).

  A1  reproduction of the published stage-1 and braking results (zero calibration)
  A2  attempted reproduction of the published Table 1 peaks
  A3  forward start-off, rigid vs segmented, archived traction characteristic
  A4  matched-manoeuvre traction demand + effective accelerated mass
  A5  manoeuvre-duration sweep
  A6  parameter sweeps (cargo mass, body mass, stiffness, guide friction, damping)
  A7  reference-trajectory shape sensitivity

Outputs -> results/archive_study/
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.archive_model import ArchiveParams, wheel_force      # noqa: E402
from tanker_dynamics.archive_studies import (Manoeuvre, compare_matched,  # noqa: E402
                                             coupling_period, demand,
                                             forward_metrics, simulate)

OUT = ROOT / "results" / "archive_study"
V_F, T_REF = 5.0, 8.0

PUB_STAGE1 = {"Terminal speed, m/s": 13.058, "Distance travelled, m": 7582.0,
              "Mean speed, m/s": 12.536, "Engine work, J": 5.337e7,
              "Mean power, W": 8.895e4}
PUB_BRAKE = {"Braking time, s": 25.054, "Braking distance, m": 163.575,
             "Mean deceleration, m/s2": -0.521, "Peak wheel force, N": 6800.0}
PUB_TABLE1 = {1250: (3.946e4, 3.853e4), 2500: (4.547e4, 4.163e4),
              3750: (5.027e4, 4.622e4), 5000: (5.375e4, 4.877e4)}


def a1_reproduction():
    p = ArchiveParams()
    t, y, _ = simulate(p, (0.0, 600.0), n=6001, rtol=1e-8, atol=1e-10)
    F = wheel_force(y, p)
    P = F * y[5]
    W = float(np.trapezoid(P, t))
    got = {"Terminal speed, m/s": float(y[5][-1]),
           "Distance travelled, m": float(y[1][-1]),
           "Mean speed, m/s": float(np.trapezoid(y[5], t) / 600.0),
           "Engine work, J": W, "Mean power, W": W / 600.0}

    y0 = y[:, -1].copy()
    tb, yb, solb = simulate(p, (600.0, 800.0), braking=True, y0=y0,
                            stop_at_rest=True, n=6001, rtol=1e-8, atol=1e-10)
    Fb = wheel_force(yb, p)
    ts = float(tb[-1] - 600.0)
    gotb = {"Braking time, s": ts,
            "Braking distance, m": float(yb[1][-1] - y0[1]),
            "Mean deceleration, m/s2": float(-y0[5] / ts),
            "Peak wheel force, N": float(np.max(np.abs(Fb)))}

    rows = []
    for src, g in ((PUB_STAGE1, got), (PUB_BRAKE, gotb)):
        for k, v in src.items():
            rows.append({"stage": "stage 1 (acceleration)" if src is PUB_STAGE1
                         else "stage 2 (braking)",
                         "quantity": k, "published": v, "archive_model": g[k],
                         "rel_error_pct": 100 * (g[k] - v) / v})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "a1_reproduction.csv", index=False)
    np.savez_compressed(OUT / "a1_timeseries.npz", t=t, v=y[5], F=F, P=P,
                        tb=tb, vb=yb[5], Fb=Fb)
    print(df.to_string(index=False))
    print("max |error| = %.2f %%" % df.rel_error_pct.abs().max())
    return df


def a2_table1():
    rows = []
    for Mv, (sr, sc) in PUB_TABLE1.items():
        out = {}
        for rigid, key in ((True, "rigid"), (False, "segmented")):
            p = ArchiveParams(Mv=Mv, rigid=rigid)
            t, y, _ = simulate(p, (0.0, 8.0), n=8001)
            out[key] = float(np.max(np.abs(wheel_force(y, p))))
        rows.append({"m_v_kg": Mv, "published_rigid_N": sr,
                     "model_rigid_N": out["rigid"],
                     "err_rigid_pct": 100 * (out["rigid"] - sr) / sr,
                     "published_compliant_N": sc,
                     "model_segmented_N": out["segmented"],
                     "err_segmented_pct": 100 * (out["segmented"] - sc) / sc,
                     "published_reduction_pct": 100 * (sr - sc) / sr,
                     "model_reduction_pct": 100 * (out["rigid"] - out["segmented"])
                     / out["rigid"]})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "a2_table1_attempt.csv", index=False)
    print(df.to_string(index=False))
    return df


def a3_forward():
    rows, ts = [], {}
    for Mv in (1250, 2500, 3750, 5000):
        for rigid, key in ((True, "rigid"), (False, "segmented")):
            p = ArchiveParams(Mv=Mv, rigid=rigid)
            t, y, _ = simulate(p, (0.0, 60.0), target_speed=V_F, n=6001)
            m = forward_metrics(p, t, y)
            rows.append({"m_v_kg": Mv, "configuration": key, **m})
            if Mv == 5000:
                ts[key] = np.vstack([t, y[5], wheel_force(y, p), y[2], y[3]])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "a3_forward.csv", index=False)
    np.savez_compressed(OUT / "a3_timeseries.npz", **ts)
    piv = df.pivot_table(index="m_v_kg", columns="configuration", values="F_peak_N")
    piv["reduction_pct"] = 100 * (piv["rigid"] - piv["segmented"]) / piv["rigid"]
    piv.to_csv(OUT / "a3_summary.csv")
    print(piv.to_string())
    return df


def a4_matched():
    p = ArchiveParams()
    man = Manoeuvre(v_f=V_F, T=T_REF)
    t = np.linspace(0.0, T_REF + 12.0, 6001)
    r = demand(p, man, t, segmented=False)
    s = demand(p, man, t, segmented=True)
    # a second, faster manoeuvre for the effective-mass figure
    man_f = Manoeuvre(v_f=V_F, T=3.0)
    tf = np.linspace(0.0, 3.0 + 12.0, 6001)
    sf = demand(p, man_f, tf, segmented=True)
    np.savez_compressed(OUT / "a4_timeseries.npz", t=t, a=r["a"], v=r["v"],
                        F_rigid=r["F"], F_seg=s["F"], P_rigid=r["P"], P_seg=s["P"],
                        m_eff=s["m_eff"], u1=s["u1"], u2=s["u2"],
                        tf=tf, m_eff_fast=sf["m_eff"], u1_fast=sf["u1"],
                        u2_fast=sf["u2"], a_fast=sf["a"])
    rows = [("Peak traction force, kN", r["F_peak"] / 1e3, s["F_peak"] / 1e3),
            ("Peak wheel torque, kN m", r["M_peak"] / 1e3, s["M_peak"] / 1e3),
            ("Peak traction power, kW", r["P_peak"] / 1e3, s["P_peak"] / 1e3),
            ("Traction work, MJ", r["W"] / 1e6, s["W"] / 1e6),
            ("Minimum effective mass, t", p.MSchkv / 1e3,
             float(np.nanmin(s["m_eff"])) / 1e3)]
    df = pd.DataFrame(rows, columns=["quantity", "rigid", "segmented"])
    df["change_pct"] = 100 * (df.rigid - df.segmented) / df.rigid
    df.to_csv(OUT / "a4_matched.csv", index=False)
    print(df.to_string(index=False))
    return df


def a5_duration():
    p = ArchiveParams()
    Tn = coupling_period(p)
    rows = []
    for T in [1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 30]:
        c = compare_matched(p, Manoeuvre(v_f=V_F, T=T))
        rows.append({"T_s": T, "T_over_Tn": T / Tn,
                     "F_rigid_kN": c["rigid"]["F_peak"] / 1e3,
                     "F_seg_kN": c["segmented"]["F_peak"] / 1e3,
                     "dF_pct": c["dF_pct"], "dP_pct": c["dP_pct"],
                     "dW_pct": c["dW_pct"], "u1_mm": c["u1_mm"],
                     "m_eff_min_t": c["m_eff_min_kg"] / 1e3})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "a5_duration.csv", index=False)
    print(f"T_n = {Tn:.3f} s; max dF = {df.dF_pct.max():.3f} %; "
          f"worst = {df.dF_pct.min():.2f} %; max dW = {df.dW_pct.max():.4f} %")
    return df


def _row(sweep, value, unit, c, extra=None):
    d = {"sweep": sweep, "value": value, "unit": unit,
         "F_rigid_kN": c["rigid"]["F_peak"] / 1e3,
         "F_seg_kN": c["segmented"]["F_peak"] / 1e3,
         "dF_pct": c["dF_pct"], "dP_pct": c["dP_pct"], "dW_pct": c["dW_pct"],
         "u1_mm": c["u1_mm"], "u2_mm": c["u2_mm"]}
    if extra:
        d.update(extra)
    return d


def a6_parametric():
    man = Manoeuvre(v_f=V_F, T=T_REF)
    rows = []
    for Mv in (1250, 2500, 3750, 5000, 7500):
        rows.append(_row("cargo mass", Mv, "kg",
                         compare_matched(ArchiveParams(Mv=Mv), man)))
    for Mkuz in (350, 700, 1400, 2800):
        rows.append(_row("body mass", Mkuz, "kg",
                         compare_matched(ArchiveParams(Mkuz=Mkuz), man)))
    for C in (2.5e4, 5e4, 1e5, 2e5, 5e5, 1e6, 5e6):
        rows.append(_row("coupling stiffness", C, "N/m",
                         compare_matched(ArchiveParams(CCHK=C, CKV=C), man)))
    for f in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40):
        rows.append(_row("guide friction", f, "-",
                         compare_matched(ArchiveParams(Alf_ch=f, Alf_k=f), man)))
    for zeta in (0.0, 0.1, 0.3, 0.7, 1.0, 1.5):
        p0 = ArchiveParams()
        c1 = 2 * zeta * np.sqrt(p0.CCHK * p0.Mkuz)
        c2 = 2 * zeta * np.sqrt(p0.CKV * p0.Mv)
        rows.append(_row("viscous damping ratio", zeta, "-",
                         compare_matched(ArchiveParams(c1=c1, c2=c2), man)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "a6_parametric.csv", index=False)
    print(df.groupby("sweep")[["dF_pct", "dW_pct"]].agg(["min", "max"]).to_string())
    return df


def a7_trajectory():
    p = ArchiveParams()
    rows = []
    for shape in ("smoothstep", "halfcosine", "quintic"):
        for T in (3.0, 8.0, 20.0):
            c = compare_matched(p, Manoeuvre(v_f=V_F, T=T, shape=shape))
            rows.append({"shape": shape, "T_s": T,
                         "F_rigid_kN": c["rigid"]["F_peak"] / 1e3,
                         "F_seg_kN": c["segmented"]["F_peak"] / 1e3,
                         "dF_pct": c["dF_pct"], "dW_pct": c["dW_pct"]})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "a7_trajectory.csv", index=False)
    print(df.to_string(index=False))
    return df


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, fn in (("A1 reproduction", a1_reproduction), ("A2 Table 1", a2_table1),
                     ("A3 forward", a3_forward), ("A4 matched", a4_matched),
                     ("A5 duration", a5_duration), ("A6 parametric", a6_parametric),
                     ("A7 trajectory", a7_trajectory)):
        print(f"\n== {name} ==", flush=True)
        fn()
    print("\n->", OUT)


if __name__ == "__main__":
    main()
