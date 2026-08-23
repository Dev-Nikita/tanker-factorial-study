"""All computational studies of the manuscript (S1-S8).

This is the paper pipeline. The 3^3 factorial scripts (02-08) belong to a separate,
future experimental study and are NOT part of it.

  S1  reproduction of the source results          -> reuses 01_reproduce_original.py
  S2  damping-term ablation (literal vs augmented source model)
  S3  forward start-off, rigid vs segmented, matched engine
  S4  inverse dynamics: matched manoeuvre, traction/torque/power demand
  S5  effective accelerated mass
  S6  manoeuvre-duration sweep -> the regime boundary
  S7  parametric sweeps: cargo mass, stiffness, damping, guide resistance
  S8  asymmetric loading

Outputs go to results/paper/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.config import load_all, build_params            # noqa: E402
from tanker_dynamics.solver import simulate, responses               # noqa: E402
from tanker_dynamics.inverse import (StartOff, compare,              # noqa: E402
                                     demand_rigid, demand_segmented)

OUT = ROOT / "results" / "paper"
V_F = 5.0            # target start-off speed, m/s (upper level of the source brief)
T_REF = 8.0          # reference manoeuvre duration, s
TAIL = 10.0


def _cfg(**over):
    c = load_all()
    for path, val in over.items():
        sect, key = path.split("__")
        grp, sub = sect.split("/")
        c[grp][sub][key] = val
    return c


def coupling_period(p) -> float:
    return 2.0 * np.pi * np.sqrt(p.m_unit_kg[0] / p.C_chv1_N_per_m)


# --------------------------------------------------------------------------- #
def s2_damping_ablation() -> pd.DataFrame:
    """Literal source model (D_kch = 0) versus the augmented dissipative model."""
    targets = {1250: (3.946e4, 3.853e4), 2500: (4.547e4, 4.163e4),
               3750: (5.027e4, 4.622e4), 5000: (5.375e4, 4.877e4)}
    rows = []
    for label, D in (("A: literal source (D_kch = 0)", 0.0),
                     ("B: augmented (D_kch identified)", None)):
        c = load_all()
        if D is not None:
            c["vehicle"]["wheels"]["D_kch_Nms_per_rad"] = D
        full = c["tanks"]["full_scale"]["m_liquid_full_kg"]
        for m_v, (src_r, src_c) in targets.items():
            p = build_params(c, m_v / full, m_v / full)
            fr = responses(simulate(p, "rigid", (0.0, 8.0), n_out=8001))["F_kch_peak_N"]
            fc = responses(simulate(p, "compliant", (0.0, 8.0), n_out=8001))["F_kch_peak_N"]
            rows.append({"model": label, "m_v_kg": m_v,
                         "F_rigid_source_N": src_r, "F_rigid_model_N": fr,
                         "err_rigid_pct": 100 * (fr - src_r) / src_r,
                         "F_compliant_source_N": src_c, "F_compliant_model_N": fc,
                         "err_compliant_pct": 100 * (fc - src_c) / src_c})
    df = pd.DataFrame(rows)
    agg = (df.groupby("model")[["err_rigid_pct", "err_compliant_pct"]]
             .apply(lambda g: np.sqrt(np.mean(np.square(g.to_numpy())))).rename("rms_error_pct"))
    df.to_csv(OUT / "s2_damping_ablation.csv", index=False)
    agg.to_csv(OUT / "s2_damping_ablation_summary.csv")
    print(agg.to_string())
    return df


def s3_forward() -> pd.DataFrame:
    """Forward start-off with the source engine characteristic, both configurations."""
    rows, ts = [], {}
    c = load_all()
    for lam in (0.25, 0.50, 0.75, 1.00):
        p = build_params(c, lam, lam)
        for cfg_name in ("rigid", "compliant"):
            sol = simulate(p, cfg_name, (0.0, 120.0), target_speed=V_F, n_out=4001)
            r = responses(sol)
            rows.append({"lambda": lam, "cargo_per_tank_kg": lam * 5000,
                         "configuration": cfg_name, **{k: v for k, v in r.items()}})
            if lam == 1.00:
                u1, u2, _, _ = sol.relative()
                ts[cfg_name] = np.vstack([sol.t, sol.v_ch, sol.wheel_force(),
                                          sol.accel(), u1, u2])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s3_forward_startoff.csv", index=False)
    np.savez_compressed(OUT / "s3_timeseries.npz", **ts)
    piv = df.pivot_table(index="cargo_per_tank_kg", columns="configuration",
                         values="F_kch_peak_N")
    piv["reduction_pct"] = 100 * (piv["rigid"] - piv["compliant"]) / piv["rigid"]
    piv.to_csv(OUT / "s3_summary.csv")
    print(piv.to_string())
    return df


def s4_inverse() -> dict:
    """Matched-manoeuvre traction demand at the reference duration."""
    p = build_params(load_all(), 1.0, 1.0)
    man = StartOff(v_f=V_F, T=T_REF)
    t = np.linspace(0.0, T_REF + TAIL, 4001)
    r, s = demand_rigid(p, man, t), demand_segmented(p, man, t)
    np.savez_compressed(OUT / "s4_inverse_timeseries.npz",
                        t=t, v=r["v"], a=r["a"],
                        F_rigid=r["F"], F_seg=s["F"],
                        P_rigid=r["P"], P_seg=s["P"],
                        m_eff_rigid=r["m_eff"], m_eff_seg=s["m_eff"],
                        u1=s["u1"], u2=s["u2"], F1=s["F1"])
    rows = [{"quantity": "peak traction force, kN", "rigid": r["F_peak"] / 1e3,
             "segmented": s["F_peak"] / 1e3},
            {"quantity": "peak wheel torque, kN m", "rigid": r["M_peak"] / 1e3,
             "segmented": s["M_peak"] / 1e3},
            {"quantity": "peak power demand, kW", "rigid": r["P_peak"] / 1e3,
             "segmented": s["P_peak"] / 1e3},
            {"quantity": "energy demand, MJ", "rigid": r["W"] / 1e6,
             "segmented": s["W"] / 1e6},
            {"quantity": "min effective mass, t", "rigid": p.m_total_kg / 1e3,
             "segmented": np.nanmin(s["m_eff"]) / 1e3}]
    df = pd.DataFrame(rows)
    df["change_pct"] = 100 * (df.rigid - df.segmented) / df.rigid
    df.to_csv(OUT / "s4_inverse_summary.csv", index=False)
    print(df.to_string(index=False))
    return {"rigid": r, "segmented": s, "params": p}


def s6_duration_sweep() -> pd.DataFrame:
    """The regime boundary: benefit versus manoeuvre duration."""
    p = build_params(load_all(), 1.0, 1.0)
    Tn = coupling_period(p)
    rows = []
    for T in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0,
              10.0, 12.0, 15.0, 20.0, 25.0, 30.0]:
        c = compare(p, StartOff(v_f=V_F, T=T), t_tail=TAIL)
        rows.append({"T_s": T, "T_over_Tn": T / Tn, "a_max_m_s2": 1.5 * V_F / T,
                     "F_rigid_kN": c["rigid"]["F_peak"] / 1e3,
                     "F_seg_kN": c["segmented"]["F_peak"] / 1e3,
                     "dF_peak_pct": c["dF_peak_pct"],
                     "dP_peak_pct": c["dP_peak_pct"],
                     "dW_pct": c["dW_pct"],
                     "u1_peak_mm": c["u1_peak_mm"], "u2_peak_mm": c["u2_peak_mm"]})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s6_duration_sweep.csv", index=False)
    pos = df[df.dF_peak_pct > 0]
    print(f"T_n = {Tn:.3f} s;  benefit appears above T/T_n ~ "
          f"{pos.T_over_Tn.min():.2f};  max benefit {df.dF_peak_pct.max():.2f} %")
    return df


def s7_parametric() -> pd.DataFrame:
    """Cargo mass, coupling stiffness, damping and guide resistance."""
    rows = []
    base = load_all()
    p0 = build_params(base, 1.0, 1.0)

    for lam in (0.25, 0.50, 0.75, 1.00):
        p = build_params(base, lam, lam)
        c = compare(p, StartOff(v_f=V_F, T=T_REF), t_tail=TAIL)
        rows.append({"sweep": "cargo mass", "value": lam * 5000, "unit": "kg/tank",
                     "T_over_Tn": T_REF / coupling_period(p), **_metrics(c)})

    for C in (2.5e4, 5e4, 1e5, 2e5, 5e5, 1e6):
        cfg = _cfg(**{"springs/compliant__C_chv1_N_per_m": C,
                      "springs/compliant__C_v1v2_N_per_m": 0.75 * C})
        p = build_params(cfg, 1.0, 1.0)
        c = compare(p, StartOff(v_f=V_F, T=T_REF), t_tail=TAIL)
        rows.append({"sweep": "coupling stiffness", "value": C, "unit": "N/m",
                     "T_over_Tn": T_REF / coupling_period(p), **_metrics(c)})

    for al in (0.0, 0.05, 0.1, 0.5, 1.0, 3.0, 6.0, 10.0, 20.0):
        cfg = _cfg(**{"springs/compliant__alpha_v1_Ns_per_m_per_kg": al,
                      "springs/compliant__alpha_v2_Ns_per_m_per_kg": 0.5 * al})
        p = build_params(cfg, 1.0, 1.0)
        zeta = al * p.m_unit_kg[0] / (2 * np.sqrt(p.C_chv1_N_per_m * p.m_unit_kg[0]))
        c = compare(p, StartOff(v_f=V_F, T=T_REF), t_tail=TAIL)
        rows.append({"sweep": "specific damping", "value": al,
                     "unit": "N s m^-1 kg^-1", "zeta": zeta,
                     "T_over_Tn": T_REF / coupling_period(p), **_metrics(c)})

    for f in (0.0, 0.02, 0.05, 0.10, 0.20, 0.40):
        cfg = _cfg(**{"springs/compliant__f_ch": f, "springs/compliant__f_k": f})
        p = build_params(cfg, 1.0, 1.0)
        c = compare(p, StartOff(v_f=V_F, T=T_REF), t_tail=TAIL)
        rows.append({"sweep": "guide resistance", "value": f, "unit": "-",
                     "T_over_Tn": T_REF / coupling_period(p), **_metrics(c)})

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s7_parametric.csv", index=False)
    print(df.groupby("sweep")[["dF_peak_pct", "u1_peak_mm"]].agg(["min", "max"]).to_string())
    return df


def s8_asymmetric() -> pd.DataFrame:
    """Asymmetric loading; the source study contains a 100 % / 60 % example."""
    rows = []
    for lam1, lam2 in ((1.0, 1.0), (1.0, 0.6), (0.6, 1.0), (1.0, 0.0), (0.0, 1.0),
                       (0.5, 0.5)):
        p = build_params(load_all(), lam1, lam2)
        c = compare(p, StartOff(v_f=V_F, T=T_REF), t_tail=TAIL)
        rows.append({"front_fill": lam1, "rear_fill": lam2,
                     "m1_kg": p.m_unit_kg[0], "m2_kg": p.m_unit_kg[1],
                     "asymmetry": abs(p.m_unit_kg[0] - p.m_unit_kg[1])
                                  / sum(p.m_unit_kg), **_metrics(c)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "s8_asymmetric.csv", index=False)
    print(df.to_string(index=False))
    return df


def _metrics(c: dict) -> dict:
    return {"F_rigid_kN": c["rigid"]["F_peak"] / 1e3,
            "F_seg_kN": c["segmented"]["F_peak"] / 1e3,
            "dF_peak_pct": c["dF_peak_pct"],
            "dP_peak_pct": c["dP_peak_pct"],
            "dW_pct": c["dW_pct"],
            "u1_peak_mm": c["u1_peak_mm"], "u2_peak_mm": c["u2_peak_mm"]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("\n== S2 damping ablation =="); s2_damping_ablation()
    print("\n== S3 forward start-off =="); s3_forward()
    print("\n== S4 inverse demand =="); s4_inverse()
    print("\n== S6 duration sweep =="); s6_duration_sweep()
    print("\n== S7 parametric =="); s7_parametric()
    print("\n== S8 asymmetric =="); s8_asymmetric()
    print("\npaper studies complete ->", OUT)


if __name__ == "__main__":
    main()
