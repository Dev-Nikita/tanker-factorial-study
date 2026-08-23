"""Identify the four model parameters that the source manuscript does not report.

SRC-A reports masses, stiffnesses, damping and resistance coefficients, but NOT

  P-1  m_tractor  - tractor mass (only the tanker chassis mass 5900 kg is given);
  P-2  r_k        - dynamic radius of the driving wheels;
  P-3  I_k        - driveline inertia reduced to the driving-wheel axes;
  P-4  D_kch      - tangential (circumferential) damping of the driving wheels.

All four are required to integrate eq. (1)-(2).  Rather than inventing them they
are identified from quantities that SRC-A *does* report, in two decoupled stages:

  Stage A  (m_tractor, r_k)  <- stage-1 integral quantities over t = 0...600 s
           (engine work, mean power, mean wheel torque, distance, mean and
            terminal speed) at m_v = 5000 kg, f_d = 0.03.
  Stage B  (I_k, D_kch)      <- the four rigid-configuration peak wheel forces of
           SRC-A Table 1, which are governed by the start-off transient only.

The staging is legitimate because the stage-1 integrals are insensitive to the
transient parameters and vice versa; this is verified by a cross-check at the end.

Identified values and residual errors -> results/tables/parameter_identification.csv
"""
from __future__ import annotations

import os
import sys
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.config import load_all, build_params      # noqa: E402
from tanker_dynamics.solver import simulate, responses         # noqa: E402

TARGETS_STAGE1 = {
    "E_engine_J": 5.337e7,
    "P_mean_W": 8.895e4,
    "M_mean_Nm": 3.835e3,
    "distance_m": 7.582e3,
    "v_mean_m_s": 12.536,
    "v_end_m_s": 13.058,
}
TARGETS_TABLE1 = {"rigid": {1250: 3.946e4, 2500: 4.547e4, 3750: 5.027e4, 5000: 5.375e4},
                  "compliant": {1250: 3.853e4, 2500: 4.163e4, 3750: 4.622e4, 5000: 4.877e4}}
T_TRANSIENT = 6.0
FAST = dict(n_out=3001, rtol=1e-6, atol=1e-8)   # identification loop
FINE = dict(n_out=8001, rtol=1e-8, atol=1e-10)  # final verification


def _set(cfg, **kw):
    m = {"m_tractor_kg": ("tractor",), "r_k_m": ("wheels",),
         "I_k_kgm2": ("wheels",), "D_kch_Nms_per_rad": ("wheels",)}
    for k, v in kw.items():
        cfg["vehicle"][m[k][0]][k] = float(v)


def eval_stage1(cfg):
    p = build_params(cfg, 1.0, 1.0)
    r = responses(simulate(p, "rigid", (0.0, 600.0), n_out=6001))
    return {k: r[k] for k in TARGETS_STAGE1}


def eval_table1(cfg, config="rigid", opts=None):
    opts = opts or FAST
    full = cfg["tanks"]["full_scale"]["m_liquid_full_kg"]
    out = {}
    for m_v in TARGETS_TABLE1["rigid"]:
        lam = m_v / full
        p = build_params(cfg, lam, lam)
        out[m_v] = responses(simulate(p, config, (0.0, T_TRANSIENT), **opts))["F_kch_peak_N"]
    return out


def cost_A(z, cfg):
    _set(cfg, m_tractor_kg=np.exp(z[0]), r_k_m=np.exp(z[1]))
    try:
        ev = eval_stage1(cfg)
    except Exception:
        return 1e6
    return float(np.mean([((ev[k] - v) / v) ** 2 for k, v in TARGETS_STAGE1.items()]))


def cost_B(z, cfg):
    _set(cfg, I_k_kgm2=np.exp(z[0]), D_kch_Nms_per_rad=np.exp(z[1]))
    try:
        e = {c: eval_table1(cfg, c) for c in ("rigid", "compliant")}
    except Exception:
        return 1e6
    err = [((e[c][m] - v) / v) ** 2
           for c, tg in TARGETS_TABLE1.items() for m, v in tg.items()]
    return float(np.mean(err))


def main():
    cfg = load_all()

    if os.environ.get("SKIP_STAGE_A"):
        m_tr = float(cfg["vehicle"]["tractor"]["m_tractor_kg"])
        r_k = float(cfg["vehicle"]["wheels"]["r_k_m"])
        print(f"Stage A: reusing stored values m_tractor={m_tr}, r_k={r_k}", flush=True)
        _set(cfg, m_tractor_kg=m_tr, r_k_m=r_k)
    else:
        print("Stage A: identifying m_tractor, r_k ...", flush=True)
        rA = minimize(cost_A, np.log([6700.0, 0.60]), args=(cfg,), method="Nelder-Mead",
                      options={"xatol": 1e-4, "fatol": 1e-9, "maxiter": 120})
        m_tr, r_k = np.exp(rA.x)
        _set(cfg, m_tractor_kg=m_tr, r_k_m=r_k)
        print(f"  m_tractor = {m_tr:.1f} kg, r_k = {r_k:.4f} m, "
              f"RMS = {np.sqrt(rA.fun):.4%}", flush=True)

    print("Stage B: identifying I_k, D_kch by bounded coarse-to-fine grid ...", flush=True)
    lo = np.log([2500.0, 1.5e3])
    hi = np.log([8000.0, 1.5e4])
    best_z, best_f = None, np.inf
    for level in range(3):
        for a in np.linspace(lo[0], hi[0], 5):
            for b in np.linspace(lo[1], hi[1], 5):
                f = cost_B(np.array([a, b]), cfg)
                if f < best_f:
                    best_f, best_z = f, np.array([a, b])
        span = (hi - lo) / 4.0
        lo, hi = best_z - span, best_z + span
        print(f"  level {level + 1}: I_k = {np.exp(best_z[0]):.0f} kg m^2, "
              f"D_kch = {np.exp(best_z[1]):.4g}, RMS = {np.sqrt(best_f):.4%}", flush=True)
    rB = SimpleNamespace(x=best_z, fun=best_f)
    I_k, D_k = np.exp(rB.x)
    _set(cfg, I_k_kgm2=I_k, D_kch_Nms_per_rad=D_k)
    print(f"  I_k = {I_k:.1f} kg m^2, D_kch = {D_k:.4g} N m s/rad, "
          f"RMS rel. err = {np.sqrt(rB.fun):.4%}")

    print("Cross-check ...")
    ev1 = eval_stage1(cfg)
    evT = {c: eval_table1(cfg, c, FINE) for c in ("rigid", "compliant")}

    rows = [{"parameter": "m_tractor_kg", "identified": m_tr, "stage": "A"},
            {"parameter": "r_k_m", "identified": r_k, "stage": "A"},
            {"parameter": "I_k_kgm2", "identified": I_k, "stage": "B"},
            {"parameter": "D_kch_Nms_per_rad", "identified": D_k, "stage": "B"}]
    for k, v in TARGETS_STAGE1.items():
        rows.append({"parameter": k, "source_value": v, "reproduced": ev1[k],
                     "rel_error_pct": 100 * (ev1[k] - v) / v, "stage": "A-target"})
    for c, tg in TARGETS_TABLE1.items():
        for m, v in tg.items():
            rows.append({"parameter": f"F_kch_peak_N {c} (m_v={m} kg)", "source_value": v,
                         "reproduced": evT[c][m], "rel_error_pct": 100 * (evT[c][m] - v) / v,
                         "stage": "B-target"})
    df = pd.DataFrame(rows)
    out = ROOT / "results" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "parameter_identification.csv", index=False)
    print(df.to_string(index=False))

    print("\n--- paste into config/vehicle.yaml ---")
    print(f"  m_tractor_kg: {m_tr:.1f}")
    print(f"  r_k_m: {r_k:.4f}")
    print(f"  I_k_kgm2: {I_k:.1f}")
    print(f"  D_kch_Nms_per_rad: {D_k:.4g}")


if __name__ == "__main__":
    main()
