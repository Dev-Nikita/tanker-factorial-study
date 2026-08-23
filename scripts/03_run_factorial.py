"""Run the 3^3 full-factorial computational experiment.

For every design point and every configuration (rigid / compliant) two manoeuvres
are integrated:

  * start-off: from rest until the chassis first reaches the speed factor x3
    (SRC-B "швидкість рушання"; see docs/MISSING_INPUTS.md M-01);
  * braking:   from the terminal state of the SRC-A stage-1 run, with the brake
    torque law of SRC-A Fig. 5.

Outputs
  results/tables/factorial_design.csv
  results/tables/factorial_responses_<scale>.csv
  results/tables/factorial_paired_<scale>.csv
  data/generated/timeseries_<scale>.npz
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.config import load_all, build_params      # noqa: E402
from tanker_dynamics.doe import factorial_design, validate     # noqa: E402
from tanker_dynamics.solver import simulate, responses         # noqa: E402

TS_RUNS = {1, 5, 14, 23, 27}   # design points whose time histories are archived


def run_scale(cfg: dict, scale: str) -> tuple:
    design = factorial_design(cfg, scale)
    validate(design)
    recs, ts = [], {}
    for _, row in design.iterrows():
        p = build_params(cfg, row.lambda_1, row.lambda_2, scale=scale)
        for config in ("rigid", "compliant"):
            t0 = time.time()
            sol = simulate(p, config, (0.0, 120.0),
                           target_speed=float(row.speed_factor_mps), n_out=4001)
            if not sol.ok:
                raise RuntimeError(f"solver failure run {row.run_id} {config}: {sol.message}")
            r = responses(sol)
            r.update({"run_id": int(row.run_id), "configuration": config, "scenario": "start_off",
                      "scale": scale, "wall_time_s": time.time() - t0})
            recs.append(r)
            if int(row.run_id) in TS_RUNS:
                u1, u2, _, _ = sol.relative()
                ts[f"{config}_{int(row.run_id)}"] = np.vstack(
                    [sol.t, sol.v_ch, sol.wheel_force(), sol.accel(), sol.jerk(), u1, u2])

        # braking stage (SRC-A scenario 2), full-scale only - rig has no 600 s run
        if scale == "full_scale":
            for config in ("rigid", "compliant"):
                acc = simulate(p, config, (0.0, 600.0), n_out=2001)
                y_end = acc.y[:, -1].copy()
                br = simulate(p, config, (600.0, 720.0), braking=True, y0=y_end,
                              stop_at_rest=True, n_out=4001)
                r = responses(br)
                v = br.v_ch
                stopped = v[-1] <= 5.0e-3
                r["braking_time_s"] = float(br.t[-1] - 600.0) if stopped else np.nan
                r["braking_distance_m"] = float(br.X_ch[-1] - br.X_ch[0]) if stopped else np.nan
                r["mean_decel_m_s2"] = (-v[0] / r["braking_time_s"]) if stopped else np.nan
                r.update({"run_id": int(row.run_id), "configuration": config,
                          "scenario": "braking", "scale": scale})
                recs.append(r)

    res = pd.DataFrame(recs)
    return design, res, ts


def paired(design: pd.DataFrame, res: pd.DataFrame, scenario="start_off") -> pd.DataFrame:
    sub = res[res.scenario == scenario]
    rg = sub[sub.configuration == "rigid"].set_index("run_id")
    cp = sub[sub.configuration == "compliant"].set_index("run_id")
    cols = [c for c in rg.columns if rg[c].dtype.kind == "f"]
    out = design.set_index("run_id").copy()
    for c in cols:
        out[f"{c}__rigid"] = rg[c]
        out[f"{c}__compliant"] = cp[c]
        with np.errstate(divide="ignore", invalid="ignore"):
            out[f"{c}__reduction_pct"] = 100 * (rg[c] - cp[c]) / rg[c].replace(0, np.nan)
    return out.reset_index()


def main():
    cfg = load_all()
    tdir = ROOT / "results" / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    (ROOT / "data" / "generated").mkdir(parents=True, exist_ok=True)

    for scale in cfg["experiment"]["scales"]:
        print(f"=== {scale} ===")
        design, res, ts = run_scale(cfg, scale)
        design.to_csv(tdir / f"factorial_design_{scale}.csv", index=False)
        res.to_csv(tdir / f"factorial_responses_{scale}.csv", index=False)
        paired(design, res).to_csv(tdir / f"factorial_paired_{scale}.csv", index=False)
        if ts:
            np.savez_compressed(ROOT / "data" / "generated" / f"timeseries_{scale}.npz", **ts)
        print(res.groupby(["scenario", "configuration"])[["F_kch_peak_N", "P_peak_W"]].mean())

    # canonical design table (identical coding for both scales)
    factorial_design(cfg, "full_scale").to_csv(tdir / "factorial_design.csv", index=False)
    print("done")


if __name__ == "__main__":
    main()
