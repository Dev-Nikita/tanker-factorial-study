"""Stage 1 - reproduce the baseline results of the source manuscript (SRC-A).

Two blocks:
  A. Stage-1 ("acceleration-motion", t = 0...600 s) integral quantities.
  B. SRC-A Table 1 - peak longitudinal forces for rigid (1) and compliant (2)
     mounting at cargo masses 1250 / 2500 / 3750 / 5000 kg.

Outputs: results/tables/original_reproduction.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.config import load_all, build_params      # noqa: E402
from tanker_dynamics.solver import simulate, responses         # noqa: E402

STAGE1 = {
    "Engine work over stage 1, J": 5.337e7,
    "Mean engine power, W": 8.895e4,
    "Mean wheel torque, N m": 3.835e3,
    "Distance travelled, m": 7.582e3,
    "Mean speed, m/s": 12.536,
    "Terminal speed, m/s": 13.058,
}
KEYS = ["E_engine_J", "P_mean_W", "M_mean_Nm", "distance_m", "v_mean_m_s", "v_end_m_s"]

TABLE1 = {  # m_v -> (rigid, compliant), N
    1250: (3.946e4, 3.853e4),
    2500: (4.547e4, 4.163e4),
    3750: (5.027e4, 4.622e4),
    5000: (5.375e4, 4.877e4),
}


def main():
    cfg = load_all()
    rows = []

    # ---- block A ---------------------------------------------------------- #
    p = build_params(cfg, 1.0, 1.0)
    r = responses(simulate(p, "rigid", (0.0, 600.0), n_out=6001))
    for (label, src), key in zip(STAGE1.items(), KEYS):
        rep = r[key]
        rows.append({"block": "Stage 1 integral quantities", "quantity": label,
                     "source_value": src, "reproduced": rep,
                     "abs_error": rep - src, "rel_error_pct": 100 * (rep - src) / src})

    # ---- block B ---------------------------------------------------------- #
    full = cfg["tanks"]["full_scale"]["m_liquid_full_kg"]
    for m_v, (src_rigid, src_comp) in TABLE1.items():
        lam = m_v / full
        pr = build_params(cfg, lam, lam)
        f_rigid = responses(simulate(pr, "rigid", (0.0, 20.0), n_out=8001))["F_kch_peak_N"]
        f_comp = responses(simulate(pr, "compliant", (0.0, 20.0), n_out=8001))["F_kch_peak_N"]
        for cfg_name, src, rep in (("rigid", src_rigid, f_rigid),
                                   ("compliant", src_comp, f_comp)):
            rows.append({"block": "Table 1 - peak wheel force F_kch, N",
                         "quantity": f"m_v = {m_v} kg, {cfg_name}",
                         "source_value": src, "reproduced": rep,
                         "abs_error": rep - src, "rel_error_pct": 100 * (rep - src) / src})
        rows.append({"block": "Table 1 - load reduction, %",
                     "quantity": f"m_v = {m_v} kg",
                     "source_value": 100 * (src_rigid - src_comp) / src_rigid,
                     "reproduced": 100 * (f_rigid - f_comp) / f_rigid,
                     "abs_error": np.nan,
                     "rel_error_pct": np.nan,
                     "note": "percentage-point comparison, see manuscript Table 4"})

    df = pd.DataFrame(rows)
    df["status"] = np.where(df.rel_error_pct.isna(), "-",
                            np.where(df.rel_error_pct.abs() <= 5, "OK (<=5 %)",
                                     np.where(df.rel_error_pct.abs() <= 15,
                                              "ACCEPTABLE (<=15 %)", "DEVIATION")))
    out = ROOT / "results" / "tables" / "original_reproduction.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
