"""Monte-Carlo uncertainty propagation at the centre point of the design.

Deterministic replicates of an ODE are meaningless, so the computational
counterpart of experimental replication is uncertainty propagation through the
physically uncertain parameters.  Uncertainty magnitudes are read from
config/uncertainty.yaml; if that file is absent the stage is skipped and the
required inputs are reported instead of being invented.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.config import load_all, build_params      # noqa: E402
from tanker_dynamics.solver import simulate, responses         # noqa: E402

UFILE = ROOT / "config" / "uncertainty.yaml"
N_MC = 250
KEYS = ["F_kch_peak_N", "F_chv1_peak_N", "a_peak_m_s2", "j_peak_m_s3",
        "P_peak_W", "P_mean_W", "u1_peak_m"]


def main():
    if not UFILE.exists():
        print("config/uncertainty.yaml not found - uncertainty stage SKIPPED.\n"
              "Required inputs: tolerance bands for m_tank_dry, liquid density, "
              "C_chv1, C_v1v2, alpha_v1, alpha_v2, f_ch, f_k, f_d, drivetrain efficiency.")
        return
    u = yaml.safe_load(UFILE.read_text(encoding="utf-8"))
    cfg = load_all()
    rng = np.random.default_rng(int(u.get("seed", 42)))
    v_t = float(cfg["experiment"]["factors"]["x3"]["levels"]["0"])

    rows = []
    for i in range(N_MC):
        c = load_all()
        for path, spec in u["parameters"].items():
            sect, key = path.split(".")
            base = float(c[sect.split("/")[0]][sect.split("/")[1]][key]) \
                if "/" in sect else float(c[sect][key])
            if spec["dist"] == "normal":
                val = base * (1.0 + rng.normal(0.0, spec["rel_sd"]))
            elif spec["dist"] == "uniform":
                val = base * (1.0 + rng.uniform(-spec["rel_halfwidth"], spec["rel_halfwidth"]))
            else:
                raise ValueError(spec["dist"])
            if "/" in sect:
                c[sect.split("/")[0]][sect.split("/")[1]][key] = val
            else:
                c[sect][key] = val
        for config in ("rigid", "compliant"):
            p = build_params(c, 0.5, 0.5)
            r = responses(simulate(p, config, (0.0, 120.0), target_speed=v_t,
                                   n_out=3001, rtol=1e-7, atol=1e-9))
            rows.append({"sample": i, "configuration": config,
                         **{k: r[k] for k in KEYS}})

    df = pd.DataFrame(rows)
    tdir = ROOT / "results" / "tables"
    df.to_csv(tdir / "uncertainty_samples.csv", index=False)

    summ = []
    for config, g in df.groupby("configuration"):
        for k in KEYS:
            s = g[k]
            summ.append({"configuration": config, "response": k,
                         "mean": s.mean(), "sd": s.std(ddof=1),
                         "cv_pct": 100 * s.std(ddof=1) / s.mean(),
                         "p2_5": s.quantile(0.025), "p97_5": s.quantile(0.975)})
    piv = df.pivot_table(index="sample", columns="configuration", values="F_kch_peak_N")
    red = 100 * (piv["rigid"] - piv["compliant"]) / piv["rigid"]
    summ.append({"configuration": "reduction", "response": "F_kch_peak_N_reduction_pct",
                 "mean": red.mean(), "sd": red.std(ddof=1),
                 "cv_pct": np.nan, "p2_5": red.quantile(0.025), "p97_5": red.quantile(0.975)})
    pd.DataFrame(summ).to_csv(tdir / "uncertainty_summary.csv", index=False)
    print(pd.DataFrame(summ).to_string(index=False))


if __name__ == "__main__":
    main()
