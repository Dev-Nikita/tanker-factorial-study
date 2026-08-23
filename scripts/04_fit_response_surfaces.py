"""Fit the second-order response surfaces and run the ANOVA for every response.

Primary model: 10 coefficients on 27 runs -> 17 residual degrees of freedom.
A saturated 27-term tensor polynomial is fitted for completeness only and is
labelled as non-inferential.

Outputs (results/tables/):
  rsm_coefficients_<scale>_<config>.csv
  rsm_anova_<scale>_<config>.csv
  rsm_summary_<scale>_<config>.csv
  rsm_saturated_note.txt
  effect_ranking_<scale>_<config>.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics import regression as rg                   # noqa: E402
from tanker_dynamics.config import load_all                    # noqa: E402

RESPONSES = [
    "P_peak_W", "P_mean_W", "E_engine_J", "M_peak_Nm", "M_mean_Nm",
    "F_kch_peak_N", "F_chv1_peak_N", "F_v1v2_peak_N",
    "a_peak_m_s2", "a_rms_m_s2", "j_peak_m_s3", "j_rms_m_s3",
    "u1_peak_m", "u2_peak_m", "t_manoeuvre_s", "distance_m",
]

TDIR = ROOT / "results" / "tables"


def analyse(design: pd.DataFrame, resp: pd.DataFrame, scale: str, config: str):
    df = design.merge(resp[resp.configuration == config], on="run_id", how="inner")
    coeffs, anovas, summaries, ranks = [], [], [], []
    for y in RESPONSES:
        if y not in df.columns or df[y].std(ddof=1) < 1e-14:
            summaries.append({"response": y, "note": "constant response - no model fitted"})
            continue
        res = rg.fit(df, y)
        coeffs.append(rg.coefficient_table(res, y))
        anovas.append(rg.anova_table(df, y))
        summaries.append(rg.summary_row(df, y))
        ct = rg.coefficient_table(res, y)
        ct = ct[ct.term != "const"].copy()
        ct["abs_std_coefficient"] = ct.std_coefficient.abs()
        ct = ct.sort_values("abs_std_coefficient", ascending=False)
        ct["rank"] = np.arange(1, len(ct) + 1)
        ranks.append(ct[["response", "term", "coefficient", "std_coefficient",
                         "p_value", "rank"]])

    tag = f"{scale}_{config}"
    if coeffs:
        pd.concat(coeffs).to_csv(TDIR / f"rsm_coefficients_{tag}.csv", index=False)
        pd.concat(anovas).to_csv(TDIR / f"rsm_anova_{tag}.csv", index=False)
        pd.concat(ranks).to_csv(TDIR / f"effect_ranking_{tag}.csv", index=False)
    pd.DataFrame(summaries).to_csv(TDIR / f"rsm_summary_{tag}.csv", index=False)
    return df, pd.DataFrame(summaries)


def saturated_check(df: pd.DataFrame, y: str) -> str:
    X = sm.add_constant(rg.saturated_matrix(df), has_constant="add")
    r = sm.OLS(df[y].to_numpy(float), X).fit()
    return (f"Saturated 27-term tensor polynomial for {y}: "
            f"R2 = {r.rsquared:.6f}, residual d.o.f. = {int(r.df_resid)}. "
            "SATURATED INTERPOLATING MODEL - NOT THE PRIMARY INFERENTIAL MODEL. "
            "With 27 coefficients on 27 design points the residual error cannot be "
            "estimated and no significance statement may be derived from it.")


def main():
    cfg = load_all()
    notes = []
    for scale in cfg["experiment"]["scales"]:
        design = pd.read_csv(TDIR / f"factorial_design_{scale}.csv")
        resp = pd.read_csv(TDIR / f"factorial_responses_{scale}.csv")
        resp = resp[resp.scenario == "start_off"]
        for config in ("rigid", "compliant"):
            df, summ = analyse(design, resp, scale, config)
            print(f"--- {scale} / {config} ---")
            print(summ.to_string(index=False))
            if scale == "full_scale" and config == "compliant":
                for y in ("P_peak_W", "F_kch_peak_N"):
                    notes.append(saturated_check(df, y))
    (TDIR / "rsm_saturated_note.txt").write_text("\n\n".join(notes), encoding="utf-8")
    print("\n".join(notes))


if __name__ == "__main__":
    main()
