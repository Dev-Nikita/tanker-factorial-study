"""Export manuscript-ready tables (CSV + Markdown + LaTeX) and the results report."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.config import load_all                     # noqa: E402

TDIR = ROOT / "results" / "tables"
MDIR = ROOT / "results" / "manuscript_tables"
SCALE = "full_scale"


def export(df: pd.DataFrame, name: str, caption: str, float_fmt="%.4g"):
    MDIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(MDIR / f"{name}.csv", index=False)
    (MDIR / f"{name}.md").write_text(f"**{caption}**\n\n" +
                                     df.to_markdown(index=False, floatfmt=".4g"),
                                     encoding="utf-8")
    (MDIR / f"{name}.tex").write_text(
        df.to_latex(index=False, float_format=float_fmt, caption=caption,
                    label=f"tab:{name}", longtable=False, escape=True),
        encoding="utf-8")
    print("table:", name)


def main():
    cfg = load_all()
    MDIR.mkdir(parents=True, exist_ok=True)

    # --- Table 1: physical parameters -------------------------------------- #
    v, s, t = cfg["vehicle"], cfg["springs"]["compliant"], cfg["tanks"]["full_scale"]
    rows = [
        ("m_ch", "tanker chassis mass", v["chassis"]["m_ch_kg"], "kg", "SRC-A"),
        ("m_tr", "tractor mass", v["tractor"]["m_tractor_kg"], "kg", "identified (P-1)"),
        ("m_k", "movable trolley mass (each)", t["m_trolley_kg"], "kg", "SRC-A"),
        ("m_v,max", "liquid cargo at 100 % fill (each tank)", t["m_liquid_full_kg"], "kg", "SRC-A"),
        ("r_k", "dynamic wheel radius", v["wheels"]["r_k_m"], "m", "identified (P-2)"),
        ("I_k", "driveline inertia reduced to wheel axes", v["wheels"]["I_k_kgm2"],
         "kg m^2", "identified (P-3)"),
        ("C_kch", "tangential wheel stiffness", v["wheels"]["C_kch_N_per_rad"],
         "N m/rad", "SRC-A"),
        ("D_kch", "tangential wheel damping", v["wheels"]["D_kch_Nms_per_rad"],
         "N m s/rad", "added term, identified (P-4)"),
        ("C_chv1", "chassis-unit 1 spring rate", s["C_chv1_N_per_m"], "N/m", "SRC-A"),
        ("C_v1v2", "unit 1-unit 2 spring rate", s["C_v1v2_N_per_m"], "N/m", "SRC-A"),
        ("alpha_v1", "specific damping, chassis-unit 1", s["alpha_v1_Ns_per_m_per_kg"],
         "N s m^-1 kg^-1", "SRC-A"),
        ("alpha_v2", "specific damping, unit 1-unit 2", s["alpha_v2_Ns_per_m_per_kg"],
         "N s m^-1 kg^-1", "SRC-A"),
        ("f_d", "road rolling-resistance coefficient", v["road"]["f_d"], "-", "SRC-A"),
        ("f_ch", "guide resistance, unit 1", s["f_ch"], "-", "SRC-A"),
        ("f_k", "guide resistance, unit 2", s["f_k"], "-", "SRC-A"),
        ("N_d", "wheel-referred plateau power",
         float(v["powertrain"]["M_peak_Nm"]) * float(v["powertrain"]["omega_peak_rad_s"]),
         "W", "SRC-A Fig. 4"),
    ]
    export(pd.DataFrame(rows, columns=["Symbol", "Quantity", "Value", "Unit", "Provenance"]),
           "table01_parameters", "Model parameters and their provenance.")

    # --- Table 2: factor definition ---------------------------------------- #
    f = cfg["experiment"]["factors"]
    rows = [{"Factor": k, "Physical quantity": f[k]["description"].strip().split("\n")[0],
             "Unit": f[k]["unit"],
             "Level -1": f[k]["levels"]["-1"], "Level 0": f[k]["levels"]["0"],
             "Level +1": f[k]["levels"]["+1"]} for k in ("x1", "x2", "x3")]
    export(pd.DataFrame(rows), "table02_factors",
           "Three-factor, three-level design definition.")

    # --- Table 3: the 27 design points ------------------------------------- #
    d = pd.read_csv(TDIR / f"factorial_design_{SCALE}.csv")
    export(d[["run_id", "randomized_run_id", "x1_coded", "x2_coded", "x3_coded",
              "front_fill_pct", "rear_fill_pct", "speed_factor_mps",
              "front_unit_mass_kg", "rear_unit_mass_kg", "total_movable_mass_kg",
              "load_asymmetry"]],
           "table03_design", "Full-factorial design matrix (3^3 = 27 runs).")

    # --- Table 4: reproduction --------------------------------------------- #
    export(pd.read_csv(TDIR / "original_reproduction.csv"),
           "table04_reproduction", "Reproduction of the source-model results.")

    # --- Table 5: identified parameters ------------------------------------ #
    export(pd.read_csv(TDIR / "parameter_identification.csv"),
           "table05_identification",
           "Parameters absent from the source and their identification residuals.")

    # --- Table 6: raw responses -------------------------------------------- #
    pr = pd.read_csv(TDIR / f"factorial_paired_{SCALE}.csv")
    keep = ["run_id", "x1_coded", "x2_coded", "x3_coded"]
    for c in ("F_kch_peak_N", "P_mean_W", "E_engine_J", "a_peak_m_s2",
              "j_peak_m_s3", "u1_peak_m"):
        keep += [f"{c}__rigid", f"{c}__compliant", f"{c}__reduction_pct"]
    export(pr[keep], "table06_responses",
           "Computational responses for all 27 design points, both configurations.")

    # --- Table 7/8: coefficients and ANOVA --------------------------------- #
    co = pd.read_csv(TDIR / f"rsm_coefficients_{SCALE}_compliant.csv")
    export(co[co.response == "F_kch_peak_N"].drop(columns=["response"]),
           "table07_coefficients_force",
           "Second-order response-surface coefficients for the peak wheel force "
           "(compliant configuration).")
    an = pd.read_csv(TDIR / f"rsm_anova_{SCALE}_compliant.csv")
    export(an[an.response.isin(["F_kch_peak_N", "P_mean_W", "a_peak_m_s2",
                                "j_peak_m_s3", "u1_peak_m"])],
           "table08_anova", "ANOVA of the second-order response surfaces.")

    # --- Table 9: model quality across responses --------------------------- #
    su = pd.read_csv(TDIR / f"rsm_summary_{SCALE}_compliant.csv")
    export(su, "table09_model_quality",
           "Response-surface quality and curvature assessment for every response.")

    # --- Table 10: effect ranking ------------------------------------------ #
    er = pd.read_csv(TDIR / f"effect_ranking_{SCALE}_compliant.csv")
    export(er[er["rank"] <= 4], "table10_effect_ranking",
           "Four strongest standardised effects per response.")

    # --- Table 11: extreme configurations ---------------------------------- #
    best = pr.loc[pr.F_kch_peak_N__compliant.idxmin()]
    worst = pr.loc[pr.F_kch_peak_N__compliant.idxmax()]
    bred = pr.loc[pr.F_kch_peak_N__reduction_pct.idxmax()]
    wred = pr.loc[pr.F_kch_peak_N__reduction_pct.idxmin()]
    ext = pd.DataFrame([
        {"case": "lowest peak load", **{k: best[k] for k in
         ("run_id", "front_fill_pct", "rear_fill_pct", "speed_factor_mps",
          "F_kch_peak_N__rigid", "F_kch_peak_N__compliant", "F_kch_peak_N__reduction_pct")}},
        {"case": "highest peak load", **{k: worst[k] for k in
         ("run_id", "front_fill_pct", "rear_fill_pct", "speed_factor_mps",
          "F_kch_peak_N__rigid", "F_kch_peak_N__compliant", "F_kch_peak_N__reduction_pct")}},
        {"case": "largest reduction", **{k: bred[k] for k in
         ("run_id", "front_fill_pct", "rear_fill_pct", "speed_factor_mps",
          "F_kch_peak_N__rigid", "F_kch_peak_N__compliant", "F_kch_peak_N__reduction_pct")}},
        {"case": "smallest reduction", **{k: wred[k] for k in
         ("run_id", "front_fill_pct", "rear_fill_pct", "speed_factor_mps",
          "F_kch_peak_N__rigid", "F_kch_peak_N__compliant", "F_kch_peak_N__reduction_pct")}},
    ])
    export(ext, "table11_extremes", "Extreme configurations of the design space.")

    # --- Table 12: uncertainty --------------------------------------------- #
    if (TDIR / "uncertainty_summary.csv").exists():
        export(pd.read_csv(TDIR / "uncertainty_summary.csv"), "table12_uncertainty",
               "Monte-Carlo uncertainty propagation at the centre point (n = 250).")

    # --- Table 13: braking -------------------------------------------------- #
    rs = pd.read_csv(TDIR / f"factorial_responses_{SCALE}.csv")
    br = rs[rs.scenario == "braking"]
    if len(br):
        g = br.groupby("configuration")[["braking_time_s", "braking_distance_m",
                                         "mean_decel_m_s2", "F_kch_peak_N",
                                         "a_peak_m_s2", "j_peak_m_s3"]].agg(["mean", "std"])
        g.columns = ["_".join(c) for c in g.columns]
        export(g.reset_index(), "table13_braking",
               "Braking stage: rigid versus compliant mounting over the 27 design points.")

    # --- summary numbers for the manuscript -------------------------------- #
    red = pr.F_kch_peak_N__reduction_pct
    facts = {
        "n_runs": len(pr),
        "reduction_min_pct": red.min(), "reduction_max_pct": red.max(),
        "reduction_mean_pct": red.mean(),
        "F_rigid_min_kN": pr.F_kch_peak_N__rigid.min() / 1e3,
        "F_rigid_max_kN": pr.F_kch_peak_N__rigid.max() / 1e3,
        "F_comp_min_kN": pr.F_kch_peak_N__compliant.min() / 1e3,
        "F_comp_max_kN": pr.F_kch_peak_N__compliant.max() / 1e3,
        "a_peak_rigid_mean": pr.a_peak_m_s2__rigid.mean(),
        "a_peak_comp_mean": pr.a_peak_m_s2__compliant.mean(),
        "j_peak_rigid_mean": pr.j_peak_m_s3__rigid.mean(),
        "j_peak_comp_mean": pr.j_peak_m_s3__compliant.mean(),
        "u1_peak_max_mm": 1e3 * pr.u1_peak_m__compliant.max(),
        "u2_peak_max_mm": 1e3 * pr.u2_peak_m__compliant.max(),
    }
    pd.Series(facts).to_csv(MDIR / "key_numbers.csv", header=["value"])
    print(pd.Series(facts).to_string())


if __name__ == "__main__":
    main()
