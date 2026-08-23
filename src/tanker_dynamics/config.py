"""Configuration loading and parameter assembly."""

from __future__ import annotations

from pathlib import Path
import yaml

from .model import TankerParams, Powertrain, Brakes

ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "config"
RESULTS = ROOT / "results"
DATA = ROOT / "data"


def load_all() -> dict:
    cfg = {}
    for name in ("vehicle", "tanks", "springs", "experiment"):
        with open(CONFIG_DIR / f"{name}.yaml", encoding="utf-8") as fh:
            cfg[name] = yaml.safe_load(fh)
    return cfg


def unit_masses(cfg: dict, lam1: float, lam2: float, scale: str = "full_scale") -> tuple:
    """m_i = m_trolley + m_tank_dry + rho * V * lambda_i  (SRC-B mass model)."""
    t = cfg["tanks"][scale]
    base = t["m_trolley_kg"] + t["m_tank_dry_kg"]
    full = t["m_liquid_full_kg"]
    return (base + full * lam1, base + full * lam2)


def build_params(cfg: dict, lam1: float, lam2: float,
                 scale: str = "full_scale", f_d: float | None = None,
                 config: str = "compliant") -> TankerParams:
    v, s = cfg["vehicle"], cfg["springs"]["compliant"]
    pw = v["powertrain"]
    return TankerParams(
        m_ch_kg=float(v["chassis"]["m_ch_kg"]),
        m_tractor_kg=float(v["tractor"]["m_tractor_kg"]),
        m_unit_kg=unit_masses(cfg, lam1, lam2, scale),
        I_k_kgm2=float(v["wheels"]["I_k_kgm2"]),
        r_k_m=float(v["wheels"]["r_k_m"]),
        C_kch_N_per_rad=float(v["wheels"]["C_kch_N_per_rad"]),
        D_kch_Nms_per_rad=float(v["wheels"].get("D_kch_Nms_per_rad", 0.0)),
        C_chv1_N_per_m=float(s["C_chv1_N_per_m"]),
        C_v1v2_N_per_m=float(s["C_v1v2_N_per_m"]),
        alpha_v1=float(s["alpha_v1_Ns_per_m_per_kg"]),
        alpha_v2=float(s["alpha_v2_Ns_per_m_per_kg"]),
        f_d=float(v["road"]["f_d"] if f_d is None else f_d),
        f_ch=float(s["f_ch"]),
        f_k=float(s["f_k"]),
        g=float(v["road"]["g"]),
        eps_v=float(v["numerics"]["sign_smoothing_eps_m_s"]),
        powertrain=Powertrain(
            M0_Nm=float(pw["M0_Nm"]),
            omega_peak_rad_s=float(pw["omega_peak_rad_s"]),
            M_peak_Nm=float(pw["M_peak_Nm"]),
            omega_ideal_end_rad_s=float(pw["omega_ideal_end_rad_s"]),
            omega_max_rad_s=float(pw["omega_max_rad_s"]),
        ),
        brakes=Brakes(
            M_inf_Nm=float(v["brakes"]["M_inf_Nm"]),
            M_drop_Nm=float(v["brakes"]["M_drop_Nm"]),
            tau_s=float(v["brakes"]["tau_s"]),
            scale=float(v["brakes"].get("scale", 1.0)),
            t_start_s=float(cfg["experiment"]["scenarios"]["braking"]["t_brake_start_s"]),
        ),
    )
