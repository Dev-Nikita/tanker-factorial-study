"""Limiting-case and conservation tests for the dynamic model."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.config import load_all, build_params, unit_masses   # noqa: E402
from tanker_dynamics.model import TankerParams, Powertrain, rhs_compliant  # noqa: E402
from tanker_dynamics.solver import simulate, responses                   # noqa: E402


@pytest.fixture(scope="module")
def cfg():
    return load_all()


# --- Test A: mass model ---------------------------------------------------- #
def test_mass_model(cfg):
    t = cfg["tanks"]["full_scale"]
    base = t["m_trolley_kg"] + t["m_tank_dry_kg"]
    assert unit_masses(cfg, 0.0, 0.0) == (base, base)
    m1, m2 = unit_masses(cfg, 1.0, 0.5)
    assert m1 == pytest.approx(base + t["m_liquid_full_kg"])
    assert m2 == pytest.approx(base + 0.5 * t["m_liquid_full_kg"])
    r = cfg["tanks"]["rig_scale"]
    assert r["m_liquid_full_kg"] == pytest.approx(
        cfg["tanks"]["liquid"]["density_kg_m3"] * r["tank_volume_m3"])


# --- Test B: symmetry ------------------------------------------------------ #
def test_symmetric_filling_gives_equal_unit_masses(cfg):
    p = build_params(cfg, 0.5, 0.5)
    assert p.m_unit_kg[0] == pytest.approx(p.m_unit_kg[1])


# --- Test C: stiff limit approaches the rigid baseline --------------------- #
def test_stiff_limit_approaches_rigid(cfg):
    import copy
    c = copy.deepcopy(cfg)
    c["springs"]["compliant"]["C_chv1_N_per_m"] = 5.0e9
    c["springs"]["compliant"]["C_v1v2_N_per_m"] = 5.0e9
    c["springs"]["compliant"]["f_ch"] = 0.0
    c["springs"]["compliant"]["f_k"] = 0.0
    p_stiff = build_params(c, 1.0, 1.0)
    p_rigid = build_params(cfg, 1.0, 1.0)
    fs = responses(simulate(p_stiff, "compliant", (0.0, 5.0), n_out=2001,
                            rtol=1e-8, atol=1e-10))["F_kch_peak_N"]
    fr = responses(simulate(p_rigid, "rigid", (0.0, 5.0), n_out=2001,
                            rtol=1e-8, atol=1e-10))["F_kch_peak_N"]
    assert abs(fs - fr) / fr < 0.02, f"stiff limit {fs:.1f} vs rigid {fr:.1f}"


# --- Test D: zero damping keeps the spring-only limit finite --------------- #
def test_zero_damping_runs(cfg):
    import copy
    c = copy.deepcopy(cfg)
    c["springs"]["compliant"]["alpha_v1_Ns_per_m_per_kg"] = 0.0
    c["springs"]["compliant"]["alpha_v2_Ns_per_m_per_kg"] = 0.0
    sol = simulate(build_params(c, 1.0, 1.0), "compliant", (0.0, 5.0), n_out=1001)
    assert sol.ok and np.all(np.isfinite(sol.y))


# --- Test E/F: energy behaviour of the unforced system --------------------- #
def _mech_energy(p: TankerParams, y: np.ndarray) -> np.ndarray:
    _, _, u1, u2, dphi, dX, du1, du2 = y
    m1, m2 = p.m_unit_kg
    m0 = p.m_tractor_kg + p.m_ch_kg
    v1, v2 = dX + du1, dX + du1 + du2
    T = 0.5 * (m0 * dX**2 + m1 * v1**2 + m2 * v2**2 + p.I_k_kgm2 * dphi**2)
    V = 0.5 * (p.C_chv1_N_per_m * u1**2 + p.C_v1v2_N_per_m * u2**2)
    return T + V


def test_energy_conserved_without_dissipation(cfg):
    p = build_params(cfg, 1.0, 1.0)
    p.alpha_v1 = p.alpha_v2 = 0.0
    p.f_d = p.f_ch = p.f_k = 0.0
    p.D_kch_Nms_per_rad = 0.0
    p.C_kch_N_per_rad = 0.0            # decouple the wheel to remove the traction path
    p.powertrain = Powertrain(M0_Nm=0.0, M_peak_Nm=0.0)
    y0 = np.array([0, 0, 0.05, -0.03, 0, 1.0, 0.0, 0.0])
    sol = simulate(p, "compliant", (0.0, 20.0), y0=y0, n_out=4001, rtol=1e-10, atol=1e-12)
    E = _mech_energy(p, sol.y)
    assert (E.max() - E.min()) / E.mean() < 1e-5


def test_energy_non_increasing_with_damping(cfg):
    p = build_params(cfg, 1.0, 1.0)
    p.f_d = p.f_ch = p.f_k = 0.0
    p.D_kch_Nms_per_rad = 0.0
    p.C_kch_N_per_rad = 0.0
    p.powertrain = Powertrain(M0_Nm=0.0, M_peak_Nm=0.0)
    y0 = np.array([0, 0, 0.05, -0.03, 0, 1.0, 0.0, 0.0])
    sol = simulate(p, "compliant", (0.0, 20.0), y0=y0, n_out=4001, rtol=1e-10, atol=1e-12)
    E = _mech_energy(p, sol.y)
    assert E[-1] <= E[0] * (1 + 1e-8)
    assert np.max(np.diff(E)) < 1e-6 * E[0]


# --- Test G: dimensional / sanity checks ----------------------------------- #
def test_static_equilibrium_is_a_fixed_point(cfg):
    p = build_params(cfg, 1.0, 1.0)
    p.powertrain = Powertrain(M0_Nm=0.0, M_peak_Nm=0.0)
    dy = rhs_compliant(0.0, np.zeros(8), p, braking=False)
    assert np.allclose(dy, 0.0, atol=1e-9)


def test_powertrain_characteristic_matches_source(cfg):
    pw = build_params(cfg, 1.0, 1.0).powertrain
    assert pw.torque(0.0) == pytest.approx(pw.M0_Nm)
    assert pw.torque(pw.omega_peak_rad_s) == pytest.approx(pw.M_peak_Nm)
    # constant-power ("ideal") region, SRC-A eq. (3)
    for w in (5.0, 12.0, 20.0):
        assert pw.power(w) == pytest.approx(pw.N_d_W, rel=1e-9)
    assert pw.torque(pw.omega_max_rad_s + 1.0) == 0.0
    assert pw.torque(0.0) > 0 and pw.N_d_W > 0


# --- Test H: internal forces cancel (Newton's third law) ------------------- #
def test_internal_forces_cancel(cfg):
    """Sum of the three body equations must leave only F_kch - R_road.

    Guide friction and the coupling forces are internal to the vehicle; if any of
    them is applied to one body without its reaction on the other, the model gains
    or loses momentum. This test previously failed: the reaction of the guide
    friction on the frame was missing.
    """
    import copy
    c = copy.deepcopy(cfg)
    c["springs"]["compliant"]["f_ch"] = 0.25      # large, so the error would be visible
    c["springs"]["compliant"]["f_k"] = 0.25
    p = build_params(c, 1.0, 0.4)
    y = np.array([0.3, 1.2, 0.02, -0.015, 4.0, 2.0, 0.30, -0.20])
    dy = rhs_compliant(0.0, y, p)
    a_ch, ddu1, ddu2 = dy[5], dy[6], dy[7]
    m1, m2 = p.m_unit_kg
    m0 = p.m_tractor_kg + p.m_ch_kg
    total = m0 * a_ch + m1 * (a_ch + ddu1) + m2 * (a_ch + ddu1 + ddu2)

    twist = y[0] - y[1] / p.r_k_m
    dtwist = y[4] - y[5] / p.r_k_m
    F_kch = (p.C_kch_N_per_rad * twist + p.D_kch_Nms_per_rad * dtwist) / p.r_k_m
    R_road = p.f_d * (m0 + m1 + m2) * p.g * np.tanh(y[5] / p.eps_v)
    assert abs(total - (F_kch - R_road)) < 1e-6 * abs(F_kch)


# --- Test I: matched-manoeuvre energy cannot decrease (passive system) ------ #
def test_matched_energy_never_decreases():
    """A passive coupling can only add dissipation, so it can never lower the
    traction work required for an identical manoeuvre. This test caught a sign
    error in the inverse-dynamics routine that reported a 26 % energy saving."""
    from tanker_dynamics.archive_model import ArchiveParams
    from tanker_dynamics.archive_studies import Manoeuvre, compare_matched
    p = ArchiveParams()
    for T in (3.0, 8.0, 20.0):
        c = compare_matched(p, Manoeuvre(v_f=5.0, T=T))
        assert c["dW_pct"] <= 1e-6, (
            f"T={T}: segmented demands {c['dW_pct']:.2f} % less work - impossible")


# --- Test J: archive scalars match the extracted worksheet values ----------- #
def test_archive_scalars():
    from tanker_dynamics.archive_model import ArchiveParams, V1, PT1
    p = ArchiveParams()
    assert p.Rk == 0.53 and p.Mch == 5900.0 and p.Mkuz == 700.0
    assert p.Mv == 5000.0 and p.Itr == 2025.0
    assert p.Alf_d == 0.06 and p.Alf_ch == 0.02 and p.Alf_k == 0.02
    assert p.CKCH == 1.2e6 and p.CCHK == 2.0e5 and p.CKV == 2.0e5
    assert p.Mgal_max == 1600.0
    assert p.MSchkv == 11600.0 and p.MSkv == 5700.0
    assert abs(p.Mchk - (1 / p.Mch + 1 / p.Mkuz)) < 1e-15
    assert abs(p.Mkv - (1 / p.Mkuz + 1 / p.Mv)) < 1e-15
    assert len(V1) == len(PT1) == 19 and PT1.max() == 65472.4
