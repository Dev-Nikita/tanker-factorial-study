"""Inverse dynamics: traction demand for a prescribed start-off manoeuvre.

The forward model answers "given this engine, how does the vehicle move?".  Because
the source powertrain characteristic holds constant power over the working range, that
question cannot reveal a change in *power demand* - the engine simply sits on its
plateau in every configuration.

The engineering claim under test is the opposite one: for the **same** start-off
manoeuvre, does segmenting the cargo inertia reduce the traction force, wheel torque
and power the tractor must supply?  That is an inverse-dynamics question and is what
this module answers.

Prescribed motion (C2, zero jerk at both ends):

    v_ref(t) = v_f * s^2 * (3 - 2 s),        s = t / T
    a_ref(t) = 6 v_f s (1 - s) / T

Rigid configuration - the whole vehicle is one body:

    F_req = m_tot a_ref + f_d m_tot g

Segmented configuration - the chassis follows the prescribed motion and the two cargo
units respond through the coupling elements; the traction force is whatever the frame
needs once the (lagging) coupling reaction is accounted for:

    m_0 a_ref = F_req - R_road + F_1     ->     F_req = m_0 a_ref + R_road - F_1

with the relative motion obtained by integrating

    m_1 (a_ref + u1'') = -F_1 - R_1 + F_2
    m_2 (a_ref + u1'' + u2'') = -F_2 - R_2
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .model import TankerParams, _sgn


# --------------------------------------------------------------------------- #
#  Prescribed manoeuvre
# --------------------------------------------------------------------------- #
@dataclass
class StartOff:
    """Smoothstep start-off from rest to ``v_f`` in ``T`` seconds."""

    v_f: float = 5.0
    T: float = 4.0

    def v(self, t):
        s = np.clip(np.asarray(t, float) / self.T, 0.0, 1.0)
        return self.v_f * s ** 2 * (3.0 - 2.0 * s)

    def a(self, t):
        s = np.clip(np.asarray(t, float) / self.T, 0.0, 1.0)
        return 6.0 * self.v_f * s * (1.0 - s) / self.T

    @property
    def a_max(self) -> float:
        return 1.5 * self.v_f / self.T


# --------------------------------------------------------------------------- #
#  Traction demand
# --------------------------------------------------------------------------- #
def demand_rigid(p: TankerParams, man: StartOff, t: np.ndarray) -> dict:
    a = man.a(t)
    v = man.v(t)
    m = p.m_total_kg
    F = m * a + p.f_d * m * p.g * np.tanh(v / p.eps_v)
    return _pack(t, v, a, F, p, m_eff=np.full_like(t, m))


def demand_segmented(p: TankerParams, man: StartOff, t: np.ndarray) -> dict:
    """Chassis motion prescribed; the two cargo units respond through the couplings."""
    m1, m2 = p.m_unit_kg
    m0 = p.m_tractor_kg + p.m_ch_kg

    def rhs(tt, y):
        u1, u2, du1, du2 = y
        a_b = float(man.a(tt))
        F1 = p.C_chv1_N_per_m * u1 + p.alpha_v1 * m1 * du1
        F2 = p.C_v1v2_N_per_m * u2 + p.alpha_v2 * m2 * du2
        R1 = p.f_ch * m1 * p.g * _sgn(du1, p.eps_v)
        R2 = p.f_k * m2 * p.g * _sgn(du2, p.eps_v)
        # guide friction is an internal pair (see model.rhs_compliant)
        # m1 (a_b + ddu1)            = -F1 - R1 + F2 + R2
        # m2 (a_b + ddu1 + ddu2)     = -F2 - R2
        ddu1 = (-F1 - R1 + F2 + R2) / m1 - a_b
        ddu2 = (-F2 - R2) / m2 - a_b - ddu1
        return [du1, du2, ddu1, ddu2]

    sol = solve_ivp(rhs, (t[0], t[-1]), np.zeros(4), method="Radau",
                    dense_output=True, rtol=1e-9, atol=1e-11,
                    max_step=(t[-1] - t[0]) / 200.0)
    u1, u2, du1, du2 = sol.sol(t)
    F1 = p.C_chv1_N_per_m * u1 + p.alpha_v1 * m1 * du1
    R1 = p.f_ch * m1 * p.g * np.tanh(du1 / p.eps_v)

    a, v = man.a(t), man.v(t)
    m_all = m0 + m1 + m2
    R_road = p.f_d * m_all * p.g * np.tanh(v / p.eps_v)
    # m0 a = F - R_road + F1 + R1
    F = m0 * a + R_road - F1 - R1

    # effective accelerated mass: inertial part of the demand divided by a_ref
    with np.errstate(divide="ignore", invalid="ignore"):
        m_eff = np.where(np.abs(a) > 1e-6, (F - R_road) / a, np.nan)
    out = _pack(t, v, a, F, p, m_eff=m_eff)
    out.update(u1=u1, u2=u2, F1=F1)
    return out


def _pack(t, v, a, F, p: TankerParams, m_eff) -> dict:
    P = F * v
    M = F * p.r_k_m
    dt = t[-1] - t[0]
    return {
        "t": t, "v": v, "a": a, "F": F, "P": P, "M": M, "m_eff": m_eff,
        "F_peak": float(np.max(F)),
        "P_peak": float(np.max(P)),
        "M_peak": float(np.max(M)),
        "P_mean": float(np.trapezoid(P, t) / dt),
        "W": float(np.trapezoid(P, t)),
        "u1_peak": float(np.max(np.abs(m_eff * 0.0))),   # replaced by caller if relevant
    }


def compare(p: TankerParams, man: StartOff, n: int = 4001, t_tail: float = 10.0) -> dict:
    """Rigid vs segmented traction demand for the same prescribed manoeuvre.

    The window is extended by ``t_tail`` beyond the ramp: after the chassis has reached
    the target speed the cargo units are still oscillating and still load the tractor,
    so truncating at t = T would flatter the segmented configuration.
    """
    t = np.linspace(0.0, man.T + t_tail, n)
    r = demand_rigid(p, man, t)
    s = demand_segmented(p, man, t)
    return {
        "rigid": r, "segmented": s,
        "dF_peak_pct": 100 * (r["F_peak"] - s["F_peak"]) / r["F_peak"],
        "dP_peak_pct": 100 * (r["P_peak"] - s["P_peak"]) / r["P_peak"],
        "dM_peak_pct": 100 * (r["M_peak"] - s["M_peak"]) / r["M_peak"],
        "dW_pct": 100 * (r["W"] - s["W"]) / r["W"],
        "u1_peak_mm": 1e3 * float(np.max(np.abs(s["u1"]))),
        "u2_peak_mm": 1e3 * float(np.max(np.abs(s["u2"]))),
    }
