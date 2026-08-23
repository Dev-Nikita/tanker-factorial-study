"""Forward and matched-manoeuvre studies on the archive-locked ETAL model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .archive_model import ArchiveParams, rhs, wheel_force, _sgn


# --------------------------------------------------------------------------- #
#  Forward integration
# --------------------------------------------------------------------------- #
def simulate(p: ArchiveParams, t_span=(0.0, 8.0), braking=False, y0=None,
             n=8001, target_speed=None, stop_at_rest=False,
             rtol=1e-9, atol=1e-11):
    y0 = np.zeros(8) if y0 is None else np.asarray(y0, float)
    events = []
    if target_speed is not None:
        def reach(t, y, *a):
            return y[5] - target_speed
        reach.terminal = True
        reach.direction = 1.0
        events.append(reach)
    if stop_at_rest:
        def rest(t, y, *a):
            return y[5] - 1.0e-3
        rest.terminal = True
        rest.direction = -1.0
        events.append(rest)
    sol = solve_ivp(rhs, t_span, y0, method="Radau", args=(p, braking),
                    rtol=rtol, atol=atol, dense_output=True,
                    events=events or None,
                    max_step=(t_span[1] - t_span[0]) / 200.0)
    t = np.linspace(t_span[0], sol.t[-1], n)
    return t, sol.sol(t), sol


def forward_metrics(p: ArchiveParams, t, y) -> dict:
    F = wheel_force(y, p)
    V = y[5]
    a = np.gradient(V, t)
    P = F * V
    dt = t[-1] - t[0]
    return {
        "F_peak_N": float(np.max(np.abs(F))),
        "a_peak_m_s2": float(np.max(np.abs(a))),
        "P_peak_W": float(np.max(P)),
        "P_mean_W": float(np.trapezoid(P, t) / dt),
        "E_J": float(np.trapezoid(P, t)),
        "u1_peak_mm": float(1e3 * np.max(np.abs(y[2]))),
        "u2_peak_mm": float(1e3 * np.max(np.abs(y[3]))),
        "t_end_s": float(t[-1]),
        "distance_m": float(y[1][-1] - y[1][0]),
        "v_end_m_s": float(V[-1]),
    }


# --------------------------------------------------------------------------- #
#  Matched manoeuvre (inverse dynamics)
# --------------------------------------------------------------------------- #
@dataclass
class Manoeuvre:
    """Prescribed chassis motion from rest to v_f in T seconds."""

    v_f: float = 5.0
    T: float = 8.0
    shape: str = "smoothstep"      # smoothstep | halfcosine | quintic

    def v(self, t):
        s = np.clip(np.asarray(t, float) / self.T, 0.0, 1.0)
        if self.shape == "smoothstep":
            return self.v_f * s ** 2 * (3 - 2 * s)
        if self.shape == "halfcosine":
            return self.v_f * 0.5 * (1 - np.cos(np.pi * s))
        if self.shape == "quintic":
            return self.v_f * s ** 3 * (10 - 15 * s + 6 * s ** 2)
        raise ValueError(self.shape)

    def a(self, t):
        s = np.clip(np.asarray(t, float) / self.T, 0.0, 1.0)
        if self.shape == "smoothstep":
            return 6 * self.v_f * s * (1 - s) / self.T
        if self.shape == "halfcosine":
            return self.v_f * 0.5 * np.pi / self.T * np.sin(np.pi * s)
        if self.shape == "quintic":
            return self.v_f * 30 * s ** 2 * (1 - s) ** 2 / self.T
        raise ValueError(self.shape)


def demand(p: ArchiveParams, man: Manoeuvre, t: np.ndarray, segmented: bool) -> dict:
    """Traction force required at the wheel-chassis interface for the prescribed motion."""
    a_ref, v_ref = man.a(t), man.v(t)
    A1 = p.Alf_d * p.MSchkv * p.g * np.tanh(v_ref / p.eps)

    if not segmented:
        F = p.MSchkv * a_ref + A1
        m_eff = np.full_like(t, p.MSchkv)
        u1 = u2 = np.zeros_like(t)
    else:
        def f(tt, z):
            u1, u2, du1, du2 = z
            ab = float(man.a(tt))
            FprCHK = -p.CCHK * u1 - p.c1 * du1
            FprKV = -p.CKV * u2 - p.c2 * du2
            A2 = p.Alf_ch * p.MSkv * p.g * _sgn(du1, p.eps)
            A3 = p.Alf_k * p.Mv * p.g * _sgn(du2, p.eps)
            B2, B3 = FprCHK - A2, FprKV - A3
            ddu1 = B2 / p.Mkuz - B3 / p.Mkuz - ab
            ddu2 = B3 / p.Mv - ab - ddu1
            return [du1, du2, ddu1, ddu2]

        s = solve_ivp(f, (t[0], t[-1]), np.zeros(4), method="Radau",
                      dense_output=True, rtol=1e-9, atol=1e-11,
                      max_step=(t[-1] - t[0]) / 400.0)
        u1, u2, du1, du2 = s.sol(t)
        FprCHK = -p.CCHK * u1 - p.c1 * du1
        A2 = p.Alf_ch * p.MSkv * p.g * np.tanh(du1 / p.eps)
        # Mch*a_ch = B1 - B2 with B1 = F - A1, B2 = FprCHK - A2
        #   ->  F = Mch*a_ch + A1 + FprCHK - A2
        F = p.Mch * a_ref + A1 + FprCHK - A2
        with np.errstate(divide="ignore", invalid="ignore"):
            m_eff = np.where(np.abs(a_ref) > 0.05, (F - A1) / a_ref, np.nan)

    P = F * v_ref
    dt = t[-1] - t[0]
    return {"t": t, "v": v_ref, "a": a_ref, "F": F, "P": P, "m_eff": m_eff,
            "u1": u1, "u2": u2,
            "F_peak": float(np.max(F)), "P_peak": float(np.max(P)),
            "M_peak": float(np.max(F) * p.Rk),
            "W": float(np.trapezoid(P, t)), "P_mean": float(np.trapezoid(P, t) / dt)}


def compare_matched(p: ArchiveParams, man: Manoeuvre, tail: float = 12.0,
                    n: int = 6001) -> dict:
    t = np.linspace(0.0, man.T + tail, n)
    r = demand(p, man, t, segmented=False)
    s = demand(p, man, t, segmented=True)
    return {
        "rigid": r, "segmented": s,
        "dF_pct": 100 * (r["F_peak"] - s["F_peak"]) / r["F_peak"],
        "dP_pct": 100 * (r["P_peak"] - s["P_peak"]) / r["P_peak"],
        "dW_pct": 100 * (r["W"] - s["W"]) / r["W"],
        "u1_mm": 1e3 * float(np.max(np.abs(s["u1"]))),
        "u2_mm": 1e3 * float(np.max(np.abs(s["u2"]))),
        "m_eff_min_kg": float(np.nanmin(s["m_eff"])),
    }


def coupling_period(p: ArchiveParams) -> float:
    """Lower coupling mode of the free chain body-cargo (the slower of the two)."""
    return 2.0 * np.pi * np.sqrt(p.Mv / p.CKV)
