"""Numerical integration and response-metric extraction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from .model import TankerParams, rhs_rigid, rhs_compliant


@dataclass
class Solution:
    t: np.ndarray
    y: np.ndarray
    config: str
    params: TankerParams
    ok: bool = True
    message: str = ""

    # ---- derived channels ------------------------------------------------ #
    @property
    def phi(self):
        return self.y[0]

    @property
    def omega(self):
        return self.y[2] if self.config == "rigid" else self.y[4]

    @property
    def X_ch(self):
        return self.y[1]

    @property
    def v_ch(self):
        return self.y[3] if self.config == "rigid" else self.y[5]

    def wheel_force(self):
        p = self.params
        twist = self.phi - self.X_ch / p.r_k_m
        dtwist = self.omega - self.v_ch / p.r_k_m
        return (p.C_kch_N_per_rad * twist + p.D_kch_Nms_per_rad * dtwist) / p.r_k_m

    def engine_torque(self):
        return np.array([self.params.powertrain.torque(w) for w in self.omega])

    def engine_power(self):
        return self.engine_torque() * np.maximum(self.omega, 0.0)

    def accel(self):
        return np.gradient(self.v_ch, self.t)

    def jerk(self):
        return np.gradient(self.accel(), self.t)

    def relative(self):
        """(u1, u2, du1, du2) - zeros for the rigid baseline."""
        if self.config == "rigid":
            z = np.zeros_like(self.t)
            return z, z, z, z
        return self.y[2], self.y[3], self.y[6], self.y[7]

    def coupling_forces(self):
        p = self.params
        u1, u2, du1, du2 = self.relative()
        m1, m2 = p.m_unit_kg
        F1 = p.C_chv1_N_per_m * u1 + p.alpha_v1 * m1 * du1
        F2 = p.C_v1v2_N_per_m * u2 + p.alpha_v2 * m2 * du2
        return F1, F2


def simulate(params: TankerParams,
             config: str = "compliant",
             t_span=(0.0, 60.0),
             braking: bool = False,
             y0: Optional[np.ndarray] = None,
             target_speed: Optional[float] = None,
             stop_at_rest: bool = False,
             n_out: int = 4001,
             rtol: float = 1e-8,
             atol: float = 1e-10) -> Solution:
    """Integrate one manoeuvre.

    If ``target_speed`` is given, integration terminates when the chassis first
    reaches it (used for the start-off scenario of the 3^3 design).
    """
    if config == "rigid":
        f, n = rhs_rigid, 4
    elif config == "compliant":
        f, n = rhs_compliant, 8
    else:
        raise ValueError(f"unknown configuration {config!r}")

    if y0 is None:
        y0 = np.zeros(n)
    y0 = np.asarray(y0, dtype=float)
    if y0.size != n:
        raise ValueError(f"state size {y0.size} != {n} for config {config!r}")

    iv = 3 if config == "rigid" else 5
    events = []
    if target_speed is not None:
        def reach(t, y, *a):
            return y[iv] - target_speed
        reach.terminal = True
        reach.direction = 1.0
        events.append(reach)
    if stop_at_rest:
        def rest(t, y, *a):
            return y[iv] - 1.0e-3
        rest.terminal = True
        rest.direction = -1.0
        events.append(rest)
    events = events or None

    sol = solve_ivp(f, t_span, y0, method="Radau", args=(params, braking),
                    rtol=rtol, atol=atol, dense_output=True, events=events,
                    max_step=(t_span[1] - t_span[0]) / 50.0)

    t_end = sol.t[-1]
    t = np.linspace(t_span[0], t_end, n_out)
    y = sol.sol(t)
    ok = sol.success and np.all(np.isfinite(y))
    return Solution(t=t, y=y, config=config, params=params, ok=ok, message=sol.message)


# --------------------------------------------------------------------------- #
#  Response metrics
# --------------------------------------------------------------------------- #
def responses(sol: Solution) -> dict:
    """Manuscript response variables Y1...Y12."""
    t = sol.t
    P = sol.engine_power()
    M = sol.engine_torque()
    F_kch = sol.wheel_force()
    F1, F2 = sol.coupling_forces()
    a = sol.accel()
    j = sol.jerk()
    u1, u2, du1, du2 = sol.relative()
    dt = t[-1] - t[0]

    def rms(x):
        return float(np.sqrt(np.mean(np.asarray(x) ** 2)))

    out = {
        # powertrain
        "P_peak_W": float(np.max(P)),
        "P_mean_W": float(np.trapezoid(P, t) / dt),
        "E_engine_J": float(np.trapezoid(P, t)),
        "M_peak_Nm": float(np.max(M)),
        "M_mean_Nm": float(np.trapezoid(M, t) / dt),
        # mechanical loads
        "F_kch_peak_N": float(np.max(np.abs(F_kch))),
        "F_chv1_peak_N": float(np.max(np.abs(F1))),
        "F_v1v2_peak_N": float(np.max(np.abs(F2))),
        # kinematics
        "a_peak_m_s2": float(np.max(np.abs(a))),
        "a_rms_m_s2": rms(a),
        "j_peak_m_s3": float(np.max(np.abs(j))),
        "j_rms_m_s3": rms(j),
        # relative motion
        "u1_peak_m": float(np.max(np.abs(u1))),
        "u2_peak_m": float(np.max(np.abs(u2))),
        "du1_peak_m_s": float(np.max(np.abs(du1))),
        "du2_peak_m_s": float(np.max(np.abs(du2))),
        # global
        "t_manoeuvre_s": float(dt),
        "distance_m": float(sol.X_ch[-1] - sol.X_ch[0]),
        "v_end_m_s": float(sol.v_ch[-1]),
        "v_mean_m_s": float(np.trapezoid(sol.v_ch, t) / dt),
    }
    out["settling_time_s"] = _settling_time(t, u1)
    return out


def _settling_time(t, u, band: float = 0.05) -> float:
    """Time after which |u| stays within ``band`` of its peak. NaN if never settles."""
    u = np.abs(np.asarray(u))
    pk = u.max()
    if pk <= 0:
        return 0.0
    thr = band * pk
    idx = np.where(u > thr)[0]
    if idx.size == 0:
        return 0.0
    return float(t[idx[-1]] - t[0])
