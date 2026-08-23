"""Longitudinal dynamic model of a tank semitrailer with compliant-damped cargo mounting.

Implements the two computational schemes of the source manuscript (SRC-A):

  * ``rigid``     - eq. (2): "tractor engine - tanker chassis together with the cargo".
                    Two generalised coordinates: driving-wheel angle ``phi`` and
                    absolute chassis coordinate ``X_ch``.
  * ``compliant`` - eq. (1): "engine - chassis - series-connected tanks (cargo)".
                    Four generalised coordinates: ``phi``, ``X_ch`` and the two
                    cargo-unit coordinates measured in the chassis-fixed frame
                    O1X1Z1 (see SRC-A Fig. 2).

Sign convention: motion in the direction of increasing X0.
Coulomb friction terms use a tanh regularisation so that the right-hand side is
Lipschitz and the stiff Radau integrator behaves.

See docs/EQUATION_TRACEABILITY.md for the mapping onto the source equations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


# --------------------------------------------------------------------------- #
#  Powertrain
# --------------------------------------------------------------------------- #
@dataclass
class Powertrain:
    """Wheel-referred engine torque characteristic, SRC-A Fig. 4 and eq. (3)."""

    M0_Nm: float = 1.50e4
    omega_peak_rad_s: float = 2.6
    M_peak_Nm: float = 3.55e4
    omega_ideal_end_rad_s: float = 22.0
    omega_max_rad_s: float = 37.0

    @property
    def N_d_W(self) -> float:
        """Plateau power of the 'ideal' region, N_d = M_peak * omega_peak."""
        return self.M_peak_Nm * self.omega_peak_rad_s

    def torque(self, omega: float) -> float:
        w = max(omega, 0.0)
        wp, we, wm = self.omega_peak_rad_s, self.omega_ideal_end_rad_s, self.omega_max_rad_s
        if w <= wp:
            # rising branch, monotone C1 blend from M0 to M_peak
            s = w / wp
            return self.M0_Nm + (self.M_peak_Nm - self.M0_Nm) * (3 * s**2 - 2 * s**3)
        if w <= we:
            return self.N_d_W / w                      # constant-power ("ideal") region
        M_e = self.N_d_W / we
        if w >= wm:
            return 0.0
        return M_e * (wm - w) / (wm - we)              # linear cut-off branch

    def power(self, omega: float) -> float:
        return self.torque(omega) * max(omega, 0.0)


@dataclass
class Brakes:
    """Total brake torque law applied to the wheels, SRC-A Fig. 5."""

    M_inf_Nm: float = 1750.0
    M_drop_Nm: float = 1200.0
    tau_s: float = 0.8
    t_start_s: float = 600.0
    scale: float = 1.0        # number of braked wheel sets (see MISSING_INPUTS P-5)

    def torque(self, t: float) -> float:
        if t < self.t_start_s:
            return 0.0
        return self.scale * (self.M_inf_Nm
                             - self.M_drop_Nm * np.exp(-(t - self.t_start_s) / self.tau_s))


# --------------------------------------------------------------------------- #
#  Parameter container
# --------------------------------------------------------------------------- #
@dataclass
class TankerParams:
    # inertias
    m_ch_kg: float = 5900.0
    m_tractor_kg: float = 6700.0
    m_unit_kg: tuple = (5700.0, 5700.0)      # trolley + tank + liquid, per movable unit
    I_k_kgm2: float = 190.0
    r_k_m: float = 0.540

    # couplings
    C_kch_N_per_rad: float = 1.2e6
    D_kch_Nms_per_rad: float = 0.0            # tangential wheel damping (see EQUATION_TRACEABILITY M-04)
    C_chv1_N_per_m: float = 2.0e5
    C_v1v2_N_per_m: float = 1.5e5
    alpha_v1: float = 0.10                   # N s / m per kg of the coupled unit
    alpha_v2: float = 0.05

    # resistances
    f_d: float = 0.03                        # road rolling resistance
    f_ch: float = 0.02                       # unit 1 along the frame guides
    f_k: float = 0.02                        # unit 2 along the guides
    g: float = 9.81

    # numerics
    eps_v: float = 1.0e-3

    powertrain: Powertrain = field(default_factory=Powertrain)
    brakes: Brakes = field(default_factory=Brakes)

    @property
    def m_total_kg(self) -> float:
        return self.m_tractor_kg + self.m_ch_kg + sum(self.m_unit_kg)


def _sgn(v: float, eps: float) -> float:
    return float(np.tanh(v / eps))


# --------------------------------------------------------------------------- #
#  Right-hand sides
# --------------------------------------------------------------------------- #
def rhs_rigid(t: float, y: np.ndarray, p: TankerParams, braking: bool = False) -> np.ndarray:
    """SRC-A eq. (2). State y = [phi, X_ch, dphi, dX_ch]."""
    phi, X_ch, dphi, dX_ch = y
    m_c = p.m_total_kg

    twist = phi - X_ch / p.r_k_m
    dtwist = dphi - dX_ch / p.r_k_m
    M_tyre = p.C_kch_N_per_rad * twist + p.D_kch_Nms_per_rad * dtwist
    F_kch = M_tyre / p.r_k_m                             # longitudinal wheel force on chassis

    M_drive = 0.0 if braking else p.powertrain.torque(dphi)
    M_brake = p.brakes.torque(t) if braking else 0.0

    ddphi = (M_drive - M_brake * _sgn(dphi, p.eps_v / p.r_k_m) - M_tyre) / p.I_k_kgm2
    ddX = (F_kch - p.f_d * m_c * p.g * _sgn(dX_ch, p.eps_v)) / m_c
    return np.array([dphi, dX_ch, ddphi, ddX])


def rhs_compliant(t: float, y: np.ndarray, p: TankerParams, braking: bool = False) -> np.ndarray:
    """SRC-A eq. (1). State y = [phi, X_ch, u1, u2, dphi, dX_ch, du1, du2].

    ``u1`` is the longitudinal displacement of cargo unit 1 relative to the chassis
    (measured in O1X1Z1, spring free length removed); ``u2`` that of unit 2
    relative to unit 1.  Absolute accelerations are recovered inside the routine.
    """
    phi, X_ch, u1, u2, dphi, dX_ch, du1, du2 = y
    m1, m2 = p.m_unit_kg
    m_ch_eff = p.m_tractor_kg + p.m_ch_kg

    twist = phi - X_ch / p.r_k_m
    dtwist = dphi - dX_ch / p.r_k_m
    M_tyre = p.C_kch_N_per_rad * twist + p.D_kch_Nms_per_rad * dtwist
    F_kch = M_tyre / p.r_k_m

    # coupling forces (positive = stretching the element)
    F1 = p.C_chv1_N_per_m * u1 + p.alpha_v1 * m1 * du1
    F2 = p.C_v1v2_N_per_m * u2 + p.alpha_v2 * m2 * du2

    # guide friction, proportional to the normal load of the mass carried
    R1 = p.f_ch * m1 * p.g * _sgn(du1, p.eps_v)
    R2 = p.f_k * m2 * p.g * _sgn(du2, p.eps_v)

    M_drive = 0.0 if braking else p.powertrain.torque(dphi)
    M_brake = p.brakes.torque(t) if braking else 0.0
    ddphi = (M_drive - M_brake * _sgn(dphi, p.eps_v / p.r_k_m) - M_tyre) / p.I_k_kgm2

    m_all = m_ch_eff + m1 + m2
    R_road = p.f_d * m_all * p.g * _sgn(dX_ch, p.eps_v)

    # Absolute accelerations: a_ch, a1 = a_ch + ddu1, a2 = a_ch + ddu1 + ddu2.
    # Guide friction is an INTERNAL force pair: -R1 on unit 1, +R1 on the frame;
    # -R2 on unit 2, +R2 on unit 1.  Summing the three equations then leaves only
    # F_kch - R_road, as it must (see tests/test_physics.py::test_internal_forces_cancel).
    #   m_ch_eff a_ch          = F_kch - R_road + F1 + R1
    #   m1 (a_ch + ddu1)       = -F1 - R1 + F2 + R2
    #   m2 (a_ch + ddu1 + ddu2)= -F2 - R2
    A = np.array([
        [m_ch_eff, 0.0, 0.0],
        [m1, m1, 0.0],
        [m2, m2, m2],
    ])
    b = np.array([
        F_kch - R_road + F1 + R1,
        -F1 - R1 + F2 + R2,
        -F2 - R2,
    ])
    a_ch, ddu1, ddu2 = np.linalg.solve(A, b)
    return np.array([dphi, dX_ch, du1, du2, ddphi, a_ch, ddu1, ddu2])
