"""Archive-locked model: a direct transcription of the ETAL Mathcad worksheet.

Source: `10 тест зад (РОЗГІН- РУХ-ГАЛЬМУВ) 23.10.09р. ETAL..xmcd`, extracted
programmatically by `tools/parse_mathcad_archive.py`. Every symbol below appears in
that worksheet; nothing is inferred. Region ids are given so each term can be checked.

State (worksheet indexing y[1..8], here 0-based):

    y0 = phi_k    driving-wheel angle                    [rad]
    y1 = X_ch     chassis longitudinal coordinate        [m]
    y2 = u1       body (kuzov) coord. relative to chassis[m]
    y3 = u2       cargo coord. relative to the body      [m]
    y4 = omega    wheel angular velocity                 [rad/s]
    y5 = V_ch     chassis velocity                       [m/s]
    y6 = du1, y7 = du2

Worksheet definitions transcribed verbatim:

    MprKCH(y) = CKCH*(y1 - y2/Rk)                                  region 4078
    FprCHK(y) = CCHK*((XCH1 - y3) - Bk - LCHK)                     region 4079
    FprKV(y)  = CKV *((XK2  - y4) - Bv - LKV)                      region 4080
    A1(y) = Alf_d *MSchkv*g*sign(y6)                               region 4066
    A2(y) = Alf_ch*MSkv  *g*sign(y7)                               region 4067
    A3(y) = Alf_k *Mv    *g*sign(y8)                               region 4068
    B1(y) = MprKCH(y)/Rk - A1(y)                                   region 4081
    B2(y) = FprCHK(y) - A2(y)                                      region 4082
    B3(y) = FprKV(y)  - A3(y)                                      region 4083
    BZ2(y) = (B1 - B2)/Mch                                         region 4084
    BZ3(y) = -B1/Mch + Mchk*B2 - B3/Mkuz                           region 4085
    BZ4(y) = -B2/Mkuz + Mkv*B3                                     region 4086
    D(t,y) = (y5, y6, y7, y8, (Mtag(y) - MprKCH(y))/Itr, BZ2, BZ3, BZ4)   region 1725
    D2(t,y) = same with Mgalm(t,y) in place of Mtag(y)             region 4406
    Mchk = (Mch+Mkuz)/(Mch*Mkuz)   MSchkv = Mch+Mkuz+Mv
    Mkv  = (Mkuz+Mv)/(Mkuz*Mv)     MSkv   = Mkuz+Mv
    Mgal(t)  = atan((t - T1_kin)*10)/1.55 * Mgal_max               region 3804
    Mgalm(t,y) = Mgal(t)*sign(y5)                                  region 3805
    Mtag(y) = interp(cspline(om1,mt1), om1, mt1, y5)               region 2106
    om1 = v1/Rk,  mt1 = pt1*Rk                                     regions 2103-2104

Note on BZ3/BZ4: with Mchk = 1/Mch + 1/Mkuz and Mkv = 1/Mkuz + 1/Mv these are exactly
the relative accelerations of a three-body chain, i.e.

    Mch *a_ch          = B1 - B2
    Mkuz*(a_ch + ddu1) = B2 - B3
    Mv  *(a_ch + ddu1 + ddu2) = B3

Summing gives B1 = MprKCH/Rk - A1, so the guide-friction terms A2 and A3 cancel
internally. The worksheet therefore confirms that the friction reaction acts on the
supporting body as well; this is checked by `tests/test_archive.py`.

The worksheet contains **no viscous damping** anywhere: the couplings are purely
elastic (CCHK, CKV) and the only dissipation is Coulomb friction (Alf_*). In
particular there is no tangential wheel damping in MprKCH.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import CubicSpline

# --------------------------------------------------------------------------- #
#  Traction characteristic, ETAL regions 2103-2104 (v1 [m/s] -> pt1 [N])
# --------------------------------------------------------------------------- #
V1 = np.array([0.0, 0.1, 0.2, 0.2932, 0.4202, 0.5473, 0.6743, 0.8014, 0.9284,
               1.0555, 1.1825, 1.3059, 1.4366, 1.5637, 2.8375, 5.0803, 7.9142,
               11.6338, 20.0])
PT1 = np.array([28000.0, 33500.0, 39000.0, 44165.9, 50622.7, 55917.3, 60049.0,
                63019.4, 64827.1, 65472.4, 64955.5, 63276.2, 60434.6, 56410.8,
                31097.6, 17369.2, 11149.6, 7584.78, 0.0])


@dataclass
class ArchiveParams:
    """ETAL worksheet scalars. Defaults are the archived values; do not edit."""

    Rk: float = 0.53           # region 2102
    g: float = 9.81            # region 2107
    Alf_d: float = 0.06        # region 2108
    Mch: float = 5900.0        # region 2109
    Mkuz: float = 700.0        # region 2110
    Mv: float = 5000.0         # region 2111
    Itr: float = 2025.0        # region 2112
    Alf_ch: float = 0.02       # region 2143
    Alf_k: float = 0.02        # region 4065
    CKCH: float = 1.2e6        # region 4069
    CCHK: float = 2.0e5        # region 4070
    CKV: float = 2.0e5         # region 4071
    Mgal_max: float = 1600.0   # region 3803
    T1_kin: float = 600.0      # region 2676

    # optional extensions, zero in the archive
    D_kch: float = 0.0         # tangential wheel damping - NOT in the worksheet
    c1: float = 0.0            # viscous damping chassis-body    - NOT in worksheet
    c2: float = 0.0            # viscous damping body-cargo      - NOT in worksheet

    eps: float = 1.0e-3        # sign() regularisation
    rigid: bool = False        # lock the couplings (baseline configuration)

    _spline: CubicSpline = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        om1, mt1 = V1 / self.Rk, PT1 * self.Rk
        self._spline = CubicSpline(om1, mt1)
        self._om_max = float(om1[-1])

    # ---- derived, ETAL regions 2113, 2140-2142 ------------------------------ #
    @property
    def MSchkv(self) -> float:
        return self.Mch + self.Mkuz + self.Mv

    @property
    def MSkv(self) -> float:
        return self.Mkuz + self.Mv

    @property
    def Mchk(self) -> float:
        return (self.Mch + self.Mkuz) / (self.Mch * self.Mkuz)

    @property
    def Mkv(self) -> float:
        return (self.Mkuz + self.Mv) / (self.Mkuz * self.Mv)

    def Mtag(self, omega: float) -> float:
        """Wheel-referred traction torque, ETAL region 2106."""
        if omega <= 0.0:
            return float(PT1[0] * self.Rk)
        if omega >= self._om_max:
            return 0.0
        return float(self._spline(omega))

    def Mgal(self, t: float) -> float:
        """Brake torque, ETAL region 3804."""
        return float(np.arctan((t - self.T1_kin) * 10.0) / 1.55 * self.Mgal_max)


def _sgn(v: float, eps: float) -> float:
    return float(np.tanh(v / eps))


# --------------------------------------------------------------------------- #
def rhs(t: float, y: np.ndarray, p: ArchiveParams, braking: bool = False) -> np.ndarray:
    phi, X, u1, u2, om, V, du1, du2 = y

    MprKCH = p.CKCH * (phi - X / p.Rk) + p.D_kch * (om - V / p.Rk)
    if p.rigid:
        FprCHK = FprKV = 0.0
    else:
        FprCHK = -p.CCHK * u1 - p.c1 * du1        # XCH1 - Bk - LCHK == 0
        FprKV = -p.CKV * u2 - p.c2 * du2          # XK2  - Bv - LKV  == 0

    A1 = p.Alf_d * p.MSchkv * p.g * _sgn(V, p.eps)
    A2 = p.Alf_ch * p.MSkv * p.g * _sgn(du1, p.eps)
    A3 = p.Alf_k * p.Mv * p.g * _sgn(du2, p.eps)

    B1 = MprKCH / p.Rk - A1
    B2 = FprCHK - A2
    B3 = FprKV - A3

    # ETAL region 1725 (traction) / 4406 (braking, Mgalm = Mgal(t)*sign(omega))
    M_wheel = -p.Mgal(t) * _sgn(om, p.eps / p.Rk) if braking else p.Mtag(om)
    dom = (M_wheel - MprKCH) / p.Itr

    if p.rigid:
        # single lumped body: Mch+Mkuz+Mv accelerated together
        a_ch = (MprKCH / p.Rk - A1) / p.MSchkv
        return np.array([om, V, 0.0, 0.0, dom, a_ch, 0.0, 0.0])

    BZ2 = (B1 - B2) / p.Mch
    BZ3 = -B1 / p.Mch + p.Mchk * B2 - B3 / p.Mkuz
    BZ4 = -B2 / p.Mkuz + p.Mkv * B3
    return np.array([om, V, du1, du2, dom, BZ2, BZ3, BZ4])


def wheel_force(y: np.ndarray, p: ArchiveParams) -> np.ndarray:
    """F_kch = MprKCH / Rk, the longitudinal force from the wheels on the chassis."""
    phi, X, _, _, om, V, _, _ = y
    return (p.CKCH * (phi - X / p.Rk) + p.D_kch * (om - V / p.Rk)) / p.Rk
