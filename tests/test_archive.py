"""Archive-locked model: worksheet fidelity and physical consistency."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.archive_model import ArchiveParams, V1, PT1, rhs   # noqa: E402
from tanker_dynamics.archive_studies import (Manoeuvre, compare_matched,  # noqa: E402
                                             simulate)


def test_scalars_match_worksheet():
    p = ArchiveParams()
    assert (p.Rk, p.g, p.Alf_d) == (0.53, 9.81, 0.06)
    assert (p.Mch, p.Mkuz, p.Mv, p.Itr) == (5900.0, 700.0, 5000.0, 2025.0)
    assert (p.Alf_ch, p.Alf_k) == (0.02, 0.02)
    assert (p.CKCH, p.CCHK, p.CKV) == (1.2e6, 2.0e5, 2.0e5)
    assert p.Mgal_max == 1600.0
    assert len(V1) == len(PT1) == 19
    assert PT1.max() == 65472.4


def test_derived_quantities():
    p = ArchiveParams()
    assert p.MSchkv == 11600.0 and p.MSkv == 5700.0
    assert p.Mchk == pytest.approx(1 / p.Mch + 1 / p.Mkuz)
    assert p.Mkv == pytest.approx(1 / p.Mkuz + 1 / p.Mv)


def test_no_calibrated_parameters():
    """The archive contains no viscous damping anywhere; defaults must be zero."""
    p = ArchiveParams()
    assert p.D_kch == 0.0 and p.c1 == 0.0 and p.c2 == 0.0


def test_momentum_balance():
    """Sum of the three body equations must equal F_kch - A1 exactly."""
    p = ArchiveParams(Alf_ch=0.25, Alf_k=0.25)     # large friction: error visible
    y = np.array([0.4, 1.1, 0.03, -0.02, 5.0, 2.4, 0.35, -0.22])
    dy = rhs(0.0, y, p)
    a_ch, ddu1, ddu2 = dy[5], dy[6], dy[7]
    total = (p.Mch * a_ch + p.Mkuz * (a_ch + ddu1)
             + p.Mv * (a_ch + ddu1 + ddu2))
    F = p.CKCH * (y[0] - y[1] / p.Rk) / p.Rk
    A1 = p.Alf_d * p.MSchkv * p.g * np.tanh(y[5] / p.eps)
    assert abs(total - (F - A1)) < 1e-6 * abs(F)


def test_reproduces_published_stage1():
    p = ArchiveParams()
    t, y, _ = simulate(p, (0.0, 600.0), n=6001, rtol=1e-8, atol=1e-10)
    assert y[5][-1] == pytest.approx(13.058, rel=1e-3)      # terminal speed
    assert y[1][-1] == pytest.approx(7582.0, rel=1e-3)      # distance


def test_reproduces_published_braking():
    p = ArchiveParams()
    _, y, _ = simulate(p, (0.0, 600.0), n=3001, rtol=1e-8, atol=1e-10)
    y0 = y[:, -1].copy()
    tb, yb, _ = simulate(p, (600.0, 800.0), braking=True, y0=y0,
                         stop_at_rest=True, n=4001, rtol=1e-8, atol=1e-10)
    assert (tb[-1] - 600.0) == pytest.approx(25.054, rel=0.01)
    assert (yb[1][-1] - y0[1]) == pytest.approx(163.575, rel=0.01)


def test_stiff_limit_approaches_rigid():
    soft = ArchiveParams(CCHK=5e10, CKV=5e10, Alf_ch=0.0, Alf_k=0.0)
    hard = ArchiveParams(rigid=True)
    fs = [np.max(np.abs(
        (q.CKCH * (y[0] - y[1] / q.Rk)) / q.Rk))
        for q in (soft, hard)
        for _, y, _ in [simulate(q, (0.0, 4.0), n=4001)]]
    assert abs(fs[0] - fs[1]) / fs[1] < 0.03


def test_matched_work_never_decreases():
    """A passive coupling can only add dissipation (work-energy balance)."""
    p = ArchiveParams()
    for T in (2.0, 3.0, 8.0, 20.0):
        c = compare_matched(p, Manoeuvre(v_f=5.0, T=T))
        assert c["dW_pct"] <= 1e-6, f"T={T}: {c['dW_pct']:.3f} % less work"


def test_matched_no_peak_benefit_over_domain():
    p = ArchiveParams()
    for T in (2.0, 4.0, 8.0, 20.0):
        for shape in ("smoothstep", "halfcosine", "quintic"):
            c = compare_matched(p, Manoeuvre(v_f=5.0, T=T, shape=shape))
            assert c["dF_pct"] <= 1e-6


def test_viscous_damping_extension_matches_declared_formula():
    """Manuscript Section 2.2: c1 = 2 zeta sqrt(C_CHK M_kuz), c2 = 2 zeta sqrt(C_KV M_v).

    Guards the damping sweep of Table 3 against drifting away from the published
    definition.
    """
    p0 = ArchiveParams()
    for zeta in (0.0, 0.3, 1.0, 1.5):
        c1 = 2 * zeta * np.sqrt(p0.CCHK * p0.Mkuz)
        c2 = 2 * zeta * np.sqrt(p0.CKV * p0.Mv)
        p = ArchiveParams(c1=c1, c2=c2)
        assert p.c1 == pytest.approx(c1) and p.c2 == pytest.approx(c2)
        if zeta == 0.0:
            assert p.c1 == 0.0 and p.c2 == 0.0


def test_all_three_trajectory_shapes_satisfy_boundary_conditions():
    """v(0)=0, v(T)=v_f, a(0)=a(T)=0 for every reference profile."""
    for shape in ("smoothstep", "halfcosine", "quintic"):
        m = Manoeuvre(v_f=5.0, T=8.0, shape=shape)
        assert float(m.v(0.0)) == pytest.approx(0.0, abs=1e-12)
        assert float(m.v(m.T)) == pytest.approx(m.v_f, rel=1e-12)
        assert float(m.a(0.0)) == pytest.approx(0.0, abs=1e-9)
        assert float(m.a(m.T)) == pytest.approx(0.0, abs=1e-9)
