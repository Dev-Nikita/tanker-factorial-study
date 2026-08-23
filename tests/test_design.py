"""Design-of-experiments integrity tests."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tanker_dynamics.config import load_all                     # noqa: E402
from tanker_dynamics.doe import factorial_design, validate      # noqa: E402
from tanker_dynamics import regression as rg                    # noqa: E402


@pytest.fixture(scope="module")
def design():
    return factorial_design(load_all(), "full_scale")


def test_27_unique_runs(design):
    validate(design)


def test_levels_map_to_physical_values(design):
    lv = design.groupby("x1_coded").front_fill_pct.unique()
    assert sorted(v[0] for v in lv) == [0.0, 50.0, 100.0]
    lv = design.groupby("x3_coded").speed_factor_mps.unique()
    assert sorted(v[0] for v in lv) == [4.0, 4.5, 5.0]


def test_design_is_balanced(design):
    for c in ("x1_coded", "x2_coded", "x3_coded"):
        assert set(design[c].value_counts()) == {9}


def test_second_order_model_has_17_residual_dof(design):
    X = rg.design_matrix(design)
    assert X.shape == (27, 9)
    assert 27 - (9 + 1) == 17


def test_saturated_model_is_saturated(design):
    X = rg.saturated_matrix(design)
    assert X.shape == (27, 26)          # + intercept = 27 coefficients
    assert 27 - 27 == 0                 # zero residual degrees of freedom


def test_randomised_order_is_a_permutation(design):
    assert sorted(design.randomized_run_id) == list(range(1, 28))
