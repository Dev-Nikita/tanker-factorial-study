"""Three-level, three-factor full factorial design (3^3 = 27), SRC-B."""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from .config import unit_masses

CODED = (-1, 0, 1)


def factorial_design(cfg: dict, scale: str = "full_scale") -> pd.DataFrame:
    """Generate the 27-run design in standard and randomised order."""
    f = cfg["experiment"]["factors"]
    lev = {k: {int(a): float(b) for a, b in f[k]["levels"].items()} for k in ("x1", "x2", "x3")}

    rows = []
    for i, (c1, c2, c3) in enumerate(itertools.product(CODED, CODED, CODED), start=1):
        lam1, lam2, v_t = lev["x1"][c1], lev["x2"][c2], lev["x3"][c3]
        m1, m2 = unit_masses(cfg, lam1, lam2, scale)
        rows.append({
            "run_id": i,
            "x1_coded": c1, "x2_coded": c2, "x3_coded": c3,
            "front_fill_pct": 100 * lam1,
            "rear_fill_pct": 100 * lam2,
            "speed_factor_mps": v_t,
            "lambda_1": lam1, "lambda_2": lam2,
            "front_unit_mass_kg": m1,
            "rear_unit_mass_kg": m2,
            "front_liquid_mass_kg": cfg["tanks"][scale]["m_liquid_full_kg"] * lam1,
            "rear_liquid_mass_kg": cfg["tanks"][scale]["m_liquid_full_kg"] * lam2,
            "scale": scale,
        })
    df = pd.DataFrame(rows)
    df["total_movable_mass_kg"] = df.front_unit_mass_kg + df.rear_unit_mass_kg
    df["load_asymmetry"] = ((df.front_unit_mass_kg - df.rear_unit_mass_kg).abs()
                            / (df.front_unit_mass_kg + df.rear_unit_mass_kg + 1e-12))

    rng = np.random.default_rng(cfg["experiment"]["randomisation"]["seed"])
    order = rng.permutation(len(df)) + 1
    df["randomized_run_id"] = order
    return df


def validate(df: pd.DataFrame) -> None:
    combos = set(zip(df.x1_coded, df.x2_coded, df.x3_coded))
    assert len(df) == 27, f"expected 27 runs, got {len(df)}"
    assert len(combos) == 27, "duplicated factor combination"
    for c in ("x1_coded", "x2_coded", "x3_coded"):
        assert set(df[c].unique()) == {-1, 0, 1}, f"{c} levels are not -1/0/+1"
    assert sorted(df.randomized_run_id) == list(range(1, 28)), "bad randomisation"
