"""Response-surface regression and ANOVA for the 3^3 design.

Primary inferential model (10 coefficients, 27 observations, 17 residual d.o.f.):

    Y = b0 + b1 x1 + b2 x2 + b3 x3
          + b12 x1x2 + b13 x1x3 + b23 x2x3
          + b11 x1^2 + b22 x2^2 + b33 x3^2 + eps

The saturated 27-term tensor polynomial is available for completeness only and is
explicitly flagged as non-inferential (zero residual degrees of freedom).
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import statsmodels.api as sm

TERMS2 = ["x1", "x2", "x3", "x1x2", "x1x3", "x2x3", "x1^2", "x2^2", "x3^2"]
TERMS1 = ["x1", "x2", "x3", "x1x2", "x1x3", "x2x3"]


def design_matrix(df: pd.DataFrame, terms=TERMS2) -> pd.DataFrame:
    x1, x2, x3 = df.x1_coded.to_numpy(float), df.x2_coded.to_numpy(float), df.x3_coded.to_numpy(float)
    cols = {
        "x1": x1, "x2": x2, "x3": x3,
        "x1x2": x1 * x2, "x1x3": x1 * x3, "x2x3": x2 * x3,
        "x1^2": x1 ** 2, "x2^2": x2 ** 2, "x3^2": x3 ** 2,
    }
    return pd.DataFrame({t: cols[t] for t in terms}, index=df.index)


def saturated_matrix(df: pd.DataFrame) -> pd.DataFrame:
    x = [df.x1_coded.to_numpy(float), df.x2_coded.to_numpy(float), df.x3_coded.to_numpy(float)]
    cols = {}
    for i, j, k in itertools.product(range(3), repeat=3):
        if i == j == k == 0:
            continue
        cols[f"x1^{i}x2^{j}x3^{k}"] = x[0] ** i * x[1] ** j * x[2] ** k
    return pd.DataFrame(cols, index=df.index)


def fit(df: pd.DataFrame, y: str, terms=TERMS2):
    X = sm.add_constant(design_matrix(df, terms), has_constant="add")
    return sm.OLS(df[y].to_numpy(float), X).fit()


def loo_r2(df: pd.DataFrame, y: str, terms=TERMS2) -> float:
    """Leave-one-out (predicted) R^2, computed from the hat matrix."""
    X = sm.add_constant(design_matrix(df, terms), has_constant="add").to_numpy(float)
    yv = df[y].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    res = yv - X @ beta
    H = X @ np.linalg.pinv(X.T @ X) @ X.T
    press = np.sum((res / (1.0 - np.diag(H))) ** 2)
    sst = np.sum((yv - yv.mean()) ** 2)
    return float(1.0 - press / sst) if sst > 0 else np.nan


def coefficient_table(res, y: str) -> pd.DataFrame:
    ci = np.asarray(res.conf_int())
    sd = res.model.exog.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    std_beta = res.params * sd / (res.model.endog.std(ddof=1) or 1.0)
    return pd.DataFrame({
        "response": y,
        "term": res.model.exog_names,
        "coefficient": res.params,
        "std_error": res.bse,
        "t_value": res.tvalues,
        "p_value": res.pvalues,
        "ci_low": ci[:, 0], "ci_high": ci[:, 1],
        "std_coefficient": std_beta,
    }).reset_index(drop=True)


def anova_table(df: pd.DataFrame, y: str, terms=TERMS2) -> pd.DataFrame:
    """Sequential (Type I) ANOVA grouped into linear / interaction / quadratic blocks."""
    yv = df[y].to_numpy(float)
    sst = float(np.sum((yv - yv.mean()) ** 2))
    blocks = [("Linear", ["x1", "x2", "x3"]),
              ("Two-factor interactions", ["x1x2", "x1x3", "x2x3"]),
              ("Quadratic (curvature)", ["x1^2", "x2^2", "x3^2"])]
    rows, used, ss_prev = [], [], 0.0
    for name, blk in blocks:
        blk = [t for t in blk if t in terms]
        if not blk:
            continue
        used += blk
        r = fit(df, y, used)
        ss = sst - float(np.sum(r.resid ** 2))
        rows.append({"source": name, "df": len(blk), "SS": ss - ss_prev})
        ss_prev = ss
    full = fit(df, y, terms)
    sse = float(np.sum(full.resid ** 2))
    dfe = int(full.df_resid)
    mse = sse / dfe if dfe else np.nan
    for r in rows:
        r["MS"] = r["SS"] / r["df"]
        r["F"] = r["MS"] / mse if mse else np.nan
    from scipy import stats
    for r in rows:
        r["p_value"] = float(1 - stats.f.cdf(r["F"], r["df"], dfe)) if np.isfinite(r["F"]) else np.nan
        r["contribution_pct"] = 100 * r["SS"] / sst if sst else np.nan
    rows.append({"source": "Residual", "df": dfe, "SS": sse, "MS": mse,
                 "F": np.nan, "p_value": np.nan, "contribution_pct": 100 * sse / sst if sst else np.nan})
    rows.append({"source": "Total", "df": len(df) - 1, "SS": sst, "MS": np.nan,
                 "F": np.nan, "p_value": np.nan, "contribution_pct": 100.0})
    out = pd.DataFrame(rows)
    out.insert(0, "response", y)
    return out


def summary_row(df: pd.DataFrame, y: str) -> dict:
    r2 = fit(df, y, TERMS2)
    r1 = fit(df, y, TERMS1)
    yv = df[y].to_numpy(float)
    n = len(df)
    rmse2 = float(np.sqrt(np.mean(r2.resid ** 2)))
    return {
        "response": y,
        "mean": float(yv.mean()), "std": float(yv.std(ddof=1)),
        "R2": float(r2.rsquared), "R2_adj": float(r2.rsquared_adj),
        "R2_pred_LOO": loo_r2(df, y, TERMS2),
        "RMSE": rmse2,
        "MAE": float(np.mean(np.abs(r2.resid))),
        "AIC": float(r2.aic), "BIC": float(r2.bic),
        "R2_first_order": float(r1.rsquared),
        "R2_pred_first_order": loo_r2(df, y, TERMS1),
        "RMSE_first_order": float(np.sqrt(np.mean(r1.resid ** 2))),
        "curvature_F": float(((np.sum(r1.resid ** 2) - np.sum(r2.resid ** 2)) / 3)
                             / (np.sum(r2.resid ** 2) / r2.df_resid)) if r2.df_resid else np.nan,
        "n_obs": n,
    }
