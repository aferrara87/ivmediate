"""
ivmediate 2.0.0: IV mediation analysis with a single instrument
===============================================================

Implementation of the IV mediation estimator of Ferrara, Gold, Heblich and
Pinto, "Identifying Causal Mediation Effects with a Single Instrument", with
diagnostics, sensitivity analysis, and support for classical/robust/clustered/
spatial standard errors and fixed effects.

Public interface (``import ivmediate``):
    ivmediate.ivmediate(data, Y, T, M, Z, ...)    estimate; returns a dict
    ivmediate.print_results(res, ...)             print the results
    ivmediate.bootstrap(data, Y, T, M, Z, ...)    bootstrap inference
    ivmediate.print_bootstrap(boot, ...)          print the bootstrap results

Authors: Andreas Ferrara (maintainer, a.ferrara@pitt.edu), Robert Gold,
Stephan Heblich, Rodrigo Pinto. License: MIT.

Citation: Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo Pinto
(2026), "Identifying Causal Mediation Effects with a Single Instrument,"
NBER Working Paper No. [number], doi:[DOI].

Requires: numpy, scipy, pandas; scikit-learn for spatial (Conley) SEs and
the LASSO; psutil (optional) to size the spatial kernel matrix.

Conventions
-----------
- Wald CIs use normal critical values z_{1−α/2}; AR sets use χ²_K.
- Degrees of freedom: absorbed FE count as the rank of the FE dummies
  (G − 1 one-way; G₁ + G₂ − c − 1 two-way, c connected components). With
  cluster SEs, FE dimensions nested in the clusters are not counted in the
  CR1 correction (reghdfe/fixest convention; with one of two dimensions
  nested, the other counts G − 1); HC1 and classical SEs count all
  absorbed FE.
- Olea–Pflueger critical values: simplified test at the estimated effective
  degrees of freedom K_eff (Montiel Olea and Pflueger 2013, eq. 4); K_eff = 1
  for one instrument (23.109 at τ = 10%, α = 5%).
- RV/XRV with K ≥ 2 instruments use the first instrument.
"""

import numpy as np
import pandas as pd
import re
import fnmatch
from scipy import stats
from scipy.linalg import solve_triangular
from scipy.optimize import brentq, minimize_scalar

__version__ = "2.0.0"

# ==============================================================================
#  Core Estimation Functions
# ==============================================================================

def lm_fit(X, Y):
    """Fast OLS via normal equations (like R's lm.fit)."""
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    resid = Y - X @ beta
    return dict(coefficients=beta, residuals=resid)


# ──────────────────────────────────────────────────────────────────────────────
#  Hajek-form weighted moments
#
#  Used for the κ-sensitivity (paper App. F.4, eq. 126). They use the
#  self-normalized (Hajek) form: Σ w_i x_i / Σ w_i, which targets the population
#  mean under both aweights (cell-mean interpretation) and pweights (inverse-
#  probability weighting, Σw ≈ population size). When w ≡ 1 these reduce to the
#  usual sample analogues (ddof = 0).
# ──────────────────────────────────────────────────────────────────────────────

def _weighted_mean(x, w):
    """Hajek-form weighted mean: Σ w_i x_i / Σ w_i."""
    return np.sum(w * x) / np.sum(w)


def _weighted_cov(x, y, w):
    """Hajek-form weighted covariance (ddof = 0):
        Σ w_i (x_i - x̄_w)(y_i - ȳ_w) / Σ w_i.
    Collapses to the MLE sample covariance when w ≡ 1.
    """
    wsum = np.sum(w)
    mx = np.sum(w * x) / wsum
    my = np.sum(w * y) / wsum
    return np.sum(w * (x - mx) * (y - my)) / wsum


def _weighted_var(x, w):
    """Hajek-form weighted variance."""
    return _weighted_cov(x, x, w)


def _weighted_std(x, w):
    """Hajek-form weighted standard deviation."""
    return np.sqrt(_weighted_var(x, w))


def _wls_fit(X, Y, w):
    """Weighted OLS on raw (un-premultiplied) data: β̂ = (X' diag(w) X)^{-1} X' diag(w) Y.

    Equivalent to OLS on √w-premultiplied data. Returns coefficients and raw-scale
    residuals (Y - X β̂), which are the usual WLS residuals — NOT premultiplied.
    """
    sw = np.sqrt(np.asarray(w, dtype=float))
    Q, R, scale = _equilibrated_qr(X * sw[:, None])
    beta = scale * solve_triangular(R, Q.T @ (Y * sw))
    resid = Y - X @ beta
    return dict(coefficients=beta, residuals=resid)


# ──────────────────────────────────────────────────────────────────────────────
#  Olea–Pflueger (2013) effective-F critical values
#
#  Source: Olea & Pflueger (2013), "A Robust Test for Weak Instruments",
#  Journal of Business & Economic Statistics 31(3), 358–369, Table 1.
#  Reject the null of weak instruments (worst-case 2SLS bias > τ × |OLS bias|)
#  if the effective F exceeds c(K, τ, α). Default τ = 0.10, α = 0.05 follow
#  OP's recommendation and are the standard values reported in applied work.
#
#  Tabulated for α = 0.05 (Wald-test size), K = 1 ... 30 (covers virtually
#  all applied IV mediation settings), τ ∈ {0.05, 0.10, 0.20, 0.30}.
#  For (τ, α) outside this grid, the lookup returns None and the printed
#  diagnostics omit the critical value with a note pointing to OP Table 1.
# ──────────────────────────────────────────────────────────────────────────────

_OP_CRIT_5PCT_WALD = {
    # K_inst -> {tau -> critical value}, alpha = 0.05.
    # OP's table column is x = 1/tau, so tau in {0.05,0.10,0.20,0.30}
    # maps to x in {20,10,5,3.33}. The canonical effective-F value
    # (K=1, 10% bias) is the x=10 cell = 23.11.
    # Source: Olea & Pflueger (2013, JBES) Table 1.
    1:  {0.05: 37.42, 0.10: 23.11, 0.20: 15.06, 0.30: 12.05},
    2:  {0.05: 32.32, 0.10: 19.29, 0.20: 12.17, 0.30:  9.57},
    3:  {0.05: 30.13, 0.10: 17.67, 0.20: 10.95, 0.30:  8.53},
    4:  {0.05: 28.85, 0.10: 16.72, 0.20: 10.23, 0.30:  7.92},
    5:  {0.05: 27.98, 0.10: 16.08, 0.20:  9.75, 0.30:  7.51},
    6:  {0.05: 27.35, 0.10: 15.62, 0.20:  9.40, 0.30:  7.21},
    7:  {0.05: 26.86, 0.10: 15.26, 0.20:  9.14, 0.30:  6.98},
    8:  {0.05: 26.47, 0.10: 14.97, 0.20:  8.92, 0.30:  6.80},
    9:  {0.05: 26.15, 0.10: 14.73, 0.20:  8.74, 0.30:  6.65},
    10: {0.05: 25.87, 0.10: 14.53, 0.20:  8.59, 0.30:  6.52},
    11: {0.05: 25.64, 0.10: 14.36, 0.20:  8.47, 0.30:  6.41},
    12: {0.05: 25.44, 0.10: 14.21, 0.20:  8.36, 0.30:  6.32},
    13: {0.05: 25.26, 0.10: 14.08, 0.20:  8.26, 0.30:  6.24},
    14: {0.05: 25.10, 0.10: 13.96, 0.20:  8.17, 0.30:  6.16},
    15: {0.05: 24.96, 0.10: 13.86, 0.20:  8.10, 0.30:  6.10},
    16: {0.05: 24.83, 0.10: 13.77, 0.20:  8.03, 0.30:  6.04},
    17: {0.05: 24.71, 0.10: 13.68, 0.20:  7.96, 0.30:  5.99},
    18: {0.05: 24.60, 0.10: 13.60, 0.20:  7.91, 0.30:  5.94},
    19: {0.05: 24.50, 0.10: 13.53, 0.20:  7.85, 0.30:  5.89},
    20: {0.05: 24.41, 0.10: 13.46, 0.20:  7.80, 0.30:  5.85},
    21: {0.05: 24.33, 0.10: 13.40, 0.20:  7.76, 0.30:  5.81},
    22: {0.05: 24.25, 0.10: 13.35, 0.20:  7.72, 0.30:  5.78},
    23: {0.05: 24.18, 0.10: 13.29, 0.20:  7.68, 0.30:  5.74},
    24: {0.05: 24.11, 0.10: 13.24, 0.20:  7.64, 0.30:  5.71},
    25: {0.05: 24.05, 0.10: 13.20, 0.20:  7.61, 0.30:  5.68},
    26: {0.05: 23.98, 0.10: 13.15, 0.20:  7.57, 0.30:  5.66},
    27: {0.05: 23.93, 0.10: 13.11, 0.20:  7.54, 0.30:  5.63},
    28: {0.05: 23.87, 0.10: 13.07, 0.20:  7.51, 0.30:  5.61},
    29: {0.05: 23.82, 0.10: 13.04, 0.20:  7.49, 0.30:  5.58},
    30: {0.05: 23.77, 0.10: 13.00, 0.20:  7.46, 0.30:  5.56},
}


# Guard: the K=1 / 10%-bias cell must equal OP's canonical effective-F value.
# (x = 1/tau, so this is the x=10 column.) Protects against a future re-keying.
assert abs(_OP_CRIT_5PCT_WALD[1][0.10] - 23.11) < 1e-9, (
    f"OP table mis-keyed: K=1, tau=0.10 should be 23.11, got "
    f"{_OP_CRIT_5PCT_WALD[1][0.10]}. Check the x = 1/tau column mapping."
)


def _op_critical_value(K, tau, alpha):
    """Olea–Pflueger (2013) simplified-test critical value at K_eff = K.

    c = χ²_{K}(x·K)⁻¹(1 − α) / K with x = 1/τ (noncentral χ² quantile); it
    reproduces Table 1 (`_OP_CRIT_5PCT_WALD`) to its two decimals at α = 0.05.
    Since 2026-10 the analysis evaluates this formula at the estimated K_eff
    (see _op_keff_critval) for any τ, α in (0, 1); this helper keeps the
    K_eff = K special case.
    """
    _check_open_unit("tau", tau)
    _check_open_unit("alpha", alpha)
    K = float(K)
    x = 1.0 / float(tau)
    return float(stats.ncx2.ppf(1.0 - float(alpha), K, x * K) / K)


def _op_effective_F(beta_Z, V_ZZ, Q_zz):
    """Olea–Pflueger (2013) effective F-statistic.

        F_eff = (β̂_Z' Q_zz β̂_Z) / tr(Q_zz · V̂_ZZ)

    where β̂_Z is the first-stage coefficient vector on the K excluded
    instruments, V̂_ZZ is its robust/cluster vcov (matching the SE regime
    used for the structural inference), and Q_zz = Z̃' Z̃ / n with Z̃ the
    instruments residualized against included exogenous.

    For K = 1 this collapses to β̂_Z² / V̂_ZZ — i.e., the robust/cluster
    Wald t² on the excluded instrument — independent of the scaling of Q_zz.
    """
    beta_Z = np.atleast_1d(np.asarray(beta_Z, dtype=float))
    V_ZZ  = np.atleast_2d(np.asarray(V_ZZ, dtype=float))
    Q_zz  = np.atleast_2d(np.asarray(Q_zz, dtype=float))
    num = float(beta_Z @ Q_zz @ beta_Z)
    den = float(np.trace(Q_zz @ V_ZZ))
    return num / den


def _is_pd(V):
    """True if V is finite and positive definite (a usable variance block).
    The spatial (Conley) sandwich can be indefinite; Wald and OP statistics
    built on such a block are reported as NaN."""
    V = np.atleast_2d(np.asarray(V, dtype=float))
    if V.size == 0 or not np.all(np.isfinite(V)):
        return False
    return bool(np.min(np.linalg.eigvalsh((V + V.T) / 2.0)) > 0)


# ══════════════════════════════════════════════════════════════════════════════
#  Fixed-Effects Handling: Singleton Dropping and Within-Group Demeaning
# ══════════════════════════════════════════════════════════════════════════════

def _drop_singletons(df, fe_col):
    """Iteratively drop singleton fixed-effect groups until none remain.

    A singleton is a group with only one observation. Dropping singletons
    can create new singletons, so we iterate. This follows the reghdfe
    convention (Correia, 2016).

    Returns: (cleaned DataFrame, number of groups, number of obs dropped)
    """
    n_before = len(df)
    while True:
        counts = df[fe_col].value_counts()
        singletons = counts[counts == 1].index
        if len(singletons) == 0:
            break
        df = df[~df[fe_col].isin(singletons)].copy()
    n_after = len(df)
    n_groups = df[fe_col].nunique()
    return df, n_groups, n_before - n_after


def _group_means(codes, w, sw, V):
    """Weighted group means of the columns of V, expanded back to rows."""
    G = sw.shape[0]
    out = np.empty_like(V)
    for j in range(V.shape[1]):
        out[:, j] = (np.bincount(codes, weights=w * V[:, j], minlength=G) / sw)[codes]
    return out


def _demean_within_groups(df, cols, fe_col, weights=None):
    """Demean the columns `cols` within groups of `fe_col`.

    With `weights`, the group means are weighted: x − Σ_g w·x / Σ_g w (as in
    reghdfe/fixest with aweights or pweights); without, ordinary means. All
    columns are demeaned on the same rows (FWL consistency).

    Returns: DataFrame with demeaned columns (original column names preserved).
    """
    df_out = df.copy()
    codes = pd.factorize(df[fe_col], sort=False)[0]
    w = np.ones(len(df)) if weights is None else np.asarray(weights, dtype=float)
    sw = np.bincount(codes, weights=w)
    V = df[cols].to_numpy(dtype=float)
    df_out[cols] = V - _group_means(codes, w, sw, V)
    return df_out


def _demean_twoway_direct(V, c1, c2, w, max_small=3000):
    """Exact weighted two-way FE residualization (Addendum A.2, O-2; the same
    algorithm as Stata's _ivm_demean2_direct).

    Model V = α_A + β_B + R by weighted least squares, with A the FE dimension
    with more levels and B the other one. Eliminating α gives the normal
    equations L β = r for the B effects, with the Schur complement
        L = diag(s_B) − Nᵀ diag(1/s_A) N,
    where s_A and s_B are the weight sums per level, N[a, b] is the weight of
    the observations in level a of A and level b of B, and r holds the B-level
    sums of w·(V − within-A means). L is singular, with one null direction per
    connected component of the FE graph, so one B level per component is fixed
    at 0. One Cholesky factorization serves all columns. The residual is
        R = (V − mean_A V) − (β_B − mean_A β_B).
    Returns None when the smaller dimension has more than `max_small` levels.
    """
    from scipy.sparse import csr_matrix, diags
    from scipy.sparse.csgraph import connected_components
    from scipy.linalg import cho_factor, cho_solve
    G1, G2 = int(c1.max()) + 1, int(c2.max()) + 1
    if G2 > G1:
        cA, cB, Ga, Gb = c2, c1, G2, G1
    else:
        cA, cB, Ga, Gb = c1, c2, G1, G2
    if Gb > max_small:
        return None
    sA = np.bincount(cA, weights=w, minlength=Ga)
    sB = np.bincount(cB, weights=w, minlength=Gb)
    N = csr_matrix((w, (cA, cB)), shape=(Ga, Gb))          # duplicates are summed
    Ns = diags(1.0 / np.sqrt(sA)) @ N
    L = np.diag(sB) - (Ns.T @ Ns).toarray()
    L = 0.5 * (L + L.T)
    # One B level per connected component of the bipartite FE graph is fixed at 0
    # (the first level of each component, in the order of the B codes).
    adj = csr_matrix((np.ones(len(cA)), (cA, Ga + cB)), shape=(Ga + Gb, Ga + Gb))
    lab = connected_components(adj, directed=False)[1]
    keep = np.ones(Gb, dtype=bool)
    keep[np.unique(lab[Ga:], return_index=True)[1]] = False
    V1 = V - _group_means(cA, w, sA, V)
    p = V.shape[1]
    rhs = np.empty((Gb, p))
    for j in range(p):
        rhs[:, j] = np.bincount(cB, weights=w * V1[:, j], minlength=Gb)
    beta = np.zeros((Gb, p))
    kb = np.flatnonzero(keep)
    if kb.size:
        Lk = L[np.ix_(kb, kb)]
        try:
            beta[kb] = cho_solve(cho_factor(Lk, check_finite=False), rhs[kb], check_finite=False)
        except np.linalg.LinAlgError:
            beta[kb] = np.linalg.lstsq(Lk, rhs[kb], rcond=None)[0]
    Bv = beta[cB]
    return V1 - (Bv - _group_means(cA, w, sA, Bv))


def _demean_twoway(df, cols, fe_cols, weights=None, tol=1e-12, max_iter=10000,
                   max_direct=3000):
    """Two-way FE demeaning (weighted).

    When the FE dimension with fewer levels has at most `max_direct` (3,000)
    levels, the exact projection of _demean_twoway_direct is used (Addendum
    A.2, O-2: alternating projections can stall on weakly connected designs,
    e.g., a worker-firm chain). Otherwise: alternating (weighted) projections
    with Irons–Tuck acceleration. One iteration demeans within fe_cols[0] and
    then within fe_cols[1]. Convergence (2026-10): the largest change in one
    iteration, divided by the standard deviation of the original column, is
    below `tol` (1e-12) for every column that varies; up to `max_iter`
    iterations. The criterion is relative, so variables measured in large
    units converge like any other. A warning is printed only on true
    non-convergence.

    Returns: DataFrame with demeaned columns.
    """
    df_out = df.copy()
    c1 = pd.factorize(df[fe_cols[0]], sort=False)[0]
    c2 = pd.factorize(df[fe_cols[1]], sort=False)[0]
    w = np.ones(len(df)) if weights is None else np.asarray(weights, dtype=float)
    V = df[cols].to_numpy(dtype=float)
    if V.shape[1] == 0:
        return df_out
    R = _demean_twoway_direct(V, c1, c2, w, max_small=max_direct)
    if R is not None:
        df_out[cols] = R
        return df_out
    sw1 = np.bincount(c1, weights=w)
    sw2 = np.bincount(c2, weights=w)
    sd = V.std(axis=0)
    active = sd > 1e-14 * np.maximum(np.abs(V).max(axis=0), 1e-300)
    scale = np.maximum(sd, 1e-300)

    def sweep(A):
        A = A - _group_means(c1, w, sw1, A)
        return A - _group_means(c2, w, sw2, A)

    converged = not active.any()
    change = 0.0
    it = 0
    while not converged and it < max_iter:
        F = sweep(V)
        it += 1
        change = float(np.max(np.abs(F - V)[:, active].max(axis=0) / scale[active]))
        if change < tol:
            V = F
            converged = True
            break
        Gs = sweep(F)
        it += 1
        D1 = F - V
        D2 = Gs - F
        DD = D2 - D1
        den = np.sum(DD * DD, axis=0)
        coef = np.where(den > 0, np.sum(D2 * DD, axis=0) / np.where(den > 0, den, 1.0), 0.0)
        V = Gs - coef[None, :] * D2          # Irons–Tuck extrapolation
    if not converged:
        print(f"  Warning: two-way FE demeaning did not converge after "
              f"{max_iter} iterations (last relative change = {change:.2e} > "
              f"tol = {tol:.0e}). Results may be imprecise; check for "
              f"nearly-collinear FE dimensions.")
    df_out[cols] = V
    return df_out


def _expand_dummies(df, dummy_vars):
    """Expand categorical variables into dummy columns (drop-first encoding).

    For each variable in dummy_vars, creates indicator columns for all but the
    first (alphabetically sorted) category. Appends columns to df and returns
    the list of new column names.

    Parameters
    ----------
    df : DataFrame (modified in place)
    dummy_vars : list of column names to expand

    Returns: list of new dummy column names
    """
    new_cols = []
    for var in dummy_vars:
        cats = sorted(df[var].dropna().unique())
        if len(cats) <= 1:
            continue
        # Drop first category (reference)
        for cat in cats[1:]:
            col_name = f"_D_{var}_{cat}"
            df[col_name] = (df[var] == cat).astype(float)
            new_cols.append(col_name)
    return new_cols


# ══════════════════════════════════════════════════════════════════════════════
#  Variance-Covariance Estimators for 2SLS
# ══════════════════════════════════════════════════════════════════════════════

#  Numerics (2026-10): every regression is solved by a QR decomposition of the
#  equilibrated design (columns scaled to unit norm), never by inverting X'X.
#  With X̂·D = Q·R (D = diag of column scalings), the 2SLS bread is
#  (X̂'X̂)^{-1} = T·T' with T = D·R^{-1}, and the score of observation i in the
#  rotated coordinates is Q_i·ê_i, so every sandwich is V = T·meat_Q·T' with
#  meat_Q built from the rows of Q. This keeps ill-conditioned but full-rank
#  designs (calendar-year polynomials, variables in large units) accurate.

def _equilibrated_qr(A):
    """QR of A with columns scaled to unit Euclidean norm.

    Returns (Q, R, scale) with A·diag(scale) = Q·R. Zero columns keep scale 1.
    Raises numpy.linalg.LinAlgError if the scaled matrix is numerically
    singular (|R_jj| < 1e-12·max|R_ii|), mirroring np.linalg.solve.
    """
    A = np.asarray(A, dtype=float)
    norms = np.sqrt(np.sum(A * A, axis=0))
    scale = np.where(norms > 0, 1.0 / np.where(norms > 0, norms, 1.0), 1.0)
    Q, R = np.linalg.qr(A * scale[None, :], mode="reduced")
    d = np.abs(np.diag(R))
    if d.size and (not np.all(np.isfinite(d)) or d.min() <= 1e-12 * max(d.max(), 1e-300)):
        raise np.linalg.LinAlgError("Singular matrix")
    return Q, R, scale


def _vcov_classical(resid, T, n, k, df_adjust=0):
    """Classical (homoskedastic) 2SLS variance σ̂²(X̂'X̂)^{-1} = σ̂²·T·T'."""
    sigma2 = np.sum(resid ** 2) / (n - k - df_adjust)
    return sigma2 * (T @ T.T)


def _vcov_robust(resid, Q, T, n, k, df_adjust=0):
    """HC1 sandwich (X̂'X̂)^{-1}(Σ ê_i² x̂_i x̂_i')(X̂'X̂)^{-1}·n/(n−k−df_adjust),
    computed as T·[(Q·ê)'(Q·ê)]·T' (see the numerics note above)."""
    S = Q * resid[:, None]
    meat = (S.T @ S) * (n / (n - k - df_adjust))
    return T @ meat @ T.T


def _vcov_cluster(resid, Q, T, n, k, cluster_ids, df_adjust=0,
                  cluster_unique=None, cluster_inverse=None):
    """Clustered (Liang–Zeger CR1) sandwich with factor
    G/(G−1)·(n−1)/(n−k−df_adjust), where S_g = Σ_{i∈g} ê_i x̂_i.

    `df_adjust` is the absorbed-FE degrees of freedom that count toward CR1:
    FE dimensions nested in the clusters are excluded by the caller
    (reghdfe/fixest convention). Cluster sums use np.bincount; precomputed
    np.unique results can be passed to avoid recomputing them.
    """
    if cluster_unique is None or cluster_inverse is None:
        cluster_unique, cluster_inverse = np.unique(cluster_ids, return_inverse=True)
    G = len(cluster_unique)
    if G < 2:
        raise ValueError("cluster SEs need at least 2 clusters.")
    S = Q * resid[:, None]
    kq = S.shape[1]
    cluster_sums = np.empty((G, kq))
    for _j in range(kq):
        cluster_sums[:, _j] = np.bincount(cluster_inverse, weights=S[:, _j], minlength=G)
    correction = (G / (G - 1)) * ((n - 1) / (n - k - df_adjust))
    meat = (cluster_sums.T @ cluster_sums) * correction
    return T @ meat @ T.T


# ──────────────────────────────────────────────────────────────────────────────
#  Spatial autocorrelation–robust vcov (Conley 1999)
#
#  For each pair of observations (i, j), the spatial sandwich estimator
#  weights the score outer-product s_i s_j' by a kernel function K(d_{ij})
#  of the distance between the two observations:
#
#      Ω̂ = Σ_i Σ_j K(d_{ij}) · ŝ_i ŝ_j'
#
#  with score ŝ_i = ê_i · x̂_i and kernel either uniform (hard cutoff at
#  `cutoff` distance) or Bartlett (linear taper to zero at `cutoff`).
#
#  Implementation uses a BallTree neighbor index (Haversine metric for
#  lat/lon in degrees, Euclidean for projected coordinates) to enumerate
#  only the within-cutoff pairs, giving O(n × N̄) complexity rather than
#  O(n²), where N̄ is the average number of within-cutoff neighbors.
#
#  Cluster-robust SEs are the special case where K(d_{ij}) = 1{cluster(i)
#  = cluster(j)}; spatial generalizes to a continuous kernel in physical
#  space rather than discrete partition into groups.
# ──────────────────────────────────────────────────────────────────────────────

# Earth radius in km (mean), used to convert Haversine angular distances
# (returned by sklearn's BallTree in radians) to physical distances in km.
_EARTH_RADIUS_KM = 6371.0


def _build_spatial_neighbors(lat, lon, cutoff_km, metric="haversine"):
    """Build a sparse list of within-cutoff neighbor pairs.

    Returns a list of (idx_i, idx_j_array, weight_j_array) tuples (one per
    observation). idx_j_array contains all j ≠ i within the cutoff (including
    i ≤ j and i > j; we keep both directions because the spatial sandwich
    sum is symmetric over the full pair set, including i = j on the diagonal).
    """
    from sklearn.neighbors import BallTree
    n = len(lat)
    coords_deg = np.column_stack([lat, lon])
    if metric == "haversine":
        coords = np.radians(coords_deg)
        cutoff_native = cutoff_km / _EARTH_RADIUS_KM
        tree = BallTree(coords, metric="haversine")
    elif metric == "euclidean":
        coords = coords_deg
        cutoff_native = cutoff_km
        tree = BallTree(coords, metric="euclidean")
    else:
        raise ValueError(f"Unknown spatial metric {metric!r}; "
                         f"use 'haversine' or 'euclidean'.")
    # query_radius returns indices and distances of neighbors within cutoff
    indices, distances = tree.query_radius(coords, r=cutoff_native,
                                           return_distance=True)
    # Convert distances back to km (only matters for Bartlett kernel)
    if metric == "haversine":
        distances = [d * _EARTH_RADIUS_KM for d in distances]
    return indices, distances


# Memory threshold (in bytes) above which the sparse adjacency build falls
# back to the loop-based per-row computation: 25% of available RAM (psutil),
# capped at 4 GB to avoid pathological cases on shared hardware; 4 GB without
# psutil.
def _spatial_sparse_memory_budget():
    """Return the memory budget (in bytes) for the sparse kernel adjacency."""
    try:
        import psutil
        avail = psutil.virtual_memory().available
        budget = int(0.25 * avail)
    except ImportError:
        # psutil not installed: use the 4 GB cap.
        budget = 4 * (1024 ** 3)   # 4 GB
    # Cap at 4 GB regardless of available RAM
    return min(budget, 4 * (1024 ** 3))


def _build_spatial_sparse_adjacency(indices, distances, kernel, cutoff_km,
                                    n, memory_budget=None):
    """Build a kernel-weighted sparse adjacency matrix from neighbor lists.

    Returns
    -------
    A : scipy.sparse.csr_matrix of shape (n, n), or None
        Sparse matrix where A[i, j] = K(d_{ij}) for j ∈ N(i), zero otherwise.
        Returns None if the predicted memory cost exceeds memory_budget; in
        that case the caller must fall back to the loop-based computation.

    The kernel-weighted sandwich Ω̂ = Σ_i Σ_{j ∈ N(i)} K(d_{ij}) ŝ_i ŝ_j'
    can then be computed in a single matrix multiplication:
        Ω̂ = scores.T @ A @ scores
    which is O(nnz(A) · k) — typically 10-100x faster than the per-row loop
    for moderate-to-large n.
    """
    from scipy.sparse import csr_matrix

    # Predict total nonzero count; this also predicts memory cost
    nnz = int(sum(len(idx) for idx in indices))
    # Memory cost: CSR stores int32 indices, int32 indptr (n+1 entries),
    # and float64 data. Roughly 16 bytes per nonzero plus 4*(n+1) for indptr.
    predicted_bytes = 16 * nnz + 4 * (n + 1)

    if memory_budget is None:
        memory_budget = _spatial_sparse_memory_budget()
    if predicted_bytes > memory_budget:
        # Caller falls back to loop-based path.
        return None

    # Build CSR directly: avoids constructing a COO matrix and converting
    # (which would temporarily double memory).
    indptr = np.empty(n + 1, dtype=np.int32)
    indptr[0] = 0
    indptr[1:] = np.cumsum([len(idx) for idx in indices])
    if kernel == "uniform":
        # Concatenate-and-fill avoids per-row Python overhead
        all_indices = np.concatenate(indices).astype(np.int32, copy=False)
        all_data = np.ones(nnz, dtype=np.float64)
    elif kernel == "bartlett":
        all_indices = np.concatenate(indices).astype(np.int32, copy=False)
        all_distances = np.concatenate(distances)
        all_data = np.maximum(1.0 - all_distances / cutoff_km, 0.0)
    else:
        raise ValueError(f"Unknown spatial kernel {kernel!r}; "
                         f"use 'uniform' or 'bartlett'.")
    return csr_matrix((all_data, all_indices, indptr), shape=(n, n))


def _vcov_spatial(resid, Q, T, n, k, lat, lon, cutoff_km, kernel,
                  metric, df_adjust=0, spatial_indices=None,
                  spatial_distances=None, spatial_adjacency=None,
                  neg_counter=None):
    """Spatial 2SLS variance (Conley 1999).

    Sandwich: (X̂'X̂)^{-1} Ω̂ (X̂'X̂)^{-1}, where
        Ω̂ = Σ_i Σ_{j in N(i)} K(d_{ij}) ŝ_i ŝ_j'
    and ŝ_i = ê_i x̂_i, N(i) is the set of observations within `cutoff_km`
    of observation i, and K is `kernel` (uniform or Bartlett). Computed in
    the rotated coordinates of the QR solve (see the numerics note above).

    Small-sample correction: (n-1)/(n-k-df_adjust). Conley's original paper
    provides no closed-form correction analogous to cluster-robust's
    G/(G-1); applied practice typically uses just (n-1)/(n-k) as a mild
    finite-sample adjustment.

    Implementation: for moderate-to-large n, builds a sparse kernel-weighted
    adjacency matrix and computes Ω̂ as a single matrix multiplication
    (scores.T @ A @ scores). Falls back to a per-row loop if the sparse
    matrix would exceed the memory budget.

    Optional cache arguments (spatial_indices, spatial_distances,
    spatial_adjacency) avoid rebuilding the BallTree and sparse adjacency
    when the same coordinates and kernel are reused across multiple vcov
    calls (typically 7 within a single ivmediate call).
    `neg_counter` (a dict) suppresses the per-call warning and counts the
    negative-variance events (diagnostic only); ivmediate warns
    once per analysis, and only for the quantities it reports.
    """
    scores = Q * resid[:, None]       # (n, k) rotated scores
    k_cols = scores.shape[1]

    # Reuse cached neighbors if provided; otherwise build.
    if spatial_indices is None or spatial_distances is None:
        indices, distances = _build_spatial_neighbors(lat, lon, cutoff_km, metric)
    else:
        indices, distances = spatial_indices, spatial_distances

    # Reuse cached adjacency if provided; otherwise build.
    if spatial_adjacency is not None:
        A = spatial_adjacency
    else:
        A = _build_spatial_sparse_adjacency(indices, distances, kernel,
                                             cutoff_km, n)

    if A is not None:
        # Vectorized path: meat = scores' A scores
        meat = scores.T @ (A @ scores)
    else:
        # Memory-pressure fallback: per-row loop (original implementation).
        # This preserves correctness on very large datasets where the sparse
        # adjacency would not fit in RAM.
        print(f"  Note: spatial vcov sparse adjacency would exceed memory "
              f"budget; using per-row loop. This is slower but uses O(n) "
              f"memory rather than O(n × avg_neighbors).")
        meat = np.zeros((k_cols, k_cols))
        for i in range(n):
            nbr_idx = indices[i]
            nbr_d = distances[i]
            if len(nbr_idx) == 0:
                continue
            if kernel == "uniform":
                w = np.ones_like(nbr_d)
            elif kernel == "bartlett":
                w = np.maximum(1.0 - nbr_d / cutoff_km, 0.0)
            else:
                raise ValueError(f"Unknown spatial kernel {kernel!r}; "
                                 f"use 'uniform' or 'bartlett'.")
            nbr_scores_weighted = (w[:, None] * scores[nbr_idx]).sum(axis=0)
            meat += np.outer(scores[i], nbr_scores_weighted)
    # Small-sample correction
    correction = (n - 1) / max(n - k - df_adjust, 1)
    meat *= correction
    vcov = T @ meat @ T.T
    # The spatial (kernel-HAC) sandwich is not guaranteed to be positive
    # semi-definite in finite samples: the uniform kernel matrix is not PSD in
    # two dimensions, whatever the cutoff. Negative variances become NaN SEs;
    # the caller prints one warning per analysis (or this function does, when
    # called without a counter).
    n_bad = int((np.diag(vcov) < 0).sum())
    if n_bad:
        if neg_counter is not None:
            neg_counter["neg_var"] = neg_counter.get("neg_var", 0) + n_bad
        else:
            print(f"  Warning: spatial vcov has {n_bad} negative diagonal entries "
                  f"(not positive semi-definite); the affected SEs are NaN.")
    return vcov


def _cross_equation_cov(res1, res2, idx1, idx2,
                        se="robust", cluster_ids=None, df_adjust=0,
                        lat=None, lon=None, spatial_cutoff=None,
                        spatial_kernel="uniform", spatial_metric="haversine",
                        vcov_cache=None):
    """Cross-equation covariance Cov(β̂_{1,idx1}, β̂_{2,idx2}) from a stacked-GMM
    treatment of two 2SLS regressions sharing observations.

    `res1` and `res2` are the dicts returned by `_2sls_with_vcov` for the two
    equations (they carry the residuals and the QR factors Q, T with
    (X̂'X̂)^{-1} = T·T' and X̂·T = Q). idx1 and idx2 select the coefficient of
    interest in each equation (e.g., γ̂ in the mediator equation, θ̂ in the
    outcome equation).

    The influence function of β̂_{j,ℓ} for observation i is
        ψ_{j,ℓ,i} = [(X̂_j' X̂_j)^{-1}]_{ℓ,:} · x̂_{j,i} · resid_{j,i}
                  = [T_j Q_{j,i}']_ℓ · resid_{j,i}

    Cross-equation covariance under each SE regime mirrors the diagonal cases:
      * classical:  σ̂_12 · [(X̂_1'X̂_1)^{-1} X̂_1'X̂_2 (X̂_2'X̂_2)^{-1}]_{idx1,idx2}
                    = σ̂_12 · [T_1 Q_1'Q_2 T_2']_{idx1,idx2}
      * robust:     Σ_i ψ_{1,idx1,i} · ψ_{2,idx2,i}     (HC1 with same n/(n-k) scale)
      * cluster:    Σ_g (Σ_i∈g ψ_{1,idx1,i})(Σ_i∈g ψ_{2,idx2,i})
      * spatial:    Σ_i Σ_j∈N(i) K(d_ij) · ψ_{1,idx1,i} · ψ_{2,idx2,j}

    The same small-sample correction factor used by the diagonal sandwich is
    applied here, with k = max(k1, k2) to be conservative.

    Used to compute the analytic SE on the indirect effect θγ via the delta
    method:
        Var(θγ) ≈ γ² Var(θ) + θ² Var(γ) + 2θγ Cov(θ, γ).
    """
    if se == "none":
        return float("nan")
    resid1, Q1, T1 = res1["resid"], res1["Q"], res1["T"]
    resid2, Q2, T2 = res2["resid"], res2["Q"], res2["T"]
    n = len(resid1)
    k_eff = max(T1.shape[0], T2.shape[0])

    # Influence-function vectors for the two coefficients of interest.
    psi1 = (Q1 @ T1[idx1, :]) * resid1            # (n,)
    psi2 = (Q2 @ T2[idx2, :]) * resid2            # (n,)

    if se == "classical":
        sigma_12 = float(resid1 @ resid2) / max(n - k_eff - df_adjust, 1)
        cross_block = T1 @ (Q1.T @ Q2) @ T2.T
        cov = sigma_12 * cross_block[idx1, idx2]

    elif se == "robust":
        meat = float(np.sum(psi1 * psi2))
        cov = meat * (n / (n - k_eff - df_adjust))

    elif se == "cluster":
        if cluster_ids is None:
            raise ValueError("cluster_ids required when se='cluster'")
        if vcov_cache and 'cluster_unique' in vcov_cache:
            cu, ci = vcov_cache['cluster_unique'], vcov_cache['cluster_inverse']
        else:
            cu, ci = np.unique(cluster_ids, return_inverse=True)
        G = len(cu)
        s1 = np.bincount(ci, weights=psi1, minlength=G)
        s2 = np.bincount(ci, weights=psi2, minlength=G)
        meat = float(np.sum(s1 * s2))
        correction = (G / (G - 1)) * ((n - 1) / (n - k_eff - df_adjust))
        cov = meat * correction

    elif se == "spatial":
        if lat is None or lon is None or spatial_cutoff is None:
            raise ValueError("se='spatial' requires lat, lon, spatial_cutoff.")
        if vcov_cache and 'spatial_indices' in vcov_cache:
            indices = vcov_cache['spatial_indices']
            distances = vcov_cache['spatial_distances']
        else:
            indices, distances = _build_spatial_neighbors(
                lat, lon, spatial_cutoff, spatial_metric)
        if vcov_cache and 'spatial_adjacency' in vcov_cache:
            A = vcov_cache['spatial_adjacency']
        else:
            A = _build_spatial_sparse_adjacency(
                indices, distances, spatial_kernel, spatial_cutoff, n)
        if A is not None:
            meat = float(psi1 @ (A @ psi2))
        else:
            print("  Note: spatial cross-equation covariance using per-row "
                  "loop (memory-budget fallback).")
            meat = 0.0
            for i in range(n):
                nbr_idx = indices[i]
                nbr_d = distances[i]
                if len(nbr_idx) == 0:
                    continue
                if spatial_kernel == "uniform":
                    w = np.ones_like(nbr_d)
                elif spatial_kernel == "bartlett":
                    w = np.maximum(1.0 - nbr_d / spatial_cutoff, 0.0)
                else:
                    raise ValueError(f"Unknown spatial kernel {spatial_kernel!r}")
                meat += float(psi1[i]) * float(np.dot(w, psi2[nbr_idx]))
        correction = (n - 1) / max(n - k_eff - df_adjust, 1)
        cov = meat * correction

    else:
        raise ValueError(f"Unknown se type: {se}")

    return float(cov)


def _2sls_with_vcov(Y, X, Z, se="robust", cluster_ids=None, df_adjust=0,
                    lat=None, lon=None, spatial_cutoff=None,
                    spatial_kernel="uniform", spatial_metric="haversine",
                    vcov_cache=None):
    """2SLS (or OLS when `Z is X`) with a flexible variance estimator.

    Point estimates come from QR decompositions of the equilibrated designs:
    X̂ = P_Z X is formed from a QR of Z, and β̂ from a QR of X̂ (see the
    numerics note above `_equilibrated_qr`).

    Parameters
    ----------
    se : {"classical", "robust", "cluster", "spatial", "none"}
        Variance estimator. "spatial" requires lat, lon, spatial_cutoff.
        "none" (internal) skips the variance and returns NaN SEs; used by
        bootstrap replicates that need point estimates only.
    cluster_ids : array or None
        Required when se="cluster".
    df_adjust : int
        Additional degrees of freedom absorbed (e.g., from fixed effects).
        Subtracted from n - k in all variance calculations.
    lat, lon : arrays or None
        Coordinates for spatial vcov (required when se="spatial").
    spatial_cutoff : float or None
        Distance cutoff for spatial vcov. Units are kilometers when
        spatial_metric="haversine", otherwise the units of lat/lon.
    spatial_kernel : {"uniform", "bartlett"}
        Kernel function for spatial weighting. Default "uniform".
    spatial_metric : {"haversine", "euclidean"}
        Distance metric. Default "haversine" (lat/lon in degrees).
    vcov_cache : dict or None
        Optional cache of precomputed quantities that don't change across
        the multiple _2sls_with_vcov calls within a single
        ivmediate invocation. Recognized keys:
          - 'cluster_unique', 'cluster_inverse': result of np.unique(
            cluster_ids, return_inverse=True). Used by cluster vcov.
          - 'spatial_indices', 'spatial_distances': result of
            _build_spatial_neighbors. Used by spatial vcov.
          - 'spatial_adjacency': sparse CSR adjacency matrix from
            _build_spatial_sparse_adjacency. Used by spatial vcov.
          - 'neg_counter': dict collecting negative spatial variances.
        Caching avoids ~7x repeated computation per analysis call.

    Returns: dict with coef, se, vcov, resid, X_hat, sigma2, df, and the QR
    factors Q (rotated design, X̂·T = Q) and T ((X̂'X̂)^{-1} = T·T').
    """
    n = len(Y)
    k = X.shape[1]
    effective_df = n - k - df_adjust

    # ── Point estimates (QR of equilibrated designs) ──
    if Z is X:
        X_hat = X
    else:
        Qz, _, _ = _equilibrated_qr(Z)
        X_hat = Qz @ (Qz.T @ X)
    Q, R, scale = _equilibrated_qr(X_hat)
    beta = scale * solve_triangular(R, Q.T @ Y)
    resid = Y - X @ beta
    sigma2 = np.sum(resid ** 2) / effective_df
    T = scale[:, None] * solve_triangular(R, np.eye(k))

    # ── Variance-covariance ──
    if se == "none":
        vcov = np.full((k, k), np.nan)
    elif se == "classical":
        vcov = _vcov_classical(resid, T, n, k, df_adjust)
    elif se == "robust":
        vcov = _vcov_robust(resid, Q, T, n, k, df_adjust)
    elif se == "cluster":
        if cluster_ids is None:
            raise ValueError("cluster_ids required when se='cluster'")
        cu = vcov_cache.get('cluster_unique') if vcov_cache else None
        ci = vcov_cache.get('cluster_inverse') if vcov_cache else None
        vcov = _vcov_cluster(resid, Q, T, n, k, cluster_ids, df_adjust,
                             cluster_unique=cu, cluster_inverse=ci)
    elif se == "spatial":
        if lat is None or lon is None or spatial_cutoff is None:
            raise ValueError("se='spatial' requires lat, lon, and spatial_cutoff.")
        si = vcov_cache.get('spatial_indices') if vcov_cache else None
        sd = vcov_cache.get('spatial_distances') if vcov_cache else None
        sa = vcov_cache.get('spatial_adjacency') if vcov_cache else None
        nc = vcov_cache.get('neg_counter') if vcov_cache else None
        vcov = _vcov_spatial(resid, Q, T, n, k, lat, lon,
                             spatial_cutoff, spatial_kernel, spatial_metric,
                             df_adjust, spatial_indices=si,
                             spatial_distances=sd, spatial_adjacency=sa,
                             neg_counter=nc)
    else:
        raise ValueError(f"Unknown se type: {se}. "
                         f"Use 'classical', 'robust', 'cluster', or 'spatial'.")

    dg = np.diag(vcov)
    with np.errstate(invalid="ignore"):
        se_vec = np.where(dg >= 0, np.sqrt(np.where(dg >= 0, dg, 0.0)), np.nan)

    return dict(coef=beta, se=se_vec, vcov=vcov, resid=resid,
                X_hat=X_hat, sigma2=sigma2, df=effective_df, Q=Q, T=T)


# ══════════════════════════════════════════════════════════════════════════════
#  Anderson–Rubin confidence sets (weak-instrument-robust inference)
# ══════════════════════════════════════════════════════════════════════════════

def _ar_confidence_set(y_dep, x_endog, W_incl, z_cols, se_kwargs, df_adjust,
                       alpha, center=0.0, scale=None,
                       grid_halfwidth=200.0, grid_points=4001):
    """Anderson–Rubin confidence set for the coefficient on a single
    endogenous regressor, by test inversion.

    Model: y_dep = beta * x_endog + (included exogenous)'delta + error, with
    the columns `z_cols` of W_incl the excluded instruments. Under
    H0: beta = b0, the AR regression of (y_dep - b0 * x_endog) on W_incl has
    zero population coefficients on the instruments. The AR statistic is the
    joint Wald statistic on those coefficients, computed under the SAME vcov
    regime as the rest of the analysis (se_kwargs), and is asymptotically
    chi2(K) under H0 REGARDLESS of first-stage strength. The confidence set
    is {b0 : AR(b0) <= chi2_{K, 1-alpha}}.

    K = 1: exact Fieller-type inversion
    -----------------------------------
    Write b(b0) for the instrument coefficient in the AR regression and V(b0)
    for its variance. OLS coefficients and residuals are linear in the
    dependent variable, so b(b0) is exactly LINEAR in b0; and every supported
    variance estimator (classical, HC-robust, cluster, spatial; weights and
    FE enter through the pre-transformed arrays) is a quadratic form in the
    residuals, so V(b0) is exactly QUADRATIC in b0. The acceptance region
    {b0 : b(b0)^2 - cv * V(b0) <= 0} is therefore the solution set of a
    scalar quadratic inequality. We recover the exact polynomial coefficients
    from THREE evaluations of _2sls_with_vcov (reusing the existing, tested
    sandwich code rather than re-deriving any variance formula), verify
    exactness at a fourth point, and solve in closed form. Possible shapes:
      'interval' : bounded [lower, upper]           (strong identification)
      'rays'     : (-inf, left_upper] U [right_lower, +inf)
      'line'     : (-inf, +inf)
    The last two are the weak-instrument signature: the data cannot bound
    the parameter at this confidence level. Degenerate (near-linear)
    quadratics give 'rays_left_only' or 'rays_right_only', and 'empty' can
    only arise from numerically degenerate inputs (the just-identified IV
    estimate attains AR = 0).

    K > 1: closed-form statistic, located on a grid and refined
    ----------------------------------------------------------
    For the same reason, b(b0) = b_y − b0·b_x is linear and the instrument
    block of the variance, V(b0), is a matrix quadratic in b0 for every
    regime, so three evaluations give AR(b0) = b(b0)'V(b0)^{-1}b(b0) exactly
    at every b0 (verified at a fourth point; otherwise the statistic is
    recomputed by brute force). The statistic is evaluated on a grid of
    `grid_points` values centered at `center` spanning +/- grid_halfwidth *
    scale to locate the accepted segments, and every boundary inside the
    grid is refined by Brent's method to a relative tolerance of 1e-10. An
    accepted interval narrower than the grid step is found by minimizing the
    exact statistic on the bracket around every interior local minimum of
    the sampled statistic. The set is reported as the union of accepted
    segments ('grid' shape), with open-endpoint flags when a segment reaches
    the grid boundary. The same path serves as a defensive fallback for
    K = 1 if the exactness check fails (which should not happen; it would
    indicate a non-quadratic vcov).

    Negative variance (spatial SEs; the kernel matrix need not be PSD)
    ------------------------------------------------------------------
    The set is computed for b0 within center ± grid_halfwidth·scale. If the
    variance of the instrument coefficients, V(b0), is not positive definite
    anywhere in that range, the AR statistic is undefined there and the set
    is reported as 'unavailable' (bounds NaN). The test at 0 is NaN when
    V(0) itself is not positive definite.

    Parameters
    ----------
    y_dep, x_endog : arrays (already weighted/demeaned as in the main flow).
    W_incl : full AR regressor matrix INCLUDING the instruments.
    z_cols : list of column indices of the instruments within W_incl.
    se_kwargs : dict forwarded to _2sls_with_vcov (se regime, cluster ids,
        spatial args, vcov cache — the cache stores only design-side
        structures, so reuse across b0 values is safe).
    df_adjust : absorbed degrees of freedom (fixed effects).
    alpha : significance level; the set has asymptotic coverage 1 - alpha.
    center, scale : centering point (2SLS estimate) and scale (its SE) used
        to place evaluation points (K = 1) or the grid (K > 1).

    Returns dict with keys:
      shape            'interval' | 'rays' | 'rays_left_only' |
                       'rays_right_only' | 'line' | 'grid' | 'empty' |
                       'unavailable' (negative variance estimate)
      lower, upper     interval endpoints (NaN unless shape == 'interval')
      left_upper       upper end of the left ray   (shape == 'rays')
      right_lower      lower end of the right ray  (shape == 'rays')
      segments         list of (lo, hi, open_left, open_right) ('grid' only)
      cv, df_chi2, alpha
      stat_at_zero, pval_at_zero   AR test of H0: beta = 0
      contains_zero    bool (None when 'unavailable')
      method           'fieller' or 'grid'
    """
    K = len(z_cols)
    cv = float(stats.chi2.ppf(1.0 - alpha, K))
    z_ix = np.ix_(z_cols, z_cols)

    def _ar_eval(b0):
        """AR regression at b0; returns (beta_Z vector, V_ZZ matrix)."""
        u = y_dep - b0 * x_endog
        res = _2sls_with_vcov(u, W_incl, W_incl, **se_kwargs,
                              df_adjust=df_adjust)
        return (np.atleast_1d(res['coef'][z_cols]).astype(float),
                np.atleast_2d(res['vcov'][z_ix]).astype(float))

    def _ar_stat(b0):
        bz, Vz = _ar_eval(b0)
        return float(bz @ np.linalg.solve(Vz, bz))

    # AR test of H0: beta = 0 (robust companion to the Wald t-test); NaN when
    # the variance estimate at 0 is not positive definite.
    bz0, Vz0 = _ar_eval(0.0)
    if _is_pd(Vz0):
        stat0 = float(bz0 @ np.linalg.solve(Vz0, bz0))
        pval0 = float(stats.chi2.sf(stat0, K))
    else:
        stat0 = pval0 = float("nan")

    out = dict(shape=None, lower=np.nan, upper=np.nan,
               left_upper=np.nan, right_lower=np.nan, segments=None,
               cv=cv, df_chi2=K, alpha=alpha,
               stat_at_zero=stat0, pval_at_zero=pval0,
               contains_zero=None, method=None)

    tc = float(center) if np.isfinite(center) else 0.0
    sc = (float(scale) if scale is not None and np.isfinite(scale)
          and scale > 0 else max(abs(tc), 1.0))
    s_hw = float(grid_halfwidth) * sc       # the range the set is computed on
    # Only the spatial (Conley) sandwich can be indefinite; elsewhere V(b0) is
    # PSD by construction and a sign check would only see rounding.
    check_neg = (se_kwargs.get('se') == 'spatial')

    use_grid = (K > 1)

    if not use_grid:
        # ── Exact Fieller inversion (K = 1) ─────────────────────────────
        d = sc  # evaluation half-step on the natural scale of the problem
        (b_m, V_m), (b_c, V_c), (b_p, V_p) = (
            _ar_eval(tc - d), _ar_eval(tc), _ar_eval(tc + d))
        bm, b0c, bp = float(b_m[0]), float(b_c[0]), float(b_p[0])
        Vm, V0, Vp = float(V_m[0, 0]), float(V_c[0, 0]), float(V_p[0, 0])

        # b(tc + s) = b0c - B*s (exact); V(tc + s) = V0 + bv*s + cvv*s^2 (exact)
        B = (bm - bp) / (2.0 * d)
        bv = (Vp - Vm) / (2.0 * d)
        cvv = (Vp - 2.0 * V0 + Vm) / (2.0 * d * d)

        # q(s) = b(s)^2 - cv * V(s)  <= 0 defines the acceptance region.
        qa = B * B - cv * cvv
        qb = -2.0 * b0c * B - cv * bv
        qc = b0c * b0c - cv * V0

        # Exactness check at a fourth point (2d): linearity of b and
        # quadraticity of V hold analytically; float error should be tiny.
        b_chk, V_chk = _ar_eval(tc + 2.0 * d)
        q_pred = qa * (2.0 * d) ** 2 + qb * (2.0 * d) + qc
        q_act = float(b_chk[0]) ** 2 - cv * float(V_chk[0, 0])
        q_scale = max(abs(qa) * (2 * d) ** 2, abs(qb) * (2 * d), abs(qc),
                      abs(q_act), 1e-300)
        if not np.isfinite(q_pred) or abs(q_pred - q_act) > 1e-6 * q_scale:
            use_grid = True  # defensive fallback; should not trigger
        else:
            out['method'] = 'fieller'
            # V(s) over |s| <= s_hw (its minimum lies at an end, at the vertex
            # or, for the center itself, at s = 0): a negative variance there
            # (non-PSD spatial kernel) leaves the set undefined. Each value is
            # compared with the size of its own terms (rounding scale).
            s_c = [-s_hw, 0.0, s_hw]
            if cvv > 0 and abs(bv / (2.0 * cvv)) <= s_hw:
                s_c.append(-bv / (2.0 * cvv))
            v_neg = check_neg and any(
                V0 + bv * s + cvv * s * s
                < -1e-12 * (abs(V0) + abs(bv * s) + abs(cvv) * s * s) for s in s_c)
            tiny = 1e-14 * max(abs(qb) * d, abs(qc), 1e-300)
            if v_neg:
                out['shape'] = 'unavailable'
            elif abs(qa) * d * d <= tiny:
                # Degenerate: effectively linear in s.
                if abs(qb) * d <= tiny:
                    out['shape'] = 'line' if qc <= 0 else 'empty'
                elif qb > 0:
                    out['shape'] = 'rays_left_only'
                    out['left_upper'] = tc + (-qc / qb)
                else:
                    out['shape'] = 'rays_right_only'
                    out['right_lower'] = tc + (-qc / qb)
            else:
                disc = qb * qb - 4.0 * qa * qc
                if disc < 0.0:
                    # No real roots: q has the sign of qa everywhere.
                    out['shape'] = 'empty' if qa > 0 else 'line'
                else:
                    sq = np.sqrt(disc)
                    # Numerically stable root pair
                    q_ = -0.5 * (qb + np.copysign(sq, qb)) if qb != 0 else 0.5 * sq
                    r1 = q_ / qa
                    r2 = qc / q_ if q_ != 0 else -qb / qa - r1
                    s_lo, s_hi = (min(r1, r2), max(r1, r2))
                    if qa > 0:
                        out['shape'] = 'interval'
                        out['lower'] = tc + s_lo
                        out['upper'] = tc + s_hi
                    else:
                        out['shape'] = 'rays'
                        out['left_upper'] = tc + s_lo
                        out['right_lower'] = tc + s_hi

    if use_grid:
        # ── Closed-form statistic on a grid, boundaries refined (K > 1) ───
        out['method'] = 'grid'
        grid = tc + np.linspace(-grid_halfwidth, grid_halfwidth,
                                int(grid_points)) * sc
        n_g = len(grid)
        d = sc
        (b_m, V_m), (b_c, V_c), (b_p, V_p) = (
            _ar_eval(tc - d), _ar_eval(tc), _ar_eval(tc + d))
        Bv = (b_m - b_p) / (2.0 * d)
        V1 = (V_p - V_m) / (2.0 * d)
        V2 = (V_p - 2.0 * V_c + V_m) / (2.0 * d * d)

        def _closed_bV(b0s):
            s_ = np.atleast_1d(np.asarray(b0s, dtype=float)) - tc
            bb = b_c[None, :] - s_[:, None] * Bv[None, :]
            VV = V_c[None] + s_[:, None, None] * V1[None] + (s_ ** 2)[:, None, None] * V2[None]
            return bb, VV

        def _closed(b0s):
            bb, VV = _closed_bV(b0s)
            with np.errstate(all="ignore"):
                try:
                    sol = np.linalg.solve(VV, bb[..., None])[..., 0]
                except np.linalg.LinAlgError:
                    return np.array([_ar_stat(b) for b in np.atleast_1d(b0s)])
            return np.sum(bb * sol, axis=1)

        chk_b = tc + 2.0 * d
        chk_act = _ar_stat(chk_b)
        chk_pred = float(_closed(chk_b)[0])
        exact = (np.all(np.isfinite([chk_act, chk_pred]))
                 and abs(chk_pred - chk_act) <= 1e-6 * max(abs(chk_act), 1.0))
        if exact:
            bb_g, VV_g = _closed_bV(grid)
            with np.errstate(all="ignore"):
                ev_min = np.linalg.eigvalsh((VV_g + np.swapaxes(VV_g, 1, 2)) / 2.0)
            stat_grid = _closed(grid)
            f_ = lambda b: float(_closed(b)[0]) - cv
        else:
            evals = [_ar_eval(b0) for b0 in grid]
            with np.errstate(all="ignore"):
                ev_min = np.array([np.linalg.eigvalsh((Vz + Vz.T) / 2.0) for _, Vz in evals])
            stat_grid = np.array([float(bz @ np.linalg.solve(Vz, bz)) if _is_pd(Vz)
                                  else np.nan for bz, Vz in evals])
            f_ = lambda b: _ar_stat(b) - cv
        # A variance block that is not positive definite anywhere on the grid
        # (non-PSD spatial kernel) leaves the set undefined. The smallest
        # eigenvalue at each grid point is compared with the size of V there.
        if check_neg:
            with np.errstate(invalid="ignore"):
                lam_lo = np.min(ev_min, axis=1)
                lam_sc = np.max(np.abs(ev_min), axis=1)
                unavailable = bool(np.any(~np.isfinite(lam_lo))
                                   or np.any(lam_lo < -1e-12 * lam_sc))
        else:
            unavailable = False
        accept = stat_grid <= cv

        def _root(a_, b_):
            fa, fb = f_(a_), f_(b_)
            if not (np.isfinite(fa) and np.isfinite(fb)) or fa * fb > 0:
                return None
            xtol = 1e-14 * max(abs(a_), abs(b_), sc)
            return float(brentq(f_, a_, b_, xtol=xtol, rtol=1e-10, maxiter=200))

        if unavailable:
            out['shape'] = 'unavailable'
        else:
            segments = []
            i = 0
            while i < n_g:
                if accept[i]:
                    j = i
                    while j + 1 < n_g and accept[j + 1]:
                        j += 1
                    lo = float(grid[i])
                    hi = float(grid[j])
                    if i > 0:
                        r_ = _root(float(grid[i - 1]), float(grid[i]))
                        lo = r_ if r_ is not None else lo
                    if j < n_g - 1:
                        r_ = _root(float(grid[j]), float(grid[j + 1]))
                        hi = r_ if r_ is not None else hi
                    segments.append((lo, hi, i == 0, j == n_g - 1))
                    i = j + 1
                else:
                    i += 1
            # An accepted interval narrower than the grid step lies between
            # two rejected grid points: minimize the statistic exactly on the
            # bracket around every interior local minimum of the sampled
            # statistic (lowest first, at most 50) and add what dips below cv.
            if n_g >= 3:
                ii = np.arange(1, n_g - 1)
                with np.errstate(invalid="ignore"):
                    is_min = ((stat_grid[ii] <= stat_grid[ii - 1])
                              & (stat_grid[ii] <= stat_grid[ii + 1])
                              & ~accept[ii] & np.isfinite(stat_grid[ii]))
                cand = ii[is_min]
                cand = cand[np.argsort(stat_grid[cand], kind="stable")][:50]
                for i_ in cand:
                    a_, b_ = float(grid[i_ - 1]), float(grid[i_ + 1])
                    xatol = 1e-12 * max(abs(a_), abs(b_), sc)
                    res_m = minimize_scalar(f_, bounds=(a_, b_), method="bounded",
                                            options=dict(xatol=xatol, maxiter=500))
                    if np.isfinite(res_m.fun) and res_m.fun < 0:
                        lo_ = _root(a_, float(res_m.x))
                        hi_ = _root(float(res_m.x), b_)
                        if lo_ is not None and hi_ is not None:
                            segments.append((lo_, hi_, False, False))
                segments.sort(key=lambda s_: s_[0])
            out['shape'] = 'grid' if segments else 'empty'
            out['segments'] = segments

    # ── Does the set contain zero? ───────────────────────────────────────
    shape = out['shape']
    if shape == 'unavailable':
        out['contains_zero'] = None
        return out
    if shape == 'interval':
        out['contains_zero'] = bool(out['lower'] <= 0.0 <= out['upper'])
    elif shape == 'rays':
        out['contains_zero'] = bool(0.0 <= out['left_upper']
                                    or 0.0 >= out['right_lower'])
    elif shape == 'rays_left_only':
        out['contains_zero'] = bool(0.0 <= out['left_upper'])
    elif shape == 'rays_right_only':
        out['contains_zero'] = bool(0.0 >= out['right_lower'])
    elif shape == 'line':
        out['contains_zero'] = True
    elif shape == 'grid':
        out['contains_zero'] = any(lo <= 0.0 <= hi
                                   for lo, hi, _, _ in out['segments'])
    else:  # 'empty'
        out['contains_zero'] = False
    # Consistency: the AR test at zero must agree with set membership
    # (up to boundary ties); prefer the direct test when they disagree.
    if np.isfinite(stat0):
        test_says_in = bool(stat0 <= cv)
        if test_says_in != out['contains_zero']:
            out['contains_zero'] = test_says_in
    return out


def _ar_unbounded_or_disconnected(ar):
    """True when an AR set is unbounded (rays, line, or a grid segment that
    reaches the grid boundary) or has more than one segment."""
    if ar is None:
        return False
    shape = ar['shape']
    if shape in ('rays', 'rays_left_only', 'rays_right_only', 'line'):
        return True
    if shape == 'grid':
        segs = ar.get('segments') or []
        return len(segs) > 1 or any(o_l or o_r for _, _, o_l, o_r in segs)
    return False


def _format_ar_set(ar, fmt="{:.4f}"):
    """Human-readable rendering of an AR confidence set dict."""
    f = fmt.format
    shape = ar['shape']
    if shape == 'interval':
        return f"[{f(ar['lower'])}, {f(ar['upper'])}]"
    if shape == 'rays':
        return (f"(-inf, {f(ar['left_upper'])}] U "
                f"[{f(ar['right_lower'])}, +inf)   <- unbounded: weak-IV signature")
    if shape == 'rays_left_only':
        return f"(-inf, {f(ar['left_upper'])}]   <- unbounded: weak-IV signature"
    if shape == 'rays_right_only':
        return f"[{f(ar['right_lower'])}, +inf)   <- unbounded: weak-IV signature"
    if shape == 'line':
        return ("(-inf, +inf)   <- AR cannot reject any value at this level "
                "(weak-IV signature)")
    if shape == 'grid':
        parts = []
        for lo, hi, open_l, open_r in ar['segments']:
            lo_s = "< " + f(lo) if open_l else f(lo)
            hi_s = "> " + f(hi) if open_r else f(hi)
            parts.append(f"[{lo_s}, {hi_s}]")
        note = ("   (grid search; '<' / '>' mark segments that reach the "
                "search grid's boundary)" if any(o_l or o_r for _, _, o_l, o_r
                                                in ar['segments']) else "")
        return " U ".join(parts) + note
    if shape == 'unavailable':
        return "unavailable: negative variance estimate"
    return "empty set (AR rejects every value: model misspecification signal)"


# ══════════════════════════════════════════════════════════════════════════════
#  Main Analysis Function
# ══════════════════════════════════════════════════════════════════════════════

# ──────────────────────────────────────────────────────────────────────────────
#  Variable-list parsing
#
#  For convenience, list-typed arguments (X, dummies) accept either a Python
#  list of strings or a single whitespace-separated string. The string form
#  also supports shell-style globs (* ? [...]) which expand against the columns
#  of the dataframe. This avoids the manual quote-and-comma overhead for long
#  control sets and lets users write e.g.,
#      X = "l_sh_routine33 l_task_outsource reg_*"
#  or, for paste-friendly multiline specs:
#      X = '''l_sh_routine33 l_task_outsource
#             reg_midatl reg_encen reg_wncen
#             reg_satl reg_escen reg_wscen reg_mount reg_pacif'''
#  Y, T, M, Z, fe, and clusterid keep their existing behavior: the user must
#  spell out each variable explicitly, which is the right default for these
#  high-stakes arguments where an accidental glob could change the structural
#  specification.
# ──────────────────────────────────────────────────────────────────────────────

_GLOB_CHARS = frozenset("*?[")


def _resolve_var_list(spec, df, *, name="<arg>", exclude=(), forbid=()):
    """Normalize a variable specification to a list of column names in df.

    Accepts:
      - None or empty string  → []
      - list/tuple of strings  → each element is a literal name or a glob
      - whitespace-separated string → tokenized the same way
    Tokens with glob characters (* ? [) are expanded against df.columns (list
    elements included); literal tokens are checked for existence.

    Role protection (2026-10): columns in `exclude` (the outcome, treatment,
    mediator, instruments and the fe/cluster/weight/coordinate columns) are
    never produced by glob expansion; a literal token that names a column in
    `forbid` (Y, T, M, Z) raises ValueError.

    Duplicates are removed preserving first occurrence. Globs that match zero
    columns and literals that aren't in df.columns raise ValueError with a
    message identifying `name`.
    Returns (resolved_list, excluded_by_glob_list) when `exclude` is given,
    otherwise just the resolved list.
    """
    exclude = set(exclude or ())
    forbid = set(forbid or ())
    excluded = []
    if spec is None:
        return ([], excluded) if exclude else []
    if isinstance(spec, str):
        s = spec.strip()
        if not s:
            return ([], excluded) if exclude else []
        tokens = re.split(r"\s+", s)
    else:
        tokens = list(spec)

    columns = list(df.columns)
    column_set = set(columns)
    resolved = []
    for tok in tokens:
        if not isinstance(tok, str):
            raise ValueError(f"{name}: variable names must be strings; got {tok!r}.")
        if any(c in tok for c in _GLOB_CHARS):
            matches = fnmatch.filter(columns, tok)
            if not matches:
                raise ValueError(
                    f"{name}: glob pattern {tok!r} matched 0 columns. "
                    f"Available columns include: "
                    f"{', '.join(columns[:8])}{'...' if len(columns) > 8 else ''}"
                )
            kept = [m for m in matches if m not in exclude]
            excluded.extend(m for m in matches if m in exclude and m not in excluded)
            if not kept:
                raise ValueError(
                    f"{name}: glob pattern {tok!r} matched only role columns "
                    f"({', '.join(matches)}), which cannot be controls."
                )
            resolved.extend(kept)
        else:
            if tok not in column_set:
                raise ValueError(
                    f"{name}: variable {tok!r} not found in data. "
                    f"Available columns include: "
                    f"{', '.join(columns[:8])}{'...' if len(columns) > 8 else ''}"
                )
            resolved.append(tok)

    bad = [v for v in resolved if v in forbid]
    if bad:
        raise ValueError(
            f"{name} must not contain the outcome, treatment, mediator or "
            f"instruments: {', '.join(dict.fromkeys(bad))}."
        )

    # Defensive: detect columns that are entirely NaN. Without this check, the
    # downstream listwise deletion would silently drop every row.
    all_nan = [v for v in resolved if df[v].isna().all()]
    if all_nan:
        raise ValueError(
            f"{name}: column(s) {all_nan!r} contain only missing values. "
            f"Either drop them from {name} or fill the missing values before "
            f"calling ivmediate."
        )

    # Deduplicate, preserving first occurrence
    out = list(dict.fromkeys(resolved))
    return (out, excluded) if exclude else out


def _summarize_resolved(name, resolved, max_show=8):
    """Format a one-line note for the resolved variable list (used in print output)."""
    if not resolved:
        return None
    n = len(resolved)
    if n <= max_show:
        body = ", ".join(resolved)
    else:
        body = ", ".join(resolved[:max_show]) + f", ... (+{n - max_show} more)"
    return f"  Resolved {name} ({n} var{'s' if n != 1 else ''}): {body}"


# ──────────────────────────────────────────────────────────────────────────────
#  Option validation, sample preparation and identification helpers (2026-10)
# ──────────────────────────────────────────────────────────────────────────────

_SE_TYPES = ("classical", "robust", "cluster", "spatial")


def _check_bool(name, value, allow_none=False):
    """Booleans must be real booleans (a string such as "false" is an error)."""
    if allow_none and value is None:
        return
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be True or False; got {value!r}.")


def _check_open_unit(name, value):
    """alpha, op_tau, op_alpha and similar must lie strictly between 0 and 1."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a number in (0, 1); got {value!r}.")
    if not (np.isfinite(v) and 0.0 < v < 1.0):
        raise ValueError(f"{name} must lie strictly between 0 and 1; got {value!r}.")


def _check_choice(name, value, choices):
    if value not in choices:
        raise ValueError(f"{name} must be one of {', '.join(repr(c) for c in choices)}; got {value!r}.")


def _level_label(alpha):
    """Confidence-level label used everywhere: 0.05 → '95%', 0.025 → '97.5%'."""
    return f"{round(100.0 * (1.0 - float(alpha)), 10):g}%"


def _normalize_Z(Z):
    """Z as a string, list or tuple of distinct column names."""
    if isinstance(Z, str):
        Zl = [Z]
    elif isinstance(Z, (list, tuple)):
        Zl = list(Z)
    else:
        raise ValueError(f"Z must be a column name or a list/tuple of names; got {Z!r}.")
    if not Zl:
        raise ValueError("Z must name at least one instrument.")
    for z in Zl:
        if not isinstance(z, str):
            raise ValueError(f"instrument names must be strings; got {z!r}.")
    dup = [z for i, z in enumerate(Zl) if z in Zl[:i]]
    if dup:
        raise ValueError(f"duplicate instrument names: {', '.join(dict.fromkeys(dup))}.")
    return Zl


def _missing_text_to_nan(series):
    """Empty and whitespace-only strings in a text id column count as missing."""
    if series.dtype == object or pd.api.types.is_string_dtype(series):
        blank = series.map(lambda v: isinstance(v, str) and v.strip() == "")
        if blank.any():
            series = series.mask(blank)
    return series


def _prepare_sample(data, Y, T, M, Z, X, dummies, fe, se, clusterid, weights,
                    lat, lon, lasso=False):
    """Resolve and validate variable roles, normalize missing values and return
    the estimation sample (before singleton removal and FE absorption).

    Steps: role checks (Y, T, M, Z distinct and present; Z without duplicates;
    X/dummies never take a role column by glob, and an explicit Y/T/M/Z among
    them is an error); text id columns ('' and whitespace-only are missing);
    listwise deletion over every used column (cluster ids only with
    se='cluster', coordinates only with se='spatial'); ±Inf in a numeric role
    is an error; negative weights are an error; zero-weight rows are dropped.

    Returns a dict with the cleaned frame `df` (used columns only, a copy),
    the boolean row mask over `data`, the resolved lists and printable notes.
    """
    notes = []
    for nm, v in (("Y", Y), ("T", T), ("M", M)):
        if not isinstance(v, str):
            raise ValueError(f"{nm} must be a column name; got {v!r}.")
    Zl = _normalize_Z(Z)
    if len({Y, T, M}) < 3:
        raise ValueError("Y, T and M must be three different columns.")
    clash = [z for z in Zl if z in (Y, T, M)]
    if clash:
        raise ValueError(f"instruments must differ from the outcome, treatment and mediator: {', '.join(clash)}.")
    fe_cols = None
    if fe is not None:
        fe_cols = [fe] if isinstance(fe, str) else list(fe)
        if len(fe_cols) == 0:
            fe_cols = None
        elif len(fe_cols) > 2:
            raise ValueError("At most two-way fixed effects are supported.")
    cl_eff = clusterid if se == "cluster" else None
    roles = [Y, T, M] + Zl
    missing_cols = [c for c in roles + (fe_cols or []) + ([cl_eff] if cl_eff else []) +
                    ([weights] if weights else []) +
                    ([lat, lon] if se == "spatial" else []) if c not in data.columns]
    if missing_cols:
        raise ValueError(f"column(s) not found in data: {', '.join(dict.fromkeys(missing_cols))}.")

    aux = [c for c in (fe_cols or []) + [clusterid, weights, lat, lon] if c]
    exclude = set(roles) | set(aux)
    X_res, ex_x = _resolve_var_list(X, data, name="X", exclude=exclude, forbid=set(roles))
    D_res, ex_d = _resolve_var_list(dummies, data, name="dummies", exclude=exclude, forbid=set(roles))
    for nm, ex in (("X", ex_x), ("dummies", ex_d)):
        if ex:
            notes.append(f"  Note: {nm} globs also matched role columns, excluded: {', '.join(ex)}.")

    numeric_roles = list(dict.fromkeys(roles + X_res + ([weights] if weights else []) +
                                       ([lat, lon] if se == "spatial" else [])))
    id_cols = list(dict.fromkeys((fe_cols or []) + ([cl_eff] if cl_eff else []) + D_res))
    used = list(dict.fromkeys(numeric_roles + id_cols))
    df = data.loc[:, used].copy()
    for c in id_cols:
        df[c] = _missing_text_to_nan(df[c])
    for c in numeric_roles:
        if not pd.api.types.is_numeric_dtype(df[c]):
            ex = df[c].dropna().astype(str).unique()[:3]
            raise ValueError(f"column {c!r} must be numeric; found values such as {', '.join(map(repr, ex))}.")
    mask = df.notna().all(axis=1).to_numpy()
    num = df.loc[mask, numeric_roles].to_numpy(dtype=float)
    inf_cols = [numeric_roles[j] for j in np.where(np.isinf(num).any(axis=0))[0]]
    if inf_cols:
        raise ValueError(f"infinite values in column(s): {', '.join(inf_cols)}.")
    n_zero = 0
    if weights:
        wv = df[weights].to_numpy(dtype=float)
        if np.any(wv[mask] < 0):
            raise ValueError("Weights must be non-negative.")
        zero = mask & (wv == 0)
        n_zero = int(zero.sum())
        if n_zero:
            if n_zero == int(mask.sum()):
                raise ValueError("All weights are zero — no observations have positive "
                                 "weight to estimate from.")
            mask = mask & ~zero
    df = df.loc[mask].copy()
    if len(df) == 0:
        raise ValueError("No observations remaining after dropping missing values.")
    return dict(df=df, mask=mask, Z=Zl, X=X_res, dummies=D_res, fe_cols=fe_cols,
                clusterid=cl_eff, notes=notes, n_zero_weights=n_zero)


def _greedy_independent(A, names, tol=1e-9, keep_first=None, ref_norms=None):
    """Columns of A that are linearly independent, tested in order.

    Each column is divided by its reference scale (`ref_norms[j]`, default its
    own norm) and kept when the norm of its residual on the previously kept
    columns exceeds `tol`, i.e., |R_jj| < 1e-9·max|R_ii| of an equilibrated QR
    marks a redundant column. With FE-demeaned columns the reference scale is
    the column's variation before the FE were absorbed, so a control that the
    FE absorb (rounding noise afterward) is recognized. Earlier columns win,
    so the later member of a collinear set is dropped (as in Stata's _rmcoll
    and R's lm). `keep_first` columns (e.g., the intercept) are always kept and
    orthogonalized first.
    Returns (kept_indices, dropped_names).
    """
    A = np.asarray(A, dtype=float)
    n, p = A.shape
    basis = []
    if keep_first is not None:
        K0 = np.asarray(keep_first, dtype=float).reshape(n, -1)
        for j in range(K0.shape[1]):
            v = K0[:, j] / max(np.linalg.norm(K0[:, j]), 1e-300)
            for b in basis:
                v = v - (b @ v) * b
            nv = np.linalg.norm(v)
            if nv > tol:
                basis.append(v / nv)
    kept, dropped = [], []
    for j in range(p):
        ref = float(ref_norms[j]) if ref_norms is not None else float(np.linalg.norm(A[:, j]))
        if not np.isfinite(ref) or ref <= 0.0:
            dropped.append(names[j])
            continue
        v = A[:, j] / ref
        for _ in range(2):                 # two passes of Gram–Schmidt for accuracy
            for b in basis:
                v = v - (b @ v) * b
        nv = np.linalg.norm(v)
        if np.isfinite(nv) and nv > tol:
            basis.append(v / nv)
            kept.append(j)
        else:
            dropped.append(names[j])
    return kept, dropped


def _fe_degrees_of_freedom(df, fe_cols, cluster_col=None):
    """Absorbed FE degrees of freedom beyond the intercept (rank of the dummies).

    One-way: G − 1. Two-way: G₁ + G₂ − c − 1, with c the number of connected
    components of the bipartite FE graph (handles nested and disconnected
    designs). With cluster SEs, an FE dimension nested in the cluster variable
    (every FE group lies within one cluster) does not count toward the CR1
    correction (reghdfe/fixest convention). With exactly one dimension nested,
    the other counts G_other − 1, whatever the connectivity of the FE graph:
    reghdfe drops the nested FE from the redundancy (mobility-group) count and
    pyfixest's default ssc counts the same.

    Returns (df_absorbed, df_for_cr1, n_groups_total, nested_flags).
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    codes = [pd.factorize(df[c], sort=False)[0] for c in fe_cols]
    G = [int(cd.max()) + 1 for cd in codes]
    nested = []
    if cluster_col is not None:
        cl = pd.factorize(df[cluster_col], sort=False)[0]
        for cd in codes:
            pairs = np.unique(np.column_stack([cd, cl]), axis=0)
            nested.append(bool(len(pairs) == len(np.unique(pairs[:, 0]))))
    else:
        nested = [False] * len(codes)
    if len(codes) == 1:
        d = G[0] - 1
        d_cr1 = 0 if nested[0] else d
    else:
        G1, G2 = G
        edges = np.unique(np.column_stack([codes[0], codes[1]]), axis=0)
        adj = csr_matrix((np.ones(len(edges)), (edges[:, 0], G1 + edges[:, 1])),
                         shape=(G1 + G2, G1 + G2))
        c = int(connected_components(adj, directed=False)[0])
        d = G1 + G2 - c - 1
        if nested[0] and nested[1]:
            d_cr1 = 0
        elif nested[0]:
            d_cr1 = G2 - 1          # reghdfe / pyfixest: one redundant level
        elif nested[1]:
            d_cr1 = G1 - 1
        else:
            d_cr1 = d               # exact count with c components (reghdfe)
    return d, d_cr1, int(sum(G)), nested


def _op_keff_critval(V_ZZ, Q_zz, tau, alpha):
    """Montiel Olea–Pflueger (2013) simplified test at the estimated K_eff.

    Ŵ₂ = Q^{1/2} V̂_π Q^{1/2}; K_eff = [tr Ŵ₂]²(1 + 2x) / (tr(Ŵ₂′Ŵ₂) +
    2x·tr(Ŵ₂)·λ_max(Ŵ₂)), x = 1/τ; critical value = noncentral-χ² quantile
    χ²_{K_eff}(x·K_eff)⁻¹(1 − α) / K_eff. K = 1 gives K_eff = 1 (23.109 at τ = 0.10).
    Returns (K_eff, critical_value); NaNs when the variance is unavailable.
    """
    V = np.atleast_2d(np.asarray(V_ZZ, dtype=float))
    Q = np.atleast_2d(np.asarray(Q_zz, dtype=float))
    x = 1.0 / float(tau)
    K = V.shape[0]
    if K == 1:
        keff = 1.0
    else:
        if not (_is_pd(V) and np.all(np.isfinite(Q))):
            return float("nan"), float("nan")
        w, U = np.linalg.eigh((Q + Q.T) / 2.0)
        Qh = (U * np.sqrt(np.maximum(w, 0.0))) @ U.T
        W2 = Qh @ V @ Qh
        W2 = (W2 + W2.T) / 2.0
        tr = float(np.trace(W2))
        lmax = float(np.max(np.linalg.eigvalsh(W2)))
        den = float(np.trace(W2 @ W2)) + 2.0 * x * tr * lmax
        if not (np.isfinite(den) and den > 0 and tr > 0):
            return float("nan"), float("nan")
        keff = tr * tr * (1.0 + 2.0 * x) / den
    crit = float(stats.ncx2.ppf(1.0 - float(alpha), keff, x * keff) / keff)
    return float(keff), crit


def _residualize(A, B):
    """Residuals of the columns of A on B (QR of the equilibrated B)."""
    if B is None or B.shape[1] == 0:
        return A.copy()
    Q, _, _ = _equilibrated_qr(B)
    return A - Q @ (Q.T @ A)


# ──────────────────────────────────────────────────────────────────────────────
#  Rich-expansion helper for post-double-selection LASSO
#
#  Builds a high-dimensional candidate set W by applying user-controlled
#  operations to the raw X variables and dummy indicators:
#    - "square", "cube" for continuous variables
#    - "interact_with_X"  pairwise interactions with continuous controls
#    - "interact_with_dummies" pairwise interactions with other dummies
#  Indicators of the same categorical variable are not interacted with each
#  other (they are mutually exclusive, so the interaction is identically zero).
#  High-cardinality FE_VAR columns are absorbed by demeaning and never enter
#  the rich expansion.
#
#  The expansion is built on raw level-form data, BEFORE FE demeaning and
#  weight pre-multiplication, because squares/cubes/interactions are
#  level-based concepts. Demeaning and weighting are applied to the resulting
#  matrix downstream.
# ──────────────────────────────────────────────────────────────────────────────

_VALID_EXPAND_OPS = frozenset(["square", "cube", "interact_with_X", "interact_with_dummies"])


def _is_continuous(values):
    """Heuristic: a variable is continuous if it has more than 2 unique values.
    Indicator variables (0/1) are correctly classified as non-continuous."""
    return pd.unique(values).size > 2


# ──────────────────────────────────────────────────────────────────────────────
#  BCH plug-in LASSO solver
#
#  Implements the iterative plug-in penalty of Belloni-Chen-Chernozhukov-Hansen
#  (2012, Econometrica) and BCH (2014, REStud). The penalty is:
#
#      λ = 2 c √n_eff Φ⁻¹(1 - α_p / (2p)),
#
#  where n_eff is the number of observations n, or the Kish effective size
#  (Σw)²/Σw² with unequal weights — also under clustering and spatial
#  dependence, where the dependence enters only through the loadings ψ_j
#  (Belloni, Chernozhukov, Hansen and Kozbur 2016, eqs. 2.4 and 2.6; lassopack).
#  The default α_p = 0.1/log(n) uses the number of observations in every
#  regime; lassopack's cluster default is 0.1/log(G).
#
#  Per-feature loadings ψ_j are computed from the current residuals û, with
#  the score w_i W_ij û_i of √w-scaled data (lassopack):
#    homosk./robust:    ψ_j = √((1/n_w) Σ_i w_i² W_ij² û_i²)
#    cluster:           ψ_j = √((1/n_w) Σ_c (Σ_{i∈c} w_i W_ij û_i)²)
#  with n_w = Σ w_i. The residuals come from a post-LASSO OLS of y on the
#  selected columns, the keep columns and an intercept (hdm/lassopack default).
#
#  Implementation (2026-10): weights are normalized to mean 1 and y and W are
#  centered at their weighted means (an unpenalized intercept), so the
#  selection does not depend on the units of the weights or on the location
#  of any variable. Then pre-multiply (y, W) by √w_i (so the squared loss is
#  weighted, loss = (1/(2n_w)) Σ w_i(y_i − W_i'β)² as in App. H eq. 179),
#  pre-divide each column of W by ψ_j (so the penalty is uniform), feed to
#  sklearn's standard Lasso, then unscale β by ψ_j to recover coefficients on
#  the original (un-rescaled, un-pre-multiplied) W.
# ──────────────────────────────────────────────────────────────────────────────

def _compute_psi(W, u, weights, cluster_ids, spatial_args=None):
    """Per-feature loadings ψ_j based on current residuals.

    Always uses raw (un-pre-multiplied) W and û with explicit weights.

    Parameters
    ----------
    W : (n, p) ndarray, raw demeaned design (un-pre-multiplied)
    u : (n,) ndarray, current residual estimate (raw scale)
    weights : (n,) ndarray of strictly-positive weights (1.0 if unweighted)
    cluster_ids : (n,) array or None
    spatial_args : dict or None
        If given, ψ_j is computed under spatial dependence using a
        kernel-HAC estimator. Required keys:
          - 'indices': list of length n; indices[i] is the array of
                       neighbor indices for observation i (within cutoff).
          - 'distances': list of length n; matching distances (km or raw).
          - 'kernel': "uniform" or "bartlett".
          - 'cutoff': float, the bandwidth.
        Mutually exclusive with cluster_ids: if both are given, spatial
        takes precedence (cluster structure is then absorbed into the
        spatial-dependence assumption via the cutoff being large enough
        to encompass each cluster).

    Returns
    -------
    psi : (p,) ndarray of loadings, all strictly positive

    Formulas (the score of observation i is w_i W_ij û_i, as with √w-scaled
    data; lassopack's weighted loadings):
        homoskedastic / robust:
            ψ_j = √((1/n_w) Σ_i w_i² W_ij² û_i²)
        cluster:
            ψ_j = √((1/n_w) Σ_c (Σ_{i∈c} w_i W_ij û_i)²)
        spatial (Conley):
            ψ_j = √((1/n_w) Σ_i Σ_{j∈N(i)} K(d_{ij}) w_i w_j W_ij W_jj' û_i û_j)
        where N(i) is the set of within-cutoff neighbors of i.
    """
    n_w = float(weights.sum())
    if spatial_args is not None:
        # Spatial: kernel-weighted score covariance.
        indices = spatial_args['indices']
        distances = spatial_args['distances']
        kernel = spatial_args['kernel']
        cutoff = spatial_args['cutoff']
        wu = weights * u                      # (n,)
        wuW = W * wu[:, None]                 # (n, p), per-obs score row
        n, p = W.shape

        # Try the vectorized (sparse-matrix) path first. If memory budget is
        # exceeded, fall back to the per-row loop. Note: the per-feature
        # diagonal psi_j² is given by the diagonal of (wuW)' A (wuW), which
        # we compute as a column-wise sum-product to avoid materializing the
        # full p×p meat matrix.
        #
        # Cache the adjacency in spatial_args: it depends only on
        # (indices, distances, kernel, cutoff), which are constant across
        # the 5-15 plug-in iterations × 4 LASSOs of a quadruple-LASSO call.
        # A 'sentinel' value of False indicates a previous attempt to build
        # the adjacency hit the memory budget; we don't retry.
        if 'adjacency' in spatial_args:
            A = spatial_args['adjacency']
        else:
            A = _build_spatial_sparse_adjacency(indices, distances, kernel,
                                                 cutoff, n)
            spatial_args['adjacency'] = A    # may be None (memory fallback)
        if A is not None:
            # Vectorized: psi_j² = Σ_i wuW[i, j] · (A @ wuW)[i, j]
            AwuW = A @ wuW                     # (n, p)
            psi2 = (wuW * AwuW).sum(axis=0)    # (p,)
        else:
            # Memory-pressure fallback (loop-based)
            psi2 = np.zeros(p)
            for i in range(n):
                nbr_idx = indices[i]
                if len(nbr_idx) == 0:
                    continue
                if kernel == "uniform":
                    w_k = np.ones(len(nbr_idx))
                elif kernel == "bartlett":
                    w_k = np.maximum(1.0 - distances[i] / cutoff, 0.0)
                else:
                    raise ValueError(f"Unknown spatial kernel {kernel!r}")
                nbr_score_sum = (w_k[:, None] * wuW[nbr_idx]).sum(axis=0)
                psi2 += wuW[i] * nbr_score_sum
        # Spatial meats can be slightly negative element-wise in finite
        # samples (PSD failure, same pathology as in _vcov_spatial). Clip
        # at zero before taking sqrt and floor below.
        psi2 = np.maximum(psi2, 0.0)
        psi = np.sqrt(psi2 / n_w)
    elif cluster_ids is None:
        # Heteroskedastic form: ψ_j = √((1/n_w) Σ_i (w_i W_ij û_i)²), the
        # same score as the cluster and spatial branches (clusters of size one
        # and a self-only spatial kernel give exactly these loadings).
        wu2 = (weights * u) ** 2
        psi = np.sqrt((W * W * wu2[:, None]).sum(axis=0) / n_w)
    else:
        # Cluster form: ψ_j = √((1/n_w) Σ_c (Σ_{i∈c} w_i W_ij û_i)²)
        wu = weights * u
        wuW = W * wu[:, None]   # (n, p), each row contributes w_i û_i W_i
        inverse, unique_clusters = pd.factorize(np.asarray(cluster_ids), sort=False)
        G = len(unique_clusters)
        cluster_sums = np.empty((G, W.shape[1]))
        for _j in range(W.shape[1]):
            cluster_sums[:, _j] = np.bincount(inverse, weights=wuW[:, _j],
                                              minlength=G)
        psi = np.sqrt((cluster_sums ** 2).sum(axis=0) / n_w)

    # Numerical safety: floor at a tiny value so column rescaling is well-defined.
    psi = np.maximum(psi, 1e-10 * float(np.abs(W).max(initial=1.0)))
    return psi


def _compute_n_eff(weights, cluster_ids, spatial_args=None):
    """Effective sample size for the BCH penalty level λ = 2c√n_eff Φ⁻¹(·).

    n_eff = n, or the Kish size (Σw)²/Σw² with unequal weights, in every
    regime. Under clustering (and spatial dependence) the dependence enters
    only through the loadings ψ_j, which are averaged over observations
    (1/n_w); pairing them with √n follows Belloni, Chernozhukov, Hansen and
    Kozbur (2016, eqs. 2.4 and 2.6) and lassopack. (Before 2026-10 the cluster
    case used √G, a penalty √(G/n) too small.) `cluster_ids` and
    `spatial_args` are accepted for interface stability and do not change n_eff.
    """
    weights = np.asarray(weights, dtype=float)
    if np.allclose(weights, weights[0]):
        return float(len(weights))
    n_w = float(weights.sum())
    return n_w * n_w / float((weights ** 2).sum())


def _bch_lambda(p, n_eff, c_slack, alpha_p):
    """BCH plug-in penalty level: λ = 2c √n_eff Φ⁻¹(1 - α_p/(2p))."""
    if p == 0:
        return 0.0
    quantile = stats.norm.ppf(1.0 - alpha_p / (2.0 * p))
    return 2.0 * c_slack * np.sqrt(n_eff) * quantile


def _post_lasso_resid(y_c, W_c, sel, sqrt_w):
    """Residuals of a weighted OLS of y on the selected columns of W (y and W
    are already centered and have the keep columns partialed out, so this is
    the post-LASSO regression on [selected, keep, intercept] by FWL)."""
    if not np.any(sel):
        return y_c.copy()
    A = W_c[:, sel] * sqrt_w[:, None]
    coef, *_ = np.linalg.lstsq(A, y_c * sqrt_w, rcond=None)
    return y_c - W_c[:, sel] @ coef


def _bch_lasso(
    y_raw, W_raw, weights, cluster_ids,
    n_iter_max=15, tol=1e-6, c_slack=1.1, alpha_p=None,
    lambda_method="plugin", cv_folds=5, cv_seed=0,
    keep_idx=None, lambda_override=None,
    spatial_args=None,
    verbose=False,
):
    """Iterative plug-in LASSO with weight + cluster + spatial support.

    See the comment block above for the formulas. Parameters:

    y_raw, W_raw : LHS / design, FE-demeaned if FE are absorbed, NOT
        √w-pre-multiplied. Both are centered at their weighted means here (an
        unpenalized intercept), so the location of any variable does not
        matter.
    weights      : (n,) positive weights; pass np.ones(n) if unweighted. They
        are normalized to mean 1 inside the LASSO, so their scale does not
        matter.
    cluster_ids  : (n,) array or None (cluster-form loadings).
    spatial_args : dict or None. If given, the loadings use a kernel-HAC
        score-variance estimator. Required keys (built once by
        _build_spatial_neighbors): indices, distances, kernel, cutoff.
        spatial_args takes precedence over cluster_ids if both are passed.
        n_eff stays n (or Kish n).
    keep_idx     : column indices that should always appear in the selection
                   (partialed out before LASSO; re-added unconditionally).
    lambda_method: "plugin" (default) or "cv".
    lambda_override : if not None, fixes λ at this value (skips plug-in/CV).

    Plug-in iterations update the loadings from post-LASSO residuals (OLS of y
    on the selected columns, the keep columns and an intercept), as hdm and
    lassopack do by default; they stop when max|Δβ| < tol or after
    n_iter_max iterations.

    Returns a dict with: beta, selected, lambda_value, psi, n_iter, converged,
    method.
    """
    from sklearn.linear_model import Lasso

    n, p = W_raw.shape
    if alpha_p is None:
        alpha_p = 0.1 / max(np.log(max(n, 2)), 1.0)

    psi = np.ones(p)
    beta = np.zeros(p)
    selected = np.zeros(p, dtype=bool)
    n_iter = 0
    converged = False

    if p == 0:
        return dict(beta=beta, selected=selected, lambda_value=0.0,
                    psi=psi, n_iter=0, converged=True, method="trivial")

    # Weight scale and location do not matter (2026-10).
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.mean()
    y_raw = y_raw - np.average(y_raw, weights=weights)
    W_raw = W_raw - np.average(W_raw, axis=0, weights=weights)
    sqrt_w = np.sqrt(weights)

    # Handle keep_idx: partial out keep cols (weighted LS) before LASSO,
    # re-add to the selection. Centered inputs stay centered.
    if keep_idx is not None and len(keep_idx) > 0:
        keep_idx = np.asarray(keep_idx, dtype=int)
        nonkeep_mask = np.ones(p, dtype=bool)
        nonkeep_mask[keep_idx] = False
        nonkeep_idx = np.where(nonkeep_mask)[0]
        K_pm = W_raw[:, keep_idx] * sqrt_w[:, None]
        coef_y, *_ = np.linalg.lstsq(K_pm, y_raw * sqrt_w, rcond=None)
        y_lasso_raw = y_raw - W_raw[:, keep_idx] @ coef_y
        W_nonkeep = W_raw[:, nonkeep_idx]
        coef_W, *_ = np.linalg.lstsq(K_pm, W_nonkeep * sqrt_w[:, None], rcond=None)
        W_lasso_raw = W_nonkeep - W_raw[:, keep_idx] @ coef_W
        p_lasso = W_lasso_raw.shape[1]
    else:
        nonkeep_idx = np.arange(p)
        y_lasso_raw = y_raw
        W_lasso_raw = W_raw
        p_lasso = p
    y_lasso_pm = y_lasso_raw * sqrt_w
    W_lasso_pm = W_lasso_raw * sqrt_w[:, None]

    n_eff = _compute_n_eff(weights, cluster_ids, spatial_args=spatial_args)

    if p_lasso == 0:
        beta = np.zeros(p)
        selected = np.zeros(p, dtype=bool)
        if keep_idx is not None:
            selected[keep_idx] = True
        return dict(beta=beta, selected=selected, lambda_value=0.0, psi=np.ones(p),
                    n_iter=0, converged=True, method="keep-only")

    if lambda_method == "cv" or lambda_override is not None:
        # CV path: compute psi once with the initial residuals (the centered y),
        # then run a single LASSO at the chosen λ.
        u0 = y_lasso_raw
        psi_lasso = _compute_psi(W_lasso_raw, u0, weights, cluster_ids, spatial_args=spatial_args)
        col_scale = 1.0 / psi_lasso
        W_scaled_pm = W_lasso_pm * col_scale[None, :]

        if lambda_override is not None:
            lam = lambda_override
            method = "override"
        else:
            if cluster_ids is not None:
                from sklearn.model_selection import GroupKFold
                cv = GroupKFold(n_splits=min(cv_folds, len(np.unique(cluster_ids))))
                splits = list(cv.split(W_scaled_pm, y_lasso_pm, groups=cluster_ids))
            else:
                from sklearn.model_selection import KFold
                rng = np.random.default_rng(cv_seed)
                cv = KFold(n_splits=cv_folds, shuffle=True,
                           random_state=int(rng.integers(0, 2**31 - 1)))
                splits = list(cv.split(W_scaled_pm, y_lasso_pm))

            alpha_max = float(np.abs(W_scaled_pm.T @ y_lasso_pm).max() / n)
            alpha_grid = alpha_max * np.logspace(-3, 0, 50)
            cv_errs = np.zeros(len(alpha_grid))
            # Capture sklearn ConvergenceWarning across the full CV grid; we
            # only want to surface a single summary if convergence failed in
            # >some fraction of fits, rather than flooding stderr.
            import warnings
            from sklearn.exceptions import ConvergenceWarning
            n_conv_fails = 0
            n_total_fits = 0
            with warnings.catch_warnings():
                warnings.simplefilter("error", ConvergenceWarning)
                for tr, te in splits:
                    Xtr, Xte = W_scaled_pm[tr], W_scaled_pm[te]
                    ytr, yte = y_lasso_pm[tr], y_lasso_pm[te]
                    for ai, a in enumerate(alpha_grid):
                        n_total_fits += 1
                        mod = Lasso(alpha=a, fit_intercept=False,
                                    max_iter=20000, tol=1e-6)
                        try:
                            mod.fit(Xtr, ytr)
                        except ConvergenceWarning:
                            n_conv_fails += 1
                            # Re-fit silently after the warning so we get coefs
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", ConvergenceWarning)
                                mod.fit(Xtr, ytr)
                        pred = Xte @ mod.coef_
                        cv_errs[ai] += float(((yte - pred) ** 2).sum())
            if n_conv_fails > 0:
                pct = 100.0 * n_conv_fails / max(n_total_fits, 1)
                if pct >= 5.0:
                    print(f"  Note: {n_conv_fails}/{n_total_fits} ({pct:.0f}%) "
                          f"CV LASSO fits hit max_iter without converging. "
                          f"Selected λ may be unstable; consider lasso_penalty='plugin' instead.")
            best = int(np.argmin(cv_errs))
            alpha_best = float(alpha_grid[best])
            lam = 2.0 * n * alpha_best
            method = "cv"

        alpha_sklearn = lam / (2.0 * n)
        mod = Lasso(alpha=alpha_sklearn, fit_intercept=False, max_iter=20000, tol=1e-6)
        mod.fit(W_scaled_pm, y_lasso_pm)
        beta_unscaled_lasso = mod.coef_ * col_scale

        beta = np.zeros(p)
        beta[nonkeep_idx] = beta_unscaled_lasso
        psi_full = np.full(p, np.nan)
        psi_full[nonkeep_idx] = psi_lasso
        selected = np.zeros(p, dtype=bool)
        selected[nonkeep_idx] = beta_unscaled_lasso != 0
        if keep_idx is not None:
            selected[keep_idx] = True
        return dict(beta=beta, selected=selected, lambda_value=lam,
                    psi=psi_full, n_iter=1, converged=True, method=method)

    # ── Iterative plug-in path ────────────────────────────────────────────
    u = y_lasso_raw.copy()                 # initial residuals: centered y
    beta_lasso = np.zeros(p_lasso)

    for it in range(n_iter_max):
        n_iter = it + 1
        psi_lasso = _compute_psi(W_lasso_raw, u, weights, cluster_ids, spatial_args=spatial_args)
        lam = _bch_lambda(p_lasso, n_eff, c_slack, alpha_p)

        col_scale = 1.0 / psi_lasso
        W_scaled_pm = W_lasso_pm * col_scale[None, :]
        alpha_sklearn = lam / (2.0 * n)

        mod = Lasso(alpha=alpha_sklearn, fit_intercept=False, max_iter=20000, tol=1e-6)
        mod.fit(W_scaled_pm, y_lasso_pm)
        beta_lasso_new = mod.coef_ * col_scale

        # Post-LASSO residuals for the next loadings (hdm/lassopack default).
        u_new = _post_lasso_resid(y_lasso_raw, W_lasso_raw, beta_lasso_new != 0, sqrt_w)
        max_change = float(np.max(np.abs(beta_lasso_new - beta_lasso)))

        if verbose:
            n_active = int((beta_lasso_new != 0).sum())
            print(f"    iter {it+1}: λ={lam:.3f}, |Δβ|_inf={max_change:.2e}, "
                  f"active={n_active}/{p_lasso}")

        beta_lasso = beta_lasso_new
        u = u_new
        if max_change < tol:
            converged = True
            break

    beta = np.zeros(p)
    beta[nonkeep_idx] = beta_lasso
    psi_full = np.full(p, np.nan)
    psi_full[nonkeep_idx] = psi_lasso
    selected = np.zeros(p, dtype=bool)
    selected[nonkeep_idx] = beta_lasso != 0
    if keep_idx is not None:
        selected[keep_idx] = True

    return dict(beta=beta, selected=selected, lambda_value=lam,
                psi=psi_full, n_iter=n_iter, converged=converged, method="plugin")


# ──────────────────────────────────────────────────────────────────────────────
#  Quadruple-LASSO orchestrator
#
#  Runs four BCH plug-in (or CV) LASSOs — one each for Y, T, M, and Z —
#  on a shared rich-expansion design W. Returns the union of selected
#  variables along with per-LASSO diagnostics (selected names, coefficients,
#  λ values, convergence flags). Seeds are deterministic offsets of a single
#  user-supplied seed: Y → seed+0, T → seed+1, M → seed+2, Z → seed+3 (per-
#  instrument for K > 1: seed+3, seed+4, ...). The seed only matters for the
#  CV path; the plug-in path is fully deterministic.
#
#  IMPORTANT: M is regressed on W ALONE (not on W + T). See the paper,
#  Appendix H (Algorithm 1 and the remark after it): including T would
#  misallocate explanatory power away from controls that matter through their
#  correlation with T, biasing the post-selection 2SLS.
# ──────────────────────────────────────────────────────────────────────────────

def _quadruple_lasso(
    y_raw, t_raw, m_raw, z_mat_raw, W_raw, col_names,
    weights, cluster_ids, keep_idx,
    lambda_method="plugin", cv_folds=5, seed=20260317,
    n_iter_max=15, tol=1e-6, c_slack=1.1, alpha_p=None,
    lambda_override=None, spatial_args=None, verbose=False,
):
    """Run four LASSOs (Y, T, M, Z[k]) on the shared rich expansion W_raw.

    Parameters
    ----------
    y_raw, t_raw, m_raw : (n,) ndarrays, FE-demeaned, NOT pre-multiplied.
    z_mat_raw : (n, K) ndarray of instruments, FE-demeaned, NOT pre-multiplied.
    W_raw : (n, p) ndarray of rich-expansion features, FE-demeaned.
    col_names : list of length p (names matching W_raw columns).
    weights : (n,) positive weights (np.ones(n) if unweighted).
    cluster_ids : (n,) array or None.
    spatial_args : dict or None. Forwarded to _bch_lasso for spatial-aware
        penalty calibration. Takes precedence over cluster_ids if both given.
    keep_idx : list/array of column indices that should be unconditionally
        selected. They are partialed out before each LASSO and re-added
        to the union selection.

    Returns
    -------
    dict with keys:
        S         : union mask, shape (p,) bool
        S_Y, S_T, S_M, S_Z_k : per-LASSO selection masks (S_Z_k is a list,
                    one mask per instrument)
        details   : per-LASSO sub-dict with names, coefs, lambda, etc.
    """
    K_inst = z_mat_raw.shape[1] if z_mat_raw.ndim == 2 else 1
    z_mat_raw = z_mat_raw.reshape(-1, K_inst) if z_mat_raw.ndim == 1 else z_mat_raw

    p = W_raw.shape[1]
    if keep_idx is None:
        keep_idx = []

    common_kwargs = dict(
        weights=weights, cluster_ids=cluster_ids, keep_idx=list(keep_idx),
        lambda_method=lambda_method, cv_folds=cv_folds,
        n_iter_max=n_iter_max, tol=tol, c_slack=c_slack, alpha_p=alpha_p,
        lambda_override=lambda_override, spatial_args=spatial_args, verbose=verbose,
    )

    if verbose:
        print("  Running quadruple-LASSO:")
        print(f"    p (rich expansion) = {p}")
        print(f"    K (instruments)    = {K_inst}")
        print(f"    keep variables     = {len(keep_idx)}")

    # Each LASSO gets a deterministic seed offset
    if verbose:
        print("  → LASSO Y on W")
    res_Y = _bch_lasso(y_raw, W_raw, cv_seed=seed + 0, **common_kwargs)
    if verbose:
        print("  → LASSO T on W")
    res_T = _bch_lasso(t_raw, W_raw, cv_seed=seed + 1, **common_kwargs)
    if verbose:
        print("  → LASSO M on W (M ~ W only — see Remark on omitting T)")
    res_M = _bch_lasso(m_raw, W_raw, cv_seed=seed + 2, **common_kwargs)
    res_Z_per_inst = []
    for k in range(K_inst):
        if verbose:
            print(f"  → LASSO Z[{k}] on W")
        res_Zk = _bch_lasso(z_mat_raw[:, k], W_raw, cv_seed=seed + 3 + k,
                            **common_kwargs)
        res_Z_per_inst.append(res_Zk)

    # Build selection masks
    S_Y = res_Y['selected']
    S_T = res_T['selected']
    S_M = res_M['selected']
    # Union over instruments for Z
    S_Z = np.zeros(p, dtype=bool)
    for r in res_Z_per_inst:
        S_Z = S_Z | r['selected']

    # Final union (also includes keep_idx, since each individual LASSO already
    # has them in its selected mask, but we add them explicitly for safety)
    S = S_Y | S_T | S_M | S_Z
    if len(keep_idx) > 0:
        S[np.asarray(keep_idx, dtype=int)] = True

    def _summarize(res, label):
        sel = res['selected']
        idx = np.where(sel)[0]
        return dict(
            label=label,
            names=[col_names[j] for j in idx],
            coefs=res['beta'][sel].tolist(),
            indices=idx.tolist(),
            lambda_value=float(res['lambda_value']),
            n_selected=int(sel.sum()),
            converged=bool(res['converged']),
            n_iter=int(res['n_iter']),
            method=res['method'],
        )

    details = dict(
        Y=_summarize(res_Y, "Y"),
        T=_summarize(res_T, "T"),
        M=_summarize(res_M, "M"),
        Z=[_summarize(r, f"Z[{k}]") for k, r in enumerate(res_Z_per_inst)],
    )

    return dict(
        S=S, S_Y=S_Y, S_T=S_T, S_M=S_M, S_Z=S_Z,
        S_Z_per_inst=[r['selected'] for r in res_Z_per_inst],
        details=details,
        union_names=[col_names[j] for j in np.where(S)[0]],
        n_selected_union=int(S.sum()),
        n_total=int(p),
    )


def _build_rich_expansion(
    df, x_names, dummy_col_names,
    expand_x=("square", "cube"),
    expand_dummies=("interact_with_X",),
    expand_override=None,
    plain_names=False,
):
    """Build the rich expansion W from raw X and dummy data.

    Parameters
    ----------
    df : DataFrame restricted to the analysis sample
    x_names : list of column names for original (continuous or categorical) X
    dummy_col_names : list of column names for already-expanded dummy indicators
    expand_x : default ops applied to each variable in x_names
    expand_dummies : default ops applied to each variable in dummy_col_names
    expand_override : dict mapping variable name → list of ops (per-variable)
    plain_names : bool, if True use ASCII names ('age_sq' instead of 'age^2')

    Returns
    -------
    W : (n, p) ndarray of expansion columns (does NOT include the originals)
    col_names : list of length p
    meta : dict with diagnostic info on what was generated
    """
    if expand_override is None:
        expand_override = {}

    # Validate operation names
    for ops in [expand_x, expand_dummies] + list(expand_override.values()):
        for op in (ops or ()):
            if op not in _VALID_EXPAND_OPS:
                raise ValueError(
                    f"Unknown expansion op {op!r}. "
                    f"Valid: {sorted(_VALID_EXPAND_OPS)}"
                )

    n = len(df)
    # Track which raw variables are continuous
    cont_x = [v for v in x_names if _is_continuous(df[v])]

    # Naming convention
    if plain_names:
        SQ, CB, INT = "_sq", "_cb", "_x_"
        def name_pow(v, p): return f"{v}{SQ}" if p == 2 else f"{v}{CB}"
        def name_int(a, b): return f"{a}{INT}{b}"
    else:
        def name_pow(v, p): return f"{v}^{p}"
        def name_int(a, b): return f"{a} × {b}"

    # Resolve effective ops for each variable
    def _ops_for(var, default):
        if var in expand_override:
            return tuple(expand_override[var])
        return tuple(default)

    cols = []
    names = []
    meta = {
        "n_originals_x": len(x_names),
        "n_originals_dummies": len(dummy_col_names),
        "n_continuous_x": len(cont_x),
        "n_dummies": len(dummy_col_names),
        "n_squares": 0,
        "n_cubes": 0,
        "n_interactions_xx": 0,
        "n_interactions_dx": 0,
        "n_interactions_dd": 0,
        "warnings": [],
    }

    # ── Originals (X variables + dummy indicators) come first ──
    # These are part of the candidate set; LASSO can keep or drop them.
    # If the user wants any to be unconditionally included, pass them via
    # `lasso_keep` (resolved to indices in this combined column ordering).
    for v in x_names:
        cols.append(df[v].values.astype(float))
        names.append(v)
    for d in dummy_col_names:
        cols.append(df[d].values.astype(float))
        names.append(d)

    # ── Squares and cubes (continuous X only) ──
    for v in x_names:
        ops = _ops_for(v, expand_x)
        is_cont = v in cont_x
        if "square" in ops:
            if not is_cont:
                meta["warnings"].append(f"'square' has no effect on indicator-like variable {v!r}; skipped.")
            else:
                cols.append(df[v].values.astype(float) ** 2)
                names.append(name_pow(v, 2))
                meta["n_squares"] += 1
        if "cube" in ops:
            if not is_cont:
                meta["warnings"].append(f"'cube' has no effect on indicator-like variable {v!r}; skipped.")
            else:
                cols.append(df[v].values.astype(float) ** 3)
                names.append(name_pow(v, 3))
                meta["n_cubes"] += 1

    # Continuous-continuous interactions (pairs from x_names)
    for v in x_names:
        ops = _ops_for(v, expand_x)
        if "interact_with_X" not in ops:
            continue
        if v not in cont_x:
            # An indicator-like x with interact_with_X: we still allow it,
            # interpreting "interact_with_X" as "interact this variable with
            # all other variables in x_names". Skip the same-name case below.
            pass
        for u in x_names:
            if u == v:
                continue
            # Avoid duplicate (v, u) and (u, v) interactions: only add (v, u)
            # when v < u in the x_names ordering.
            if x_names.index(u) <= x_names.index(v):
                continue
            cols.append(df[v].values.astype(float) * df[u].values.astype(float))
            names.append(name_int(v, u))
            meta["n_interactions_xx"] += 1

    # Dummy × X interactions (dummy_col_names × x_names)
    for d in dummy_col_names:
        ops = _ops_for(d, expand_dummies)
        # Silently skip square/cube on dummies (mechanically pointless)
        if "square" in ops or "cube" in ops:
            meta["warnings"].append(
                f"'square'/'cube' on dummy variable {d!r} has no effect; skipped."
            )
        if "interact_with_X" in ops:
            for u in x_names:
                cols.append(df[d].values.astype(float) * df[u].values.astype(float))
                names.append(name_int(d, u))
                meta["n_interactions_dx"] += 1

    # Dummy × dummy interactions (dummy_col_names × dummy_col_names)
    # Skip (a, b) where a and b are indicators of the same source categorical
    # (their product is identically zero).
    for d in dummy_col_names:
        ops = _ops_for(d, expand_dummies)
        if "interact_with_dummies" not in ops:
            continue
        for e in dummy_col_names:
            if e == d:
                continue
            if dummy_col_names.index(e) <= dummy_col_names.index(d):
                continue
            # Heuristic: indicators of same source share a "_D_<source>_" prefix
            shared_prefix = None
            for src in set(d.rsplit("_", 1)[:1]) & set(e.rsplit("_", 1)[:1]):
                shared_prefix = src
            if d.startswith("_D_") and e.startswith("_D_"):
                # Reconstruct source name conservatively: source ends just
                # before the last underscore. This is heuristic; with the
                # naming "_D_<source>_<level>" we strip "_<level>" by finding
                # the last "_" position.
                src_d = d.rsplit("_", 1)[0]
                src_e = e.rsplit("_", 1)[0]
                if src_d == src_e:
                    continue  # same source → product is 0
            cols.append(df[d].values.astype(float) * df[e].values.astype(float))
            names.append(name_int(d, e))
            meta["n_interactions_dd"] += 1

    # Build the matrix. If empty, return a (n, 0) array.
    if not cols:
        W = np.zeros((n, 0))
    else:
        W = np.column_stack(cols).astype(float)

    # Detect duplicate column names. This catches the case where a user passes
    # X = ["age", "age^2"] (already squared in their data) AND requests
    # expand_x=("square",), which would generate two columns named "age^2".
    # Since the LASSO selection identifies columns by name, duplicates create
    # ambiguity. Raise rather than silently keeping or dropping.
    if len(names) != len(set(names)):
        from collections import Counter
        counts = Counter(names)
        dups = sorted(k for k, v in counts.items() if v > 1)
        raise ValueError(
            "Rich expansion produced duplicate column names: "
            f"{dups!r}. This usually means a user-specified variable was "
            "already an expansion-style term (e.g., 'age^2') and the same "
            "expansion would generate it again. Either remove the duplicate "
            "from X, exclude the operation via expand_override, or rename the "
            "user-supplied column."
        )

    # Drop near-constant columns (variance below tolerance)
    if W.shape[1] > 0:
        col_var = W.var(axis=0)
        col_mean_abs = np.maximum(np.abs(W.mean(axis=0)), 1.0)
        rel_tol = 1e-12
        keep_mask = col_var > (rel_tol * col_mean_abs ** 2)
        n_dropped = (~keep_mask).sum()
        if n_dropped > 0:
            meta["warnings"].append(f"Dropped {n_dropped} near-constant columns from rich expansion.")
        W = W[:, keep_mask]
        names = [names[j] for j in range(len(names)) if keep_mask[j]]

    # Emit warnings (only first few to avoid flooding the console)
    for msg in meta["warnings"][:3]:
        print(f"  Note: {msg}")
    if len(meta["warnings"]) > 3:
        print(f"  Note: ... (+{len(meta['warnings']) - 3} more expansion warnings)")

    meta["p"] = W.shape[1]
    return W, names, meta


def ivmediate(data, Y, T, M, Z, X=None, dummies=None, fe=None,
                          se="robust", clusterid=None, alpha=0.05,
                          weights=None, wtype="aweight",
                          op_tau=0.10, op_alpha=0.05,
                          # ── Anderson–Rubin confidence sets ──
                          ar_ci=True,
                          ar_grid_halfwidth=200.0,
                          ar_grid_points=4001,
                          # ── Spatial (Conley) standard errors ──
                          lat=None, lon=None,
                          spatial_cutoff=None,
                          spatial_kernel=None,
                          spatial_metric=None,
                          # ── Post-double-selection LASSO ──
                          lasso=False,
                          lasso_penalty="plugin",
                          lasso_keep=None,
                          lasso_seed=20260317,
                          lasso_cv_folds=5,
                          expand_x=("square", "cube"),
                          expand_dummies=("interact_with_X",),
                          expand_override=None,
                          show_selection=False,
                          lasso_c_slack=1.1,
                          lasso_alpha_p=None,
                          lasso_n_iter_max=15,
                          lasso_tol=1e-6,
                          lasso_lambda_override=None,
                          lasso_verbose=False,
                          _point_only=False):
    """
    Run the full IV mediation analysis with diagnostics.

    Parameters
    ----------
    data : pd.DataFrame
        The dataset (never modified).
    Y : str
        Name of the outcome variable.
    T : str
        Name of the treatment variable.
    M : str
        Name of the mediator variable.
    Z : str, list of str or tuple of str
        Name(s) of instrument variable(s), typically one. They must differ
        from Y, T and M, and may not repeat.
    X : list of str, str, or None
        Control variables (covariates). Accepts a list of column names or a
        single whitespace-separated string; list elements and string tokens
        may be shell-style globs (e.g., "reg_*") resolved against
        data.columns. Multiline strings are fine. Examples:
            X = ["age", "income", "reg_a", "reg_b"]            # explicit list
            X = "age income reg_*"                              # one-liner with glob
            X = '''age income
                   reg_a reg_b reg_c'''                        # multiline paste
        Globs never pick up Y, T, M, Z or the fe/cluster/weight/coordinate
        columns (a note lists the excluded matches); naming Y, T, M or Z
        explicitly is an error. Exactly collinear controls (after absorbing
        the FE) are dropped with a note "dropped (collinear): …", keeping the
        earlier column of a collinear set; df count only the kept columns.
        Default: no controls.
    dummies : list of str, str, or None
        Categorical variables to expand into indicator dummies (drop-first).
        Same parsing rules and role protection as X. Empty or whitespace-only
        text values count as missing. For variables with few levels (< ~50).
        Default: none.
    fe : str or list of str or None
        Name(s) of fixed-effect variable(s), absorbed by within-group
        demeaning (FWL); one-way or two-way (the exact projection when the
        FE dimension with fewer levels has at most 3,000 levels, otherwise
        alternating projections with Irons–Tuck acceleration, relative
        tolerance 1e-12). With weights the
        group means are weighted. Singleton groups are dropped iteratively.
        Absorbed df = rank of the FE dummies beyond the intercept.
        Empty or whitespace-only text values count as missing. Default: none.
    se : str
        Standard error type: 'classical', 'robust' (HC1), 'cluster' (CR1) or
        'spatial' (Conley). Default: 'robust'.
    clusterid : str or None
        Name of the cluster variable. Required if se='cluster' and used only
        then (otherwise a note says it is ignored). With FE nested in the
        clusters, those FE are not counted in the CR1 small-sample factor.
    alpha : float in (0, 1)
        Significance level for CIs, AR sets and the RV/XRV threshold.
        Default: 0.05.
    weights : str or None
        Name of the weight column. Negative weights are an error; zero-weight
        rows are dropped (with a note) before singleton detection and FE
        demeaning. Default: None (unweighted).
    wtype : str
        Weight type: 'aweight' (analytical) or 'pweight' (probability).
        'pweight' with se='classical' is upgraded to 'robust'. Default: 'aweight'.
    lat, lon : str or None
        Names of latitude and longitude columns (in degrees, by default).
        Required when se='spatial' and an error otherwise. Rows with missing
        coordinates are dropped. The columns may also serve as projected
        (Cartesian) coordinates if spatial_metric='euclidean'.
    spatial_cutoff : float or None
        Distance cutoff defining the kernel support for spatial vcov.
        Required when se='spatial' (an error otherwise). Units are kilometers
        under spatial_metric='haversine' (the default for lat/lon in degrees);
        otherwise the units of lat/lon. Conley SE results are sensitive to
        this choice — applied practice is to report results at multiple
        cutoffs (e.g., 100, 250, 500 km) as a robustness check. A uniform-
        kernel cutoff that covers every pair is an error (the Conley meat is
        0 by construction); a neighbor share above 50% prints a warning.
    spatial_kernel : {"uniform", "bartlett"} or None
        Kernel function K(d) for spatial weighting. "uniform" gives equal
        weight 1 to all neighbors within cutoff (hard cutoff); "bartlett"
        gives a linear taper K(d) = max(1 - d/cutoff, 0). None means
        "uniform" with se='spatial'; giving it with another se is an error.
    spatial_metric : {"haversine", "euclidean"} or None
        Distance metric. "haversine" interprets lat/lon as degrees and
        computes great-circle distance (in km, given Earth radius 6371 km);
        "euclidean" computes raw Euclidean distance in the units of lat/lon
        (use this if your coordinates are projected). None means "haversine"
        with se='spatial'; giving it with another se is an error.
    op_tau : float in (0, 1)
        Worst-case 2SLS bias tolerance, as a fraction of the OLS bias
        benchmark, for the Olea–Pflueger (2013) critical value of the
        first-stage F-statistics. Default: 0.10 (OP's recommendation).
    op_alpha : float in (0, 1)
        Size of the OP weak-instrument test. Default: 0.05. Critical values
        come from the simplified test at the estimated effective degrees of
        freedom K_eff (OP eq. 4; K_eff = 1 for one instrument), computed for
        any (op_tau, op_alpha). Distinct from `alpha` (CI level).
    ar_ci : bool
        If True (default), invert the Anderson–Rubin test to produce
        weak-instrument-robust confidence sets for θ (mediator effect,
        instrumenting M with Z conditional on T), γ (treatment → mediator),
        and τ_total. The inversion uses the same SE regime and `alpha` as
        the analytic inference. With one instrument the set is computed in
        closed form (Fieller inversion; shapes: bounded interval, union of
        two rays, or the whole line); with K > 1 instruments the joint
        chi2(K) AR statistic is computed in closed form, located on a search
        grid and its boundaries refined by Brent's method. With spatial SEs,
        a set whose variance estimate turns negative within the search range
        (the kernel matrix is not PSD) is reported as unavailable. No AR set
        is reported for the direct effect τ: with M endogenous, τ is the
        coefficient on an included regressor and requires subset-AR
        methods. Set to False to skip.
    ar_grid_halfwidth : float
        Half-width of the AR search grid in units of the parameter's
        analytic SE (K > 1, and the K = 1 fallback; also the range checked
        for negative spatial variances). Default: 200.
    ar_grid_points : int
        Number of AR search-grid points (K > 1, and the K = 1 fallback).
        Default: 4001.
    lasso : bool
        If True, run post-double-selection LASSO (the "quadruple LASSO"
        procedure, paper App. H) before the structural 2SLS. Four LASSO
        regressions are run on a rich expansion of (X, dummies): one each
        for Y, T, M, and each instrument in Z. The union of selected
        variables forms the post-selection control set. Each LASSO has an
        unpenalized intercept and uses weights normalized to mean 1.
        Default: False.
    lasso_penalty : {"plugin", "cv"}
        Method for choosing the LASSO penalty λ. "plugin" uses the rigorous
        iterative BCH plug-in formula with post-LASSO residual loadings
        (recommended; deterministic given the data). "cv" uses K-fold
        cross-validation (folds drawn with `lasso_seed`; with se='cluster'
        the folds are whole clusters). Any other value is an error.
        Default: "plugin".
    lasso_keep : list of str, str, or None
        Variable names that should always appear in the post-selection
        control set, regardless of LASSO selection. These are partialed out
        before each LASSO and re-added unconditionally. Names are resolved
        against the *rich expansion*, so they can refer to:
            - originals (e.g., "X1", "region")
            - powers (e.g., "X1^2", "income^3")
            - interactions (e.g., "X1 × X2", "_D_region_2 × age")
        Globs (* ? [...]) are supported and expand against rich-expansion
        names: e.g., "X1*" matches X1, X1^2, X1^3. Note: interaction names
        contain spaces around the × symbol, so to force an interaction use
        the list form (e.g., lasso_keep=["X1 × X2"]) — the string form would
        split on the embedded whitespace.
        Default: no forced inclusions.
    lasso_seed : int
        Random seed for the CV path (irrelevant for "plugin"). The four LASSOs
        get deterministic offsets: seed+0 (Y), seed+1 (T), seed+2 (M),
        seed+3 (first instrument), etc. Default: 20260317.
    lasso_cv_folds : int
        Number of folds for cross-validation (>= 2). Default: 5.
    expand_x : tuple of str
        Default expansion ops for the variables in X: "square" and "cube"
        add powers of variables with more than two distinct values (others
        are skipped with a note), and "interact_with_X" adds the products of
        pairs of X variables. "interact_with_dummies" is accepted but has no
        effect on X. Default: ("square", "cube").
    expand_dummies : tuple of str
        Default expansion ops for the indicators made from `dummies`:
        "interact_with_X" adds their products with every X variable, and
        "interact_with_dummies" the products of indicators of different
        variables. Powers have no effect on indicators (skipped with a
        note). Default: ("interact_with_X",).
    expand_override : dict or None
        Per-variable override of expand_x / expand_dummies. Keys are names
        in X or indicator names such as "_D_region_2" (other keys are
        ignored); values are tuples of ops. Example:
        {"age": ("square",)} (age gets only square, no cube). Default: None.
    show_selection : bool
        If True, print_results prints the LASSO selection details (sets,
        penalties, convergence). Default: False (the selections are in the
        result dict either way).
    lasso_c_slack : float > 0
        Slack constant c in the BCH penalty formula λ = 2c√n_eff Φ⁻¹(...),
        with n_eff = n (Kish n with weights) in every SE regime.
        Default: 1.1 (BCH's recommendation).
    lasso_alpha_p : float in (0, 1) or None
        Tail probability in the Bonferroni-style penalty quantile. If None,
        defaults to 0.1/log(n) (BCH's data-dependent default), with n the
        number of observations also under clustering (lassopack's cluster
        default uses 0.1/log(G) instead). Default: None.
    lasso_n_iter_max : int
        Maximum iterations for the plug-in loadings. Default: 15.
    lasso_tol : float
        Convergence tolerance on |Δβ|_∞ across plug-in iterations. Default: 1e-6.
    lasso_lambda_override : float or None
        Optional fixed λ value, bypasses plug-in/CV. For replicability
        studies that need to match a specific published λ. Default: None.
    lasso_verbose : bool
        Print per-iteration progress for each LASSO. Default: False.

    Booleans must be True or False (a string such as "false" is an error).

    Returns
    -------
    dict with all estimation results, diagnostics, and sensitivity measures:
    the five estimates with SEs and Wald CIs (theta_hat, tau_hat, gamma_hat,
    total_hat, indirect; normal critical values), pi1/se_pi1 (first stage of T
    on Z, lists for K ≥ 2), n, n_clusters, df_fe_absorbed,
    fe_nested_in_cluster, dropped_collinear, the F statistics with K_eff_TZ,
    K_eff_MZ_T and the OP critical values, the AR sets, the κ block
    (c_TM_hat, kappa_star, kappa_star_CI, kappa_star_CI_informative,
    rho_TY_star, kappa_adj_table, F_MZ_T_above_crit), the R² block
    (F_reduced_form, F_first_stage_cond — both t², with deprecated aliases
    f_reduced_form, f_first_stage_cond — RV_theta, XRV_theta and components),
    the LASSO output and the residuals eps_T_hat, eps_M_hat, eps_Y_hat.
    """

    # ── Validate options (2026-10: enumerations, ranges and booleans) ──────
    _check_choice("se", se, _SE_TYPES)
    _check_choice("wtype", wtype, ("aweight", "pweight"))
    _check_open_unit("alpha", alpha)
    _check_open_unit("op_tau", op_tau)
    _check_open_unit("op_alpha", op_alpha)
    for _nm, _v in (("ar_ci", ar_ci), ("lasso", lasso), ("show_selection", show_selection),
                    ("lasso_verbose", lasso_verbose), ("_point_only", _point_only)):
        _check_bool(_nm, _v)
    _check_choice("lasso_penalty", lasso_penalty, ("plugin", "cv"))
    try:
        _hw, _gp = float(ar_grid_halfwidth), float(ar_grid_points)
    except (TypeError, ValueError):
        raise ValueError("ar_grid_halfwidth and ar_grid_points must be numbers.")
    if not (np.isfinite(_hw) and _hw > 0):
        raise ValueError(f"ar_grid_halfwidth must be positive; got {ar_grid_halfwidth!r}.")
    if not (_gp.is_integer() and _gp >= 3):
        raise ValueError(f"ar_grid_points must be an integer >= 3; got {ar_grid_points!r}.")
    if lasso:
        if not (np.isfinite(float(lasso_c_slack)) and float(lasso_c_slack) > 0):
            raise ValueError(f"lasso_c_slack must be positive; got {lasso_c_slack!r}.")
        if lasso_alpha_p is not None:
            _check_open_unit("lasso_alpha_p", lasso_alpha_p)
        if not (float(lasso_n_iter_max).is_integer() and lasso_n_iter_max >= 1):
            raise ValueError(f"lasso_n_iter_max must be a positive integer; got {lasso_n_iter_max!r}.")
        if not (float(lasso_tol) > 0):
            raise ValueError(f"lasso_tol must be positive; got {lasso_tol!r}.")
        if not (float(lasso_cv_folds).is_integer() and lasso_cv_folds >= 2):
            raise ValueError(f"lasso_cv_folds must be an integer >= 2; got {lasso_cv_folds!r}.")
        if lasso_lambda_override is not None and not (float(lasso_lambda_override) >= 0):
            raise ValueError(f"lasso_lambda_override must be >= 0; got {lasso_lambda_override!r}.")

    # Spatial arguments are used only with se='spatial' (2026-10: an error otherwise).
    spatial_given = [nm for nm, v in (("lat", lat), ("lon", lon), ("spatial_cutoff", spatial_cutoff),
                                      ("spatial_kernel", spatial_kernel),
                                      ("spatial_metric", spatial_metric)) if v is not None]
    if se == "spatial":
        if lat is None or lon is None:
            raise ValueError("se='spatial' requires lat and lon column names.")
        if spatial_cutoff is None:
            raise ValueError(
                "se='spatial' requires spatial_cutoff (distance bandwidth). "
                "Conley (1999) results are sensitive to this choice; report "
                "results at multiple cutoffs as a robustness check."
            )
        try:
            _cut = float(spatial_cutoff)
        except (TypeError, ValueError):
            raise ValueError(f"spatial_cutoff must be a positive number; got {spatial_cutoff!r}.")
        if not (np.isfinite(_cut) and _cut > 0):
            raise ValueError(f"spatial_cutoff must be positive; got {spatial_cutoff}.")
        spatial_kernel = "uniform" if spatial_kernel is None else spatial_kernel
        spatial_metric = "haversine" if spatial_metric is None else spatial_metric
        _check_choice("spatial_kernel", spatial_kernel, ("uniform", "bartlett"))
        _check_choice("spatial_metric", spatial_metric, ("haversine", "euclidean"))
    elif spatial_given:
        raise ValueError(
            f"{', '.join(spatial_given)} given but se={se!r}: the spatial arguments "
            f"are used only with se='spatial'."
        )

    # Cluster ids are used only with cluster SEs (D3).
    if se == "cluster" and clusterid is None:
        raise ValueError("se='cluster' requires clusterid (the cluster variable).")
    if clusterid is not None and se != "cluster":
        print("  Note: clusterid ignored (se != 'cluster').")

    # ── Weight-type / SE-regime consistency ──────────────────────────────
    # Under pweights (w_i = 1/π_i, survey sampling), classical homoskedastic SEs
    # are invalid: the weighting introduces mechanical heteroskedasticity that is
    # not a nuisance. Force robust. Under aweights (cell-mean weights), classical
    # SEs are valid if the aweight specification is correctly calibrated to the
    # underlying variance structure; robust/cluster SEs are a defensive default.
    if weights is not None and wtype == "pweight" and se == "classical":
        print("  Note: pweight + classical SE is inconsistent. "
              "Upgrading SE to 'robust' (HC1 sandwich).")
        se = "robust"
    elif weights is not None and wtype == "aweight" and se == "classical":
        print("  Note: aweight + classical SE assumes the aweight specification "
              "is correctly calibrated. Consider 'robust' or 'cluster' as a "
              "defensive default.")
    se_v = "none" if _point_only else se           # variance regime actually computed
    if _point_only:
        ar_ci = False

    # ── Roles, missing values, listwise deletion, zero weights ───────────
    # X and dummies accept a list/tuple or a whitespace-separated string with
    # optional shell-style globs (resolved against data.columns); globs never
    # pick up Y, T, M, Z or the fe/cluster/weight/coordinate columns, and an
    # explicit Y, T, M or Z among the controls is an error. Text id columns
    # (fe, cluster, dummies) treat '' and whitespace-only values as missing;
    # ±Inf in a numeric role is an error. Zero-weight rows are dropped here,
    # before singleton detection and FE demeaning.
    prep = _prepare_sample(data, Y, T, M, Z, X, dummies, fe, se, clusterid, weights,
                           lat, lon, lasso=lasso)
    Z, X, dummies = prep["Z"], prep["X"], prep["dummies"]
    fe_cols, clusterid, df = prep["fe_cols"], prep["clusterid"], prep["df"]
    note_X = _summarize_resolved("X", X)
    note_D = _summarize_resolved("dummies", dummies)
    if note_X is not None:
        print(note_X)
    if note_D is not None:
        print(note_D)
    for line in prep["notes"]:
        print(line)
    is_weighted = weights is not None
    if prep["n_zero_weights"]:
        print(f"  Note: dropping {prep['n_zero_weights']} observation(s) with zero weight.")

    # Expand dummy variables into indicator columns (on the estimation sample)
    dummy_col_names = []
    if len(dummies) > 0:
        dummy_col_names = _expand_dummies(df, dummies)
        if not lasso:
            # Standard path: dummy indicators get appended to X and treated as
            # ordinary controls. Under LASSO, they are kept separate and fed
            # to _build_rich_expansion (which knows that indicators don't get
            # squared/cubed and applies dummy-specific interaction rules).
            X = X + dummy_col_names

    # Drop singleton FE groups (iteratively). With two-way FE, dropping a
    # singleton in one dimension can create a new singleton in the other, so
    # we alternate over the dimensions until a full pass drops nothing
    # (reghdfe convention, Correia 2016). Group counts are taken AFTER
    # convergence so the FE degrees-of-freedom adjustment reflects the final
    # estimation sample.
    n_dropped_singletons = 0
    n_fe_groups = 0
    if fe_cols:
        while True:
            dropped_this_pass = 0
            for fe_col in fe_cols:
                df, _, n_drop = _drop_singletons(df, fe_col)
                dropped_this_pass += n_drop
            n_dropped_singletons += dropped_this_pass
            if dropped_this_pass == 0:
                break
        if n_dropped_singletons > 0:
            print(f"  Note: dropped {n_dropped_singletons} singleton observations "
                  f"across {len(fe_cols)} FE dimension(s).")

    n = len(df)
    if n == 0:
        raise ValueError("No observations remaining after dropping missing values and singletons.")
    n_clusters = int(df[clusterid].nunique()) if clusterid else None
    if clusterid and n_clusters < 2:
        raise ValueError("cluster SEs need at least 2 clusters.")

    # Absorbed FE degrees of freedom (rank of the FE dummies beyond the
    # intercept: G−1 one-way, G₁+G₂−c−1 two-way with c connected components).
    # With cluster SEs, FE dimensions nested in the clusters do not count
    # toward the CR1 correction (reghdfe/fixest convention).
    df_fe_absorbed = 0
    df_vcov = 0
    fe_nested_in_cluster = None
    if fe_cols:
        df_fe_absorbed, df_fe_cr1, n_fe_groups, nested = _fe_degrees_of_freedom(
            df, fe_cols, clusterid if se == "cluster" else None)
        if se == "cluster":
            fe_nested_in_cluster = nested
            df_vcov = df_fe_cr1
        else:
            df_vcov = df_fe_absorbed

    # ── LASSO: build rich expansion BEFORE FE demeaning ──────────────────
    # The expansion uses level-form X and dummy indicators, then is demeaned
    # alongside the structural variables. This preserves the level-based
    # interpretation of squares/cubes/interactions while letting downstream
    # FWL handle the FE absorption uniformly.
    rich_col_names = []
    rich_meta = None
    if lasso:
        W_rich_levels, rich_col_names, rich_meta = _build_rich_expansion(
            df, x_names=X, dummy_col_names=dummy_col_names,
            expand_x=expand_x, expand_dummies=expand_dummies,
            expand_override=expand_override,
        )
        # Store rich expansion columns in df with the names we'll use
        # downstream, in a single pd.concat (no per-column inserts). A name that
        # clashes with an existing column is stored as _RX_{j}_{name}.
        rich_in_df = []
        existing_cols = set(df.columns)
        for j, nm in enumerate(rich_col_names):
            rich_in_df.append(f"_RX_{j}_{nm}" if nm in existing_cols else nm)
        if rich_in_df:
            _saved_index = df.index
            df = df.reset_index(drop=True)
            df = pd.concat(
                [df, pd.DataFrame(W_rich_levels, columns=rich_in_df)], axis=1
            )
            df.index = _saved_index
        # Resolve lasso_keep AGAINST THE RICH-EXPANSION NAMES, not data.columns.
        # This lets users force any candidate term — original ("X1") OR derived
        # ("X1^2", "X1 × X2") — into the post-selection set. Globs expand against
        # rich_col_names too, e.g., lasso_keep = "f_pruss*" would force every
        # term involving f_pruss (the variable, its powers, and its interactions).
        if lasso_keep is None:
            lasso_keep_resolved = []
        elif isinstance(lasso_keep, str) and not lasso_keep.strip():
            lasso_keep_resolved = []
        else:
            if isinstance(lasso_keep, str):
                tokens_keep = re.split(r"\s+", lasso_keep.strip())
            else:
                tokens_keep = list(lasso_keep)
            lasso_keep_resolved = []
            seen_keep = set()
            for tok in tokens_keep:
                if any(c in tok for c in _GLOB_CHARS):
                    matches = fnmatch.filter(rich_col_names, tok)
                    if not matches:
                        raise ValueError(
                            f"lasso_keep: glob pattern {tok!r} matched 0 names "
                            f"in the rich expansion. Available names include: "
                            f"{', '.join(rich_col_names[:8])}"
                            f"{'...' if len(rich_col_names) > 8 else ''}"
                        )
                    for m_ in matches:
                        if m_ not in seen_keep:
                            seen_keep.add(m_); lasso_keep_resolved.append(m_)
                else:
                    if tok not in rich_col_names:
                        raise ValueError(
                            f"lasso_keep variable {tok!r} is not in the rich "
                            f"expansion. lasso_keep accepts originals (e.g., 'X1'), "
                            f"powers (e.g., 'X1^2'), and interactions (e.g., 'X1 × X2'). "
                            f"Available names include: "
                            f"{', '.join(rich_col_names[:8])}"
                            f"{'...' if len(rich_col_names) > 8 else ''}"
                        )
                    if tok not in seen_keep:
                        seen_keep.add(tok); lasso_keep_resolved.append(tok)
        keep_idx_lasso = [rich_col_names.index(nm) for nm in lasso_keep_resolved]
        if rich_meta and rich_meta.get('p', 0) > 0:
            print(f"  Rich expansion built: p = {rich_meta['p']} "
                  f"(originals_x={rich_meta['n_originals_x']}, "
                  f"originals_dummies={rich_meta['n_originals_dummies']}, "
                  f"sq={rich_meta['n_squares']}, cb={rich_meta['n_cubes']}, "
                  f"x×x={rich_meta['n_interactions_xx']}, "
                  f"d×x={rich_meta['n_interactions_dx']}, "
                  f"d×d={rich_meta['n_interactions_dd']})")
        if lasso_keep_resolved:
            print(f"  LASSO keep ({len(lasso_keep_resolved)}): "
                  f"{', '.join(lasso_keep_resolved[:8])}"
                  f"{'...' if len(lasso_keep_resolved) > 8 else ''}")
        control_cols = rich_in_df
    else:
        rich_in_df = []
        keep_idx_lasso = []
        control_cols = list(X)

    # Weights (all positive after the zero-weight drop above)
    w_vec = df[weights].to_numpy(dtype=float) if is_weighted else np.ones(n)
    sqrt_w = np.sqrt(w_vec)

    # Reference scales for the collinearity and identification checks: the
    # weighted variation of each column around its weighted mean BEFORE the
    # FE are absorbed.
    def _ref_norm(v):
        v = np.asarray(v, dtype=float)
        return float(np.sqrt(np.sum(w_vec * (v - np.average(v, weights=w_vec)) ** 2)))
    ref_norm_ctrl = {c: _ref_norm(df[c].to_numpy(dtype=float)) for c in control_cols}
    ref_norm_role = {c: _ref_norm(df[c].to_numpy(dtype=float)) for c in [T, M] + Z}

    # Demean within FE groups (weighted group means under weights; FWL
    # projection onto the fixed effects). All analysis variables are demeaned
    # on the same estimation sample.
    analysis_cols = [Y, T, M] + Z + control_cols
    if fe_cols:
        wts = w_vec if is_weighted else None
        if len(fe_cols) == 1:
            df = _demean_within_groups(df, analysis_cols, fe_cols[0], weights=wts)
        else:
            df = _demean_twoway(df, analysis_cols, fe_cols, weights=wts)

    # After demeaning, the intercept is absorbed into the FE. An intercept is
    # still included in every regression (it is ≈ 0 with FE).
    y = df[Y].to_numpy(dtype=float)
    t = df[T].to_numpy(dtype=float)
    m = df[M].to_numpy(dtype=float)
    z_mat = df[Z].to_numpy(dtype=float).reshape(n, -1)
    x_mat = df[control_cols].to_numpy(dtype=float).reshape(n, -1) if control_cols else None
    cluster_ids = df[clusterid].to_numpy() if clusterid else None

    # Coordinates for spatial SEs (rows with missing coordinates were dropped
    # by the listwise deletion above).
    if se == "spatial":
        lat_arr = df[lat].to_numpy(dtype=float)
        lon_arr = df[lon].to_numpy(dtype=float)
        if spatial_metric == "haversine":
            if np.any(np.abs(lat_arr) > 90.0) or np.any(np.abs(lon_arr) > 180.0):
                raise ValueError(
                    "Lat/lon values are outside the (-90, 90) × (-180, 180) "
                    "degree range. Check that you have lat/lon in degrees, or "
                    "switch to spatial_metric='euclidean' for projected "
                    "coordinates."
                )
    else:
        lat_arr = None
        lon_arr = None

    # Centralized dict of SE-related kwargs threaded into every call to
    # _2sls_with_vcov, so that the structural coefficients, first-stage F's,
    # RV f-statistics and AR statistics use the SAME variance estimator.
    se_kwargs = dict(
        se=se_v,
        cluster_ids=cluster_ids,
        lat=lat_arr,
        lon=lon_arr,
        spatial_cutoff=spatial_cutoff,
        spatial_kernel=spatial_kernel,
        spatial_metric=spatial_metric,
    )

    # ── Vcov cache: quantities that do not change across the regressions ──
    # Cluster: factorized ids. Spatial: neighbor lists and the kernel matrix
    # (built once; reused by the LASSO loadings).
    vcov_cache = {"neg_counter": {}}
    sp_indices = sp_distances = None
    if se == "spatial" and (not _point_only or lasso):
        sp_indices, sp_distances = _build_spatial_neighbors(lat_arr, lon_arr, spatial_cutoff,
                                                            spatial_metric)
        nnz = int(sum(len(ix) for ix in sp_indices))
        share = nnz / float(n) / float(n)
        if spatial_kernel == "uniform" and nnz >= n * n:
            raise ValueError("the spatial cutoff covers all pairs; Conley SEs are degenerate "
                             "(0 by construction). Choose a smaller spatial_cutoff.")
        if share > 0.5:
            print(f"  Warning: the spatial cutoff puts {100 * share:.0f}% of all pairs within "
                  f"the cutoff (average neighbor share > 50%); Conley SEs are unreliable "
                  f"in this regime. Consider a smaller cutoff.")
    if se_v == "cluster":
        ci_, cu_ = pd.factorize(cluster_ids, sort=False)
        vcov_cache['cluster_unique'] = cu_
        vcov_cache['cluster_inverse'] = ci_
    elif se_v == "spatial":
        vcov_cache['spatial_indices'] = sp_indices
        vcov_cache['spatial_distances'] = sp_distances
        vcov_cache['spatial_adjacency'] = _build_spatial_sparse_adjacency(
            sp_indices, sp_distances, spatial_kernel, spatial_cutoff, n
        )
    se_kwargs['vcov_cache'] = vcov_cache

    # Raw (un-premultiplied) arrays are retained alongside the √w-transformed
    # versions. The √w-transformed arrays flow through 2SLS for point estimates
    # and the sandwich vcov; the raw arrays + w_vec feed the sensitivity
    # machinery and the LASSO solver, which need Hajek-form weighted moments.
    y_raw = y.copy()
    t_raw = t.copy()
    m_raw = m.copy()
    z_mat_raw = z_mat.copy()
    x_mat_raw = x_mat.copy() if x_mat is not None else None

    # ── LASSO step (run on raw-scale, FE-demeaned data) ────────────────────
    lasso_info = None
    if lasso:
        if x_mat_raw is None or x_mat_raw.shape[1] == 0:
            # Nothing to select from. Treat as no controls.
            lasso_info = dict(
                S=np.zeros(0, dtype=bool),
                S_Y=np.zeros(0, dtype=bool),
                S_T=np.zeros(0, dtype=bool),
                S_M=np.zeros(0, dtype=bool),
                S_Z=np.zeros(0, dtype=bool),
                S_Z_per_inst=[],
                details=dict(Y=None, T=None, M=None, Z=[]),
                union_names=[], n_selected_union=0, n_total=0,
            )
            x_mat_raw = None
            X = []
        else:
            if se == "spatial":
                # The LASSO loadings must reflect the spatial dependence; reuse
                # the neighbor lists built above.
                lasso_spatial_args = dict(
                    indices=sp_indices,
                    distances=sp_distances,
                    kernel=spatial_kernel,
                    cutoff=spatial_cutoff,
                )
            else:
                lasso_spatial_args = None
            lasso_info = _quadruple_lasso(
                y_raw=y_raw, t_raw=t_raw, m_raw=m_raw, z_mat_raw=z_mat_raw,
                W_raw=x_mat_raw, col_names=rich_col_names,
                weights=w_vec, cluster_ids=cluster_ids,
                keep_idx=keep_idx_lasso,
                lambda_method=lasso_penalty,
                cv_folds=lasso_cv_folds,
                seed=lasso_seed,
                n_iter_max=lasso_n_iter_max,
                tol=lasso_tol,
                c_slack=lasso_c_slack,
                alpha_p=lasso_alpha_p,
                lambda_override=lasso_lambda_override,
                spatial_args=lasso_spatial_args,
                verbose=lasso_verbose,
            )
            sel_idx = np.where(lasso_info['S'])[0]
            if len(sel_idx) == 0:
                x_mat_raw = None
                X = []
            else:
                x_mat_raw = x_mat_raw[:, sel_idx]
                X = [rich_col_names[j] for j in sel_idx]
                control_cols = [rich_in_df[j] for j in sel_idx]
            print(f"  LASSO selected {lasso_info['n_selected_union']} of "
                  f"{lasso_info['n_total']} candidate variables "
                  f"({lasso_penalty} penalty).")

    # ── Collinear controls (D9): drop exactly redundant columns of [1, X] ──
    # Columns are tested in order on the weighted, FE-demeaned design; a column
    # whose residual on [1, earlier columns] is below 1e-9 of its variation
    # before the FE were absorbed is dropped with a note. df count only the
    # kept columns.
    dropped_collinear = []
    if x_mat_raw is not None and x_mat_raw.shape[1] > 0:
        A = x_mat_raw * sqrt_w[:, None]
        refs = np.array([ref_norm_ctrl[c] for c in control_cols])
        kept, dropped_collinear = _greedy_independent(
            A, list(X), keep_first=sqrt_w, ref_norms=refs)
        if dropped_collinear:
            print(f"  Note: dropped (collinear): {', '.join(dropped_collinear)}")
            x_mat_raw = x_mat_raw[:, kept] if kept else None
            X = [X[j] for j in kept]
            control_cols = [control_cols[j] for j in kept]
    x_mat = x_mat_raw.copy() if x_mat_raw is not None else None

    # ── √w premultiplication of every analysis array ──
    y = y_raw * sqrt_w
    t = t_raw * sqrt_w
    m = m_raw * sqrt_w
    z_mat = z_mat_raw * sqrt_w[:, None]
    ones = sqrt_w.copy()
    if x_mat is not None:
        x_mat = x_mat * sqrt_w[:, None]
        exog = np.column_stack([ones, x_mat])              # [1, X]
        Z_base = np.column_stack([ones, z_mat, x_mat])     # [1, Z, X]
    else:
        exog = ones.reshape(-1, 1)                         # [1]
        Z_base = np.column_stack([ones, z_mat])            # [1, Z]

    z_crit = stats.norm.ppf(1 - alpha / 2)
    K_inst = z_mat.shape[1]  # number of instruments
    n_exog = exog.shape[1]

    # ── Identification checks (D9; paper Theorem 1 rank conditions, eq. 26) ──
    # After partialing out [1, X] (and the FE), every instrument, the
    # treatment and the mediator must keep variation, and the instruments must
    # not be collinear.
    def _rel_resid_var(v_w, ref):
        r_ = _residualize(v_w.reshape(n, -1), exog)[:, 0]
        return float(r_ @ r_) / ref ** 2 if ref > 0 else 0.0
    for j_, zc in enumerate(Z):
        if _rel_resid_var(z_mat[:, j_], ref_norm_role[zc]) < 1e-12:
            raise ValueError(f"instrument {zc} has no variation after the controls/fixed "
                             f"effects (σ_ZZ = 0, paper eq. 26).")
    if K_inst > 1:
        Zt_ = _residualize(z_mat, exog)
        refs_z = np.array([ref_norm_role[zc] for zc in Z])
        _, dz = _greedy_independent(Zt_, list(Z), ref_norms=refs_z)
        if dz:
            raise ValueError("instruments are collinear after the controls/fixed effects "
                             f"({', '.join(dz)}).")
    for vv, nm_, role_ in ((t, T, "treatment"), (m, M, "mediator")):
        if _rel_resid_var(vv, ref_norm_role[nm_]) < 1e-12:
            raise ValueError(f"{role_} {nm_} has no variation after the controls/fixed effects.")
    # [T̃, M̃] must have rank 2: a mediator that is a linear function of the
    # treatment and the controls leaves no conditional first stage for θ.
    r_mt = _residualize(m.reshape(n, -1), np.column_stack([exog, t]))[:, 0]
    if float(r_mt @ r_mt) / ref_norm_role[M] ** 2 < 1e-12:
        raise ValueError("the mediator is collinear with the treatment after the "
                         "controls/fixed effects.")

    # Not enough observations for the largest regression ([1, X, Z, T])?
    k_max = n_exog + K_inst + 1
    if n - k_max - df_fe_absorbed < 1:
        raise ValueError(
            f"not enough observations: n = {n}, but the largest regression has {k_max} "
            f"coefficients plus {df_fe_absorbed} absorbed fixed-effect degrees of freedom "
            f"(n must exceed their sum)."
        )

    # ══════════════════════════════════════════════════════════════════════
    # 1. FIRST STAGE: T ~ Z + X
    #    π̂₁ (paper eq. 27) and F_{T, Z|X}: joint Wald F for the EXCLUDED
    #    INSTRUMENTS only, partialing out X (and, via demeaning, any FE), under
    #    the SAME SE regime as θ̂ and τ̂. For K = 1 this is the Olea–Pflueger
    #    (2013) effective F; for K > 1 the OP rescaling is applied via
    #    _op_effective_F. The OP critical value comes from the simplified test
    #    at the estimated K_eff (OP 2013, eq. 4) at (op_tau, op_alpha).
    # ══════════════════════════════════════════════════════════════════════
    X_fs = Z_base  # [1, Z, X]   (Z columns 1 .. K_inst, X after)
    res_fs = _2sls_with_vcov(t, X_fs, X_fs, **se_kwargs, df_adjust=df_vcov)
    z_idx_fs = list(range(1, 1 + K_inst))
    beta_Z_fs = res_fs['coef'][z_idx_fs]
    V_ZZ_fs = res_fs['vcov'][np.ix_(z_idx_fs, z_idx_fs)]
    se_Z_fs = res_fs['se'][z_idx_fs]
    pi1 = float(beta_Z_fs[0]) if K_inst == 1 else [float(v) for v in beta_Z_fs]
    se_pi1 = float(se_Z_fs[0]) if K_inst == 1 else [float(v) for v in se_Z_fs]

    def _wald_over_K(b, V):
        # NaN unless the variance block is positive definite (an indefinite
        # spatial sandwich would give a negative or meaningless F).
        if not _is_pd(V):
            return float("nan")
        return float(b @ np.linalg.solve(V, b) / len(b))

    F_TZ = _wald_over_K(beta_Z_fs, V_ZZ_fs)
    z_resid_fs = _residualize(z_mat, exog)
    Q_zz_fs = z_resid_fs.T @ z_resid_fs / n
    if K_inst == 1:
        F_TZ_op = F_TZ
    else:
        F_TZ_op = (_op_effective_F(beta_Z_fs, V_ZZ_fs, Q_zz_fs)
                   if _is_pd(V_ZZ_fs) else float("nan"))
    K_eff_TZ, F_TZ_op_critval = _op_keff_critval(V_ZZ_fs, Q_zz_fs, op_tau, op_alpha)

    # ══════════════════════════════════════════════════════════════════════
    # 2. TREATMENT-MEDIATOR: M ~ T | Z + X (2SLS, instrument T with Z)
    # ══════════════════════════════════════════════════════════════════════
    X_gamma = np.column_stack([exog, t])       # [1, X, T]
    Z_gamma = Z_base                            # [1, Z, X]
    res_gamma = _2sls_with_vcov(m, X_gamma, Z_gamma, **se_kwargs, df_adjust=df_vcov)
    idx_T_in_gamma = X_gamma.shape[1] - 1  # T is appended last
    gamma_hat = res_gamma['coef'][idx_T_in_gamma]
    se_gamma = res_gamma['se'][idx_T_in_gamma]

    # ══════════════════════════════════════════════════════════════════════
    # 3. CONDITIONAL 2SLS: Y ~ M + T | Z + T + X  (paper eqs. 32–34)
    #    (instrument M with Z, conditioning on T as exogenous)
    # ══════════════════════════════════════════════════════════════════════
    X_theta = np.column_stack([exog, m, t])         # [1, X, M, T]
    Z_theta = np.column_stack([Z_base, t])           # [1, Z, X, T]
    res_theta = _2sls_with_vcov(y, X_theta, Z_theta, **se_kwargs, df_adjust=df_vcov)
    idx_M = n_exog      # M column index in X_theta
    idx_T = n_exog + 1  # T column index in X_theta
    theta_hat = res_theta['coef'][idx_M]
    tau_hat = res_theta['coef'][idx_T]
    se_theta = res_theta['se'][idx_M]
    se_tau = res_theta['se'][idx_T]

    # ══════════════════════════════════════════════════════════════════════
    # 4. TOTAL EFFECT: Y ~ T | Z + X (standard 2SLS)
    # ══════════════════════════════════════════════════════════════════════
    X_total = np.column_stack([exog, t])   # [1, X, T]
    Z_total = Z_base                        # [1, Z, X]
    res_total = _2sls_with_vcov(y, X_total, Z_total, **se_kwargs, df_adjust=df_vcov)
    idx_T_total = X_total.shape[1] - 1
    total_hat = res_total['coef'][idx_T_total]
    se_total = res_total['se'][idx_T_total]

    # ══════════════════════════════════════════════════════════════════════
    # 5. CONDITIONAL FIRST-STAGE F: M ~ Z + T + X  (paper eq. 32, Remark 3)
    #    F_{M, Z|T, X}: joint Wald F for the excluded instruments, partialing
    #    out T, X and the FE, under the same SE regime; OP effective F and the
    #    K_eff critical value with its own Ŵ₂.
    # ══════════════════════════════════════════════════════════════════════
    X_cfs = np.column_stack([exog, z_mat, t])  # [1, X, Z, T]
    res_cfs = _2sls_with_vcov(m, X_cfs, X_cfs, **se_kwargs, df_adjust=df_vcov)
    z_idx_cfs = list(range(n_exog, n_exog + K_inst))
    beta_Z_cfs = res_cfs['coef'][z_idx_cfs]
    V_ZZ_cfs = res_cfs['vcov'][np.ix_(z_idx_cfs, z_idx_cfs)]
    F_MZ_T = _wald_over_K(beta_Z_cfs, V_ZZ_cfs)
    # Per-instrument t² values (supplementary diagnostics, always stored;
    # NaN where the variance is not positive)
    F_MZ_T_individual = [
        float(res_cfs['coef'][i] ** 2 / res_cfs['vcov'][i, i])
        if res_cfs['vcov'][i, i] > 0 else float("nan") for i in z_idx_cfs
    ]
    incl_exog_cfs = np.column_stack([exog, t])   # [1, X, T]
    z_resid_cfs = _residualize(z_mat, incl_exog_cfs)
    Q_zz_cfs = z_resid_cfs.T @ z_resid_cfs / n
    if K_inst == 1:
        F_MZ_T_op = F_MZ_T
    else:
        F_MZ_T_op = (_op_effective_F(beta_Z_cfs, V_ZZ_cfs, Q_zz_cfs)
                     if _is_pd(V_ZZ_cfs) else float("nan"))
    K_eff_MZ_T, F_MZ_T_op_critval = _op_keff_critval(V_ZZ_cfs, Q_zz_cfs, op_tau, op_alpha)

    # ══════════════════════════════════════════════════════════════════════
    # 6. DECOMPOSITION CHECK: tau + theta * gamma = total  (paper eq. 30)
    # ══════════════════════════════════════════════════════════════════════
    indirect = theta_hat * gamma_hat
    decomp_check = tau_hat + indirect
    decomp_error = abs(decomp_check - total_hat)
    decomp_error_rel = decomp_error / max(abs(total_hat), abs(tau_hat) + abs(indirect), 1e-300)

    # ══════════════════════════════════════════════════════════════════════
    # 6b. ANALYTIC SE FOR THE INDIRECT EFFECT θγ (delta method)
    # ══════════════════════════════════════════════════════════════════════
    # The mediator and outcome equations share observations, so γ̂ (from
    # res_gamma) and θ̂ (from res_theta) have a non-zero cross-equation
    # covariance, computed under the active SE regime by a stacked-GMM
    # treatment (see _cross_equation_cov). The indirect-effect variance
    # follows from the delta method:
    #     Var(θγ) ≈ γ² Var(θ) + θ² Var(γ) + 2θγ Cov(θ, γ).
    cov_gamma_theta = _cross_equation_cov(
        res_gamma, res_theta,
        idx1=idx_T_in_gamma,        # γ is the T coefficient in the mediator eq.
        idx2=idx_M,                  # θ is the M coefficient in the outcome eq.
        df_adjust=df_vcov, **se_kwargs,
    )
    var_indirect = (gamma_hat ** 2) * (se_theta ** 2) \
                   + (theta_hat ** 2) * (se_gamma ** 2) \
                   + 2.0 * theta_hat * gamma_hat * cov_gamma_theta
    # Numerical safeguard: in pathological samples (e.g., poorly identified
    # first stage with cluster errors) the delta-method variance can come
    # out slightly negative due to cancellation. Floor at zero and warn.
    if var_indirect < 0:
        if var_indirect > -1e-8 * max((gamma_hat * se_theta) ** 2, 1e-300):
            var_indirect = 0.0
        else:
            print(f"  Warning: delta-method variance for the indirect effect "
                  f"is negative ({var_indirect:.2e}). This typically reflects "
                  f"a poorly identified first stage. Reporting NaN; consider "
                  f"ivmediate.bootstrap() for inference.")
            var_indirect = float('nan')
    se_indirect = float(np.sqrt(var_indirect)) if var_indirect == var_indirect else float('nan')

    # ══════════════════════════════════════════════════════════════════════
    # 7. κ-BASED SENSITIVITY (paper Section 5.1.1, Proposition 4, eqs. 36–42;
    #    App. F.1 algorithm and App. F.4 step 3)
    # ══════════════════════════════════════════════════════════════════════
    # Under weights, κ = Cov(ε_T, ε_Y) is a *population* covariance, estimable
    # by Hajek-form weighted moments (App. F.4 eq. 126). We work with the RAW
    # (un-pre-multiplied) arrays together with w_vec.
    #
    # Proposition 4(b) is stated for variables already projected onto the
    # control space [1, X]; we partial X out of T, M, Z (weighted FWL) before
    # computing the first-stage / IV residuals that feed κ̂*. This applies to
    # the no-LASSO path (user controls) and the LASSO path (selected set S).
    if x_mat_raw is not None and x_mat_raw.shape[1] > 0:
        X_partial = np.column_stack([np.ones(n), x_mat_raw])
    else:
        X_partial = np.ones((n, 1))

    def _partial_out(v_raw):
        return _wls_fit(X_partial, v_raw, w_vec)['residuals']

    t_orth = _partial_out(t_raw)
    m_orth = _partial_out(m_raw)
    z_mat_orth = np.column_stack([_partial_out(z_mat_raw[:, k]) for k in range(K_inst)])

    # First-stage residuals: ε̂_T from WLS of T on [1, Z] — RAW (orthogonalized) scale
    X_fs_kappa = np.column_stack([np.ones(n), z_mat_orth])
    eps_T_hat = _wls_fit(X_fs_kappa, t_orth, w_vec)['residuals']

    # Mediator IV residuals: ε̂_M = M − γ̂·T with the structural γ̂ (for K = 1
    # the IV ratio cov_w(Z, M)/cov_w(Z, T) on X-orthogonalized variables).
    eps_M_hat = m_orth - gamma_hat * t_orth

    # Outcome residuals: raw-scale ε̂_Y from the conditional 2SLS (the stored
    # residual is √w·ε̂_Y when weighted).
    eps_Y_hat = res_theta['resid'] / sqrt_w

    # Hajek-form weighted covariances/variances (paper eq. 126). c_TM_hat is
    # the COVARIANCE Cov(ε_T, ε_M) — c_TM in the paper (ρ_TM is the correlation
    # of Remark 3). Known as rho_TM_hat pre-2026-07; kept as an alias.
    c_TM_hat = _weighted_cov(eps_T_hat, eps_M_hat, w_vec)
    sigma_eT_hat = _weighted_std(eps_T_hat, w_vec)
    sigma_eY_hat = _weighted_std(eps_Y_hat, w_vec)

    # Breakeven κ* (eq. 37): the value of Cov(ε_T, ε_Y) that sets θ_adj to 0
    kappa_star = theta_hat * c_TM_hat

    # Inference breakeven κ*_CI (D1): the smallest κ, in the direction of κ*,
    # at which the κ-adjusted CI first includes 0 (eq. 38 with the near bound):
    #     κ*_CI = (θ̂ − sign(θ̂)·z·SE(θ̂))·ĉ_TM.
    # Uninformative (flagged) when the CI already contains 0; θ̂ = 0 gives 0.
    kappa_star_CI = (theta_hat - np.sign(theta_hat) * z_crit * se_theta) * c_TM_hat
    kappa_star_CI_informative = (bool(abs(theta_hat) > z_crit * se_theta)
                                 if np.isfinite(se_theta) else None)

    # Breakeven as a correlation (Remark 8, eq. 42)
    if sigma_eT_hat > 0 and sigma_eY_hat > 0:
        rho_TY_star = kappa_star / (sigma_eT_hat * sigma_eY_hat)
    else:
        rho_TY_star = np.nan

    # κ-adjusted estimates for κ₀ = f·κ̂* (App. F.4 step 3(e), eq. 127): at
    # f = 0 we have θ̂; at f = 1 the breakeven, θ̂_adj = 0.
    kappa_frac_grid = np.array([0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00])
    if np.isnan(kappa_star) or kappa_star == 0:
        # No meaningful breakeven (e.g., ĉ_TM ≈ 0). Show baseline only.
        kappa_grid = np.array([0.0])
        kappa_frac_used = np.array([0.0])
    else:
        kappa_grid = kappa_frac_grid * kappa_star
        kappa_frac_used = kappa_frac_grid
    kappa_adj_table = []
    for k0, frac in zip(kappa_grid, kappa_frac_used):
        th_adj = theta_hat - k0 / c_TM_hat
        ci_lo = th_adj - z_crit * se_theta
        ci_hi = th_adj + z_crit * se_theta
        kappa_adj_table.append(dict(
            kappa=k0, kappa_frac=frac,
            theta_adj=th_adj, ci_lo=ci_lo, ci_hi=ci_hi,
        ))
    # Strength of the conditional first stage (App. F.4 step 3): True iff the
    # effective F_MZ|T reaches the OP critical value at the user's op_tau; None
    # when either is not finite. print_results notes it below the κ table.
    if np.isfinite(F_MZ_T_op) and np.isfinite(F_MZ_T_op_critval):
        F_MZ_T_above_crit = bool(F_MZ_T_op >= F_MZ_T_op_critval)
    else:
        F_MZ_T_above_crit = None

    # ══════════════════════════════════════════════════════════════════════
    # 8. R²-BASED SENSITIVITY (paper App. F.2–F.3, eqs. 111–125)
    # ══════════════════════════════════════════════════════════════════════
    # Residualize w.r.t. [1, X, T] via FWL, then regress the residualized Y on
    # the residualized Z (conditional reduced form, eq. 111) and M on Z
    # (conditional first stage, eq. 112). The t-statistics use the SAME
    # inference regime as θ̂ (paper footnote 52). With K ≥ 2 only the first
    # instrument is used.
    W_fwl = np.column_stack([exog, t])  # [1, X, T], √w-premultiplied if weighted
    Y_tilde = lm_fit(W_fwl, y)['residuals']
    M_tilde = lm_fit(W_fwl, m)['residuals']
    Z_tilde = lm_fit(W_fwl, z_mat[:, 0])['residuals']   # first instrument

    # After FWL we regress on [1, Z_tilde]; df_adjust accounts for the columns
    # absorbed by FWL (everything in W_fwl except the intercept) plus the FE
    # degrees of freedom that count for this SE regime.
    X_Z = np.column_stack([ones, Z_tilde])
    df_adj_rv = (W_fwl.shape[1] - 1) + df_vcov

    res_rf = _2sls_with_vcov(Y_tilde, X_Z, X_Z, **se_kwargs, df_adjust=df_adj_rv)
    t_rf = res_rf['coef'][1] / res_rf['se'][1]
    F_reduced_form = t_rf ** 2

    res_fs_cond = _2sls_with_vcov(M_tilde, X_Z, X_Z, **se_kwargs, df_adjust=df_adj_rv)
    t_fs_cond = res_fs_cond['coef'][1] / res_fs_cond['se'][1]
    F_first_stage_cond = t_fs_cond ** 2

    # Degrees of freedom for the RV formulas: the residual df of the regression
    # of Y (or M) on [1, X, T, Z] (sensemakr's model$df.residual).
    df_rv = n - W_fwl.shape[1] - 1 - df_fe_absorbed

    # RV and XRV closed forms (App. F.3, eqs. 123–125; Cinelli–Hazlett):
    #   partial Cohen's f:  f_beta = |t| / sqrt(df)
    #   critical threshold: f_star = |t*_{alpha, df-1}| / sqrt(df - 1)
    #   XRV_{1,alpha} = max{0, (f_beta^2 - f_star^2) / (1 + f_beta^2)}
    #   RV_{1,alpha}  = 0                                  if f_beta - f_star < 0
    #                 = (1/2)(sqrt(fqa^4 + 4 fqa^2) - fqa^2),
    #                   with fqa = f_beta - f_star,        if f_beta < 1/f_star
    #                 = XRV_{1,alpha}                      otherwise (slack case)
    t_star_rv = abs(stats.t.ppf(alpha / 2.0, df_rv - 1))
    f_star_rv = t_star_rv / np.sqrt(df_rv - 1)

    def _xrv_ols(f_stat):
        """XRV_{1,alpha} for an individual OLS coefficient (f_stat = t^2)."""
        f_beta2 = abs(f_stat) / df_rv          # squared partial Cohen's f
        return max((f_beta2 - f_star_rv ** 2) / (1.0 + f_beta2), 0.0)

    def _rv_ols(f_stat):
        """RV_{1,alpha} for an individual OLS coefficient (f_stat = t^2)."""
        f_beta = np.sqrt(abs(f_stat) / df_rv)  # partial Cohen's f
        fqa = f_beta - f_star_rv
        if fqa < 0.0:
            return 0.0
        if f_beta < 1.0 / f_star_rv:
            return 0.5 * (np.sqrt(fqa ** 4 + 4.0 * fqa ** 2) - fqa ** 2)
        return _xrv_ols(f_stat)

    if np.isfinite(F_reduced_form) and np.isfinite(F_first_stage_cond):
        # Decomposition (eqs. 117–118): RV(θ) = min(RV_rf, RV_fs)
        xrv_rf = _xrv_ols(F_reduced_form)
        xrv_fs = _xrv_ols(F_first_stage_cond)
        rv_rf = _rv_ols(F_reduced_form)
        rv_fs = _rv_ols(F_first_stage_cond)
        XRV_theta = min(xrv_rf, xrv_fs)
        RV_theta = min(rv_rf, rv_fs)
    else:
        xrv_rf = xrv_fs = rv_rf = rv_fs = XRV_theta = RV_theta = float("nan")

    # ══════════════════════════════════════════════════════════════════════
    # 9. CONFIDENCE INTERVALS (Wald, normal critical values)
    # ══════════════════════════════════════════════════════════════════════
    ci_theta = (theta_hat - z_crit * se_theta, theta_hat + z_crit * se_theta)
    ci_tau = (tau_hat - z_crit * se_tau, tau_hat + z_crit * se_tau)
    ci_gamma = (gamma_hat - z_crit * se_gamma, gamma_hat + z_crit * se_gamma)
    ci_total = (total_hat - z_crit * se_total, total_hat + z_crit * se_total)
    if np.isfinite(se_indirect):
        ci_indirect = (indirect - z_crit * se_indirect,
                       indirect + z_crit * se_indirect)
    else:
        ci_indirect = (float('nan'), float('nan'))

    # ══════════════════════════════════════════════════════════════════════
    # 9b. ANDERSON–RUBIN (WEAK-INSTRUMENT-ROBUST) CONFIDENCE SETS
    # ══════════════════════════════════════════════════════════════════════
    # Test inversion under the SAME SE regime as the analytic inference.
    #   θ:       AR regression of (y − θ₀·m) on [1, X, T, Z]; Z block tested.
    #   γ:       AR regression of (m − γ₀·t) on [1, Z, X];    Z block tested.
    #   τ_total: AR regression of (y − τ₀·t) on [1, Z, X];    Z block tested.
    # K = 1: exact Fieller inversion; K ≥ 2: closed-form statistic located on
    # the grid and refined by Brent's method.
    ar_theta = ar_gamma = ar_total = None
    if ar_ci:
        W_ar_theta = np.column_stack([exog, t, z_mat])   # [1, X, T, Z]
        z_cols_theta = list(range(exog.shape[1] + 1,
                                  exog.shape[1] + 1 + K_inst))
        ar_theta = _ar_confidence_set(
            y, m, W_ar_theta, z_cols_theta, se_kwargs, df_vcov, alpha,
            center=theta_hat, scale=se_theta,
            grid_halfwidth=ar_grid_halfwidth, grid_points=ar_grid_points)
        z_cols_zb = list(range(1, 1 + K_inst))           # Z inside [1, Z, X]
        ar_gamma = _ar_confidence_set(
            m, t, Z_base, z_cols_zb, se_kwargs, df_vcov, alpha,
            center=gamma_hat, scale=se_gamma,
            grid_halfwidth=ar_grid_halfwidth, grid_points=ar_grid_points)
        ar_total = _ar_confidence_set(
            y, t, Z_base, z_cols_zb, se_kwargs, df_vcov, alpha,
            center=total_hat, scale=se_total,
            grid_halfwidth=ar_grid_halfwidth, grid_points=ar_grid_points)

    # Negative spatial variances (the kernel matrix is not PSD): one warning
    # per analysis, and only for reported quantities (θ, τ, γ, τ_total, π̂₁,
    # the first-stage F statistics and the RV t statistics). Intercepts,
    # controls and the AR regressions do not count; an AR set whose variance
    # turns negative is reported as unavailable instead.
    if se == "spatial" and not _point_only:
        rep_var = [res_gamma['vcov'][idx_T_in_gamma, idx_T_in_gamma],
                   res_theta['vcov'][idx_M, idx_M], res_theta['vcov'][idx_T, idx_T],
                   res_total['vcov'][idx_T_total, idx_T_total],
                   res_rf['vcov'][1, 1], res_fs_cond['vcov'][1, 1]]
        if any(v < 0 for v in rep_var) or not (_is_pd(V_ZZ_fs) and _is_pd(V_ZZ_cfs)):
            print("  Warning: the spatial (Conley) variance was negative for some "
                  "reported coefficients (the kernel matrix is not positive "
                  "semi-definite); their SEs and the statistics that use them "
                  "are reported as NaN.")

    # ══════════════════════════════════════════════════════════════════════
    # 10. ASSEMBLE RESULTS
    # ══════════════════════════════════════════════════════════════════════
    results = dict(
        version=__version__,
        # ── Sample info ──
        n=n,
        n_instruments=K_inst,
        n_covariates=len(X),
        n_dummies=len(dummy_col_names) if dummies else 0,
        se_type=se,
        alpha=alpha,
        fe=fe_cols,
        n_fe_groups=n_fe_groups,
        n_dropped_singletons=n_dropped_singletons,
        df_fe_absorbed=df_fe_absorbed,
        fe_nested_in_cluster=fe_nested_in_cluster,
        n_clusters=n_clusters,
        dropped_collinear=dropped_collinear,
        weighted=is_weighted,
        wtype=wtype if is_weighted else None,
        # Spatial-SE settings (None when se != 'spatial')
        spatial_cutoff=spatial_cutoff if se == 'spatial' else None,
        spatial_kernel=spatial_kernel if se == 'spatial' else None,
        spatial_metric=spatial_metric if se == 'spatial' else None,

        # ── Point estimates ──
        theta_hat=theta_hat, se_theta=se_theta, ci_theta=ci_theta,
        tau_hat=tau_hat, se_tau=se_tau, ci_tau=ci_tau,
        gamma_hat=gamma_hat, se_gamma=se_gamma, ci_gamma=ci_gamma,
        total_hat=total_hat, se_total=se_total, ci_total=ci_total,
        indirect=indirect, se_indirect=se_indirect, ci_indirect=ci_indirect,
        # Cross-equation covariance from stacked-GMM treatment (used by
        # the delta-method SE on θγ; exposed in case users want to inspect it)
        cov_gamma_theta=cov_gamma_theta,
        # First stage of T on Z (paper eq. 27): π̂₁ and its SE (lists for K ≥ 2)
        pi1=pi1, se_pi1=se_pi1,

        # ── Decomposition check ──
        decomp_check=decomp_check,
        decomp_error=decomp_error,
        decomp_error_rel=decomp_error_rel,

        # ── First-stage diagnostics ──
        F_TZ=F_TZ,
        F_MZ_T=F_MZ_T,
        F_MZ_T_individual=F_MZ_T_individual,
        # Olea–Pflueger effective F, estimated K_eff and the simplified-test
        # critical value at (op_tau, op_alpha)
        F_TZ_op=F_TZ_op,
        F_MZ_T_op=F_MZ_T_op,
        K_eff_TZ=K_eff_TZ,
        K_eff_MZ_T=K_eff_MZ_T,
        F_TZ_op_critval=F_TZ_op_critval,
        F_MZ_T_op_critval=F_MZ_T_op_critval,
        op_tau=op_tau,
        op_alpha=op_alpha,

        # ── Anderson–Rubin confidence sets (None when ar_ci=False) ──
        ar_ci=ar_ci,
        ar_theta=ar_theta,
        ar_gamma=ar_gamma,
        ar_total=ar_total,

        # ── κ-based sensitivity ──
        c_TM_hat=c_TM_hat,
        rho_TM_hat=c_TM_hat,   # deprecated alias (pre-2026-07 name); the
                               # paper's notation for Cov(ε̂_T, ε̂_M) is c_TM
        sigma_eT=sigma_eT_hat,
        sigma_eY=sigma_eY_hat,
        kappa_star=kappa_star,
        kappa_star_CI=kappa_star_CI,
        kappa_star_CI_informative=kappa_star_CI_informative,
        rho_TY_star=rho_TY_star,
        kappa_adj_table=kappa_adj_table,
        F_MZ_T_above_crit=F_MZ_T_above_crit,

        # ── R²-based sensitivity ──
        F_reduced_form=F_reduced_form,
        F_first_stage_cond=F_first_stage_cond,
        f_reduced_form=F_reduced_form,          # deprecated alias (= t², not Cohen's f)
        f_first_stage_cond=F_first_stage_cond,  # deprecated alias (= t², not Cohen's f)
        RV_theta=RV_theta,
        XRV_theta=XRV_theta,
        rv_rf=rv_rf,
        rv_fs=rv_fs,
        xrv_rf=xrv_rf,
        xrv_fs=xrv_fs,

        # ── Post-double-selection LASSO ──
        lasso=lasso,
        lasso_penalty=lasso_penalty if lasso else None,
        lasso_seed=lasso_seed if (lasso and lasso_penalty == "cv") else None,
        lasso_info=lasso_info,           # full structured output of _quadruple_lasso
        lasso_expansion_meta=rich_meta,  # rich-expansion metadata
        lasso_show_selection=show_selection if lasso else False,
        rich_col_names=list(rich_col_names) if lasso else None,

        # ── Raw-scale residuals (un-pre-multiplied), for bootstrap consumers ──
        # Each is on the X-orthogonalized scale, post-FE-demean. eps_T_hat is
        # from T regressed on [1, Z]; eps_M_hat is the IV residual M - γ̂T;
        # eps_Y_hat is from the conditional 2SLS.
        eps_T_hat=eps_T_hat,
        eps_M_hat=eps_M_hat,
        eps_Y_hat=eps_Y_hat,
    )

    return results


# ══════════════════════════════════════════════════════════════════════════════
#  Print Function
# ══════════════════════════════════════════════════════════════════════════════

def _validate_star_cutoffs(cutoffs):
    """Normalize and validate the (c***, c**, c*) cutoffs.

    Accepts any iterable of three numbers in (0, 1); returns them sorted
    ascending so that result[0] is the *** cutoff (tightest), result[1] is
    **, result[2] is *. Raises ValueError otherwise.
    """
    try:
        cuts = tuple(float(c) for c in cutoffs)
    except (TypeError, ValueError):
        raise ValueError(
            f"star_cutoffs must be an iterable of 3 numbers; got {cutoffs!r}"
        )
    if len(cuts) != 3:
        raise ValueError(
            f"star_cutoffs must have exactly 3 entries; got {len(cuts)}"
        )
    if not all(np.isfinite(c) and 0 < c < 1 for c in cuts):
        raise ValueError(
            f"star_cutoffs must all be strictly between 0 and 1; got {cuts}"
        )
    return tuple(sorted(cuts))


def _wald_pvalue(est, se):
    """Two-sided Wald p-value from a z-statistic = est/se.

    Returns np.nan if est, se are not both finite or se is non-positive.
    Uses scipy.stats.norm.sf for numerical stability in the tail.
    """
    if est is None or se is None:
        return np.nan
    if not (np.isfinite(est) and np.isfinite(se)) or se <= 0:
        return np.nan
    return 2.0 * stats.norm.sf(abs(est / se))


def _sig_stars(pval, cutoffs):
    """Return '***'/'**'/'*'/'' for a p-value given strictness cutoffs.

    cutoffs is the (already-sorted) output of _validate_star_cutoffs:
    cutoffs[0] is the *** threshold, cutoffs[1] is **, cutoffs[2] is *.
    """
    if pval is None or not np.isfinite(pval):
        return ""
    if pval < cutoffs[0]:
        return "***"
    if pval < cutoffs[1]:
        return "**"
    if pval < cutoffs[2]:
        return "*"
    return ""


def _format_cut(c):
    """Format a cutoff for the significance legend (no trailing zeros)."""
    return f"{c:g}"


def print_results(r, title="IV Mediation Analysis", boot=None,
                  stars=True, star_cutoffs=(0.01, 0.05, 0.10),
                  kappa_table=False, r2_sensitivity=False):
    """Print a formatted summary of results.

    Parameters
    ----------
    r : dict
        Output of ivmediate().
    title : str
    boot : dict or None
        If a bootstrap result (output of bootstrap()) is provided,
        the coefficient table reports bootstrap SEs/CIs for all five rows
        (τ_total, τ, θ, γ, θγ), the κ-inference breakeven uses the bootstrap
        SE on θ, and the bootstrap CIs of κ* and ρ*_TY are shown. Otherwise the
        analytic SEs from r are used. First-stage F-statistics, Olea–Pflueger
        critical values, AR sets and the RV/XRV diagnostics always use the
        analytic SE regime (classical / robust / cluster / spatial).
    stars : bool, default True
        Append significance stars (***/**/*) to the five coefficient rows
        of the MEDIATION ESTIMATES table, based on two-sided Wald p-values
        computed from estimate / SE (using whichever SE is in the row —
        analytic or bootstrap). Rows without a finite SE are left
        unstarred. Set False to turn stars off entirely. Stars are not
        applied to the F-statistic, RV/XRV, or κ-frame diagnostics.
    star_cutoffs : tuple of 3 floats, default (0.01, 0.05, 0.10)
        Significance thresholds for ***, **, * respectively. Any three
        numbers in (0, 1) are accepted; ordering does not matter (the
        strictest cutoff is always assigned ***).
    kappa_table : bool, default False
        Also print the κ-adjusted estimates (θ and its CI with κ_0 set to 0,
        10, 25, 50, 75, 90 and 100% of κ*). When the effective F_{M,Z|T,X} is
        below the Olea–Pflueger critical value (F_MZ_T_above_crit is False), a
        note below the table says that c_TM is then estimated imprecisely. The
        κ section always prints c_TM, κ*, κ*_CI and ρ*_TY.
    r2_sensitivity : bool, default False
        Also print the R²-based sensitivity section: the reduced-form F
        (= t²), RV(θ), XRV(θ) and their reduced-form and conditional
        first-stage components. The values are always in the results.
    """
    _check_bool("stars", stars)
    _check_bool("kappa_table", kappa_table)
    _check_bool("r2_sensitivity", r2_sensitivity)
    if stars:
        star_cutoffs = _validate_star_cutoffs(star_cutoffs)
    lvl = _level_label(r['alpha'])
    using_boot = boot is not None
    z_crit = stats.norm.ppf(1 - r['alpha'] / 2)

    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")
    info_line = (f"  n = {r['n']:,}  |  Instruments: {r['n_instruments']}  |  "
                 f"Covariates: {r['n_covariates']}  |  SE: {r['se_type']}")
    if r.get('n_clusters') is not None:
        info_line += f" (G = {r['n_clusters']} clusters)"
    if using_boot:
        info_line += f"  +  Bootstrap (B={boot['B']}, {boot['method']})"
    print(info_line)
    if r.get('weighted'):
        print(f"  Weighted: {r['wtype']}")
    if r.get('n_dummies', 0) > 0:
        print(f"  Dummy indicators: {r['n_dummies']}")
    if r.get('dropped_collinear'):
        print(f"  Dropped (collinear): {', '.join(r['dropped_collinear'])}")
    if r.get('fe'):
        line = (f"  FE: {', '.join(r['fe'])}  |  Groups: {r['n_fe_groups']}  |  "
                f"Singletons dropped: {r['n_dropped_singletons']}  |  "
                f"df absorbed: {r['df_fe_absorbed']}")
        nested = r.get('fe_nested_in_cluster')
        if nested and any(nested):
            line += ("  |  nested in clusters: " +
                     ", ".join(fe_ for fe_, nb in zip(r['fe'], nested) if nb) +
                     " (not counted in the CR1 df)")
        print(line)
    if r.get('se_type') == 'spatial':
        print(f"  Spatial: cutoff = {r.get('spatial_cutoff', '?')} "
              f"(metric={r.get('spatial_metric', '?')}, "
              f"kernel={r.get('spatial_kernel', '?')})")
    if r.get('lasso'):
        info = r.get('lasso_info') or {}
        n_sel = info.get('n_selected_union', 0)
        n_tot = info.get('n_total', 0)
        print(f"  LASSO: post-double selection, penalty={r['lasso_penalty']!r}, "
              f"selected {n_sel}/{n_tot} candidates")
    print(f"{'─'*70}")

    print(f"\n  MEDIATION ESTIMATES")
    ci_header = f"[{lvl} CI]"
    if using_boot:
        se_header = "Boot SE"
        ci_header = f"[Boot {lvl} CI]"
    else:
        se_header = "SE"
    print(f"  {'Parameter':<20} {'Estimate':>10} {se_header:>10} "
          f"{ci_header:>24}")
    print(f"  {'─'*64}")
    if using_boot:
        # Estimates come from r so the printed point estimates match the
        # structural table exactly (bootstrap point == structural point).
        rows = [
            ("τ_total (Wald)",  r['total_hat'],  boot['boot_se']['total_hat'],   boot['boot_ci']['total_hat']),
            ("τ (direct T→Y)",  r['tau_hat'],    boot['boot_se']['tau_hat'],     boot['boot_ci']['tau_hat']),
            ("θ (mediator→Y)",  r['theta_hat'],  boot['boot_se']['theta_hat'],   boot['boot_ci']['theta_hat']),
            ("γ (T→mediator)",  r['gamma_hat'],  boot['boot_se']['gamma_hat'],   boot['boot_ci']['gamma_hat']),
            ("θγ (indirect)",   r['indirect'],   boot['boot_se']['indirect_hat'], boot['boot_ci']['indirect_hat']),
        ]
    else:
        rows = [
            ("τ_total (Wald)",  r['total_hat'], r['se_total'],    r['ci_total']),
            ("τ (direct T→Y)",  r['tau_hat'],   r['se_tau'],      r['ci_tau']),
            ("θ (mediator→Y)",  r['theta_hat'], r['se_theta'],    r['ci_theta']),
            ("γ (T→mediator)",  r['gamma_hat'], r['se_gamma'],    r['ci_gamma']),
            ("θγ (indirect)",   r['indirect'],  r['se_indirect'], r['ci_indirect']),
        ]
    for name, est, se_val, ci in rows:
        if se_val is not None and np.isfinite(se_val):
            if stars:
                pval = _wald_pvalue(est, se_val)
                star_str = _sig_stars(pval, star_cutoffs)
            else:
                star_str = ""
            suffix = f" {star_str}" if star_str else ""
            print(f"  {name:<20} {est:10.4f} {se_val:10.4f} "
                  f"  [{ci[0]:9.4f}, {ci[1]:9.4f}]{suffix}")
        else:
            print(f"  {name:<20} {est:10.4f}")
    if stars:
        c3, c2, c1 = star_cutoffs  # tightest → loosest
        print(f"  Sig. codes: *** p<{_format_cut(c3)}, "
              f"** p<{_format_cut(c2)}, * p<{_format_cut(c1)}")

    # ── Decomposition check (silent unless violated) ───────────────────
    # For K=1 the identity τ + θγ = τ_total holds exactly (to rounding), so a
    # relative error above 1e-8 signals a genuine inconsistency. For K>1 the
    # objects are estimated equation-by-equation and coincide only
    # asymptotically; the note fires above 1e-2 relative.
    rel = r.get('decomp_error_rel')
    if rel is None:
        rel = r['decomp_error'] / max(abs(r['total_hat']),
                                      abs(r['tau_hat']) + abs(r['indirect']), 1e-300)
    if r['n_instruments'] == 1:
        if rel > 1e-8:
            print(f"\n  WARNING: decomposition identity violated.")
            print(f"  τ + θγ = {r['decomp_check']:.6f}  vs  "
                  f"τ_total = {r['total_hat']:.6f}")
            print(f"  relative error = {rel:.2e}  (threshold = 1e-08)")
            print(f"  For a single instrument this identity is exact; a nonzero "
                  f"error indicates an internal inconsistency. Please report.")
    else:
        if rel > 1e-2:
            print(f"\n  Note: with K>1 instruments, τ + θγ = {r['decomp_check']:.4f} "
                  f"and the Wald")
            print(f"  τ_total = {r['total_hat']:.4f} differ by {rel:.1e} (relative). "
                  f"This is expected — the")
            print(f"  objects are estimated equation-by-equation (not by joint GMM) "
                  f"and coincide")
            print(f"  only asymptotically. Not an error.")

    print(f"\n  FIRST-STAGE DIAGNOSTICS  (SE: {r['se_type']})")
    pi1, se_pi1 = r.get('pi1'), r.get('se_pi1')
    if pi1 is not None:
        if isinstance(pi1, list):
            body = ", ".join(f"{p_:.4f} ({s_:.4f})" for p_, s_ in zip(pi1, se_pi1))
            print(f"  π_1 (T on Z | X)        = {body}")
        else:
            print(f"  π_1 (T on Z | X)      = {pi1:.4f} ({se_pi1:.4f})")

    def _op_str(crit, keff):
        if crit is None or not np.isfinite(crit):
            return "  |  OP critical value: unavailable"
        return (f"  |  OP critical value = {crit:.2f}"
                f"  (K_eff={keff:.2f}, τ={r['op_tau']}, α={r['op_alpha']})")

    op_str_TZ = _op_str(r.get('F_TZ_op_critval'), r.get('K_eff_TZ', float('nan')))
    if r['n_instruments'] > 1:
        # K > 1: the label column is two characters wider (the widest label
        # is "F_{M, Z|T, X} (OP eff.)")
        print(f"  F_{{T, Z|X}}              = {r['F_TZ']:.2f}  (joint Wald)")
        print(f"  F_{{T, Z|X}} (OP eff.)    = {r['F_TZ_op']:.2f}{op_str_TZ}")
    else:
        # For K=1, joint Wald and OP F_eff coincide — single line is cleaner
        print(f"  F_{{T, Z|X}}            = {r['F_TZ']:.2f}{op_str_TZ}")

    op_str_MZ = _op_str(r.get('F_MZ_T_op_critval'), r.get('K_eff_MZ_T', float('nan')))
    if r['n_instruments'] > 1:
        print(f"  F_{{M, Z|T, X}}           = {r['F_MZ_T']:.2f}  (joint Wald)")
        print(f"  F_{{M, Z|T, X}} (OP eff.) = {r['F_MZ_T_op']:.2f}{op_str_MZ}")
        per_inst = "  ".join(f"{v:.1f}" for v in r['F_MZ_T_individual'])
        print(f"    per-instrument t²:      {per_inst}")
    else:
        print(f"  F_{{M, Z|T, X}}         = {r['F_MZ_T']:.2f}{op_str_MZ}")

    # ── Anderson–Rubin confidence sets ──────────────────────────────────
    if r.get('ar_ci') and r.get('ar_theta') is not None:
        print(f"\n  ANDERSON–RUBIN {lvl} CONFIDENCE SETS  "
              f"(weak-IV-robust, SE: {r['se_type']})")
        any_unbounded = False
        for label, key in (("θ (M → Y | T)", 'ar_theta'),
                           ("γ (T → M)", 'ar_gamma'),
                           ("τ_total", 'ar_total')):
            ar = r[key]
            if ar is None:
                continue
            print(f"  {label:<14} {_format_ar_set(ar)}")
            pv0 = ar['pval_at_zero']
            pv0_s = f"{pv0:.4f}" if np.isfinite(pv0) else "n/a"
            print(f"  {'':<14}   AR p-value at 0 = {pv0_s}"
                  f"   ({'exact Fieller inversion' if ar['method'] == 'fieller' else 'exact statistic, grid search'})")
            if _ar_unbounded_or_disconnected(ar):
                any_unbounded = True
        if any_unbounded:
            print(f"  Note: an unbounded or disconnected AR set is the expected "
                  f"outcome under weak")
            print(f"  identification — the data cannot bound the parameter at "
                  f"this confidence level.")
        print(f"  AR sets remain valid regardless of first-stage strength; "
              f"compare with the Wald")
        print(f"  CIs above, which require the corresponding effective F to "
              f"clear the OP threshold.")

    print(f"\n  κ-BASED SENSITIVITY (MCA diagnostics)")
    if r.get('weighted'):
        print(f"  (Hajek-weighted moments, wtype = {r['wtype']})")
    if using_boot:
        print(f"  (κ-frame CIs use bootstrap SE on θ)")
    print(f"  c_TM = Cov(ε_T, ε_M) = {r['c_TM_hat']:.4f}")
    print(f"  κ* (breakeven)       = {r['kappa_star']:.4f}")
    # κ̂*_CI = (θ̂ − sign(θ̂)·z·SE(θ̂))·ĉ_TM (D1); with a bootstrap, SE(θ̂) is the
    # bootstrap SE. Uninformative when the CI already contains 0.
    if using_boot:
        se_theta_use = boot['boot_se']['theta_hat']
    else:
        se_theta_use = r['se_theta']
    th = r['theta_hat']
    kappa_star_CI_use = (th - np.sign(th) * z_crit * se_theta_use) * r['c_TM_hat']
    informative = bool(np.isfinite(se_theta_use) and abs(th) > z_crit * se_theta_use)
    flag = "" if informative else "  (uninformative: the CI already contains 0)"
    print(f"  κ*_CI (inf. break.)  = {kappa_star_CI_use:.4f}{flag}")
    rho = r['rho_TY_star']
    rho_note = ("  (|ρ*| > 1: no feasible error correlation overturns θ)"
                if np.isfinite(rho) and abs(rho) > 1 else "")
    print(f"  ρ*_TY (as corr.)     = {rho:.4f}{rho_note}")
    if using_boot:
        for key, lab in (('kappa_star', 'κ*'), ('rho_TY_star', 'ρ*_TY')):
            ci_b = boot['boot_ci'].get(key)
            if ci_b is not None:
                print(f"  {lab} bootstrap {_level_label(boot['alpha'])} CI ({boot['ci_type']}) = "
                      f"[{ci_b[0]:.4f}, {ci_b[1]:.4f}]")
    if kappa_table:
        print(f"\n  κ-adjusted estimates  (κ_0 swept as fraction of breakeven κ*):")
        print(f"  {'% κ*':>6} {'κ_0':>12} {'θ_adj':>10} {'CI_lo':>10} {'CI_hi':>10}")
        for row in r['kappa_adj_table']:
            pct = row.get('kappa_frac', 0.0) * 100.0
            # With a bootstrap only the CI half-width changes (z·SE_boot(θ̂)).
            if using_boot:
                ci_lo_use = row['theta_adj'] - z_crit * se_theta_use
                ci_hi_use = row['theta_adj'] + z_crit * se_theta_use
            else:
                ci_lo_use = row['ci_lo']
                ci_hi_use = row['ci_hi']
            print(f"  {pct:>5.0f}% {row['kappa'] + 0.0:>12.4f} {row['theta_adj']:>10.4f} "
                  f"{ci_lo_use:>10.4f} {ci_hi_use:>10.4f}")
        # A conditional first stage below the OP critical value: ĉ_TM, and so
        # the κ quantities, are estimated imprecisely (App. F.4 step 3).
        if r.get('F_MZ_T_above_crit') is False:
            f_lab = "F_{M, Z|T, X}" if r['n_instruments'] == 1 else "F_{M, Z|T, X} (OP eff.)"
            ar_shown = bool(r.get('ar_ci') and r.get('ar_theta') is not None)
            print(f"  Note: {f_lab} = {r['F_MZ_T_op']:.2f} is below the OP critical value "
                  f"({r['F_MZ_T_op_critval']:.2f}, τ = {r['op_tau']}), so c_TM is estimated")
            print("  imprecisely; the AR set for θ above is robust to weak instruments."
                  if ar_shown else "  imprecisely.")

    if r2_sensitivity:
        print(f"\n  R²-BASED SENSITIVITY (reduced-form robustness)")
        print(f"  (t-statistics computed under {r['se_type']} SEs"
              + ("; with K > 1, RV/XRV use the first instrument)" if r['n_instruments'] > 1 else ")"))
        F_rf = r.get('F_reduced_form', r.get('f_reduced_form'))
        print(f"  F (reduced form, = t²) = {F_rf:.1f}")
        print(f"  RV(θ)                  = {r['RV_theta']:.4f}")
        print(f"  XRV(θ)                 = {r['XRV_theta']:.4f}")
        print(f"    (RV_rf={r['rv_rf']:.4f}, RV_fs={r['rv_fs']:.4f})")

    if r.get('lasso') and r.get('lasso_show_selection') and r.get('lasso_info'):
        info = r['lasso_info']
        details = info.get('details', {})
        print(f"\n  LASSO SELECTION DETAILS")
        print(f"  Penalty method: {r['lasso_penalty']}")
        if r.get('lasso_seed') is not None:
            print(f"  CV seed: {r['lasso_seed']}")

        def _print_set(label, d, max_show=20):
            if d is None:
                return
            n_sel = d.get('n_selected', 0)
            k_it = d.get('n_iter', 0)
            conv = (f"converged in {k_it} iteration{'s' if k_it != 1 else ''}"
                    if d.get('converged') else
                    f"did not converge after {k_it} iteration{'s' if k_it != 1 else ''}")
            print(f"  S_{label}  ({n_sel:>3} selected, λ = {d.get('lambda_value', float('nan')):.3f}, "
                  f"{conv})")
            if n_sel > 0:
                names = d['names'][:max_show]
                more = d['n_selected'] - len(names)
                line = ", ".join(names) + (f", ... (+{more} more)" if more > 0 else "")
                print(f"    {line}")
        _print_set("Y", details.get('Y'))
        _print_set("T", details.get('T'))
        _print_set("M", details.get('M'))
        for k, dz in enumerate(details.get('Z', []) or []):
            _print_set(f"Z[{k}]", dz)

        # Union summary
        union_names = info.get('union_names', [])
        n_union = info.get('n_selected_union', 0)
        n_total = info.get('n_total', 0)
        print(f"  Union S = S_Y ∪ S_T ∪ S_M ∪ S_Z: {n_union} of {n_total} candidates")
        if union_names:
            shown = union_names[:20]
            more = n_union - len(shown)
            line = ", ".join(shown) + (f", ... (+{more} more)" if more > 0 else "")
            print(f"    {line}")

    print(f"\n{'='*70}\n")


# ══════════════════════════════════════════════════════════════════════════════
#  Bootstrap Inference
#
#  Two methods are exposed:
#    'pairs' — resample observations with replacement. Valid for
#              cross-sectional iid data only.
#    'block' — resample whole clusters with replacement. Valid for
#              clustered data; requires clusterid.
#
#  Both methods refit the entire mediation analysis on each replicate and
#  produce bootstrap CIs for the structural parameters (θ̂, τ̂, γ̂, τ_total)
#  and the κ-frame quantities (κ̂*, ρ̂*_TY).
#
#  Bootstrap is intentionally NOT supported with lasso=True. Naive bootstrap
#  ignores selection variability and produces under-coverage. Use the
#  analytic post-double-selection inference from ivmediate,
#  which is asymptotically valid under Belloni-Chernozhukov-Hansen (2014).
#
#  For spatial dependence, the recommended approach is grid-cell clustering
#  per Bester, Conley, and Hansen (2011): partition observations into a
#  sufficient number of (otherwise arbitrary) spatial cells and use those
#  as the clusterid. The block bootstrap then provides honest inference
#  under spatial dependence without committing to a specific cutoff.
# ══════════════════════════════════════════════════════════════════════════════

def bootstrap(
    data, Y, T, M, Z, X=None, dummies=None, fe=None,
    se='robust', clusterid=None, alpha=0.05,
    weights=None, wtype='aweight',
    op_tau=0.10, op_alpha=0.05,
    lat=None, lon=None,
    spatial_cutoff=None, spatial_kernel=None, spatial_metric=None,
    # Bootstrap-specific
    method=None,
    B=999,
    seed=20260317,
    ci_type='basic',
    return_boot_dist=False,
    bootstrap_diagnostics=False,
    progress=True,
    verbose=False,
    _return_indices=False,
):
    """Bootstrap inference for IV mediation.

    Refits the mediation analysis on B resampled datasets and produces
    bootstrap CIs for each structural parameter and for κ* and ρ*_TY.

    Resampling is done on the estimation sample: rows with missing values in
    any used column, zero-weight rows and (with spatial SEs) rows without
    coordinates are removed first, so every 'pairs' replicate has n rows and
    every 'block' replicate has the G clusters of the estimation sample.
    Singleton FE groups are re-detected inside each replicate.

    Parameters
    ----------
    data, Y, T, M, Z, X, dummies, fe, se, clusterid, alpha, weights, wtype,
    op_tau, op_alpha, lat, lon, spatial_cutoff, spatial_kernel,
    spatial_metric :
        Forwarded to ivmediate. The arguments must specify the
        same model that bootstrap inference is targeted at.
    method : {'pairs', 'block', None}
        Resampling scheme:
          'pairs' — resample observations with replacement.
                    Valid for cross-sectional iid data ONLY.
          'block' — resample clusters (of `clusterid`) with replacement.
                    Valid under cluster dependence; requires clusterid.
                    Each drawn copy is its own cluster, and fixed effects
                    nested in the clusters get their own levels per copy
                    (as Stata's idcluster()).
        If None (default), uses 'block' when se='cluster', else 'pairs'.
    B : int
        Number of bootstrap replicates (at least 50). Default: 999.
    seed : int
        Master random seed. Replicate b draws from the b-th independent child
        stream of np.random.SeedSequence(seed) (SeedSequence.spawn), so
        different seeds share no replicates. The global NumPy random state is
        never used or changed. Default: 20260317.
    ci_type : {'basic', 'percentile', 'studentized'}
        Type of bootstrap confidence interval to construct.
          'basic'        — twice the estimate minus the replicate quantiles
                           q_{1-α/2} and q_{α/2}. Default; more robust to
                           skewness in the bootstrap distribution (which arises
                           naturally for products like θγ and ratios like κ*).
          'percentile'   — direct quantiles of the bootstrap distribution.
          'studentized'  — bootstrap-t with each replicate's analytic SE and the
                           original analytic SE. Falls back to 'basic' for κ*,
                           ρ*_TY, RV and XRV (no analytic SE).
        For 'basic' and 'percentile' without bootstrap_diagnostics the
        replicates compute point estimates only (no sandwich and no spatial
        neighbor search); 'studentized' and the diagnostics refit fully.
    return_boot_dist : bool
        If True, include the full B-vector for each bootstrapped parameter
        in the output dict under key 'boot_dist'. Default: False.
    bootstrap_diagnostics : bool
        If True, also bootstrap RV(θ) and XRV(θ). Default: False.
    progress : bool
        Print a progress line every ~5% of replicates. Default: True.
    verbose : bool
        Print extra diagnostic information. Default: False.

    Returns
    -------
    dict with keys:
        point         : dict of original-data point estimates (total_hat, tau_hat,
                       theta_hat, gamma_hat, indirect_hat, kappa_star,
                       rho_TY_star; plus RV_theta, XRV_theta with diagnostics)
        boot_se       : dict of bootstrap-estimated SEs (SD of the replicates)
        boot_ci       : dict of bootstrap CIs (of the chosen ci_type)
        analytic_se   : dict of analytic SEs from the original fit (for comparison)
        analytic_ci   : dict of analytic CIs from the original fit
        n_failed      : number of bootstrap replicates that failed
        n_est         : size of the estimation sample that is resampled
        method, B, seed, ci_type, alpha : bookkeeping for reproducibility
        boot_dist     : dict of B-vectors per parameter (only if return_boot_dist=True)

    Raises
    ------
    ValueError
        For invalid options (e.g., method='block' without clusterid, B < 50).
    RuntimeError
        If all bootstrap replicates fail.
    TypeError
        Implicit: lasso=True is rejected at the function-signature level
        because it is not in the kwarg list. Use ivmediate
        directly for LASSO-based analyses.

    Examples
    --------
    Basic usage with clustered data:
        boot = ivmediate.bootstrap(
            df, 'Y', 'T', 'M', ['Z'], X=['X1', 'X2'],
            se='cluster', clusterid='state',
        )
        ivmediate.print_bootstrap(boot)

    Pairs bootstrap on cross-sectional data:
        boot = ivmediate.bootstrap(
            df, 'Y', 'T', 'M', ['Z'], X=['X1'], method='pairs',
        )

    Returning the full bootstrap distribution for plotting:
        boot = ivmediate.bootstrap(
            df, 'Y', 'T', 'M', ['Z'], se='cluster', clusterid='state',
            return_boot_dist=True,
        )
        plt.hist(boot['boot_dist']['theta_hat'], bins=50)
    """
    import contextlib, io as _io

    # ── Options ──
    for _nm, _v in (("return_boot_dist", return_boot_dist),
                    ("bootstrap_diagnostics", bootstrap_diagnostics),
                    ("progress", progress), ("verbose", verbose),
                    ("_return_indices", _return_indices)):
        _check_bool(_nm, _v)
    _check_choice("se", se, _SE_TYPES)
    _check_open_unit("alpha", alpha)
    if method is None:
        method = 'block' if (se == 'cluster' and clusterid is not None) else 'pairs'
    _check_choice("method", method, ('pairs', 'block'))
    if method == 'block' and clusterid is None:
        raise ValueError(
            "method='block' (cluster bootstrap) requires clusterid. "
            "Provide clusterid or switch to method='pairs' for iid data."
        )
    if method == 'pairs' and clusterid is not None and se == 'cluster':
        print(
            f"  Note: method='pairs' was selected but clusterid={clusterid!r} "
            f"is provided. Pairs bootstrap is valid only under cross-sectional "
            f"independence; if your data are clustered, use method='block'."
        )
    if method == 'block' and se != 'cluster':
        print(f"  Note: block resampling by {clusterid!r} with se={se!r} "
              f"(the replicates do not use cluster SEs).")
    _check_choice("ci_type", ci_type, ('basic', 'percentile', 'studentized'))
    if not (isinstance(B, (int, np.integer)) and not isinstance(B, (bool, np.bool_))):
        raise ValueError(f"B must be an integer; got {B!r}.")
    if B < 50:
        raise ValueError(f"B should be at least 50 (got {B}); 999 is recommended.")
    if not (isinstance(seed, (int, np.integer)) and not isinstance(seed, (bool, np.bool_))) or seed < 0:
        raise ValueError(f"seed must be a non-negative integer; got {seed!r}.")

    # ── Estimation sample: resample only rows the analysis would use ──
    prep = _prepare_sample(data, Y, T, M, Z, X, dummies, fe, se, clusterid, weights,
                           lat, lon)
    mask = prep["mask"].copy()
    if method == 'block':
        if clusterid not in data.columns:
            raise ValueError(f"column(s) not found in data: {clusterid}.")
        cl_ok = _missing_text_to_nan(data[clusterid]).notna().to_numpy()
        mask &= cl_ok
    est_rows = np.flatnonzero(mask)
    data = data.loc[mask].reset_index(drop=True)
    if len(data) == 0:
        raise ValueError("No observations remaining after dropping missing values.")

    # ── Common kwargs forwarded to ivmediate ──
    # ar_ci=False: AR sets are analytic test inversions, not bootstrap targets.
    # X and dummies go in as the lists resolved on the original data, so a
    # glob cannot pick up the replicate's helper columns or the original
    # cluster id (which is not the cluster variable inside a block replicate).
    iv_kwargs = dict(
        Y=Y, T=T, M=M, Z=Z, X=list(prep["X"]), dummies=list(prep["dummies"]), fe=fe,
        se=se, clusterid=clusterid if se == 'cluster' else None, alpha=alpha,
        weights=weights, wtype=wtype,
        op_tau=op_tau, op_alpha=op_alpha,
        ar_ci=False,
        lat=lat, lon=lon, spatial_cutoff=spatial_cutoff,
        spatial_kernel=spatial_kernel, spatial_metric=spatial_metric,
    )

    # ── Original-data fit (point estimates, analytic SEs and CIs) ──
    if verbose:
        print(f"  Fitting original sample...")
    r_original = ivmediate(data=data, **iv_kwargs)

    # Order matches the structural-output coefficient table:
    # total → direct → mediator coefs → indirect → MCA diagnostics.
    params = ['total_hat', 'tau_hat', 'theta_hat', 'gamma_hat',
              'indirect_hat', 'kappa_star', 'rho_TY_star']
    if bootstrap_diagnostics:
        params += ['RV_theta', 'XRV_theta']

    # Some bootstrap-tracked parameters live under different keys in the
    # ivmediate result dict.
    result_key = {
        'theta_hat':    'theta_hat',
        'tau_hat':      'tau_hat',
        'gamma_hat':    'gamma_hat',
        'total_hat':    'total_hat',
        'indirect_hat': 'indirect',         # iv_mediation stores it as 'indirect'
        'kappa_star':   'kappa_star',
        'rho_TY_star':  'rho_TY_star',
        'RV_theta':     'RV_theta',
        'XRV_theta':    'XRV_theta',
    }

    point = {p: float(r_original[result_key[p]]) for p in params}
    _se_ind0 = float(r_original.get('se_indirect', float('nan')))
    analytic_se = {
        'theta_hat': float(r_original['se_theta']),
        'tau_hat':   float(r_original['se_tau']),
        'gamma_hat': float(r_original['se_gamma']),
        'total_hat': float(r_original['se_total']),
        # Delta-method SE for θγ (may be NaN in pathological samples).
        'indirect_hat': _se_ind0 if np.isfinite(_se_ind0) else None,
        # κ-frame and RV/XRV have no closed-form analytic SEs here.
        'kappa_star':  None,
        'rho_TY_star': None,
    }
    if bootstrap_diagnostics:
        analytic_se['RV_theta']  = None
        analytic_se['XRV_theta'] = None
    analytic_ci = {
        'theta_hat': list(r_original['ci_theta']),
        'tau_hat':   list(r_original['ci_tau']),
        'gamma_hat': list(r_original['ci_gamma']),
        'total_hat': list(r_original['ci_total']),
        'indirect_hat': (list(r_original['ci_indirect'])
                         if np.all(np.isfinite(r_original['ci_indirect']))
                         else None),
        'kappa_star':  None,
        'rho_TY_star': None,
    }
    if bootstrap_diagnostics:
        analytic_ci['RV_theta']  = None
        analytic_ci['XRV_theta'] = None

    # Replicates need their own SEs only for studentized CIs or RV/XRV.
    need_se = (ci_type == 'studentized') or bootstrap_diagnostics

    # ── Set up resampling (on the estimation sample) ──
    n = len(data)
    if method == 'block':
        cl_codes, unique_clusters = pd.factorize(data[clusterid], sort=False)
        G = len(unique_clusters)
        order = np.argsort(cl_codes, kind="stable")
        starts = np.searchsorted(cl_codes[order], np.arange(G))
        ends = np.searchsorted(cl_codes[order], np.arange(G), side="right")
        cluster_rows = [order[starts[g]:ends[g]] for g in range(G)]
        cluster_sizes = ends - starts
        # FE nested in the clusters get one level per drawn copy (as Stata's
        # bootstrap, cluster() idcluster() does): they stay nested in the
        # replicate's clusters, so replicate SEs use the original fit's df
        # convention, and a singleton group stays a singleton in each copy.
        nested_fe = {}
        for c in (prep["fe_cols"] or []):
            fe_codes = pd.factorize(data[c], sort=False)[0].astype(np.int64)
            pairs = np.unique(np.column_stack([fe_codes, cl_codes]), axis=0)
            if len(pairs) == len(np.unique(pairs[:, 0])):
                nested_fe[c] = (fe_codes, int(fe_codes.max()) + 1)
    child_seeds = np.random.SeedSequence(int(seed)).spawn(B)
    idx_record = [] if _return_indices else None

    # ── Bootstrap loop ──
    boot_dist = {p: np.full(B, np.nan) for p in params}
    if ci_type == 'studentized':
        # Per-replicate analytic SEs (theta/tau/gamma/total plus the
        # delta-method SE for the indirect effect).
        boot_se_replicate = {
            p: np.full(B, np.nan) for p in
            ['theta_hat', 'tau_hat', 'gamma_hat', 'total_hat', 'indirect_hat']
        }
    n_failed = 0
    progress_step = max(B // 20, 1) if progress else B + 1

    if progress:
        print(f"  Bootstrap ({method}, B={B})...")

    for b in range(B):
        rng_b = np.random.default_rng(child_seeds[b])

        # Construct the b-th replicate dataset
        if method == 'pairs':
            idx_b = rng_b.integers(0, n, size=n)
            if idx_record is not None:
                idx_record.append(idx_b)
            data_b = data.iloc[idx_b].reset_index(drop=True)
            iv_kwargs_b = iv_kwargs
        else:  # block
            g_b = rng_b.integers(0, G, size=G)
            if idx_record is not None:
                idx_record.append(g_b)
            concat_idx = np.concatenate([cluster_rows[g] for g in g_b])
            data_b = data.iloc[concat_idx].reset_index(drop=True)
            copy_id = np.repeat(np.arange(G, dtype=np.int64), cluster_sizes[g_b])
            data_b['_boot_cluster'] = copy_id
            for c, (codes_c, n_lev) in nested_fe.items():
                data_b[c] = copy_id * n_lev + codes_c[concat_idx]
            # Each drawn cluster is a separate cluster in the replicate
            iv_kwargs_b = dict(iv_kwargs)
            if se == 'cluster':
                iv_kwargs_b['clusterid'] = '_boot_cluster'

        try:
            with contextlib.redirect_stdout(_io.StringIO()):
                r_b = ivmediate(data=data_b, **iv_kwargs_b,
                                            _point_only=not need_se)
            for p in params:
                boot_dist[p][b] = float(r_b[result_key[p]])
            if ci_type == 'studentized':
                boot_se_replicate['theta_hat'][b] = float(r_b['se_theta'])
                boot_se_replicate['tau_hat'][b]   = float(r_b['se_tau'])
                boot_se_replicate['gamma_hat'][b] = float(r_b['se_gamma'])
                boot_se_replicate['total_hat'][b] = float(r_b['se_total'])
                boot_se_replicate['indirect_hat'][b] = float(r_b['se_indirect'])
        except Exception as e:
            n_failed += 1
            if verbose:
                print(f"    Replicate {b}: failed ({type(e).__name__}: {e})")

        if progress and (b + 1) % progress_step == 0:
            successful = b + 1 - n_failed
            print(f"    {b+1:>5}/{B}  ({100*(b+1)/B:>4.1f}%, "
                  f"{successful} successful)")

    if n_failed >= B:
        raise RuntimeError(
            f"All {B} bootstrap replicates failed. "
            "Check that the model specification is consistent with the data."
        )
    if n_failed > 0.10 * B:
        print(f"  Warning: {n_failed}/{B} bootstrap replicates failed "
              f"({100*n_failed/B:.1f}%). Bootstrap CIs may be unreliable.")

    # ── Aggregate to bootstrap SE and CI ──
    boot_se = {}
    boot_ci = {}
    for p in params:
        boot = boot_dist[p][~np.isnan(boot_dist[p])]
        boot_se[p] = float(np.std(boot, ddof=1)) if len(boot) > 1 else float('nan')
        if ci_type == 'percentile':
            q_lo, q_hi = np.quantile(boot, [alpha/2, 1 - alpha/2])
            boot_ci[p] = [float(q_lo), float(q_hi)]
        elif ci_type == 'basic':
            q_lo, q_hi = np.quantile(boot, [alpha/2, 1 - alpha/2])
            boot_ci[p] = [float(2 * point[p] - q_hi),
                          float(2 * point[p] - q_lo)]
        elif ci_type == 'studentized':
            # t_b = (θ̂_b - θ̂) / SE_b, then quantile of t_b → CI for θ̂
            if p in boot_se_replicate:
                mask_b = (~np.isnan(boot_dist[p])) & (~np.isnan(boot_se_replicate[p]))
                if mask_b.sum() < 10:
                    boot_ci[p] = [float('nan'), float('nan')]
                    continue
                t_stats = (boot_dist[p][mask_b] - point[p]) / boot_se_replicate[p][mask_b]
                t_stats = t_stats[np.isfinite(t_stats)]
                if len(t_stats) < 10:
                    boot_ci[p] = [float('nan'), float('nan')]
                    continue
                q_lo, q_hi = np.quantile(t_stats, [alpha/2, 1 - alpha/2])
                # Use the analytic SE from the original fit for the studentization
                point_se = analytic_se.get(p)
                if point_se is None or not np.isfinite(point_se):
                    boot_ci[p] = [float('nan'), float('nan')]
                else:
                    boot_ci[p] = [float(point[p] - q_hi * point_se),
                                  float(point[p] - q_lo * point_se)]
            else:
                # κ-frame and diagnostics don't have analytic SE; fall back to basic
                q_lo, q_hi = np.quantile(boot, [alpha/2, 1 - alpha/2])
                boot_ci[p] = [float(2 * point[p] - q_hi),
                              float(2 * point[p] - q_lo)]

    out = dict(
        point=point,
        boot_se=boot_se,
        boot_ci=boot_ci,
        analytic_se=analytic_se,
        analytic_ci=analytic_ci,
        n_failed=n_failed,
        n_est=n,
        method=method,
        B=B,
        seed=seed,
        ci_type=ci_type,
        alpha=alpha,
    )
    if return_boot_dist:
        out['boot_dist'] = boot_dist
    if _return_indices:
        # Draws (0-based: rows of the estimation sample for pairs, clusters in
        # order of first appearance for block) and the estimation sample's
        # 0-based row positions in the input data.
        out['_indices'] = np.vstack(idx_record)
        out['_rows'] = est_rows
    return out


def print_bootstrap(boot, title="Bootstrap Mediation",
                    stars=True, star_cutoffs=(0.01, 0.05, 0.10)):
    """Print a compact bootstrap-results table.

    Parameters
    ----------
    boot : dict
        Output of bootstrap().
    title : str
    stars, star_cutoffs
        Same semantics as in print_results. Stars are applied only to the
        five coefficient rows (θ, τ, γ, τ_total, θγ); the κ* / ρ*_TY /
        RV(θ) / XRV(θ) rows are sensitivity statistics, not coefficients,
        and remain unstarred regardless of their bootstrap distribution.
    """
    _check_bool("stars", stars)
    if stars:
        star_cutoffs = _validate_star_cutoffs(star_cutoffs)
    # Set of bootstrap keys that ARE coefficients (eligible for stars).
    coef_keys = {'theta_hat', 'tau_hat', 'gamma_hat',
                 'total_hat', 'indirect_hat'}

    print(f"\n{'='*70}")
    print(f"  {title}  (method={boot['method']}, B={boot['B']}, "
          f"seed={boot['seed']})")
    print(f"{'='*70}")
    print(f"  Bootstrap CI type: {boot['ci_type']}  "
          f"({_level_label(boot['alpha'])}, α = {boot['alpha']})")
    failed_str = f", {boot['n_failed']} failed" if boot['n_failed'] > 0 else ""
    print(f"  {boot['B'] - boot['n_failed']} successful replicates{failed_str}.")
    print(f"  {'─'*68}")

    label_map = {
        'theta_hat':    'θ (mediator→Y)',
        'tau_hat':      'τ (direct T→Y)',
        'gamma_hat':    'γ (T→mediator)',
        'total_hat':    'τ_total',
        'indirect_hat': 'θγ (indirect)',
        'kappa_star':  'κ*',
        'rho_TY_star': 'ρ*_TY',
        'RV_theta':    'RV(θ)',
        'XRV_theta':   'XRV(θ)',
    }
    print(f"  {'Parameter':<18} {'Estimate':>10} "
          f"{'Boot SE':>10} {'Boot CI':>26}")
    for p in boot['point']:
        label = label_map.get(p, p)
        est = boot['point'][p]
        sd  = boot['boot_se'][p]
        ci  = boot['boot_ci'][p]
        ci_str = f"[{ci[0]:>10.4f}, {ci[1]:>10.4f}]"
        if stars and p in coef_keys:
            pval = _wald_pvalue(est, sd)
            star_str = _sig_stars(pval, star_cutoffs)
        else:
            star_str = ""
        suffix = f" {star_str}" if star_str else ""
        print(f"  {label:<18} {est:>10.4f} {sd:>10.4f} {ci_str:>26}{suffix}")
    if stars:
        c3, c2, c1 = star_cutoffs
        print(f"  Sig. codes: *** p<{_format_cut(c3)}, "
              f"** p<{_format_cut(c2)}, * p<{_format_cut(c1)}")
    print(f"{'='*70}\n")
