# ivmediate for Python

The Python package `ivmediate` estimates causal mediation effects with a single instrument by the conditional 2SLS estimator of Ferrara, Gold, Heblich and Pinto (2026), with weak-instrument diagnostics, Anderson–Rubin confidence sets and sensitivity analysis.

[Main page](https://github.com/aferrara87/ivmediate) · [Python](https://github.com/aferrara87/ivmediate/blob/main/python/README.md) · [R](https://github.com/aferrara87/ivmediate/blob/main/r/README.md) · [Stata](https://github.com/aferrara87/ivmediate/blob/main/stata/README.md) · [MATLAB](https://github.com/aferrara87/ivmediate/blob/main/matlab/README.md)

The [main page](https://github.com/aferrara87/ivmediate#the-estimator) explains the estimator; this page covers the Python package. Paper references are to sections, remarks, propositions and appendices of Ferrara, Gold, Heblich and Pinto (2026).

## Installation

```bash
pip install "git+https://github.com/aferrara87/ivmediate#subdirectory=python"
```

pip needs git for this. Requirements: Python ≥ 3.9, numpy ≥ 1.21, scipy ≥ 1.7, pandas ≥ 1.3 and scikit-learn ≥ 1.0; pip installs any that are missing. psutil is optional; with it, spatial (Conley) standard errors adapt their memory use to the free memory.

Check the installed version:

```python
import ivmediate
print(ivmediate.__version__)
```

```text
2.0.0
```

Each result also records the version that produced it, in `r["version"]`.

## Quick start

Becker and Woessmann (2009): the effect of the share of Protestants (T) on income tax revenue per capita (Y) in 426 Prussian counties, mediated by literacy (M), with the distance to Wittenberg (Z) as the instrument and 11 county controls. The data are in [examples/becker_woessmann](https://github.com/aferrara87/ivmediate/tree/main/examples/becker_woessmann); run the code from the repository root, or run [bw.py](https://github.com/aferrara87/ivmediate/blob/main/examples/becker_woessmann/bw.py) from anywhere.

```python
import pandas as pd
import ivmediate

bw = pd.read_csv("examples/becker_woessmann/becker_woessmann_2009.csv")
r = ivmediate.ivmediate(
    bw, Y="inctaxpc", T="f_prot", M="f_rw", Z="kmwittenberg",
    X="f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop",
    se="robust",
)
ivmediate.print_results(r, title="Becker and Woessmann (2009)")
```

```text
  Resolved X (11 vars): f_young, f_jew, f_fem, f_ortsgeb, f_pruss, f_blind, f_deaf, f_dumb, ... (+3 more)

======================================================================
  Becker and Woessmann (2009)
======================================================================
  n = 426  |  Instruments: 1  |  Covariates: 11  |  SE: robust
──────────────────────────────────────────────────────────────────────

  MEDIATION ESTIMATES
  Parameter              Estimate         SE                 [95% CI]
  ────────────────────────────────────────────────────────────────
  τ_total (Wald)           0.5865     0.2425   [   0.1112,    1.0617] **
  τ (direct T→Y)          -0.2992     0.2141   [  -0.7187,    0.1204]
  θ (mediator→Y)           4.7091     2.1638   [   0.4681,    8.9502] **
  γ (T→mediator)           0.1881     0.0269   [   0.1354,    0.2407] ***
  θγ (indirect)            0.8857     0.4076   [   0.0869,    1.6845] **
  Sig. codes: *** p<0.01, ** p<0.05, * p<0.1

  FIRST-STAGE DIAGNOSTICS  (SE: robust)
  π_1 (T on Z | X)      = -0.0984 (0.0113)
  F_{T, Z|X}            = 75.69  |  OP critical value = 23.11  (K_eff=1.00, τ=0.1, α=0.05)
  F_{M, Z|T, X}         = 16.36  |  OP critical value = 23.11  (K_eff=1.00, τ=0.1, α=0.05)

  ANDERSON–RUBIN 95% CONFIDENCE SETS  (weak-IV-robust, SE: robust)
  θ (M → Y | T)  [0.4311, 10.2039]
                   AR p-value at 0 = 0.0333   (exact Fieller inversion)
  γ (T → M)      [0.1400, 0.2488]
                   AR p-value at 0 = 0.0000   (exact Fieller inversion)
  τ_total        [0.1258, 1.1030]
                   AR p-value at 0 = 0.0126   (exact Fieller inversion)
  AR sets remain valid regardless of first-stage strength; compare with the Wald
  CIs above, which require the corresponding effective F to clear the OP threshold.

  κ-BASED SENSITIVITY (MCA diagnostics)
  c_TM = Cov(ε_T, ε_M) = -89.2747
  κ* (breakeven)       = -420.4053
  κ*_CI (inf. break.)  = -41.7858
  ρ*_TY (as corr.)     = -0.2567

======================================================================
```

These are the numbers of the paper's Table 6, column 1; the [main page](https://github.com/aferrara87/ivmediate#example-becker-and-woessmann-2009) explains how to read them.

## Syntax

```python
ivmediate.ivmediate(data, Y, T, M, Z, X=None, dummies=None, fe=None,
                    se="robust", clusterid=None, alpha=0.05,
                    weights=None, wtype="aweight", op_tau=0.1, op_alpha=0.05,
                    ar_ci=True, ar_grid_halfwidth=200.0, ar_grid_points=4001,
                    lat=None, lon=None, spatial_cutoff=None,
                    spatial_kernel=None, spatial_metric=None,
                    lasso=False, lasso_penalty="plugin", lasso_keep=None,
                    lasso_seed=20260317, lasso_cv_folds=5,
                    expand_x=("square", "cube"), expand_dummies=("interact_with_X",),
                    expand_override=None, show_selection=False,
                    lasso_c_slack=1.1, lasso_alpha_p=None, lasso_n_iter_max=15,
                    lasso_tol=1e-6, lasso_lambda_override=None, lasso_verbose=False)

ivmediate.print_results(r, title="IV Mediation Analysis", boot=None, stars=True,
                        star_cutoffs=(0.01, 0.05, 0.1), kappa_table=False,
                        r2_sensitivity=False)

ivmediate.bootstrap(data, Y, T, M, Z, X=None, dummies=None, fe=None,
                    se="robust", clusterid=None, alpha=0.05,
                    weights=None, wtype="aweight", op_tau=0.1, op_alpha=0.05,
                    lat=None, lon=None, spatial_cutoff=None,
                    spatial_kernel=None, spatial_metric=None,
                    method=None, B=999, seed=20260317, ci_type="basic",
                    return_boot_dist=False, bootstrap_diagnostics=False,
                    progress=True, verbose=False)

ivmediate.print_bootstrap(boot, title="Bootstrap Mediation", stars=True,
                          star_cutoffs=(0.01, 0.05, 0.1))
```

- `ivmediate()` estimates the model and returns a dict (see [Output](#output)). While it runs, it prints short notes: the resolved list of controls, dropped observations and controls, ignored options.
- `print_results()` prints the tables for that dict.
- `bootstrap()` refits the model on resampled data and returns a dict; `print_bootstrap()` prints it, and `print_results(r, boot=boot)` puts the bootstrap standard errors and CIs into the main table.
- `data` is a pandas DataFrame and all variables are given by column name. `Z` is one name or a list of names (`Z=["z1", "z2"]`); unlike `X`, a string is never split.

## Options

Every option after `Z` has a default; pass it by name. Booleans must be `True` or `False` (a string such as `"false"` is an error).

### Model and controls

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `data` | required | pandas DataFrame, one row per observation; it is not modified. | |
| `Y`, `T`, `M` | required | Column names of the outcome, the treatment and the mediator; three different columns. | Section 2 |
| `X` | `None` | Controls: a list of column names or one string of names separated by spaces or line breaks. Names may be globs (`*`, `?`, `[...]`), e.g., `"age income reg_*"`. Globs skip the Y, T, M, Z, FE, cluster, weight and coordinate columns (a note lists them); naming Y, T, M or Z is an error, and so is a glob that matches nothing. Exactly collinear controls are dropped with a note, keeping the earlier one. | |
| `dummies` | `None` | Categorical variables turned into indicators, leaving out the first category in sorted order; given like `X`. The indicators are named `_D_<variable>_<category>`, e.g., `_D_region_2`. Meant for variables with few categories (up to about 50). | |

Rows with a missing value in any column the model uses are dropped. Empty or whitespace-only text in FE, cluster and `dummies` columns counts as missing. Infinite values are an error, and so is text in a numeric role.

### Fixed effects and weights

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `fe` | `None` | One or two fixed-effect columns (`fe="firm"` or `fe=["firm", "year"]`), absorbed by within-group demeaning. Two-way effects use the exact projection when the dimension with fewer levels has at most 3,000 levels, and alternating projections with Irons–Tuck acceleration (relative tolerance 1e-12) otherwise. Singleton groups are dropped repeatedly until none are left (a note gives the count). The absorbed degrees of freedom are the rank of the FE dummies: G − 1 one-way, G_1 + G_2 − c − 1 two-way with c connected components. | |
| `weights` | `None` | Column of non-negative weights. Rows with zero weight are dropped (with a note) before singletons are detected; negative weights are an error. With FE, the group means are weighted. | |
| `wtype` | `"aweight"` | `"aweight"` (analytic) or `"pweight"` (probability). Both give the same results except with `se="classical"`: `"pweight"` then switches to robust standard errors (with a note), and `"aweight"` keeps classical ones and prints a reminder that robust or cluster standard errors are the safer choice. | |

### Standard errors

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `se` | `"robust"` | `"classical"`, `"robust"` (HC1), `"cluster"` (CR1, needs `clusterid`) or `"spatial"` (Conley, needs `lat`, `lon` and `spatial_cutoff`). The same regime applies to the Wald CIs, the F statistics and their critical values, and the AR sets. | Remark 4; Appendix D.1 |
| `clusterid` | `None` | Cluster column; required with `se="cluster"` and ignored (with a note) otherwise. At least 2 clusters. FE nested in the clusters do not count in the CR1 small-sample factor, as in reghdfe and fixest. | |
| `alpha` | `0.05` | Significance level of the Wald CIs, the AR sets, the κ statistics and the RV/XRV threshold; strictly between 0 and 1. | |
| `lat`, `lon` | `None` | Coordinate columns, in degrees by default. Required with `se="spatial"` and an error with any other `se`. Rows with missing coordinates are dropped. | |
| `spatial_cutoff` | `None` | Kernel bandwidth, required with `se="spatial"`: kilometers under the haversine metric, otherwise the units of `lat` and `lon`. With the uniform kernel, a cutoff that covers every pair of observations is an error; with either kernel, one that covers more than half of the pairs gives a warning. Conley standard errors depend on the cutoff, so report several. | |
| `spatial_kernel` | `None` | `"uniform"` (the default with `se="spatial"`) or `"bartlett"`, which tapers linearly to zero at the cutoff. | |
| `spatial_metric` | `None` | `"haversine"` (the default with `se="spatial"`; great-circle distance in km from degrees) or `"euclidean"` (projected coordinates). | |

The spatial arguments are an error with any other `se`. A Conley variance estimate need not be positive: a negative variance gives a NaN standard error, and an AR set whose variance turns negative anywhere in its search range is reported as unavailable.

### Instruments and Anderson–Rubin sets

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `Z` | required | The instrument's column name, or a list or tuple of names for several instruments. The names must differ from Y, T and M and from each other. | Section 3.3, Proposition 2; Section 4.1, Theorem 1 |
| `ar_ci` | `True` | Anderson–Rubin confidence sets for θ, γ and τ_total, valid however weak the instrument, with the same standard-error regime and `alpha` as the Wald CIs. With one instrument the set is computed exactly (Fieller inversion); with several, the exact statistic is scanned on a grid and the boundaries are refined by Brent's method. There is no AR set for τ, which would need subset-AR methods. `False` skips them. | Remark 3 discussion (Section 4.2); Table 6 notes (Section 8) |
| `ar_grid_halfwidth` | `200.0` | Half-width of the AR search grid, in analytic standard errors of the parameter around its estimate (several instruments; also the range checked for negative Conley variances). | |
| `ar_grid_points` | `4001` | Number of grid points; an integer of at least 3. | |

### Weak-instrument diagnostics

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `op_tau` | `0.1` | Largest 2SLS bias the Olea–Pflueger test tolerates, as a share of the OLS bias benchmark (0.1 = 10%). It sets the critical values for both first stages: T on Z, and M on Z given T (the conditional first stage). | Remark 3 (Section 4.2); Appendix F.4, step 2 |
| `op_alpha` | `0.05` | Size of the Olea–Pflueger test. The critical values come from the simplified test at the estimated effective degrees of freedom K_eff (1 with one instrument), computed for any `op_tau` and `op_alpha` in (0, 1). Separate from `alpha`. | Remark 3 (Section 4.2); Appendix F.4, step 2 |

### Sensitivity and printing

The κ and R² sensitivity statistics are always computed and stored in the result; these `print_results()` options control what is printed.

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `r` | required | The dict returned by `ivmediate()`. | |
| `title` | `"IV Mediation Analysis"` | Title line of the printout. | |
| `boot` | `None` | A `bootstrap()` result. The coefficient table then shows bootstrap standard errors and CIs, κ*_CI uses the bootstrap standard error of θ, and the bootstrap CIs of κ* and ρ*_TY are added. The F statistics, critical values, AR sets and RV/XRV stay analytic. | |
| `stars` | `True` | Significance stars on the five coefficient rows, from two-sided Wald p-values with the standard error shown in the row. Rows without a finite standard error get none. | |
| `star_cutoffs` | `(0.01, 0.05, 0.1)` | p-value thresholds for `***`, `**` and `*`: any three numbers in (0, 1), in any order. | |
| `kappa_table` | `False` | Adds the κ-adjusted estimates: θ and its CI with κ_0 at 0, 10, 25, 50, 75, 90 and 100% of κ*. When the conditional first stage is below its Olea–Pflueger critical value (`F_MZ_T_above_crit` is `False`), a note under the table says that c_TM is estimated imprecisely. | Section 5.1.1, Proposition 5; Appendix F.1 and F.4 |
| `r2_sensitivity` | `False` | Adds the R²-based section: the reduced-form F (= t²), RV(θ), XRV(θ) and their reduced-form and first-stage components. | Appendix F.2 (Definitions 2–3) and F.3 |

### LASSO

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `lasso` | `False` | Post-double-selection ("quadruple") LASSO before the 2SLS: one LASSO each for Y, T, M and every instrument on a rich expansion of `X` and `dummies`; the union of the selections becomes the set of controls. Each LASSO has an unpenalized intercept, weights are normalized to mean 1, and the penalty loadings allow for heteroskedasticity, and for clustering or spatial correlation under those standard errors. Not available in `bootstrap()`. | Remark 6 (Section 4.2); Appendix H (Algorithm 1 in H.2) |
| `lasso_penalty` | `"plugin"` | `"plugin"`: the iterated plug-in penalty of Belloni, Chernozhukov and Hansen with post-LASSO residual loadings (deterministic). `"cv"`: K-fold cross-validation (folds by cluster with `se="cluster"`). | |
| `lasso_keep` | `None` | Names that are always kept in the controls (partialed out before each LASSO). They refer to the rich expansion: originals (`"x1"`), powers (`"x1^2"`, `"x1^3"`) and interactions (`"x1 × x2"`, `"_D_region_2 × x1"`). Globs work (`"x1*"` matches x1, x1^2, x1^3 and the interactions that start with x1). Interaction names contain spaces, so give them in a list. | |
| `lasso_seed` | `20260317` | Seed of the CV folds; the four LASSOs use seed + 0 (Y), + 1 (T), + 2 (M) and + 3 + k (instrument k, from 0). Not used by the plug-in penalty. | |
| `lasso_cv_folds` | `5` | Number of CV folds; an integer of at least 2. | |
| `expand_x` | `("square", "cube")` | Terms added for the `X` variables: `"square"` and `"cube"` (variables with more than two distinct values only; others are skipped with a note) and `"interact_with_X"` (products of pairs of `X` variables). | Appendix H |
| `expand_dummies` | `("interact_with_X",)` | Terms added for the `dummies` indicators: `"interact_with_X"` (products with every `X` variable) and `"interact_with_dummies"` (products of indicators of different variables). Powers have no effect on indicators. | Appendix H |
| `expand_override` | `None` | Terms for single variables, overriding the two options above, e.g., `{"age": ("square",)}`. Keys are `X` names or indicator names such as `"_D_region_2"`; other keys are ignored. | |
| `show_selection` | `False` | `print_results()` also prints the LASSO selections (sets, penalties, convergence). They are in the result either way. | Appendix H.4 |
| `lasso_c_slack` | `1.1` | Slack c of the plug-in penalty λ = 2c√n Φ⁻¹(1 − α_p/(2p)), where n is the number of observations (with unequal weights, the Kish effective sample size). | |
| `lasso_alpha_p` | `None` | Tail probability α_p of the plug-in penalty; `None` means 0.1/log(n), with n the number of observations also under clustering (lassopack's cluster default is 0.1/log(G)). | |
| `lasso_n_iter_max` | `15` | Maximum number of plug-in loading iterations. | |
| `lasso_tol` | `1e-6` | The iterations stop once the largest coefficient change is below this. | |
| `lasso_lambda_override` | `None` | A fixed λ that bypasses the plug-in and CV rules, e.g., to replicate a published λ. | |
| `lasso_verbose` | `False` | Prints each LASSO's iterations. | |

### Bootstrap

`bootstrap()` takes the model arguments of `ivmediate()` (`data` to `spatial_metric`, except `ar_ci` and the AR grid; no LASSO) and these:

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `method` | `None` | `"pairs"` resamples observations (independent data only); `"block"` resamples the clusters of `clusterid`, treats each drawn copy as its own cluster and gives FE nested in the clusters new levels in each copy. `None` means `"block"` with `se="cluster"` and `"pairs"` otherwise. | |
| `B` | `999` | Number of replicates; an integer of at least 50. | |
| `seed` | `20260317` | Non-negative integer. Replicate b draws from the b-th child stream of numpy's `SeedSequence(seed)`; numpy's global random state is never used or changed. | |
| `ci_type` | `"basic"` | `"basic"` (twice the estimate minus the replicate quantiles), `"percentile"` (the quantiles) or `"studentized"` (bootstrap-t with each replicate's analytic standard error; κ*, ρ*_TY, RV and XRV fall back to basic). | |
| `return_boot_dist` | `False` | Adds the B replicate values of every parameter under `"boot_dist"`. | |
| `bootstrap_diagnostics` | `False` | Also bootstraps RV(θ) and XRV(θ). | Appendix F.2 |
| `progress` | `True` | Prints a progress line every 5% of the replicates. | |
| `verbose` | `False` | Prints extra details, such as why a replicate failed. | |

The bootstrap resamples the estimation sample: rows with missing values, zero weights or (with `"block"`) a missing cluster id are removed first. Singleton FE groups are found again in each replicate. With `"basic"` or `"percentile"` CIs and no diagnostics, the replicates compute point estimates only, which is fast. Failed replicates are left out and counted (with a warning above 10%); if all fail, `bootstrap()` raises an error. `print_bootstrap(boot, title="Bootstrap Mediation", stars=True, star_cutoffs=(0.01, 0.05, 0.1))` prints the result; its stars work as in `print_results()` and go on the five coefficient rows only.

## Output

### The printout

| Block | What it shows | Paper |
|---|---|---|
| Header | n, number of instruments and controls, standard-error type (and number of clusters); weights, indicators, dropped collinear controls, FE (groups, singletons, nesting in clusters), spatial and LASSO settings when used. | |
| MEDIATION ESTIMATES | τ_total, τ, θ, γ and θγ with standard errors, Wald CIs and stars; the standard error of θγ comes from the delta method. | Section 4.2; Remark 4, Appendix D.1 |
| FIRST-STAGE DIAGNOSTICS | π_1 (T on Z given X), F_{T, Z\|X} and the conditional F_{M, Z\|T, X} with their Olea–Pflueger critical values and K_eff. With several instruments: the joint Wald F, the effective F and each instrument's t². | Remark 3 (Section 4.2); Appendix F.4, step 2 |
| ANDERSON–RUBIN CONFIDENCE SETS | The AR sets for θ, γ and τ_total and the AR p-value at 0; unbounded sets are flagged as the sign of a weak instrument. | Remark 3 discussion (Section 4.2) |
| κ-BASED SENSITIVITY | c_TM = Cov(ε_T, ε_M), the breakeven κ*, the inference breakeven κ*_CI (flagged as uninformative when the CI for θ already contains 0) and κ* as a correlation, ρ*_TY. | Section 5.1.1, Proposition 5; Remark 8 |
| κ-adjusted estimates | Only with `kappa_table=True`. | Appendix F.4 |
| R²-BASED SENSITIVITY | Only with `r2_sensitivity=True`. | Appendix F.2 and F.3 |
| LASSO SELECTION DETAILS | Only with `lasso=True, show_selection=True`: each LASSO's selected controls, λ and convergence, and their union. | Appendix H.4 |

With `kappa_table=True` and `r2_sensitivity=True`, the quick-start printout gains these two blocks after the κ section:

```python
ivmediate.print_results(r, title="Becker and Woessmann (2009)",
                        kappa_table=True, r2_sensitivity=True)
```

```text
  κ-adjusted estimates  (κ_0 swept as fraction of breakeven κ*):
    % κ*          κ_0      θ_adj      CI_lo      CI_hi
      0%       0.0000     4.7091     0.4681     8.9502
     10%     -42.0405     4.2382    -0.0029     8.4793
     25%    -105.1013     3.5318    -0.7092     7.7729
     50%    -210.2027     2.3546    -1.8865     6.5956
     75%    -315.3040     1.1773    -3.0638     5.4183
     90%    -378.3648     0.4709    -3.7701     4.7120
    100%    -420.4053     0.0000    -4.2411     4.2411
  Note: F_{M, Z|T, X} = 16.36 is below the OP critical value (23.11, τ = 0.1), so c_TM is estimated
  imprecisely; the AR set for θ above is robust to weak instruments.

  R²-BASED SENSITIVITY (reduced-form robustness)
  (t-statistics computed under robust SEs)
  F (reduced form, = t²) = 4.5
  RV(θ)                  = 0.0078
  XRV(θ)                 = 0.0016
    (RV_rf=0.0078, RV_fs=0.0972)
```

### The returned dict

```python
ar = r["ar_theta"]
print(f"theta = {r['theta_hat']:.4f} (SE {r['se_theta']:.4f})")
print(f"AR set for theta: {ar['shape']} [{ar['lower']:.4f}, {ar['upper']:.4f}]")
print(f"F_MZ_T = {r['F_MZ_T_op']:.2f}, OP critical value = {r['F_MZ_T_op_critval']:.2f}, "
      f"above: {r['F_MZ_T_above_crit']}")
```

```text
theta = 4.7091 (SE 2.1638)
AR set for theta: interval [0.4311, 10.2039]
F_MZ_T = 16.36, OP critical value = 23.11, above: False
```

| Keys | What they hold |
|---|---|
| `total_hat`, `tau_hat`, `theta_hat`, `gamma_hat`, `indirect` | The estimates of τ_total, τ, θ, γ and θγ. |
| `se_total`, `se_tau`, `se_theta`, `se_gamma`, `se_indirect` | Their standard errors. |
| `ci_total`, `ci_tau`, `ci_theta`, `ci_gamma`, `ci_indirect` | Their Wald CIs, as (lower, upper). |
| `pi1`, `se_pi1` | The first-stage coefficient of T on Z given X and its standard error (lists with several instruments). |
| `F_TZ`, `F_MZ_T` | First-stage F statistics of T on Z and of M on Z given T (joint Wald with several instruments). |
| `F_TZ_op`, `F_MZ_T_op`, `F_TZ_op_critval`, `F_MZ_T_op_critval`, `K_eff_TZ`, `K_eff_MZ_T` | Olea–Pflueger effective F statistics, their critical values and effective degrees of freedom. |
| `F_MZ_T_above_crit` | `True` if the conditional first stage reaches its Olea–Pflueger critical value, `False` if not (then the κ table is imprecise), `None` if either value is not finite. |
| `ar_theta`, `ar_gamma`, `ar_total` | The AR sets (below); `None` with `ar_ci=False`. |
| `c_TM_hat`, `kappa_star`, `kappa_star_CI`, `rho_TY_star` | c_TM, κ*, κ*_CI and ρ*_TY. |
| `kappa_star_CI_informative` | `True` if the Wald CI for θ excludes 0, so that κ*_CI is informative; `False` if it contains 0; `None` without a finite standard error. |
| `kappa_adj_table` | The κ table: a list of dicts with `kappa`, `kappa_frac`, `theta_adj`, `ci_lo` and `ci_hi`. |
| `F_reduced_form`, `F_first_stage_cond`, `RV_theta`, `XRV_theta`, `rv_rf`, `rv_fs`, `xrv_rf`, `xrv_fs` | The R²-based sensitivity statistics and their components. |
| `n`, `n_clusters`, `n_instruments`, `n_covariates`, `n_dummies` | Observations used, clusters (`None` unless `se="cluster"`), instruments, controls used (after collinear drops; with LASSO, those selected) and indicators. |
| `n_fe_groups`, `n_dropped_singletons`, `df_fe_absorbed` | FE groups (summed over dimensions), dropped singletons and absorbed degrees of freedom. |
| `fe_nested_in_cluster` | One `True`/`False` per FE dimension: whether it is nested in the clusters (and so left out of the CR1 factor). `None` without FE or without cluster standard errors. |
| `dropped_collinear` | Controls dropped as exactly collinear. |
| `lasso_info` | With `lasso=True`: `union_names` (the selected controls), `n_selected_union`, `n_total` (candidates), and `details` with each LASSO's selection, λ and convergence. `None` otherwise. |
| `eps_T_hat`, `eps_M_hat`, `eps_Y_hat` | The residuals behind the κ statistics, as numpy arrays over the estimation sample. |
| `version` | The package version. |

Each AR set is a dict. Its `shape` says how to read it:

| `shape` | The set | Fields |
|---|---|---|
| `"interval"` | [lower, upper] | `lower`, `upper` |
| `"rays"` | (−∞, left_upper] ∪ [right_lower, +∞) | `left_upper`, `right_lower` |
| `"rays_left_only"`, `"rays_right_only"` | a single ray | `left_upper` or `right_lower` |
| `"line"` | the whole real line | |
| `"grid"` | a union of segments (several instruments) | `segments`: a list of (lo, hi, open_left, open_right); an open end reaches the edge of the search grid |
| `"empty"` | no value is accepted | |
| `"unavailable"` | the Conley variance is negative within the search range | |

The other fields are `pval_at_zero` and `stat_at_zero` (the AR test of 0), `contains_zero`, `cv` (the χ² critical value), `df_chi2`, `alpha` and `method` (`"fieller"` or `"grid"`).

### The bootstrap result

```python
boot = ivmediate.bootstrap(
    bw, Y="inctaxpc", T="f_prot", M="f_rw", Z="kmwittenberg",
    X="f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop",
    B=999, seed=20260317, progress=False,
)
ivmediate.print_bootstrap(boot, title="Becker and Woessmann (2009)")
```

```text
  Resolved X (11 vars): f_young, f_jew, f_fem, f_ortsgeb, f_pruss, f_blind, f_deaf, f_dumb, ... (+3 more)

======================================================================
  Becker and Woessmann (2009)  (method=pairs, B=999, seed=20260317)
======================================================================
  Bootstrap CI type: basic  (95%, α = 0.05)
  999 successful replicates.
  ────────────────────────────────────────────────────────────────────
  Parameter            Estimate    Boot SE                    Boot CI
  τ_total                0.5865     0.2592   [    0.0571,     1.0800] **
  τ (direct T→Y)        -0.2992     0.4009   [   -0.7490,     0.2334]
  θ (mediator→Y)         4.7091     4.4634   [   -0.8436,     9.2793]
  γ (T→mediator)         0.1881     0.0283   [    0.1274,     0.2349] ***
  θγ (indirect)          0.8857     0.5497   [   -0.0258,     1.7411]
  κ*                  -420.4053   213.8456   [ -827.8824,    -4.1637]
  ρ*_TY                 -0.2567     0.1261   [   -0.5056,    -0.0042]
  Sig. codes: *** p<0.01, ** p<0.05, * p<0.1
======================================================================
```

| Key | What it holds |
|---|---|
| `point` | The estimates on the original data: `total_hat`, `tau_hat`, `theta_hat`, `gamma_hat`, `indirect_hat`, `kappa_star`, `rho_TY_star`, and `RV_theta`, `XRV_theta` with `bootstrap_diagnostics=True`. |
| `boot_se` | Bootstrap standard errors (the standard deviation of the replicates), with the same keys. |
| `boot_ci` | Bootstrap CIs of type `ci_type`, as [lower, upper]. |
| `analytic_se`, `analytic_ci` | The analytic standard errors and CIs of the original fit, for comparison (`None` for κ*, ρ*_TY, RV and XRV). |
| `n_failed`, `n_est` | Failed replicates and the size of the resampled estimation sample. |
| `method`, `B`, `seed`, `ci_type`, `alpha` | The settings, for reproducibility. |
| `boot_dist` | The B replicate values of each parameter (only with `return_boot_dist=True`). |

The indirect effect is `indirect_hat` here and `indirect` in the result of `ivmediate()`.

## Differences from the other languages

- Python is the reference implementation: the R, Stata and MATLAB versions are checked against its results on a common set of validation runs ([validation](https://github.com/aferrara87/ivmediate/tree/main/validation)), to a relative tolerance of 1e-8.
- On Windows, glob patterns in `X`, `dummies` and `lasso_keep` ignore case, because Python's `fnmatch` follows the operating system. On macOS and Linux, and in the other languages, they are case-sensitive.
- `dummies` on a pandas float column (for example an integer code with missing values) gives indicator names like `_D_v_2.0` where R gives `_D_v_2`. The indicators are the same; only the names differ, which matters for `lasso_keep` and `expand_override`.
- Bootstrap draws come from Python's own random-number stream. With the same seed, the bootstrap results therefore differ from the other languages' by simulation noise, which can be large for a parameter with heavy tails: for θ in the Becker–Woessmann example, whose conditional first stage is only moderately strong, the bootstrap SE ranges from 2.4 to 4.5 across the four languages (B = 999), while the bootstrap CIs differ much less.
- The cross-validated LASSO penalty (`lasso_penalty="cv"`) exists only in Python; the other languages offer the plug-in penalty.
- The printed output is identical to R's.
- Speed: see the [comparison on the main page](https://github.com/aferrara87/ivmediate#speed).

## Help inside Python

```python
import ivmediate
help(ivmediate)                  # overview and conventions
help(ivmediate.ivmediate)        # all estimation options and the returned dict
help(ivmediate.print_results)
help(ivmediate.bootstrap)
help(ivmediate.print_bootstrap)
```

In IPython and Jupyter, `ivmediate.ivmediate?` shows the same text.
