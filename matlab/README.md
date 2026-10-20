# ivmediate for MATLAB

Causal mediation analysis with a single instrument (Ferrara, Gold, Heblich and Pinto 2026) in base
MATLAB, with the same results as the Python, R and Stata versions.

[Main page](../README.md) · [Python](../python/README.md) · [R](../r/README.md) · [Stata](../stata/README.md) · [MATLAB](../matlab/README.md)

## Installation

Download the repository from <https://github.com/aferrara87/ivmediate> (Code → Download ZIP) and unzip
it, or clone it with `git clone https://github.com/aferrara87/ivmediate.git`. Then add its `matlab`
folder to the MATLAB path:

```matlab
addpath('<repo>/matlab')   % <repo>: the folder you unzipped or cloned
savepath                   % optional: keep the path for later sessions
```

`savepath` writes MATLAB's `pathdef.m`. If it cannot (for example without write access to the MATLAB
installation folder), put the `addpath` line in your `startup.m` instead.

- **Requirements:** MATLAB R2016b or later, base MATLAB only: no toolboxes. Tested on R2026b.
- **GNU Octave:** not tested.
- **Version:** `ver` lists the package (`ivmediate: IV mediation analysis with a single ...`) with
  `Version 2.0.0`, and every result carries it: `res.version` is `'2.0.0'`. `which ivmediate` shows
  which copy MATLAB uses.
- **Check the installation:** `run('<repo>/matlab/tests/test_bw.m')` reproduces the example below and
  prints `ivmediate 2.0.0: Becker-Woessmann Table 6, column 1 reproduced (16 values)`;
  `runtests('<repo>/matlab/tests')` runs the same test in MATLAB's unit-test framework.

## Quick start

The Becker and Woessmann (2009) example of the main page: column 1 of Table 6 in the paper (Section 8).
[`examples/becker_woessmann/bw.m`](../examples/becker_woessmann/bw.m) runs the same example.

```matlab
data = ivmediate_read_csv('<repo>/examples/becker_woessmann/becker_woessmann_2009.csv');
res = ivmediate(data, 'inctaxpc', 'f_prot', 'f_rw', {'kmwittenberg'}, ...
    'X', 'f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop');
ivmediate_print(res, 'title', 'Becker and Woessmann (2009)');
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

The numbers are those of the main page, which explains the output block by block:
[Example: Becker and Woessmann (2009)](../README.md#example-becker-and-woessmann-2009). In
`matlab -batch` on Windows the console cannot show Greek letters, so the output then spells them out
(`theta`, `gamma`, `tau_total`, `kappa*`, ...); the option `'ascii'` (`true` or `false`) sets the
mode in any session.

## Syntax

```text
res  = ivmediate(data, Y, T, M, Z)
res  = ivmediate(data, Y, T, M, Z, Name, Value, ...)
ivmediate_print(res)
ivmediate_print(res, Name, Value, ...)
boot = ivmediate_bootstrap(data, Y, T, M, Z)
boot = ivmediate_bootstrap(data, Y, T, M, Z, Name, Value, ...)
ivmediate_print_bootstrap(boot)
ivmediate_print_bootstrap(boot, Name, Value, ...)
data = ivmediate_read_csv(path)
```

| Argument | Meaning |
|---|---|
| `data` | A table, a scalar struct of equal-length column vectors, or the output of `ivmediate_read_csv`. |
| `Y`, `T`, `M` | Names of the outcome, treatment and mediator columns. |
| `Z` | The instrument: a name, or a cell or string array of names for several instruments. |
| `res` | The results struct of `ivmediate` (see [Output](#output)). |
| `boot` | The result of `ivmediate_bootstrap`; `ivmediate_print(res, 'boot', boot)` puts its SEs and CIs in the table. |

- Names are text (char or string). Option names are case-insensitive; choices such as `'robust'` are
  case-sensitive. Logical options take `true`/`false` or `1`/`0`.
- `Y`, `T`, `M`, `Z`, the controls and the weights must be numeric; dummies, fixed effects and cluster
  ids may be numbers or text.
- Missing values are `NaN` in numeric columns and `''`, blanks, `<missing>` or `<undefined>` in text
  columns. Rows with a missing value in any column the model uses are dropped.

## Options

`[]` means "none" or, where stated, "automatic". The paper column gives where the paper discusses the
quantity.

### Model and controls

| Option | Default | What it does | Paper |
|---|---|---|---|
| `X` | `[]` | Controls: a cell or string array of names, or one space-separated string. The globs `*`, `?` and `[...]` expand against the column names (case-sensitive). A glob never picks up Y, T, M, Z or the fixed-effect, cluster, weight and coordinate columns, and naming Y, T, M or Z is an error. Collinear controls are dropped and listed in `res.dropped_collinear`. | Section 4.1 |
| `dummies` | `[]` | Categorical variables, expanded to indicators for every level but the first (in sorted order); the same input forms as `X`. The indicators are named `_D_<var>_<level>`. | |

### Fixed effects and weights

| Option | Default | What it does | Paper |
|---|---|---|---|
| `fe` | `[]` | One or two fixed-effect columns (a name, or a cell array of two names), absorbed. Singleton groups are dropped, repeatedly, until none is left. Two-way: an exact projection when the smaller dimension has at most 3,000 levels; otherwise alternating projections (Irons–Tuck), with a warning if they do not converge. The absorbed degrees of freedom are the rank of the fixed-effect dummies; with cluster SEs, a dimension nested in the clusters is not counted (as in reghdfe). | |
| `weights` | `[]` | A column of non-negative weights. Rows with zero weight are dropped. | |
| `wtype` | `'aweight'` | `'aweight'` or `'pweight'`. The two differ only with `se = 'classical'`: `'pweight'` then switches to robust SEs, with a note. | |

### Standard errors

| Option | Default | What it does | Paper |
|---|---|---|---|
| `se` | `'robust'` | `'classical'`, `'robust'` (HC1), `'cluster'` (one-way, CR1) or `'spatial'` (Conley 1999). The small-sample factors are n/(n − k) for HC1 and G/(G − 1) · (n − 1)/(n − k) for CR1, where k counts the regressors and the absorbed fixed effects (for CR1, without the dimensions nested in the clusters). | Remark 4; Appendix D.1 |
| `clusterid` | `[]` | The cluster column. Required with `se = 'cluster'`; ignored, with a note, otherwise. | |
| `lat`, `lon` | `[]` | Coordinate columns. Required with `se = 'spatial'`. | |
| `spatial_cutoff` | `[]` | The distance bandwidth. Required with `se = 'spatial'`; in km with the haversine metric. | |
| `spatial_kernel` | `'uniform'` | `'uniform'` or `'bartlett'`. | |
| `spatial_metric` | `'haversine'` | `'haversine'` (coordinates in degrees, distances in km) or `'euclidean'` (in the units of the coordinates). | |
| `alpha` | `0.05` | Level of the Wald CIs (normal critical values) and the AR sets: 0.05 gives 95% sets. | |

The spatial options (`lat`, `lon`, `spatial_cutoff`, `spatial_kernel`, `spatial_metric`) are an
error with any other `se`.

### Instruments and Anderson–Rubin sets

| Option | Default | What it does | Paper |
|---|---|---|---|
| `Z` (argument) | | One instrument, or several. With several, the diagnostics use the Olea–Pflueger effective F, and the R²-based sensitivity uses the first instrument. | Section 4.2 |
| `ar_ci` | `true` | Anderson–Rubin confidence sets for θ, γ and τ_total, valid however weak the instruments are. With one instrument they are exact (Fieller inversion); with several, the AR statistic is searched on a grid and the boundaries are refined. | Remark 3 (Section 4.2); Table 6 notes (Section 8) |
| `ar_grid_halfwidth` | `200` | With several instruments: half-width of the search grid, in standard errors around the estimate. | |
| `ar_grid_points` | `4001` | With several instruments: number of grid points (an integer of at least 3). | |

### Weak-instrument diagnostics

| Option | Default | What it does | Paper |
|---|---|---|---|
| `op_tau` | `0.10` | Olea–Pflueger worst-case bias tolerance for the critical values of the first stage (F_TZ) and the conditional first stage (F_MZ_T); any value in (0, 1). | Remark 3 (Section 4.2); Appendix F.4, step 2 |
| `op_alpha` | `0.05` | Size of the Olea–Pflueger test, in (0, 1). | |

The critical values come from the simplified Olea–Pflueger test at the estimated effective degrees of
freedom K_eff (Montiel Olea and Pflueger 2013); with one instrument and the defaults the critical
value is 23.11. `res.F_MZ_T_above_crit` records whether the conditional first stage reaches it.

### Sensitivity and printing

The κ summary (c_TM, κ\*, κ\*_CI and ρ\*_TY) is always printed, and every value below is in `res`
whatever is printed. The options of `ivmediate_print`:

| Option | Default | What it does | Paper |
|---|---|---|---|
| `kappa_table` | `false` | Also print the κ-adjusted estimates: θ and its CI with κ_0 at 0, 10, 25, 50, 75, 90 and 100% of κ\*. When `res.F_MZ_T_above_crit` is false, a note below the table says that c_TM is then estimated imprecisely. | Section 5.1.1, Proposition 5; Appendix F.4 |
| `r2_sensitivity` | `false` | Also print the R²-based sensitivity: the reduced-form F (= t²), the robustness values RV and XRV for θ, and their components. | Appendix F.2–F.3 |
| `boot` | `[]` | A result of `ivmediate_bootstrap`. The table then shows the bootstrap SEs and CIs, κ\*_CI uses the bootstrap SE of θ, and the bootstrap CIs of κ\* and ρ\*_TY are added. | |
| `title` | `'IV Mediation Analysis'` | The title line. | |
| `stars` | `true` | Significance stars, from the Wald p-values (with `boot`, from the bootstrap SEs). | |
| `star_cutoffs` | `[0.01 0.05 0.10]` | The p-value cutoffs for `***`, `**` and `*`. | |
| `ascii` | `[]` | Plain-ASCII names (`theta`, `gamma`, ...) instead of Greek letters. `[]`: on in `matlab -batch` on Windows, off otherwise. | |

`ivmediate_print_bootstrap` takes `title` (default `'Bootstrap Mediation'`), `stars`, `star_cutoffs`
and `ascii`. `ivmediate` and `ivmediate_bootstrap` take `ascii` for their notes, and `quiet`
(default `false`) to suppress their notes, warnings and progress lines.

### LASSO

| Option | Default | What it does | Paper |
|---|---|---|---|
| `lasso` | `false` | Post-double-selection LASSO (Belloni, Chernozhukov and Hansen 2014). The candidate controls are `X`, the `dummies` and the terms built from them (below); a LASSO for Y, T, M and each instrument selects among them, and the union of the selections enters the estimation. Not available in `ivmediate_bootstrap`. | Remark 6 (Section 4.2); Appendix H |
| `lasso_penalty` | `'plugin'` | The plug-in penalty. `'cv'` is an error in MATLAB (see [Differences](#differences-from-the-other-languages)). | Appendix H.2 |
| `expand_x` | `{'square', 'cube'}` | Terms built from each continuous control: `'square'`, `'cube'`, `'interact_with_X'`. | |
| `expand_dummies` | `{'interact_with_X'}` | Terms built from the dummy indicators: `'interact_with_X'`, `'interact_with_dummies'`. | |
| `expand_override` | `[]` | Terms per variable: a `containers.Map` (name → cell array of terms), an N×2 cell `{name, {terms}}`, or a struct. | |
| `lasso_keep` | `[]` | Candidates that are always kept (names or globs). | |
| `show_selection` | `false` | `ivmediate_print` then lists the controls selected in each equation. | Appendix H.4 |
| `lasso_c_slack` | `1.1` | The plug-in constant c (> 0). | |
| `lasso_alpha_p` | `[]` | The plug-in α_p, in (0, 1); `[]` = 0.1/log(n). | |
| `lasso_n_iter_max` | `15` | Iterations for the plug-in penalty loadings. | |
| `lasso_tol` | `1e-6` | Convergence tolerance of those iterations. | |
| `lasso_lambda_override` | `[]` | A fixed penalty level λ (≥ 0) instead of the plug-in. | |
| `lasso_verbose` | `false` | Print the plug-in iterations. | |
| `lasso_seed` | `20260317` | Seed for cross-validation; unused by the plug-in penalty. | |
| `lasso_cv_folds` | `5` | Cross-validation folds; checked, but cross-validation is not available in MATLAB. | |

### Bootstrap

`ivmediate_bootstrap` takes the model options of `ivmediate` (`X`, `dummies`, `fe`, `se`,
`clusterid`, `alpha`, `weights`, `wtype`, `op_tau`, `op_alpha` and the spatial options; not the
LASSO), and these:

| Option | Default | What it does |
|---|---|---|
| `method` | `[]` | `'pairs'` (resample observations) or `'block'` (resample clusters). `[]`: `'block'` with `se = 'cluster'`, `'pairs'` otherwise. |
| `B` | `999` | Number of replications, an integer of at least 50. |
| `seed` | `20260317` | A non-negative integer; the same seed gives the same draws. |
| `ci_type` | `'basic'` | `'basic'`, `'percentile'` or `'studentized'`. |
| `bootstrap_diagnostics` | `false` | Also bootstrap RV and XRV. |
| `return_boot_dist` | `false` | Add the replicate estimates (`boot.boot_dist`). |
| `progress` | `true` | Print progress lines. |
| `verbose` | `false` | Report failed replicates. |
| `quiet`, `ascii` | `false`, `[]` | As in `ivmediate`. |
| `boot_indices`, `return_indices` | `[]`, `false` | Testing hooks: replay given draws; return the draws. |

The bootstrap resamples the estimation sample: rows with missing values, rows with zero weight and,
for `'block'`, rows without a cluster id are removed first. Each drawn cluster counts as a new cluster,
and fixed effects nested in the clusters get one level per drawn copy. The AR sets are not
bootstrapped. With `'basic'` or `'percentile'` intervals (and no `bootstrap_diagnostics`), the
replicates compute point estimates only, which is faster.

## Output

### The printout

- **Header:** n, the number of instruments and covariates, the SE type (and the number of clusters).
  When they apply, lines for the weights, the dummy indicators, dropped collinear controls, the fixed
  effects (groups, singletons dropped, absorbed degrees of freedom, dimensions nested in the clusters),
  the Conley settings and the LASSO selection.
- **Mediation estimates:** τ_total, τ, θ, γ and θγ, with SEs, CIs and stars (bootstrap SEs and CIs with
  `'boot'`). A warning appears if τ + θγ differs from τ_total (Corollary 1); with several
  instruments the two agree only asymptotically, and a note appears when they differ by more than 1%.
- **First-stage diagnostics:** π_1 (T on Z given X), the first-stage F_{T, Z|X} and the conditional
  F_{M, Z|T, X}, each with its Olea–Pflueger critical value (K_eff, τ, α). With several instruments,
  also the joint Wald F, the effective F and, for the conditional first stage, each instrument's t².
- **Anderson–Rubin sets** for θ, γ and τ_total, each with the AR p-value at zero. Unbounded sets are
  marked as the signature of weak instruments.
- **κ-based sensitivity:** c_TM, κ\* (the breakeven), κ\*_CI (flagged when the CI for θ already
  contains zero) and ρ\*_TY.
- **On request:** the κ-adjusted table (`kappa_table`), the R² section (`r2_sensitivity`) and the LASSO
  selection (`show_selection`).

### The results struct

| Field | Content |
|---|---|
| `theta_hat`, `se_theta`, `ci_theta` | The estimate of θ, its SE and its Wald CI (1×2); likewise `tau_*`, `gamma_*` and `total_*`, and `indirect`, `se_indirect`, `ci_indirect` for θγ. |
| `n`, `n_instruments`, `n_covariates`, `n_dummies`, `n_clusters` | Sizes; `n_clusters` is `[]` without cluster SEs. |
| `se_type`, `alpha`, `weighted`, `wtype` | The settings used (`se_type` after a `'pweight'` switch). |
| `pi1`, `se_pi1` | The first-stage coefficients of Z for T (1×K). |
| `F_TZ`, `F_MZ_T` | The Wald F of the instruments in the first stage and in the conditional first stage (joint with several instruments; `F_MZ_T_individual` has each instrument's t²). |
| `F_TZ_op`, `F_MZ_T_op`, `K_eff_TZ`, `K_eff_MZ_T`, `F_TZ_op_critval`, `F_MZ_T_op_critval` | The Olea–Pflueger effective F, the effective degrees of freedom and the critical values. |
| `F_MZ_T_above_crit` | `true` if `F_MZ_T_op` reaches `F_MZ_T_op_critval` (at `op_tau`), `false` if not, `[]` if either is not finite. |
| `ar_theta`, `ar_gamma`, `ar_total` | The AR sets (structs, below); `[]` with `'ar_ci', false`. |
| `c_TM_hat`, `kappa_star`, `kappa_star_CI`, `rho_TY_star` | c_TM, κ\*, κ\*_CI and ρ\*_TY. |
| `kappa_star_CI_informative` | `true` if the Wald CI for θ excludes zero, so that κ\*_CI is informative; `false` if not; `[]` if the SE of θ is not finite. |
| `kappa_adj_table` | The κ table as a struct array with fields `kappa`, `kappa_frac`, `theta_adj`, `ci_lo` and `ci_hi`, one element per fraction of κ\*. |
| `F_reduced_form`, `RV_theta`, `XRV_theta`, `rv_rf`, `rv_fs`, `xrv_rf`, `xrv_fs` | The R²-based sensitivity. |
| `fe`, `n_fe_groups`, `n_dropped_singletons`, `df_fe_absorbed` | The fixed effects: names, groups (summed over dimensions), singletons dropped, absorbed degrees of freedom. |
| `fe_nested_in_cluster` | With cluster SEs, one logical per fixed-effect dimension: `true` if every group lies within one cluster (the dimension is then not counted in the CR1 degrees of freedom). `[]` without cluster SEs or without fixed effects. |
| `dropped_collinear` | Controls dropped as collinear. |
| `lasso_info`, `rich_col_names` | The LASSO selection (`lasso_info.union_names`, `n_selected_union`, `n_total`, and per equation `details`) and the candidate names. |
| `eps_T_hat`, `eps_M_hat`, `eps_Y_hat` | The residuals behind the κ analysis. |
| `version` | `'2.0.0'`. |

Each AR set is a struct:

| Field | Content |
|---|---|
| `shape` | `'interval'` [`lower`, `upper`]; `'rays'` (−∞, `left_upper`] ∪ [`right_lower`, ∞); `'rays_left_only'` (−∞, `left_upper`]; `'rays_right_only'` [`right_lower`, ∞); `'line'`, the whole real line; `'empty'`; `'grid'` with several instruments (see `segments`); `'unavailable'` when a Conley variance is not positive. |
| `lower`, `upper`, `left_upper`, `right_lower` | The bounds; `NaN` where the shape has none. |
| `segments` | With several instruments, a struct array with fields `lo`, `hi`, `open_left` and `open_right` (`true` when the segment reaches the edge of the search grid). |
| `contains_zero`, `stat_at_zero`, `pval_at_zero` | The AR test of the value zero. |
| `cv`, `df_chi2`, `alpha`, `method` | The χ² critical value, its degrees of freedom (the number of instruments), the level, and `'fieller'` or `'grid'`. |

### The bootstrap result

```matlab
boot = ivmediate_bootstrap(data, 'inctaxpc', 'f_prot', 'f_rw', {'kmwittenberg'}, ...
    'X', 'f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop', ...
    'progress', false);
ivmediate_print_bootstrap(boot, 'title', 'Becker and Woessmann (2009): bootstrap');
```

```text
  Resolved X (11 vars): f_young, f_jew, f_fem, f_ortsgeb, f_pruss, f_blind, f_deaf, f_dumb, ... (+3 more)

======================================================================
  Becker and Woessmann (2009): bootstrap  (method=pairs, B=999, seed=20260317)
======================================================================
  Bootstrap CI type: basic  (95%, α = 0.05)
  999 successful replicates.
  ────────────────────────────────────────────────────────────────────
  Parameter            Estimate    Boot SE                    Boot CI
  τ_total                0.5865     0.2483   [    0.0491,     1.0210] **
  τ (direct T→Y)        -0.2992     0.2409   [   -0.7246,     0.1832]
  θ (mediator→Y)         4.7091     2.4320   [   -0.2680,     8.9602] *
  γ (T→mediator)         0.1881     0.0268   [    0.1335,     0.2353] ***
  θγ (indirect)          0.8857     0.4379   [    0.0083,     1.6934] **
  κ*                  -420.4053   205.6373   [ -806.9583,   -10.3533]
  ρ*_TY                 -0.2567     0.1190   [   -0.4930,    -0.0225]
  Sig. codes: *** p<0.01, ** p<0.05, * p<0.1
======================================================================
```

The bootstrap results differ from those of the other languages, which draw different samples
(see [Differences](#differences-from-the-other-languages)). `boot` holds:

| Field | Content |
|---|---|
| `point`, `boot_se`, `boot_ci` | Structs with the fields `total_hat`, `tau_hat`, `theta_hat`, `gamma_hat`, `indirect_hat`, `kappa_star` and `rho_TY_star` (and `RV_theta`, `XRV_theta` with `bootstrap_diagnostics`): the estimates, the bootstrap SEs and the bootstrap CIs (1×2). |
| `analytic_se`, `analytic_ci` | The analytic SEs and CIs of the same quantities (`[]` where there are none). |
| `n_failed`, `n_est` | Failed replicates, and the size of the estimation sample. |
| `method`, `B`, `seed`, `ci_type`, `alpha` | The settings. |
| `boot_dist` | With `return_boot_dist`: the replicate estimates (B×1 per quantity). |

## Differences from the other languages

- **Data input.** `data` is a table, a struct of columns, or the result of `ivmediate_read_csv`. That
  reader treats empty fields and the tokens `NA`, `N/A`, `NULL`, `NaN`, `nan` and `#N/A` as missing in
  every column. A table from `readtable` works as well; for a file that codes missing values as `NA`,
  use `readtable(file, 'TreatAsMissing', 'NA')`.
- **ASCII output in `matlab -batch`** on Windows, as in the [Quick start](#quick-start); the option
  `'ascii'` overrides it.
- **Python's `None` is `[]` here**, for example `n_clusters` without cluster SEs, the AR sets with
  `'ar_ci', false`, and `F_MZ_T_above_crit` when the F or its critical value is not finite.
- **Bootstrap draws.** The package draws from its own generator (splitmix64, keyed by the seed and the
  replicate number), not from MATLAB's `rng`, so your random-number state is left untouched and the
  same seed gives the same draws on every machine. With the same seed, the bootstrap results therefore
  differ from the other languages' by simulation noise, which can be large for a parameter with heavy
  tails: for θ in the Becker–Woessmann example, whose conditional first stage is only moderately strong,
  the bootstrap SE ranges from 2.4 to 4.5 across the four languages (B = 999), while the bootstrap CIs
  differ much less. The point estimates and the analytic results agree.
- **Not available in MATLAB:** `lasso_penalty = 'cv'` (cross-validation), which only the Python package
  offers.
- **Small conventions:** glob matching is case-sensitive (Python's is not on Windows), and integer
  levels of a numeric dummy are labeled `_D_v_2` (pandas: `_D_v_2.0`).
- **Speed.** MATLAB uses Intel's multithreaded linear algebra on all cores by default. From N = 10,000
  on it is the fastest of the four in the [speed table](../README.md#speed); at N = 1,000 every
  language takes about a tenth of a second.

The estimator itself is described on the main page: [The estimator](../README.md#the-estimator).

## Help inside MATLAB

```matlab
help ivmediate
help ivmediate_print
help ivmediate_bootstrap
help ivmediate_print_bootstrap
help ivmediate_read_csv
```

Each help text gives the syntax and every option with its default; `doc ivmediate` shows the same text
in the Help browser.
