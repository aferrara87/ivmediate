# ivmediate for R

The R package `ivmediate` (version 2.0.0): causal mediation analysis with a single instrument.

[Main page](../README.md) · [Python](../python/README.md) · [R](../r/README.md) · [Stata](../stata/README.md) · [MATLAB](../matlab/README.md)

## Installation

```r
# install.packages("remotes")
remotes::install_github("aferrara87/ivmediate", subdir = "r")
```

From a clone of the repository, run `R CMD INSTALL r` in its root folder.

**Requirements:**
- R ≥ 4.2;
- the package Matrix, which ships with R;
- parallel, stats and utils, which are part of base R.

`remotes` is needed only for `install_github()`. On Windows it may print a message asking you to install
Rtools. The package has no compiled code and installs without it.

**Version.** Check it with:

```r
packageVersion("ivmediate")
#> [1] '2.0.0'
```

Every result also stores it, in `res$version`.

**Tip: linear algebra.** R on Windows uses its built-in, single-threaded reference BLAS by default. The
[speed table](../README.md#speed) was run with it, and R falls behind Python and MATLAB there for large
models; they use Intel's multithreaded MKL. An optimized BLAS can speed up large models. To see which BLAS your R uses, run `sessionInfo()`. On Windows,
"Matrix products: default" with no BLAS line means the reference BLAS. How to install another BLAS depends
on the platform; see the section on linear algebra in the *R Installation and Administration* manual.

## Quick start

The Becker and Woessmann (2009) data ship with the package. This reproduces column 1 of Table 6 in the
paper (Section 8):

```r
library(ivmediate)

bw <- read.csv(system.file("extdata", "becker_woessmann_2009.csv", package = "ivmediate"))
res <- ivmediate(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw", Z = "kmwittenberg",
                 X = "f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop")
print(res, title = "Becker and Woessmann (2009)")
```

```
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

The main page explains these numbers: [Example: Becker and Woessmann (2009)](../README.md#example-becker-and-woessmann-2009).
The estimator itself is described in [The estimator](../README.md#the-estimator). The same data and a
script, `bw.R`, are in [`examples/becker_woessmann/`](../examples/becker_woessmann/).

## Syntax

```r
res <- ivmediate(data, Y, T, M, Z, X = NULL, dummies = NULL, fe = NULL,
                 se = "robust", clusterid = NULL, alpha = 0.05,
                 weights = NULL, wtype = "aweight",
                 op_tau = 0.10, op_alpha = 0.05,
                 ar_ci = TRUE, ar_grid_halfwidth = 200, ar_grid_points = 4001,
                 lat = NULL, lon = NULL, spatial_cutoff = NULL,
                 spatial_kernel = NULL, spatial_metric = NULL,
                 lasso = FALSE, lasso_penalty = "plugin", lasso_keep = NULL,
                 lasso_seed = 20260317, lasso_cv_folds = 5,
                 expand_x = c("square", "cube"), expand_dummies = "interact_with_X",
                 expand_override = NULL, show_selection = FALSE,
                 lasso_c_slack = 1.1, lasso_alpha_p = NULL, lasso_n_iter_max = 15,
                 lasso_tol = 1e-6, lasso_lambda_override = NULL, lasso_verbose = FALSE)

print(res, title = "IV Mediation Analysis", boot = NULL, stars = TRUE,
      star_cutoffs = c(0.01, 0.05, 0.10), kappa_table = FALSE, r2_sensitivity = FALSE)

boot <- ivmediate_bootstrap(data, Y, T, M, Z, X = NULL, dummies = NULL, fe = NULL,
                            se = "robust", clusterid = NULL, alpha = 0.05,
                            weights = NULL, wtype = "aweight", op_tau = 0.10, op_alpha = 0.05,
                            lat = NULL, lon = NULL, spatial_cutoff = NULL,
                            spatial_kernel = NULL, spatial_metric = NULL,
                            method = NULL, B = 999, seed = 20260317, ci_type = "basic",
                            return_boot_dist = FALSE, bootstrap_diagnostics = FALSE,
                            progress = TRUE, verbose = FALSE)

print(boot, title = "Bootstrap Mediation", stars = TRUE, star_cutoffs = c(0.01, 0.05, 0.10))
```

- `ivmediate()` returns a list of class `"ivmediate"`. Typing `res` prints it with the default title.
- `print(res, boot = boot)` puts the bootstrap SEs and CIs into the main table.
- `ivmediate_bootstrap()` takes the model arguments of `ivmediate()` except the AR-set and LASSO options.
- Logical options must be `TRUE` or `FALSE`; text such as `"false"` is an error.

## Options

Paper references are to sections, remarks, propositions and appendices of Ferrara, Gold, Heblich and Pinto
(2026).

### Model and controls

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `data` | required | A data frame. data.table and tibble inputs are accepted. It is never modified. | |
| `Y`, `T`, `M` | required | Names of the outcome, treatment and mediator columns. | Section 2 |
| `Z` | required | Name(s) of the instrument(s): one or several, different from `Y`, `T` and `M`, without repeats. | Section 3.3, Proposition 2 |
| `X` | `NULL` | Controls: a character vector or one whitespace-separated string (several lines are fine). Shell-style globs (`*`, `?`, `[...]`) are matched against `names(data)`. Globs never match `Y`, `T`, `M`, `Z` or the fixed-effect, cluster, weight and coordinate columns, and a note lists what they skipped. Naming `Y`, `T`, `M` or `Z` in `X` is an error. Exactly collinear controls are dropped, with a note. | |
| `dummies` | `NULL` | Categorical variables, expanded into indicators with the first category dropped: the first in sorted order for text columns, the first level for factors. Same parsing and protection as `X`. Empty or whitespace-only text counts as missing. | |

Rows with a missing value in any column the model uses are dropped. A factor or character column in a
numeric role is converted through its labels when every label is a number, with a note; otherwise it is
an error.

### Fixed effects and weights

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `fe` | `NULL` | One or two fixed-effect variables, absorbed by demeaning. For two-way FE the projection is exact when the FE with fewer levels has at most 3,000 levels; otherwise it uses alternating projections (Irons–Tuck acceleration, relative tolerance 1e-12). Singleton groups are dropped iteratively. The absorbed degrees of freedom are the rank of the FE dummies; with cluster SEs, FE nested in the clusters are not counted in the small-sample factor. | |
| `weights` | `NULL` | Name of the weight column. Negative weights are an error; zero-weight rows are dropped, with a note. Fixed effects use weighted group means. | |
| `wtype` | `"aweight"` | `"aweight"` (analytic) or `"pweight"` (probability). With `se = "classical"`, `"pweight"` switches to robust SEs, with a note. | |

### Standard errors

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `se` | `"robust"` | `"classical"`, `"robust"` (HC1), `"cluster"` (CR1) or `"spatial"` (Conley). | |
| `clusterid` | `NULL` | Cluster variable. It is required with `se = "cluster"` and ignored otherwise, with a note. Rows without a cluster id are dropped, and at least 2 clusters are needed. | |
| `lat`, `lon` | `NULL` | Coordinate columns for `se = "spatial"`: degrees for `"haversine"`, projected units for `"euclidean"`. Giving them with another `se` is an error. | |
| `spatial_cutoff` | `NULL` | Kernel cutoff, required with `se = "spatial"`: km with `"haversine"`, coordinate units with `"euclidean"`. A uniform cutoff that covers every pair is an error. A neighbor share above 50% prints a warning. Conley SEs depend on the cutoff, so report several. | |
| `spatial_kernel` | `NULL` (= `"uniform"`) | `"uniform"` (weight 1 within the cutoff) or `"bartlett"` (linear taper). | |
| `spatial_metric` | `NULL` (= `"haversine"`) | `"haversine"` (great-circle km, Earth radius 6371 km) or `"euclidean"`. | |
| `alpha` | `0.05` | Level for the confidence intervals, the AR sets and the RV/XRV threshold. | |

### Instruments and Anderson–Rubin sets

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `ar_ci` | `TRUE` | Anderson–Rubin confidence sets for θ, γ and τ_total, with the same SE type and `alpha` as the Wald intervals. They stay valid when the instruments are weak. With one instrument, each set is computed exactly (Fieller): an interval, two rays or the whole line. With several, the joint χ²(K) statistic is solved on a grid, with exact boundaries. No set is reported for τ. | Remark 3 (Section 4.2); Table 6 notes (Section 8) |
| `ar_grid_halfwidth` | `200` | Half-width of the AR search grid, in analytic SEs of the parameter. | |
| `ar_grid_points` | `4001` | Number of grid points (an integer of at least 3). | |

### Weak-instrument diagnostics

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `op_tau` | `0.10` | Worst-case bias tolerance for the Olea–Pflueger critical values of both first-stage F statistics: 0.10 is a bias of 10% of the OLS bias. | Remark 3; Appendix F.4, step 2 |
| `op_alpha` | `0.05` | Size of the Olea–Pflueger test. Critical values come from the simplified test at the estimated effective degrees of freedom K_eff (K_eff = 1 with one instrument). | Remark 3 |

### Sensitivity and printing

These are arguments of `print(res, ...)`, except `show_selection`. The values behind every section are
always in `res`.

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `kappa_table` | `FALSE` | Also print the κ-adjusted estimates: θ and its CI with κ_0 at 0, 10, 25, 50, 75, 90 and 100% of the breakeven κ*. If the conditional first stage is below the Olea–Pflueger critical value, a note follows the table. | Section 5.1.1, Proposition 5; Appendix F.1 and F.4 |
| `r2_sensitivity` | `FALSE` | Also print the R²-based section: the reduced-form F (= t²), RV(θ) and XRV(θ), with their components. | Appendix F.2–F.3 |
| `title` | `"IV Mediation Analysis"` | Title line. | |
| `boot` | `NULL` | A result of `ivmediate_bootstrap()` for the same model. The table then shows bootstrap SEs and CIs, and the κ section adds bootstrap CIs of κ* and ρ*_TY. | |
| `stars` | `TRUE` | Significance stars, from two-sided normal p-values of estimate/SE. | |
| `star_cutoffs` | `c(0.01, 0.05, 0.10)` | Three p-value cutoffs for `***`, `**` and `*`, in any order. | |
| `show_selection` | `FALSE` | An argument of `ivmediate()`: print the LASSO selections (variables, penalty, convergence). | |

The κ section is always printed: c_TM, κ*, κ*_CI and ρ*_TY (Section 5.1.1; Remark 8 for ρ*_TY). The
main page shows the two optional sections for the BW example; in R:
`print(res, kappa_table = TRUE, r2_sensitivity = TRUE)`.

### LASSO

| Argument | Default | What it does | Paper |
|---|---|---|---|
| `lasso` | `FALSE` | Post-double-selection LASSO before the 2SLS. One plug-in LASSO is run for each of `Y`, `T`, `M` and every instrument, on an expansion of `X` and `dummies`; the union of the selections becomes the control set. Each LASSO has an unpenalized intercept. | Remark 6; Appendix H (Algorithm 1 in H.2) |
| `lasso_penalty` | `"plugin"` | The iterative plug-in penalty with post-LASSO loadings. `"cv"` is not available in R. | Appendix H.2 |
| `lasso_keep` | `NULL` | Names that always enter the control set: originals (`"x1"`), powers (`"x1^2"`) or interactions, with globs allowed. They are partialed out before each LASSO. | |
| `expand_x` | `c("square", "cube")` | Candidates from `X`: `"square"` and `"cube"` add powers of continuous variables, and `"interact_with_X"` adds products of pairs. | |
| `expand_dummies` | `"interact_with_X"` | Candidates from `dummies`: `"interact_with_X"` adds products with each `X` variable, and `"interact_with_dummies"` adds products of indicators from different variables. | |
| `expand_override` | `NULL` | Per-variable operations, as a named list, for example `list(age = "square")`. | |
| `lasso_c_slack` | `1.1` | Constant c of the penalty λ = 2c·√n·Φ⁻¹(1 − α_p/(2p)), with n the number of observations (Kish's effective n with weights) and p the number of candidates. | Appendix H.2 |
| `lasso_alpha_p` | `NULL` (= 0.1/log n) | Tail probability α_p of the penalty, with n the number of observations (also with clusters). | Appendix H.2 |
| `lasso_n_iter_max` | `15` | Maximum number of penalty-loading iterations. | |
| `lasso_tol` | `1e-6` | Convergence tolerance on the largest coefficient change. | |
| `lasso_lambda_override` | `NULL` | A fixed λ instead of the plug-in rule. | |
| `lasso_verbose` | `FALSE` | Print the iterations of each LASSO. | |
| `lasso_seed`, `lasso_cv_folds` | `20260317`, `5` | Used only by the cross-validated penalty, which R does not have. Accepted for compatibility. | |

### Bootstrap

`ivmediate_bootstrap()` refits the model on resampled data. Resampling starts from the estimation sample.
Rows the analysis would drop are removed first: those with missing values, zero weights or no cluster id,
and with Conley SEs those without coordinates.

| Argument | Default | What it does |
|---|---|---|
| `method` | `NULL` | `"pairs"` (resample observations) or `"block"` (resample clusters of `clusterid`). `NULL` means `"block"` with `se = "cluster"` and `"pairs"` otherwise. In the block bootstrap each drawn copy of a cluster is its own cluster, and FE nested in the clusters get new levels per copy. |
| `B` | `999` | Number of replicates, an integer of at least 50. |
| `seed` | `20260317` | Seed of the replicates' random-number streams. |
| `ci_type` | `"basic"` | `"basic"`, `"percentile"` or `"studentized"`. The studentized intervals use each replicate's analytic SE, and fall back to `"basic"` for κ*, ρ*_TY, RV and XRV, which have none. |
| `return_boot_dist` | `FALSE` | Also return all replicates (`boot_dist`). |
| `bootstrap_diagnostics` | `FALSE` | Also bootstrap RV(θ) and XRV(θ). |
| `progress` | `TRUE` | Print a progress line every 5% of the replicates. |
| `verbose` | `FALSE` | Print extra information. |

With `"basic"` or `"percentile"` intervals and no diagnostics, the replicates compute point estimates only,
which is faster. Three arguments are internal: `.point_only` in `ivmediate()`, which the replicates use,
and `boot_indices` and `.return_indices`, which the validation harness uses to replay the same draws in
every language.

## Output

### The printout

- **Header:** n, the number of instruments and covariates, the SE type (with the number of clusters), and
  when they apply, the weights, the dummy indicators, dropped collinear controls, the FE (groups,
  singletons dropped, absorbed df), the Conley settings and the LASSO selection count.
- **Mediation estimates:** τ_total, τ, θ, γ and θγ, with SEs, Wald CIs (normal critical values) and
  stars. The SE of θγ uses the delta method (Remark 4; Appendix D.1). With one instrument,
  τ + θγ = τ_total exactly (Corollary 1), and a warning appears only if that fails. With several
  instruments, a note appears when the two differ by more than 1%. They are then estimated equation by
  equation and agree only asymptotically.
- **First-stage diagnostics:** π_1, the coefficient of `Z` in the first stage, with its SE. Then
  F_{T,Z|X} for the usual first stage and F_{M,Z|T,X} for the conditional one, each with its
  Olea–Pflueger critical value, K_eff, τ and α (Remark 3). With several instruments there are separate
  lines for the joint Wald F and the effective F, and the conditional first stage also shows each
  instrument's t².
- **Anderson–Rubin sets:** for θ, γ and τ_total, the set and the AR p-value at zero. A note appears when a
  set is unbounded or disconnected.
- **κ-based sensitivity:**
  - c_TM, the covariance of the treatment and mediator errors;
  - κ*, the breakeven;
  - κ*_CI, marked "(uninformative)" when the CI for θ already contains zero;
  - ρ*_TY, the breakeven as a correlation, with a note when its absolute value exceeds 1;
  - with `boot`, the bootstrap CIs of κ* and ρ*_TY.
- **On request:** the κ-adjusted table with its note, the R²-based section, and the LASSO selections.

### The returned list

`ivmediate()` returns a list of class `c("ivmediate", "list")`. The elements you are most likely to need:

| Element | Content |
|---|---|
| `theta_hat`, `se_theta`, `ci_theta` | θ, its SE and Wald CI. The same holds for `tau_*`, `gamma_*`, `total_*` and for `indirect`, `se_indirect`, `ci_indirect` (θγ). |
| `cov_gamma_theta` | The covariance of the γ and θ estimates, used for the SE of θγ. |
| `pi1`, `se_pi1` | First-stage coefficient(s) of `Z` and SE(s); vectors with several instruments. |
| `F_TZ`, `F_MZ_T` | Wald F statistics of the two first stages (t² with one instrument). |
| `F_TZ_op`, `F_MZ_T_op`, `F_TZ_op_critval`, `F_MZ_T_op_critval`, `K_eff_TZ`, `K_eff_MZ_T` | Olea–Pflueger effective F statistics, critical values and effective degrees of freedom. |
| `F_MZ_T_above_crit` | `TRUE` if `F_MZ_T_op` reaches its critical value, `FALSE` if not, `NA` if either is not finite. |
| `ar_theta`, `ar_gamma`, `ar_total` | The AR sets, described below (`NULL` with `ar_ci = FALSE`). |
| `c_TM_hat`, `kappa_star`, `kappa_star_CI`, `rho_TY_star` | The κ quantities. |
| `kappa_star_CI_informative` | `TRUE` if the Wald CI for θ excludes zero, so that κ*_CI is informative. |
| `kappa_adj_table` | The κ table: a list of rows with `kappa`, `kappa_frac`, `theta_adj`, `ci_lo` and `ci_hi`. |
| `F_reduced_form`, `F_first_stage_cond`, `RV_theta`, `XRV_theta`, `rv_rf`, `rv_fs`, `xrv_rf`, `xrv_fs` | The R²-based quantities. |
| `n`, `n_clusters`, `n_instruments`, `n_covariates`, `n_dummies`, `se_type`, `alpha` | Sample and settings. |
| `fe`, `n_fe_groups`, `n_dropped_singletons`, `df_fe_absorbed` | Fixed effects. |
| `fe_nested_in_cluster` | With `se = "cluster"` and FE, a logical per FE variable: `TRUE` if the FE is nested in the clusters and therefore not counted in the CR1 factor. `NULL` otherwise. |
| `dropped_collinear` | Names of the controls dropped as collinear. |
| `lasso_info` | With `lasso = TRUE`: `union_names`, `n_selected_union`, `n_total`, and `details` per LASSO (`names`, `coefs`, `lambda_value`, `n_selected`, `converged`, `n_iter`). |
| `eps_T_hat`, `eps_M_hat`, `eps_Y_hat` | Residuals of the treatment, mediator and outcome equations. |
| `version` | The package version. |

**AR sets.** Each of `ar_theta`, `ar_gamma` and `ar_total` is a list:

| Field | Content |
|---|---|
| `shape` | `"interval"`, `"rays"`, `"rays_left_only"`, `"rays_right_only"`, `"line"`, `"empty"`, `"grid"` (several instruments) or `"unavailable"` (a negative Conley variance). |
| `lower`, `upper` | Bounds of an `"interval"`. |
| `left_upper`, `right_lower` | For two rays, (−∞, `left_upper`] ∪ [`right_lower`, ∞); one-sided rays use one of the two. |
| `segments` | For `"grid"`, a list of segments with `lo`, `hi`, `open_left` and `open_right`. |
| `pval_at_zero`, `stat_at_zero`, `contains_zero` | The AR test of the value 0. |
| `cv`, `df_chi2`, `alpha`, `method` | Critical value, χ² degrees of freedom, level, and `"fieller"` or `"grid"`. |

For the BW example, `res$ar_theta$shape` is `"interval"`, with `lower` 0.431 and `upper` 10.2.

### The bootstrap

```r
boot <- ivmediate_bootstrap(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw", Z = "kmwittenberg",
                            X = "f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop",
                            B = 999, progress = FALSE)
print(boot, title = "Becker and Woessmann (2009): bootstrap")
```

```
  Resolved X (11 vars): f_young, f_jew, f_fem, f_ortsgeb, f_pruss, f_blind, f_deaf, f_dumb, ... (+3 more)

======================================================================
  Becker and Woessmann (2009): bootstrap  (method=pairs, B=999, seed=20260317)
======================================================================
  Bootstrap CI type: basic  (95%, α = 0.05)
  999 successful replicates.
  ────────────────────────────────────────────────────────────────────
  Parameter            Estimate    Boot SE                    Boot CI
  τ_total                0.5865     0.2389   [    0.1144,     1.0572] **
  τ (direct T→Y)        -0.2992     0.2840   [   -0.7202,     0.2676]
  θ (mediator→Y)         4.7091     2.8339   [   -0.9453,     9.0492] *
  γ (T→mediator)         0.1881     0.0275   [    0.1305,     0.2380] ***
  θγ (indirect)          0.8857     0.4542   [    0.0284,     1.6970] *
  κ*                  -420.4053   201.2135   [ -806.5613,   -24.8422]
  ρ*_TY                 -0.2567     0.1182   [   -0.4926,    -0.0360]
  Sig. codes: *** p<0.01, ** p<0.05, * p<0.1
======================================================================
```

`ivmediate_bootstrap()` returns a list of class `c("ivmediate_bootstrap", "list")` with these elements:
- `point`, the original estimates;
- `boot_se`, the standard deviations of the replicates;
- `boot_ci`, the bootstrap CIs;
- `analytic_se` and `analytic_ci`, from the original fit;
- `n_failed`, the number of failed replicates;
- `n_est`, the size of the resampled estimation sample;
- `method`, `B`, `seed`, `ci_type` and `alpha`;
- with `return_boot_dist = TRUE`, `boot_dist`.

The parameters inside these elements are named `total_hat`, `tau_hat`, `theta_hat`, `gamma_hat`,
`indirect_hat`, `kappa_star` and `rho_TY_star`. With `bootstrap_diagnostics = TRUE`, `RV_theta` and
`XRV_theta` are added. Note that the indirect effect is `indirect_hat` here but `indirect` in the result of
`ivmediate()`.

## Differences from the other languages

- **Printing.** The printout is the same as Python's, character for character; it is checked on the BW
  example.
- **Globs.** Matching is case-sensitive in R. Python's matching is case-insensitive on Windows.
- **Dummy labels.** If pandas stores a numeric category column as floats (for example because it has
  missing values), its indicators are named like `_D_v_2.0` in Python and `_D_v_2` in R. This matters only
  when you name them in `lasso_keep`.
- **Bootstrap draws.** R draws from its own L'Ecuyer-CMRG random-number streams, one stream per
  replicate from `seed` (`parallel::nextRNGStream`). With the same seed, the bootstrap results therefore
  differ from the other languages' by simulation noise, which can be large for a parameter with heavy
  tails: for θ in the Becker–Woessmann example, whose conditional first stage is only moderately strong,
  the bootstrap SE ranges from 2.4 to 4.5 across the four languages (B = 999), while the bootstrap CIs
  differ much less. The point estimates agree. Your random-number state (`.Random.seed` and `RNGkind()`)
  is restored afterward.
- **Cross-validated LASSO penalty.** `lasso_penalty = "cv"` is not available in R. Use the default
  plug-in penalty, or fix λ with `lasso_lambda_override`.
- **R data types.** data.table and tibble inputs are accepted. Factor and character columns in a numeric
  role are converted when their labels are numbers. Factor dummies drop their first level, rather than the
  first label in sorted order.
- **Speed.** See the [speed table](../README.md#speed) and the BLAS tip under Installation.

## Help inside R

```r
?ivmediate                      # the analysis and all its arguments
?print.ivmediate                # the printout and its options
?ivmediate_bootstrap            # bootstrap inference
?print.ivmediate_bootstrap
package?ivmediate               # overview, with the citation
example(ivmediate)              # runs the Becker and Woessmann example
```
