# ivmediate for Stata

Causal mediation analysis for linear IV models: `ivmediate` splits the effect of a treatment into a
direct effect and an indirect effect through a mediator, using instruments for the treatment and none
for the mediator.

[Main page](../README.md) · [Python](../python/README.md) · [R](../r/README.md) · [Stata](../stata/README.md) · [MATLAB](../matlab/README.md)

## Installation

```stata
net install ivmediate, from("https://raw.githubusercontent.com/aferrara87/ivmediate/main/stata") replace
```

- **Requirements:** Stata 16 or later. No other packages are needed; the computations run in Mata.
- **Replaces `ivmediate` 1.0.4** (SSC, 2020). Version 2.0 keeps the command name and the basic
  syntax, but the estimator, the defaults (robust standard errors and Anderson–Rubin sets) and the
  output change, and the `full` option is gone. Re-run results obtained with 1.0.4. The section
  "What's new in 2.0" of `help ivmediate` lists the changes. `replace` overwrites an installed 1.0.4.
- **Check the version:** `which ivmediate` prints the path of the file Stata uses and its header,
  which starts with `*! ivmediate 2.0.0  17 Oct 2026`. After an estimation, `e(version)` is `2.0.0`.

## Quick start

Becker and Woessmann (2009): 426 Prussian counties, with income tax revenue per capita as the outcome
(`inctaxpc`), the share of Protestants as the treatment (`f_prot`), the literacy rate as the mediator
(`f_rw`), the distance to Wittenberg as the instrument (`kmwittenberg`) and 11 county controls. This is
Table 6, column 1, of the paper. The data are in
[`examples/becker_woessmann/`](../examples/becker_woessmann/); run this from that folder:

```stata
import delimited using "becker_woessmann_2009.csv", clear
* The CSV codes missing values as NA, which import delimited reads as text.
destring inctaxpc f_prot f_rw kmwittenberg f_young f_jew f_fem f_ortsgeb f_pruss ///
    f_blind f_deaf f_dumb hhsize lnpop gpop, replace force

ivmediate inctaxpc f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop, ///
    treatment(f_prot) mediator(f_rw) instrument(kmwittenberg)
```

```
------------------------------------------------------------------------------
Linear IV Mediation Analysis                             ivmediate 2.0.0
------------------------------------------------------------------------------
Outcome:   inctaxpc                            Obs           =       426
Treatment: f_prot                              Instruments   =         1
Mediator:  f_rw                                Covariates    =        11
SE:        robust (HC1)

Mediation estimates
------------------------------------------------------------------------------
Parameter                 Estimate  Std. err.   [95% conf. interval]
------------------------------------------------------------------------------
tau_total (Wald)            0.5865     0.2425      0.1112     1.0617 **
tau (direct T->Y)          -0.2992     0.2141     -0.7187     0.1204
theta (M->Y | T)            4.7091     2.1638      0.4681     8.9502 **
gamma (T->M)                0.1881     0.0269      0.1354     0.2407 ***
theta*gamma (indirect)      0.8857     0.4076      0.0869     1.6845 **
------------------------------------------------------------------------------
*** p<0.01, ** p<0.05, * p<0.1 (two-sided Wald, normal critical values)
  Mediator explains 151.01% of the total effect.

First-stage diagnostics  (SE: robust)
------------------------------------------------------------------------------
  pi1 (T on Z | X)         = -0.0984 (0.0113)
  F(T, Z|X)                =      75.69   OP critical value = 23.11 (K_eff = 1.00, tau = 0.1, alpha = 0.05)
  F(M, Z|T,X)              =      16.36   OP critical value = 23.11 (K_eff = 1.00, tau = 0.1, alpha = 0.05)

Anderson-Rubin 95% confidence sets  (weak-IV-robust, SE: robust)
------------------------------------------------------------------------------
  theta (M->Y | T)  [0.4311, 10.2039]
                      AR p-value at 0 = 0.0333   (exact Fieller inversion)
  gamma (T->M)      [0.1400, 0.2488]
                      AR p-value at 0 = 0.0000   (exact Fieller inversion)
  tau_total         [0.1258, 1.1030]
                      AR p-value at 0 = 0.0126   (exact Fieller inversion)
  AR sets (chi2(1) critical value 3.841) remain valid whatever the first-stage strength;
  the Wald CIs above require the effective F to clear the OP critical value.

kappa-based sensitivity (MCA diagnostics)
------------------------------------------------------------------------------
  c_TM = Cov(eps_T, eps_M)   =     -89.2747
  kappa* (breakeven)         =    -420.4053
  kappa*_CI (inf. breakeven) =     -41.7858
  rho*_TY (as correlation)   =      -0.2567
------------------------------------------------------------------------------
Excluded instrument(s): kmwittenberg
```

These are the numbers of the paper: τ_total 0.586 (0.242), τ −0.299 (0.214), θ 4.709 (2.164) with
the AR set [0.43, 10.20], γ 0.188 (0.027), θγ 0.886 (0.408). The main page explains
[how to read them](../README.md#example-becker-and-woessmann-2009) and
[the estimator](../README.md#the-estimator). In a clone of the repository,
[`bw.do`](../examples/becker_woessmann/bw.do) runs the same example with the repository's `stata/`
folder on the ado-path, without installing.

## Syntax

```
ivmediate depvar [indepvars] [if] [in] [weight], treatment(varname) mediator(varname)
          instrument(varlist) [options]
```

*depvar* is the outcome *Y* and *indepvars* are the controls *X*. `treatment()`, `mediator()` and
`instrument()` are required. `aweight`s and `pweight`s are allowed.

## Options

Each table gives the shortest abbreviation, the default and, where there is one, a reference to the
paper (Ferrara, Gold, Heblich and Pinto 2026; section, remark and proposition numbers of the current
draft).

### Model and controls

| Option | Short | Default | What it does | Paper |
|---|---|---|---|---|
| *depvar* | | | The outcome *Y*. | Section 2 |
| `treatment(varname)` | `treat()` | required | The treatment *T*. | Section 2 |
| `mediator(varname)` | `med()` | required | The mediator *M*. The outcome, treatment and mediator must be three different variables. | Section 2 |
| `instrument(varlist)` | `inst()` | required | One or more excluded instruments *Z* for the treatment; the mediator needs no instrument of its own. With several instruments each equation is estimated by overidentified 2SLS. | Section 3.3, Proposition 2; Section 4.2 |
| *indepvars* | | none | The controls *X*, partialed out and not reported. Factor variables (`i.region`) enter as dummies without the base level; time-series operators are allowed. Wildcards and ranges (`x*`, `x1-x8`) never pick up the outcome, treatment, mediator, instruments or the absorb, cluster, weight and coordinate variables; such matches are dropped with a note. Exactly collinear controls are dropped with a note. | |
| `if`, `in` | | | Restrict the sample. The sample is the observations with no missing value in any variable the model uses and, with weights, a positive weight. | |

### Fixed effects and weights

| Option | Short | Default | What it does | Paper |
|---|---|---|---|---|
| `absorb(varlist)` | `a()` | none | One or two categorical variables (numeric or string) absorbed as fixed effects. Singleton groups are dropped, iteratively across the dimensions. Two-way fixed effects are removed exactly when the smaller dimension has at most 3,000 levels, and by accelerated alternating projections otherwise. | |
| `[aweight=`*w*`]`, `[pweight=`*w*`]` | `[aw=`*w*`]`, `[pw=`*w*`]` | none | Analytic or sampling weights; the results do not depend on their scale. Observations with zero weight are dropped and negative weights are an error. With `pweight`s, `vce(unadjusted)` becomes `vce(robust)`. The κ block uses weighted (Hajek) moments. | Appendix F.4, step 3 |

### Standard errors

| Option | Short | Default | What it does | Paper |
|---|---|---|---|---|
| `vce(robust)` | `vce(r)` | default | Heteroskedasticity-robust (HC1). | |
| `vce(unadjusted)` | `vce(un)`, `vce(ols)`, `vce(classical)` | | Homoskedastic. | |
| `vce(cluster` *clustvar*`)` | `vce(cl` *clustvar*`)` | | One-way cluster-robust (Liang–Zeger CR1). The cluster variable may be numeric or string; at least 2 clusters. Fixed effects nested in the clusters do not count in the small-sample factor (as in `reghdfe`). | |
| `vce(conley)` | `vce(spatial)` | | Conley spatial standard errors. Needs `latitude()`, `longitude()` and `cutoff()`. | |
| `latitude(varname)` | `lat()` | | Latitude in degrees, or the y-coordinate under `metric(euclidean)`. | |
| `longitude(varname)` | `lon()` | | Longitude in degrees, or the x-coordinate under `metric(euclidean)`. | |
| `cutoff(#)` | `cut()` | required with `vce(conley)` | Distance beyond which observations are independent: kilometers under `haversine`, coordinate units under `euclidean`. Report results for several cutoffs. | |
| `kernel(string)` | `kern()` | `uniform` | `uniform` or `bartlett`. | |
| `metric(string)` | `metr()` | `haversine` | `haversine` (great-circle distance) or `euclidean`. | |

The spatial options are an error without `vce(conley)`. Neither kernel guarantees a positive
variance in two dimensions; a negative variance is reported as "SE unavailable" (`help ivmediate`,
Remarks).

### Instruments and Anderson–Rubin sets

| Option | Short | Default | What it does | Paper |
|---|---|---|---|---|
| `noar` | | sets computed | Skips the Anderson–Rubin sets for θ, γ and the total effect. `ar` is accepted and changes nothing. | Remark 3 (Section 4.2); Section 8 |
| `arpoints(#)` | `arp()` | 4001 | Grid points of the AR search with more than one instrument (at least 3). With one instrument the set is exact (Fieller inversion) and needs no grid. | |
| `arwidth(#)` | `arw()` | 200 | Half-width of the AR search grid, in standard errors around the estimate; also the range on which a Conley variance is checked for negativity. | |

The AR sets use the confidence level of `level()` and the same type of standard errors as the Wald
intervals.

### Weak-instrument diagnostics

| Option | Short | Default | What it does | Paper |
|---|---|---|---|---|
| `optau(#)` | `opt()` | 0.10 | Bias tolerance τ of the Olea–Pflueger test, in (0, 1). | Remark 3 (Section 4.2); Appendix F.4, step 2; Section 7.3 |
| `opalpha(#)` | `opa()` | 0.05 | Size α of the Olea–Pflueger test, in (0, 1). | |

Both first stages are tested: *T* on *Z*, and the conditional first stage of *M* on *Z* given *T*. The
critical value is that of the simplified test at the estimated effective degrees of freedom *K_eff*
(23.11 with one instrument at the defaults).

### Sensitivity and display

| Option | Short | Default | What it does | Paper |
|---|---|---|---|---|
| `kappatable` | `kappa` | off | Prints the κ-adjusted table: θ adjusted for κ0 = 0, 10, 25, 50, 75, 90 and 100% of κ*, with Wald intervals, and a note when the conditional first stage is weak. The table is always stored in `e(kappa_adj)`. | Section 5.1.1, Proposition 5; Appendix F.4, step 3 |
| `r2sensitivity` | `r2` | off | Prints the R²-based section: the reduced-form F, RV(θ), XRV(θ) and their components. The values are always stored. | Appendix F.2 (Definitions 2–3) and F.3 |
| `level(#)` | `l()` | 95, or as set by `set level` | Confidence level of the Wald, AR and bootstrap intervals and of the RV threshold. | |
| `nostars` | `nostar` | stars shown | Suppresses the significance stars. | |
| `starcutoffs(# # #)` | `starcut()` | 0.01 0.05 0.10 | The three p-value thresholds of the stars, in any order. | |
| `residuals(`*stub*`)` | `resid()` | none | Saves the residuals behind the κ block as new variables *stub*`_eT`, *stub*`_eM` and *stub*`_eY`. | Appendix F.1 |

The κ summary (c_TM, κ*, κ*_CI and ρ*_TY) is always printed; ρ*_TY is Remark 8 (Section 5.1.1).

### LASSO

| Option | Short | Default | What it does | Paper |
|---|---|---|---|---|
| `lasso` | | off | Post-double-selection LASSO: selects controls for *Y*, *T*, *M* and each instrument from a basis expansion of the controls, then estimates on the union of the selected terms. *indepvars* may contain continuous variables and `i.` factors but no interactions (the expansion builds them). Not available with `bootstrap`. | Remark 6 (Section 4.2); Appendix H (Algorithm 1 in H.2) |
| `lpenalty(plugin)` | | `plugin` | The penalty rule; the plug-in penalty is the only one offered. | |
| `expandx(string)` | | `square cube` | Expansion of the continuous controls: any of `square`, `cube` and `interact_with_X` (pairwise products), or `none`. | |
| `expanddum(string)` | | `interact_with_X` | Expansion with the factor dummies: any of `interact_with_X` (dummy × control) and `interact_with_dummies` (dummy × dummy of different factors), or `none`. | |
| `lkeep(patterns)` | | none | Expansion terms that are always kept, matched by name with `strmatch()` patterns, for example `turn^2`, `2.rep78#turn` or `c1*`. A pattern that matches no term is an error. | |
| `showsel` | | off | Prints the terms each LASSO selected, with its penalty and convergence. | |
| `lcslack(#)` | | 1.1 | Slack *c* of the plug-in penalty λ = 2*c*√*n* Φ⁻¹(1 − α_p/(2*p*)), with *p* candidate terms. | |
| `lalphap(#)` | | 0.1/ln(*n*) | Tail probability α_p of the penalty, in (0, 1). | |
| `lniter(#)` | | 15 | Maximum number of iterations of the penalty loadings. | |
| `ltol(#)` | | 1e-6 | Convergence tolerance of the loadings. | |
| `llambda(#)` | | plug-in | Fixes λ instead of iterating the plug-in penalty. | |

### Bootstrap

| Option | Short | Default | What it does | Paper |
|---|---|---|---|---|
| `bootstrap` | `b` | off | Adds bootstrap standard errors and intervals for the five effects, κ* and ρ*_TY, printed in a separate block. The estimates table, `e(b)` and `e(V)` stay analytic. Not available with `lasso`. | |
| `breps(#)` | `bre()` | 999 | Number of replications (at least 50). | |
| `bseed(#)` | `bse()` | 20260317 | Random-number seed. The command saves `c(rngstate)` and restores it afterward, so your random-number stream is unchanged. | |
| `bmethod(string)` | `bm()` | `block` with `vce(cluster)`, else `pairs` | `pairs` resamples observations; `block` resamples clusters and needs `vce(cluster)`. | |
| `bcitype(string)` | `bci()` | `basic` | `basic`, `percentile` or `studentized` (bootstrap-t). | |
| `bootidx(matname)` | | none | Supplies the resamples instead of drawing them: one row per replicate, holding 1-based positions in the estimation sample (`pairs`) or cluster numbers in order of first appearance (`block`). | |
| `bdiagnostics` | `bdiag` | off | Also bootstraps RV(θ) and XRV(θ). | |

## Output

### Printed blocks

The printout has the same blocks, in the same order, as in the other languages, with ASCII labels in
place of the Greek letters. The main page walks through them with the Python output; this is how the
labels map:

| Block | Stata | Main page (Python) |
|---|---|---|
| Mediation estimates | `tau_total (Wald)`, `tau (direct T->Y)`, `theta (M->Y \| T)`, `gamma (T->M)`, `theta*gamma (indirect)` | τ_total, τ, θ, γ, θγ |
| First-stage diagnostics | `pi1 (T on Z \| X)`, `F(T, Z\|X)`, `F(M, Z\|T,X)` | π_1, F_{T,Z\|X}, F_{M,Z\|T,X} |
| Anderson–Rubin sets | `theta`, `gamma`, `tau_total` | θ, γ, τ_total |
| κ-based sensitivity | `c_TM`, `kappa*`, `kappa*_CI`, `rho*_TY` | c_TM, κ*, κ*_CI, ρ*_TY |

- **Header:** the outcome, treatment and mediator, the number of observations, instruments and
  controls, and the type of standard errors.
- **Mediation estimates:** the five effects with Wald standard errors, normal-based intervals and
  stars. The standard error of θγ is the delta-method one, with the covariance of γ and θ across the
  two equations (Remark 4; Appendix D.1). Below the table, Stata adds the share of the total effect
  that runs through the mediator (`e(mepct)`).
- **First-stage diagnostics:** π_1, the coefficient of *Z* in the first stage of *T*; the F statistic
  of that first stage, and of the conditional first stage of *M* on *Z* given *T*, each with its
  Olea–Pflueger critical value. With one instrument F is the effective F (the squared t statistic).
  With several, Stata prints the joint Wald F and the OP effective F on separate lines, and the
  squared t statistic of each instrument in the conditional first stage.
- **Anderson–Rubin sets** for θ, γ and τ_total, each with the AR p-value at 0. Unbounded or
  disconnected sets are printed as such, with a note.
- **κ-based sensitivity:** c_TM = Cov(ε_T, ε_M), the breakeven κ*, the inference breakeven κ*_CI
  (marked "uninformative" when the interval for θ already contains 0) and ρ*_TY, the breakeven as a
  correlation.
- **Footer:** the excluded instruments.

With `kappatable r2sensitivity`, the κ block gets the table, and the R² section follows (Becker and
Woessmann; the other blocks are as above):

```
kappa-based sensitivity (MCA diagnostics)
------------------------------------------------------------------------------
  c_TM = Cov(eps_T, eps_M)   =     -89.2747
  kappa* (breakeven)         =    -420.4053
  kappa*_CI (inf. breakeven) =     -41.7858
  rho*_TY (as correlation)   =      -0.2567

  kappa-adjusted theta (kappa0 swept as a fraction of kappa*):
  % kappa*       kappa0    theta_adj        CI_lo        CI_hi
        0%       0.0000       4.7091       0.4681       8.9502
       10%     -42.0405       4.2382      -0.0029       8.4793
       25%    -105.1013       3.5318      -0.7092       7.7729
       50%    -210.2027       2.3546      -1.8865       6.5956
       75%    -315.3040       1.1773      -3.0638       5.4183
       90%    -378.3648       0.4709      -3.7701       4.7120
      100%    -420.4053       0.0000      -4.2411       4.2411
  Note: F(M, Z|T,X) = 16.36 is below the OP critical value (23.11, tau = 0.1), so c_TM is estimated imprecisely; the
  AR set for theta above is robust to weak instruments.

R^2-based sensitivity (reduced-form robustness)
------------------------------------------------------------------------------
  (t-statistics under robust SEs)
  F (reduced form, = t^2) =       4.53
  RV(theta)               =     0.0078
  XRV(theta)              =     0.0016
    (RV_rf = 0.0078, RV_fs = 0.0972)
------------------------------------------------------------------------------
```

With `bootstrap`, a bootstrap block comes last (Becker and Woessmann, default settings):

```
Bootstrap  (method = pairs, B = 999, seed = 20260317)
------------------------------------------------------------------------------
  CI type: basic (95%);  999 successful replicates
  Parameter                 Estimate    Boot SE          [Boot 95% CI]
  tau_total (Wald)            0.5865     0.2481      0.0108     1.0134 **
  tau (direct T->Y)          -0.2992     0.2546     -0.7512     0.2185
  theta (M->Y | T)            4.7091     2.4924     -0.4906     8.9944 *
  gamma (T->M)                0.1881     0.0274      0.1263     0.2346 ***
  theta*gamma (indirect)      0.8857     0.4448     -0.0669     1.6868 **
  kappa*                   -420.4053   208.0278   -806.1157     9.8988
  rho*_TY                    -0.2567     0.1203     -0.4908    -0.0327
  kappa*_CI with the bootstrap SE of theta = 15.7093  (uninformative: the CI already contains 0)
  The mediation estimates table, e(b) and e(V) keep the analytic SEs.
------------------------------------------------------------------------------
```

With `lasso`, the header says how many expansion terms were selected, and a LASSO block after the
sensitivity sections lists them; `showsel` adds the selection of each LASSO, with its penalty and
convergence.

### Stored results

`ivmediate` is an e-class command. These are the results most users need; `help ivmediate` lists all
of them.

**Estimates and diagnostics (scalars)**

| Result | Content |
|---|---|
| `e(N)` | Observations used |
| `e(total)`, `e(se_total)` | Total effect τ_total and its standard error |
| `e(tau)`, `e(se_tau)` | Direct effect τ |
| `e(theta)`, `e(se_theta)` | Mediator effect θ |
| `e(gamma)`, `e(se_gamma)` | Effect of the treatment on the mediator, γ |
| `e(indirect)`, `e(se_indirect)` | Indirect effect θγ (delta-method standard error) |
| `e(cov_gamma_theta)` | Covariance of γ and θ |
| `e(mepct)` | Indirect effect as a percentage of the total effect |
| `e(F_TZ)`, `e(F_TZ_op)`, `e(K_eff_TZ)`, `e(F_TZ_crit)` | First stage of *T*: F (joint Wald), OP effective F, *K_eff* and the OP critical value |
| `e(F_MZ_T)`, `e(F_MZ_T_op)`, `e(K_eff_MZ_T)`, `e(F_MZ_T_crit)` | The same for the conditional first stage of *M* given *T* |
| `e(F_MZ_T_above_crit)` | 1 if the effective F of the conditional first stage reaches its OP critical value, 0 if not, missing if either is unavailable |
| `e(c_TM)` | Cov(ε_T, ε_M) |
| `e(kappa_star)`, `e(kappa_star_CI)` | Breakeven κ* and inference breakeven κ*_CI |
| `e(kappa_star_CI_informative)` | 1 if the interval for θ excludes 0, 0 if it contains 0 (then κ*_CI is uninformative) |
| `e(rho_TY_star)` | Breakeven as a correlation, ρ*_TY |
| `e(RV_theta)`, `e(XRV_theta)` | Robustness values RV(θ) and XRV(θ) |
| `e(V_unavailable)` | 1 if a variance is unavailable; `e(V)` is then not posted |
| `e(k_inst)`, `e(k_cov)` | Number of instruments and of controls used |
| `e(N_clust)`, `e(df_fe)`, `e(N_singleton)` | Clusters, absorbed degrees of freedom, singleton observations dropped |

**Anderson–Rubin sets** (*x* is `theta`, `gamma` or `total`)

| Result | Content |
|---|---|
| `e(ar_`*x*`_shape)` (macro) | `interval`, `rays`, `rays_left_only`, `rays_right_only`, `line`, `grid` (more than one instrument), `empty` or `unavailable` |
| `e(ar_`*x*`_lower)`, `e(ar_`*x*`_upper)` | Endpoints of an interval |
| `e(ar_`*x*`_left_upper)`, `e(ar_`*x*`_right_lower)` | Inner endpoints of a pair of rays |
| `e(ar_`*x*`_seg)` (matrix) | With more than one instrument: one row per segment, with columns `lower`, `upper`, `open_left` and `open_right` (1 if the segment reaches the edge of the search grid) |
| `e(ar_`*x*`_p)`, `e(ar_`*x*`_stat0)`, `e(ar_`*x*`_contains0)` | AR p-value and statistic at 0; 1 if the set contains 0 |
| `e(ar_cv)`, `e(ar_df)`, `e(ar_level)` | χ² critical value, its degrees of freedom (the number of instruments) and the level |

**Bootstrap** (with `bootstrap`)

| Result | Content |
|---|---|
| `e(boot_point)`, `e(boot_se)`, `e(boot_ci)` (matrices) | Estimate, bootstrap standard error and interval (`lower`, `upper`) of `total`, `tau`, `theta`, `gamma`, `indirect`, `kappa_star` and `rho_TY_star`, plus RV and XRV with `bdiagnostics` |
| `e(boot_analytic_se)` (matrix) | Analytic standard errors of the same quantities, missing where none exists |
| `e(boot_dist)` (matrix) | The replicates, one row per replication |
| `e(boot_se_total)` … `e(boot_se_indirect)` | Bootstrap standard errors of the five effects |
| `e(boot_B)`, `e(boot_n_failed)`, `e(boot_seed)` | Replications, failed replications, seed (missing with `bootidx()`) |
| `e(boot_method)`, `e(boot_citype)` (macros) | Resampling scheme and interval type |

**Macros and matrices**

| Result | Content |
|---|---|
| `e(cmd)`, `e(version)`, `e(cmdline)` | `ivmediate`, `2.0.0`, and the command as typed |
| `e(depvar)`, `e(treat)`, `e(med)`, `e(inst)`, `e(controls)` | Variable names; `e(controls)` lists the controls used |
| `e(vce)`, `e(se_type)` | Type of standard errors: Stata names (`unadjusted`, `robust`, `cluster`, `conley`) and the names of the other languages (`classical`, `robust`, `cluster`, `spatial`) |
| `e(clustvar)`, `e(absorb)`, `e(wtype)`, `e(wexp)` | Cluster variable, absorbed variables, weights |
| `e(b)`, `e(V)` | Total, direct and indirect effect, and a diagonal matrix of their variances (not a joint covariance matrix: the three come from different equations) |
| `e(pi1)`, `e(se_pi1)` | First stage of *T* on *Z* (1 × *K*) |
| `e(F_MZ_T_individual)` | Squared t statistic of each instrument in the conditional first stage (*K* × 1) |
| `e(kappa_adj)` | The κ-adjusted table (columns `kappa0`, `fraction`, `theta_adj`, `ci_lo`, `ci_hi`) |
| `e(sample)` | Marks the estimation sample |

## Differences from the other languages

Apart from the bootstrap draws (below), the numbers are the same as in Python, R and MATLAB: the
[validation harness](../validation/) checks this to a relative tolerance of 1e-8. A user notices
these differences:

- **Layout and names.** Stata prints its own layout, with ASCII names such as `theta*gamma` and
  `kappa*`, and adds the share of the total effect through the mediator.
- **One command.** Estimation and printing are one command, and the bootstrap is the `bootstrap`
  option. The other languages have separate functions for printing and for the bootstrap.
- **Option names.** Python, R and MATLAB use the same argument names, which map to Stata as follows:

  | Stata | Python, R and MATLAB |
  |---|---|
  | *depvar*, *indepvars*; `treatment()`, `mediator()`, `instrument()` | `Y`, `X`; `T`, `M`, `Z` |
  | `i.` factor variables among the controls | `dummies` |
  | `absorb()` | `fe` |
  | `vce(unadjusted)`, `vce(robust)`, `vce(cluster` *clustvar*`)`, `vce(conley)` | `se = "classical"`, `"robust"`, `"cluster"` with `clusterid`, `"spatial"` |
  | `[aweight=`*w*`]`, `[pweight=`*w*`]` | `weights` with `wtype = "aweight"` or `"pweight"` |
  | `level(95)` | `alpha = 0.05` |
  | `latitude()`, `longitude()`, `cutoff()`, `kernel()`, `metric()` | `lat`, `lon`, `spatial_cutoff`, `spatial_kernel`, `spatial_metric` |
  | `optau()`, `opalpha()` | `op_tau`, `op_alpha` |
  | `noar`, `arpoints()`, `arwidth()` | `ar_ci` set to false, `ar_grid_points`, `ar_grid_halfwidth` |
  | `kappatable`, `r2sensitivity`, `nostars`, `starcutoffs()` | the print function's `kappa_table`, `r2_sensitivity`, `stars`, `star_cutoffs` |
  | `lasso`, `lpenalty()`, `lkeep()`, `expandx()`, `expanddum()`, `showsel` | `lasso`, `lasso_penalty`, `lasso_keep`, `expand_x`, `expand_dummies`, `show_selection` |
  | `lcslack()`, `lalphap()`, `lniter()`, `ltol()`, `llambda()` | `lasso_c_slack`, `lasso_alpha_p`, `lasso_n_iter_max`, `lasso_tol`, `lasso_lambda_override` |
  | `breps()`, `bseed()`, `bmethod()`, `bcitype()`, `bdiagnostics` | the bootstrap function's `B`, `seed`, `method`, `ci_type`, `bootstrap_diagnostics` |
  | `residuals()` | the returned residuals `eps_T_hat`, `eps_M_hat`, `eps_Y_hat` |
  | `bootidx()` | `boot_indices` in R and MATLAB (a testing hook with 0-based draws); none in Python |

- **No per-variable expansion rules.** The other languages' `expand_override` has no Stata
  counterpart: `expandx()` and `expanddum()` apply to all controls.
- **No infinite values.** Stata variables cannot hold ±∞, so the other languages' error for infinite
  values does not arise. Unbounded AR sets have no infinite endpoints either: an endpoint that would
  be infinite is missing (`.`), and with more than one instrument `open_left` and `open_right` in
  `e(ar_`*x*`_seg)` flag the segments that reach the edge of the search grid. `e(ar_`*x*`_shape)` says
  what kind of set it is.
- **Variable patterns.** Controls take Stata's wildcards and ranges (`x*`, `x?`, `x1-x8`) rather than
  shell-style globs, so a character class such as `x[12]` is not available.
- **Bootstrap draws.** The replicates come from Stata's own random-number generator, seeded with
  `bseed()`, and your `c(rngstate)` is restored afterward. With the same seed, the bootstrap results
  therefore differ from the other languages' by simulation noise, which can be large for a parameter
  with heavy tails: for θ in the Becker–Woessmann example, whose conditional first stage is only
  moderately strong, the bootstrap SE ranges from 2.4 to 4.5 across the four languages (B = 999), while
  the bootstrap CIs differ much less. `bootidx()` replays a given set of resamples exactly.
- **Speed.** Stata/SE runs on one core; Stata/MP may be faster. Timings are on the
  [main page](../README.md#speed).

## Help inside Stata

```stata
help ivmediate
```

The help file documents every option and stored result, the conventions (degrees of freedom and
critical values) and "What's new in 2.0".
