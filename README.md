# ivmediate

**Causal mediation analysis with a single instrument, in Python, R, Stata and MATLAB.**

`ivmediate` implements the estimator of Ferrara, Gold, Heblich and Pinto (2026), *Identifying Causal
Mediation Effects with a Single Instrument* (NBER Working Paper No. [number]). With one instrument for the
treatment, it splits the treatment's effect on an outcome into a direct effect and an indirect effect
through a mediator. It also reports the weak-instrument diagnostics and the sensitivity analyses that the
paper develops. The four implementations give the same results; they are tested against each other in
[`validation/`](validation/).

**Documentation:** [Python](python/README.md) · [R](r/README.md) · [Stata](stata/README.md) ·
[MATLAB](matlab/README.md)

**On this page:** [The estimator](#the-estimator) ·
[Installation and documentation](#installation-and-documentation) ·
[Example: Becker and Woessmann (2009)](#example-becker-and-woessmann-2009) · [Speed](#speed) ·
[Citation](#citation)

## The estimator

In the linear model of the paper (Section 2),

$$T = \pi_1 Z + X'\beta_T + \varepsilon_T,\qquad M = \gamma T + X'\beta_M + \varepsilon_M,\qquad
Y = \tau T + \theta M + X'\beta_Y + \varepsilon_Y,$$

- the instrument $Z$ shifts the treatment $T$;
- the treatment affects the mediator $M$ ($\gamma$) and the outcome $Y$ directly ($\tau$);
- the mediator affects the outcome ($\theta$).

The total effect of $T$ on $Y$ is $\tau_{\text{total}} = \tau + \theta\gamma$, and $\theta\gamma$ is the
indirect effect through the mediator. A standard IV regression identifies $\tau_{\text{total}}$ and
$\gamma$, but not $\theta$ and $\tau$ separately (Section 2.4).

**Identification.** The paper's identifying assumption is *mediated confounding* (MCA, Section 3). It lets
the mediator be endogenous, but restricts how the unobservables of the treatment and the outcome are
linked. In the linear model, $\varepsilon_T$ and $\varepsilon_Y$ are uncorrelated, so any dependence
between them runs through the mediator. When the treatment's and the mediator's unobservables are
correlated (*mediated relevance*, MRA), the instrument predicts the mediator given the treatment. It then
serves as an instrument for $M$ conditional on $T$ (Section 3.3, Proposition 2).

**Estimation.**
- $\theta$ and $\tau$: two-stage least squares of $Y$ on $M$ and $T$, instrumenting $M$ with $Z$
  (conditional 2SLS, Section 4.2).
- $\gamma$: the IV regression of $M$ on $T$; $\tau_{\text{total}}$: the IV regression of $Y$ on $T$.
- The indirect effect $\theta\gamma$ comes with a delta-method standard error (Remark 4 and
  Appendix D.1).

**Diagnostics.** Two first stages matter: the usual one, of $T$ on $Z$, and the *conditional* first
stage, of $M$ on $Z$ given $T$. The conditional first stage identifies $\theta$ and can be weak even when
the usual one is strong (Remark 3). `ivmediate` reports both as Olea–Pflueger effective F-statistics,
with their critical values. It also reports Anderson–Rubin confidence sets for $\theta$, $\gamma$ and
$\tau_{\text{total}}$, which remain valid when the instruments are weak.

**Sensitivity.** Mediated confounding cannot be tested with instruments for the treatment alone
(Proposition 4). The $\kappa$-based sensitivity analysis instead asks how large a violation,
$\kappa = \mathrm{Cov}(\varepsilon_T, \varepsilon_Y)$, would have to be to overturn the estimate of
$\theta$ (Section 5.1.1). It reports:
- the breakeven $\kappa^*$;
- its confidence-interval version $\kappa^*_{CI}$;
- the same breakeven as a correlation, $\rho^*_{TY}$ (Remark 8).

On request, `ivmediate` also prints a table of $\kappa$-adjusted estimates, and robustness values for
confounders of the reduced form (Appendix F.2–F.4).

**Options.**
- **Controls:** categorical dummies, and post-double-selection LASSO when there are many candidate
  controls (Appendix H).
- **Fixed effects and weights:** one- and two-way fixed effects, and weights.
- **Standard errors:** classical, robust, cluster and spatial (Conley).
- **Instruments:** one or several.
- **Bootstrap.**

The [language pages](#installation-and-documentation) document each option.

## Installation and documentation

Each language has its own documentation page, linked below. It covers the installation, the syntax,
every option with references to the paper, the output, and what differs from the other languages.

**[Python documentation](python/README.md)** (Python 3.9 or later)

```bash
pip install "git+https://github.com/aferrara87/ivmediate#subdirectory=python"
```

**[R documentation](r/README.md)** (R 4.2 or later)

```r
# install.packages("remotes")
remotes::install_github("aferrara87/ivmediate", subdir = "r")
```

**[Stata documentation](stata/README.md)** (Stata 16 or later). The command replaces `ivmediate` 1.0.4
on SSC; its help file lists what is new in version 2.0.

```stata
net install ivmediate, from("https://raw.githubusercontent.com/aferrara87/ivmediate/main/stata") replace
```

**[MATLAB documentation](matlab/README.md)** (R2016b or later, no toolboxes). Download this repository
(Code → Download ZIP) or clone it, then add its `matlab` folder to the MATLAB path:

```matlab
addpath('<repo>/matlab')   % <repo>: the folder you unzipped or cloned
```

## Example: Becker and Woessmann (2009)

Becker and Woessmann (2009) ask whether Protestantism made Prussian counties more prosperous through
literacy. The variables:
- $Y$: income tax revenue per capita in 1877, in Pfennig (hundredths of a Mark; Becker and Woessmann
  report it in Marks);
- $T$: the share of Protestants (in percent);
- $M$: the literacy rate (in percent);
- $Z$: the distance to Wittenberg (in km);
- 11 county-level controls, for 426 counties.

The example reproduces column 1 of Table 6 in the paper (Section 8). The data come from Becker, Sascha O.,
and Ludger Woessmann (2009), "Was Weber Wrong? A Human Capital Theory of Protestant Economic History,"
*Quarterly Journal of Economics* 124 (2), 531–596. The folder
[`examples/becker_woessmann/`](examples/becker_woessmann/) holds the data and the same example in R, Stata
and MATLAB.

```python
import pandas as pd
import ivmediate

bw = pd.read_csv("examples/becker_woessmann/becker_woessmann_2009.csv")
res = ivmediate.ivmediate(
    bw, Y="inctaxpc", T="f_prot", M="f_rw", Z="kmwittenberg",
    X="f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop",
)
ivmediate.print_results(res, title="Becker and Woessmann (2009)")
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

**Reading the output, block by block.** One percentage point of the Protestant share is a small change, so
the effects are also given for the move from an all-Catholic to an all-Protestant county (100 points), the
comparison Becker and Woessmann use.

**Mediation estimates.**
- **Total effect.** A one-percentage-point higher Protestant share raises income tax revenue per capita
  by 0.59 Pfennig ($\tau_{\text{total}}$, the usual IV estimate). An all-Protestant county thus collects
  0.59 Marks more per head than an all-Catholic county, about 30 percent of the average of 1.98 Marks;
  Becker and Woessmann (2009, p. 569) report the same figure.
- **Mediator effect.** Holding the Protestant share fixed, one percentage point of literacy raises income
  tax revenue per capita by 4.71 Pfennig ($\theta$), so ten points raise it by 0.47 Marks.
- **Effect on the mediator.** A one-percentage-point higher Protestant share raises the literacy rate by
  0.19 percentage points ($\gamma$): an all-Protestant county is 19 points more literate than an
  all-Catholic one (18.9 points in Becker and Woessmann 2009, p. 558).
- **Indirect and direct effects.** The indirect effect through literacy, $\theta\gamma = 0.89$ Pfennig per
  percentage point (0.89 Marks from an all-Catholic to an all-Protestant county), is larger than the total
  effect. The direct effect, holding literacy fixed, is negative and not significant ($\tau = {-0.30}$
  Pfennig per point).
- **Decomposition.** The estimates add up: $\tau + \theta\gamma = \tau_{\text{total}}$
  (Corollary 1). Stars come from the Wald p-values.

**First-stage diagnostics.**
- **Instrument.** Each kilometer from Wittenberg lowers the Protestant share by 0.098 percentage points
  ($\pi_1$), about 10 points per 100 km.
- **Usual first stage.** It is strong: $F_{T,Z|X} = 75.69$, well above the Olea–Pflueger critical value
  of 23.11 (10% worst-case bias, 5% level).
- **Conditional first stage.** It identifies $\theta$ and is weaker: $F_{M,Z|T,X} = 16.36$. That is
  below 23.11 but above the 15.06 for a 20% bias tolerance (Remark 3). The Wald interval for $\theta$
  should therefore be read together with the Anderson–Rubin set.

**Anderson–Rubin sets.** These sets are valid however weak the instruments are.
- The set for $\theta$, $[0.43, 10.20]$, excludes zero ($p = 0.033$ at zero). It is wider and more
  skewed than the Wald interval, which reflects the moderate conditional first stage.
- The sets for $\gamma$ and $\tau_{\text{total}}$ are close to their Wald intervals.

**κ-based sensitivity.**
- **$c_{TM}$.** This is the estimated covariance of the treatment and mediator unobservables. It is
  the source of identification (mediated relevance).
- **$\kappa^{\ast} = {-420.4}$.** The covariance between the treatment and outcome unobservables that would
  move $\theta$ to zero; the negative sign gives the direction of that confounding.
- **$\kappa^{\ast}_{CI} = {-41.8}$.** The covariance at which the 95% interval for $\theta$ starts to include
  zero.
- **$\rho^{\ast}_{TY} = {-0.257}$.** The breakeven expressed as a correlation. Overturning the mediator
  effect would take a correlation of about −0.26 between the unobserved drivers of Protestantism and of
  prosperity (Remark 8): a sizable violation of mediated confounding.

**Optional output.** Two more sections are printed on request. Their values are always stored in the
results.

```python
ivmediate.print_results(res, title="Becker and Woessmann (2009)", kappa_table=True, r2_sensitivity=True)
```

```
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

**The κ table.** It shows how $\theta$ and its interval move if mediated confounding is violated by
a fraction of the breakeven $\kappa^{\ast}$ (Appendix F.4). When the conditional first stage is below the
critical value, a note says so.

**The robustness values.** RV and XRV (Appendix F.2–F.3) measure how strong an omitted confounder of the
*reduced form* would have to be to make $\theta$ insignificant. This is a different threat from a
violation of mediated confounding.

## Speed

How long one estimation takes in each language, on the same simulated data: one instrument,
50 covariates, and a fixed effect with 100 groups. The comparison is meant to help choose a platform;
absolute times depend on the machine. Scripts and full results are in [`benchmarks/`](benchmarks/).

<table>
<thead>
<tr><th rowspan="2">N</th><th colspan="4">Defaults</th><th colspan="4">Defaults + post-double-selection LASSO</th></tr>
<tr><th>Python</th><th>R</th><th>Stata</th><th>MATLAB</th><th>Python</th><th>R</th><th>Stata</th><th>MATLAB</th></tr>
</thead>
<tbody>
<tr><td align="right">1,000</td><td align="right">0.08</td><td align="right">0.11</td><td align="right">0.07</td><td align="right">0.07</td><td align="right">0.08</td><td align="right">0.11</td><td align="right">0.11</td><td align="right">0.05</td></tr>
<tr><td align="right">10,000</td><td align="right">0.74</td><td align="right">1.12</td><td align="right">0.78</td><td align="right">0.32</td><td align="right">2.54</td><td align="right">6.33</td><td align="right">4.37</td><td align="right">1.19</td></tr>
<tr><td align="right">100,000</td><td align="right">5.36</td><td align="right">15.69</td><td align="right">10.94</td><td align="right">1.84</td><td align="right">20.55</td><td align="right">92.86</td><td align="right">127</td><td align="right">10.39</td></tr>
</tbody>
</table>

*Note: Median wall-clock seconds of three runs of the estimation call, after a warm-up run. Loading the
data and printing are excluded.*
- *Defaults: robust standard errors, Anderson–Rubin sets and the κ diagnostics, with the 50 covariates
  and the fixed effect absorbed.*
- *The LASSO selects the controls from the 50 covariates and their squares and cubes (200 candidates).*
- *Machine: Intel Core Ultra 9 285K (24 cores), 64 GB RAM, Windows 11.*
- *Software: Python 3.13 with numpy 2.3 (Intel MKL); R 4.5.2 (default reference BLAS); StataNow 19.5 SE
  (one core); MATLAB R2026b (Intel MKL, 24 threads).*
- *Point estimates: all four languages give the same estimates, to a relative difference of $`10^{-12}`$.
  The true values are $`\theta = 0.5`$ and $`\theta\gamma = 0.2`$.*
  - *Defaults: 0.882, 0.521 and 0.486 for $`\theta`$, and 0.417, 0.205 and 0.196 for $`\theta\gamma`$
    (N = 1,000, 10,000 and 100,000).*
  - *LASSO: 1.197, 0.521 and 0.486 for $`\theta`$, and 0.515, 0.205 and 0.195 for $`\theta\gamma`$.*
- *Why the LASSO differs at N = 1,000: in this design each of the 50 covariates has a small effect, and
  with 1,000 observations the LASSO keeps only 8 controls. From N = 10,000 on, it keeps all 50, and the
  estimates are close to the defaults.*

What the table shows:
- From N = 10,000 on, MATLAB is the fastest, followed by Python. Both use Intel's multithreaded linear
  algebra.
- R on Windows uses its single-threaded reference BLAS by default and is slowest without the LASSO.
- Stata/SE runs on one core and is slowest with the LASSO at N = 100,000.
- At N = 1,000, every run takes about a tenth of a second.

## Citation

If you use `ivmediate`, please cite the paper:

> Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo Pinto (2026). "Identifying Causal Mediation
> Effects with a Single Instrument." NBER Working Paper No. [number]. doi:[DOI]

```bibtex
@techreport{FerraraGoldHeblichPinto2026,
  author      = {Ferrara, Andreas and Gold, Robert and Heblich, Stephan and Pinto, Rodrigo},
  title       = {Identifying Causal Mediation Effects with a Single Instrument},
  institution = {National Bureau of Economic Research},
  type        = {Working Paper},
  number      = {[number]},
  year        = {2026},
  doi         = {[DOI]}
}
```

GitHub's "Cite this repository" button gives the same reference.

## License and contact

MIT License, see [LICENSE](LICENSE). Copyright (c) 2026 Andreas Ferrara, Robert Gold, Stephan Heblich,
Rodrigo Pinto. Maintainer: Andreas Ferrara (a.ferrara@pitt.edu). Please report problems through the
[issue tracker](https://github.com/aferrara87/ivmediate/issues).

**Use of AI.** This page and the implementations in Python, R, Stata and MATLAB were developed and tested
with Claude Fable 5.1 and Claude Opus 5.5 (Anthropic).
