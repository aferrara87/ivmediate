# Cross-language validation

`ivmediate` is implemented four times. This folder checks that the R, Stata and MATLAB versions reproduce
the Python version, which is the reference, on 75 configurations. Numbers must agree to a relative
difference of 1e-8.

## What is compared

Each configuration runs the estimation, or the bootstrap, and compares every reported number: estimates,
standard errors, confidence intervals, F statistics, Olea–Pflueger critical values, Anderson–Rubin sets,
κ and RV/XRV sensitivity, and LASSO selections. The configurations cover:

- **Standard errors:** classical, robust (HC1), cluster (CR1) and Conley (uniform and Bartlett kernels;
  haversine and euclidean distances).
- **Weights:** analytic and probability weights, including zero weights.
- **Fixed effects:**
  - one- and two-way fixed effects;
  - fixed effects nested in clusters;
  - disconnected two-way fixed effects;
  - a weakly connected worker–firm design;
  - singleton groups.
- **Controls:** categorical dummies (text, integer, float and with blanks), the different ways of listing
  controls (lists, strings, globs), and collinear and ill-conditioned controls.
- **Instruments:** one or two, homo- and heteroskedastic errors, a negative θ̂, and weak instruments with
  unbounded Anderson–Rubin sets.
- **LASSO:** the plug-in penalty, a fixed penalty, kept controls, per-variable expansion rules,
  interactions with dummies, fixed effects, clusters, weights, Conley loadings and two instruments, plus
  invariance to shifts in the variables.
- **Bootstrap:** pairs and block, basic, percentile and studentized intervals, missing values, and fixed
  effects nested in the clusters. Every language replays the same draws, so the replicates match one for
  one.
- **Degenerate inputs:** a negative Conley variance, and 13 inputs that must stop with an error of the
  same kind in every language.
- **The Olea–Pflueger table:** all 120 tabulated critical values, and critical values computed off the
  table.

## Running it

From this folder:

```sh
python gen_data.py          # optional: regenerates the data, configs.json and fn_expected.json (all included)
python run_py.py            # Python results -> out_py.json, plus the bootstrap draws bootidx_*.csv
Rscript run_r.R             # -> out_r.json
matlab -batch "run_matlab"  # -> out_matlab.json
python compare.py r         # compares out_r.json with out_py.json; exit code 0 if everything agrees
python compare.py matlab
python compare.py stata
```

For Stata, run `do run_stata.do` from the `stata/` subfolder; it writes `out_stata.json` here. The
runners use the packages of this repository (`../python/src`, `../r/R`, `../matlab`, `../stata`), not
installed copies. R also needs the package `jsonlite`. `stress_r.R` and `glob_check.R` are additional R
checks.

## The comparison

- **Format.** Each runner writes one flat JSON object per configuration, with every value as text:
  numbers to 17 digits, lists as `key.1`, `key.2`, …
- **Tolerance.** Numbers agree within a relative difference of 1e-8 or an absolute one of 1e-10, and text
  must match exactly.
- **Errors.** A configuration that stops with an error must stop in every language with an error of the
  same class (`configs.json`, `error_classes`); the messages may differ.
- **Bootstrap draws.** R replays the draws that `run_py.py` writes here. MATLAB and Stata replay the copies
  in `reference_2026-10/`, which also holds the Python reference output of version 2.0.0.
- **Expected result.** Python, R and MATLAB agree on all 76 entries (75 configurations and the static
  critical-value table). Stata agrees on 74: two configurations cannot be expressed in Stata, an infinite
  value and a per-variable LASSO expansion rule.

## Data

| File | n | Purpose |
|---|---|---|
| `data_main.csv`, `data_weak.csv` | 900, 400 | base designs, and a weak instrument |
| `data_ext.csv` | 1,200 | the remaining scenarios: weights, fixed effects and clusters of every kind, text ids, LASSO candidates, shifted copies, ill-conditioned and degenerate inputs |
| `data_negvar.csv` | 144 | a lattice on which the Conley variance is negative |
| `data_chain.csv` | 4,000 | a weakly connected worker–firm design for two-way fixed effects |

## Known differences in labels

- **Dummy labels.** When pandas reads a numeric category column as floats, Python names the dummies like
  `_D_v_2.0`, where R names them `_D_v_2`. The values agree.
- **Glob matching.** It is case-insensitive in Python on Windows, and case-sensitive elsewhere and in the
  other languages.
