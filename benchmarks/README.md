# Speed comparison across languages

This folder reproduces the run-time table on the main page. It compares the
four implementations of `ivmediate` on the same simulated data. The goal is a
relative comparison, to help choose a platform; the absolute times depend on
the machine.

## Design

- **Data** (`gen_data.py`): N = 1,000, 10,000 and 100,000. One instrument, a
  treatment, a mediator and an outcome, with 50 covariates and a fixed effect
  with 100 groups. The data-generating process is described at the top of the
  script. The CSV files go to `data/` and are not part of the repository.
- **Specifications:**
  - **Defaults:** robust standard errors, Anderson–Rubin sets and the κ
    diagnostics, with the 50 covariates and the fixed effect absorbed.
  - **Defaults plus the post-double-selection LASSO:** the controls are selected
    from the 50 covariates and their squares and cubes (200 candidates,
    plug-in penalty). Everything else is unchanged.
- **Timing:** each runner times only the estimation call, excluding loading the
  data and printing. It makes one untimed warm-up run, then three timed runs,
  and reports the median.
- **Check:** `make_table.py` stops if the four languages do not give the same
  point estimates (relative difference below 1e-8).

## Running it

From this folder, one language at a time, with the machine otherwise idle:

```sh
python gen_data.py
python run_python.py
Rscript run_r.R
stata-se -b do run_stata.do           # or: do run_stata.do, from Stata
matlab -batch "run_matlab"
python make_table.py                  # writes results/table.md
```

Each runner uses the installed `ivmediate` package if there is one, and the
package in this repository otherwise. Each runner writes
`results/<language>.csv`, with all three times and the point estimates, and
`results/env_<language>.txt`, with the program versions.
