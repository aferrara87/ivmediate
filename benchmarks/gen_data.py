"""Simulated data for the speed comparison (benchmarks/README.md).

Writes data/bench_<N>.csv for N = 1,000, 10,000 and 100,000: one instrument z,
treatment t, mediator m, outcome y, 50 covariates x1..x50 and a fixed effect g
with 100 equal groups.

    T = 0.5 Z + X'b + a_g + e_T
    M = 0.4 T + X'b + a_g + e_M
    Y = 0.3 T + 0.5 M + X'b + a_g + e_Y

with b_j = 0.1, a_g ~ N(0, 1), Z and X_j ~ N(0, 1), and (e_T, e_M, e_Y) normal
with unit variances, corr(e_T, e_M) = 0.5 (mediated relevance),
corr(e_T, e_Y) = 0 (mediated confounding) and corr(e_M, e_Y) = 0.3.

Usage:  python gen_data.py            (all three sizes)
        python gen_data.py 1000       (selected sizes)
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SIZES = [int(a) for a in sys.argv[1:]] or [1_000, 10_000, 100_000]
SEED = 20261017
N_X, N_GROUPS = 50, 100
COV = np.array([[1.0, 0.5, 0.0], [0.5, 1.0, 0.3], [0.0, 0.3, 1.0]])


def simulate(n, seed=SEED):
    rng = np.random.default_rng([seed, n])
    X = rng.standard_normal((n, N_X))
    g = np.arange(n) % N_GROUPS + 1
    a = rng.standard_normal(N_GROUPS)[g - 1]
    z = rng.standard_normal(n)
    e = rng.multivariate_normal(np.zeros(3), COV, size=n)
    xb = X @ np.full(N_X, 0.1)
    t = 0.5 * z + xb + a + e[:, 0]
    m = 0.4 * t + xb + a + e[:, 1]
    y = 0.3 * t + 0.5 * m + xb + a + e[:, 2]
    df = pd.DataFrame({"y": y, "t": t, "m": m, "z": z, "g": g})
    for j in range(N_X):
        df[f"x{j + 1}"] = X[:, j]
    return df


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    for n in SIZES:
        path = os.path.join(HERE, "data", f"bench_{n}.csv")
        simulate(n).to_csv(path, index=False, float_format="%.17g")
        print(f"wrote {path}")
