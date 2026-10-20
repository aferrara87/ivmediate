"""Speed comparison, Python (benchmarks/README.md).

Times ivmediate() on data/bench_<N>.csv with the defaults and with the
post-double-selection LASSO: one untimed warm-up, then 3 timed runs (only the
estimation call; loading the data and printing are excluded). Writes
results/python.csv and results/env_python.txt.

Usage:  python run_python.py            (all sizes)
        python run_python.py 1000       (selected sizes)
"""
import contextlib
import io
import os
import platform
import statistics
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
try:
    import ivmediate
except ImportError:  # a checkout without the installed package
    sys.path.insert(0, os.path.join(HERE, "..", "python", "src"))
    import ivmediate

SIZES = [int(a) for a in sys.argv[1:]] or [1_000, 10_000, 100_000]
X = " ".join(f"x{j}" for j in range(1, 51))
SPECS = {"default": {}, "lasso": {"lasso": True}}
REPS = 3


def fit(df, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return ivmediate.ivmediate(df, Y="y", T="t", M="m", Z="z", X=X, fe="g", **kw)


rows = []
for n in SIZES:
    df = pd.read_csv(os.path.join(HERE, "data", f"bench_{n}.csv"), float_precision="round_trip")
    for spec, kw in SPECS.items():
        fit(df, **kw)                                   # warm-up
        times = []
        for _ in range(REPS):
            t0 = time.perf_counter()
            r = fit(df, **kw)
            times.append(time.perf_counter() - t0)
        rows.append(dict(language="Python", n=n, spec=spec,
                         **{f"t{i + 1}": t for i, t in enumerate(times)},
                         median=statistics.median(times),
                         theta=r["theta_hat"], gamma=r["gamma_hat"], tau=r["tau_hat"],
                         total=r["total_hat"], se_type=r["se_type"]))
        print(f"Python n={n:>6} {spec:7s} median {statistics.median(times):8.3f} s", flush=True)

os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
sizes = "_".join(str(n) for n in SIZES) if sys.argv[1:] else ""
pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", f"python{sizes and '_' + sizes}.csv"),
                          index=False, float_format="%.17g")
blas = np.show_config(mode="dicts")["Build Dependencies"]["blas"]
with open(os.path.join(HERE, "results", "env_python.txt"), "w") as f:
    f.write(f"Python {platform.python_version()}; ivmediate {ivmediate.__version__}; "
            f"numpy {np.__version__} ({blas.get('name')} {blas.get('version')}); "
            f"pandas {pd.__version__}\n")
