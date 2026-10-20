"""Smoke test: Becker and Woessmann (2009), the paper's Table 6, column 1.

Runs under pytest (``pytest python/tests``) and as a plain script
(``python python/tests/test_bw.py``). The data are the Becker-Woessmann extract in
``examples/becker_woessmann/`` of the repository.
"""
import os
import sys

import pandas as pd

try:
    import ivmediate
except ImportError:  # running from a checkout without installing the package
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
    import ivmediate

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "examples", "becker_woessmann", "becker_woessmann_2009.csv")
X = "f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop"

# Table 6, column 1, as printed in the paper (robust SEs, one instrument)
EXPECTED = [
    ("n", 426, 0),
    ("total_hat", 0.586, 3), ("se_total", 0.242, 3),
    ("tau_hat", -0.299, 3), ("se_tau", 0.214, 3),
    ("theta_hat", 4.709, 3), ("se_theta", 2.164, 3),
    ("gamma_hat", 0.188, 3), ("se_gamma", 0.027, 3),
    ("indirect", 0.886, 3), ("se_indirect", 0.408, 3),
    ("F_TZ", 75.69, 2), ("F_MZ_T", 16.36, 2),
    ("rho_TY_star", -0.257, 3),
]
AR_THETA = (0.43, 10.20)


def run_bw():
    data = pd.read_csv(DATA)
    return ivmediate.ivmediate(data, Y="inctaxpc", T="f_prot", M="f_rw", Z="kmwittenberg", X=X)


def test_bw_table6_col1():
    res = run_bw()
    for key, value, digits in EXPECTED:
        assert round(float(res[key]), digits) == value, (key, res[key], value)
    ar = res["ar_theta"]
    assert ar["shape"] == "interval"
    assert (round(ar["lower"], 2), round(ar["upper"], 2)) == AR_THETA, (ar["lower"], ar["upper"])


def test_version():
    assert ivmediate.__version__ == "2.0.0"


if __name__ == "__main__":
    test_bw_table6_col1()
    test_version()
    print(f"ivmediate {ivmediate.__version__}: Becker-Woessmann Table 6, column 1 reproduced "
          f"({len(EXPECTED) + 2} values)")
