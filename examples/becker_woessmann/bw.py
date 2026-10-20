"""Becker and Woessmann (2009): IV mediation estimates (paper, Table 6, column 1).

Install the package first (pip install ./python from the repository root), then
run from anywhere:  python -B bw.py
"""
import os
import sys

import pandas as pd

try:
    import ivmediate
except ImportError:  # running from a checkout without installing the package
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python", "src"))
    import ivmediate

HERE = os.path.dirname(os.path.abspath(__file__))
bw = pd.read_csv(os.path.join(HERE, "becker_woessmann_2009.csv"))
r = ivmediate.ivmediate(
    bw, Y="inctaxpc", T="f_prot", M="f_rw", Z="kmwittenberg",
    X="f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop",
    se="robust",
)
ivmediate.print_results(r, title="Becker and Woessmann (2009)")
