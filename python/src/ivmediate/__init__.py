"""ivmediate: IV mediation analysis with a single instrument.

Implements the estimator of Ferrara, Gold, Heblich and Pinto, "Identifying
Causal Mediation Effects with a Single Instrument", with weak-instrument
diagnostics, Anderson-Rubin confidence sets, sensitivity analysis, fixed
effects, weights, LASSO controls and classical, robust, cluster and spatial
(Conley) standard errors.

Usage::

    import ivmediate
    res = ivmediate.ivmediate(df, Y="y", T="t", M="m", Z="z", X="x1 x2")
    ivmediate.print_results(res)
    boot = ivmediate.bootstrap(df, Y="y", T="t", M="m", Z="z", X="x1 x2")
    ivmediate.print_bootstrap(boot)

Authors: Andreas Ferrara (maintainer, a.ferrara@pitt.edu), Robert Gold,
Stephan Heblich, Rodrigo Pinto. License: MIT.

Citation: Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo Pinto
(2026), "Identifying Causal Mediation Effects with a Single Instrument,"
NBER Working Paper No. [number], doi:[DOI].
"""
from ._core import __version__, ivmediate, print_results, bootstrap, print_bootstrap

__all__ = ["ivmediate", "print_results", "bootstrap", "print_bootstrap", "__version__"]
