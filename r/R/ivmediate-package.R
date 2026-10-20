#' ivmediate: causal mediation analysis with a single instrument
#'
#' Estimates how much of the effect of a treatment on an outcome runs through a
#' mediator when one instrument for the treatment is available, following
#' Ferrara, Gold, Heblich and Pinto (2026). The estimates come with
#' classical, robust, clustered or spatial (Conley) standard errors, absorbed
#' fixed effects and weights, Olea-Pflueger first-stage diagnostics,
#' Anderson-Rubin confidence sets, kappa- and R^2-based sensitivity analyses,
#' an optional post-double-selection LASSO and bootstrap inference.
#'
#' The public functions are [ivmediate()] (the analysis), its print method
#' [print.ivmediate()], [ivmediate_bootstrap()] (bootstrap inference) and its
#' print method [print.ivmediate_bootstrap()]. The same estimator is available
#' in Python, MATLAB and Stata, with the same argument names; the four
#' implementations agree to rounding error on a common validation harness.
#'
#' @section Example data:
#' `system.file("extdata", "becker_woessmann_2009.csv", package = "ivmediate")`
#' is the county-level data of Becker and Woessmann (2009) for Prussia in 1871
#' (public use), with which the paper's Table 6, column 1 is estimated; see the
#' examples in [ivmediate()].
#'
#' @references
#' Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo Pinto (2026),
#' "Identifying Causal Mediation Effects with a Single Instrument," NBER
#' Working Paper No. \[number\], doi:\[DOI\].
#'
#' Becker, Sascha O., and Ludger Woessmann (2009), "Was Weber Wrong? A Human
#' Capital Theory of Protestant Economic History," Quarterly Journal of
#' Economics 124(2), 531-596.
#'
#' @aliases ivmediate-package
"_PACKAGE"
