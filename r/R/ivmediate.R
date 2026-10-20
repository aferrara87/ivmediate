# =============================================================================
# ivmediate 2.0.0 -- IV mediation analysis with a single instrument (R)
# =============================================================================
#
# The IV mediation estimator of Ferrara, Gold, Heblich and Pinto, "Identifying
# Causal Mediation Effects with a Single Instrument": estimates, diagnostics,
# sensitivity analysis, classical/robust/clustered/spatial standard errors,
# absorbed fixed effects, weights, the bootstrap and the post-double-selection
# LASSO. Public functions: ivmediate(), its print method, ivmediate_bootstrap()
# and its print method; everything else is internal.
#
# Copyright (c) 2026 Andreas Ferrara, Robert Gold, Stephan Heblich, Rodrigo
# Pinto. MIT license (see the LICENSE file of the package).
#
# Requires R (>= 4.2 for UTF-8 output of the Greek labels on Windows), the
# 'Matrix' package (sparse kernels for spatial SEs and the two-way FE
# projection) and the base package 'parallel' (the bootstrap's random-number
# streams).
#
# Conventions
# -----------
# - Wald CIs use normal critical values z_{1-alpha/2}; AR sets use chi^2_K.
# - Degrees of freedom: absorbed FE count as the rank of the FE dummies
#   (G - 1 one-way; G1 + G2 - c - 1 two-way, c connected components). With
#   cluster SEs, FE dimensions nested in the clusters are not counted in the
#   CR1 correction (reghdfe/fixest convention; with one of two dimensions
#   nested, the other counts G - 1); HC1 and classical SEs count all absorbed
#   FE.
# - Two-way FE are absorbed by an exact projection when the FE dimension with
#   fewer levels has at most 3,000 levels, otherwise by alternating
#   projections with Irons-Tuck acceleration.
# - Olea-Pflueger critical values: simplified test at the estimated effective
#   degrees of freedom K_eff (Montiel Olea and Pflueger 2013, eq. 4); K_eff = 1
#   for one instrument (23.109 at tau = 10%, alpha = 5%).
# - RV/XRV with K >= 2 instruments use the first instrument.
# - RV/XRV use the residual df of the OLS regressions behind theta-hat
#   (sensemakr's dof); their t statistics use the SE regime of the analysis.
#
# R specifics
# -----------
# - data.table and tibble inputs are converted to data.frame.
# - Text id columns (cluster, fe, dummies): "", whitespace-only and NA are
#   missing.
# - Factor or character columns in a numeric role are converted through their
#   labels when every non-missing label is a number (with a note); otherwise an
#   error names the column.
# - Returned scalars are unnamed numerics.
# - ivmediate_bootstrap() never changes the caller's random-number state: the
#   replicates use L'Ecuyer-CMRG streams (parallel::nextRNGStream) from one
#   seed, and the caller's RNG kind and .Random.seed are restored on exit.
#
# Not available in R: lasso_penalty = "cv" (it depends on scikit-learn's fold
# assignment, which cannot be reproduced in R); it stops with a clear message.
# Use the default plug-in penalty or lasso_lambda_override.
#
# Validation: validation/run_r.R in the repository runs every harness config
# and validation/compare.py compares the output with the Python implementation
# (relative 1e-8, absolute 1e-10, equal error classes).
#
# Known cosmetic divergences from Python:
#   - numeric dummy-category labels: "_D_v_2" here vs "_D_v_2.0" for float
#     columns in pandas (indicator values and counts agree); use the R names
#     in lasso_keep;
#   - glob matching is case-sensitive (Python's fnmatch is case-insensitive
#     on Windows only).
# =============================================================================

IVMED_VERSION <- "2.0.0"


# =============================================================================
#  Hajek-form weighted moments
# =============================================================================

# Notes and warnings go through message(), so suppressMessages() silences them; the text is the same
# as with cat().
.inform <- function(..., sep = " ") message(paste(..., sep = sep), appendLF = FALSE)

.weighted_cov <- function(x, y, w) {
  wsum <- sum(w)
  mx <- sum(w * x) / wsum
  my <- sum(w * y) / wsum
  sum(w * (x - mx) * (y - my)) / wsum
}

.weighted_var <- function(x, w) .weighted_cov(x, x, w)

.weighted_std <- function(x, w) sqrt(.weighted_var(x, w))


# =============================================================================
#  Numerically stable least squares (WP1.5)
# =============================================================================

# QR of the equilibrated design A diag(scale) (columns scaled to unit norm).
# Householder QR is invariant to positive column scaling, so A itself is
# factorized and the columns of R are scaled: the equilibrated factorization
# is (Q, R diag(scale)). Every regression is solved from this factorization
# instead of inverting X'X. Singular (or numerically singular: min |R_jj| <=
# 1e-12 * max |R_jj| of the equilibrated R) designs stop with an error. The
# thin Q (n x k) is formed only if need_Q; otherwise the Householder
# reflections in $qr are applied directly (qr.qty, qr.resid, qr.fitted).
.eq_qr <- function(A, need_Q = TRUE) {
  A <- as.matrix(A)
  storage.mode(A) <- "double"
  k <- ncol(A)
  norms <- sqrt(colSums(A * A))
  scale <- ifelse(norms > 0, 1 / ifelse(norms > 0, norms, 1), 1)
  q <- qr(A, tol = 1e-14)
  R <- qr.R(q) * rep(scale, each = k)
  d <- abs(diag(R))
  if (q$rank < k || any(q$pivot != seq_len(k)) ||
      (length(d) > 0L && (!all(is.finite(d)) || min(d) <= 1e-12 * max(max(d), 1e-300))))
    stop("Singular matrix")
  list(qr = q, Q = if (need_Q) qr.Q(q) else NULL, R = R, scale = scale)
}

# Residuals of A (n x p) after projecting on the columns of B.
.residualize <- function(A, B) {
  A <- as.matrix(A)
  if (is.null(B) || ncol(B) == 0L) return(A)
  qr.resid(.eq_qr(B, need_Q = FALSE)$qr, A)
}

# Equilibrated QR with a per-analysis memo (an environment): a design that is
# identical to one factorized before in the same analysis reuses that
# factorization (Q is added when first needed).
.eq_qr_memo <- function(A, need_Q = TRUE, memo = NULL) {
  if (!is.null(memo)) {
    for (i in seq_along(memo$entries)) {
      e <- memo$entries[[i]]
      if (identical(e$A, A)) {
        if (need_Q && is.null(e$f$Q)) { e$f$Q <- qr.Q(e$f$qr); memo$entries[[i]] <- e }
        return(e$f)
      }
    }
  }
  f <- .eq_qr(A, need_Q = need_Q)
  if (!is.null(memo)) memo$entries <- c(utils::tail(memo$entries, 5L), list(list(A = A, f = f)))
  f
}

# Weighted least squares on raw data: residuals are on the raw scale. Y may be
# a matrix (one factorization for several outcomes).
.wls_fit <- function(X, Y, w) {
  sw <- sqrt(as.numeric(w))
  f <- .eq_qr(X * sw, need_Q = FALSE)
  qty <- qr.qty(f$qr, as.matrix(Y * sw))[seq_len(ncol(f$R)), , drop = FALSE]
  beta <- f$scale * backsolve(f$R, qty)
  if (is.matrix(Y)) return(list(coefficients = beta, residuals = Y - X %*% beta))
  list(coefficients = as.numeric(beta), residuals = as.numeric(Y - X %*% beta))
}


# =============================================================================
#  Olea-Pflueger (2013) effective F and critical values (WP1.7)
# =============================================================================

# Published Table 1 (alpha = 0.05), kept for reference: rows K = 1..30,
# columns tau in {0.05, 0.10, 0.20, 0.30}. The critical values the analysis
# uses are computed (see .op_critical_value), which reproduces this table to
# its rounding at K_eff = K.
.OP_CRIT_5PCT_WALD <- matrix(c(
  37.42, 23.11, 15.06, 12.05,
  32.32, 19.29, 12.17,  9.57,
  30.13, 17.67, 10.95,  8.53,
  28.85, 16.72, 10.23,  7.92,
  27.98, 16.08,  9.75,  7.51,
  27.35, 15.62,  9.40,  7.21,
  26.86, 15.26,  9.14,  6.98,
  26.47, 14.97,  8.92,  6.80,
  26.15, 14.73,  8.74,  6.65,
  25.87, 14.53,  8.59,  6.52,
  25.64, 14.36,  8.47,  6.41,
  25.44, 14.21,  8.36,  6.32,
  25.26, 14.08,  8.26,  6.24,
  25.10, 13.96,  8.17,  6.16,
  24.96, 13.86,  8.10,  6.10,
  24.83, 13.77,  8.03,  6.04,
  24.71, 13.68,  7.96,  5.99,
  24.60, 13.60,  7.91,  5.94,
  24.50, 13.53,  7.85,  5.89,
  24.41, 13.46,  7.80,  5.85,
  24.33, 13.40,  7.76,  5.81,
  24.25, 13.35,  7.72,  5.78,
  24.18, 13.29,  7.68,  5.74,
  24.11, 13.24,  7.64,  5.71,
  24.05, 13.20,  7.61,  5.68,
  23.98, 13.15,  7.57,  5.66,
  23.93, 13.11,  7.54,  5.63,
  23.87, 13.07,  7.51,  5.61,
  23.82, 13.04,  7.49,  5.58,
  23.77, 13.00,  7.46,  5.56
), nrow = 30, byrow = TRUE,
  dimnames = list(as.character(1:30), c("0.05", "0.10", "0.20", "0.30")))

stopifnot(abs(.OP_CRIT_5PCT_WALD["1", "0.10"] - 23.11) < 1e-9)

# Simplified OP critical value c(K, tau, alpha) for (possibly non-integer) K:
# the non-central chi^2 quantile chi2^-1_{K, x K}(1 - alpha) / K with x = 1/tau.
.op_critical_value <- function(K, tau, alpha) {
  .check_open_unit("tau", tau)
  .check_open_unit("alpha", alpha)
  K <- as.numeric(K)
  x <- 1 / as.numeric(tau)
  as.numeric(stats::qchisq(1 - as.numeric(alpha), K, ncp = x * K) / K)
}

# Olea-Pflueger effective F: (beta_Z' Q_zz beta_Z) / tr(Q_zz V_ZZ).
.op_effective_F <- function(beta_Z, V_ZZ, Q_zz) {
  beta_Z <- as.numeric(beta_Z)
  V_ZZ <- as.matrix(V_ZZ); Q_zz <- as.matrix(Q_zz)
  as.numeric(t(beta_Z) %*% Q_zz %*% beta_Z) / sum(diag(Q_zz %*% V_ZZ))
}

# Positive definite (finite, smallest eigenvalue > 0)?
.is_pd <- function(V) {
  V <- as.matrix(V)
  if (length(V) == 0L || !all(is.finite(V))) return(FALSE)
  min(eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values) > 0
}

# K_eff (OP 2013, eq. 4) and the simplified critical value at K_eff. With
# W2 = Q^{1/2} V Q^{1/2}: K_eff = (tr W2)^2 (1 + 2x) / (tr(W2'W2) + 2x tr(W2) lmax(W2)).
.op_keff_critval <- function(V_ZZ, Q_zz, tau, alpha) {
  V <- as.matrix(V_ZZ); Q <- as.matrix(Q_zz)
  x <- 1 / as.numeric(tau)
  K <- nrow(V)
  if (K == 1L) {
    keff <- 1
  } else {
    if (!(.is_pd(V) && all(is.finite(Q)))) return(c(NaN, NaN))
    e <- eigen((Q + t(Q)) / 2, symmetric = TRUE)
    Qh <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), K) %*% t(e$vectors)
    W2 <- Qh %*% V %*% Qh
    W2 <- (W2 + t(W2)) / 2
    tr <- sum(diag(W2))
    lmax <- max(eigen(W2, symmetric = TRUE, only.values = TRUE)$values)
    den <- sum(diag(W2 %*% W2)) + 2 * x * tr * lmax
    if (!(is.finite(den) && den > 0 && tr > 0)) return(c(NaN, NaN))
    keff <- tr * tr * (1 + 2 * x) / den
  }
  crit <- as.numeric(stats::qchisq(1 - as.numeric(alpha), keff, ncp = x * keff) / keff)
  c(keff, crit)
}


# =============================================================================
#  Option validation (WP5.5)
# =============================================================================

.SE_TYPES <- c("classical", "robust", "cluster", "spatial")

.show_value <- function(v) paste(utils::capture.output(dput(v)), collapse = "")

.check_bool <- function(name, value) {
  if (!(is.logical(value) && length(value) == 1L && !is.na(value)))
    stop(sprintf("%s must be TRUE or FALSE; got %s.", name, .show_value(value)), call. = FALSE)
}

.check_open_unit <- function(name, value) {
  ok <- is.numeric(value) && length(value) == 1L && is.finite(value)
  if (!ok) {
    if (is.numeric(value) && length(value) == 1L)
      stop(sprintf("%s must lie strictly between 0 and 1; got %s.", name, .show_value(value)), call. = FALSE)
    stop(sprintf("%s must be a number in (0, 1); got %s.", name, .show_value(value)), call. = FALSE)
  }
  if (!(value > 0 && value < 1))
    stop(sprintf("%s must lie strictly between 0 and 1; got %s.", name, .show_value(value)), call. = FALSE)
}

.check_choice <- function(name, value, choices) {
  if (!(is.character(value) && length(value) == 1L && !is.na(value) && value %in% choices))
    stop(sprintf("%s must be one of %s; got %s.", name,
                 paste(sprintf("'%s'", choices), collapse = ", "), .show_value(value)), call. = FALSE)
}

.is_whole <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v == round(v)

# One formatter for every confidence-level label: 0.05 -> "95%", 0.025 -> "97.5%".
.level_label <- function(alpha) sprintf("%g%%", round(100 * (1 - as.numeric(alpha)), 10))


# =============================================================================
#  Variable-list parsing (lists / whitespace strings / shell-style globs)
# =============================================================================

# Translate a shell glob to an anchored regex, following Python's
# fnmatch.translate: '*' -> '.*', '?' -> '.', '[seq]' -> character class
# ('[!seq]' negates; a leading '^' is literal; an unterminated '[' is a
# literal bracket). Matching is case-sensitive (Python's fnmatch is
# case-insensitive on Windows only -- a platform quirk we do not replicate).
.glob_to_regex <- function(pat) {
  chars <- strsplit(pat, "")[[1]]
  n <- length(chars); i <- 1L; out <- character(0)
  while (i <= n) {
    ch <- chars[i]
    if (ch == "*") {
      out <- c(out, ".*")
    } else if (ch == "?") {
      out <- c(out, ".")
    } else if (ch == "[") {
      j <- i + 1L
      if (j <= n && chars[j] == "!") j <- j + 1L
      if (j <= n && chars[j] == "]") j <- j + 1L
      while (j <= n && chars[j] != "]") j <- j + 1L
      if (j > n) {
        out <- c(out, "\\[")                       # unterminated: literal '['
      } else {
        inner <- paste(chars[(i + 1L):(j - 1L)], collapse = "")
        if (startsWith(inner, "!")) inner <- paste0("^", substring(inner, 2))
        else if (startsWith(inner, "^")) inner <- paste0("\\", inner)
        out <- c(out, "[", inner, "]")
        i <- j
      }
    } else {
      if (grepl("[][{}()*+?.\\\\^$|]", ch)) out <- c(out, paste0("\\", ch))
      else out <- c(out, ch)
    }
    i <- i + 1L
  }
  paste0("^", paste(out, collapse = ""), "$")
}

.has_glob <- function(tok) grepl("[*?[]", tok)

# Resolve a variable spec (NULL / whitespace string / character vector, globs
# allowed) against the column names. Glob matches never include the columns in
# `exclude` (the role and auxiliary columns); naming a `forbid` column
# explicitly is an error. Returns list(vars, excluded).
.resolve_var_list <- function(spec, df, name = "<arg>", exclude = character(0),
                              forbid = character(0)) {
  excluded <- character(0)
  if (is.null(spec)) return(list(vars = character(0), excluded = excluded))
  if (is.list(spec)) spec <- unlist(spec, use.names = FALSE)
  if (!is.character(spec))
    stop(sprintf("%s: variable names must be strings; got %s.", name, .show_value(spec)), call. = FALSE)
  if (length(spec) == 1L) {
    s <- trimws(spec)
    if (is.na(s) || !nzchar(s)) return(list(vars = character(0), excluded = excluded))
    tokens <- strsplit(s, "\\s+")[[1]]
  } else {
    tokens <- spec
  }
  cols <- colnames(df)
  avail <- paste0(paste(utils::head(cols, 8), collapse = ", "), if (length(cols) > 8) "..." else "")
  resolved <- character(0)
  for (tok in tokens) {
    if (is.na(tok)) stop(sprintf("%s: variable names must not be NA.", name), call. = FALSE)
    if (.has_glob(tok)) {
      matches <- grep(.glob_to_regex(tok), cols, value = TRUE)
      if (length(matches) == 0L)
        stop(sprintf("%s: glob pattern '%s' matched 0 columns. Available columns include: %s",
                     name, tok, avail), call. = FALSE)
      kept <- matches[!(matches %in% exclude)]
      excluded <- c(excluded, setdiff(matches[matches %in% exclude], excluded))
      if (length(kept) == 0L)
        stop(sprintf("%s: glob pattern '%s' matched only role columns (%s), which cannot be controls.",
                     name, tok, paste(matches, collapse = ", ")), call. = FALSE)
      resolved <- c(resolved, kept)
    } else {
      if (!(tok %in% cols))
        stop(sprintf("%s: variable '%s' not found in data. Available columns include: %s",
                     name, tok, avail), call. = FALSE)
      resolved <- c(resolved, tok)
    }
  }
  bad <- unique(resolved[resolved %in% forbid])
  if (length(bad) > 0L)
    stop(sprintf("%s must not contain the outcome, treatment, mediator or instruments: %s.",
                 name, paste(bad, collapse = ", ")), call. = FALSE)
  all_na <- resolved[vapply(resolved, function(v) all(is.na(df[[v]])), logical(1))]
  if (length(all_na) > 0L)
    stop(sprintf("%s: column(s) %s contain only missing values. Either drop them from %s or fill the missing values before calling ivmediate.",
                 name, paste(unique(all_na), collapse = ", "), name), call. = FALSE)
  list(vars = unique(resolved), excluded = excluded)
}

.summarize_resolved <- function(name, resolved, max_show = 8L) {
  if (length(resolved) == 0L) return(NULL)
  n <- length(resolved)
  body <- if (n <= max_show) paste(resolved, collapse = ", ")
          else paste0(paste(resolved[1:max_show], collapse = ", "), sprintf(", ... (+%d more)", n - max_show))
  sprintf("  Resolved %s (%d var%s): %s", name, n, if (n != 1L) "s" else "", body)
}

# Instruments: one name or a character vector of names, without duplicates.
.normalize_Z <- function(Z) {
  if (is.list(Z)) Z <- unlist(Z, use.names = FALSE)
  if (!is.character(Z))
    stop(sprintf("Z must be a column name or a character vector of names; got %s.", .show_value(Z)), call. = FALSE)
  if (length(Z) == 0L) stop("Z must name at least one instrument.", call. = FALSE)
  if (anyNA(Z)) stop("instrument names must not be NA.", call. = FALSE)
  dup <- unique(Z[duplicated(Z)])
  if (length(dup) > 0L)
    stop(sprintf("duplicate instrument names: %s.", paste(dup, collapse = ", ")), call. = FALSE)
  Z
}

# Text ids: "", whitespace-only and NA are missing (WP5.1).
.missing_text_to_na <- function(v) {
  if (is.character(v) || is.factor(v)) {
    blank <- !is.na(v) & trimws(as.character(v)) == ""
    if (any(blank)) v[blank] <- NA
  }
  v
}

# A column in a numeric role: numeric and logical columns pass; factor or
# character columns are converted through their labels when every non-missing
# label is a number (WP5.2), otherwise an error names the column.
.as_numeric_role <- function(v, name) {
  if (is.numeric(v)) return(list(v = as.numeric(v), note = NULL))
  if (is.logical(v)) return(list(v = as.numeric(v), note = NULL))
  if (is.factor(v) || is.character(v)) {
    lab <- trimws(as.character(v))
    lab[!is.na(lab) & lab == ""] <- NA
    num <- suppressWarnings(as.numeric(lab))
    bad <- !is.na(lab) & is.na(num)
    if (any(bad)) {
      ex <- utils::head(unique(lab[bad]), 3)
      stop(sprintf("column '%s' must be numeric; found values such as %s.", name,
                   paste(sprintf("'%s'", ex), collapse = ", ")), call. = FALSE)
    }
    kind <- if (is.factor(v)) "factor" else "character"
    return(list(v = num, note = sprintf("  Note: column '%s' is a %s; converted to numbers through its labels.", name, kind)))
  }
  stop(sprintf("column '%s' must be numeric; found a column of class %s.", name, class(v)[1]), call. = FALSE)
}

# Validate the roles and build the estimation sample: listwise deletion over
# the used columns, +-Inf check, zero weights dropped (before singleton
# detection), clusterid used only with se = "cluster" (D3).
.prepare_sample <- function(data, Y, T, M, Z, X, dummies, fe, se, clusterid, weights, lat, lon) {
  for (nm in c("Y", "T", "M")) {
    v <- get(nm)
    if (!(is.character(v) && length(v) == 1L && !is.na(v)))
      stop(sprintf("%s must be a column name; got %s.", nm, .show_value(v)), call. = FALSE)
  }
  Zl <- .normalize_Z(Z)
  if (length(unique(c(Y, T, M))) < 3L)
    stop("Y, T and M must be three different columns.", call. = FALSE)
  clash <- Zl[Zl %in% c(Y, T, M)]
  if (length(clash) > 0L)
    stop(sprintf("instruments must differ from the outcome, treatment and mediator: %s.",
                 paste(clash, collapse = ", ")), call. = FALSE)
  fe_cols <- NULL
  if (!is.null(fe)) {
    if (!is.character(fe)) stop("fe must be a column name or a character vector of names.", call. = FALSE)
    fe_cols <- as.character(fe)
    if (length(fe_cols) == 0L) fe_cols <- NULL
    else if (length(fe_cols) > 2L) stop("At most two-way fixed effects are supported.", call. = FALSE)
  }
  cl_eff <- if (identical(se, "cluster")) clusterid else NULL
  roles <- c(Y, T, M, Zl)
  need <- c(roles, fe_cols, cl_eff, weights, if (identical(se, "spatial")) c(lat, lon))
  miss_cols <- unique(need[!(need %in% colnames(data))])
  if (length(miss_cols) > 0L)
    stop(sprintf("column(s) not found in data: %s.", paste(miss_cols, collapse = ", ")), call. = FALSE)
  aux <- c(fe_cols, clusterid, weights, lat, lon)
  exclude <- unique(c(roles, aux))
  rx <- .resolve_var_list(X, data, name = "X", exclude = exclude, forbid = roles)
  rd <- .resolve_var_list(dummies, data, name = "dummies", exclude = exclude, forbid = roles)
  notes <- character(0)
  if (length(rx$excluded)) notes <- c(notes, sprintf("  Note: X globs also matched role columns, excluded: %s.", paste(rx$excluded, collapse = ", ")))
  if (length(rd$excluded)) notes <- c(notes, sprintf("  Note: dummies globs also matched role columns, excluded: %s.", paste(rd$excluded, collapse = ", ")))
  X_res <- rx$vars; D_res <- rd$vars

  numeric_roles <- unique(c(roles, X_res, weights, if (identical(se, "spatial")) c(lat, lon)))
  id_cols <- unique(c(fe_cols, cl_eff, D_res))
  used <- unique(c(numeric_roles, id_cols))
  df <- data[, used, drop = FALSE]
  for (cc in id_cols) df[[cc]] <- .missing_text_to_na(df[[cc]])
  for (cc in numeric_roles) {
    cv <- .as_numeric_role(df[[cc]], cc)
    df[[cc]] <- cv$v
    if (!is.null(cv$note)) notes <- c(notes, cv$note)
  }
  mask <- stats::complete.cases(df)
  num <- as.matrix(df[mask, numeric_roles, drop = FALSE])
  inf_cols <- numeric_roles[colSums(is.infinite(num)) > 0]
  if (length(inf_cols) > 0L)
    stop(sprintf("infinite values in column(s): %s.", paste(inf_cols, collapse = ", ")), call. = FALSE)
  n_zero <- 0L
  if (!is.null(weights)) {
    wv <- df[[weights]]
    if (any(wv[mask] < 0)) stop("Weights must be non-negative.", call. = FALSE)
    zero <- mask & wv == 0
    zero[is.na(zero)] <- FALSE
    n_zero <- sum(zero)
    if (n_zero > 0L) {
      if (n_zero == sum(mask))
        stop("All weights are zero \u2014 no observations have positive weight to estimate from.", call. = FALSE)
      mask <- mask & !zero
    }
  }
  df <- df[mask, , drop = FALSE]
  if (nrow(df) == 0L) stop("No observations remaining after dropping missing values.", call. = FALSE)
  list(df = df, mask = mask, Z = Zl, X = X_res, dummies = D_res, fe_cols = fe_cols,
       clusterid = cl_eff, notes = notes, n_zero_weights = n_zero)
}

# Order-respecting greedy Gram-Schmidt: keep a column if its component
# orthogonal to the kept ones (and to keep_first) exceeds tol relative to its
# reference norm. Returns 1-based kept indices and the dropped names.
# The basis is preallocated (unused columns are zero) and each projection is
# a block Gram-Schmidt step, two passes for the tested columns.
.greedy_independent <- function(A, names, tol = 1e-9, keep_first = NULL, ref_norms = NULL) {
  A <- as.matrix(A); n <- nrow(A); p <- ncol(A)
  K0 <- if (!is.null(keep_first)) matrix(as.numeric(keep_first), nrow = n) else matrix(0, n, 0)
  basis <- matrix(0, n, ncol(K0) + p); nb <- 0L
  for (j in seq_len(ncol(K0))) {
    v <- K0[, j] / max(sqrt(sum(K0[, j]^2)), 1e-300)
    if (nb > 0L) v <- as.numeric(v - basis %*% crossprod(basis, v))
    nv <- sqrt(sum(v^2))
    if (nv > tol) { nb <- nb + 1L; basis[, nb] <- v / nv }
  }
  kept <- integer(0); dropped <- character(0)
  for (j in seq_len(p)) {
    ref <- if (!is.null(ref_norms)) as.numeric(ref_norms[j]) else sqrt(sum(A[, j]^2))
    if (!is.finite(ref) || ref <= 0) { dropped <- c(dropped, names[j]); next }
    v <- A[, j] / ref
    if (nb > 0L) for (pass in 1:2) v <- as.numeric(v - basis %*% crossprod(basis, v))
    nv <- sqrt(sum(v^2))
    if (is.finite(nv) && nv > tol) { nb <- nb + 1L; basis[, nb] <- v / nv; kept <- c(kept, j) }
    else dropped <- c(dropped, names[j])
  }
  list(kept = kept, dropped = dropped)
}


# =============================================================================
#  Fixed effects: singletons, degrees of freedom, (weighted) demeaning
# =============================================================================

# Integer codes 1..G in order of first appearance (pandas.factorize, sort=False).
.codes <- function(x) match(x, unique(x))

# Iteratively drop singleton FE groups until none remain (reghdfe convention).
.drop_singletons <- function(df, fe_col) {
  n_before <- nrow(df)
  repeat {
    g <- as.character(df[[fe_col]])
    counts <- table(g)
    singles <- names(counts)[counts == 1L]
    if (length(singles) == 0L) break
    df <- df[!(g %in% singles), , drop = FALSE]
  }
  list(df = df, n_dropped = n_before - nrow(df))
}

# Union-find on the bipartite FE graph: number of connected components.
.n_components <- function(c1, c2, G1, G2) length(unique(.fe_component_labels(c1, c2, G1, G2)))

# Connected-component label (a root node id) of every level of the bipartite FE
# graph: nodes 1..G1 are the levels of the first dimension, G1 + 1..G1 + G2
# those of the second; every observation links its two levels (union-find).
.fe_component_labels <- function(c1, c2, G1, G2) {
  parent <- seq_len(G1 + G2)
  find <- function(i) { while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
  e <- unique(cbind(c1, G1 + c2))
  for (r in seq_len(nrow(e))) {
    a <- find(e[r, 1]); b <- find(e[r, 2])
    if (a != b) parent[a] <- b
  }
  vapply(seq_len(G1 + G2), find, integer(1))
}

# Absorbed df = rank of the FE dummies beyond the intercept (WP1.3), and the df
# used by the CR1 correction (FE nested in the clusters not counted; with one of
# two dimensions nested, the other counts G - 1: reghdfe / fixest / pyfixest).
.fe_degrees_of_freedom <- function(df, fe_cols, cluster_col = NULL) {
  codes <- lapply(fe_cols, function(cc) .codes(df[[cc]]))
  G <- vapply(codes, max, integer(1))
  if (!is.null(cluster_col)) {
    cl <- .codes(df[[cluster_col]])
    nested <- vapply(codes, function(cd) nrow(unique(cbind(cd, cl))) == length(unique(cd)), logical(1))
  } else nested <- rep(FALSE, length(codes))
  if (length(codes) == 1L) {
    d <- G[1] - 1L
    d_cr1 <- if (nested[1]) 0L else d
  } else {
    cc <- .n_components(codes[[1]], codes[[2]], G[1], G[2])
    d <- G[1] + G[2] - cc - 1L
    d_cr1 <- if (nested[1] && nested[2]) 0L else if (nested[1]) G[2] - 1L else if (nested[2]) G[1] - 1L else d
  }
  list(d = as.integer(d), d_cr1 = as.integer(d_cr1), n_groups = as.integer(sum(G)), nested = nested)
}

# Weighted group means of the columns of V, expanded to the rows.
.group_means <- function(codes, w, sw, V) {
  (rowsum(V * w, codes, reorder = TRUE) / sw)[codes, , drop = FALSE]
}

# One-way demeaning with weighted group means (WP1.2).
.demean_within_groups <- function(df, cols, fe_col, weights = NULL) {
  codes <- .codes(df[[fe_col]])
  w <- if (is.null(weights)) rep(1, nrow(df)) else as.numeric(weights)
  sw <- as.numeric(rowsum(w, codes, reorder = TRUE))
  V <- as.matrix(df[, cols, drop = FALSE]); storage.mode(V) <- "double"
  df[, cols] <- V - .group_means(codes, w, sw, V)
  df
}

# Minimum-norm least-squares solution of A X = B (as numpy's lstsq with its
# default rcond: singular values below eps * max(dim) * max(d) count as zero).
.lstsq_min_norm <- function(A, B) {
  s <- svd(A)
  tol <- max(dim(A)) * .Machine$double.eps * max(s$d, 0)
  dinv <- ifelse(s$d > tol, 1 / s$d, 0)
  s$v %*% (dinv * crossprod(s$u, B))
}

# Exact weighted two-way FE residualization (Addendum A.2, O-2; the algorithm of
# Stata's _ivm_demean2_direct and Python's _demean_twoway_direct). Model
# V = alpha_A + beta_B + R by weighted least squares, with A the FE dimension
# with more levels and B the other one. Eliminating alpha gives L beta = r for
# the B effects, with the Schur complement L = diag(s_B) - N' diag(1/s_A) N,
# where s_A and s_B are the weight sums per level, N[a, b] is the weight of the
# observations in level a of A and level b of B, and r holds the B-level sums of
# w * (V - within-A means). L is singular, with one null direction per
# connected component of the FE graph, so one B level per component (the first
# in B code order) is fixed at 0. One Cholesky factorization serves all columns.
# The residual is R = (V - mean_A V) - (beta_B - mean_A beta_B). Returns NULL
# when the smaller dimension has more than max_small levels.
.demean_twoway_direct <- function(V, c1, c2, w, max_small = 3000L) {
  G1 <- max(c1); G2 <- max(c2)
  if (G2 > G1) { cA <- c2; cB <- c1; Ga <- G2; Gb <- G1 } else { cA <- c1; cB <- c2; Ga <- G1; Gb <- G2 }
  if (Gb > max_small) return(NULL)
  sA <- as.numeric(rowsum(w, cA, reorder = TRUE))
  sB <- as.numeric(rowsum(w, cB, reorder = TRUE))
  N <- Matrix::sparseMatrix(i = cA, j = cB, x = w, dims = c(Ga, Gb))      # duplicates are summed
  Ns <- Matrix::Diagonal(x = 1 / sqrt(sA)) %*% N
  L <- diag(sB, nrow = Gb) - as.matrix(Matrix::crossprod(Ns))
  L <- 0.5 * (L + t(L))
  lab <- .fe_component_labels(cA, cB, Ga, Gb)[Ga + seq_len(Gb)]
  keep <- duplicated(lab)                       # the first B level of each component is fixed
  V1 <- V - .group_means(cA, w, sA, V)
  rhs <- rowsum(w * V1, cB, reorder = TRUE)
  beta <- matrix(0, Gb, ncol(V))
  kb <- which(keep)
  if (length(kb) > 0L) {
    Lk <- L[kb, kb, drop = FALSE]; rk <- rhs[kb, , drop = FALSE]
    Rc <- tryCatch(chol(Lk), error = function(e) NULL)
    beta[kb, ] <- if (!is.null(Rc)) backsolve(Rc, backsolve(Rc, rk, transpose = TRUE)) else .lstsq_min_norm(Lk, rk)
  }
  Bv <- beta[cB, , drop = FALSE]
  V1 - (Bv - .group_means(cA, w, sA, Bv))
}

# Two-way demeaning (weighted). When the FE dimension with fewer levels has at
# most max_direct (3,000) levels: the exact projection above (A.2, O-2;
# alternating projections can stall on weakly connected designs such as a
# worker-firm chain). Otherwise alternating weighted projections with
# Irons-Tuck extrapolation, relative criterion max|change| / sd(column) < tol
# (WP1.4), up to max_iter iterations, with a warning only on true
# non-convergence.
.demean_twoway <- function(df, cols, fe_cols, weights = NULL, tol = 1e-12, max_iter = 10000L,
                           max_direct = 3000L) {
  c1 <- .codes(df[[fe_cols[1]]]); c2 <- .codes(df[[fe_cols[2]]])
  w <- if (is.null(weights)) rep(1, nrow(df)) else as.numeric(weights)
  V <- as.matrix(df[, cols, drop = FALSE]); storage.mode(V) <- "double"
  if (ncol(V) == 0L) return(df)
  R <- .demean_twoway_direct(V, c1, c2, w, max_small = max_direct)
  if (!is.null(R)) { df[, cols] <- R; return(df) }
  sw1 <- as.numeric(rowsum(w, c1, reorder = TRUE)); sw2 <- as.numeric(rowsum(w, c2, reorder = TRUE))
  sd <- sqrt(colMeans(sweep(V, 2, colMeans(V))^2))
  active <- sd > 1e-14 * pmax(apply(abs(V), 2, max), 1e-300)
  scale <- pmax(sd, 1e-300)
  sweep2 <- function(A) { A <- A - .group_means(c1, w, sw1, A); A - .group_means(c2, w, sw2, A) }
  converged <- !any(active); change <- 0; it <- 0L
  while (!converged && it < max_iter) {
    F <- sweep2(V); it <- it + 1L
    change <- max(apply(abs(F - V)[, active, drop = FALSE], 2, max) / scale[active])
    if (change < tol) { V <- F; converged <- TRUE; break }
    Gs <- sweep2(F); it <- it + 1L
    D1 <- F - V; D2 <- Gs - F; DD <- D2 - D1
    den <- colSums(DD * DD)
    coef <- ifelse(den > 0, colSums(D2 * DD) / ifelse(den > 0, den, 1), 0)
    V <- Gs - D2 * rep(coef, each = nrow(D2))            # Irons-Tuck extrapolation
  }
  if (!converged)
    .inform(sprintf("  Warning: two-way FE demeaning did not converge after %d iterations (last relative change = %.2e > tol = %.0e). Results may be imprecise; check for nearly-collinear FE dimensions.\n",
                as.integer(max_iter), change, tol))
  df[, cols] <- V
  df
}

# Expand categorical variables into drop-first indicator columns. Categories
# are sorted with method = "radix" (C-locale byte order == Python's codepoint
# sort for text; factors keep their level order).
.expand_dummies <- function(df, dummy_vars) {
  new_cols <- character(0)
  for (var in dummy_vars) {
    cats <- sort(unique(df[[var]][!is.na(df[[var]])]), method = "radix")
    if (length(cats) <= 1L) next
    for (cat in cats[-1L]) {                       # drop first category (reference)
      col_name <- sprintf("_D_%s_%s", var, as.character(cat))
      df[[col_name]] <- as.numeric(df[[var]] == cat)
      new_cols <- c(new_cols, col_name)
    }
  }
  list(df = df, new_cols = new_cols)
}


# =============================================================================
#  Variance-covariance estimators for 2SLS (sandwiches from the QR solves)
# =============================================================================
#
# With X_hat = Q R diag(1/scale), (X_hat'X_hat)^{-1} = Tm Tm' where
# Tm = diag(scale) R^{-1}, and the rotated scores are Q * resid, so every
# sandwich is Tm meat Tm'.

# Number of negative diagonal entries met by the spatial sandwich in one
# analysis (counted, then reported once per analysis).
.new_neg_counter <- function() { e <- new.env(); e$neg_var <- 0L; e }

.vcov_classical <- function(resid, Tm, n, k, df_adjust = 0) {
  sigma2 <- sum(resid^2) / (n - k - df_adjust)
  sigma2 * tcrossprod(Tm)
}

.vcov_robust <- function(resid, Q, Tm, n, k, df_adjust = 0) {
  S <- Q * resid
  meat <- crossprod(S) * (n / (n - k - df_adjust))
  Tm %*% meat %*% t(Tm)
}

# cl = list(codes = integer codes 1..G, G = number of clusters).
.vcov_cluster <- function(resid, Q, Tm, n, k, cl, df_adjust = 0) {
  G <- cl$G
  if (G < 2L) stop("cluster SEs need at least 2 clusters.", call. = FALSE)
  sums <- rowsum(Q * resid, cl$codes, reorder = TRUE)
  correction <- (G / (G - 1)) * ((n - 1) / (n - k - df_adjust))
  meat <- crossprod(sums) * correction
  Tm %*% meat %*% t(Tm)
}

.vcov_spatial <- function(resid, Q, Tm, n, k, A, df_adjust = 0, neg = NULL) {
  if (is.null(A)) stop("internal error: spatial adjacency not built.")
  S <- Q * resid
  meat <- crossprod(S, as.matrix(A %*% S)) * ((n - 1) / max(n - k - df_adjust, 1))
  vcov <- Tm %*% meat %*% t(Tm)
  n_bad <- sum(diag(vcov) < 0)
  if (n_bad > 0L) {
    if (!is.null(neg)) neg$neg_var <- neg$neg_var + n_bad
    else .inform(sprintf("  Warning: spatial vcov has %d negative diagonal entries (not positive semi-definite); the affected SEs are NaN.\n", n_bad))
  }
  vcov
}

# 2SLS with the variance under se in {"none", "classical", "robust", "cluster",
# "spatial"}. Passing the same matrix as X and Z gives OLS. vc = list(cl, A, neg).
# memo: optional per-analysis environment for reusing factorizations.
.twosls_with_vcov <- function(Y, X, Z, se = "robust", vc = NULL, df_adjust = 0, memo = NULL) {
  Y <- as.numeric(Y)
  n <- length(Y); k <- ncol(X)
  if (identical(Z, X)) {
    X_hat <- X
  } else {
    # X_hat = P_Z X. Columns of X that are also columns of Z (the included
    # exogenous regressors) project onto themselves and are copied; candidates
    # are found by their column sums and confirmed exactly.
    qz <- .eq_qr_memo(Z, need_Q = FALSE, memo = memo)$qr
    sx <- colSums(X); sz <- colSums(Z)
    in_Z <- vapply(seq_len(k), function(j) {
      cand <- which(sz == sx[j])
      if (length(cand) == 0L) return(FALSE)
      x <- as.numeric(X[, j])
      for (l in cand) if (identical(x, as.numeric(Z[, l]))) return(TRUE)
      FALSE
    }, logical(1))
    X_hat <- X
    if (!all(in_Z)) X_hat[, !in_Z] <- qr.fitted(qz, X[, !in_Z, drop = FALSE])
  }
  # Q (n x k) is needed only for the sandwich meat and the cross-equation
  # covariance; point estimates use the reflections directly.
  f <- .eq_qr_memo(X_hat, need_Q = !identical(se, "none"), memo = memo)
  beta <- as.numeric(f$scale * backsolve(f$R, qr.qty(f$qr, Y)[seq_len(k)]))
  resid <- as.numeric(Y - X %*% beta)
  Tm <- f$scale * backsolve(f$R, diag(k))            # diag(scale) R^{-1}
  vcov <- switch(se,
    none      = matrix(NaN, k, k),
    classical = .vcov_classical(resid, Tm, n, k, df_adjust),
    robust    = .vcov_robust(resid, f$Q, Tm, n, k, df_adjust),
    cluster   = {
      if (is.null(vc$cl)) stop("cluster_ids required when se='cluster'")
      .vcov_cluster(resid, f$Q, Tm, n, k, vc$cl, df_adjust)
    },
    spatial   = .vcov_spatial(resid, f$Q, Tm, n, k, vc$A, df_adjust, vc$neg),
    stop(sprintf("Unknown se type: %s. Use 'classical', 'robust', 'cluster', or 'spatial'.", se))
  )
  dg <- diag(vcov)
  se_vec <- ifelse(!is.na(dg) & dg >= 0, sqrt(pmax(dg, 0)), NaN)
  list(coef = beta, se = as.numeric(se_vec), vcov = vcov, resid = resid,
       X_hat = X_hat, df = n - k - df_adjust, Q = f$Q, Tm = Tm)
}

# Cross-equation covariance Cov(beta_{1,idx1}, beta_{2,idx2}) of two 2SLS fits
# on the same observations (stacked-GMM sandwich), for the delta-method SE of
# theta * gamma. Small-sample factors use k_eff = max(k1, k2).
.cross_equation_cov <- function(res1, res2, idx1, idx2, se = "robust", vc = NULL, df_adjust = 0) {
  if (identical(se, "none")) return(NaN)
  n <- length(res1$resid)
  k_eff <- max(nrow(res1$Tm), nrow(res2$Tm))
  psi1 <- as.numeric(res1$Q %*% res1$Tm[idx1, ]) * res1$resid
  psi2 <- as.numeric(res2$Q %*% res2$Tm[idx2, ]) * res2$resid
  cov <- switch(se,
    classical = {
      sigma_12 <- sum(res1$resid * res2$resid) / max(n - k_eff - df_adjust, 1)
      cross <- res1$Tm %*% crossprod(res1$Q, res2$Q) %*% t(res2$Tm)
      sigma_12 * cross[idx1, idx2]
    },
    robust = sum(psi1 * psi2) * (n / (n - k_eff - df_adjust)),
    cluster = {
      G <- vc$cl$G
      s1 <- as.numeric(rowsum(psi1, vc$cl$codes, reorder = TRUE))
      s2 <- as.numeric(rowsum(psi2, vc$cl$codes, reorder = TRUE))
      sum(s1 * s2) * (G / (G - 1)) * ((n - 1) / (n - k_eff - df_adjust))
    },
    spatial = sum(psi1 * as.numeric(vc$A %*% psi2)) * ((n - 1) / max(n - k_eff - df_adjust, 1)),
    stop(sprintf("Unknown se type: %s", se))
  )
  as.numeric(cov)
}


# =============================================================================
#  Spatial (Conley 1999) kernel adjacency
# =============================================================================
#
# A[i, j] = K(d_ij) for every pair within `cutoff` (self included), with the
# haversine distance in km (Earth radius 6371) for lat/lon in degrees, or the
# Euclidean distance for projected coordinates. The search is exact but only
# compares points whose first coordinate lies within the cutoff band: the
# great-circle distance is at least R * |dlat| (radians), and the Euclidean
# distance at least |dlat|. Rows are processed in latitude order in blocks of
# about 5e6 distance cells.
.EARTH_RADIUS_KM <- 6371.0

.spatial_adjacency <- function(lat, lon, cutoff, kernel = "uniform", metric = "haversine") {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("se='spatial' needs the 'Matrix' package (ships with R): install.packages('Matrix').")
  n <- length(lat)
  is_hav <- identical(metric, "haversine")
  ord <- order(lat); lat_s <- lat[ord]; lon_s <- lon[ord]
  band <- if (is_hav) cutoff / .EARTH_RADIUS_KM * (180 / pi) else cutoff
  band <- band * (1 + 1e-9) + 1e-12                       # safety margin; distances decide
  if (is_hav) { latr <- lat_s * (pi / 180); lonr <- lon_s * (pi / 180); coslat <- cos(latr) }
  ii <- vector("list", 0); jj <- vector("list", 0); kk <- vector("list", 0)
  start <- 1L
  while (start <= n) {
    lo_j <- findInterval(lat_s[start] - band, lat_s, left.open = TRUE) + 1L
    # rows per block so that rows x candidates stays near 5e6 cells
    span_hi <- findInterval(lat_s[start] + band, lat_s)
    width <- max(span_hi - lo_j + 1L, 1L)
    nb <- max(1L, min(n - start + 1L, as.integer(5e6 / width)))
    idx <- start:(start + nb - 1L)
    hi_j <- findInterval(lat_s[idx[length(idx)]] + band, lat_s)
    cand <- lo_j:hi_j
    if (is_hav) {
      dlat <- outer(latr[idx], latr[cand], "-")
      dlon <- outer(lonr[idx], lonr[cand], "-")
      a <- sin(dlat / 2)^2 + outer(coslat[idx], coslat[cand], "*") * sin(dlon / 2)^2
      a[a < 0] <- 0; a[a > 1] <- 1
      D <- 2 * .EARTH_RADIUS_KM * asin(sqrt(a))
    } else {
      dx <- outer(lat_s[idx], lat_s[cand], "-"); dy <- outer(lon_s[idx], lon_s[cand], "-")
      D <- sqrt(dx * dx + dy * dy)
    }
    D <- matrix(D, nrow = length(idx))
    hit <- which(D <= cutoff, arr.ind = TRUE)
    if (nrow(hit) > 0L) {
      d <- D[hit]
      ii[[length(ii) + 1L]] <- ord[idx[hit[, 1L]]]
      jj[[length(jj) + 1L]] <- ord[cand[hit[, 2L]]]
      kk[[length(kk) + 1L]] <- if (kernel == "uniform") rep(1.0, length(d)) else pmax(1.0 - d / cutoff, 0.0)
    }
    start <- start + nb
  }
  ii <- unlist(ii); jj <- unlist(jj); kk <- unlist(kk)
  list(A = Matrix::sparseMatrix(i = ii, j = jj, x = kk, dims = c(n, n)), nnz = length(ii))
}


# =============================================================================
#  Anderson-Rubin confidence sets (weak-instrument-robust inference, WP4)
# =============================================================================
#
# Confidence set for the coefficient on one endogenous regressor by test
# inversion: under H0: beta = b0 the regression of (y_dep - b0 x_endog) on
# W_incl has zero coefficients on the excluded instruments z_cols; the AR
# statistic is their joint Wald statistic under the analysis' SE regime, and
# the set is {b0 : AR(b0) <= chi2_{K, 1-alpha}}.
#
# The instrument coefficients are linear in b0, b(b0) = b_c - s B (s = b0 - tc),
# and every variance (classical, HC1, CR1, Conley; weights and FE enter through
# the transformed arrays) is quadratic, V(b0) = V_c + s V1 + s^2 V2; both are
# recovered from three refits and checked at a fourth.
#   K = 1: exact Fieller-type inversion of a scalar quadratic.
#   K >= 2: AR(b0) = b' V^{-1} b in closed form on a search grid that only
#     brackets sign changes; boundaries by Brent's method; interior local
#     minima between grid points are checked so narrow intervals are found.
# Under spatial SEs, a negative variance on the computed range makes the set
# "unavailable" (NaN bounds) rather than empty.
.ar_confidence_set <- function(y_dep, x_endog, W_incl, z_cols, se, vc, df_adjust, alpha,
                               center = 0.0, scale = NULL,
                               grid_halfwidth = 200.0, grid_points = 4001L) {
  K <- length(z_cols)
  cv <- as.numeric(stats::qchisq(1 - alpha, K))

  ar_eval <- function(b0) {
    u <- y_dep - b0 * x_endog
    res <- .twosls_with_vcov(u, W_incl, W_incl, se = se, vc = vc, df_adjust = df_adjust)
    list(bz = as.numeric(res$coef[z_cols]), Vz = matrix(res$vcov[z_cols, z_cols], K, K))
  }
  ar_stat <- function(b0) { e <- ar_eval(b0); as.numeric(t(e$bz) %*% solve(e$Vz, e$bz)) }

  e0 <- ar_eval(0)
  if (.is_pd(e0$Vz)) {
    stat0 <- as.numeric(t(e0$bz) %*% solve(e0$Vz, e0$bz))
    pval0 <- as.numeric(stats::pchisq(stat0, K, lower.tail = FALSE))
  } else { stat0 <- NaN; pval0 <- NaN }

  out <- list(shape = NA_character_, lower = NaN, upper = NaN, left_upper = NaN, right_lower = NaN,
              segments = NULL, cv = cv, df_chi2 = K, alpha = alpha,
              stat_at_zero = stat0, pval_at_zero = pval0, contains_zero = NA, method = NA_character_)

  tc <- if (is.finite(center)) center else 0.0
  sc <- if (!is.null(scale) && is.finite(scale) && scale > 0) scale else max(abs(tc), 1.0)
  s_hw <- as.numeric(grid_halfwidth) * sc                # the range the set is computed on
  check_neg <- identical(se, "spatial")
  use_grid <- (K > 1L)

  if (!use_grid) {
    d <- sc
    em <- ar_eval(tc - d); ec <- ar_eval(tc); ep <- ar_eval(tc + d)
    bm <- em$bz[1]; b0c <- ec$bz[1]; bp <- ep$bz[1]
    Vm <- em$Vz[1, 1]; V0 <- ec$Vz[1, 1]; Vp <- ep$Vz[1, 1]
    B <- (bm - bp) / (2 * d)
    bv <- (Vp - Vm) / (2 * d)
    cvv <- (Vp - 2 * V0 + Vm) / (2 * d * d)
    qa <- B * B - cv * cvv
    qb <- -2 * b0c * B - cv * bv
    qc <- b0c * b0c - cv * V0
    echk <- ar_eval(tc + 2 * d)
    q_pred <- qa * (2 * d)^2 + qb * (2 * d) + qc
    q_act <- echk$bz[1]^2 - cv * echk$Vz[1, 1]
    q_scale <- max(abs(qa) * (2 * d)^2, abs(qb) * (2 * d), abs(qc), abs(q_act), 1e-300)
    if (!is.finite(q_pred) || abs(q_pred - q_act) > 1e-6 * q_scale) {
      use_grid <- TRUE  # defensive fallback; should not trigger
    } else {
      out$method <- "fieller"
      s_c <- c(-s_hw, 0, s_hw)
      if (cvv > 0 && abs(bv / (2 * cvv)) <= s_hw) s_c <- c(s_c, -bv / (2 * cvv))
      v_neg <- check_neg && any(V0 + bv * s_c + cvv * s_c * s_c <
                                  -1e-12 * (abs(V0) + abs(bv * s_c) + abs(cvv) * s_c * s_c))
      tiny <- 1e-14 * max(abs(qb) * d, abs(qc), 1e-300)
      if (v_neg) {
        out$shape <- "unavailable"
      } else if (abs(qa) * d * d <= tiny) {
        if (abs(qb) * d <= tiny) {
          out$shape <- if (qc <= 0) "line" else "empty"
        } else if (qb > 0) {
          out$shape <- "rays_left_only"; out$left_upper <- tc + (-qc / qb)
        } else {
          out$shape <- "rays_right_only"; out$right_lower <- tc + (-qc / qb)
        }
      } else {
        disc <- qb * qb - 4 * qa * qc
        if (disc < 0) {
          out$shape <- if (qa > 0) "empty" else "line"
        } else {
          sq <- sqrt(disc)
          q_ <- if (qb != 0) -0.5 * (qb + sign(qb) * sq) else 0.5 * sq
          r1 <- q_ / qa
          r2 <- if (q_ != 0) qc / q_ else -qb / qa - r1
          s_lo <- min(r1, r2); s_hi <- max(r1, r2)
          if (qa > 0) {
            out$shape <- "interval"; out$lower <- tc + s_lo; out$upper <- tc + s_hi
          } else {
            out$shape <- "rays"; out$left_upper <- tc + s_lo; out$right_lower <- tc + s_hi
          }
        }
      }
    }
  }

  if (use_grid) {
    out$method <- "grid"
    grid <- tc + seq(-grid_halfwidth, grid_halfwidth, length.out = as.integer(grid_points)) * sc
    n_g <- length(grid)
    d <- sc
    em <- ar_eval(tc - d); ec <- ar_eval(tc); ep <- ar_eval(tc + d)
    Bv <- (em$bz - ep$bz) / (2 * d)
    V1 <- (ep$Vz - em$Vz) / (2 * d)
    V2 <- (ep$Vz - 2 * ec$Vz + em$Vz) / (2 * d * d)
    closed_one <- function(b0) {
      s <- b0 - tc
      bb <- ec$bz - s * Bv
      VV <- ec$Vz + s * V1 + s * s * V2
      r <- tryCatch(as.numeric(t(bb) %*% solve(VV, bb)), error = function(e) NA_real_)
      if (is.na(r)) r <- tryCatch(ar_stat(b0), error = function(e) NaN)
      r
    }
    chk_b <- tc + 2 * d
    chk_act <- tryCatch(ar_stat(chk_b), error = function(e) NaN)
    chk_pred <- closed_one(chk_b)
    exact <- all(is.finite(c(chk_act, chk_pred))) && abs(chk_pred - chk_act) <= 1e-6 * max(abs(chk_act), 1)
    if (exact) {
      stat_grid <- vapply(grid, closed_one, numeric(1))
      f_ <- function(b) closed_one(b) - cv
      Vfun <- function(b0) { s <- b0 - tc; ec$Vz + s * V1 + s * s * V2 }
    } else {
      evals <- lapply(grid, ar_eval)
      stat_grid <- vapply(evals, function(e) if (.is_pd(e$Vz)) as.numeric(t(e$bz) %*% solve(e$Vz, e$bz)) else NaN, numeric(1))
      f_ <- function(b) ar_stat(b) - cv
      Vfun <- NULL
    }
    unavailable <- FALSE
    if (check_neg) {
      ev <- if (exact) lapply(grid, function(b0) { V <- Vfun(b0); eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values })
            else lapply(evals, function(e) eigen((e$Vz + t(e$Vz)) / 2, symmetric = TRUE, only.values = TRUE)$values)
      lam_lo <- vapply(ev, min, numeric(1)); lam_sc <- vapply(ev, function(x) max(abs(x)), numeric(1))
      unavailable <- any(!is.finite(lam_lo)) || any(lam_lo < -1e-12 * lam_sc)
    }
    accept <- !is.na(stat_grid) & stat_grid <= cv
    root <- function(a_, b_) {
      fa <- f_(a_); fb <- f_(b_)
      if (!(is.finite(fa) && is.finite(fb)) || fa * fb > 0) return(NULL)
      xtol <- 1e-14 * max(abs(a_), abs(b_), sc)
      stats::uniroot(f_, c(a_, b_), f.lower = fa, f.upper = fb, tol = xtol, maxiter = 1000L)$root
    }
    if (unavailable) {
      out$shape <- "unavailable"
    } else {
      segments <- list()
      i <- 1L
      while (i <= n_g) {
        if (accept[i]) {
          j <- i
          while (j + 1L <= n_g && accept[j + 1L]) j <- j + 1L
          lo <- grid[i]; hi <- grid[j]
          if (i > 1L) { r_ <- root(grid[i - 1L], grid[i]); if (!is.null(r_)) lo <- r_ }
          if (j < n_g) { r_ <- root(grid[j], grid[j + 1L]); if (!is.null(r_)) hi <- r_ }
          segments[[length(segments) + 1L]] <- list(lo = lo, hi = hi, open_left = (i == 1L), open_right = (j == n_g))
          i <- j + 1L
        } else i <- i + 1L
      }
      if (n_g >= 3L) {
        ii <- 2:(n_g - 1L)
        is_min <- (stat_grid[ii] <= stat_grid[ii - 1L]) & (stat_grid[ii] <= stat_grid[ii + 1L]) &
          !accept[ii] & is.finite(stat_grid[ii])
        is_min[is.na(is_min)] <- FALSE
        cand <- ii[is_min]
        cand <- utils::head(cand[order(stat_grid[cand], method = "radix")], 50)
        for (i_ in cand) {
          a_ <- grid[i_ - 1L]; b_ <- grid[i_ + 1L]
          xatol <- 1e-12 * max(abs(a_), abs(b_), sc)
          om <- stats::optimize(f_, c(a_, b_), tol = xatol)
          if (is.finite(om$objective) && om$objective < 0) {
            lo_ <- root(a_, om$minimum); hi_ <- root(om$minimum, b_)
            if (!is.null(lo_) && !is.null(hi_))
              segments[[length(segments) + 1L]] <- list(lo = lo_, hi = hi_, open_left = FALSE, open_right = FALSE)
          }
        }
        if (length(segments) > 1L)
          segments <- segments[order(vapply(segments, function(s) s$lo, numeric(1)), method = "radix")]
      }
      out$shape <- if (length(segments)) "grid" else "empty"
      out$segments <- segments
    }
  }

  shape <- out$shape
  if (identical(shape, "unavailable")) { out$contains_zero <- NA; return(out) }
  out$contains_zero <- switch(shape,
    interval = (out$lower <= 0 && 0 <= out$upper),
    rays = (0 <= out$left_upper || 0 >= out$right_lower),
    rays_left_only = (0 <= out$left_upper),
    rays_right_only = (0 >= out$right_lower),
    line = TRUE,
    grid = any(vapply(out$segments, function(s) s$lo <= 0 && 0 <= s$hi, logical(1))),
    FALSE)
  if (is.finite(stat0)) {
    test_says_in <- stat0 <= cv
    if (test_says_in != out$contains_zero) out$contains_zero <- test_says_in
  }
  out
}

# "Unbounded or disconnected": rays, the whole line, or a grid set with more
# than one segment or a segment reaching the search boundary.
.ar_unbounded_or_disconnected <- function(ar) {
  if (is.null(ar)) return(FALSE)
  shape <- ar$shape
  if (shape %in% c("rays", "rays_left_only", "rays_right_only", "line")) return(TRUE)
  if (identical(shape, "grid")) {
    segs <- ar$segments
    return(length(segs) > 1L || any(vapply(segs, function(s) isTRUE(s$open_left) || isTRUE(s$open_right), logical(1))))
  }
  FALSE
}

.format_ar_set <- function(ar, fmt = "%.4f") {
  f <- function(x) sprintf(fmt, x)
  shape <- ar$shape
  if (identical(shape, "interval")) return(sprintf("[%s, %s]", f(ar$lower), f(ar$upper)))
  if (identical(shape, "rays"))
    return(sprintf("(-inf, %s] U [%s, +inf)   <- unbounded: weak-IV signature", f(ar$left_upper), f(ar$right_lower)))
  if (identical(shape, "rays_left_only")) return(sprintf("(-inf, %s]   <- unbounded: weak-IV signature", f(ar$left_upper)))
  if (identical(shape, "rays_right_only")) return(sprintf("[%s, +inf)   <- unbounded: weak-IV signature", f(ar$right_lower)))
  if (identical(shape, "line")) return("(-inf, +inf)   <- AR cannot reject any value at this level (weak-IV signature)")
  if (identical(shape, "grid")) {
    parts <- vapply(ar$segments, function(s) {
      lo_s <- if (isTRUE(s$open_left)) paste0("< ", f(s$lo)) else f(s$lo)
      hi_s <- if (isTRUE(s$open_right)) paste0("> ", f(s$hi)) else f(s$hi)
      sprintf("[%s, %s]", lo_s, hi_s)
    }, character(1))
    any_open <- any(vapply(ar$segments, function(s) isTRUE(s$open_left) || isTRUE(s$open_right), logical(1)))
    note <- if (any_open) "   (grid search; '<' / '>' mark segments that reach the search grid's boundary)" else ""
    return(paste0(paste(parts, collapse = " U "), note))
  }
  if (identical(shape, "unavailable")) return("unavailable: negative variance estimate")
  "empty set (AR rejects every value: model misspecification signal)"
}


# =============================================================================
#  Post-double-selection LASSO (BCH plug-in; sklearn-replicating coordinate
#  descent so the SELECTED SET matches Python; WP2)
# =============================================================================

.VALID_EXPAND_OPS <- c("square", "cube", "interact_with_X", "interact_with_dummies")

# Continuous iff more than 2 unique values (matches pd.unique(x).size > 2).
.is_continuous <- function(values) length(unique(values)) > 2L

# Build the rich expansion W. Column ORDER mirrors the Python builder:
# originals (X, then dummies) -> squares/cubes (per var: sq then cb) -> x*x ->
# d*x -> d*d. Returns list(W, col_names, meta).
.build_rich_expansion <- function(df, x_names, dummy_col_names,
                                  expand_x = c("square", "cube"),
                                  expand_dummies = c("interact_with_X"),
                                  expand_override = NULL) {
  if (is.null(expand_override)) expand_override <- list()
  for (ops in c(list(expand_x, expand_dummies), expand_override))
    for (op in ops) if (!(op %in% .VALID_EXPAND_OPS))
      stop(sprintf("Unknown expansion op '%s'. Valid: %s", op, paste(sort(.VALID_EXPAND_OPS), collapse = ", ")))
  n <- nrow(df)
  cont_x <- x_names[vapply(x_names, function(v) .is_continuous(df[[v]]), logical(1))]
  name_pow <- function(v, p) sprintf("%s^%d", v, p)
  name_int <- function(a, b) sprintf("%s \u00d7 %s", a, b)
  ops_for <- function(var, default) if (var %in% names(expand_override)) expand_override[[var]] else default

  cols <- list(); names_v <- character(0)
  meta <- list(n_originals_x = length(x_names), n_originals_dummies = length(dummy_col_names),
               n_continuous_x = length(cont_x), n_dummies = length(dummy_col_names),
               n_squares = 0L, n_cubes = 0L, n_interactions_xx = 0L,
               n_interactions_dx = 0L, n_interactions_dd = 0L, warnings = character(0))
  add <- function(vec, nm) { cols[[length(cols) + 1L]] <<- as.numeric(vec); names_v <<- c(names_v, nm) }

  for (v in x_names)        add(df[[v]], v)
  for (d in dummy_col_names) add(df[[d]], d)
  for (v in x_names) {
    ops <- ops_for(v, expand_x); is_cont <- v %in% cont_x
    if ("square" %in% ops) { if (!is_cont) meta$warnings <- c(meta$warnings, sprintf("'square' has no effect on indicator-like variable '%s'; skipped.", v))
                             else { add(df[[v]]^2, name_pow(v, 2)); meta$n_squares <- meta$n_squares + 1L } }
    if ("cube"   %in% ops) { if (!is_cont) meta$warnings <- c(meta$warnings, sprintf("'cube' has no effect on indicator-like variable '%s'; skipped.", v))
                             else { add(df[[v]]^3, name_pow(v, 3)); meta$n_cubes <- meta$n_cubes + 1L } }
  }
  for (vi in seq_along(x_names)) {
    v <- x_names[vi]; ops <- ops_for(v, expand_x)
    if (!("interact_with_X" %in% ops)) next
    for (ui in seq_along(x_names)) {
      if (ui <= vi) next
      u <- x_names[ui]
      add(df[[v]] * df[[u]], name_int(v, u)); meta$n_interactions_xx <- meta$n_interactions_xx + 1L
    }
  }
  for (d in dummy_col_names) {
    ops <- ops_for(d, expand_dummies)
    if ("square" %in% ops || "cube" %in% ops)
      meta$warnings <- c(meta$warnings, sprintf("'square'/'cube' on dummy variable '%s' has no effect; skipped.", d))
    if ("interact_with_X" %in% ops)
      for (u in x_names) { add(df[[d]] * df[[u]], name_int(d, u)); meta$n_interactions_dx <- meta$n_interactions_dx + 1L }
  }
  for (di in seq_along(dummy_col_names)) {
    d <- dummy_col_names[di]; ops <- ops_for(d, expand_dummies)
    if (!("interact_with_dummies" %in% ops)) next
    for (ei in seq_along(dummy_col_names)) {
      if (ei <= di) next
      e <- dummy_col_names[ei]
      if (startsWith(d, "_D_") && startsWith(e, "_D_")) {
        src_d <- sub("_[^_]*$", "", d); src_e <- sub("_[^_]*$", "", e)
        if (src_d == src_e) next   # same source -> product is 0
      }
      add(df[[d]] * df[[e]], name_int(d, e)); meta$n_interactions_dd <- meta$n_interactions_dd + 1L
    }
  }

  if (length(cols) == 0L) W <- matrix(0.0, nrow = n, ncol = 0L)
  else { W <- do.call(cbind, cols); storage.mode(W) <- "double" }
  if (length(names_v) != length(unique(names_v))) {
    dups <- sort(unique(names_v[duplicated(names_v)]))
    stop(sprintf("Rich expansion produced duplicate column names: %s. Remove the duplicate from X, exclude the op via expand_override, or rename the column.",
                 paste(dups, collapse = ", ")))
  }
  # Drop near-constant columns (POPULATION variance, as numpy's W.var(axis=0)).
  if (ncol(W) > 0L) {
    mu <- colMeans(W)
    col_var <- colSums(sweep(W, 2, mu, "-")^2) / nrow(W)
    keep_mask <- col_var > (1e-12 * pmax(abs(mu), 1.0)^2)
    n_dropped <- sum(!keep_mask)
    if (n_dropped > 0L) meta$warnings <- c(meta$warnings, sprintf("Dropped %d near-constant columns from rich expansion.", n_dropped))
    W <- W[, keep_mask, drop = FALSE]; names_v <- names_v[keep_mask]
  }
  for (msg in utils::head(meta$warnings, 3)) .inform(sprintf("  Note: %s\n", msg))
  if (length(meta$warnings) > 3) .inform(sprintf("  Note: ... (+%d more expansion warnings)\n", length(meta$warnings) - 3))
  meta$p <- ncol(W)
  list(W = W, col_names = names_v, meta = meta)
}

# Resolve lasso_keep (NULL / string / vector) against the rich-expansion names; globs allowed.
.resolve_lasso_keep <- function(lasso_keep, rich_col_names) {
  if (is.null(lasso_keep)) return(character(0))
  if (is.character(lasso_keep) && length(lasso_keep) == 1L && !nzchar(trimws(lasso_keep))) return(character(0))
  tokens <- if (is.character(lasso_keep) && length(lasso_keep) == 1L) strsplit(trimws(lasso_keep), "\\s+")[[1]] else as.character(lasso_keep)
  avail <- paste0(paste(utils::head(rich_col_names, 8), collapse = ", "), if (length(rich_col_names) > 8) "..." else "")
  resolved <- character(0)
  for (tok in tokens) {
    if (.has_glob(tok)) {
      matches <- rich_col_names[grepl(.glob_to_regex(tok), rich_col_names)]
      if (length(matches) == 0L)
        stop(sprintf("lasso_keep: glob pattern '%s' matched 0 names in the rich expansion. Available names include: %s", tok, avail))
      resolved <- c(resolved, setdiff(matches, resolved))
    } else {
      if (!(tok %in% rich_col_names))
        stop(sprintf("lasso_keep variable '%s' is not in the rich expansion. lasso_keep accepts originals (e.g., 'X1'), powers (e.g., 'X1^2'), and interactions (e.g., 'X1 \u00d7 X2'). Available names include: %s", tok, avail))
      if (!(tok %in% resolved)) resolved <- c(resolved, tok)
    }
  }
  resolved
}

# Per-feature BCH loadings psi_j from the current residuals u (W and u centered,
# weights normalized to mean 1): the score of the sqrt(w)-scaled data in every
# form (A.1-3): heteroskedastic sum(w^2 W^2 u^2)/sum(w), cluster
# sum_g (sum_{i in g} w u W)^2 / sum(w), spatial (w u W)' A (w u W) / sum(w).
.compute_psi <- function(W, u, weights, cluster_codes = NULL, spatial_A = NULL) {
  n_w <- sum(weights)
  wu <- weights * u
  if (!is.null(spatial_A)) {
    wuW <- W * wu
    psi2 <- pmax(colSums(wuW * as.matrix(spatial_A %*% wuW)), 0)
    psi <- sqrt(psi2 / n_w)
  } else if (is.null(cluster_codes)) {
    psi <- sqrt(colSums(W * W * (wu * wu)) / n_w)
  } else {
    cs <- rowsum(W * wu, cluster_codes, reorder = TRUE)
    psi <- sqrt(colSums(cs^2) / n_w)
  }
  pmax(psi, 1e-10 * max(abs(W), 1.0))
}

# Effective sample size for the penalty: n with equal weights, else the Kish n
# (also under clustering, D5).
.compute_n_eff <- function(weights) {
  if (all(abs(weights - weights[1]) <= 1e-8 + 1e-5 * abs(weights[1]))) return(as.numeric(length(weights)))
  sum(weights)^2 / sum(weights^2)
}

# BCH plug-in penalty level: lambda = 2 c sqrt(n_eff) Phi^{-1}(1 - alpha_p/(2p)).
# The default alpha_p = 0.1/log(n) uses the number of observations in every SE
# regime (A.1-8); lassopack's cluster default is 0.1/log(G) instead.
.bch_lambda <- function(p, n_eff, c_slack, alpha_p) {
  if (p == 0L) return(0.0)
  2.0 * c_slack * sqrt(n_eff) * stats::qnorm(1.0 - alpha_p / (2.0 * p))
}

# Weighted least-squares residuals of y on the columns of W (raw scale).
.wls_resid <- function(y, W, sqrt_w) {
  if (is.null(W) || ncol(W) == 0L) return(y)
  q <- qr(W * sqrt_w, tol = 1e-12)
  as.numeric(qr.resid(q, y * sqrt_w)) / sqrt_w
}

# Residuals of the post-LASSO OLS of y on the selected columns (WP2.4, D7).
.post_lasso_resid <- function(y_c, W_c, sel, sqrt_w) {
  if (!any(sel)) return(y_c)
  .wls_resid(y_c, W_c[, sel, drop = FALSE], sqrt_w)
}

# Cyclic coordinate descent replicating sklearn's enet_coordinate_descent for
# Lasso (l1_ratio=1, fit_intercept=FALSE): minimize (1/(2n))||y-Xb||^2 + alpha*||b||_1.
.lasso_cd <- function(X, y, alpha, max_iter = 20000L, tol = 1e-6) {
  n <- nrow(X); p <- ncol(X)
  beta <- numeric(p)
  norm_cols <- colSums(X * X)
  R <- y - as.numeric(X %*% beta)
  alpha_n <- alpha * n
  tol_gap <- tol * sum(y * y)
  for (it in seq_len(max_iter)) {
    w_max <- 0.0; d_w_max <- 0.0
    for (j in seq_len(p)) {
      nj <- norm_cols[j]
      if (nj == 0) next
      w_old <- beta[j]
      tmp <- sum(X[, j] * R) + nj * w_old
      bj <- sign(tmp) * max(abs(tmp) - alpha_n, 0) / nj
      if (bj != w_old) R <- R - X[, j] * (bj - w_old)
      beta[j] <- bj
      dw <- abs(bj - w_old); if (dw > d_w_max) d_w_max <- dw
      if (abs(bj) > w_max) w_max <- abs(bj)
    }
    if (w_max == 0.0 || d_w_max / w_max < tol || it == max_iter) {
      XtR <- as.numeric(crossprod(X, R))
      dual_norm_XtA <- max(abs(XtR))
      R_norm2 <- sum(R * R)
      if (dual_norm_XtA > alpha_n) { const <- alpha_n / dual_norm_XtA; gap <- 0.5 * (R_norm2 + R_norm2 * const^2) }
      else                        { const <- 1.0;                       gap <- R_norm2 }
      gap <- gap + alpha_n * sum(abs(beta)) - const * sum(R * y)
      if (gap < tol_gap) break
    }
  }
  beta
}

# Iterative plug-in LASSO of y on W with an unpenalized intercept (weighted
# centering, WP2.1), weights normalized to mean 1 (WP2.2), keep columns
# partialed out, and loadings updated from post-LASSO residuals (WP2.4).
.bch_lasso <- function(y_raw, W_raw, weights, cluster_codes = NULL, n_iter_max = 15L, tol = 1e-6,
                       c_slack = 1.1, alpha_p = NULL, lambda_method = "plugin",
                       keep_idx = NULL, lambda_override = NULL, spatial_A = NULL, verbose = FALSE) {
  n <- nrow(W_raw); p <- ncol(W_raw)
  if (is.null(alpha_p)) alpha_p <- 0.1 / max(log(max(n, 2)), 1.0)
  if (p == 0L) return(list(beta = numeric(0), selected = logical(0), lambda_value = 0.0,
                           n_iter = 0L, converged = TRUE, method = "trivial"))
  if (identical(lambda_method, "cv") && is.null(lambda_override))
    stop("lasso_penalty='cv' is not available in the R port: sklearn's KFold/GroupKFold ",
         "fold assignment cannot be reproduced in R, so the CV path cannot be validated ",
         "against the Python output. Use lasso_penalty='plugin' (the default, and the ",
         "paper's choice), or pass lasso_lambda_override to fix lambda directly.", call. = FALSE)
  weights <- as.numeric(weights); weights <- weights / mean(weights)
  y_raw <- y_raw - sum(weights * y_raw) / sum(weights)
  W_raw <- sweep(W_raw, 2, colSums(W_raw * weights) / sum(weights), "-")
  sqrt_w <- sqrt(weights)
  if (is.null(keep_idx)) keep_idx <- integer(0)
  keep_idx <- as.integer(keep_idx)

  if (length(keep_idx) > 0L) {
    nonkeep_idx <- setdiff(seq_len(p), keep_idx)
    qk <- qr(W_raw[, keep_idx, drop = FALSE] * sqrt_w, tol = 1e-12)
    y_lasso_raw <- as.numeric(qr.resid(qk, y_raw * sqrt_w)) / sqrt_w
    W_lasso_raw <- if (length(nonkeep_idx)) qr.resid(qk, W_raw[, nonkeep_idx, drop = FALSE] * sqrt_w) / sqrt_w
                   else matrix(0, n, 0)
    W_lasso_raw <- matrix(W_lasso_raw, nrow = n)
    p_lasso <- ncol(W_lasso_raw)
  } else {
    nonkeep_idx <- seq_len(p)
    y_lasso_raw <- y_raw; W_lasso_raw <- W_raw
    p_lasso <- p
  }
  y_lasso_pm <- y_lasso_raw * sqrt_w
  W_lasso_pm <- W_lasso_raw * sqrt_w
  n_eff <- .compute_n_eff(weights)

  if (p_lasso == 0L) {
    selected <- logical(p); if (length(keep_idx) > 0L) selected[keep_idx] <- TRUE
    return(list(beta = numeric(p), selected = selected, lambda_value = 0.0,
                n_iter = 0L, converged = TRUE, method = "keep-only"))
  }

  if (!is.null(lambda_override)) {
    psi_lasso <- .compute_psi(W_lasso_raw, y_lasso_raw, weights, cluster_codes, spatial_A)
    col_scale <- 1.0 / psi_lasso
    lam <- lambda_override
    beta_cd <- .lasso_cd(sweep(W_lasso_pm, 2, col_scale, "*"), y_lasso_pm, lam / (2.0 * n))
    beta_unscaled <- beta_cd * col_scale
    beta <- numeric(p); beta[nonkeep_idx] <- beta_unscaled
    selected <- logical(p); selected[nonkeep_idx] <- (beta_unscaled != 0)
    if (length(keep_idx) > 0L) selected[keep_idx] <- TRUE
    return(list(beta = beta, selected = selected, lambda_value = lam, n_iter = 1L,
                converged = TRUE, method = "override"))
  }

  u <- y_lasso_raw
  beta_lasso <- numeric(p_lasso); n_iter <- 0L; converged <- FALSE; lam <- 0.0
  for (it in seq_len(n_iter_max)) {
    n_iter <- it
    psi_lasso <- .compute_psi(W_lasso_raw, u, weights, cluster_codes, spatial_A)
    lam <- .bch_lambda(p_lasso, n_eff, c_slack, alpha_p)
    col_scale <- 1.0 / psi_lasso
    beta_cd <- .lasso_cd(sweep(W_lasso_pm, 2, col_scale, "*"), y_lasso_pm, lam / (2.0 * n))
    beta_new <- beta_cd * col_scale
    u <- .post_lasso_resid(y_lasso_raw, W_lasso_raw, beta_new != 0, sqrt_w)
    max_change <- max(abs(beta_new - beta_lasso))
    if (verbose) cat(sprintf("    iter %d: lambda=%.3f, |dbeta|_inf=%.2e, active=%d/%d\n", it, lam, max_change, sum(beta_new != 0), p_lasso))
    beta_lasso <- beta_new
    if (max_change < tol) { converged <- TRUE; break }
  }
  beta <- numeric(p); beta[nonkeep_idx] <- beta_lasso
  selected <- logical(p); selected[nonkeep_idx] <- (beta_lasso != 0)
  if (length(keep_idx) > 0L) selected[keep_idx] <- TRUE
  list(beta = beta, selected = selected, lambda_value = lam, n_iter = n_iter,
       converged = converged, method = "plugin")
}

# Run four LASSOs (Y, T, M, Z[k]) on the shared rich expansion; union selection.
.quadruple_lasso <- function(y_raw, t_raw, m_raw, z_mat_raw, W_raw, col_names, weights, cluster_codes,
                             keep_idx, lambda_method = "plugin", n_iter_max = 15L, tol = 1e-6,
                             c_slack = 1.1, alpha_p = NULL, lambda_override = NULL,
                             spatial_A = NULL, verbose = FALSE) {
  if (is.null(dim(z_mat_raw))) z_mat_raw <- matrix(z_mat_raw, ncol = 1L)
  K_inst <- ncol(z_mat_raw); p <- ncol(W_raw)
  if (is.null(keep_idx)) keep_idx <- integer(0)
  call_bch <- function(y) .bch_lasso(y, W_raw, weights, cluster_codes, n_iter_max = n_iter_max, tol = tol,
                                     c_slack = c_slack, alpha_p = alpha_p, lambda_method = lambda_method,
                                     keep_idx = keep_idx, lambda_override = lambda_override,
                                     spatial_A = spatial_A, verbose = verbose)
  res_Y <- call_bch(y_raw); res_T <- call_bch(t_raw); res_M <- call_bch(m_raw)
  res_Z_per <- lapply(seq_len(K_inst), function(k) call_bch(z_mat_raw[, k]))
  S_Y <- res_Y$selected; S_T <- res_T$selected; S_M <- res_M$selected
  S_Z <- logical(p); for (r in res_Z_per) S_Z <- S_Z | r$selected
  S <- S_Y | S_T | S_M | S_Z
  if (length(keep_idx) > 0L) S[as.integer(keep_idx)] <- TRUE
  summarize <- function(res, label) {
    idx <- which(res$selected)
    list(label = label, names = col_names[idx], coefs = res$beta[res$selected], indices = idx,
         lambda_value = as.numeric(res$lambda_value), n_selected = sum(res$selected),
         converged = res$converged, n_iter = res$n_iter, method = res$method)
  }
  list(S = S, S_Y = S_Y, S_T = S_T, S_M = S_M, S_Z = S_Z,
       S_Z_per_inst = lapply(res_Z_per, function(r) r$selected),
       details = list(Y = summarize(res_Y, "Y"), T = summarize(res_T, "T"), M = summarize(res_M, "M"),
                      Z = lapply(seq_along(res_Z_per), function(k) summarize(res_Z_per[[k]], sprintf("Z[%d]", k - 1L)))),
       union_names = col_names[which(S)], n_selected_union = sum(S), n_total = p)
}


# =============================================================================
#  Main analysis function
# =============================================================================

#' IV mediation analysis with a single instrument
#'
#' Estimates how much of the effect of a treatment `T` on an outcome `Y` runs
#' through a mediator `M` (the indirect effect, theta times gamma) and how much
#' runs around it (the direct effect, tau), with one instrument `Z` for the
#' treatment (or several), as in Ferrara, Gold, Heblich and Pinto (2026). The
#' function also reports the first-stage diagnostics (Olea-Pflueger effective F
#' statistics and critical values), Anderson-Rubin confidence sets, the kappa-
#' and R^2-based sensitivity analyses and, optionally, a post-double-selection
#' LASSO for the controls. Print the result with [print.ivmediate()].
#'
#' @param data A data frame (data.table and tibble inputs are accepted). It is
#'   never modified.
#' @param Y Name of the outcome variable.
#' @param T Name of the treatment variable.
#' @param M Name of the mediator variable.
#' @param Z Name(s) of the instrument variable(s), typically one. They must
#'   differ from `Y`, `T` and `M` and may not repeat.
#' @param X Control variables: a character vector of column names or one
#'   whitespace-separated string (multi-line strings are fine). Elements and
#'   tokens may be shell-style globs such as `"reg_*"`, resolved against
#'   `names(data)`. Globs never pick up `Y`, `T`, `M`, `Z` or the fixed-effect,
#'   cluster, weight and coordinate columns (a note lists the excluded
#'   matches); naming `Y`, `T`, `M` or `Z` explicitly is an error. Exactly
#'   collinear controls (after absorbing the fixed effects) are dropped with a
#'   note, keeping the earlier column of a collinear set; the degrees of
#'   freedom count only the kept columns. Default: no controls.
#' @param dummies Categorical variables to expand into indicator dummies
#'   (drop-first). Same parsing rules and role protection as `X`. Empty or
#'   whitespace-only text values count as missing. Default: none.
#' @param fe Name(s) of one or two fixed-effect variables, absorbed by
#'   within-group demeaning (with weighted group means when weights are given).
#'   Two-way FE use an exact projection when the FE dimension with fewer levels
#'   has at most 3,000 levels, otherwise alternating projections with
#'   Irons-Tuck acceleration (relative tolerance 1e-12). Singleton groups are
#'   dropped iteratively. Absorbed degrees of freedom: the rank of the FE
#'   dummies beyond the intercept. Empty or whitespace-only text values count
#'   as missing. Default: none.
#' @param se Standard errors: `"classical"`, `"robust"` (HC1), `"cluster"` (CR1)
#'   or `"spatial"` (Conley). Default: `"robust"`.
#' @param clusterid Name of the cluster variable. Required with
#'   `se = "cluster"` and used only then (otherwise a note says it is
#'   ignored). Fixed effects nested in the clusters are not counted in the CR1
#'   small-sample factor.
#' @param alpha Significance level for the confidence intervals, the AR sets and
#'   the RV/XRV threshold, in (0, 1). Default: 0.05.
#' @param weights Name of the weight column. Negative weights are an error;
#'   zero-weight rows are dropped (with a note) before singleton detection and
#'   FE demeaning. Default: unweighted.
#' @param wtype Weight type: `"aweight"` (analytic) or `"pweight"` (probability;
#'   with `se = "classical"` the SEs are upgraded to `"robust"`). Default:
#'   `"aweight"`.
#' @param op_tau Worst-case 2SLS bias tolerance, as a fraction of the OLS bias
#'   benchmark, for the Olea-Pflueger (2013) critical values of the
#'   first-stage F statistics, in (0, 1). Default: 0.10.
#' @param op_alpha Size of the Olea-Pflueger weak-instrument test, in (0, 1).
#'   The critical values come from the simplified test at the estimated
#'   effective degrees of freedom K_eff (K_eff = 1 for one instrument),
#'   computed for any `op_tau` and `op_alpha`. Default: 0.05.
#' @param ar_ci If `TRUE` (default), invert the Anderson-Rubin test to obtain
#'   weak-instrument-robust confidence sets for theta (the mediator effect),
#'   gamma (treatment on mediator) and tau_total, with the SE regime and
#'   `alpha` of the analysis. With one instrument the set is computed in closed
#'   form (Fieller inversion: a bounded interval, two rays or the whole line);
#'   with several instruments the joint chi^2(K) statistic is computed in
#'   closed form, located on a search grid and its boundaries refined by
#'   Brent's method. With spatial SEs, a set whose variance estimate turns
#'   negative is reported as unavailable. No AR set is reported for the direct
#'   effect tau.
#' @param ar_grid_halfwidth Half-width of the AR search grid, in units of the
#'   parameter's analytic SE. Default: 200.
#' @param ar_grid_points Number of AR search-grid points (an integer of at
#'   least 3). Default: 4001.
#' @param lat,lon Names of the latitude and longitude columns (in degrees, or
#'   projected coordinates with `spatial_metric = "euclidean"`). Required with
#'   `se = "spatial"` and an error otherwise. Rows with missing coordinates are
#'   dropped.
#' @param spatial_cutoff Distance cutoff of the spatial kernel, required with
#'   `se = "spatial"`: kilometers with the haversine metric, otherwise the
#'   units of `lat` and `lon`. Conley SEs are sensitive to the cutoff; report
#'   several. A uniform-kernel cutoff that covers every pair is an error, and a
#'   neighbor share above one half prints a warning.
#' @param spatial_kernel `"uniform"` (weight 1 within the cutoff) or
#'   `"bartlett"` (linear taper); `NULL` means `"uniform"` with
#'   `se = "spatial"`. Giving it with another `se` is an error.
#' @param spatial_metric `"haversine"` (great-circle distance in km for
#'   coordinates in degrees, Earth radius 6371 km) or `"euclidean"`; `NULL`
#'   means `"haversine"` with `se = "spatial"`. Giving it with another `se` is
#'   an error.
#' @param lasso If `TRUE`, run the post-double-selection ("quadruple") LASSO
#'   before the structural 2SLS: one plug-in LASSO each for `Y`, `T`, `M` and
#'   every instrument on a rich expansion of `X` and `dummies`; the union of
#'   the selections forms the control set. Each LASSO has an unpenalized
#'   intercept and uses weights normalized to mean 1. Default: `FALSE`.
#' @param lasso_penalty `"plugin"`: the iterative BCH plug-in penalty with
#'   post-LASSO residual loadings (deterministic). `"cv"` (cross-validation) is
#'   not available in R and stops with a message. Default: `"plugin"`.
#' @param lasso_keep Names that always enter the control set, resolved against
#'   the rich expansion: originals (`"X1"`), powers (`"X1^2"`) and interactions.
#'   Globs are supported. These columns are partialed out before each LASSO.
#'   Default: none.
#' @param lasso_seed,lasso_cv_folds Used only by the cross-validated penalty,
#'   which is not available in R; accepted for compatibility with the Python
#'   implementation.
#' @param expand_x Expansion operations for the variables of `X` in the LASSO's
#'   candidate set: `"square"` and `"cube"` add powers of the continuous
#'   variables (indicator-like variables are skipped with a warning), and
#'   `"interact_with_X"` adds the products of pairs of `X` variables. Default:
#'   `c("square", "cube")`.
#' @param expand_dummies Expansion operations for the indicators made from
#'   `dummies`: `"interact_with_X"` adds their products with every `X`
#'   variable, and `"interact_with_dummies"` the products of indicators of
#'   different variables. Powers have no effect on indicators. Default:
#'   `"interact_with_X"`. In both arguments the valid names are `"square"`,
#'   `"cube"`, `"interact_with_X"` and `"interact_with_dummies"`; any other name
#'   is an error.
#' @param expand_override Per-variable override of `expand_x` and
#'   `expand_dummies`: a named list of operation vectors, for example
#'   `list(age = "square")`. Default: `NULL`.
#' @param show_selection If `TRUE`, the print method shows the LASSO selections
#'   (sets, penalties, convergence). The selections are in the results either
#'   way. Default: `FALSE`.
#' @param lasso_c_slack Slack constant c of the plug-in penalty
#'   lambda = 2 c sqrt(n) qnorm(1 - alpha_p / (2p)). Default: 1.1.
#' @param lasso_alpha_p Tail probability alpha_p of the plug-in penalty; `NULL`
#'   means 0.1 / log(n), with n the number of observations also under
#'   clustering (lassopack's cluster default is 0.1 / log(G)).
#' @param lasso_n_iter_max Maximum number of plug-in loading iterations.
#'   Default: 15.
#' @param lasso_tol Convergence tolerance on the largest coefficient change
#'   between plug-in iterations. Default: 1e-6.
#' @param lasso_lambda_override A fixed penalty level that bypasses the
#'   plug-in rule, for example to replicate a published lambda. Default:
#'   `NULL`.
#' @param lasso_verbose Print per-iteration progress for each LASSO. Default:
#'   `FALSE`.
#' @param .point_only Internal: compute point estimates only (used by
#'   [ivmediate_bootstrap()] for its replicates). Leave at `FALSE`.
#'
#' @details
#' Logical options must be `TRUE` or `FALSE`; text such as `"false"` is an
#' error. Wald confidence intervals use normal critical values; AR sets use
#' the chi^2 distribution with K degrees of freedom. The R implementation
#' agrees with the Python implementation to rounding error on the
#' cross-language validation harness of the repository.
#'
#' @return A list of class `c("ivmediate", "list")` with the five estimates,
#'   their SEs and Wald CIs (`theta_hat`, `tau_hat`, `gamma_hat`, `total_hat`,
#'   `indirect` with `se_*` and `ci_*`); the first stage of `T` on `Z` (`pi1`,
#'   `se_pi1`); sample information (`n`, `n_clusters`, `df_fe_absorbed`,
#'   `fe_nested_in_cluster`, `dropped_collinear`); the F statistics with
#'   `K_eff_TZ`, `K_eff_MZ_T` and the Olea-Pflueger critical values; the AR sets
#'   (`ar_theta`, `ar_gamma`, `ar_total`); the kappa block (`c_TM_hat`,
#'   `kappa_star`, `kappa_star_CI`, `kappa_star_CI_informative`, `rho_TY_star`,
#'   `kappa_adj_table`, `F_MZ_T_above_crit`); the R^2 block (`F_reduced_form`,
#'   `F_first_stage_cond`, `RV_theta`, `XRV_theta` and their components); the
#'   LASSO output; the residuals `eps_T_hat`, `eps_M_hat` and `eps_Y_hat`; and
#'   the package `version`.
#'
#' @references
#' Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo Pinto (2026),
#' "Identifying Causal Mediation Effects with a Single Instrument," NBER
#' Working Paper No. \[number\], doi:\[DOI\].
#'
#' Montiel Olea, Jose Luis, and Carolin Pflueger (2013), "A Robust Test for
#' Weak Instruments," Journal of Business and Economic Statistics 31(3),
#' 358-369.
#'
#' @seealso [print.ivmediate()] for the printout, [ivmediate_bootstrap()] for
#'   bootstrap inference.
#'
#' @examples
#' # Becker and Woessmann (2009): Protestantism, literacy and income
#' # (Table 6, column 1 of the paper)
#' bw <- read.csv(system.file("extdata", "becker_woessmann_2009.csv",
#'                            package = "ivmediate"))
#' res <- ivmediate(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw",
#'                  Z = "kmwittenberg",
#'                  X = "f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf
#'                       f_dumb hhsize lnpop gpop")
#' print(res, title = "Becker and Woessmann (2009)")
#' res$theta_hat
#'
#' # The optional kappa-adjusted table and R^2-based section
#' print(res, kappa_table = TRUE, r2_sensitivity = TRUE)
#' @export
ivmediate <- function(data, Y, T, M, Z, X = NULL, dummies = NULL, fe = NULL,
                      se = "robust", clusterid = NULL, alpha = 0.05,
                      weights = NULL, wtype = "aweight",
                      op_tau = 0.10, op_alpha = 0.05,
                      ar_ci = TRUE, ar_grid_halfwidth = 200.0, ar_grid_points = 4001,
                      lat = NULL, lon = NULL, spatial_cutoff = NULL,
                      spatial_kernel = NULL, spatial_metric = NULL,
                      lasso = FALSE, lasso_penalty = "plugin",
                      lasso_keep = NULL, lasso_seed = 20260317, lasso_cv_folds = 5,
                      expand_x = c("square", "cube"),
                      expand_dummies = c("interact_with_X"),
                      expand_override = NULL, show_selection = FALSE,
                      lasso_c_slack = 1.1, lasso_alpha_p = NULL,
                      lasso_n_iter_max = 15, lasso_tol = 1e-6,
                      lasso_lambda_override = NULL, lasso_verbose = FALSE,
                      .point_only = FALSE) {
  # ---- inputs and options (WP5.5) ----
  if (!is.data.frame(data)) stop("data must be a data.frame (data.table and tibble inputs are accepted).", call. = FALSE)
  data <- as.data.frame(data)
  .check_choice("se", se, .SE_TYPES)
  .check_choice("wtype", wtype, c("aweight", "pweight"))
  .check_open_unit("alpha", alpha)
  .check_open_unit("op_tau", op_tau)
  .check_open_unit("op_alpha", op_alpha)
  .check_bool("ar_ci", ar_ci); .check_bool("lasso", lasso); .check_bool("show_selection", show_selection)
  .check_bool("lasso_verbose", lasso_verbose); .check_bool(".point_only", .point_only)
  .check_choice("lasso_penalty", lasso_penalty, c("plugin", "cv"))
  if (!(is.numeric(ar_grid_halfwidth) && length(ar_grid_halfwidth) == 1L && is.finite(ar_grid_halfwidth) && ar_grid_halfwidth > 0))
    stop(sprintf("ar_grid_halfwidth must be positive; got %s.", .show_value(ar_grid_halfwidth)), call. = FALSE)
  if (!(.is_whole(ar_grid_points) && ar_grid_points >= 3))
    stop(sprintf("ar_grid_points must be an integer >= 3; got %s.", .show_value(ar_grid_points)), call. = FALSE)
  if (lasso) {
    if (!(is.numeric(lasso_c_slack) && length(lasso_c_slack) == 1L && is.finite(lasso_c_slack) && lasso_c_slack > 0))
      stop(sprintf("lasso_c_slack must be positive; got %s.", .show_value(lasso_c_slack)), call. = FALSE)
    if (!is.null(lasso_alpha_p)) .check_open_unit("lasso_alpha_p", lasso_alpha_p)
    if (!(.is_whole(lasso_n_iter_max) && lasso_n_iter_max >= 1))
      stop(sprintf("lasso_n_iter_max must be a positive integer; got %s.", .show_value(lasso_n_iter_max)), call. = FALSE)
    if (!(is.numeric(lasso_tol) && length(lasso_tol) == 1L && isTRUE(lasso_tol > 0)))
      stop(sprintf("lasso_tol must be positive; got %s.", .show_value(lasso_tol)), call. = FALSE)
    if (!(.is_whole(lasso_cv_folds) && lasso_cv_folds >= 2))
      stop(sprintf("lasso_cv_folds must be an integer >= 2; got %s.", .show_value(lasso_cv_folds)), call. = FALSE)
    if (!is.null(lasso_lambda_override) && !(is.numeric(lasso_lambda_override) && length(lasso_lambda_override) == 1L &&
                                             isTRUE(lasso_lambda_override >= 0)))
      stop(sprintf("lasso_lambda_override must be >= 0; got %s.", .show_value(lasso_lambda_override)), call. = FALSE)
  }

  # ---- spatial arguments only with se = "spatial" (D3) ----
  spatial_given <- c("lat", "lon", "spatial_cutoff", "spatial_kernel", "spatial_metric")[
    !vapply(list(lat, lon, spatial_cutoff, spatial_kernel, spatial_metric), is.null, logical(1))]
  if (se == "spatial") {
    if (is.null(lat) || is.null(lon)) stop("se='spatial' requires lat and lon column names.", call. = FALSE)
    if (is.null(spatial_cutoff))
      stop("se='spatial' requires spatial_cutoff (distance bandwidth). Conley (1999) results are sensitive to this choice; report results at multiple cutoffs as a robustness check.", call. = FALSE)
    if (!(is.numeric(spatial_cutoff) && length(spatial_cutoff) == 1L))
      stop(sprintf("spatial_cutoff must be a positive number; got %s.", .show_value(spatial_cutoff)), call. = FALSE)
    if (!(is.finite(spatial_cutoff) && spatial_cutoff > 0))
      stop(sprintf("spatial_cutoff must be positive; got %s.", .show_value(spatial_cutoff)), call. = FALSE)
    if (is.null(spatial_kernel)) spatial_kernel <- "uniform"
    if (is.null(spatial_metric)) spatial_metric <- "haversine"
    .check_choice("spatial_kernel", spatial_kernel, c("uniform", "bartlett"))
    .check_choice("spatial_metric", spatial_metric, c("haversine", "euclidean"))
  } else if (length(spatial_given) > 0L) {
    stop(sprintf("%s given but se='%s': the spatial arguments are used only with se='spatial'.",
                 paste(spatial_given, collapse = ", "), se), call. = FALSE)
  }

  # ---- cluster information only with cluster SEs (D3) ----
  if (se == "cluster" && is.null(clusterid))
    stop("se='cluster' requires clusterid (the cluster variable).", call. = FALSE)
  if (!is.null(clusterid) && se != "cluster") .inform("  Note: clusterid ignored (se != 'cluster').\n")

  # ---- weight-type / SE-regime consistency ----
  if (!is.null(weights) && wtype == "pweight" && se == "classical") {
    .inform("  Note: pweight + classical SE is inconsistent. Upgrading SE to 'robust' (HC1 sandwich).\n")
    se <- "robust"
  } else if (!is.null(weights) && wtype == "aweight" && se == "classical") {
    .inform("  Note: aweight + classical SE assumes the aweight specification is correctly calibrated. Consider 'robust' or 'cluster' as a defensive default.\n")
  }
  se_v <- if (.point_only) "none" else se          # variance regime actually computed
  if (.point_only) ar_ci <- FALSE

  # ---- estimation sample (WP5.1, WP1.6) ----
  prep <- .prepare_sample(data, Y, T, M, Z, X, dummies, fe, se, clusterid, weights, lat, lon)
  Z <- prep$Z; X <- prep$X; dummies <- prep$dummies
  fe_cols <- prep$fe_cols; clusterid <- prep$clusterid; df <- prep$df
  note_X <- .summarize_resolved("X", X); note_D <- .summarize_resolved("dummies", dummies)
  if (!is.null(note_X)) .inform(note_X, "\n", sep = "")
  if (!is.null(note_D)) .inform(note_D, "\n", sep = "")
  for (line in prep$notes) .inform(line, "\n", sep = "")
  is_weighted <- !is.null(weights)
  if (prep$n_zero_weights > 0L) .inform(sprintf("  Note: dropping %d observation(s) with zero weight.\n", prep$n_zero_weights))

  # ---- dummy indicators ----
  dummy_col_names <- character(0)
  if (length(dummies) > 0L) {
    ed <- .expand_dummies(df, dummies); df <- ed$df; dummy_col_names <- ed$new_cols
    if (!lasso) X <- c(X, dummy_col_names)
  }

  # ---- singletons (alternating over the FE dimensions until none remain) ----
  n_dropped_singletons <- 0L; n_fe_groups <- 0L
  if (!is.null(fe_cols)) {
    repeat {
      dropped_this_pass <- 0L
      for (fe_col in fe_cols) { ds <- .drop_singletons(df, fe_col); df <- ds$df; dropped_this_pass <- dropped_this_pass + ds$n_dropped }
      n_dropped_singletons <- n_dropped_singletons + dropped_this_pass
      if (dropped_this_pass == 0L) break
    }
    if (n_dropped_singletons > 0L)
      .inform(sprintf("  Note: dropped %d singleton observations across %d FE dimension(s).\n", n_dropped_singletons, length(fe_cols)))
  }
  n <- nrow(df)
  if (n == 0L) stop("No observations remaining after dropping missing values and singletons.", call. = FALSE)
  n_clusters <- if (!is.null(clusterid)) length(unique(df[[clusterid]])) else NULL
  if (!is.null(clusterid) && n_clusters < 2L) stop("cluster SEs need at least 2 clusters.", call. = FALSE)

  # ---- absorbed degrees of freedom (WP1.3) ----
  df_fe_absorbed <- 0L; df_vcov <- 0L; fe_nested_in_cluster <- NULL
  if (!is.null(fe_cols)) {
    fd <- .fe_degrees_of_freedom(df, fe_cols, if (se == "cluster") clusterid else NULL)
    df_fe_absorbed <- fd$d; n_fe_groups <- fd$n_groups
    if (se == "cluster") { fe_nested_in_cluster <- fd$nested; df_vcov <- fd$d_cr1 } else df_vcov <- df_fe_absorbed
  }

  # ---- LASSO: rich expansion built on levels, before FE demeaning ----
  rich_col_names <- character(0); rich_meta <- NULL; rich_in_df <- character(0); keep_idx_lasso <- integer(0)
  if (lasso) {
    rex <- .build_rich_expansion(df, x_names = X, dummy_col_names = dummy_col_names, expand_x = expand_x,
                                 expand_dummies = expand_dummies, expand_override = expand_override)
    rich_col_names <- rex$col_names; rich_meta <- rex$meta
    existing <- colnames(df)
    rich_in_df <- vapply(seq_along(rich_col_names), function(j) {
      nm <- rich_col_names[j]; if (nm %in% existing) sprintf("_RX_%d_%s", j - 1L, nm) else nm }, character(1))
    if (length(rich_in_df) > 0L) {
      Wdf <- as.data.frame(rex$W, stringsAsFactors = FALSE); colnames(Wdf) <- rich_in_df
      rn <- rownames(df); df <- cbind(df, Wdf); rownames(df) <- rn
    }
    keep_resolved <- .resolve_lasso_keep(lasso_keep, rich_col_names)
    keep_idx_lasso <- match(keep_resolved, rich_col_names)
    if (!is.null(rich_meta) && rich_meta$p > 0L)
      .inform(sprintf("  Rich expansion built: p = %d (originals_x=%d, originals_dummies=%d, sq=%d, cb=%d, x\u00d7x=%d, d\u00d7x=%d, d\u00d7d=%d)\n",
                  rich_meta$p, rich_meta$n_originals_x, rich_meta$n_originals_dummies, rich_meta$n_squares,
                  rich_meta$n_cubes, rich_meta$n_interactions_xx, rich_meta$n_interactions_dx, rich_meta$n_interactions_dd))
    if (length(keep_resolved) > 0L)
      .inform(sprintf("  LASSO keep (%d): %s%s\n", length(keep_resolved), paste(utils::head(keep_resolved, 8), collapse = ", "),
                  if (length(keep_resolved) > 8) "..." else ""))
    control_cols <- rich_in_df
  } else {
    control_cols <- X
  }

  w_vec <- if (is_weighted) as.numeric(df[[weights]]) else rep(1, n)
  sqrt_w <- sqrt(w_vec)
  # Reference norms (weighted, centered; before the FE are absorbed) for the
  # collinearity and identification checks.
  ref_norm <- function(v) { v <- as.numeric(v); sqrt(sum(w_vec * (v - sum(w_vec * v) / sum(w_vec))^2)) }
  ref_norm_ctrl <- vapply(control_cols, function(cc) ref_norm(df[[cc]]), numeric(1))
  ref_norm_role <- vapply(c(T, M, Z), function(cc) ref_norm(df[[cc]]), numeric(1))
  names(ref_norm_role) <- c(T, M, Z)

  # ---- absorb the FE (weighted group means, WP1.2) ----
  analysis_cols <- c(Y, T, M, Z, control_cols)
  if (!is.null(fe_cols)) {
    wts <- if (is_weighted) w_vec else NULL
    df <- if (length(fe_cols) == 1L) .demean_within_groups(df, analysis_cols, fe_cols[1], wts)
          else .demean_twoway(df, analysis_cols, fe_cols, wts)
  }

  y <- as.numeric(df[[Y]]); t <- as.numeric(df[[T]]); m <- as.numeric(df[[M]])
  z_mat <- as.matrix(df[, Z, drop = FALSE]); storage.mode(z_mat) <- "double"
  x_mat <- if (length(control_cols)) { xm <- as.matrix(df[, control_cols, drop = FALSE]); storage.mode(xm) <- "double"; xm } else NULL
  cl <- NULL
  if (!is.null(clusterid)) { cc <- .codes(df[[clusterid]]); cl <- list(codes = cc, G = max(cc)) }
  lat_arr <- lon_arr <- NULL
  if (se == "spatial") {
    lat_arr <- as.numeric(df[[lat]]); lon_arr <- as.numeric(df[[lon]])
    if (spatial_metric == "haversine" && (any(abs(lat_arr) > 90) || any(abs(lon_arr) > 180)))
      stop("Lat/lon values are outside the (-90, 90) x (-180, 180) degree range. Check that you have lat/lon in degrees, or switch to spatial_metric='euclidean' for projected coordinates.", call. = FALSE)
  }

  # ---- spatial kernel (degenerate cutoffs, item 21) ----
  neg <- .new_neg_counter()
  spatial_A <- NULL
  if (se == "spatial" && (!.point_only || lasso)) {
    sa <- .spatial_adjacency(lat_arr, lon_arr, spatial_cutoff, kernel = spatial_kernel, metric = spatial_metric)
    spatial_A <- sa$A
    share <- sa$nnz / n / n
    if (spatial_kernel == "uniform" && sa$nnz >= n * n)
      stop("the spatial cutoff covers all pairs; Conley SEs are degenerate (0 by construction). Choose a smaller spatial_cutoff.", call. = FALSE)
    if (share > 0.5)
      .inform(sprintf("  Warning: the spatial cutoff puts %.0f%% of all pairs within the cutoff (average neighbor share > 50%%); Conley SEs are unreliable in this regime. Consider a smaller cutoff.\n", 100 * share))
  }
  vc <- list(cl = if (se_v == "cluster") cl else NULL, A = if (se_v == "spatial") spatial_A else NULL, neg = neg)

  y_raw <- y; t_raw <- t; m_raw <- m; z_mat_raw <- z_mat; x_mat_raw <- x_mat

  # ---- LASSO selection (WP2) ----
  lasso_info <- NULL
  if (lasso) {
    if (is.null(x_mat_raw) || ncol(x_mat_raw) == 0L) {
      lasso_info <- list(S = logical(0), S_Y = logical(0), S_T = logical(0), S_M = logical(0), S_Z = logical(0),
                         S_Z_per_inst = list(), details = list(Y = NULL, T = NULL, M = NULL, Z = list()),
                         union_names = character(0), n_selected_union = 0L, n_total = 0L)
      x_mat_raw <- NULL; X <- character(0); control_cols <- character(0)
    } else {
      lasso_info <- .quadruple_lasso(y_raw, t_raw, m_raw, z_mat_raw, W_raw = x_mat_raw, col_names = rich_col_names,
                                     weights = w_vec, cluster_codes = if (!is.null(cl)) cl$codes else NULL,
                                     keep_idx = keep_idx_lasso, lambda_method = lasso_penalty,
                                     n_iter_max = as.integer(lasso_n_iter_max), tol = lasso_tol, c_slack = lasso_c_slack,
                                     alpha_p = lasso_alpha_p, lambda_override = lasso_lambda_override,
                                     spatial_A = if (se == "spatial") spatial_A else NULL, verbose = lasso_verbose)
      sel_idx <- which(lasso_info$S)
      if (length(sel_idx) == 0L) { x_mat_raw <- NULL; X <- character(0); control_cols <- character(0) }
      else { x_mat_raw <- x_mat_raw[, sel_idx, drop = FALSE]; X <- rich_col_names[sel_idx]; control_cols <- rich_in_df[sel_idx] }
      .inform(sprintf("  LASSO selected %d of %d candidate variables (%s penalty).\n",
                  lasso_info$n_selected_union, lasso_info$n_total, lasso_penalty))
    }
  }

  # ---- collinear controls (WP1.6) ----
  dropped_collinear <- character(0)
  if (!is.null(x_mat_raw) && ncol(x_mat_raw) > 0L) {
    gi <- .greedy_independent(x_mat_raw * sqrt_w, X, keep_first = sqrt_w, ref_norms = ref_norm_ctrl[control_cols])
    dropped_collinear <- gi$dropped
    if (length(dropped_collinear) > 0L) {
      .inform(sprintf("  Note: dropped (collinear): %s\n", paste(dropped_collinear, collapse = ", ")))
      x_mat_raw <- if (length(gi$kept)) x_mat_raw[, gi$kept, drop = FALSE] else NULL
      X <- X[gi$kept]; control_cols <- control_cols[gi$kept]
    }
  }

  # ---- sqrt(w) premultiplication (WLS as OLS on transformed data) ----
  y <- y_raw * sqrt_w; t <- t_raw * sqrt_w; m <- m_raw * sqrt_w
  z_mat <- z_mat_raw * sqrt_w
  ones <- sqrt_w
  if (!is.null(x_mat_raw)) {
    x_mat <- x_mat_raw * sqrt_w
    exog <- cbind(ones, x_mat); Z_base <- cbind(ones, z_mat, x_mat)
  } else {
    x_mat <- NULL; exog <- matrix(ones, ncol = 1L); Z_base <- cbind(ones, z_mat)
  }
  z_crit <- stats::qnorm(1 - alpha / 2)
  K_inst <- ncol(z_mat)
  n_exog <- ncol(exog)

  # ---- identification checks (paper eq. 26; A.1-6) ----
  # [1, X] and [1, X, T] are factorized once here and reused below (Q_zz of
  # both first stages, the RV regressions); qr_memo serves the 2SLS fits.
  qr_memo <- new.env(); qr_memo$entries <- list()
  q_exog <- .eq_qr(exog, need_Q = FALSE)$qr
  z_res_exog <- qr.resid(q_exog, z_mat)
  rel_resid_var <- function(r_, ref) if (ref > 0) sum(r_^2) / ref^2 else 0
  for (j_ in seq_len(K_inst))
    if (rel_resid_var(z_res_exog[, j_], ref_norm_role[[Z[j_]]]) < 1e-12)
      stop(sprintf("instrument %s has no variation after the controls/fixed effects (\u03c3_ZZ = 0, paper eq. 26).", Z[j_]), call. = FALSE)
  if (K_inst > 1L) {
    dz <- .greedy_independent(z_res_exog, Z, ref_norms = ref_norm_role[Z])$dropped
    if (length(dz) > 0L)
      stop(sprintf("instruments are collinear after the controls/fixed effects (%s).", paste(dz, collapse = ", ")), call. = FALSE)
  }
  if (rel_resid_var(qr.resid(q_exog, t), ref_norm_role[[T]]) < 1e-12)
    stop(sprintf("treatment %s has no variation after the controls/fixed effects.", T), call. = FALSE)
  if (rel_resid_var(qr.resid(q_exog, m), ref_norm_role[[M]]) < 1e-12)
    stop(sprintf("mediator %s has no variation after the controls/fixed effects.", M), call. = FALSE)
  q_exog_t <- .eq_qr(cbind(exog, t), need_Q = FALSE)$qr
  r_mt <- qr.resid(q_exog_t, m)
  if (sum(r_mt^2) / ref_norm_role[[M]]^2 < 1e-12)
    stop("the mediator is collinear with the treatment after the controls/fixed effects.", call. = FALSE)
  k_max <- n_exog + K_inst + 1L
  if (n - k_max - df_fe_absorbed < 1)
    stop(sprintf("not enough observations: n = %d, but the largest regression has %d coefficients plus %d absorbed fixed-effect degrees of freedom (n must exceed their sum).",
                 n, k_max, df_fe_absorbed), call. = FALSE)

  wald_over_K <- function(b, V) if (!.is_pd(V)) NaN else as.numeric(t(b) %*% solve(V, b)) / length(b)

  # 1. FIRST STAGE: T ~ Z + X (OLS)
  res_fs <- .twosls_with_vcov(t, Z_base, Z_base, se = se_v, vc = vc, df_adjust = df_vcov, memo = qr_memo)
  z_idx_fs <- 2:(1L + K_inst)
  beta_Z_fs <- res_fs$coef[z_idx_fs]
  V_ZZ_fs <- res_fs$vcov[z_idx_fs, z_idx_fs, drop = FALSE]
  se_Z_fs <- res_fs$se[z_idx_fs]
  pi1 <- as.numeric(beta_Z_fs); se_pi1 <- as.numeric(se_Z_fs)
  F_TZ <- wald_over_K(beta_Z_fs, V_ZZ_fs)
  Q_zz_fs <- crossprod(z_res_exog) / n
  F_TZ_op <- if (K_inst == 1L) F_TZ else if (.is_pd(V_ZZ_fs)) .op_effective_F(beta_Z_fs, V_ZZ_fs, Q_zz_fs) else NaN
  kc <- .op_keff_critval(V_ZZ_fs, Q_zz_fs, op_tau, op_alpha); K_eff_TZ <- kc[1]; F_TZ_op_critval <- kc[2]

  # 2. TREATMENT-MEDIATOR: M ~ T | Z + X
  X_gamma <- cbind(exog, t)
  res_gamma <- .twosls_with_vcov(m, X_gamma, Z_base, se = se_v, vc = vc, df_adjust = df_vcov, memo = qr_memo)
  idx_T_in_gamma <- ncol(X_gamma)
  gamma_hat <- res_gamma$coef[idx_T_in_gamma]; se_gamma <- res_gamma$se[idx_T_in_gamma]

  # 3. CONDITIONAL 2SLS: Y ~ M + T | Z + T + X
  X_theta <- cbind(exog, m, t); Z_theta <- cbind(Z_base, t)
  res_theta <- .twosls_with_vcov(y, X_theta, Z_theta, se = se_v, vc = vc, df_adjust = df_vcov, memo = qr_memo)
  idx_M <- n_exog + 1L; idx_T <- n_exog + 2L
  theta_hat <- res_theta$coef[idx_M]; tau_hat <- res_theta$coef[idx_T]
  se_theta <- res_theta$se[idx_M]; se_tau <- res_theta$se[idx_T]

  # 4. TOTAL EFFECT: Y ~ T | Z + X
  X_total <- cbind(exog, t)
  res_total <- .twosls_with_vcov(y, X_total, Z_base, se = se_v, vc = vc, df_adjust = df_vcov, memo = qr_memo)
  idx_T_total <- ncol(X_total)
  total_hat <- res_total$coef[idx_T_total]; se_total <- res_total$se[idx_T_total]

  # 5. CONDITIONAL FIRST STAGE: M ~ Z + T + X
  X_cfs <- cbind(exog, z_mat, t)
  res_cfs <- .twosls_with_vcov(m, X_cfs, X_cfs, se = se_v, vc = vc, df_adjust = df_vcov, memo = qr_memo)
  z_idx_cfs <- (n_exog + 1L):(n_exog + K_inst)
  beta_Z_cfs <- res_cfs$coef[z_idx_cfs]
  V_ZZ_cfs <- res_cfs$vcov[z_idx_cfs, z_idx_cfs, drop = FALSE]
  F_MZ_T <- wald_over_K(beta_Z_cfs, V_ZZ_cfs)
  F_MZ_T_individual <- vapply(z_idx_cfs, function(i) {
    v <- res_cfs$vcov[i, i]; if (isTRUE(v > 0)) res_cfs$coef[i]^2 / v else NaN }, numeric(1))
  Q_zz_cfs <- crossprod(qr.resid(q_exog_t, z_mat)) / n
  F_MZ_T_op <- if (K_inst == 1L) F_MZ_T else if (.is_pd(V_ZZ_cfs)) .op_effective_F(beta_Z_cfs, V_ZZ_cfs, Q_zz_cfs) else NaN
  kc <- .op_keff_critval(V_ZZ_cfs, Q_zz_cfs, op_tau, op_alpha); K_eff_MZ_T <- kc[1]; F_MZ_T_op_critval <- kc[2]

  # 6. Decomposition check and the delta-method SE of theta * gamma
  indirect <- theta_hat * gamma_hat
  decomp_check <- tau_hat + indirect
  decomp_error <- abs(decomp_check - total_hat)
  decomp_error_rel <- decomp_error / max(abs(total_hat), abs(tau_hat) + abs(indirect), 1e-300)
  cov_gamma_theta <- .cross_equation_cov(res_gamma, res_theta, idx_T_in_gamma, idx_M, se = se_v, vc = vc, df_adjust = df_vcov)
  var_indirect <- gamma_hat^2 * se_theta^2 + theta_hat^2 * se_gamma^2 + 2 * theta_hat * gamma_hat * cov_gamma_theta
  if (!is.na(var_indirect) && var_indirect < 0) {
    if (var_indirect > -1e-8 * max((gamma_hat * se_theta)^2, 1e-300)) var_indirect <- 0
    else {
      .inform(sprintf("  Warning: delta-method variance for the indirect effect is negative (%.2e). This typically reflects a poorly identified first stage. Reporting NaN; consider ivmediate_bootstrap for inference.\n", var_indirect))
      var_indirect <- NaN
    }
  }
  se_indirect <- if (!is.na(var_indirect)) sqrt(var_indirect) else NaN

  # 7. kappa-based sensitivity (D1, D4)
  X_partial <- if (!is.null(x_mat_raw) && ncol(x_mat_raw) > 0L) cbind(1, x_mat_raw) else matrix(1, n, 1)
  orth <- .wls_fit(X_partial, cbind(t_raw, m_raw, z_mat_raw), w_vec)$residuals   # one factorization
  t_orth <- orth[, 1]; m_orth <- orth[, 2]
  z_mat_orth <- orth[, 2L + seq_len(K_inst), drop = FALSE]
  eps_T_hat <- .wls_fit(cbind(1, z_mat_orth), t_orth, w_vec)$residuals
  eps_M_hat <- m_orth - gamma_hat * t_orth
  eps_Y_hat <- res_theta$resid / sqrt_w
  c_TM_hat <- .weighted_cov(eps_T_hat, eps_M_hat, w_vec)
  sigma_eT_hat <- .weighted_std(eps_T_hat, w_vec)
  sigma_eY_hat <- .weighted_std(eps_Y_hat, w_vec)
  kappa_star <- theta_hat * c_TM_hat
  kappa_star_CI <- (theta_hat - sign(theta_hat) * z_crit * se_theta) * c_TM_hat
  kappa_star_CI_informative <- if (is.finite(se_theta)) abs(theta_hat) > z_crit * se_theta else NA
  rho_TY_star <- if (isTRUE(sigma_eT_hat > 0 && sigma_eY_hat > 0)) kappa_star / (sigma_eT_hat * sigma_eY_hat) else NaN
  kappa_frac_grid <- c(0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00)
  if (is.na(kappa_star) || kappa_star == 0) { kappa_grid <- 0; kappa_frac_used <- 0 }
  else { kappa_grid <- kappa_frac_grid * kappa_star; kappa_frac_used <- kappa_frac_grid }
  kappa_adj_table <- lapply(seq_along(kappa_grid), function(ii) {
    k0 <- kappa_grid[ii]; th_adj <- theta_hat - k0 / c_TM_hat
    list(kappa = k0, kappa_frac = kappa_frac_used[ii], theta_adj = th_adj,
         ci_lo = th_adj - z_crit * se_theta, ci_hi = th_adj + z_crit * se_theta)
  })
  # Strength of the conditional first stage (App. F.4 step 3): TRUE iff the
  # effective F_MZ|T reaches the OP critical value at the user's op_tau; NA
  # when either is not finite. print.ivmediate notes it below the kappa table.
  F_MZ_T_above_crit <- if (is.finite(F_MZ_T_op) && is.finite(F_MZ_T_op_critval)) F_MZ_T_op >= F_MZ_T_op_critval else NA

  # 8. R^2-based sensitivity (App. F.3): t statistics under the SE regime,
  #    df = residual df of the OLS regressions (sensemakr's dof, item 24)
  W_fwl <- cbind(exog, t)
  rv_res <- qr.resid(q_exog_t, cbind(y, m, z_mat[, 1]))                # W_fwl = [1, X, T]; first instrument
  Y_tilde <- rv_res[, 1]; M_tilde <- rv_res[, 2]; Z_tilde <- rv_res[, 3]
  X_Z <- cbind(ones, Z_tilde)
  df_adj_rv <- (ncol(W_fwl) - 1L) + df_vcov
  res_rf <- .twosls_with_vcov(Y_tilde, X_Z, X_Z, se = se_v, vc = vc, df_adjust = df_adj_rv, memo = qr_memo)
  F_reduced_form <- (res_rf$coef[2] / res_rf$se[2])^2
  res_fs_cond <- .twosls_with_vcov(M_tilde, X_Z, X_Z, se = se_v, vc = vc, df_adjust = df_adj_rv, memo = qr_memo)
  F_first_stage_cond <- (res_fs_cond$coef[2] / res_fs_cond$se[2])^2
  df_rv <- n - ncol(W_fwl) - 1L - df_fe_absorbed
  # df_rv = 1 (n at its minimum) leaves no df for t*: RV and XRV are NaN, as in the Python reference
  t_star_rv <- if (df_rv > 1) abs(stats::qt(alpha / 2, df_rv - 1)) else NaN
  f_star_rv <- if (df_rv > 1) t_star_rv / sqrt(df_rv - 1) else NaN
  xrv_ols <- function(f_stat) { f_beta2 <- abs(f_stat) / df_rv; max((f_beta2 - f_star_rv^2) / (1 + f_beta2), 0) }
  rv_ols <- function(f_stat) {
    f_beta <- sqrt(abs(f_stat) / df_rv); fqa <- f_beta - f_star_rv
    if (!is.na(fqa) && fqa < 0) return(0)
    if (!is.na(fqa) && f_beta < 1 / f_star_rv) 0.5 * (sqrt(fqa^4 + 4 * fqa^2) - fqa^2) else xrv_ols(f_stat)
  }
  if (is.finite(F_reduced_form) && is.finite(F_first_stage_cond)) {
    xrv_rf <- xrv_ols(F_reduced_form); xrv_fs <- xrv_ols(F_first_stage_cond)
    rv_rf <- rv_ols(F_reduced_form); rv_fs <- rv_ols(F_first_stage_cond)
    XRV_theta <- min(xrv_rf, xrv_fs); RV_theta <- min(rv_rf, rv_fs)
  } else {
    xrv_rf <- xrv_fs <- rv_rf <- rv_fs <- XRV_theta <- RV_theta <- NaN
  }

  # 9. Confidence intervals
  ci_theta <- c(theta_hat - z_crit * se_theta, theta_hat + z_crit * se_theta)
  ci_tau   <- c(tau_hat   - z_crit * se_tau,   tau_hat   + z_crit * se_tau)
  ci_gamma <- c(gamma_hat - z_crit * se_gamma, gamma_hat + z_crit * se_gamma)
  ci_total <- c(total_hat - z_crit * se_total, total_hat + z_crit * se_total)
  ci_indirect <- if (is.finite(se_indirect)) c(indirect - z_crit * se_indirect, indirect + z_crit * se_indirect) else c(NaN, NaN)

  # 9b. Anderson-Rubin sets (theta: [1, X, T, Z]; gamma and tau_total: [1, Z, X])
  ar_theta <- ar_gamma <- ar_total <- NULL
  if (isTRUE(ar_ci)) {
    W_ar_theta <- cbind(exog, t, z_mat)
    z_cols_theta <- (n_exog + 2L):(n_exog + 1L + K_inst)
    ar_theta <- .ar_confidence_set(y, m, W_ar_theta, z_cols_theta, se_v, vc, df_vcov, alpha,
                                   center = theta_hat, scale = se_theta,
                                   grid_halfwidth = ar_grid_halfwidth, grid_points = ar_grid_points)
    z_cols_zb <- 2:(1L + K_inst)
    ar_gamma <- .ar_confidence_set(m, t, Z_base, z_cols_zb, se_v, vc, df_vcov, alpha,
                                   center = gamma_hat, scale = se_gamma,
                                   grid_halfwidth = ar_grid_halfwidth, grid_points = ar_grid_points)
    ar_total <- .ar_confidence_set(y, t, Z_base, z_cols_zb, se_v, vc, df_vcov, alpha,
                                   center = total_hat, scale = se_total,
                                   grid_halfwidth = ar_grid_halfwidth, grid_points = ar_grid_points)
  }

  # One warning per analysis for negative Conley variances of reported quantities (A.1-5)
  if (se == "spatial" && !.point_only) {
    rep_var <- c(res_gamma$vcov[idx_T_in_gamma, idx_T_in_gamma], res_theta$vcov[idx_M, idx_M],
                 res_theta$vcov[idx_T, idx_T], res_total$vcov[idx_T_total, idx_T_total],
                 res_rf$vcov[2, 2], res_fs_cond$vcov[2, 2])
    if (any(rep_var < 0) || !(.is_pd(V_ZZ_fs) && .is_pd(V_ZZ_cfs)))
      .inform("  Warning: the spatial (Conley) variance was negative for some reported coefficients (the kernel matrix is not positive semi-definite); their SEs and the statistics that use them are reported as NaN.\n")
  }

  u <- function(v) as.numeric(unname(v))
  structure(list(
    version = IVMED_VERSION,
    n = n, n_instruments = K_inst, n_covariates = length(X),
    n_dummies = if (length(dummies) > 0L) length(dummy_col_names) else 0L,
    se_type = se, alpha = alpha, fe = fe_cols, n_fe_groups = n_fe_groups,
    n_dropped_singletons = n_dropped_singletons, df_fe_absorbed = df_fe_absorbed,
    fe_nested_in_cluster = fe_nested_in_cluster, n_clusters = n_clusters,
    dropped_collinear = dropped_collinear,
    weighted = is_weighted, wtype = if (is_weighted) wtype else NULL,
    spatial_cutoff = if (se == "spatial") spatial_cutoff else NULL,
    spatial_kernel = if (se == "spatial") spatial_kernel else NULL,
    spatial_metric = if (se == "spatial") spatial_metric else NULL,
    theta_hat = u(theta_hat), se_theta = u(se_theta), ci_theta = u(ci_theta),
    tau_hat = u(tau_hat), se_tau = u(se_tau), ci_tau = u(ci_tau),
    gamma_hat = u(gamma_hat), se_gamma = u(se_gamma), ci_gamma = u(ci_gamma),
    total_hat = u(total_hat), se_total = u(se_total), ci_total = u(ci_total),
    indirect = u(indirect), se_indirect = u(se_indirect), ci_indirect = u(ci_indirect),
    cov_gamma_theta = u(cov_gamma_theta),
    pi1 = u(pi1), se_pi1 = u(se_pi1),
    decomp_check = u(decomp_check), decomp_error = u(decomp_error), decomp_error_rel = u(decomp_error_rel),
    F_TZ = u(F_TZ), F_MZ_T = u(F_MZ_T), F_MZ_T_individual = u(F_MZ_T_individual),
    F_TZ_op = u(F_TZ_op), F_MZ_T_op = u(F_MZ_T_op), K_eff_TZ = u(K_eff_TZ), K_eff_MZ_T = u(K_eff_MZ_T),
    F_TZ_op_critval = u(F_TZ_op_critval), F_MZ_T_op_critval = u(F_MZ_T_op_critval),
    op_tau = op_tau, op_alpha = op_alpha,
    ar_ci = ar_ci, ar_theta = ar_theta, ar_gamma = ar_gamma, ar_total = ar_total,
    c_TM_hat = u(c_TM_hat),
    rho_TM_hat = u(c_TM_hat),     # deprecated alias (pre-2026-07 name for Cov(eps_T, eps_M))
    sigma_eT = u(sigma_eT_hat), sigma_eY = u(sigma_eY_hat),
    kappa_star = u(kappa_star), kappa_star_CI = u(kappa_star_CI),
    kappa_star_CI_informative = kappa_star_CI_informative,
    rho_TY_star = u(rho_TY_star),
    kappa_adj_table = kappa_adj_table, F_MZ_T_above_crit = F_MZ_T_above_crit,
    F_reduced_form = u(F_reduced_form), F_first_stage_cond = u(F_first_stage_cond),
    f_reduced_form = u(F_reduced_form),            # deprecated alias (= t^2, not Cohen's f)
    f_first_stage_cond = u(F_first_stage_cond),    # deprecated alias (= t^2, not Cohen's f)
    RV_theta = u(RV_theta), XRV_theta = u(XRV_theta),
    rv_rf = u(rv_rf), rv_fs = u(rv_fs), xrv_rf = u(xrv_rf), xrv_fs = u(xrv_fs),
    lasso = lasso, lasso_penalty = if (lasso) lasso_penalty else NULL,
    lasso_seed = NULL,                             # the CV path is not ported
    lasso_info = lasso_info, lasso_expansion_meta = rich_meta,
    lasso_show_selection = if (lasso) show_selection else FALSE,
    rich_col_names = if (lasso) rich_col_names else NULL,
    eps_T_hat = u(eps_T_hat), eps_M_hat = u(eps_M_hat), eps_Y_hat = u(eps_Y_hat)
  ), class = c("ivmediate", "list"))
}


# =============================================================================
#  Significance stars and printing
#  (Greek letters and arrows are written as \uXXXX escapes so the source file
#   stays pure ASCII.)
# =============================================================================

.validate_star_cutoffs <- function(cutoffs) {
  cuts <- suppressWarnings(as.numeric(cutoffs))
  if (length(cuts) != 3L || any(is.na(cuts)))
    stop(sprintf("star_cutoffs must be 3 numbers; got %s", paste(cutoffs, collapse = ", ")), call. = FALSE)
  if (!all(is.finite(cuts) & cuts > 0 & cuts < 1))
    stop(sprintf("star_cutoffs must all be strictly between 0 and 1; got %s", paste(cuts, collapse = ", ")), call. = FALSE)
  sort(cuts)   # ascending: [1] = *** , [2] = ** , [3] = *
}

# Two-sided Wald p-value from z = est/se. NA if not both finite or se <= 0.
.wald_pvalue <- function(est, se) {
  if (is.null(est) || is.null(se)) return(NA_real_)
  if (!(is.finite(est) && is.finite(se)) || se <= 0) return(NA_real_)
  2.0 * stats::pnorm(abs(est / se), lower.tail = FALSE)
}

.sig_stars <- function(pval, cutoffs) {
  if (is.null(pval) || !is.finite(pval)) return("")
  if (pval < cutoffs[1]) return("***")
  if (pval < cutoffs[2]) return("**")
  if (pval < cutoffs[3]) return("*")
  ""
}

.format_cut <- function(c) formatC(c, format = "g")

# Python-style %g of a number (used for tau and alpha in the output).
.fmt_g <- function(x) sprintf("%g", x)

# Character-count padding (sprintf pads by bytes; multi-byte labels would be short).
.pad_left_align <- function(s, w) { d <- w - nchar(s, type = "chars"); if (d > 0) paste0(s, strrep(" ", d)) else s }
.pad_right_align <- function(s, w) { d <- w - nchar(s, type = "chars"); if (d > 0) paste0(strrep(" ", d), s) else s }

# kappa_table = TRUE also prints the kappa-adjusted estimates (theta and its CI
# with kappa_0 at 0, 10, 25, 50, 75, 90 and 100% of the breakeven kappa*); when
# the effective F_{M,Z|T,X} is below the Olea-Pflueger critical value
# (F_MZ_T_above_crit FALSE), a note below the table says that c_TM is then
# estimated imprecisely. r2_sensitivity = TRUE also prints the R^2-based
# section (reduced-form F, RV, XRV and their components). The values are in
# the results either way (Addendum A.2).
#' Print an IV mediation analysis
#'
#' Prints the mediation estimates, the first-stage diagnostics, the
#' Anderson-Rubin confidence sets and the kappa-based sensitivity summary of an
#' [ivmediate()] result and, on request, the kappa-adjusted table and the
#' R^2-based sensitivity section. The values are in the result either way.
#'
#' @param x An object of class `ivmediate`, the result of [ivmediate()].
#' @param title Title line of the printout. Default: `"IV Mediation Analysis"`.
#' @param boot Optional result of [ivmediate_bootstrap()] for the same model.
#'   The coefficient table then shows the bootstrap SEs and CIs, and the kappa
#'   section adds the bootstrap CIs of kappa* and rho*_TY. Default: `NULL`.
#' @param stars Print significance stars. Default: `TRUE`.
#' @param star_cutoffs Three p-value cutoffs in (0, 1) for `***`, `**` and `*`,
#'   in any order. Default: `c(0.01, 0.05, 0.10)`.
#' @param kappa_table Also print the kappa-adjusted estimates: theta and its CI
#'   with kappa_0 at 0, 10, 25, 50, 75, 90 and 100 percent of the breakeven
#'   kappa*. When the effective conditional first-stage F is below the
#'   Olea-Pflueger critical value (`F_MZ_T_above_crit` is `FALSE`), a note below
#'   the table says that c_TM is then estimated imprecisely. Default: `FALSE`.
#' @param r2_sensitivity Also print the R^2-based sensitivity section: the
#'   reduced-form F (= t^2), RV(theta) and XRV(theta) with their reduced-form
#'   and conditional first-stage components. Default: `FALSE`.
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#' @seealso [ivmediate()], [ivmediate_bootstrap()]
#' @examples
#' bw <- read.csv(system.file("extdata", "becker_woessmann_2009.csv",
#'                            package = "ivmediate"))
#' res <- ivmediate(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw",
#'                  Z = "kmwittenberg",
#'                  X = "f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf
#'                       f_dumb hhsize lnpop gpop")
#' print(res, title = "Becker and Woessmann (2009)", kappa_table = TRUE)
#' @export
print.ivmediate <- function(x, title = "IV Mediation Analysis", boot = NULL,
                            stars = TRUE, star_cutoffs = c(0.01, 0.05, 0.10),
                            kappa_table = FALSE, r2_sensitivity = FALSE, ...) {
  r <- x
  .check_bool("stars", stars)
  .check_bool("kappa_table", kappa_table)
  .check_bool("r2_sensitivity", r2_sensitivity)
  if (stars) star_cutoffs <- .validate_star_cutoffs(star_cutoffs)
  lvl <- .level_label(r$alpha)
  using_boot <- !is.null(boot)
  z_crit <- stats::qnorm(1 - r$alpha / 2)

  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf("  %s\n", title))
  cat(sprintf("%s\n", strrep("=", 70)))
  info_line <- sprintf("  n = %s  |  Instruments: %d  |  Covariates: %d  |  SE: %s",
                       formatC(r$n, format = "d", big.mark = ","), r$n_instruments, r$n_covariates, r$se_type)
  if (!is.null(r$n_clusters)) info_line <- paste0(info_line, sprintf(" (G = %d clusters)", r$n_clusters))
  if (using_boot) info_line <- paste0(info_line, sprintf("  +  Bootstrap (B=%d, %s)", as.integer(boot$B), boot$method))
  cat(info_line, "\n", sep = "")
  if (isTRUE(r$weighted)) cat(sprintf("  Weighted: %s\n", r$wtype))
  if (!is.null(r$n_dummies) && r$n_dummies > 0) cat(sprintf("  Dummy indicators: %d\n", r$n_dummies))
  if (length(r$dropped_collinear) > 0L) cat(sprintf("  Dropped (collinear): %s\n", paste(r$dropped_collinear, collapse = ", ")))
  if (!is.null(r$fe)) {
    line <- sprintf("  FE: %s  |  Groups: %d  |  Singletons dropped: %d  |  df absorbed: %d",
                    paste(r$fe, collapse = ", "), r$n_fe_groups, r$n_dropped_singletons, r$df_fe_absorbed)
    nested <- r$fe_nested_in_cluster
    if (!is.null(nested) && any(nested))
      line <- paste0(line, "  |  nested in clusters: ", paste(r$fe[nested], collapse = ", "), " (not counted in the CR1 df)")
    cat(line, "\n", sep = "")
  }
  if (identical(r$se_type, "spatial"))
    cat(sprintf("  Spatial: cutoff = %s (metric=%s, kernel=%s)\n", as.character(r$spatial_cutoff),
                as.character(r$spatial_metric), as.character(r$spatial_kernel)))
  if (isTRUE(r$lasso)) {
    n_sel <- if (!is.null(r$lasso_info)) r$lasso_info$n_selected_union else 0L
    n_tot <- if (!is.null(r$lasso_info)) r$lasso_info$n_total else 0L
    cat(sprintf("  LASSO: post-double selection, penalty='%s', selected %d/%d candidates\n", r$lasso_penalty, n_sel, n_tot))
  }
  cat(sprintf("%s\n", strrep("\u2500", 70)))

  # ---- MEDIATION ESTIMATES ----
  cat("\n  MEDIATION ESTIMATES\n")
  if (using_boot) { se_header <- "Boot SE"; ci_header <- sprintf("[Boot %s CI]", lvl) }
  else { se_header <- "SE"; ci_header <- sprintf("[%s CI]", lvl) }
  cat(sprintf("  %-20s %10s %10s %24s\n", "Parameter", "Estimate", se_header, ci_header))
  cat(sprintf("  %s\n", strrep("\u2500", 64)))
  if (using_boot) {
    rows <- list(
      list("\u03c4_total (Wald)",      r$total_hat, boot$boot_se[["total_hat"]],    boot$boot_ci[["total_hat"]]),
      list("\u03c4 (direct T\u2192Y)", r$tau_hat,   boot$boot_se[["tau_hat"]],      boot$boot_ci[["tau_hat"]]),
      list("\u03b8 (mediator\u2192Y)", r$theta_hat, boot$boot_se[["theta_hat"]],    boot$boot_ci[["theta_hat"]]),
      list("\u03b3 (T\u2192mediator)", r$gamma_hat, boot$boot_se[["gamma_hat"]],    boot$boot_ci[["gamma_hat"]]),
      list("\u03b8\u03b3 (indirect)",  r$indirect,  boot$boot_se[["indirect_hat"]], boot$boot_ci[["indirect_hat"]]))
  } else {
    rows <- list(
      list("\u03c4_total (Wald)",      r$total_hat, r$se_total,    r$ci_total),
      list("\u03c4 (direct T\u2192Y)", r$tau_hat,   r$se_tau,      r$ci_tau),
      list("\u03b8 (mediator\u2192Y)", r$theta_hat, r$se_theta,    r$ci_theta),
      list("\u03b3 (T\u2192mediator)", r$gamma_hat, r$se_gamma,    r$ci_gamma),
      list("\u03b8\u03b3 (indirect)",  r$indirect,  r$se_indirect, r$ci_indirect))
  }
  for (row in rows) {
    name <- row[[1]]; est <- row[[2]]; se_val <- row[[3]]; ci <- row[[4]]
    if (!is.null(se_val) && is.finite(se_val)) {
      star_str <- if (stars) .sig_stars(.wald_pvalue(est, se_val), star_cutoffs) else ""
      suffix <- if (nzchar(star_str)) paste0(" ", star_str) else ""
      cat(sprintf("  %s %10.4f %10.4f   [%9.4f, %9.4f]%s\n", .pad_left_align(name, 20), est, se_val, ci[1], ci[2], suffix))
    } else {
      cat(sprintf("  %s %10.4f\n", .pad_left_align(name, 20), est))
    }
  }
  if (stars)
    cat(sprintf("  Sig. codes: *** p<%s, ** p<%s, * p<%s\n",
                .format_cut(star_cutoffs[1]), .format_cut(star_cutoffs[2]), .format_cut(star_cutoffs[3])))

  # ---- decomposition check (silent unless violated; relative, item 26) ----
  rel <- r$decomp_error_rel
  if (is.null(rel)) rel <- r$decomp_error / max(abs(r$total_hat), abs(r$tau_hat) + abs(r$indirect), 1e-300)
  if (r$n_instruments == 1L) {
    if (isTRUE(rel > 1e-8)) {
      cat("\n  WARNING: decomposition identity violated.\n")
      cat(sprintf("  \u03c4 + \u03b8\u03b3 = %.6f  vs  \u03c4_total = %.6f\n", r$decomp_check, r$total_hat))
      cat(sprintf("  relative error = %.2e  (threshold = 1e-08)\n", rel))
      cat("  For a single instrument this identity is exact; a nonzero error indicates an internal inconsistency. Please report.\n")
    }
  } else if (isTRUE(rel > 1e-2)) {
    cat(sprintf("\n  Note: with K>1 instruments, \u03c4 + \u03b8\u03b3 = %.4f and the Wald\n", r$decomp_check))
    cat(sprintf("  \u03c4_total = %.4f differ by %.1e (relative). This is expected \u2014 the\n", r$total_hat, rel))
    cat("  objects are estimated equation-by-equation (not by joint GMM) and coincide\n")
    cat("  only asymptotically. Not an error.\n")
  }

  # ---- FIRST-STAGE DIAGNOSTICS ----
  cat(sprintf("\n  FIRST-STAGE DIAGNOSTICS  (SE: %s)\n", r$se_type))
  if (!is.null(r$pi1)) {
    body <- paste(sprintf("%.4f (%.4f)", r$pi1, r$se_pi1), collapse = ", ")
    # K > 1: the label column is two characters wider (the widest label is
    # "F_{M, Z|T, X} (OP eff.)")
    pad <- if (r$n_instruments > 1L) "        " else "      "
    cat(sprintf("  \u03c0_1 (T on Z | X)%s= %s\n", pad, body))
  }
  op_str <- function(crit, keff) {
    if (is.null(crit) || !is.finite(crit)) return("  |  OP critical value: unavailable")
    sprintf("  |  OP critical value = %.2f  (K_eff=%.2f, \u03c4=%s, \u03b1=%s)", crit, keff, .fmt_g(r$op_tau), .fmt_g(r$op_alpha))
  }
  op_TZ <- op_str(r$F_TZ_op_critval, r$K_eff_TZ)
  if (r$n_instruments > 1L) {
    cat(sprintf("  F_{T, Z|X}              = %.2f  (joint Wald)\n", r$F_TZ))
    cat(sprintf("  F_{T, Z|X} (OP eff.)    = %.2f%s\n", r$F_TZ_op, op_TZ))
  } else cat(sprintf("  F_{T, Z|X}            = %.2f%s\n", r$F_TZ, op_TZ))
  op_MZ <- op_str(r$F_MZ_T_op_critval, r$K_eff_MZ_T)
  if (r$n_instruments > 1L) {
    cat(sprintf("  F_{M, Z|T, X}           = %.2f  (joint Wald)\n", r$F_MZ_T))
    cat(sprintf("  F_{M, Z|T, X} (OP eff.) = %.2f%s\n", r$F_MZ_T_op, op_MZ))
    cat(sprintf("    per-instrument t\u00b2:      %s\n", paste(sprintf("%.1f", r$F_MZ_T_individual), collapse = "  ")))
  } else cat(sprintf("  F_{M, Z|T, X}         = %.2f%s\n", r$F_MZ_T, op_MZ))

  # ---- Anderson-Rubin confidence sets ----
  if (isTRUE(r$ar_ci) && !is.null(r$ar_theta)) {
    cat(sprintf("\n  ANDERSON\u2013RUBIN %s CONFIDENCE SETS  (weak-IV-robust, SE: %s)\n", lvl, r$se_type))
    any_unbounded <- FALSE
    for (row in list(list("\u03b8 (M \u2192 Y | T)", r$ar_theta), list("\u03b3 (T \u2192 M)", r$ar_gamma),
                     list("\u03c4_total", r$ar_total))) {
      ar <- row[[2]]
      if (is.null(ar)) next
      cat(sprintf("  %s %s\n", .pad_left_align(row[[1]], 14), .format_ar_set(ar)))
      pv0 <- ar$pval_at_zero
      pv0_s <- if (is.finite(pv0)) sprintf("%.4f", pv0) else "n/a"
      cat(sprintf("  %-14s   AR p-value at 0 = %s   (%s)\n", "", pv0_s,
                  if (identical(ar$method, "fieller")) "exact Fieller inversion" else "exact statistic, grid search"))
      if (.ar_unbounded_or_disconnected(ar)) any_unbounded <- TRUE
    }
    if (any_unbounded) {
      cat("  Note: an unbounded or disconnected AR set is the expected outcome under weak\n")
      cat("  identification \u2014 the data cannot bound the parameter at this confidence level.\n")
    }
    cat("  AR sets remain valid regardless of first-stage strength; compare with the Wald\n")
    cat("  CIs above, which require the corresponding effective F to clear the OP threshold.\n")
  }

  # ---- kappa-BASED SENSITIVITY ----
  cat("\n  \u03ba-BASED SENSITIVITY (MCA diagnostics)\n")
  if (isTRUE(r$weighted)) cat(sprintf("  (Hajek-weighted moments, wtype = %s)\n", r$wtype))
  if (using_boot) cat("  (\u03ba-frame CIs use bootstrap SE on \u03b8)\n")
  cat(sprintf("  c_TM = Cov(\u03b5_T, \u03b5_M) = %.4f\n", r$c_TM_hat))
  cat(sprintf("  \u03ba* (breakeven)       = %.4f\n", r$kappa_star))
  se_theta_use <- if (using_boot) boot$boot_se[["theta_hat"]] else r$se_theta
  th <- r$theta_hat
  kappa_star_CI_use <- (th - sign(th) * z_crit * se_theta_use) * r$c_TM_hat
  informative <- is.finite(se_theta_use) && abs(th) > z_crit * se_theta_use
  flag <- if (informative) "" else "  (uninformative: the CI already contains 0)"
  cat(sprintf("  \u03ba*_CI (inf. break.)  = %.4f%s\n", kappa_star_CI_use, flag))
  rho <- r$rho_TY_star
  rho_note <- if (is.finite(rho) && abs(rho) > 1) "  (|\u03c1*| > 1: no feasible error correlation overturns \u03b8)" else ""
  cat(sprintf("  \u03c1*_TY (as corr.)     = %.4f%s\n", rho, rho_note))
  if (using_boot) {
    for (kl in list(c("kappa_star", "\u03ba*"), c("rho_TY_star", "\u03c1*_TY"))) {
      ci_b <- boot$boot_ci[[kl[1]]]
      if (!is.null(ci_b))
        cat(sprintf("  %s bootstrap %s CI (%s) = [%.4f, %.4f]\n", kl[2], .level_label(boot$alpha), boot$ci_type, ci_b[1], ci_b[2]))
    }
  }
  if (kappa_table) {
    cat("\n  \u03ba-adjusted estimates  (\u03ba_0 swept as fraction of breakeven \u03ba*):\n")
    cat(sprintf("  %s %s %s %s %s\n", .pad_right_align("% \u03ba*", 6), .pad_right_align("\u03ba_0", 12),
                .pad_right_align("\u03b8_adj", 10), .pad_right_align("CI_lo", 10), .pad_right_align("CI_hi", 10)))
    for (row in r$kappa_adj_table) {
      pct <- (if (!is.null(row$kappa_frac)) row$kappa_frac else 0) * 100
      if (using_boot) { ci_lo_use <- row$theta_adj - z_crit * se_theta_use; ci_hi_use <- row$theta_adj + z_crit * se_theta_use }
      else { ci_lo_use <- row$ci_lo; ci_hi_use <- row$ci_hi }
      cat(sprintf("  %5.0f%% %12.4f %10.4f %10.4f %10.4f\n", pct, row$kappa + 0, row$theta_adj, ci_lo_use, ci_hi_use))
    }
    # A conditional first stage below the OP critical value: c_TM, and so the
    # kappa quantities, are estimated imprecisely (App. F.4 step 3).
    if (identical(r$F_MZ_T_above_crit, FALSE)) {
      f_lab <- if (r$n_instruments == 1L) "F_{M, Z|T, X}" else "F_{M, Z|T, X} (OP eff.)"
      ar_shown <- isTRUE(r$ar_ci) && !is.null(r$ar_theta)
      cat(sprintf("  Note: %s = %.2f is below the OP critical value (%.2f, \u03c4 = %s), so c_TM is estimated\n",
                  f_lab, r$F_MZ_T_op, r$F_MZ_T_op_critval, .fmt_g(r$op_tau)))
      cat(if (ar_shown) "  imprecisely; the AR set for \u03b8 above is robust to weak instruments.\n" else "  imprecisely.\n")
    }
  }

  # ---- R^2-BASED SENSITIVITY (on request) ----
  if (r2_sensitivity) {
    cat("\n  R\u00b2-BASED SENSITIVITY (reduced-form robustness)\n")
    cat(sprintf("  (t-statistics computed under %s SEs%s\n", r$se_type,
                if (r$n_instruments > 1L) "; with K > 1, RV/XRV use the first instrument)" else ")"))
    F_rf <- if (!is.null(r$F_reduced_form)) r$F_reduced_form else r$f_reduced_form
    cat(sprintf("  F (reduced form, = t\u00b2) = %.1f\n", F_rf))
    cat(sprintf("  RV(\u03b8)                  = %.4f\n", r$RV_theta))
    cat(sprintf("  XRV(\u03b8)                 = %.4f\n", r$XRV_theta))
    cat(sprintf("    (RV_rf=%.4f, RV_fs=%.4f)\n", r$rv_rf, r$rv_fs))
  }

  # ---- LASSO SELECTION DETAILS ----
  if (isTRUE(r$lasso) && isTRUE(r$lasso_show_selection) && !is.null(r$lasso_info)) {
    info <- r$lasso_info; details <- info$details
    cat("\n  LASSO SELECTION DETAILS\n")
    cat(sprintf("  Penalty method: %s\n", r$lasso_penalty))
    print_set <- function(label, d, max_show = 20L) {
      if (is.null(d)) return(invisible(NULL))
      n_sel <- d$n_selected; k_it <- d$n_iter
      conv <- if (isTRUE(d$converged)) sprintf("converged in %d iteration%s", k_it, if (k_it != 1) "s" else "")
              else sprintf("did not converge after %d iteration%s", k_it, if (k_it != 1) "s" else "")
      cat(sprintf("  S_%s  (%3d selected, \u03bb = %.3f, %s)\n", label, n_sel, d$lambda_value, conv))
      if (n_sel > 0L) {
        nms <- utils::head(d$names, max_show); more <- n_sel - length(nms)
        line <- paste(nms, collapse = ", "); if (more > 0L) line <- paste0(line, sprintf(", ... (+%d more)", more))
        cat(sprintf("    %s\n", line))
      }
    }
    print_set("Y", details$Y); print_set("T", details$T); print_set("M", details$M)
    if (!is.null(details$Z)) for (k in seq_along(details$Z)) print_set(sprintf("Z[%d]", k - 1L), details$Z[[k]])
    n_union <- info$n_selected_union; n_total <- info$n_total; union_names <- info$union_names
    cat(sprintf("  Union S = S_Y \u222a S_T \u222a S_M \u222a S_Z: %d of %d candidates\n", n_union, n_total))
    if (length(union_names) > 0L) {
      shown <- utils::head(union_names, 20L); more <- n_union - length(shown)
      line <- paste(shown, collapse = ", "); if (more > 0L) line <- paste0(line, sprintf(", ... (+%d more)", more))
      cat(sprintf("    %s\n", line))
    }
  }
  cat(sprintf("\n%s\n\n", strrep("=", 70)))
  invisible(x)
}

#' Print a bootstrap of an IV mediation analysis
#'
#' Prints the bootstrap point estimates, SEs and confidence intervals of the
#' five effects and of kappa* and rho*_TY, next to the analytic SEs and CIs.
#'
#' @param x An object of class `ivmediate_bootstrap`, the result of
#'   [ivmediate_bootstrap()].
#' @param title Title line of the printout. Default: `"Bootstrap Mediation"`.
#' @param stars,star_cutoffs Significance stars and their cutoffs, as in
#'   [print.ivmediate()].
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#' @seealso [ivmediate_bootstrap()]
#' @export
print.ivmediate_bootstrap <- function(x, title = "Bootstrap Mediation", stars = TRUE,
                                      star_cutoffs = c(0.01, 0.05, 0.10), ...) {
  boot <- x
  .check_bool("stars", stars)
  if (stars) star_cutoffs <- .validate_star_cutoffs(star_cutoffs)
  coef_keys <- c("theta_hat", "tau_hat", "gamma_hat", "total_hat", "indirect_hat")
  cat(sprintf("\n%s\n", strrep("=", 70)))
  cat(sprintf("  %s  (method=%s, B=%d, seed=%s)\n", title, boot$method, as.integer(boot$B), format(boot$seed, scientific = FALSE)))
  cat(sprintf("%s\n", strrep("=", 70)))
  cat(sprintf("  Bootstrap CI type: %s  (%s, \u03b1 = %s)\n", boot$ci_type, .level_label(boot$alpha), .fmt_g(boot$alpha)))
  failed_str <- if (boot$n_failed > 0) sprintf(", %d failed", boot$n_failed) else ""
  cat(sprintf("  %d successful replicates%s.\n", as.integer(boot$B - boot$n_failed), failed_str))
  cat(sprintf("  %s\n", strrep("\u2500", 68)))
  label_map <- c(theta_hat = "\u03b8 (mediator\u2192Y)", tau_hat = "\u03c4 (direct T\u2192Y)",
                 gamma_hat = "\u03b3 (T\u2192mediator)", total_hat = "\u03c4_total",
                 indirect_hat = "\u03b8\u03b3 (indirect)", kappa_star = "\u03ba*",
                 rho_TY_star = "\u03c1*_TY", RV_theta = "RV(\u03b8)", XRV_theta = "XRV(\u03b8)")
  cat(sprintf("  %-18s %10s %10s %26s\n", "Parameter", "Estimate", "Boot SE", "Boot CI"))
  for (p in names(boot$point)) {
    label <- if (p %in% names(label_map)) label_map[[p]] else p
    est <- boot$point[[p]]; sdv <- boot$boot_se[[p]]; ci <- boot$boot_ci[[p]]
    ci_str <- sprintf("[%10.4f, %10.4f]", ci[1], ci[2])
    star_str <- if (stars && (p %in% coef_keys)) .sig_stars(.wald_pvalue(est, sdv), star_cutoffs) else ""
    suffix <- if (nzchar(star_str)) sprintf(" %s", star_str) else ""
    cat(sprintf("  %s %10.4f %10.4f %26s%s\n", .pad_left_align(label, 18), est, sdv, ci_str, suffix))
  }
  if (stars)
    cat(sprintf("  Sig. codes: *** p<%s, ** p<%s, * p<%s\n",
                .format_cut(star_cutoffs[1]), .format_cut(star_cutoffs[2]), .format_cut(star_cutoffs[3])))
  cat(sprintf("%s\n\n", strrep("=", 70)))
  invisible(x)
}


# =============================================================================
#  Bootstrap inference (WP3)
# =============================================================================
#
# Resamples the estimation sample (rows for "pairs", clusters present in the
# estimation sample for "block"). Replicate streams come from one seed through
# L'Ecuyer-CMRG sub-streams (parallel::nextRNGStream); the caller's RNG kind and
# .Random.seed are restored on exit, so the user's random numbers are never
# changed. For basic and percentile CIs without diagnostics the replicates
# compute point estimates only. In the block bootstrap, each drawn copy of a
# cluster gets its own levels of every FE nested in the clusters.
#' Bootstrap inference for IV mediation
#'
#' Refits the mediation analysis of [ivmediate()] on `B` resampled data sets and
#' returns bootstrap SEs and confidence intervals for the five effects and for
#' kappa* and rho*_TY. Print the result with [print.ivmediate_bootstrap()], or
#' pass it to [print.ivmediate()] as `boot` to show the bootstrap inference in
#' the main table.
#'
#' Resampling uses the estimation sample: rows with missing values in any used
#' column, zero-weight rows and (with spatial SEs) rows without coordinates are
#' removed first, so every `"pairs"` replicate has n rows and every `"block"`
#' replicate has the G clusters of the estimation sample. Singleton FE groups
#' are detected again in each replicate. Replicates reuse the control lists
#' resolved on the original data. For `"basic"` and `"percentile"` intervals
#' without `bootstrap_diagnostics`, the replicates compute point estimates
#' only.
#'
#' @inheritParams ivmediate
#' @param method Resampling scheme. `"pairs"` resamples observations with
#'   replacement (valid for independent observations only). `"block"`
#'   resamples clusters of `clusterid` with replacement (valid under cluster
#'   dependence; requires `clusterid`); each drawn copy is its own cluster,
#'   and fixed effects nested in the clusters get their own levels per copy,
#'   as with Stata's `idcluster()`. `NULL` (default) means `"block"` when
#'   `se = "cluster"` and `"pairs"` otherwise.
#' @param B Number of bootstrap replicates (an integer of at least 50).
#'   Default: 999.
#' @param seed Random seed. Replicate b draws from the b-th L'Ecuyer-CMRG
#'   stream of `seed` (`parallel::nextRNGStream`), so different seeds share no
#'   replicates. The caller's random-number generator kind and `.Random.seed`
#'   are restored on exit. Default: 20260317.
#' @param ci_type `"basic"` (2 theta-hat minus the bootstrap quantiles; the
#'   default, robust to skewness in products such as theta-hat times gamma-hat
#'   and in ratios), `"percentile"` (the bootstrap quantiles) or
#'   `"studentized"` (bootstrap-t with each replicate's analytic SE; kappa*,
#'   rho*_TY, RV and XRV, which have no analytic SE, use `"basic"`).
#' @param return_boot_dist If `TRUE`, the result includes the replicates of
#'   every bootstrapped quantity (`boot_dist`). Default: `FALSE`.
#' @param bootstrap_diagnostics If `TRUE`, RV(theta) and XRV(theta) are also
#'   bootstrapped. Default: `FALSE`.
#' @param progress Print a progress line every five percent of the replicates.
#'   Default: `TRUE`.
#' @param verbose Print extra diagnostic information. Default: `FALSE`.
#' @param boot_indices Testing hook: a `B` x n (`"pairs"`) or `B` x G
#'   (`"block"`) matrix of 0-based draws (rows of the estimation sample, or
#'   clusters in order of first appearance in it). The random-number generator
#'   is not used. The cross-language validation harness uses it to replay the
#'   same draws in every language. Default: `NULL`.
#' @param .return_indices Internal: also return the draws and the
#'   estimation-sample rows. Default: `FALSE`.
#'
#' @return A list of class `c("ivmediate_bootstrap", "list")` with `point`
#'   (the original-data estimates), `boot_se` (standard deviations of the
#'   replicates), `boot_ci` (bootstrap CIs of type `ci_type`), `analytic_se`
#'   and `analytic_ci` (from the original fit), `n_failed` (failed
#'   replicates), `n_est` (size of the resampled estimation sample),
#'   `method`, `B`, `seed`, `ci_type`, `alpha` and, with
#'   `return_boot_dist = TRUE`, `boot_dist`.
#'
#' @references
#' Ferrara, Andreas, Robert Gold, Stephan Heblich, and Rodrigo Pinto (2026),
#' "Identifying Causal Mediation Effects with a Single Instrument," NBER
#' Working Paper No. \[number\], doi:\[DOI\].
#'
#' @seealso [ivmediate()], [print.ivmediate_bootstrap()]
#' @examples
#' bw <- read.csv(system.file("extdata", "becker_woessmann_2009.csv",
#'                            package = "ivmediate"))
#' xs <- "f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize
#'        lnpop gpop"
#' boot <- ivmediate_bootstrap(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw",
#'                             Z = "kmwittenberg", X = xs, B = 99,
#'                             progress = FALSE)
#' print(boot)
#' res <- ivmediate(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw",
#'                  Z = "kmwittenberg", X = xs)
#' print(res, boot = boot)
#' @export
ivmediate_bootstrap <- function(data, Y, T, M, Z, X = NULL, dummies = NULL, fe = NULL,
                                se = "robust", clusterid = NULL, alpha = 0.05,
                                weights = NULL, wtype = "aweight",
                                op_tau = 0.10, op_alpha = 0.05,
                                lat = NULL, lon = NULL, spatial_cutoff = NULL,
                                spatial_kernel = NULL, spatial_metric = NULL,
                                method = NULL, B = 999, seed = 20260317,
                                ci_type = "basic", return_boot_dist = FALSE,
                                bootstrap_diagnostics = FALSE, progress = TRUE, verbose = FALSE,
                                # Testing hook: a B x (n_est or G) matrix of 0-based draws (rows of
                                # the estimation sample for "pairs", clusters in order of first
                                # appearance in it for "block"); the RNG is bypassed (harness).
                                boot_indices = NULL, .return_indices = FALSE) {
  if (!is.data.frame(data)) stop("data must be a data.frame (data.table and tibble inputs are accepted).", call. = FALSE)
  data <- as.data.frame(data)
  for (nm in c("return_boot_dist", "bootstrap_diagnostics", "progress", "verbose", ".return_indices"))
    .check_bool(nm, get(nm))
  .check_choice("se", se, .SE_TYPES)
  .check_open_unit("alpha", alpha)
  if (is.null(method)) method <- if (se == "cluster" && !is.null(clusterid)) "block" else "pairs"
  .check_choice("method", method, c("pairs", "block"))
  if (method == "block" && is.null(clusterid))
    stop("method='block' (cluster bootstrap) requires clusterid. Provide clusterid or switch to method='pairs' for iid data.", call. = FALSE)
  if (method == "pairs" && !is.null(clusterid) && se == "cluster")
    .inform(sprintf("  Note: method='pairs' was selected but clusterid='%s' is provided. Pairs bootstrap is valid only under cross-sectional independence; if your data are clustered, use method='block'.\n", clusterid))
  if (method == "block" && se != "cluster")
    .inform(sprintf("  Note: block resampling by '%s' with se='%s' (the replicates do not use cluster SEs).\n", clusterid, se))
  .check_choice("ci_type", ci_type, c("basic", "percentile", "studentized"))
  if (!(.is_whole(B))) stop(sprintf("B must be an integer; got %s.", .show_value(B)), call. = FALSE)
  if (B < 50) stop(sprintf("B should be at least 50 (got %d); 999 is recommended.", as.integer(B)), call. = FALSE)
  if (!(.is_whole(seed) && seed >= 0)) stop(sprintf("seed must be a non-negative integer; got %s.", .show_value(seed)), call. = FALSE)
  B <- as.integer(B)

  # ---- the estimation sample is what is resampled (item 16) ----
  prep <- .prepare_sample(data, Y, T, M, Z, X, dummies, fe, se, clusterid, weights, lat, lon)
  mask <- prep$mask
  if (method == "block") {
    if (!(clusterid %in% colnames(data))) stop(sprintf("column(s) not found in data: %s.", clusterid), call. = FALSE)
    mask <- mask & !is.na(.missing_text_to_na(data[[clusterid]]))
  }
  est_rows <- which(mask) - 1L
  data <- data[mask, , drop = FALSE]; rownames(data) <- NULL
  if (nrow(data) == 0L) stop("No observations remaining after dropping missing values.", call. = FALSE)

  # Replicates reuse the control lists resolved on the original data (A.1-9).
  iv_args <- list(Y = Y, T = T, M = M, Z = Z, X = prep$X, dummies = prep$dummies, fe = fe,
                  se = se, clusterid = if (se == "cluster") clusterid else NULL, alpha = alpha,
                  weights = weights, wtype = wtype, op_tau = op_tau, op_alpha = op_alpha, ar_ci = FALSE,
                  lat = lat, lon = lon, spatial_cutoff = spatial_cutoff,
                  spatial_kernel = spatial_kernel, spatial_metric = spatial_metric)
  if (verbose) cat("  Fitting original sample...\n")
  r_original <- do.call(ivmediate, c(list(data = data), iv_args))

  params <- c("total_hat", "tau_hat", "theta_hat", "gamma_hat", "indirect_hat", "kappa_star", "rho_TY_star")
  if (bootstrap_diagnostics) params <- c(params, "RV_theta", "XRV_theta")
  result_key <- c(theta_hat = "theta_hat", tau_hat = "tau_hat", gamma_hat = "gamma_hat",
                  total_hat = "total_hat", indirect_hat = "indirect", kappa_star = "kappa_star",
                  rho_TY_star = "rho_TY_star", RV_theta = "RV_theta", XRV_theta = "XRV_theta")
  point <- vapply(params, function(p) as.numeric(r_original[[result_key[[p]]]]), numeric(1))
  se_ind0 <- as.numeric(r_original$se_indirect)
  analytic_se <- list(theta_hat = r_original$se_theta, tau_hat = r_original$se_tau,
                      gamma_hat = r_original$se_gamma, total_hat = r_original$se_total,
                      indirect_hat = if (is.finite(se_ind0)) se_ind0 else NA_real_,
                      kappa_star = NA_real_, rho_TY_star = NA_real_)
  if (bootstrap_diagnostics) { analytic_se$RV_theta <- NA_real_; analytic_se$XRV_theta <- NA_real_ }
  analytic_ci <- list(theta_hat = r_original$ci_theta, tau_hat = r_original$ci_tau,
                      gamma_hat = r_original$ci_gamma, total_hat = r_original$ci_total,
                      indirect_hat = if (all(is.finite(r_original$ci_indirect))) r_original$ci_indirect else NULL,
                      kappa_star = NULL, rho_TY_star = NULL)
  need_se <- identical(ci_type, "studentized") || bootstrap_diagnostics

  n <- nrow(data)
  if (method == "block") {
    cl_codes <- .codes(data[[clusterid]])
    G <- max(cl_codes)
    cluster_rows <- split(seq_len(n), factor(cl_codes, levels = seq_len(G)))
    cluster_sizes <- vapply(cluster_rows, length, integer(1))
    nested_fe <- list()
    for (cc in prep$fe_cols) {
      fe_codes <- .codes(data[[cc]])
      if (nrow(unique(cbind(fe_codes, cl_codes))) == length(unique(fe_codes)))
        nested_fe[[cc]] <- list(codes = fe_codes, n_lev = max(fe_codes))
    }
  }

  # ---- random-number streams (D8): the caller's RNG state is restored ----
  if (is.null(boot_indices)) {
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
    on.exit({
      RNGkind(old_kind[1], old_kind[2], old_kind[3])
      if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
    stream <- get(".Random.seed", envir = globalenv())
  } else {
    boot_indices <- as.matrix(boot_indices)
    want <- if (method == "pairs") n else G
    if (nrow(boot_indices) < B || ncol(boot_indices) != want)
      stop(sprintf("boot_indices must have at least B = %d rows and %d columns (%s).", B, want,
                   if (method == "pairs") "rows of the estimation sample" else "clusters"), call. = FALSE)
  }
  idx_record <- if (.return_indices) vector("list", B) else NULL

  boot_dist <- matrix(NA_real_, nrow = B, ncol = length(params), dimnames = list(NULL, params))
  studentized <- identical(ci_type, "studentized")
  st_params <- c("theta_hat", "tau_hat", "gamma_hat", "total_hat", "indirect_hat")
  if (studentized) boot_se_replicate <- matrix(NA_real_, nrow = B, ncol = length(st_params), dimnames = list(NULL, st_params))
  n_failed <- 0L
  progress_step <- if (progress) max(B %/% 20L, 1L) else B + 1L
  if (progress) cat(sprintf("  Bootstrap (%s, B=%d)...\n", method, B))
  run_quiet <- function(dat, args) {
    sink(nullfile()); on.exit(sink(), add = TRUE)
    suppressWarnings(do.call(ivmediate, c(list(data = dat), args, list(.point_only = !need_se))))
  }

  for (b in seq_len(B)) {
    if (is.null(boot_indices)) {
      stream <- parallel::nextRNGStream(stream)
      assign(".Random.seed", stream, envir = globalenv())
      draw <- if (method == "pairs") sample.int(n, n, replace = TRUE) else sample.int(G, G, replace = TRUE)
    } else {
      draw <- as.integer(boot_indices[b, ]) + 1L
    }
    if (!is.null(idx_record)) idx_record[[b]] <- draw - 1L
    if (method == "pairs") {
      data_b <- data[draw, , drop = FALSE]; rownames(data_b) <- NULL
      args_b <- iv_args
    } else {
      concat_idx <- unlist(cluster_rows[draw], use.names = FALSE)
      data_b <- data[concat_idx, , drop = FALSE]; rownames(data_b) <- NULL
      copy_id <- rep(seq_len(G), times = cluster_sizes[draw])
      data_b[["_boot_cluster"]] <- copy_id
      for (cc in names(nested_fe))
        data_b[[cc]] <- (copy_id - 1L) * nested_fe[[cc]]$n_lev + nested_fe[[cc]]$codes[concat_idx]
      args_b <- iv_args
      if (se == "cluster") args_b$clusterid <- "_boot_cluster"
    }
    rb <- tryCatch(run_quiet(data_b, args_b), error = function(e) e)
    if (!inherits(rb, "error")) {
      for (p in params) boot_dist[b, p] <- as.numeric(rb[[result_key[[p]]]])
      if (studentized) {
        boot_se_replicate[b, "theta_hat"] <- rb$se_theta; boot_se_replicate[b, "tau_hat"] <- rb$se_tau
        boot_se_replicate[b, "gamma_hat"] <- rb$se_gamma; boot_se_replicate[b, "total_hat"] <- rb$se_total
        boot_se_replicate[b, "indirect_hat"] <- rb$se_indirect
      }
    } else {
      n_failed <- n_failed + 1L
      if (verbose) cat(sprintf("    Replicate %d: failed (%s)\n", b, conditionMessage(rb)))
    }
    if (progress && (b %% progress_step == 0L))
      cat(sprintf("    %5d/%d  (%4.1f%%, %d successful)\n", b, B, 100 * b / B, b - n_failed))
  }
  if (n_failed >= B)
    stop(sprintf("All %d bootstrap replicates failed. Check that the model specification is consistent with the data.", B), call. = FALSE)
  if (n_failed > 0.10 * B)
    .inform(sprintf("  Warning: %d/%d bootstrap replicates failed (%.1f%%). Bootstrap CIs may be unreliable.\n", n_failed, B, 100 * n_failed / B))

  boot_se <- stats::setNames(rep(NA_real_, length(params)), params)
  boot_ci <- stats::setNames(vector("list", length(params)), params)
  qq <- function(x) as.numeric(stats::quantile(x, c(alpha / 2, 1 - alpha / 2), type = 7, names = FALSE))
  for (p in params) {
    bd <- boot_dist[, p]; bv <- bd[!is.na(bd)]
    boot_se[p] <- if (length(bv) > 1L) stats::sd(bv) else NA_real_
    if (ci_type == "percentile") {
      boot_ci[[p]] <- qq(bv)
    } else if (ci_type == "basic") {
      q <- qq(bv); boot_ci[[p]] <- c(2 * point[[p]] - q[2], 2 * point[[p]] - q[1])
    } else if (p %in% st_params) {
      seb <- boot_se_replicate[, p]
      mk <- !is.na(bd) & !is.na(seb)
      if (sum(mk) < 10L) { boot_ci[[p]] <- c(NaN, NaN); next }
      t_stats <- (bd[mk] - point[[p]]) / seb[mk]
      t_stats <- t_stats[is.finite(t_stats)]
      if (length(t_stats) < 10L) { boot_ci[[p]] <- c(NaN, NaN); next }
      q <- qq(t_stats)
      point_se <- analytic_se[[p]]
      boot_ci[[p]] <- if (is.null(point_se) || !is.finite(point_se)) c(NaN, NaN)
                      else c(point[[p]] - q[2] * point_se, point[[p]] - q[1] * point_se)
    } else {
      q <- qq(bv); boot_ci[[p]] <- c(2 * point[[p]] - q[2], 2 * point[[p]] - q[1])
    }
  }
  out <- list(point = point, boot_se = boot_se, boot_ci = boot_ci, analytic_se = analytic_se,
              analytic_ci = analytic_ci, n_failed = n_failed, n_est = n, method = method, B = B,
              seed = seed, ci_type = ci_type, alpha = alpha)
  class(out) <- c("ivmediate_bootstrap", "list")
  if (return_boot_dist) out$boot_dist <- boot_dist
  if (.return_indices) { out$.indices <- do.call(rbind, idx_record); out$.rows <- est_rows }
  out
}
