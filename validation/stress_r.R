# R-only stress tests: error paths, input-style equivalence, invariants, edges.
#
#   Rscript stress_r.R   (from any directory; uses the data files next to this script)
#
# Loads the code of the R package (r/R/*.R); the expectations are those of version 2.0.0.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (.f in list.files(file.path(HERE, "..", "r", "R"), pattern = "[.]R$", full.names = TRUE))
  sys.source(.f, envir = globalenv())

D  <- utils::read.csv(file.path(HERE, "data_main.csv"), check.names = FALSE, stringsAsFactors = FALSE)
DW <- utils::read.csv(file.path(HERE, "data_weak.csv"), check.names = FALSE, stringsAsFactors = FALSE)

PASS <- 0L; FAIL <- 0L; FAILED <- character(0)
ok <- function(label, cond) {
  if (isTRUE(cond)) { PASS <<- PASS + 1L; cat(sprintf("  PASS  %s\n", label)) }
  else { FAIL <<- FAIL + 1L; FAILED <<- c(FAILED, label); cat(sprintf("  FAIL  %s\n", label)) }
}
# The package's notes and warnings go through message(); capture them along with the printed output.
capture_all <- function(expr) {
  msgs <- character()
  out <- withCallingHandlers(utils::capture.output(invisible(expr)),
                             message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  c(out, unlist(strsplit(paste(msgs, collapse = ""), "\n", fixed = TRUE)))
}
expect_error <- function(label, expr, pattern = NULL) {
  err <- tryCatch({ capture_all(force(expr)); NULL }, error = function(e) conditionMessage(e))
  if (is.null(err)) ok(paste0(label, " [expected error, got none]"), FALSE)
  else if (!is.null(pattern) && !grepl(pattern, err, fixed = TRUE))
    ok(sprintf("%s [error text mismatch: %s]", label, substr(err, 1, 100)), FALSE)
  else ok(label, TRUE)
}
run_quiet <- function(expr) { out <- capture_all(r <- force(expr)); attr(r, "console") <- out; r }
iv <- function(...) run_quiet(ivmediate(D, Y = "y", T = "tvar", M = "mvar", ...))
XB <- c("x1", "x2", "x3")
releq <- function(a, b, tol = 1e-8) {
  a <- as.numeric(a); b <- as.numeric(b)
  all(abs(a - b) <= tol * pmax(abs(a), abs(b), 1))
}

cat("== A. error paths ==\n")
expect_error("unknown se", iv(Z = "z1", X = XB, se = "bogus"), "se must be one of")                        # Phase B wording
expect_error("cluster without clusterid", iv(Z = "z1", X = XB, se = "cluster"), "requires clusterid")     # Phase B wording
expect_error("spatial missing lat/lon", iv(Z = "z1", X = XB, se = "spatial", spatial_cutoff = 100), "requires lat and lon")
expect_error("spatial missing cutoff", iv(Z = "z1", X = XB, se = "spatial", lat = "lat", lon = "lon"), "requires spatial_cutoff")
expect_error("spatial cutoff <= 0", iv(Z = "z1", X = XB, se = "spatial", lat = "lat", lon = "lon", spatial_cutoff = -5), "must be positive")
expect_error("bad spatial kernel", iv(Z = "z1", X = XB, se = "spatial", lat = "lat", lon = "lon", spatial_cutoff = 100, spatial_kernel = "gauss"), "uniform")
expect_error("bad spatial metric", iv(Z = "z1", X = XB, se = "spatial", lat = "lat", lon = "lon", spatial_cutoff = 100, spatial_metric = "manhattan"), "haversine")
expect_error("lat col missing", iv(Z = "z1", X = XB, se = "spatial", lat = "nope", lon = "lon", spatial_cutoff = 100), "not found")
D_badcoord <- D; D_badcoord$lat <- D_badcoord$lat * 10
expect_error("lat/lon out of range stops (fix 4)",
             run_quiet(ivmediate(D_badcoord, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB,
                                             se = "spatial", lat = "lat", lon = "lon", spatial_cutoff = 100)),
             "degree range")
expect_error("unknown wtype", iv(Z = "z1", X = XB, weights = "w", wtype = "fweight"), "wtype must be one of")   # Phase B wording
D_negw <- D; D_negw$w[5] <- -1
expect_error("negative weights",
             run_quiet(ivmediate(D_negw, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, weights = "w")),
             "non-negative")
D_zw <- D; D_zw$w <- 0
expect_error("all-zero weights",
             run_quiet(ivmediate(D_zw, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, weights = "w")),
             "All weights are zero")
expect_error("X var not found", iv(Z = "z1", X = "nonexistent_var"), "not found in data")
expect_error("glob matches 0 cols", iv(Z = "z1", X = "qq_*"), "matched 0 columns")
D_allna <- D; D_allna$x1 <- NA_real_
expect_error("all-NA column",
             run_quiet(ivmediate(D_allna, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB)),
             "only missing values")
expect_error("three-way FE", iv(Z = "z1", X = XB, fe = c("fe1", "fe2", "clu")), "two-way")
expect_error("lasso cv not ported", iv(Z = "z1", X = XB, lasso = TRUE, lasso_penalty = "cv"), "not available in the R port")
expect_error("lasso_keep unknown", iv(Z = "z1", X = XB, lasso = TRUE, lasso_keep = "zzz"), "not in the rich expansion")
expect_error("lasso_keep glob 0", iv(Z = "z1", X = XB, lasso = TRUE, lasso_keep = "zzz*"), "matched 0 names")
expect_error("unknown expand op", iv(Z = "z1", X = XB, lasso = TRUE, expand_x = c("quartic")), "Unknown expansion op")
D_dup <- D; D_dup[["x1^2"]] <- D_dup$x1^2
expect_error("duplicate rich names",
             run_quiet(ivmediate(D_dup, Y = "y", T = "tvar", M = "mvar", Z = "z1",
                                             X = c("x1", "x1^2"), lasso = TRUE)),
             "duplicate column names")
expect_error("boot B<50", ivmediate_bootstrap(D, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, B = 10), "at least 50")
expect_error("boot bad ci_type", ivmediate_bootstrap(D, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, ci_type = "bca"), "ci_type")
expect_error("boot block no clusterid", ivmediate_bootstrap(D, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, method = "block"), "requires clusterid")
expect_error("boot bad method", ivmediate_bootstrap(D, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, method = "wild"), "pairs")
ok("Phase D: the stability stubs are removed (Python-only, left out of 2.0.0)",
   !exists("lasso_selection_stability") && !exists("print_stability"))
expect_error("bad star cutoffs (2 values)", .validate_star_cutoffs(c(0.01, 0.05)), "3 numbers")
expect_error("bad star cutoffs (>1)", .validate_star_cutoffs(c(0.01, 0.05, 1.5)), "strictly between")
# All-NA Y is not caught by the X/dummies all-NA guard (same as Python);
# it falls through to the complete-cases drop -> "No observations remaining".
D_empty <- D; D_empty$y[] <- NA_real_
expect_error("no obs remaining",
             run_quiet(ivmediate(D_empty, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB)),
             "No observations remaining")

cat("== B. input-style equivalence ==\n")
r_vec   <- iv(Z = "z1", X = XB, se = "robust")
r_str   <- iv(Z = "z1", X = "x1 x2 x3", se = "robust")
r_multi <- iv(Z = "z1", X = "x1
                             x2   x3", se = "robust")
r_glob  <- iv(Z = "z1", X = c("x[12]", "x3"), se = "robust")   # glob inside vector
ok("X string == vector", releq(r_str$theta_hat, r_vec$theta_hat, 0) && releq(r_str$se_theta, r_vec$se_theta, 0))
ok("X multiline == vector", releq(r_multi$theta_hat, r_vec$theta_hat, 0))
ok("X glob-in-vector == vector", releq(r_glob$theta_hat, r_vec$theta_hat, 0) && r_glob$n_covariates == 3)
r_z_str <- iv(Z = "z1", X = XB)                                 # Z as bare string
ok("Z bare string works", releq(r_z_str$theta_hat, r_vec$theta_hat, 0))
r_x_null  <- iv(Z = "z1")
r_x_empty <- iv(Z = "z1", X = character(0))
r_x_blank <- iv(Z = "z1", X = "  ")
ok("X NULL == character(0) == blank string",
   releq(r_x_null$theta_hat, r_x_empty$theta_hat, 0) && releq(r_x_null$theta_hat, r_x_blank$theta_hat, 0) &&
   r_x_null$n_covariates == 0)
r_fe_str <- iv(Z = "z1", X = XB, fe = "fe1", se = "classical")
r_fe_vec <- iv(Z = "z1", X = XB, fe = c("fe1"), se = "classical")
ok("fe string == fe vector", releq(r_fe_str$theta_hat, r_fe_vec$theta_hat, 0))
r_dum_str <- iv(Z = "z1", X = XB, dummies = "dcat dnum")
r_dum_vec <- iv(Z = "z1", X = XB, dummies = c("dcat", "dnum"))
ok("dummies string == vector", releq(r_dum_str$theta_hat, r_dum_vec$theta_hat, 0) && r_dum_str$n_dummies == r_dum_vec$n_dummies)

cat("== C. invariants ==\n")
# C1: uniform weights == unweighted
D_u <- D; D_u$wu <- 2.0
r_unw <- iv(Z = "z1", X = XB, se = "robust")
r_uw  <- run_quiet(ivmediate(D_u, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, se = "robust",
                                         weights = "wu", wtype = "aweight"))
ok("uniform aweights == unweighted (theta, se, F, kappa*, RV)",
   releq(r_uw$theta_hat, r_unw$theta_hat, 1e-10) && releq(r_uw$se_theta, r_unw$se_theta, 1e-10) &&
   releq(r_uw$F_MZ_T, r_unw$F_MZ_T, 1e-10) && releq(r_uw$kappa_star, r_unw$kappa_star, 1e-10) &&
   releq(r_uw$RV_theta, r_unw$RV_theta, 1e-10))
# C2: one-way FE demeaning == explicit group dummies (drop singleton groups first)
tabs <- table(D$fe1); keepg <- names(tabs)[tabs >= 2]
D_nos <- D[D$fe1 %in% keepg, ]
r_fe  <- run_quiet(ivmediate(D_nos, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, fe = "fe1", se = "classical"))
r_dm  <- run_quiet(ivmediate(D_nos, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, dummies = "fe1", se = "classical"))
ok("FE demeaning == explicit dummies (theta)", releq(r_fe$theta_hat, r_dm$theta_hat, 1e-7))
ok("FE demeaning == explicit dummies (se_theta, classical)", releq(r_fe$se_theta, r_dm$se_theta, 1e-7))
ok("FE demeaning == explicit dummies (gamma, tau, total)",
   releq(r_fe$gamma_hat, r_dm$gamma_hat, 1e-7) && releq(r_fe$tau_hat, r_dm$tau_hat, 1e-7) &&
   releq(r_fe$total_hat, r_dm$total_hat, 1e-7))
# C3: K=1 decomposition identity is exact
r1 <- iv(Z = "z1", X = XB, se = "cluster", clusterid = "clu")
ok("K=1 decomposition exact", r1$decomp_error < 1e-10)
# C4: two-way cascade accounting matches Python (n=878, dropped=4)
r2w <- iv(Z = "z1", X = XB, fe = c("fe1", "fe2"))
ok("two-way FE cascade: 4 singletons dropped (fix 1)", r2w$n_dropped_singletons == 4 && r2w$n == 878)
r1w <- iv(Z = "z1", X = XB, fe = "fe1")
ok("one-way FE: exactly the direct singleton dropped", r1w$n_dropped_singletons == 1)
# C5: AR set consistent with Wald under strong ID (contains estimate; similar width)
ar_w <- r_vec$ar_theta$upper - r_vec$ar_theta$lower
wald_w <- r_vec$ci_theta[2] - r_vec$ci_theta[1]
ok("AR interval ~ Wald CI under strong ID",
   r_vec$ar_theta$shape == "interval" &&
   r_vec$ar_theta$lower < r_vec$theta_hat && r_vec$theta_hat < r_vec$ar_theta$upper &&
   ar_w / wald_w > 0.85 && ar_w / wald_w < 1.15)
# C6: ar_ci=FALSE -> NULL sets, print still works
r_noar <- iv(Z = "z1", X = XB, ar_ci = FALSE)
pr <- capture.output(print(r_noar, title = "no AR"))
ok("ar_ci=FALSE: NULL sets + printable", is.null(r_noar$ar_theta) && !any(grepl("ANDERSON", pr)))
# C7: alpha scales CI width by z-ratio
r05 <- iv(Z = "z1", X = XB, alpha = 0.05); r10 <- iv(Z = "z1", X = XB, alpha = 0.10)
ratio <- (r10$ci_theta[2] - r10$ci_theta[1]) / (r05$ci_theta[2] - r05$ci_theta[1])
ok("alpha=0.10 CI width ratio = z ratio", releq(ratio, qnorm(0.95) / qnorm(0.975), 1e-10))
# C8: weak data -> rays; contains_zero must be consistent with the AR test at 0
rw <- run_quiet(ivmediate(DW, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = "x1", se = "robust"))
ok("weak data: theta AR set is rays; zero-membership consistent with AR p-value",
   rw$ar_theta$shape == "rays" &&
   identical(isTRUE(rw$ar_theta$contains_zero), rw$ar_theta$pval_at_zero > 0.05))
prw <- capture.output(print(rw, title = "weak"))
ok("weak-IV print shows unbounded note", any(grepl("weak-IV signature", prw)) && any(grepl("unbounded or disconnected", prw)))

cat("== D. printing / cosmetics ==\n")
r_l <- iv(Z = "z1", X = XB, lasso = TRUE, show_selection = TRUE)
pl <- capture.output(print(r_l, title = "lasso print"))
ok("print shows LASSO header line (fix 5)", any(grepl("LASSO: post-double selection", pl)))
ok("print shows LASSO selection details", any(grepl("LASSO SELECTION DETAILS", pl)))
pr1 <- capture.output(print(r_vec, title = "t", stars = FALSE))
ok("stars=FALSE removes sig codes", !any(grepl("Sig. codes", pr1)))
pr2 <- capture.output(print(r_vec, title = "t", star_cutoffs = c(0.10, 0.01, 0.05)))
ok("star cutoffs auto-sorted", any(grepl("\\*\\*\\* p<0.01", pr2)))
# two-way demeaning non-convergence warning (fix 3): force with max_iter=1. A.2 (O-2): designs whose smaller FE
# dimension has <= 3,000 levels use the exact projection, so Irons-Tuck is forced with max_direct = 0 here.
dd <- data.frame(a = rnorm(50), g1 = rep(letters[1:5], 10), g2 = rep(LETTERS[1:10], each = 5))
wtxt <- capture_all(.demean_twoway(dd, "a", c("g1", "g2"), tol = 1e-16, max_iter = 1L, max_direct = 0L))
ok("two-way demeaning warns on non-convergence (fix 3; Irons-Tuck forced)", any(grepl("did not converge", wtxt)))
etxt <- capture_all(ex <- .demean_twoway(dd, "a", c("g1", "g2")))
ok("A.2: the default (exact) two-way projection equals the dummy-variable residuals, no warning",
   !any(grepl("did not converge", etxt)) && max(abs(ex$a - residuals(lm(a ~ g1 + g2, dd)))) < 1e-12)

cat("== E. LASSO edge cases ==\n")
r_le <- iv(Z = "z1", lasso = TRUE)               # no X, no dummies: empty candidate set
ok("LASSO with empty candidate set runs", r_le$n_covariates == 0 && r_le$lasso_info$n_total == 0)
r_lo <- iv(Z = "z1", X = XB, lasso = TRUE, lasso_lambda_override = 1e9)   # nothing survives
ok("LASSO huge lambda: selects nothing, still estimates", r_lo$lasso_info$n_selected_union == 0 && is.finite(r_lo$theta_hat))
r_lk <- iv(Z = "z1", X = XB, lasso = TRUE, lasso_keep = c("x1", "x2"), lasso_lambda_override = 1e9)
ok("LASSO keep survives huge lambda", r_lk$lasso_info$n_selected_union == 2 &&
     all(c("x1", "x2") %in% r_lk$lasso_info$union_names))
D_one <- D; D_one$dsingle <- "only"
r_d1 <- run_quiet(ivmediate(D_one, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, dummies = "dsingle"))
ok("single-category dummy skipped", r_d1$n_dummies == 0)

cat("== F. bootstrap (R RNG path) ==\n")
set.seed(1)
bt <- run_quiet(ivmediate_bootstrap(D, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB,
                                    se = "cluster", clusterid = "clu", B = 60, seed = 99,
                                    ci_type = "basic", progress = FALSE))
ok("boot block auto-selected", bt$method == "block")
ok("boot SEs positive & CIs ordered",
   all(is.finite(unlist(bt$boot_se[c("theta_hat", "gamma_hat", "total_hat")]))) &&
   bt$boot_ci$theta_hat[1] < bt$boot_ci$theta_hat[2] && bt$n_failed == 0)
ok("boot point == analytic point", releq(bt$point[["theta_hat"]], r1$theta_hat, 1e-12))
# Phase B (D3): the pairs note needs se = "cluster"; with other SEs the bootstrap does not use clusterid and
# prints nothing about it (as the Python reference).
nt <- capture_all(bt2 <- ivmediate_bootstrap(D, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB, se = "cluster",
                                                method = "pairs", clusterid = "clu", B = 60, seed = 99,
                                                ci_type = "percentile", progress = FALSE))
ok("pairs+clusterid note printed (se = 'cluster')", any(grepl("valid only under cross-sectional independence", nt)))
nt3 <- capture_all(bt3 <- ivmediate_bootstrap(D, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB,
                                                 method = "pairs", clusterid = "clu", B = 60, seed = 99,
                                                 ci_type = "percentile", progress = FALSE))
bt4 <- run_quiet(ivmediate_bootstrap(D, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = XB,
                                     method = "pairs", B = 60, seed = 99, ci_type = "percentile", progress = FALSE))
ok("pairs+clusterid with se = 'robust': no pairs note, same bootstrap as without clusterid",
   !any(grepl("valid only under cross-sectional independence", nt3)) && identical(bt3$boot_ci, bt4$boot_ci))
pb <- capture.output(print(bt, title = "boot"))
ok("print runs with stars", any(grepl("Sig. codes", pb)))

cat("== G. spatial extremes ==\n")
# Cutoff >> data extent: the uniform-kernel meat degenerates to ~0 (normal
# equations), diagonals can go negative -> NaN SEs. The run must COMPLETE
# (Python parity: NaN propagates through RV/XRV instead of crashing) and the
# PSD-failure warning must be surfaced whenever SEs did go NaN.
# Phase B (PY-13): a uniform cutoff that covers all pairs is an error (also with LASSO); a Bartlett kernel
# over all pairs runs with the > 50% neighbor-share warning (NaN SEs allowed, surfaced once).
expect_error("huge uniform cutoff (all pairs) is an error", iv(Z = "z1", X = XB, se = "spatial", lat = "lat", lon = "lon",
                                                                spatial_cutoff = 5000), "covers all pairs")
expect_error("huge uniform cutoff with LASSO is an error", iv(Z = "z1", X = XB, se = "spatial", lat = "lat", lon = "lon",
                                                               spatial_cutoff = 5000, lasso = TRUE), "covers all pairs")
r_sp_big <- iv(Z = "z1", X = XB, se = "spatial", lat = "lat", lon = "lon", spatial_cutoff = 5000, spatial_kernel = "bartlett")
sp_console <- attr(r_sp_big, "console")
ok("huge Bartlett cutoff completes with the > 50% warning (NaN SEs allowed, warning surfaced)",
   is.finite(r_sp_big$theta_hat) && any(grepl("> 50%", sp_console)) &&
   (is.finite(r_sp_big$se_theta) || any(grepl("variance was negative", sp_console))))
r_sp_tiny <- iv(Z = "z1", X = XB, se = "spatial", lat = "lat", lon = "lon", spatial_cutoff = 1e-6)
ok("tiny cutoff (self-only) ~ HC0-style robust", releq(r_sp_tiny$theta_hat, r_vec$theta_hat, 1e-10) &&
     abs(r_sp_tiny$se_theta / r_vec$se_theta - 1) < 0.05)

cat("== H. kitchen sink + moderate-n timing ==\n")
t0 <- proc.time()[3]
r_ks <- iv(Z = "z1", X = "x1 x2 x3 reg_*", dummies = c("dcat", "dnum"), fe = "fe1",
           se = "cluster", clusterid = "clu", weights = "w", wtype = "pweight", alpha = 0.10)
ok("kitchen sink (glob+dummies+FE+cluster+pweight+AR) runs",
   is.finite(r_ks$theta_hat) && is.finite(r_ks$se_theta) && !is.null(r_ks$ar_theta))
nbig <- 20000L
set.seed(42)
zb <- rnorm(nbig); eb <- matrix(rnorm(nbig * 3), ncol = 3)
tb <- 0.8 * zb + eb[, 1]
mb <- 0.6 * tb + 0.5 * eb[, 1] + sqrt(0.75) * eb[, 2]   # corr(e_T, e_M) -> strong F_{M,Z|T}
yb <- 0.3 * tb + 0.5 * mb + eb[, 3]
DB <- data.frame(y = yb, tvar = tb, mvar = mb, z1 = zb, x1 = rnorm(nbig),
                 clu = sample(sprintf("c%03d", 1:200), nbig, TRUE))
rb <- run_quiet(ivmediate(DB, Y = "y", T = "tvar", M = "mvar", Z = "z1", X = "x1",
                                      se = "cluster", clusterid = "clu"))
el <- proc.time()[3] - t0
ok(sprintf("n=20k cluster+AR runs fast (%.1fs total incl. kitchen sink)", el), el < 60)
ok("n=20k estimates sane", abs(rb$theta_hat - 0.5) < 0.05 && rb$ar_theta$shape == "interval")

cat(sprintf("\n==== STRESS SUMMARY: %d passed, %d failed ====\n", PASS, FAIL))
if (FAIL > 0) { cat("Failed:\n"); for (f in FAILED) cat("  -", f, "\n") }
