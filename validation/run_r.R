# Run the R implementation across the config grid; flatten results to JSON.
#
#   Rscript run_r.R            (from any directory; writes out_r.json next to this file)
#
# Mirrors run_py.py. configs.json holds three lists, each run against data_<data>.csv:
#   configs        analysis configs (ivmediate + print)
#   boot_configs   bootstrap configs (ivmediate_bootstrap + print), replaying the
#                  Python draws in bootidx_<name>.csv (B x n_est for pairs, B x G for block)
#   error_configs  degenerate inputs that must stop with an error of class "expect"
# A config may name Y, T and M (default y, tvar, mvar); every other field except name, data
# and expect is passed to the function. Output: one flat object of strings per config, with
# metadata _status, _error_class, _ERROR and _stdout_lines (see README.md).
suppressMessages(library(jsonlite))

HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

# The code of the R package (r/R/*.R), loaded into the global environment (the harness also
# calls internal helpers such as .op_critical_value).
for (f in list.files(file.path(HERE, "..", "r", "R"), pattern = "[.]R$", full.names = TRUE))
  sys.source(f, envir = globalenv())

CFG <- fromJSON(file.path(HERE, "configs.json"),
                simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)

# read.csv parses every float exactly; empty text fields stay "" (missing for the port).
read_data <- function(nm) utils::read.csv(file.path(HERE, sprintf("data_%s.csv", nm)),
                                          check.names = FALSE, stringsAsFactors = FALSE)
DATASETS <- lapply(c(main = "main", weak = "weak", ext = "ext", negvar = "negvar", chain = "chain"), read_data)
ROLE_DEFAULTS <- list(Y = "y", T = "tvar", M = "mvar")
ERROR_CLASSES <- CFG$error_classes

error_class <- function(msg) {
  low <- tolower(msg)
  for (ec in ERROR_CLASSES) {
    phrases <- unlist(ec[[2]])
    if (any(vapply(phrases, function(ph) grepl(ph, low, fixed = TRUE), logical(1)))) return(ec[[1]])
  }
  "other"
}

fnum <- function(v) {
  v <- as.numeric(v)
  if (is.nan(v)) return("NaN")
  if (is.na(v)) return("NULL")          # R's NA plays Python's None
  sprintf("%.17g", v)
}

put <- function(out, path, v) {
  s <- if (is.null(v) || length(v) == 0) "NULL"
  else if (is.logical(v)) { if (is.na(v)) "NULL" else if (v) "TRUE" else "FALSE" }
  else if (is.numeric(v)) { if (is.integer(v)) { if (is.na(v)) "NULL" else as.character(v) } else fnum(v) }
  else as.character(v)
  assign(path, s, envir = out)
}

# Scalars as `path`, vectors (K >= 2) as `path.1`, `path.2`, ...
put_vec <- function(out, path, v) {
  if (length(v) > 1L) for (i in seq_along(v)) put(out, sprintf("%s.%d", path, i), v[i])
  else put(out, path, v)
}

flat_ar <- function(out, pre, ar) {
  if (is.null(ar)) { put(out, pre, "NULL"); return(invisible(NULL)) }
  for (k in c("shape", "method", "contains_zero", "cv", "df_chi2", "alpha",
              "stat_at_zero", "pval_at_zero", "lower", "upper",
              "left_upper", "right_lower")) put(out, paste0(pre, ".", k), ar[[k]])
  if (!is.null(ar$segments)) {
    put(out, paste0(pre, ".n_segments"), length(ar$segments))
    for (i in seq_along(ar$segments)) {
      s <- ar$segments[[i]]
      put(out, sprintf("%s.seg%d.lo", pre, i), s$lo)
      put(out, sprintf("%s.seg%d.hi", pre, i), s$hi)
      put(out, sprintf("%s.seg%d.open_left", pre, i), isTRUE(s$open_left))
      put(out, sprintf("%s.seg%d.open_right", pre, i), isTRUE(s$open_right))
    }
  }
}

flat_eps <- function(out, pre, e) {
  e <- as.numeric(e)
  put(out, paste0(pre, ".mean"), mean(e))
  put(out, paste0(pre, ".sd"), sqrt(mean((e - mean(e))^2)))
  for (i in 1:3) put(out, sprintf("%s.v%d", pre, i), e[i])
}

flat_lasso_detail <- function(out, pre, d) {
  if (is.null(d)) { put(out, pre, "NULL"); return(invisible(NULL)) }
  put(out, paste0(pre, ".n_selected"), as.integer(d$n_selected))
  put(out, paste0(pre, ".lambda_value"), d$lambda_value)
  put(out, paste0(pre, ".converged"), isTRUE(d$converged))
  put(out, paste0(pre, ".n_iter"), as.integer(d$n_iter))
  put(out, paste0(pre, ".method"), d$method)
  put(out, paste0(pre, ".names"), paste(d$names, collapse = "|"))
  put(out, paste0(pre, ".coefs"), paste(sprintf("%.10g", d$coefs), collapse = "|"))
  put(out, paste0(pre, ".indices"), paste(as.integer(d$indices), collapse = "|"))   # 1-based
}

flatten_result <- function(r) {
  out <- new.env()
  for (k in c("n", "n_instruments", "n_covariates", "n_dummies", "se_type", "alpha",
              "n_fe_groups", "n_dropped_singletons", "df_fe_absorbed", "weighted", "wtype",
              "spatial_cutoff", "spatial_kernel", "spatial_metric",
              "theta_hat", "se_theta", "tau_hat", "se_tau", "gamma_hat", "se_gamma",
              "total_hat", "se_total", "indirect", "se_indirect", "cov_gamma_theta",
              "decomp_check", "decomp_error",
              "F_TZ", "F_MZ_T", "F_TZ_op", "F_MZ_T_op",
              "F_TZ_op_critval", "F_MZ_T_op_critval", "op_tau", "op_alpha",
              "c_TM_hat", "rho_TM_hat", "sigma_eT", "sigma_eY",
              "kappa_star", "kappa_star_CI", "rho_TY_star",
              "f_reduced_form", "f_first_stage_cond",
              "RV_theta", "XRV_theta", "rv_rf", "rv_fs", "xrv_rf", "xrv_fs",
              "lasso", "lasso_penalty",
              # 2026-10 (WP8)
              "n_clusters", "K_eff_TZ", "K_eff_MZ_T",
              "kappa_star_CI_informative", "F_MZ_T_above_crit")) put(out, k, r[[k]])
  put_vec(out, "pi1", r$pi1)
  put_vec(out, "se_pi1", r$se_pi1)
  nested <- r$fe_nested_in_cluster
  put(out, "fe_nested_in_cluster", if (is.null(nested)) "NULL" else paste(ifelse(nested, "TRUE", "FALSE"), collapse = ","))
  put(out, "fe", if (is.null(r$fe)) "NULL" else paste(r$fe, collapse = ","))
  for (nm in c("ci_theta", "ci_tau", "ci_gamma", "ci_total", "ci_indirect")) {
    put(out, paste0(nm, ".1"), r[[nm]][1]); put(out, paste0(nm, ".2"), r[[nm]][2])
  }
  for (i in seq_along(r$F_MZ_T_individual))
    put(out, sprintf("F_MZ_T_individual.%d", i), r$F_MZ_T_individual[i])
  for (i in seq_along(r$kappa_adj_table)) {
    row <- r$kappa_adj_table[[i]]
    for (k in c("kappa", "kappa_frac", "theta_adj", "ci_lo", "ci_hi"))
      put(out, sprintf("kappa_adj.%d.%s", i, k), row[[k]])
  }
  for (nm in c("ar_theta", "ar_gamma", "ar_total")) flat_ar(out, nm, r[[nm]])
  for (nm in c("eps_T_hat", "eps_M_hat", "eps_Y_hat")) flat_eps(out, nm, r[[nm]])
  li <- r$lasso_info
  if (!is.null(li)) {
    put(out, "lasso_info.n_selected_union", as.integer(li$n_selected_union))
    put(out, "lasso_info.n_total", as.integer(li$n_total))
    put(out, "lasso_info.union_names", paste(li$union_names, collapse = "|"))
    det <- li$details
    flat_lasso_detail(out, "lasso_info.Y", det$Y)
    flat_lasso_detail(out, "lasso_info.T", det$T)
    flat_lasso_detail(out, "lasso_info.M", det$M)
    if (!is.null(det$Z)) for (k in seq_along(det$Z))
      flat_lasso_detail(out, sprintf("lasso_info.Z%d", k), det$Z[[k]])
  }
  as.list(out)
}

flatten_boot <- function(bt) {
  out <- new.env()
  for (k in c("n_failed", "method", "B", "ci_type", "alpha", "n_est")) put(out, k, bt[[k]])
  for (p in names(bt$point)) put(out, paste0("point.", p), bt$point[[p]])
  for (p in names(bt$boot_se)) put(out, paste0("boot_se.", p), bt$boot_se[[p]])
  for (p in names(bt$boot_ci)) {
    put(out, paste0("boot_ci.", p, ".1"), bt$boot_ci[[p]][1])
    put(out, paste0("boot_ci.", p, ".2"), bt$boot_ci[[p]][2])
  }
  for (p in names(bt$analytic_se)) put(out, paste0("analytic_se.", p), bt$analytic_se[[p]])
  if (!is.null(bt$boot_dist)) {
    for (p in colnames(bt$boot_dist)) {
      vec <- bt$boot_dist[, p]; ok <- vec[!is.na(vec)]
      put(out, sprintf("boot_dist.%s.mean", p), mean(ok))
      put(out, sprintf("boot_dist.%s.sd", p), stats::sd(ok))
      for (i in 1:3) put(out, sprintf("boot_dist.%s.v%d", p, i), vec[i])
    }
  }
  as.list(out)
}

split_cfg <- function(cfg) {
  name <- cfg$name; data <- DATASETS[[cfg$data]]; expect <- cfg$expect
  roles <- lapply(names(ROLE_DEFAULTS), function(k) if (!is.null(cfg[[k]])) cfg[[k]] else ROLE_DEFAULTS[[k]])
  names(roles) <- names(ROLE_DEFAULTS)
  kw <- cfg[setdiff(names(cfg), c("name", "data", "expect", names(ROLE_DEFAULTS)))]
  list(name = name, data = data, roles = roles, kw = kw, expect = expect)
}

error_entry <- function(msg) list(`_status` = "error", `_error_class` = error_class(msg), `_ERROR` = paste0("Error: ", msg))

ALL <- list()
n_unexpected <- 0L
report <- function(name, out, expect) {
  got <- if (identical(out[["_status"]], "error")) out[["_error_class"]] else "ok"
  want <- if (is.null(expect)) "ok" else expect
  msg <- if (got == "ok") "ok" else sprintf("ERROR [%s] %s", got, out[["_ERROR"]])
  if (got != want) { n_unexpected <<- n_unexpected + 1L; msg <- paste0(msg, sprintf("   <-- UNEXPECTED (expected %s)", want)) }
  cat(sprintf("[R] %s: %s\n", name, msg))
}

# ---- static: version, OP table, computed critical values ----
st <- new.env()
put(st, "version", IVMED_VERSION)
for (K in 1:30) for (tau in c("0.05", "0.10", "0.20", "0.30"))
  put(st, sprintf("op.%d.%s", K, tau), .OP_CRIT_5PCT_WALD[as.character(K), tau])
cases <- list(c(1, 0.10, 0.05), c(5, 0.30, 0.05), c(30, 0.05, 0.05),
              c(31, 0.10, 0.05), c(1, 0.15, 0.05), c(1, 0.10, 0.01))
for (cc in cases)
  put(st, sprintf("opcv.%d.%s.%s", as.integer(cc[1]), format(cc[2]), format(cc[3])),
      .op_critical_value(cc[1], cc[2], cc[3]))
ALL[["_static"]] <- as.list(st)

# The package's notes and warnings go through message(); capture them along with the printed output.
capture_all <- function(expr) {
  msgs <- character()
  out <- withCallingHandlers(utils::capture.output(invisible(expr)),
                             message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  c(out, unlist(strsplit(paste(msgs, collapse = ""), "\n", fixed = TRUE)))
}

# ---- analysis configs and degenerate inputs ----
for (cfg in c(CFG$configs, CFG$error_configs)) {
  sc <- split_cfg(cfg)
  res <- tryCatch({
    txt_out <- capture_all({
      r <- do.call(ivmediate, c(list(data = sc$data), sc$roles, sc$kw))
      print(r, title = sc$name, kappa_table = TRUE, r2_sensitivity = TRUE)   # smoke test, optional sections on (A.2)
    })
    fl <- flatten_result(r)
    fl[["_status"]] <- "ok"
    fl[["_stdout_lines"]] <- as.character(length(txt_out))
    fl
  }, error = function(e) error_entry(conditionMessage(e)))
  report(sc$name, res, sc$expect)
  ALL[[sc$name]] <- res
}

# ---- bootstrap configs (replay the Python draws on the same estimation sample) ----
for (bcfg in CFG$boot_configs) {
  sc <- split_cfg(bcfg)
  idx <- as.matrix(utils::read.csv(file.path(HERE, sprintf("bootidx_%s.csv", sc$name)), header = FALSE))
  rows_py <- scan(file.path(HERE, sprintf("bootidx_%s_rows.csv", sc$name)), quiet = TRUE)
  res <- tryCatch({
    txt_out <- capture_all({
      bt <- do.call(ivmediate_bootstrap,
                    c(list(data = sc$data), sc$roles,
                      list(return_boot_dist = TRUE, progress = FALSE, boot_indices = idx, .return_indices = TRUE),
                      sc$kw))
      print(bt, title = sc$name)
    })
    if (!identical(as.numeric(bt$.rows), as.numeric(rows_py)))
      cat(sprintf("[R] %s: WARNING estimation sample differs from bootidx_%s_rows.csv (%d vs %d rows)\n",
                  sc$name, sc$name, length(bt$.rows), length(rows_py)))
    fl <- flatten_boot(bt)
    fl[["_status"]] <- "ok"
    fl
  }, error = function(e) error_entry(conditionMessage(e)))
  report(sc$name, res, sc$expect)
  ALL[[sc$name]] <- res
}

writeLines(toJSON(ALL, auto_unbox = TRUE, null = "null"), file.path(HERE, "out_r.json"))
cat(sprintf("wrote out_r.json (%d configs; %d with an unexpected status)\n", length(ALL) - 1L, n_unexpected))
quit(status = if (n_unexpected > 0L) 1L else 0L)
