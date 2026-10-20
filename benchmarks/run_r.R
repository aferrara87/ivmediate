# Speed comparison, R (benchmarks/README.md).
#
# Times ivmediate() on data/bench_<N>.csv with the defaults and with the
# post-double-selection LASSO: one untimed warm-up, then 3 timed runs (only the
# estimation call; loading the data and printing are excluded). Writes
# results/r.csv and results/env_r.txt.
#
# Usage:  Rscript run_r.R            (all sizes)
#         Rscript run_r.R 1000       (selected sizes)

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
if (requireNamespace("ivmediate", quietly = TRUE)) {
  library(ivmediate)
} else {  # a checkout without the installed package
  for (f in list.files(file.path(here, "..", "r", "R"), pattern = "[.]R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}

user_sizes <- as.integer(commandArgs(trailingOnly = TRUE))
sizes <- if (length(user_sizes)) user_sizes else c(1000L, 10000L, 100000L)
X <- paste0("x", 1:50, collapse = " ")
specs <- list(default = FALSE, lasso = TRUE)
reps <- 3L

fit <- function(df, lasso) {
  out <- NULL
  invisible(capture.output(suppressMessages(
    out <- ivmediate(df, Y = "y", T = "t", M = "m", Z = "z", X = X, fe = "g", lasso = lasso))))
  out
}

rows <- list()
for (n in sizes) {
  df <- utils::read.csv(file.path(here, "data", sprintf("bench_%d.csv", n)))
  for (spec in names(specs)) {
    fit(df, specs[[spec]])                                   # warm-up
    times <- numeric(reps)
    for (i in seq_len(reps)) {
      t0 <- proc.time()[["elapsed"]]
      r <- fit(df, specs[[spec]])
      times[i] <- proc.time()[["elapsed"]] - t0
    }
    rows[[length(rows) + 1L]] <- data.frame(
      language = "R", n = n, spec = spec, t1 = times[1], t2 = times[2], t3 = times[3],
      median = stats::median(times), theta = r$theta_hat, gamma = r$gamma_hat,
      tau = r$tau_hat, total = r$total_hat, se_type = r$se_type)
    cat(sprintf("R      n=%6d %-7s median %8.3f s\n", n, spec, stats::median(times)))
  }
}

dir.create(file.path(here, "results"), showWarnings = FALSE)
suffix <- if (length(user_sizes)) paste0("_", paste(sizes, collapse = "_")) else ""
utils::write.csv(do.call(rbind, rows), file.path(here, "results", paste0("r", suffix, ".csv")),
                 row.names = FALSE)
blas <- basename(extSoftVersion()[["BLAS"]])
writeLines(sprintf("%s; ivmediate %s; BLAS %s; LAPACK %s", R.version.string, r$version,
                   if (nzchar(blas)) blas else "(R's reference BLAS)", La_version()),
           file.path(here, "results", "env_r.txt"))
