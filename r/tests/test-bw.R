# Smoke test: Becker and Woessmann (2009) with the default options reproduces Table 6, column 1 of
# Ferrara, Gold, Heblich and Pinto (2026), rounded as printed. Base R only.
library(ivmediate)

bw <- read.csv(system.file("extdata", "becker_woessmann_2009.csv", package = "ivmediate"))
xs <- "f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop"
res <- ivmediate(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw", Z = "kmwittenberg", X = xs)

table6 <- list(                          # value, digits printed in the paper
  n           = list(res$n, 426, 0),
  total       = list(c(res$total_hat, res$se_total), c(0.586, 0.242), 3),
  tau         = list(c(res$tau_hat, res$se_tau), c(-0.299, 0.214), 3),
  theta       = list(c(res$theta_hat, res$se_theta), c(4.709, 2.164), 3),
  gamma       = list(c(res$gamma_hat, res$se_gamma), c(0.188, 0.027), 3),
  indirect    = list(c(res$indirect, res$se_indirect), c(0.886, 0.408), 3),
  ar_theta    = list(c(res$ar_theta$lower, res$ar_theta$upper), c(0.43, 10.20), 2),
  F_TZ        = list(res$F_TZ, 75.69, 2),
  F_MZ_T      = list(res$F_MZ_T, 16.36, 2),
  rho_TY_star = list(res$rho_TY_star, -0.257, 3))
for (nm in names(table6)) {
  v <- table6[[nm]]
  ok <- all(abs(round(as.numeric(v[[1]]), v[[3]]) - v[[2]]) < 0.5 * 10^-v[[3]] + 1e-12)
  cat(sprintf("%-12s %-30s Table 6: %s\n", nm, paste(format(as.numeric(v[[1]]), digits = 7), collapse = ", "),
              paste(v[[2]], collapse = ", ")))
  if (!ok) stop("Table 6, column 1 not reproduced: ", nm)
}

stopifnot(inherits(res, "ivmediate"), is.list(res), identical(res$se_type, "robust"),
          identical(res$version, as.character(utils::packageVersion("ivmediate"))),
          identical(res$F_MZ_T_above_crit, FALSE))

# The print method: default output, and the optional sections on request
p_default <- capture.output(print(res, title = "Becker and Woessmann (2009)"))
p_full <- capture.output(print(res, title = "Becker and Woessmann (2009)", kappa_table = TRUE, r2_sensitivity = TRUE))
stopifnot(any(grepl("Becker and Woessmann (2009)", p_default, fixed = TRUE)),
          !any(grepl("SENSITIVITY (reduced-form", p_default, fixed = TRUE)),
          any(grepl("SENSITIVITY (reduced-form", p_full, fixed = TRUE)),
          length(p_full) > length(p_default))

# A short bootstrap: class, print method, and the bootstrap passed to print()
boot <- ivmediate_bootstrap(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw", Z = "kmwittenberg", X = xs,
                            B = 50, progress = FALSE)
stopifnot(inherits(boot, "ivmediate_bootstrap"), boot$B == 50, boot$n_failed == 0, boot$n_est == 426,
          all(is.finite(unlist(boot$boot_se[c("theta_hat", "gamma_hat", "indirect_hat")]))))
invisible(capture.output(print(boot), print(res, boot = boot)))

print(res, title = "Becker and Woessmann (2009)")
cat("test-bw.R: Table 6, column 1 reproduced\n")
