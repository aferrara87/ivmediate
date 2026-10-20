# Becker and Woessmann (2009): IV mediation estimates (paper, Table 6, column 1).
#
# Install the package first (R CMD INSTALL r from the repository root, or
# remotes::install_github("aferrara87/ivmediate", subdir = "r")), then run from anywhere:
#   Rscript bw.R

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
if (requireNamespace("ivmediate", quietly = TRUE)) {
  library(ivmediate)
} else {  # running from a checkout without installing the package
  for (f in list.files(file.path(here, "..", "..", "r", "R"), pattern = "[.]R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}

bw <- utils::read.csv(file.path(here, "becker_woessmann_2009.csv"), check.names = FALSE)
res <- ivmediate(bw, Y = "inctaxpc", T = "f_prot", M = "f_rw", Z = "kmwittenberg",
                 X = "f_young f_jew f_fem f_ortsgeb f_pruss f_blind f_deaf f_dumb hhsize lnpop gpop",
                 se = "robust")
print(res, title = "Becker and Woessmann (2009)")
