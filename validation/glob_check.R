suppressMessages(library(jsonlite))
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
# The code of the R package (r/R/*.R); .glob_to_regex is internal.
for (.f in list.files(file.path(HERE, "..", "r", "R"), pattern = "[.]R$", full.names = TRUE))
  sys.source(.f, envir = globalenv())
names <- c('x1','x2','x3','x22','reg_a','reg_b','y1','x[','a.b1','aXb1','x1^2','x1^23','x1 Q x2')
exp <- fromJSON(file.path(HERE, "fn_expected.json"))
nfail <- 0
for (p in names(exp)) {
  got <- sort(names[grepl(.glob_to_regex(p), names)])
  want <- sort(unlist(exp[[p]]))
  if (!identical(got, if (length(want)) want else character(0))) {
    nfail <- nfail + 1
    cat(sprintf("MISMATCH %-10s  R={%s}  PY={%s}\n", p, paste(got, collapse=","), paste(want, collapse=",")))
  }
}
cat(sprintf("glob battery: %d/%d patterns agree with fnmatch\n", length(exp) - nfail, length(exp)))
