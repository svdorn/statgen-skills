#!/usr/bin/env Rscript
# sbayesrc_run.R --- companion to prs.py (the `prs` skill orchestrator).
#
# Reads a COJO-format sumstats file + an SBayesRC LD-eigendecomposition
# folder, runs the canonical SBayesRC 3-step pipeline, and writes
# per-SNP posterior effect sizes (PRS weights) to <out>.snpRes.
#
# Usage (called from prs.py; not intended for direct invocation):
#   Rscript sbayesrc_run.R \
#     --cojo <file.cojo.tsv> \
#     --ld-dir <unzipped LD folder> \
#     --annot-cache <path to baseline2.2 annotation file> \
#     --out <output prefix>

suppressPackageStartupMessages({
  ok <- requireNamespace("optparse", quietly = TRUE)
  if (!ok) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  library(optparse)
})

op <- OptionParser(option_list = list(
  make_option("--cojo", type = "character"),
  make_option("--ld-dir", type = "character"),
  make_option("--annot-cache", type = "character"),
  make_option("--out", type = "character")
))
args <- parse_args(op)

stopifnot(file.exists(args$cojo))
stopifnot(dir.exists(args$`ld-dir`))

log_path <- paste0(args$out, ".sbayesrc.log")
log_fp <- file(log_path, open = "w")
log_line <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), " ",
                paste(..., collapse = " "))
  writeLines(msg, log_fp); flush(log_fp)
  message(msg)
}

log_line("prs/sbayesrc_run.R starting")
log_line("cojo =", args$cojo)
log_line("ld-dir =", args$`ld-dir`)
log_line("out =", args$out)

# ----- Install SBayesRC if needed -----------------------------------------
if (!requireNamespace("SBayesRC", quietly = TRUE)) {
  log_line("SBayesRC R package not installed; installing from GitHub...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("zhilizheng/SBayesRC", upgrade = "never")
}
suppressPackageStartupMessages(library(SBayesRC))
log_line("SBayesRC version:",
          as.character(utils::packageVersion("SBayesRC")))

# ----- Cache the BaselineLD annotation file -------------------------------
# The zip at the GCTB host unpacks to a single ~2 GB TSV. SBayesRC's
# `sbayesrc(annot = ...)` parameter wants that TSV path directly.
annot_path <- args$`annot-cache`
annot_url <- "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/annot_baseline2.2.zip"
if (!file.exists(annot_path) || file.size(annot_path) < 1e8) {
  log_line("downloading BaselineLD 2.2 annotation zip ->", annot_path)
  dir.create(dirname(annot_path), showWarnings = FALSE, recursive = TRUE)
  zip_path <- paste0(annot_path, ".zip")
  utils::download.file(annot_url, destfile = zip_path, mode = "wb")
  utils::unzip(zip_path, exdir = dirname(annot_path))
}
log_line("annotation file:", annot_path)

# ----- Tidy / harmonise sumstats to the LD reference's SNP scaffold -------
out_prefix <- args$out
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

# SBayesRC's tidy/impute write to <output> directly (no .ma suffix);
# subsequent steps consume the same path that tidy wrote.
tidy_path <- paste0(out_prefix, ".tidy")
log_line("step 1/3: SBayesRC::tidy()")
SBayesRC::tidy(
  mafile = args$cojo,
  LDdir = args$`ld-dir`,
  output = tidy_path,
  log2file = TRUE
)

# ----- Impute any missing SNPs from the LD reference ----------------------
imputed_path <- paste0(out_prefix, ".imputed")
log_line("step 2/3: SBayesRC::impute()")
SBayesRC::impute(
  mafile = tidy_path,
  LDdir = args$`ld-dir`,
  output = imputed_path,
  log2file = TRUE
)

# ----- Run the SBayesRC sampler -------------------------------------------
log_line("step 3/3: SBayesRC::sbayesrc()")
fit <- SBayesRC::sbayesrc(
  mafile = imputed_path,
  LDdir = args$`ld-dir`,
  annot = annot_path,
  outPrefix = out_prefix,
  log2file = TRUE
)

# ----- Surface key fit stats to the log for the orchestrator to parse -----
if (is.list(fit)) {
  if (!is.null(fit$hsq))   log_line(sprintf("hsq = %.6g", fit$hsq))
  if (!is.null(fit$hsqSE)) log_line(sprintf("hsq_se = %.6g", fit$hsqSE))
  if (!is.null(fit$Pi))    log_line(sprintf("Pi = %s",
                                              paste(signif(fit$Pi, 4),
                                                    collapse = ",")))
  if (!is.null(fit$niter)) log_line(sprintf("niter = %d", fit$niter))
  if (!is.null(fit$nburn)) log_line(sprintf("nburn = %d", fit$nburn))
}

log_line("done; weights at", paste0(out_prefix, ".snpRes"))
close(log_fp)
