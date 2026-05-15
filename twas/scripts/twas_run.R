#!/usr/bin/env Rscript
# twas_run.R --- thin wrapper for FUSION.assoc_test.R + dep installer.
#
# Currently the Python orchestrator calls FUSION.assoc_test.R directly
# (no R logic of our own required); this stub exists to (1) install the
# R packages FUSION needs on first use, and (2) act as a future hook
# for post-processing FUSION output in R if we ever add it.
#
# Usage:
#   Rscript twas_run.R install-deps
#   Rscript twas_run.R --help

suppressPackageStartupMessages({
  ensure_cran <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("[twas] installing CRAN %s", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org",
                        quiet = TRUE)
    }
  }
  ensure_remotes <- function() {
    ensure_cran("remotes")
    library(remotes)
  }
  ensure_plink2R <- function() {
    if (!requireNamespace("plink2R", quietly = TRUE)) {
      message("[twas] installing plink2R from github (gabraham/plink2R)")
      ensure_remotes()
      remotes::install_github("gabraham/plink2R/plink2R", upgrade = "never",
                                quiet = TRUE)
    }
  }
})

args <- commandArgs(trailingOnly = TRUE)
op <- if (length(args)) args[[1]] else "install-deps"

if (op == "install-deps") {
  ensure_cran("optparse")
  ensure_cran("glmnet")
  ensure_cran("methods")
  ensure_cran("Rcpp")
  ensure_cran("RcppEigen")
  ensure_cran("here")     # FUSION.assoc_test.R uses here::here() for paths
  ensure_cran("magrittr") # commonly transitive but pin explicitly
  ensure_cran("data.table")
  ensure_cran("coloc")    # for --coloc_P option
  ensure_plink2R()
  message("[twas] dependencies installed successfully")
  quit(status = 0)
}

if (op == "--help" || op == "-h") {
  cat("Usage:\n",
      "  Rscript twas_run.R install-deps   # install FUSION R deps\n",
      "  Rscript twas_run.R --help\n",
      sep = "")
  quit(status = 0)
}

message("twas_run.R: unknown subcommand ", op)
quit(status = 1)
