#!/usr/bin/env Rscript
# locuszoom_run.R --- companion to locuszoom.py (the `locuszoom` skill orchestrator).
#
# Reads a normalised sumstats TSV (rsid, chrom, pos, p, beta, se), builds a
# locuszoomr locus object for the requested gene/region/lead-SNP window,
# annotates LD via LDlink, and writes <out>.pdf + <out>.png.
#
# Auto-installs locuszoomr, LDlinkR, ensembldb, EnsDb.Hsapiens.v75/v86,
# optparse from CRAN/Bioconductor on first use.
#
# Exit codes:
#   0  -> success
#   3  -> LDLINK_TOKEN_INVALID (HTTP 401/403 from LDlink)
#   4  -> BAD_LOCUS (gene not found / region empty)
#   5  -> SPARSE_LOCUS (< 50 SNPs in window without --force-sparse)

suppressPackageStartupMessages({
  ensure <- function(pkg, src = "cran") {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("[locuszoom] installing %s (%s)", pkg, src))
      if (src == "cran") {
        install.packages(pkg, repos = "https://cloud.r-project.org",
                          quiet = TRUE)
      } else if (src == "bioc") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org",
                            quiet = TRUE)
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
      }
    }
  }
  ensure("optparse", "cran")
  ensure("locuszoomr", "cran")
  ensure("LDlinkR", "cran")
  ensure("ensembldb", "bioc")
  library(optparse)
})

op <- OptionParser(option_list = list(
  make_option("--sumstats", type = "character"),
  make_option("--build", type = "character", default = "hg19"),
  make_option("--ld-pop", type = "character", default = "EUR"),
  make_option("--ldlink-token", type = "character"),
  make_option("--out", type = "character"),
  make_option("--gene", type = "character", default = NA_character_),
  make_option("--chrom", type = "character", default = NA_character_),
  make_option("--start", type = "integer", default = NA_integer_),
  make_option("--end", type = "integer", default = NA_integer_),
  make_option("--lead-snp", type = "character", default = NA_character_),
  make_option("--flank", type = "integer", default = 100000L),
  make_option("--force-sparse", type = "character", default = "0"),
  make_option("--finemap-pip", type = "character", default = NA_character_)
))
args <- parse_args(op)

stopifnot(file.exists(args$sumstats))
out_pdf <- paste0(args$out, ".pdf")
out_png <- paste0(args$out, ".png")
dir.create(dirname(args$out), showWarnings = FALSE, recursive = TRUE)

ensdb_pkg <- if (args$build == "hg38") "EnsDb.Hsapiens.v86" else "EnsDb.Hsapiens.v75"
if (!requireNamespace(ensdb_pkg, quietly = TRUE)) {
  message(sprintf("[locuszoom] installing %s from Bioconductor", ensdb_pkg))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org",
                      quiet = TRUE)
  }
  BiocManager::install(ensdb_pkg, ask = FALSE, update = FALSE, quiet = TRUE)
}
suppressPackageStartupMessages(library(locuszoomr))
suppressPackageStartupMessages(library(LDlinkR))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(ensdb_pkg, character.only = TRUE))

ss <- read.table(args$sumstats, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE,
                  na.strings = c("NA", "", "."))
# Coerce types
ss$pos  <- suppressWarnings(as.integer(ss$pos))
ss$p    <- suppressWarnings(as.numeric(ss$p))
ss$beta <- suppressWarnings(as.numeric(ss$beta))
ss$se   <- suppressWarnings(as.numeric(ss$se))
ss <- ss[!is.na(ss$pos) & !is.na(ss$p) & !is.na(ss$chrom) & ss$chrom != "", ]
ss$chrom <- sub("^chr", "", ss$chrom)

message(sprintf("[locuszoom] sumstats: %d rows after coercion", nrow(ss)))

# --- Build the locus object ---------------------------------------------
mode <- NA_character_
if (!is.na(args$gene)) {
  mode <- "gene"
} else if (!is.na(args$chrom)) {
  mode <- "region"
} else if (!is.na(args$`lead-snp`)) {
  mode <- "lead_snp"
}
if (is.na(mode)) {
  message("BAD_LOCUS: pass --gene, --chrom/--start/--end, or --lead-snp"); quit(status = 4)
}

ens_db <- get(ensdb_pkg)
loc <- tryCatch({
  if (mode == "gene") {
    locus(data = ss, gene = args$gene, flank = args$flank,
          ens_db = ens_db,
          chrom = "chrom", pos = "pos", p = "p", labs = "rsid")
  } else if (mode == "region") {
    locus(data = ss, seqname = args$chrom,
          xrange = c(args$start, args$end),
          ens_db = ens_db,
          chrom = "chrom", pos = "pos", p = "p", labs = "rsid")
  } else {
    # lead_snp: find row, derive a window
    snp_row <- ss[ss$rsid == args$`lead-snp`, , drop = FALSE]
    if (nrow(snp_row) == 0) {
      message(sprintf("BAD_LOCUS: lead SNP %s not in sumstats", args$`lead-snp`))
      quit(status = 4)
    }
    snp_chr <- as.character(snp_row$chrom[1])
    snp_pos <- as.integer(snp_row$pos[1])
    locus(data = ss, seqname = snp_chr,
          xrange = c(snp_pos - args$flank, snp_pos + args$flank),
          ens_db = ens_db,
          chrom = "chrom", pos = "pos", p = "p", labs = "rsid")
  }
}, error = function(e) {
  message(sprintf("BAD_LOCUS: locus() failed: %s", conditionMessage(e)))
  quit(status = 4)
})

n_window <- nrow(loc$data)
message(sprintf("[locuszoom] %d SNPs in window", n_window))
if (n_window < 50 && args$`force-sparse` != "1") {
  message(sprintf("SPARSE_LOCUS: only %d SNPs; pass --force-sparse to override",
                   n_window))
  quit(status = 5)
}

# Lead SNP (lowest p within window)
lead_idx <- which.min(loc$data$p)
lead_snp <- loc$data$rsid[lead_idx]
lead_p   <- loc$data$p[lead_idx]
lead_chr <- loc$data$chrom[lead_idx]
lead_pos <- loc$data$pos[lead_idx]
message(sprintf("[locuszoom] lead SNP %s @ %s:%s, p=%.2g",
                 lead_snp, lead_chr, lead_pos, lead_p))

# --- Fine-map PIP / CS overlay (optional) ---------------------------------
has_finemap <- !is.na(args$`finemap-pip`) && nzchar(args$`finemap-pip`) &&
  file.exists(args$`finemap-pip`)
finemap_summary <- list(n_variants_with_pip = NA_integer_,
                         max_pip = NA_real_, pip_lead_snp = NA_character_,
                         n_credible_sets = NA_integer_,
                         credible_set_sizes = NULL)
if (has_finemap) {
  fp <- args$`finemap-pip`
  message(sprintf("[locuszoom] reading finemap PIP file: %s", fp))
  pip_df <- tryCatch(
    read.table(fp, header = TRUE, sep = "\t",
                stringsAsFactors = FALSE, check.names = FALSE,
                na.strings = c("NA", "", ".")),
    error = function(e) {
      message(sprintf("[locuszoom] failed to read PIP file (%s); skipping overlay",
                       conditionMessage(e))); NULL
    })
  if (!is.null(pip_df) && nrow(pip_df) > 0) {
    cols_lc <- tolower(colnames(pip_df))
    pick <- function(...) {
      for (cand in tolower(c(...))) {
        i <- match(cand, cols_lc)
        if (!is.na(i)) return(colnames(pip_df)[i])
      }
      NA_character_
    }
    snp_col <- pick("snp", "rsid", "marker", "variant")
    pip_col <- pick("sushie_pip_all", "pip_all", "pip")
    cs_col  <- pick("csindex", "cs", "credible_set", "cs_index")
    if (is.na(snp_col) || is.na(pip_col)) {
      message(sprintf("[locuszoom] PIP file missing snp/pip column; cols = %s",
                       paste(colnames(pip_df), collapse = ", ")))
    } else {
      pip_df$.pip <- suppressWarnings(as.numeric(pip_df[[pip_col]]))
      pip_df$.cs  <- if (!is.na(cs_col)) as.character(pip_df[[cs_col]]) else NA_character_
      pip_df$.snp <- as.character(pip_df[[snp_col]])
      pip_df <- pip_df[!is.na(pip_df$.pip), ]
      m <- match(loc$data$rsid, pip_df$.snp)
      loc$data$pip <- pip_df$.pip[m]
      loc$data$cs  <- pip_df$.cs[m]
      n_pip <- sum(!is.na(loc$data$pip))
      message(sprintf("[locuszoom] joined %d/%d variants on rsid -> PIP",
                       n_pip, nrow(loc$data)))
      if (n_pip > 0) {
        pip_lead_i <- which.max(loc$data$pip)
        finemap_summary$n_variants_with_pip <- n_pip
        finemap_summary$max_pip <- max(loc$data$pip, na.rm = TRUE)
        finemap_summary$pip_lead_snp <- loc$data$rsid[pip_lead_i]
        cs_vec <- loc$data$cs[!is.na(loc$data$cs) & loc$data$cs != "" &
                                loc$data$cs != "NA"]
        if (length(cs_vec) > 0) {
          tbl <- table(cs_vec)
          finemap_summary$n_credible_sets <- length(tbl)
          finemap_summary$credible_set_sizes <- as.integer(tbl)
        } else {
          finemap_summary$n_credible_sets <- 0L
        }
      } else {
        has_finemap <- FALSE
        message("[locuszoom] no rsid overlap with sumstats; PIP overlay skipped")
      }
    }
  } else {
    has_finemap <- FALSE
  }
}

# --- LD annotation via LDlink ----------------------------------------------
n_ld_pairs <- NA_integer_
ld_status  <- "skipped"
loc_ld <- tryCatch({
  link_LD(loc, token = args$`ldlink-token`, pop = args$`ld-pop`,
          method = "matrix")
}, error = function(e) {
  msg <- conditionMessage(e)
  message(sprintf("[locuszoom] link_LD error: %s", msg))
  # If it's an auth error, return code 3
  if (grepl("401|403|invalid token|Unauthorized|authentication",
            msg, ignore.case = TRUE)) {
    message("LDLINK_TOKEN_INVALID")
    quit(status = 3)
  }
  NULL
})
if (!is.null(loc_ld)) {
  loc <- loc_ld
  if (!is.null(loc$data$ld)) {
    n_ld_pairs <- sum(!is.na(loc$data$ld))
    ld_status  <- "ok"
    message(sprintf("[locuszoom] %d SNPs annotated with LDlink r^2",
                     n_ld_pairs))
  } else {
    ld_status <- "no_ld_column"
  }
}

# --- Render plot to PDF + PNG ---------------------------------------------
# 2-panel (default): p-value scatter + gene tracks via locus_plot()
# 3-panel (--finemap-pip given): p-value scatter / PIP scatter / gene tracks
render_plot <- function() {
  if (has_finemap) {
    # Hand-built 3-panel layout via locuszoomr's atomic plot functions.
    layout(matrix(1:3, ncol = 1), heights = c(2, 2, 1.4))
    op <- par(mar = c(2.5, 4.2, 1.5, 1.0), no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    # Panel 1: -log10(p) scatter (LD-coloured if available)
    tryCatch(scatter_plot(loc, labels = c("index"), legend_pos = "topleft"),
             error = function(e) {
               message(sprintf("[locuszoom] p-scatter warn: %s",
                                conditionMessage(e)))
             })
    # Panel 2: PIP scatter, coloured by credible set if `cs` column present
    has_cs <- any(!is.na(loc$data$cs) & loc$data$cs != "" & loc$data$cs != "NA")
    pip_col <- rep("grey60", nrow(loc$data))
    if (has_cs) {
      cs_levels <- sort(unique(loc$data$cs[!is.na(loc$data$cs) &
                                              loc$data$cs != "" &
                                              loc$data$cs != "NA"]))
      palette_cs <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                       "#FF7F00", "#A65628", "#F781BF", "#999999")
      for (i in seq_along(cs_levels)) {
        m <- !is.na(loc$data$cs) & loc$data$cs == cs_levels[i]
        pip_col[m] <- palette_cs[((i - 1) %% length(palette_cs)) + 1]
      }
    }
    pip_pos <- loc$data$pos
    pip_y   <- loc$data$pip
    plot(pip_pos, pip_y, type = "n",
         xlab = "", ylab = "PIP", ylim = c(0, 1.02),
         xlim = c(loc$xrange[1], loc$xrange[2]), xaxt = "n", bty = "l")
    abline(h = c(0.5, 0.95), lty = 3, col = "grey80")
    # Plot NA-PIP variants first (background), then with-PIP on top.
    na_mask <- is.na(pip_y)
    points(pip_pos[na_mask], rep(0, sum(na_mask)),
           pch = 20, col = "grey85", cex = 0.5)
    points(pip_pos[!na_mask], pip_y[!na_mask],
           pch = 21, bg = pip_col[!na_mask], col = "black", cex = 1.1)
    # Annotate the top-PIP variant
    if (!is.na(finemap_summary$pip_lead_snp)) {
      pl <- which(loc$data$rsid == finemap_summary$pip_lead_snp)
      if (length(pl)) text(loc$data$pos[pl[1]], loc$data$pip[pl[1]],
                            labels = finemap_summary$pip_lead_snp,
                            pos = 3, cex = 0.8, offset = 0.3)
    }
    if (has_cs) {
      legend("topright",
             legend = paste0("CS ", cs_levels),
             pt.bg = palette_cs[seq_along(cs_levels)],
             pch = 21, bty = "n", cex = 0.8)
    }
    # Panel 3: gene tracks
    par(mar = c(3.5, 4.2, 0.5, 1.0))
    tryCatch(genetracks(loc),
             error = function(e) {
               message(sprintf("[locuszoom] genetracks warn: %s",
                                conditionMessage(e)))
             })
  } else {
    tryCatch(locus_plot(loc, labels = c("index"), legend_pos = "topleft"),
             error = function(e) {
               message(sprintf("[locuszoom] locus_plot warn: %s",
                                conditionMessage(e)))
               par(mfrow = c(2, 1)); scatter_plot(loc); genetracks(loc)
             })
  }
}

# Heights / widths grow when there's a PIP panel
plot_h <- if (has_finemap) 8 else 6
plot_w <- 8

pdf(out_pdf, width = plot_w, height = plot_h); render_plot(); dev.off()
message(sprintf("[locuszoom] wrote %s", out_pdf))

png(out_png, width = plot_w, height = plot_h, units = "in", res = 300)
render_plot(); dev.off()
message(sprintf("[locuszoom] wrote %s", out_png))

# --- Sidecar summary for the Python orchestrator --------------------------
summary_path <- Sys.getenv("LOCUSZOOM_R_SUMMARY", "")
if (nzchar(summary_path)) {
  summary <- list(
    locus = list(
      seqname = as.character(loc$seqname),
      xrange_start = as.integer(loc$xrange[1]),
      xrange_end   = as.integer(loc$xrange[2]),
      gene_hit = if (!is.na(args$gene)) args$gene else NA_character_
    ),
    n_snps_in_window = n_window,
    lead_snp = lead_snp,
    lead_p   = lead_p,
    lead_chr = as.character(lead_chr),
    lead_pos = as.integer(lead_pos),
    ld_status = ld_status,
    n_ld_pairs = n_ld_pairs,
    finemap = finemap_summary,
    locuszoomr_version = as.character(utils::packageVersion("locuszoomr"))
  )
  jsonlite_ok <- requireNamespace("jsonlite", quietly = TRUE)
  if (!jsonlite_ok) {
    install.packages("jsonlite", repos = "https://cloud.r-project.org",
                      quiet = TRUE)
  }
  writeLines(jsonlite::toJSON(summary, auto_unbox = TRUE, na = "null", pretty = TRUE),
              summary_path)
}

quit(status = 0)
