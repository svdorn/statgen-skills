# locuszoomr + LDlink reference

The `locuszoom` skill wraps the [`locuszoomr`](https://github.com/myles-lewis/locuszoomr) R package (Lewis & Knight). This file documents the function surface used by `scripts/locuszoom_run.R` and the LDlink REST endpoint it hits.

## locuszoomr functions (the four we use)

### `locus(...)` — build a locus object

```r
loc <- locus(
  data    = sumstats_df,      # data.frame with chrom, pos, p columns
  gene    = "APOE",           # OR `seqname = "19", xrange = c(start, end)`
  flank   = 1e5,              # window each side of gene; ignored if xrange given
  ens_db  = EnsDb.Hsapiens.v75,  # Ensembl DB for the right build
  chrom   = "chrom",          # column name for chromosome (auto-detected often)
  pos     = "pos",            # column name for base-pair position
  p       = "p",              # column name for p-value
  labs    = "rsid"            # column name for variant id (used in labels)
)
```

Returns an `S3` object with `$data` (sumstats subset), `$seqname`, `$xrange`, `$gene`, `$genes` (Ensembl gene records), and (after `link_LD`) `$data$ld` r² column.

**Common errors:**
- `gene = "XXX"` not in Ensembl DB → throws "gene not found".
- All variants on a different chromosome than the gene → empty `$data`.

### `link_LD(loc, token, pop, method)` — annotate r²

```r
loc <- link_LD(loc,
               token  = "<LDlink API token>",
               pop    = "EUR",            # 1000G population
               method = "matrix")         # "matrix" pulls all-pairs r²; default is single-pair
```

Hits the LDlink REST API at `https://ldlink.nih.gov/LDlinkRest/`. Adds a `ld` column to `loc$data` with r² values against the lead SNP. On HTTP 401/403, raises an error containing "401" or "Unauthorized" — the R script detects this and exits with code `3` (LDLINK_TOKEN_INVALID).

### `locus_plot(loc, ...)` — the canonical plot

```r
locus_plot(loc,
           labels    = c("index"),        # label the lead SNP
           legend_pos = "topleft",
           cex       = 0.8)
```

Renders a two-panel base-graphics plot: top = -log10(p) scatter coloured by LD r² to the lead SNP, bottom = Ensembl gene tracks within the window. The R script writes this twice — once to PDF via `pdf()` and once to PNG via `png()` at 300 dpi.

### `scatter_plot(loc)` + `genetracks(loc)` — fallback

If `locus_plot()` fails (some edge cases with sparse LD or missing gene records), the R script falls back to calling `scatter_plot()` and `genetracks()` separately in a 2-row `par(mfrow=c(2,1))` layout. Same visual content, fewer auto-decisions.

## LDlink API

Endpoint: `https://ldlink.nih.gov/LDlinkRest/ldproxy` (the matrix variant uses `/ldmatrix`).

- **Auth**: required token via `?token=<token>` query param. Get one at https://ldlink.nih.gov/?tab=apiaccess (free; takes seconds).
- **Populations**: `EUR` (default), `AFR`, `EAS`, `AMR`, `SAS`, plus 1000G sub-population codes like `CEU`, `YRI`, `CHB`, etc.
- **Rate limits**: generous for the matrix endpoint — a single regional plot is one request. Long-running batch use should add `Sys.sleep(0.5)` between loci.
- **Errors**: HTTP 401 (bad token), HTTP 400 (bad SNP / pop), HTTP 500 (transient). The R script exits with code 3 on auth failures and propagates other errors with a clear stderr line.

## Ensembl DB packages

The skill needs the right Ensembl DB for the input build:

| Build | Bioconductor package | Approx. download |
|---|---|---|
| hg19 / GRCh37 | `EnsDb.Hsapiens.v75` | ~150 MB |
| hg38 / GRCh38 | `EnsDb.Hsapiens.v86` | ~150 MB |

These install once via `BiocManager::install` and are cached in the system R library (typically `~/Library/R/<version>/library/` on macOS).

## Sumstats column normalisation

The Python orchestrator pre-normalises the input to a 6-column TSV:

```
rsid    chrom    pos    p    beta    se
```

It auto-detects from four common shapes:

| Input shape | rsid | chrom | pos | p | beta | se |
|---|---|---|---|---|---|---|
| GWAS Catalog harmonised (`*.h.tsv.gz`) | `hm_rsid` | `hm_chrom` or `chromosome` | `hm_pos` or `base_pair_location` | `p_value` | `hm_beta` or `beta` | `standard_error` |
| LDSC munged (`.sumstats.gz`) | `SNP` | (none — chr/pos derived from a separate map) | (none) | `P` | (Z-derived) | (Z-derived) |
| COJO / `.ma` | `SNP` | (none) | (none) | `p` | `b` | `se` |
| GWAS-SSF v1.0 | `rsid` | `chromosome` | `base_pair_location` | `p_value` | `beta` | `standard_error` |

For LDSC munged and COJO inputs that lack chrom/pos, the skill currently refuses if the user requested `--region` (it can't filter without coordinates). For `--gene` and `--lead-snp` modes the rsid is enough since `locus()` resolves coordinates from the Ensembl DB. But you'll generally get better results from the raw harmonised file — pass that when you have it.

## Output

- **`<prefix>.pdf`** — vector format, ~50–500 kB depending on SNP count. Good for figures in papers.
- **`<prefix>.png`** — 300 dpi raster, ~1–5 MB. Good for inline previews.
- **`<prefix>.locuszoom.json`** — sidecar with the locus selection, LD population, r² coverage stats, OKG node IDs, and locuszoomr version.

## Fine-mapping overlay (3-panel mode)

When `--finemap-pip <path>` or `--finemap-sidecar <path>` is supplied, the R worker:

1. Reads the PIP TSV with `read.table(header=TRUE, sep="\t")`.
2. Auto-detects columns:
   - **snp**: `snp`, `rsid`, `SNP`, `marker`, `variant`
   - **pip**: `sushie_pip_all`, `pip_all`, `pip`, `PIP`
   - **cs** (optional): `CSIndex`, `cs`, `credible_set`, `cs_index`
3. Joins to `loc$data` by rsid → fills `loc$data$pip` and `loc$data$cs`.
4. Renders three vertically-stacked panels via base-R `layout(matrix(1:3, ncol=1), heights=c(2,2,1.4))`:
   - Top: `locuszoomr::scatter_plot(loc)` — -log10(p) scatter, LD-coloured if `link_LD` succeeded.
   - Middle: a hand-built PIP scatter (`plot(pos, pip, ...)`) with credible-set colour mapping, top-PIP label, and reference lines at PIP=0.5 and PIP=0.95.
   - Bottom: `locuszoomr::genetracks(loc)` — Ensembl gene tracks.
5. Reports join coverage and credible-set summary to stderr, and writes them into the sidecar JSON via `LOCUSZOOM_R_SUMMARY`.

The credible-set palette uses an 8-colour distinguishable set (`#E41A1C`, `#377EB8`, `#4DAF4A`, `#984EA3`, `#FF7F00`, `#A65628`, `#F781BF`, `#999999`) and recycles modulo 8 for fine-mapping runs with >8 credible sets — rare in practice.

### finemap-skill sidecar handshake

Compatible with the `mancusolab/sushie`-based `finemap` skill's `.finemap.json` output. The Python orchestrator inspects these keys in order to locate the PIP TSV:

1. `cs_path` — preferred (has per-CS rows with `CSIndex`, `snp`, `pip_all`)
2. `weights_path` / `weight_path` — fallback (has per-SNP `sushie_pip_all` or `PIP`)
3. `pip_file` — generic key for future fine-mappers

And these OKG node IDs are propagated forward (when present) into the locuszoom sidecar under `okg_node_ids.finemap_*`:

- `okg_node_ids.method` → `finemap_method`
- `okg_node_ids.software` → `finemap_software`
- `okg_node_ids.paper` → `finemap_paper`
- `okg_node_ids.ld_panel` → `finemap_ld_panel`

This keeps the provenance chain intact: sumstats → fine-mapping (with its own LD panel and method) → plot (with LDlink-1000G LD for the visual).

## Caveats

- **LD population must match the GWAS ancestry.** Plotting an EUR GWAS with `--ld-pop EAS` makes the colour gradient meaningless.
- **The chrom column** must be unprefixed integers (`19`, not `chr19`). The Python normaliser strips `chr` prefixes.
- **`locus(gene=...)` is fuzzy** — if your gene symbol matches multiple Ensembl records, `locuszoomr` picks the first; pass `--region` explicitly to be deterministic.
- **`link_LD(method="matrix")`** loads all-pairs r², which is what we want for the colour gradient. Single-pair mode (`method = "default"`) only compares against the lead SNP; the plot looks similar but the underlying data is less rich.
