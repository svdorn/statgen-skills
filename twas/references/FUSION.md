# TWAS-FUSION reference

The `twas` skill wraps **TWAS-FUSION** (Gusev et al. 2016 *Nat Genet*, PMID 26854917, doi 10.1038/ng.3506; repo: <https://github.com/gusevlab/fusion_twas>; project page: <https://gusevlab.org/projects/fusion/>).

## What FUSION does

FUSION imputes the genetic component of gene expression from cis-SNP weights pre-computed in a reference panel (e.g., GTEx tissues), then tests the imputed expression against GWAS summary statistics for trait association — producing a per-gene TWAS statistic.

## Pipeline (the `assoc` subcommand)

1. **`FUSION.assoc_test.R`** — the primary script. Per chromosome:
   - Reads the sumstats (rsid + A1 + A2 + Z).
   - Loads the panel's `*.pos` manifest, which lists every gene's weight file (`*.wgt.RDat`) and physical position.
   - For each gene: aligns GWAS SNPs to the prediction-model SNPs, imputes GWAS Z-scores for any reference SNPs missing in the GWAS (IMPG algorithm), computes the TWAS statistic.
   - Writes a `.dat` table with one row per gene tested.

The `twas` skill loops over every chromosome present in the panel's `.pos`, concatenates the per-chromosome `.dat` files into one `<out>.assoc.tsv`, and parses the result for the sidecar.

## Sumstats format expected by FUSION

Whitespace-delimited (FUSION accepts space or tab) with at minimum:

| Column | Meaning |
|---|---|
| `SNP` | variant identifier (rsid) |
| `A1`  | first allele (effect allele) |
| `A2`  | second allele (other allele) |
| `Z`   | Z-score, sign with respect to `A1` |

Additional columns are allowed and ignored. The `twas` skill auto-detects from harmonised GWAS Catalog, LDSC munged `.sumstats.gz`, COJO `.ma`, or generic; if `Z` is missing, it computes `Z = beta/SE` (or `Z = log(OR)/SE`).

## `FUSION.assoc_test.R` output schema

The per-chromosome `.dat` table columns (in order):

| Column | Meaning |
|---|---|
| `FILE` | Path to the gene's weight `.wgt.RDat` |
| `ID` | Gene identifier (ENSG or HGNC symbol depending on panel) |
| `CHR` | Chromosome |
| `P0` | Gene start position |
| `P1` | Gene end position |
| `HSQ` | cis-heritability estimate for this gene's expression |
| `BEST.GWAS.ID` | rsID of the most significant GWAS SNP in the locus |
| `BEST.GWAS.Z` | Z-score of that top GWAS SNP |
| `EQTL.ID` | rsID of the best eQTL in the locus (per the FUSION model) |
| `EQTL.R2` | eQTL r² for the best eQTL |
| `EQTL.Z` | Z-score of the best eQTL on GWAS |
| `EQTL.GWAS.Z` | GWAS Z-score for the best eQTL SNP |
| `NSNP` | Number of SNPs used in the model |
| `MODEL` | Best-performing prediction model (BLUP / LASSO / TopSNP / Enet) |
| `MODELCV.R2` | Cross-validated R² of the prediction model |
| `MODELCV.PV` | P-value of the cross-validated R² |
| **`TWAS.Z`** | **TWAS Z-score (primary statistic of interest)** |
| **`TWAS.P`** | **TWAS P-value** |

The `twas` skill summarises this into `n_genes_tested`, `n_significant_bonferroni` (at 0.05 / n_genes), and the top 10 hits (by TWAS.P).

## GTEx v8 weight panels

Each FUSION GTEx v8 panel is one `.tar.gz` containing:

- One `<Tissue>.pos` file (filtered — significantly-heritable genes only; recommended)
- One `<Tissue>.no_filter.pos` file (all genes)
- A `WEIGHTS/` subdir with one `<gene>.wgt.RDat` per gene

Pass `--no-filter` to switch from the filtered `.pos` (default) to `.no_filter.pos`.

URL pattern (49 tissues × 2 ancestries = 98 panels):

```
https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/{ALL,EUR}/GTExv8.{ALL,EUR}.<Tissue>.tar.gz
```

Examples:

```
GTExv8.EUR.Brain_Cortex.tar.gz
GTExv8.ALL.Whole_Blood.tar.gz
GTExv8.EUR.Adipose_Subcutaneous.tar.gz
```

The OKG registers all 98 panels as `dataset:fusion_gtex_v8_{eur,all}:<tissue_slug>` nodes with `source_url` pointing at the S3 URL.

## LD reference

FUSION expects the canonical 1000G EUR LDREF (hg19, PLINK1 bfile per chromosome) at `--ref_ld_chr ./LDREF/1000G.EUR.`:

```
https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
```

The `twas` skill downloads + caches this at `~/.cache/twas/ldref/LDREF/` on first use. OKG node: `ld_panel:fusion_1000g_eur`.

## R dependencies

Installed automatically by `scripts/twas_run.R install-deps` on first use:

- `optparse`, `glmnet`, `methods`, `Rcpp`, `RcppEigen` (CRAN)
- `plink2R` (via `remotes::install_github("gabraham/plink2R/plink2R")`)

## Sidecar manifest

`<out>.twas.json` cites:

- `okg_node_ids.method` = `method:fusion_twas`
- `okg_node_ids.software` = `software:fusion`
- `okg_node_ids.paper` = `paper:fusion_2016`
- `okg_node_ids.ld_panel` = `ld_panel:fusion_1000g_eur`
- `okg_node_ids.tissue` = `tissue:<slug>` (the GTEx v8 tissue used)
- `okg_node_ids.cohort` = `cohort:gtex_v8`
- `okg_node_ids.dataset_panel` = `dataset:fusion_gtex_v8_{eur,all}:<slug>`
- `okg_node_ids.dataset_gwas` = the user-supplied GWAS dataset (optional)

Plus reproducibility metadata: FUSION repo commit, panel URL + SHA-256, LDREF URL.

## Caveats

- **`MODELCV.PV` significance** doesn't mean the gene has a useful prediction model. Some genes pass cis-heritability filters but have low cross-validated R²; their TWAS Z is then noisy. Filter by both `MODELCV.PV < 0.05` AND `TWAS.P < 0.05/n_genes` for confident hits.
- **MHC region** (chr6:25–34 Mb) often dominates immune-trait TWAS scans. Inspect carefully or use FUSION's `--no_MHC` flag (forwarded via `extra`).
- **TWAS associations don't prove colocalization** — use `--coloc-p` (calls FUSION's COLOC integration) or follow up with the `finemap` + `locuszoom` skills to assess whether the TWAS signal is driven by a single shared causal variant.
