# Summary-statistics fine-mapping (`sumstats` subcommand)

The `sumstats` subcommand calls `sushie.infer_ss.infer_sushie_ss` programmatically. It works for any K ≥ 1:

- **K=1** → SuSiE-RSS on z-scores + LD (single-ancestry).
- **K≥2** → SuShiE sumstats mode (multi-ancestry).

This is the path to use when you have GWAS summary statistics and an LD reference panel but no individual-level genotypes.

## Inputs

| Flag           | Meaning                                                                                      |
|----------------|----------------------------------------------------------------------------------------------|
| `--z FILE...`  | One z-score TSV per ancestry. SNPs must be in the same row order as the LD matrix.           |
| `--ld FILE...` | One LD matrix per ancestry. `.npy` (numpy save) or whitespace TSV. Must be symmetric.        |
| `--n INT...`   | Sample size (one int per ancestry).                                                          |
| `--out PREFIX` | Output prefix.                                                                               |
| `--L INT`      | Max single effects (default 10, the SuSiE/SuShiE paper default).                             |
| `--max-iter`, `--min-tol`, `--threshold`, `--purity`, `--min-snps` | sushie's standard inference knobs. |

### Z-score file format

Whitespace, tab, or comma-separated, **with a header**. The loader tries (in order):

1. A column named `Z`, `z`, `zscore`, or `z_score` → used directly.
2. Otherwise: columns `BETA`/`beta`/`b` **and** `SE`/`se`/`standard_error`/`stderr` → computes `z = beta / se`.

An optional `SNP` / `rsid` / `variant_id` / `SNPID` column is preserved into outputs but not required (row index is used as the fallback ID).

```text
SNP        BETA    SE     P
rs1001    -0.024  0.011  0.029
rs1002     0.103  0.014  2.3e-13
...
```

### LD matrix format

`.npy` is preferred (faster + lossless). Whitespace TSV also works.

Per-ancestry LD must be:
- square (m × m) with **the same m across ancestries**,
- in the **same SNP order** as the z-score file,
- a correlation matrix (sushie internally treats it as such).

Generate from PLINK with `plink --bfile <prefix> --r square --out <ld>` and save the lower-triangle output, or compute directly:

```python
import numpy as np
# X is samples x m centered+scaled genotype matrix
LD = np.corrcoef(X, rowvar=False)
np.save("locus.ld.npy", LD)
```

## Examples

### Single-ancestry (K=1, SuSiE-RSS)

```bash
python3 scripts/finemap.py sumstats \
    --z  GWAS/sumstats/MSH3.z.tsv \
    --ld LD/MSH3.ld.npy \
    --n  361194 \
    --out results/finemap/MSH3
```

### Multi-ancestry (K=3, SuShiE)

```bash
python3 scripts/finemap.py sumstats \
    --z  GWAS/sumstats/EUR_MSH3.z.tsv  GWAS/sumstats/AFR_MSH3.z.tsv  GWAS/sumstats/EAS_MSH3.z.tsv \
    --ld LD/EUR_MSH3.ld.npy            LD/AFR_MSH3.ld.npy            LD/EAS_MSH3.ld.npy \
    --n  361194                         48000                          92000 \
    --out results/finemap/MSH3_3anc
```

## Outputs

| File                     | Contents                                                              |
|--------------------------|-----------------------------------------------------------------------|
| `<out>.cs.tsv`           | Credible sets (from `result.cs`), plus a resolved `SNP` column.       |
| `<out>.weight.tsv`       | `SNPIndex`, `SNP`, `PIP` for every SNP (from `result.pip`).           |
| `<out>.finemap.json`     | Sidecar with sushie commit, parsed CS/PIP summary, OKG node IDs.      |

The sidecar `okg_node_ids` use the SuShiE node set when K ≥ 2 (`method:sushie`, `software:sushie`, `paper:sushie_2025`) and the SuSiE node set when K = 1 (`method:susie_finemapping`, `software:sushie`, `software:susie`, `paper:susie_2020`).

## Notes

- The sumstats path is more practical than `susie`/`sushie` for GWAS-scale data, since individual genotypes are not needed.
- The skill does **not** harmonize alleles or check that the LD reference panel matches the GWAS population. That responsibility is yours — mismatched LD will silently distort PIPs. For UKB-scale sumstats consider GCTB's `ldm13M` panel (`ld_panel:gctb_ukb_ldm13m` in the KG).
- For genome-build mismatches between the GWAS and the LD panel, use the [`liftover` skill](../../liftover/SKILL.md) on the sumstats first.
- For pulling raw sumstats from the GWAS Catalog, see the [`gwas-fetch` skill](../../gwas-fetch/SKILL.md).
