---
name: finemap
description: Fine-mapping wrapper around the mancusolab/sushie Python package. Four modes: `susie` (single-ancestry, individual genotypes), `sushie` (multi-ancestry, individual genotypes), `region` (single-locus from GWAS sumstats + a reference VCF/PLINK; sushie computes LD internally), `sumstats` (precomputed Z + LD matrices, K=1 SuSiE-RSS or K≥2 SuShiE via `infer_sushie_ss`). Auto-installs sushie from GitHub on first use and runs the bundled 3-ancestry tutorial as a verification step. Records OKG provenance (method:sushie or method:susie_finemapping, software:sushie, paper:sushie_2025, plus ld_panel + dataset when supplied) in a per-run .finemap.json sidecar when $OKG_REPO is set. Use when a user asks to fine-map a locus with SuSiE or to do multi-ancestry fine-mapping with SuShiE.
license: MIT
compatibility: Requires Python 3.9+, git, and ~500 MB of disk for the sushie clone + tutorial data. Auto-installs sushie via `pip install .` (or `uv pip install`). OKG provenance optional, gated on $OKG_REPO.
metadata:
  author: stephendorn
  version: "0.1"
  source: agentskills.io/specification
  sushie_repo: https://github.com/mancusolab/sushie
  sushie_docs: https://mancusolab.github.io/sushie/
  sushie_paper_doi: 10.1038/s41588-025-02262-7
---

## What this skill does

Wraps the `sushie finemap` CLI behind a thin orchestrator that:

1. **Installs sushie** on first use (clones `mancusolab/sushie` to `~/.cache/sushie/repo`, runs `pip install .`).
2. **Verifies the install** by running the bundled 3-ancestry tutorial against `<repo>/data/`. Once-only; results recorded under `~/.cache/sushie/verify_install/`.
3. **Runs fine-mapping** in one of four modes:
   - **`susie`**: single-population, individual-level (one VCF + one phenotype file).
   - **`sushie`**: multi-population, individual-level (K VCFs + K phenotype files).
   - **`region`**: single-locus from GWAS sumstats + reference genotypes. Pass `--gwas-sumstats`, `--chrom/--start/--end`, `--N`, and one of `--ref-vcf` / `--ref-plink` / `--ref-bgen`; sushie computes the LD matrix internally. This is the simplest mode when you have a GWAS sumstats file and an ancestry-matched genotype reference panel (e.g. 1000G EUR).
   - **`sumstats`**: precomputed Z-score files + precomputed LD matrices per ancestry (K=1 SuSiE-RSS, K≥2 SuShiE) via `infer_sushie_ss` — use when you've already computed LD yourself.
4. **Records provenance** in a `.finemap.json` sidecar with sushie version (git SHA), bundled CLI flags, parsed CS / PIP summary, and OKG node IDs (method, software, paper, ld_panel, dataset where applicable).

### Choosing between `region` and `sumstats`

| You have… | Use `region` | Use `sumstats` |
|---|---|---|
| GWAS sumstats + reference VCF / PLINK / BGEN | ✓ (sushie computes LD) | requires precomputed LD |
| Z-scores + precomputed LD `.npy` per locus | n/a | ✓ |
| Multiple ancestries simultaneously (SuShiE) | not yet (use `sumstats`) | ✓ (pass K files per flag) |

For LD-source options + how to get a 1000G EUR (hg19/hg38) reference VCF, see [references/LD_OPTIONS.md](references/LD_OPTIONS.md).

## First-run setup (verify-install)

The first time the skill is invoked (or any time the user passes `--verify-install`), the skill runs:

```bash
cd ~/.cache/sushie/repo/data/
sushie finemap \
    --pheno EUR.pheno AFR.pheno EAS.pheno \
    --vcf vcf/EUR.vcf vcf/AFR.vcf vcf/EAS.vcf \
    --covar EUR.covar AFR.covar EAS.covar \
    --output ~/.cache/sushie/verify_install/test_result
```

This is the canonical tutorial from https://mancusolab.github.io/sushie/manual.html — 3 ancestries (489 EUR, 639 AFR, 481 EAS), 123/129/113 HapMap SNPs.

Success criteria:
- `test_result.cs.tsv` exists and contains at least one credible set
- `test_result.weight.tsv` exists
- No tracebacks in stderr

A `~/.cache/sushie/verify_install/.verified` sentinel file is written on success; future invocations skip the verification unless `--re-verify` is passed.

## Inputs you may need to elicit

For `susie` (single-ancestry):
- `--vcf <path>` — genotype VCF
- `--pheno <path>` — phenotype file (tab-separated, sushie format)
- `--covar <path>` (optional) — covariates
- `--out <prefix>` — output prefix

For `sushie` (multi-ancestry, K ≥ 2):
- `--vcf <eur.vcf> <afr.vcf> [...]` — K VCFs (order matters, must match pheno order)
- `--pheno <eur.pheno> <afr.pheno> [...]` — K phenotype files
- `--covar <eur.covar> <afr.covar> [...]` (optional)
- `--out <prefix>`

For `sumstats` (GWAS sumstats + LD; K=1 SuSiE-RSS, K≥2 SuShiE):
- `--z <z1.tsv> [<z2.tsv> ...]` — z-score TSV per ancestry (or BETA+SE columns; computed on the fly)
- `--ld <ld1.npy> [<ld2.npy> ...]` — LD matrix per ancestry (`.npy` or whitespace TSV, square, same SNP order as `--z`)
- `--n <n1> [<n2> ...]` — sample size per ancestry
- `--out <prefix>`
- See [references/SUMMARY_STATS.md](references/SUMMARY_STATS.md) for file formats and examples.

## How to execute

1. **Run the script**: `python3 scripts/finemap.py <subcommand> <flags>`.
2. **First-run install + verify** are automatic on the first invocation; subsequent runs skip them.
3. **Report back**:
   - Output prefix and files generated (`<out>.cs.tsv`, `<out>.weight.tsv`, `<out>.cv.tsv`, `<out>.log`)
   - Number of credible sets, mean CS size, top-PIP variant per CS
   - For `sushie`: the estimated cross-ancestry effect-size correlation matrix from `<out>.corr.tsv`
   - Path to `.finemap.json` sidecar with full provenance

## Examples

```bash
# Verify install (run the 3-ancestry tutorial once)
python3 scripts/finemap.py --verify-install

# Single-ancestry SuSiE on one locus
python3 scripts/finemap.py susie \
    --vcf GWAS/genotypes/MSH3.vcf \
    --pheno GWAS/phenos/T01.pheno \
    --out results/finemap/MSH3_T01

# Multi-ancestry SuShiE across EUR + AFR + EAS
python3 scripts/finemap.py sushie \
    --vcf GWAS/genotypes/EUR_MSH3.vcf GWAS/genotypes/AFR_MSH3.vcf GWAS/genotypes/EAS_MSH3.vcf \
    --pheno GWAS/phenos/EUR_T01.pheno GWAS/phenos/AFR_T01.pheno GWAS/phenos/EAS_T01.pheno \
    --out results/finemap/MSH3_T01_3anc

# Sumstats fine-mapping (single ancestry, SuSiE-RSS)
python3 scripts/finemap.py sumstats \
    --z  GWAS/sumstats/MSH3.z.tsv \
    --ld LD/MSH3.ld.npy \
    --n  361194 \
    --out results/finemap/MSH3_sumstats

# Sumstats fine-mapping (multi-ancestry SuShiE)
python3 scripts/finemap.py sumstats \
    --z  EUR_MSH3.z.tsv  AFR_MSH3.z.tsv  EAS_MSH3.z.tsv \
    --ld EUR_MSH3.ld.npy AFR_MSH3.ld.npy EAS_MSH3.ld.npy \
    --n  361194           48000           92000 \
    --out results/finemap/MSH3_3anc_sumstats

# Pin to a specific sushie commit (reproducibility)
python3 scripts/finemap.py sushie --sushie-commit <sha> --vcf ... --pheno ... --out ...
```

## Sidecar manifest schema

`<out>.finemap.json`:

```json
{
  "subcommand": "sushie",
  "output_prefix": "results/finemap/MSH3_T01_3anc",
  "sushie_repo": "https://github.com/mancusolab/sushie",
  "sushie_commit": "<git sha>",
  "n_ancestries": 3,
  "inputs": {
    "vcf":   ["EUR.vcf", "AFR.vcf", "EAS.vcf"],
    "pheno": ["EUR.pheno", "AFR.pheno", "EAS.pheno"],
    "covar": ["EUR.covar", "AFR.covar", "EAS.covar"]
  },
  "summary": {
    "n_cs": 3,
    "mean_cs_size": 5.7,
    "max_pip": 0.98,
    "cross_ancestry_correlation": [[1.0, 0.85, 0.74], [0.85, 1.0, 0.79], [0.74, 0.79, 1.0]]
  },
  "okg_node_ids": {
    "method":   "method:sushie",
    "software": "software:sushie",
    "paper":    "paper:sushie_2025"
  },
  "captured_at": "2026-05-13T..."
}
```

For `susie` mode, `okg_node_ids.method` is `method:susie_finemapping` (the canonical SuSiE method), with `software:sushie` still as the operational implementation — the same dual-provenance pattern used by the `ldsc` skill.

## Notes and edge cases

- **sushie's CLI is `sushie finemap`**; this skill wraps it with `python3 finemap.py {susie,sushie}` for the individual-level path, and calls `sushie.infer_ss.infer_sushie_ss` directly for the `sumstats` path.
- **Summary-statistics inputs**: pass per-ancestry z-score TSVs, LD matrices, and sample sizes via `--z`, `--ld`, `--n`. See [references/SUMMARY_STATS.md](references/SUMMARY_STATS.md) for file formats, an alignment guide, and notes on LD-panel pairing.
- **VCF format**: sushie accepts plink/VCF/BGEN. The skill's flags use `--vcf` for simplicity; pass `--plink <prefix>` to use plink format instead.
- **Single-ancestry mode** is sushie's degenerate K=1 case; results are equivalent to SuSiE-RSS at numerical tolerance.

For details on sushie's bundled test data and the verification step, see [references/EXAMPLES.md](references/EXAMPLES.md).
