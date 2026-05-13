---
name: ldsc
description: Run LD Score Regression (LDSC) for sumstats munging, SNP heritability (h2), and genetic correlation (rg). Auto-installs LDSC from the CBIIT Python 3 / Mac-compatible fork and caches the canonical EUR LD-score reference files on first use. Resolves method/software/paper provenance from the OKG when $OKG_REPO is set; otherwise runs without OKG provenance. Use when a user asks to munge GWAS sumstats, estimate heritability from sumstats, or compute genetic correlation between two GWAS traits.
license: MIT
compatibility: Requires Python 3.9+, git, and an internet connection on first run to clone CBIIT/ldsc and download the EUR LD-score reference (~250 MB extracted). Auto-installs the LDSC Python deps (numpy, pandas, scipy, bitarray) into the active environment via uv pip / pip --user.
metadata:
  author: stephendorn
  version: "0.1"
  source: agentskills.io/specification
  ldsc_fork: https://github.com/CBIIT/ldsc
  ldsc_upstream: https://github.com/bulik/ldsc
---

## What this skill does

Wraps the three core LDSC workflows behind a single CLI:

| Subcommand | Wraps | Purpose |
|---|---|---|
| `munge` | `munge_sumstats.py` | Standardize a raw GWAS TSV into LDSC's `.sumstats.gz` format. |
| `h2` | `ldsc.py --h2` | Estimate SNP heritability from one set of munged sumstats. |
| `rg` | `ldsc.py --rg` | Estimate genetic correlation between two sets of munged sumstats. |

On first use the skill:
1. Clones `https://github.com/CBIIT/ldsc` to `~/.cache/ldsc/repo/` (override via `--repo-cache`).
2. Installs LDSC's Python deps into the active env (`numpy`, `pandas`, `scipy`, `bitarray`).
3. Downloads `eur_w_ld_chr.tar.bz2` from the Broad and extracts to `~/.cache/ldsc/ld_scores/eur_w_ld_chr/` (override via `--ld-scores-dir`).

Subsequent runs reuse the cache; no network needed unless the user passes `--refresh`.

## Inputs you may need to elicit

Ask via `AskUserQuestion` when missing:

- **`munge`**: input TSV path, output prefix, and (if not standard) the sumstats column names (`--snp`, `--a1`, `--a2`, `--p`, `--n`, etc. mirror LDSC's `munge_sumstats.py` flags).
- **`h2`**: input `.sumstats.gz` and (optionally) a partitioned LD-score reference if the user wants stratified h2. Default reference is `eur_w_ld_chr/`.
- **`rg`**: two `.sumstats.gz` files (`--in1`, `--in2`) and the LD-score reference.

## How to execute

1. **Run the script**: `python3 scripts/ldsc.py <subcommand> <flags>`. The script handles install, cache, deps, sidecar manifest.

2. **Report back**:
   - Output paths (`.log`, `.sumstats.gz`, `.results`, etc. per subcommand)
   - Key numbers parsed from the LDSC `.log`:
     - munge: number of variants retained / dropped, mean chi-squared
     - h2: SNP h² estimate ± SE, intercept ± SE, ratio
     - rg: rg estimate ± SE, p-value, h² of each trait
   - Path to `<output>.ldsc.json` sidecar manifest with OKG provenance + reference-file SHA-256s + LDSC version (git commit).

3. **OKG provenance**:
   - If `$OKG_REPO` set: query `get_node("method:ldsc")` and `get_node("software:ldsc")` at the pinned generation; embed the node IDs in the sidecar manifest.
   - If `$OKG_REPO` unset: the skill still runs; manifest just omits `okg_node_ids`.

4. **Refusal triggers**:
   - munge: input has fewer than 100,000 variants after standardization (LDSC needs ≥200k for stable estimates; flag a warning at 200k, refuse below 100k).
   - h2: mean chi-squared < 1.02 (the GWAS is unpowered for h² estimation — flag, but don't refuse).
   - rg: either input's h² < 0 or SE > estimate (the input is too noisy for genetic correlation).

## Examples

```bash
# 1. Munge raw sumstats into LDSC format
python3 scripts/ldsc.py munge \
    --in GWAS/raw/T01.tsv.gz \
    --out GWAS/munged/T01 \
    --N 419013   # or auto-read from a column

# 2. Estimate heritability
python3 scripts/ldsc.py h2 \
    --in GWAS/munged/T01.sumstats.gz \
    --out h2/T01

# 3. Genetic correlation between two traits
python3 scripts/ldsc.py rg \
    --in1 GWAS/munged/T01.sumstats.gz \
    --in2 GWAS/munged/T02.sumstats.gz \
    --out rg/T01_T02

# 4. With OKG provenance (records method:ldsc, software:ldsc, paper:ldsc_2015 in the manifest)
OKG_REPO=~/Lab/KG/okg python3 scripts/ldsc.py h2 --in T01.sumstats.gz --out h2/T01
```

## Sidecar manifest schema

Every run writes `<out_prefix>.ldsc.json` with:

```json
{
  "subcommand": "h2",
  "input": "...", "output_prefix": "...",
  "ldsc_repo": "https://github.com/CBIIT/ldsc",
  "ldsc_commit": "<git sha>",
  "ld_scores_dir": "~/.cache/ldsc/ld_scores/eur_w_ld_chr",
  "ld_scores_sha256": "<sha256 of the tarball>",
  "key_results": {
    "h2": 0.123, "h2_se": 0.014,
    "intercept": 1.02, "intercept_se": 0.008,
    "ratio": 0.05, "mean_chi2": 1.34
  },
  "okg_node_ids": {"method": "method:ldsc", "software": "software:ldsc", "paper": "paper:ldsc_2015"},
  "captured_at": "2026-05-13T..."
}
```

## Notes and edge cases

- **CBIIT fork vs bulik/ldsc**: the OKG currently records `bulik/ldsc` as the canonical software. The CBIIT fork is functionally equivalent for the workflows here but Python-3 / macOS-arm64 compatible. Override via `--repo-url <url>` if you want the upstream or another fork.
- **LD-score files** are downloaded once from the Broad's public mirror. The default (EUR) covers EUR-ancestry GWAS; for trans-ancestry / non-EUR use, pass `--ld-scores-dir` to a panel matching your sumstats' ancestry. The skill flags ancestry-mismatch warnings when the input's `population` column doesn't match the LD-scores' assumed ancestry.
- **Allele alignment** during munge: LDSC handles allele flips automatically when an effect-allele frequency column is present. The skill prints how many SNPs were ambiguous and dropped.
- **Reproducibility**: pin LDSC by passing `--repo-commit <sha>` to clone-checkout a specific revision. Default = `main` HEAD.

See [references/LDSC_REFS.md](references/LDSC_REFS.md) for the canonical EUR/EAS/AFR LD-score-file URLs and citation guidance.
