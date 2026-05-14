# COJO sumstats format

COJO = **"Conditional and Joint"**, the GWAS sumstats convention invented for
[GCTA-COJO](https://yanglab.westlake.edu.cn/software/gcta/#COJO) by the
Yang lab. SBayesRC (and most other GCTB-family tools) consume sumstats in
this format.

A COJO file is a plain TSV with exactly these 8 columns, in this order:

| Col | Meaning | Notes |
|---|---|---|
| `SNP`  | variant identifier (rsid or `chr:pos:A1:A2`) | must match the LD reference's variant set |
| `A1`   | effect allele                                  | a.k.a. tested allele |
| `A2`   | other allele                                   | a.k.a. reference allele |
| `freq` | frequency of A1 in the GWAS / reference        | use `NA` if unknown — SBayesRC tolerates it |
| `b`    | signed effect size (log OR for binary traits)  | NOT z-score |
| `se`   | standard error of `b`                          | |
| `p`    | association p-value                            | numeric, not `1.0e-300` strings |
| `N`    | per-variant sample size                        | for case/control use `n_cases + n_controls` |

Header line uses those exact lower/mixed-case names. Tab-separated.

## How the `prs` skill produces this

The orchestrator `scripts/prs.py` auto-detects the input shape and converts
it to COJO before handing to SBayesRC. Currently supported inputs:

| Input shape | Detection heuristic |
|---|---|
| GWAS Catalog harmonised TSV (`*_h.tsv.gz`) | columns `hm_rsid`, `hm_effect_allele`, `hm_other_allele`, `hm_beta`, `p_value` |
| LDSC-munged `.sumstats.gz` | columns `SNP A1 A2 Z N` (signed effect derived as `b = Z * SE`; needs SE) |
| GWAS-SSF v1.0 bare-name | columns `rsid effect_allele other_allele beta standard_error p_value` |
| Pre-converted COJO TSV | already has the 8 columns; the converter is a no-op |

If a row is missing any required field (or its effect/p is NA), it's
dropped. The skill prints how many rows survive and refuses if fewer than
50,000 (SBayesRC needs the full HapMap3 scaffold).

## When you can skip the conversion

If you've already got a clean COJO TSV (e.g. from a previous SBayesRC run
or an upstream harmonization pipeline), the converter sees the 8 columns
and just rewrites them. No information lost. So passing
`--gwas-sumstats <file>.cojo.tsv` is fine and idempotent.
