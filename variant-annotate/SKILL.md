---
name: variant-annotate
description: Per-rsID variant annotation against gnomAD (v4, GRCh38) and GTEx (v8). For a list of rsIDs, fetches the most-severe coding consequence + gene + per-population allele frequency from gnomAD, plus the top single-tissue eQTL and sQTL from GTEx. Designed for the post-fine-mapping "what is this SNP doing" lookup that researchers normally do by hand at gnomad.broadinstitute.org and gtexportal.org. Caches per-rsID JSON under ~/.cache/variant-annotate/ so repeated runs are free. Use when a user has a credible set, a GWAS top-hit list, or any rsID list and wants coding/regulatory context attached.
license: MIT
compatibility: Pure Python 3.9+, stdlib only (urllib + json). Network access to gnomad.broadinstitute.org and gtexportal.org required on cache miss.
metadata:
  author: stephendorn
  version: "0.1"
  source: agentskills.io/specification
  gnomad_api: https://gnomad.broadinstitute.org/api
  gtex_api: https://gtexportal.org/api/v2
---

## What this skill does

Given a list of rsIDs (CLI arg, file, or column of a TSV), enriches each
with two annotation blocks:

**gnomAD v4 (GRCh38)** — one POST against the public GraphQL endpoint:
- `variant_id` (chr-pos-ref-alt)
- most-severe canonical coding consequence (`missense_variant`,
  `synonymous_variant`, `intron_variant`, ...)
- gene symbol on the canonical transcript
- HGVS protein notation when coding
- per-population allele frequency (computed as `ac/an`) for nfe, afr,
  eas, sas, amr, fin, asj, mid, remaining

**GTEx v8** — two-step REST:
1. `/api/v2/dataset/variant?snpId=<rsid>` → resolves to a GTEx variantId
   like `chr19_44908684_T_C_b38`
2. `/api/v2/association/singleTissueEqtl?variantId=...` and
   `singleTissueSqtl?variantId=...` → all significant tissue/gene
   associations
- Reports the top eQTL and top sQTL by p-value (tissue, gene symbol,
  NES, p), plus the number of tissues with a reported association.

Outputs a single TSV with one row per input rsID and columns for both
annotation blocks. Empty cells when the rsID isn't found in either
source.

## Why a separate skill

This is the "I just got my credible set, now what?" lookup. It's the
same pattern whether the SNP list came from fine-mapping (sushie CS),
GWAS top hits (gwas-fetch), PRS scoring (prs), or replication
candidates — none of those workflows should own the gnomAD/GTEx
plumbing.

The downstream `finemap` skill calls into this one for its
`--annotate` table flag; it does not duplicate the helpers.

## Quick start

```bash
# Single rsID, both sources, print as a markdown table:
python scripts/annotate.py --rsid rs429358 --sources gnomad,gtex \
    --format markdown

# Many rsIDs from a file, gnomAD only, TSV output:
python scripts/annotate.py --rsid-file my_cs.rsids.txt \
    --sources gnomad --out annotated.tsv

# From a TSV column (joins back to your original rows):
python scripts/annotate.py --tsv my_cs.tsv --rsid-col rsID \
    --sources both --out my_cs.annotated.tsv
```

## CLI

```
python scripts/annotate.py
    [--rsid rs429358 ...           # one or more rsIDs inline
     | --rsid-file path.txt        # one rsID per line
     | --tsv path.tsv --rsid-col rsID]  # join onto an existing TSV
    --sources {gnomad|gtex|both}   # default: both
    --format {tsv|markdown|json}   # default: tsv
    --out <path>                   # default: stdout
    [--cache-dir ~/.cache/variant-annotate]
    [--no-cache]                   # bypass cached JSON, refetch
```

### Output columns

When `--sources both`:

| Column | Source | Description |
|---|---|---|
| `rsid` | input | The rsID queried |
| `gnomad_variant_id` | gnomAD | chr-pos-ref-alt (e.g. 19-44908684-T-C) |
| `consequence` | gnomAD | most-severe canonical consequence term |
| `gene` | gnomAD | gene symbol on the canonical transcript |
| `hgvsp` | gnomAD | HGVS protein change (or HGVS coding if non-protein) |
| `af_nfe`, `af_afr`, ... | gnomAD | per-population allele frequency |
| `gtex_variant_id` | GTEx | chr19_44908684_T_C_b38-style ID |
| `top_eqtl_tissue` | GTEx | tissue with smallest eQTL p-value |
| `top_eqtl_gene` | GTEx | eGene at that tissue |
| `top_eqtl_nes` | GTEx | normalized effect size |
| `top_eqtl_p` | GTEx | eQTL p-value |
| `top_sqtl_tissue` | GTEx | tissue with smallest sQTL p-value |
| `top_sqtl_gene` | GTEx | sGene at that tissue |
| `top_sqtl_p` | GTEx | sQTL p-value |
| `n_eqtl_tissues` | GTEx | # tissues with reported eQTL |
| `n_sqtl_tissues` | GTEx | # tissues with reported sQTL |

## Caching

Per-rsID JSON is written to:

- `~/.cache/variant-annotate/gnomad/<rsid>.json`
- `~/.cache/variant-annotate/gtex/<rsid>.json`

The cache is keyed strictly on rsID. Re-running the same list is free.
Pass `--no-cache` to force a refetch (e.g. after a gnomAD release).

## Limits, notes, and pitfalls

- **Build**: gnomAD r4 is GRCh38. GTEx v8 is GRCh38. The rsID is
  build-agnostic, so this skill works regardless of what build your
  upstream sumstats are in.
- **Multi-allelic rsIDs**: `gnomad.variant(rsid:...)` returns one
  variant. If an rsID has been split into multiple ALT entries, we
  get the first match — not all of them.
- **Rate**: queries are sequential by default to be nice to both APIs.
  For lists >100 rsIDs, expect ~30 s on a cold cache. We do not
  parallelize.
- **Missing tissues**: if GTEx reports no significant association for
  the variant, the `top_eqtl_*` and `top_sqtl_*` columns are empty;
  the variant's `gtex_variant_id` is still recorded so you know the
  rsID resolved.
- **NES sign**: GTEx NES is on the alt allele. Compare against your
  GWAS effect-allele convention before interpreting direction.

## Related

- `finemap` calls this skill for its `--annotate` flag on the `region`
  subcommand to enrich the credible-set markdown table.
- A future OKG-side change should ingest gnomAD/GTEx coverage as
  proper substrate (eqtl edges, sqtl edges, molQTL edges); this skill
  is intentionally read-only and stateless w.r.t. the graph.
