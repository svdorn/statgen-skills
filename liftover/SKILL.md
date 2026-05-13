---
name: liftover
description: KG-aware genomic liftover for GWAS summary statistics. Resolves source/target genome build (hg18/hg19/hg38) via the statgen-analysis OKG when available, falling back to the GWAS Catalog REST API. Refuses to lift unknown-provenance data and emits a coverage-gap proposal stub. Supports `--gcst <accession>` to fetch+lift directly from the GWAS Catalog. Use when a user asks to lift a sumstats file, asks whether a lift is needed for a GWAS/LD pair, or asks to fetch a GWAS Catalog study and align its build to an LD panel.
license: MIT
compatibility: Requires Python 3.9+ and uv (or pip). Installs pyliftover on first use. OKG resolution requires the user to set $OKG_REPO (or pass --okg-repo) to a local clone of an `okg` repo with a deployment exposing an MCP server at deployments/<deployment>/server.py (e.g. statgen-analysis). The skill refuses any --okg-* / --gcst invocation when $OKG_REPO is unset.
metadata:
  author: stephendorn
  version: "0.2"
  source: agentskills.io/specification
---

## What this skill does

Takes a GWAS sumstats TSV (or `.tsv.gz`) with chromosome + position columns and produces a build-lifted copy. Two entry points:

1. **Local file**: user supplies the input path; the skill resolves source/target builds (OKG-first, explicit fallback) and lifts.
2. **GWAS Catalog fetch**: user supplies `--gcst <accession>`; the skill resolves build via OKG, falls back to the GWAS Catalog REST API for build + FTP URL, downloads to a cache dir, then lifts.

The skill is **OKG-aware**: every successful lift cites the OKG node IDs it consulted. If the OKG can't resolve the dataset, the skill emits a paired-OpenSpec coverage-gap stub instead of guessing — the same provenance pattern used by the statgen-analysis benchmark harness.

## Inputs you may need to elicit

If the user's invocation is incomplete, ask via `AskUserQuestion`:

1. **Input** — either `--in <path>` to a local TSV/TSV.gz, or `--gcst <accession>` for a GWAS Catalog study (e.g. `GCST90704615`).
2. **Target build** — either `--target hg19|hg38`, or `--okg-ld-panel-id <ld_panel_node>` (the skill will resolve target from the panel's `attrs.genome_build`).
3. **Source build (optional)** — `--source hg19|hg38` to bypass OKG/API resolution, or `--okg-dataset-id <dataset_node>` to resolve via OKG.

## How to execute

1. **Run the script**: `python3 scripts/lift_sumstats.py <args>` — handles install of pyliftover on first use, OKG lookup, optional GWAS Catalog fetch, the actual lift, and sidecar manifest writing.

2. **Report back to the user**:
   - Output file path
   - Lift summary (`n_input → n_lifted`, drop rate)
   - If drop rate > 5%, flag as a warning (suggests source build is wrong)
   - Path to the sidecar `<output>.lift.json` for provenance review
   - If the skill refused: read the emitted coverage-gap stub at `okg-coverage-stubs/add-dataset-<slug>/proposal.md` and summarize what the user needs to add to the OKG.

3. **Choosing the source-resolution path**:
   - If the user gave `--source` explicitly: use it; skip OKG.
   - If the user gave `--okg-dataset-id`: query OKG; refuse with stub if missing.
   - If the user gave `--gcst`: try OKG alias-search for the GCST; if hit, use OKG's build; if miss, fall back to the GWAS Catalog REST API.
   - Otherwise: prompt the user.

4. **Refusal triggers** (the script handles these; surface them clearly):
   - OKG returned no node for the dataset AND no `--source` was supplied
   - Source and target builds disagree but pyliftover doesn't ship a chain for that pair (e.g. some hg17 paths)
   - Input has fewer than 1,000 variants (likely wrong file)
   - More than 20% of variants drop during lift (suggests source-build inference is wrong)

   In every refusal case, the script writes `<output>.lift.json` with `status: "refused"` and `reason: "<code>"`, plus a coverage-gap stub if appropriate.

## Examples

```bash
# Lift a local file to hg19, source build derived from OKG ld_panel
python3 scripts/lift_sumstats.py \
    --in GWAS/raw/T01.tsv.gz \
    --source hg38 \
    --okg-ld-panel-id ld_panel:gctb_ukb_ldm13m

# Fetch from GWAS Catalog (build inferred from REST API), lift to hg19
python3 scripts/lift_sumstats.py \
    --gcst GCST90704615 \
    --target hg19 \
    --cache-dir ~/.cache/gwas-catalog

# Explicit overrides (bypass OKG, bypass GWAS Catalog API)
python3 scripts/lift_sumstats.py \
    --in T01.tsv.gz --source hg38 --target hg19 --out T01.hg19.tsv
```

## Sidecar manifest schema

Every successful lift writes `<output>.lift.json` next to the output:

```json
{
  "status": "lifted",
  "input": "...", "output": "...",
  "source_build": "hg38", "target_build": "hg19",
  "chain_file_source": "pyliftover (UCSC)",
  "n_input": 6874, "n_lifted": 6840, "n_dropped": 34,
  "drop_rate": 0.0049,
  "dropped_reasons": {"unmapped": 30, "ambiguous": 3, "missing_field": 1},
  "chrom_col": "chromosome", "pos_col": "base_pair_location",
  "okg_node_ids": {"dataset": "...", "ld_panel": "..."},
  "captured_at": "2026-05-13T..."
}
```

Refusal manifests use `"status": "refused"` and a `"reason"` field.

For `--gcst` fetches, a `<cached>.fetch.json` is also written with the download URL, SHA-256, and the OKG-or-API provenance trail.

## Notes and edge cases

- **No conda needed** — pyliftover is a pure-Python wheel; `uv pip install pyliftover` or `pip install --user pyliftover` is enough.
- **Chain files** are downloaded by pyliftover from UCSC on first use; subsequent runs are offline.
- **Column-name overrides** — pass `--chrom-col <name>` and `--pos-col <name>` if the input uses different headers (default: `chromosome`, `base_pair_location` — the GWAS Catalog harmonized convention).
- **`source == target`** — writes a build-suffixed copy with a no-op manifest so file naming stays uniform across a pipeline.
- **`gzip` input** — auto-detected via `.gz` extension; output is plain `.tsv` (re-gzip downstream if needed).

For the full coverage-gap workflow and how the skill integrates with the paired-OpenSpec pattern used in the statgen-analysis OKG, see [references/COVERAGE_GAPS.md](references/COVERAGE_GAPS.md).
