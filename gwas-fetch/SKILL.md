---
name: gwas-fetch
description: Fetch harmonised GWAS summary statistics from the GWAS Catalog by accession (GCST*). Resolves the study's metadata (genome build, FTP URL, sample size, ancestry) via the OKG when available, falling back to the GWAS Catalog REST API. Caches downloads to a configurable directory and writes a sidecar provenance manifest with SHA-256 + KG node IDs. Use when a user asks to download a GWAS Catalog study, fetch sumstats for a GCST accession, or pull sumstats for an analysis pipeline. Does NOT lift coordinates — chain into the `liftover` skill for that.
license: MIT
compatibility: Requires Python 3.9+ and network access (to GWAS Catalog FTP / REST API). OKG provenance optional, gated on $OKG_REPO; without it the skill falls back to the GWAS Catalog REST API for build/URL metadata.
metadata:
  author: stephendorn
  version: "0.1"
  source: agentskills.io/specification
  upstream_catalog: https://www.ebi.ac.uk/gwas/
  rest_api_base: https://www.ebi.ac.uk/gwas/rest/api
---

## What this skill does

Resolves a GWAS Catalog accession (`GCST<7-digit>`) into its harmonised sumstats file and downloads it to a cache dir with a provenance sidecar. Two-tier metadata resolution:

1. **OKG-first**: if `$OKG_REPO` is set, query the OKG via the MCP server (`search(method="alias", query=<accession>)`). Reads `attrs.genome_build` + `attrs.source_url` from the matched `dataset_metadata` or `paper` node.
2. **GWAS Catalog REST API fallback**: if the OKG has no node, hit `https://www.ebi.ac.uk/gwas/rest/api/studies/<accession>` for metadata + construct the harmonised FTP URL (`https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<author>_<pubmed>/<accession>/harmonised/...`).

Output: the downloaded `.h.tsv.gz` (or whatever the catalog provides) plus `<file>.fetch.json` recording the resolution path, SHA-256, and OKG node IDs if any.

## Inputs you may need to elicit

- **`--gcst <accession>`** (required). Examples: `GCST90704615`, `GCST006907`.
- **`--cache-dir <path>`** — default `~/.cache/gwas-catalog/`. Override for a project-local cache.
- **`--refresh`** — re-download even if cached.

## How to execute

1. Run: `python3 scripts/gwas_fetch.py --gcst <accession>` with any overrides.
2. The script:
   a. Queries OKG (if `$OKG_REPO` set) → falls back to REST API on miss.
   b. Resolves the harmonised sumstats URL.
   c. Downloads with `urllib.request.urlretrieve` (skip if cached).
   d. Computes SHA-256.
   e. Writes `<cached>.fetch.json` provenance manifest.
3. Report back:
   - Cached file path
   - Genome build, sample size, ancestry (from OKG or REST API)
   - Path to the sidecar `.fetch.json`
   - Suggested next step: if the user's target LD panel is on a different build, chain into the `liftover` skill via `python3 ~/.claude/skills/liftover/scripts/lift_sumstats.py --in <cached> --source <build> --target <other-build>`.

## Refusal triggers

- Neither OKG nor REST API yields a `genome_build`. The script writes a coverage-gap proposal stub at `okg-coverage-stubs/add-dataset-<slug>/proposal.md` (same pattern as the `liftover` skill).
- REST API returns 404 / non-JSON. Cite the URL in the error.
- Cached file present but SHA-256 mismatches a previously recorded value (signals a silent re-host or corruption).

## Examples

```bash
# Fetch with OKG provenance
OKG_REPO=~/Lab/KG/okg python3 scripts/gwas_fetch.py --gcst GCST90704615

# Fetch without OKG (REST API only)
python3 scripts/gwas_fetch.py --gcst GCST006907 --cache-dir ./GWAS/raw

# Fetch then lift in one shot (chain with liftover skill)
python3 scripts/gwas_fetch.py --gcst GCST90704615 --cache-dir ./GWAS/raw
python3 ~/.claude/skills/liftover/scripts/lift_sumstats.py \
    --in ./GWAS/raw/<file>.tsv.gz \
    --okg-ld-panel-id ld_panel:gctb_ukb_ldm13m
```

## Sidecar manifest schema

`<cached>.fetch.json`:

```json
{
  "accession": "GCST90704615",
  "download_url": "https://ftp.ebi.ac.uk/pub/databases/gwas/...",
  "cached_path": "...",
  "sha256": "...",
  "genome_build": "hg38",
  "okg_node_ids": {"dataset": "...", "paper": "..."},
  "provenance": {"source": "okg" | "gwas_catalog_api", "...": "..."},
  "captured_at": "..."
}
```

## Notes

- The skill is intentionally **fetch-only**; it does not lift coordinates. For build alignment use the `liftover` skill. Keeping fetch and lift as separate skills means each is composable in batch pipelines.
- The harmonised sumstats convention (chrom + base_pair_location + GRCh38 by default for newer GWAS Catalog studies) is the GWAS Catalog's choice; downstream tools should not assume the schema beyond what the file actually contains.
- For very large studies, downloads may take many minutes; the script prints progress.

See [references/CATALOG_API.md](references/CATALOG_API.md) for the GWAS Catalog REST API endpoints the skill consults.
