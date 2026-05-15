---
name: locuszoom
description: Produce LocusZoom-style regional association plots from a GWAS sumstats file using the `locuszoomr` R package (Lewis & Knight, github.com/myles-lewis/locuszoomr). Selects a window by gene symbol, region (chr:start-end), or lead-SNP rsid, fetches pairwise LD from the LDlink REST API (1000G, configurable population), and renders a multi-track PDF + PNG with the association scatter on top and Ensembl gene tracks below. With `--finemap-pip` or `--finemap-sidecar`, adds a middle panel of per-variant PIPs coloured by credible-set membership — the canonical 3-panel fine-mapped-locus figure. Requires the user's LDlink API token; the skill prompts for it on first use and caches it under `~/.cache/locuszoom/ldlink_token`. OKG-aware - cites software:locuszoomr, ld_panel:thousand_genomes_<pop>, and the input dataset_metadata node in the sidecar manifest when `$OKG_REPO` is set; propagates method/software/paper/ld_panel OKG node IDs from a `--finemap-sidecar` into the locuszoom manifest. Use when a user asks for a locus plot, a regional association plot, a Manhattan close-up, a fine-mapped locus plot, or wants to visualize a specific gene / SNP / region in their GWAS sumstats (with or without overlaid fine-mapping results).
license: MIT
compatibility: Requires Python 3.9+ and R 4.0+. The `locuszoomr`, `LDlinkR`, `ensembldb`, `EnsDb.Hsapiens.v75` (hg19), and `EnsDb.Hsapiens.v86` (hg38) R packages auto-install from CRAN/Bioconductor on first use. LDlink API access requires a free token from https://ldlink.nih.gov/?tab=apiaccess. OKG provenance optional, gated on $OKG_REPO.
metadata:
  author: stephendorn
  version: "0.1"
  source: agentskills.io/specification
  locuszoomr_repo: https://github.com/myles-lewis/locuszoomr
  ldlink_api_url: https://ldlink.nih.gov/?tab=apiaccess
  okg_software_node: software:locuszoomr
---

## What this skill does

Renders a LocusZoom-style regional association plot from a GWAS sumstats file. Workflow:

1. **Resolve the LDlink API token** in this order: `--ldlink-token <token>` (and cache it), then `$LDLINK_TOKEN`, then `~/.cache/locuszoom/ldlink_token`. If none found, the script exits with `LDLINK_TOKEN_MISSING`; the agent should then ask the user for their token and re-invoke with `--ldlink-token`.
2. **Normalize sumstats** to a small data frame with columns `rsid`, `chrom`, `pos`, `p`, `beta`, `se` (auto-detected from GWAS Catalog harmonised, LDSC munged `.sumstats.gz`, COJO `.ma`, or GWAS-SSF v1.0).
3. **Select the locus window** — by `--gene <SYMBOL>`, `--region <chr:start-end>`, or `--lead-snp <rsid>` (+ `--flank <bp>`).
4. **Build the locus object** via `locuszoomr::locus(...)` with the appropriate Ensembl DB (`EnsDb.Hsapiens.v75` for hg19, `EnsDb.Hsapiens.v86` for hg38).
5. **Annotate LD** by calling `locuszoomr::link_LD(loc, token = <api_token>, pop = <pop>)`, which queries LDlink's `/LDproxy` endpoint for r² vs the lead SNP.
6. **Render** the plot to `<out>.pdf` (vector) and `<out>.png` (raster, 300 dpi).
7. **Write a sidecar manifest** at `<out>.locuszoom.json` with the locus selection, LD population, LDlink request count, OKG node IDs, output SHA-256s, and timestamps.

## Inputs you may need to elicit

| Flag | Required | Notes |
|---|---|---|
| `--sumstats <path>` | yes | Raw harmonised TSV (`*.h.tsv.gz`), LDSC munged `.sumstats.gz`, COJO `.ma`, or GWAS-SSF v1.0. |
| `--gene <SYMBOL>` | one of these three | Locus by gene symbol (e.g. `APOE`, `IRF5`). |
| `--region <chr:start-end>` | | Locus by explicit region (e.g. `19:45000000-45500000`). |
| `--lead-snp <rsid>` | | Locus centred on a SNP (e.g. `rs429358`). Combine with `--flank`. |
| `--flank <bp>` | optional | Window each side of gene/SNP (default `100000`, i.e. 100 kb). Ignored if `--region` is given. |
| `--ld-pop EUR\|AFR\|EAS\|AMR\|SAS` | optional | LDlink 1000G population for r² calc (default `EUR`). |
| `--build hg19\|hg38` | optional | Genome build of the sumstats — picks the Ensembl DB. Default `hg19`. |
| `--ldlink-token <token>` | first use only | LDlink API token. Cached on first use; later runs read from cache. |
| `--out <prefix>` | yes | Output prefix; writes `<prefix>.pdf`, `<prefix>.png`, `<prefix>.locuszoom.json`. |
| `--okg-dataset-id <dataset:...>` | optional | OKG dataset_metadata node ID; recorded in the sidecar. |
| `--finemap-pip <path>` | optional | Per-variant PIP TSV. Column auto-detection: `snp`/`rsid`/`SNP` for variant id; `pip`/`PIP`/`pip_all`/`sushie_pip_all` for PIP; `cs`/`CSIndex`/`credible_set` for credible-set membership. Adds a middle PIP panel; if `cs` is present, points are coloured by credible set. |
| `--finemap-sidecar <path>` | optional | A `.finemap.json` sidecar from the `finemap` skill; the PIP TSV is derived from its `cs_path` / `weights_path` field. Carries method/software/paper/ld_panel OKG node IDs forward into the locuszoom sidecar. |

## First-time setup: the LDlink API token

The `locuszoomr` package queries the NIH LDlink REST API for pairwise r² with the lead SNP. LDlink requires a free per-user token from https://ldlink.nih.gov/?tab=apiaccess (sign up with an email; takes seconds).

**On first invocation:**

1. Try the script with no token flag. If the user already cached one earlier, it just runs.
2. If the script exits with `LDLINK_TOKEN_MISSING`, ask the user (via `AskUserQuestion` or a plain prompt) for their LDlink token. Link them to https://ldlink.nih.gov/?tab=apiaccess if they don't have one yet.
3. Re-invoke the script with `--ldlink-token <answer>`. The script writes it to `~/.cache/locuszoom/ldlink_token` (chmod 0600) and proceeds.
4. **All subsequent runs** in this and future sessions read from the cached file silently — do not prompt again unless the user explicitly says "refresh my token" or the script reports `LDLINK_TOKEN_INVALID`.

To force a refresh: pass `--ldlink-token <new_token>` again (overwrites the cache) or delete `~/.cache/locuszoom/ldlink_token` and re-prompt.

## How to execute

1. Run: `python3 scripts/locuszoom.py --sumstats <path> --gene <SYMBOL> --out <prefix>` with any overrides.
2. Report back:
   - Output PDF + PNG paths
   - Number of SNPs in the window, lead SNP rsid + p-value
   - Number of LD pairs returned by LDlink
   - Path to the `.locuszoom.json` sidecar
   - If LDlink returned fewer than ~10 pairs (sparse coverage), flag it as a warning — the plot will still render but the colour gradient will be uninformative.

## Refusal triggers

- LDlink token missing AND `--ldlink-token` not supplied → exit code 2 with stderr `LDLINK_TOKEN_MISSING`. Agent prompts user; re-invokes.
- LDlink token rejected (HTTP 401/403 from LDlink) → exit with `LDLINK_TOKEN_INVALID`. Agent informs user, asks for a fresh token.
- Locus selection ambiguous: more than one gene symbol matches in the Ensembl DB, or `--region` is malformed.
- Sumstats window has fewer than 50 SNPs → not enough for a meaningful plot; refuses unless `--force-sparse` passed.
- All three of `--gene`, `--region`, `--lead-snp` missing → script asks the user which one to use.
- `--finemap-pip` file missing required columns (`snp`/`rsid` + `pip`/`PIP`) → skill warns, skips the PIP overlay, falls back to the 2-panel plot. Not a hard refusal.
- `--finemap-sidecar` references a `cs_path` / `weights_path` that doesn't exist on disk → hard error before the R worker runs.

## Examples

```bash
# Plot the APOE locus from a Kunkle 2019 AD sumstats file (hg19)
OKG_REPO=~/Lab/KG/okg python3 scripts/locuszoom.py \
    --sumstats ~/Lab/KG/skills-test/kg-skills/GWAS/raw/Kunkle_etal_Stage1_results.txt \
    --gene APOE --build hg19 --ld-pop EUR \
    --okg-dataset-id dataset:gcst007511_ad \
    --out ~/Lab/KG/skills-test/kg-skills/locus/AD_APOE

# Region-based selection on a GRCh38 Height sumstats (Yengo 2022)
python3 scripts/locuszoom.py \
    --sumstats ~/Lab/KG/skills-test/kg-skills/GWAS/raw/GCST90245992.h.tsv.gz \
    --region 19:45000000-45500000 --build hg38 \
    --out plots/height_APOE_region

# Lead-SNP centred (rs429358 = APOE-e4), 200kb flank
python3 scripts/locuszoom.py \
    --sumstats sumstats.tsv \
    --lead-snp rs429358 --flank 200000 --build hg19 \
    --out plots/locus_rs429358

# Fine-mapped 3-panel plot — pull PIPs + CS directly from a finemap-skill sidecar
OKG_REPO=~/Lab/KG/okg python3 scripts/locuszoom.py \
    --sumstats sumstats.tsv \
    --gene APOE --build hg19 \
    --finemap-sidecar ~/Lab/KG/skills-test/kg-skills/finemap/AD_APOE.finemap.json \
    --okg-dataset-id dataset:gcst007511_ad \
    --out plots/AD_APOE_finemapped

# Or pass a bare PIP TSV (works with any fine-mapper output that has snp+pip columns)
python3 scripts/locuszoom.py \
    --sumstats sumstats.tsv \
    --gene APOE --build hg19 \
    --finemap-pip susie_results.tsv \
    --out plots/AD_APOE_finemapped
```

## Fine-mapping overlay (3-panel mode)

When `--finemap-pip` or `--finemap-sidecar` is supplied, the rendering switches to a stacked 3-panel layout:

1. **Top** — the standard -log10(p) scatter, LD-coloured to the lead SNP (same as 2-panel mode).
2. **Middle** — per-variant PIP scatter. y-axis is `[0, 1]`. Points are coloured by credible-set membership (CS1, CS2, …) if a `cs`/`CSIndex`/`credible_set` column was found in the PIP file; variants outside any credible set are grey. The top-PIP variant is labelled. Reference lines at PIP=0.5 and PIP=0.95 are drawn dotted.
3. **Bottom** — Ensembl gene tracks (same as 2-panel mode).

Column auto-detection (case-insensitive) on the PIP TSV:

| Field | Header names tried (in order) |
|---|---|
| variant id | `snp`, `rsid`, `SNP`, `marker`, `variant` |
| PIP | `sushie_pip_all`, `pip_all`, `pip`, `PIP` |
| credible set (optional) | `CSIndex`, `cs`, `credible_set`, `cs_index` |

The PIP file is joined to the sumstats window **by rsid**. Variants in the window with no PIP entry are plotted at y=0 in grey so the absence is visible. The skill reports the join coverage (`joined N/M variants on rsid -> PIP`) — if it's much smaller than the window size, you likely have an rsid-vs-coordinate naming mismatch.

When `--finemap-sidecar` is used, the locuszoom sidecar inherits the fine-mapper's OKG node IDs under `okg_node_ids.finemap_method` / `finemap_software` / `finemap_paper` / `finemap_ld_panel`, preserving the full provenance chain from sumstats → fine-mapping → plot.

## Sidecar manifest schema

`<out>.locuszoom.json`:

```json
{
  "sumstats_input": "...",
  "sumstats_sha256": "...",
  "output_pdf": "<prefix>.pdf",
  "output_png": "<prefix>.png",
  "locus_selection": {
    "mode": "gene",
    "gene": "APOE",
    "region_chr": "19",
    "region_start": 45313000,
    "region_end": 45813000,
    "flank_bp": 100000,
    "lead_snp": "rs429358",
    "lead_p": 1.0e-300
  },
  "build": "hg19",
  "ensembl_db": "EnsDb.Hsapiens.v75",
  "ld": {
    "source": "LDlink",
    "endpoint": "https://ldlink.nih.gov/LDlinkRest/ldproxy",
    "population": "EUR",
    "lead_snp": "rs429358",
    "n_ld_pairs": 487
  },
  "n_snps_in_window": 1284,
  "okg_node_ids": {
    "software": "software:locuszoomr",
    "ld_resource": "external:ldlink_1000g_EUR",
    "dataset": "dataset:gcst007511_ad"
  },
  "locuszoomr_version": "0.3.7",
  "captured_at": "2026-05-15T..."
}
```

## Notes and edge cases

- **LD pop should match the GWAS ancestry**. Plotting an EUR GWAS with `--ld-pop EAS` gives r² values from a population that doesn't match your discovery sumstats — the colour gradient won't reflect actual LD in the GWAS. Default `EUR`; override for non-EUR studies.
- **`locuszoomr` queries LDlink one request per locus**, so quota burn is gentle (a few requests per session). The free token has generous limits.
- **Ensembl DB packages are heavy** (~150 MB each). The R script downloads them via `BiocManager::install` on first use and caches them in the system R library. Subsequent runs use the cache.
- **No conda required** — pure R + Bioconductor. The R script self-installs missing packages.
- **OKG provenance is optional**; the skill works without `$OKG_REPO`. When set, the sidecar cites `software:locuszoomr` if such a node exists; otherwise records the GitHub URL as a free-form provenance field.

See [references/LOCUSZOOMR.md](references/LOCUSZOOMR.md) for the locuszoomr function reference and LDlink endpoint details.
