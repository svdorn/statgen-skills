---
name: twas
description: Run a transcriptome-wide association study (TWAS) on a GWAS sumstats file using TWAS-FUSION (Gusev et al. 2016, Nat Genet) with pre-computed GTEx v8 multi-tissue expression weights. Resolves the weight panel via the OKG (`--tissue Brain_Cortex --ancestry EUR` → `dataset:fusion_gtex_v8_eur:brain_cortex`), downloads + caches the panel + the FUSION 1000G LDREF on first use, runs `FUSION.assoc_test.R` per chromosome, concatenates per-gene results, and writes a `.twas.json` sidecar citing the method/software/paper/ld_panel/dataset/tissue/cohort OKG nodes. OKG-aware - if `$OKG_REPO` is set, can also auto-suggest tissues for a given trait via the `tissue → trait relevant_to` mapping (e.g., AD → Brain_Cortex / Brain_Hippocampus / Brain_Frontal_Cortex_BA9 / Whole_Blood). Use when a user asks to run TWAS, compute gene-level association from GWAS sumstats, integrate expression weights with a GWAS, or fine-map a GWAS locus to candidate causal genes.
license: MIT
compatibility: Requires Python 3.9+ and R 4.0+ (FUSION is R). Installs the `gusevlab/fusion_twas` repo + the FUSION 1000G LDREF (~120 MB) + the requested GTEx v8 weight panel (~50-200 MB each) on first use, cached under `~/.cache/twas/`. R deps `plink2R`, `glmnet`, `methods`, `optparse` auto-install on first use. PLINK1.9 binary required (skill auto-installs from the FUSION-shipped utilities or `brew install plink` on macOS). OKG provenance optional, gated on `$OKG_REPO`.
metadata:
  author: stephendorn
  version: "0.1"
  source: agentskills.io/specification
  fusion_repo: https://github.com/gusevlab/fusion_twas
  fusion_paper_doi: 10.1038/ng.3506
  fusion_paper_pmid: "26854917"
  okg_method_node: method:fusion_twas
  okg_software_node: software:fusion
  okg_paper_node: paper:fusion_2016
  okg_ld_panel_node: ld_panel:fusion_1000g_eur
  okg_cohort_node: cohort:gtex_v8
---

## What this skill does

Wraps `FUSION.assoc_test.R` behind a single CLI. Workflow:

1. **Resolve the FUSION weight panel** for a (`--tissue`, `--ancestry`) pair via the OKG: `--tissue Brain_Cortex --ancestry EUR` → `dataset:fusion_gtex_v8_eur:brain_cortex` → reads the node's `source_url` (the S3 tar.gz). User can also pass `--okg-dataset-id` directly to bypass the convenience lookup.
2. **Download + unpack** the panel to `~/.cache/twas/fusion_gtex_v8_<anc>/<tissue>/` on first use; cached afterwards.
3. **Download + unpack** the FUSION 1000G EUR LDREF (`ld_panel:fusion_1000g_eur`) to `~/.cache/twas/ldref/` on first use; cached afterwards.
4. **Normalize sumstats** to the FUSION-expected `SNP A1 A2 Z` layout (auto-detects from harmonised TSV, LDSC munged, COJO `.ma`, or generic; computes `Z` from `beta/SE` or `OR` when needed).
5. **Run `FUSION.assoc_test.R`** for each chromosome that has weights in the panel's `*.pos` file (auto-iterates 1..22 + X if present). Concatenates the per-chromosome `.dat` outputs into a single `<out>.assoc.tsv`.
6. **Write a sidecar manifest** at `<out>.twas.json` with the OKG node IDs (method/software/paper/ld_panel/dataset/tissue/cohort), the panel SHA-256, parsed key stats (n_genes_tested, n_significant_at_bonferroni, top hits), and reproducibility metadata (FUSION repo commit, panel release).

## Inputs you may need to elicit

| Flag | Required | Notes |
|---|---|---|
| `--sumstats <path>` | yes | GWAS sumstats: harmonised `.h.tsv.gz`, LDSC munged `.sumstats.gz`, COJO `.ma`, or a TSV with at minimum `SNP/A1/A2/Z` or `SNP/A1/A2/BETA/SE`. The skill normalises to FUSION's `SNP A1 A2 Z` format internally. |
| `--tissue <Name>` | yes (or `--okg-dataset-id`) | FUSION-style tissue name (e.g. `Brain_Cortex`, `Whole_Blood`, `Adipose_Subcutaneous`). Resolves to the GTEx v8 panel via OKG when `--ancestry` is given. |
| `--ancestry EUR\|ALL` | yes (or `--okg-dataset-id`) | Picks `dataset:fusion_gtex_v8_eur:<tissue>` (EUR-only) or `dataset:fusion_gtex_v8_all:<tissue>` (multi-ancestry). EUR is the default match for the rest of the panel. |
| `--okg-dataset-id <dataset:...>` | alt | Explicit OKG panel node ID (bypasses the tissue/ancestry heuristic). |
| `--okg-trait-id <trait:...>` | optional | If passed without `--tissue`, the skill queries the OKG for `tissue → trait relevant_to` edges and runs TWAS over every relevant tissue (~3-4 per trait), concatenating results into one table. |
| `--use-filtered-pos` | optional (default true) | Use the panel's `*.pos` filtered for genes with significant heritability (the FUSION-recommended default) rather than `no_filter.pos` (all genes). |
| `--chr <list>` | optional | Comma-separated chromosomes to run (default: every chromosome with weights in the panel). |
| `--coloc-p <float>` | optional | If set, runs the COLOC sub-test for any gene with TWAS.P below this threshold. |
| `--perm <int>` | optional | Permutation count for the FUSION permutation test (`--perm` flag in `FUSION.assoc_test.R`). |
| `--out <prefix>` | yes | Writes `<prefix>.assoc.tsv`, `<prefix>.twas.json`, `<prefix>.log`. |
| `--okg-repo <path>` | optional | OKG repo path (honors `$OKG_REPO`). Required if you want OKG-resolution of the panel. |

## How to execute

1. Run: `python3 scripts/twas.py --sumstats <path> --tissue <Name> --ancestry EUR --out <prefix>` with any overrides.
2. Report back:
   - Path to `<prefix>.assoc.tsv` (the concatenated per-gene TWAS results).
   - `n_genes_tested`, `n_significant` (Bonferroni-significant at 0.05 / n_genes), and the top 5 hits with their TWAS.Z + TWAS.P.
   - Path to `<prefix>.twas.json` sidecar.
   - If `--okg-trait-id` was used, the table is per-(gene × tissue), and the report breaks down the top hits per tissue.

## Refusal triggers

- Sumstats missing both Z-score and beta/SE columns → can't construct Z; refuses.
- Neither `--tissue + --ancestry` nor `--okg-dataset-id` supplied → can't pick a panel; refuses.
- Panel resolves to a node that has no `source_url` attr → refuses with a coverage-gap stub at `okg-coverage-stubs/add-twas-panel-<slug>/proposal.md`.
- After running, fewer than 100 genes were testable in the panel (too small for a meaningful TWAS); flags a warning but doesn't refuse.

## Examples

```bash
# Kunkle 2019 AD × GTEx v8 Brain Cortex (EUR), OKG-resolved
OKG_REPO=~/Lab/KG/okg python3 scripts/twas.py \
    --sumstats ~/Lab/KG/skills-test/kg-skills/GWAS/raw/Kunkle_etal_Stage1_results.txt \
    --tissue Brain_Cortex --ancestry EUR \
    --okg-dataset-id dataset:gcst007511_ad \
    --out ~/Lab/KG/skills-test/kg-skills/twas/AD_Brain_Cortex

# Yengo 2022 Height × all relevant tissues (auto-resolved from the trait→tissue mapping)
python3 scripts/twas.py \
    --sumstats height.h.tsv.gz \
    --okg-trait-id trait:height \
    --out twas/height_multi_tissue

# Explicit OKG-dataset pin (skip tissue heuristic)
python3 scripts/twas.py \
    --sumstats t01.sumstats.gz \
    --okg-dataset-id dataset:fusion_gtex_v8_eur:liver \
    --out plots/T01_liver
```

## Sidecar manifest schema

`<out>.twas.json`:

```json
{
  "sumstats_input": "...",
  "sumstats_sha256": "...",
  "output_assoc": "<prefix>.assoc.tsv",
  "output_log": "<prefix>.log",
  "tissue": "Brain_Cortex",
  "ancestry": "EUR",
  "fusion_repo": "https://github.com/gusevlab/fusion_twas",
  "fusion_commit": "<git sha>",
  "ldref_sha256": "...",
  "panel_sha256": "...",
  "n_genes_tested": 4783,
  "n_significant_bonferroni": 12,
  "top_hits": [
    {"id": "APOE", "chr": 19, "twas_z": 18.4, "twas_p": 1.2e-75, "best_eqtl": "rs429358"}
  ],
  "okg_node_ids": {
    "method": "method:fusion_twas",
    "software": "software:fusion",
    "paper": "paper:fusion_2016",
    "ld_panel": "ld_panel:fusion_1000g_eur",
    "tissue": "tissue:brain_cortex",
    "cohort": "cohort:gtex_v8",
    "dataset_panel": "dataset:fusion_gtex_v8_eur:brain_cortex",
    "dataset_gwas": "dataset:gcst007511_ad"
  },
  "captured_at": "2026-05-15T..."
}
```

## Notes and edge cases

- **Panel build is GRCh38**; FUSION's LDREF is hg19. FUSION handles the build mismatch internally via the panel's `pos` file. If your sumstats are on hg19 (e.g., LDSC munged), no manual liftover is needed for the TWAS step — but if you pre-liftover to hg38, that's fine too, FUSION joins by rsid.
- **EUR vs ALL panels**: EUR-only weights are stronger for EUR-ancestry GWAS (which is most of our panel); ALL is the safe choice for trans-ancestry studies. Default is EUR if neither `--ancestry` nor `--okg-dataset-id` is specified.
- **Filtered vs no_filter pos**: each archive ships two `*.pos` files — one for genes with significant cis-heritability (FUSION-recommended; the default) and one for all genes (`no_filter` suffix). Pass `--no-filter` to switch.
- **Heavy first run**: downloads the FUSION repo + the LDREF (~120 MB) + the requested panel (~50-200 MB). Subsequent runs reuse the cache.

See [references/FUSION.md](references/FUSION.md) for the FUSION command-line reference, sumstats column convention, and the assoc_test output schema.
