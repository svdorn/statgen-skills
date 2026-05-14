---
name: prs
description: Build a polygenic risk score from GWAS summary statistics using a parameterised PRS method. Currently implements SBayesRC (Zheng et al. 2024, Nat Genet) via the zhilizheng/SBayesRC R package; LDpred2, PRS-CS, and PRS-CSx can be added behind the same `--method` flag. OKG-aware - resolves the ancestry-matched LD eigendecomposition reference (HapMap3 EUR/EAS/AFR) via `ld_panel:sbayesrc_hm3_<anc>` and pulls sample size from the input dataset_metadata node. Writes per-SNP PRS weights plus a sidecar manifest citing the method/software/paper/LD-panel node IDs. Use when a user asks to compute PRS weights, run SBayesRC, or build a polygenic score from sumstats.
license: MIT
compatibility: Requires Python 3.9+, R 4.0+ with `remotes`, and `unzip`. The SBayesRC R package is installed from github.com/zhilizheng/SBayesRC on first use. LD eigendecomposition references (~3-5 GB per ancestry) are downloaded from `gctbhub.cloud.edu.au` to `~/.cache/sbayesrc/` and cached. OKG provenance optional, gated on $OKG_REPO.
metadata:
  author: stephendorn
  version: "0.1"
  source: agentskills.io/specification
  sbayesrc_repo: https://github.com/zhilizheng/SBayesRC
  sbayesrc_paper_doi: 10.1038/s41588-024-01704-y
  okg_method_node: method:sbayesrc
  okg_ld_panels:
    - ld_panel:sbayesrc_hm3_eur
    - ld_panel:sbayesrc_hm3_eas
    - ld_panel:sbayesrc_hm3_afr
---

## What this skill does

Computes per-SNP polygenic-score weights from a GWAS sumstats file using a pluggable PRS method, defaulting to SBayesRC. Workflow:

1. **Resolve LD panel** — either by `--ancestry eur|eas|afr` (maps to the matching `ld_panel:sbayesrc_hm3_<anc>` node) or by explicit `--okg-ld-panel-id`. Reads the OKG node's `source_url` + `local_path_hint`; downloads + unzips the eigendecomposition on first use.
2. **Convert sumstats to COJO format** — `SNP A1 A2 freq b se p N`. Accepts (a) raw GWAS Catalog harmonised files with `hm_*` columns, (b) LDSC munged `.sumstats.gz`, or (c) pre-converted COJO TSV. N can be auto-resolved from `--okg-dataset-id` if not given.
3. **Run the chosen method** — for `--method sbayesrc`, calls the companion R script `scripts/sbayesrc_run.R` which loads the SBayesRC R package, points it at the LD folder + a baseline-LD annotation file, and writes per-SNP weights.
4. **Write provenance sidecar** — `<out>.prs.json` cites method/software/paper/ld_panel/source-dataset node IDs, parameters, key fit stats from the R log, and SHA-256 of the LD reference.

## Inputs you may need to elicit

| Flag | Required | Notes |
|---|---|---|
| `--method` | yes (default `sbayesrc`) | Currently only `sbayesrc`. Future: `ldpred2`, `prscs`, `prscsx`. |
| `--gwas-sumstats <path>` | yes | Raw harmonised TSV, LDSC munged `.sumstats.gz`, or COJO TSV. |
| `--ancestry eur\|eas\|afr` | yes (or `--okg-ld-panel-id`) | Picks the matching SBayesRC HM3 LD eigendecomposition. |
| `--okg-ld-panel-id <ld_panel:...>` | alt | Explicit OKG panel node (bypasses the ancestry heuristic). |
| `--okg-dataset-id <dataset:...>` | optional | Auto-resolves `N` for COJO from `n_cases + n_controls`, or `n_samples` if quantitative. |
| `--N <int>` | optional | Constant total sample size (overrides OKG-resolved). |
| `--out <prefix>` | yes | Output prefix; the weights file lands at `<prefix>.snpRes` (SBayesRC's native suffix). |

## How to execute

1. Run: `python3 scripts/prs.py --method sbayesrc --gwas-sumstats <path> --ancestry eur --out <prefix>` with any overrides.
2. Report back:
   - Path to the PRS weights file (`<prefix>.snpRes`)
   - Number of SNPs retained, mean/max effect size
   - Convergence indicator (MCMC iterations, final hsq estimate, polygenicity Pi)
   - Path to the `.prs.json` sidecar
   - For binary traits the `manifest.method_specific` block carries the SBayesRC variance components

## Refusal triggers

- Method is unknown (e.g. `--method ldpred2` until that backend is implemented). Skill suggests opening a feature request and shows the supported set.
- No `--ancestry` or `--okg-ld-panel-id` given; skill can't pick an LD reference.
- Sumstats has fewer than 50,000 SNPs after HM3 intersection (SBayesRC needs the full HM3 scaffold).
- LD download fails AND no `local_path_hint` already populated → skill writes a coverage-gap stub at `okg-coverage-stubs/add-ld-panel-<slug>/proposal.md` and refuses.

## Examples

```bash
# Default: SBayesRC on a EUR IBD GWAS, OKG-resolved everything
OKG_REPO=~/Lab/KG/okg python3 scripts/prs.py \
    --method sbayesrc \
    --gwas-sumstats ~/Lab/KG/skills-test/GWAS/28067908-GCST004131-EFO_0003767.h.tsv.gz \
    --ancestry eur \
    --okg-dataset-id dataset:gcst004131_ibd \
    --out ~/Lab/KG/skills-test/prs/GCST004131_sbayesrc

# Explicit LD-panel pin (skip ancestry heuristic)
python3 scripts/prs.py \
    --method sbayesrc \
    --gwas-sumstats T01.sumstats.gz \
    --okg-ld-panel-id ld_panel:sbayesrc_hm3_eas \
    --N 200000 \
    --out prs/T01_sbayesrc_eas
```

## Sidecar manifest schema

`<out>.prs.json`:

```json
{
  "method": "sbayesrc",
  "gwas_sumstats_input": "...",
  "output_prefix": "...",
  "weights_path": "<prefix>.snpRes",
  "n_snps_retained": 1142318,
  "okg_node_ids": {
    "method": "method:sbayesrc",
    "software": "software:sbayesrc_r",
    "paper": "paper:zheng_2024_sbayesrc",
    "ld_panel": "ld_panel:sbayesrc_hm3_eur",
    "dataset": "dataset:gcst004131_ibd"
  },
  "ld_reference": {
    "local_path": "~/.cache/sbayesrc/ukbEUR_HM3/",
    "source_url": "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/ukbEUR_HM3.zip",
    "sha256_zip": "..."
  },
  "method_specific": {
    "hsq": 0.121, "hsq_se": 0.011,
    "polygenicity_pi": 0.014,
    "n_mcmc_iter": 3000, "n_burnin": 1000
  },
  "captured_at": "2026-05-14T..."
}
```

## Adding a new PRS method

The orchestrator picks an implementation via `args.method`. To add `--method ldpred2`:
1. Add a new branch in `scripts/prs.py` (a small `run_ldpred2(args, ...)` function).
2. Add an OKG `software:` node + `method:` node + (if needed) an LD-panel node.
3. Document the new method's input quirks in `references/<METHOD>.md`.

See [references/SBAYESRC.md](references/SBAYESRC.md) for the SBayesRC R-side details (R package install, COJO conversion, annotation file requirement).

## Notes

- **x86_64-only Apptainer**: SBayesRC's official container `docker://zhiliz/sbayesrc` is x86_64-only. This skill defaults to the **pure-R backend** for portability (works on macOS arm64). Pass `--backend apptainer` to use the container instead — only valid on x86_64 Linux hosts (typically a cluster). The OKG records both `software:sbayesrc` (GCTB C++) and `software:sbayesrc_r` (R/Apptainer); the manifest cites whichever was used.
- **Baseline-LD annotation file** is required by SBayesRC for the functional-prior step. The R script downloads it from the GCTB host on first use (~50 MB) and caches at `~/.cache/sbayesrc/annot_baseline2.2.txt`.
- **LD reference sizes**: EUR 3.1 GB, EAS 2.4 GB, AFR 5.0 GB (zipped). Unzipped, each expands to ~2× that (per-block eigen-decomposition matrices). Plan disk accordingly.
