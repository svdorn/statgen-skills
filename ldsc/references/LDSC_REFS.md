# LDSC reference files

## EUR (default)
- `eur_w_ld_chr.tar.bz2` — https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
  - 1000G Phase 3 EUR samples, HapMap3 SNPs
  - Used as both `--ref-ld-chr` and `--w-ld-chr` for h2/rg

## EAS
- `eas_ldscores.tar.bz2` — https://data.broadinstitute.org/alkesgroup/LDSCORE/eas_ldscores.tar.bz2
  - 1000G Phase 3 EAS samples; pass with `--ld-scores-dir <extracted-dir>`

## AFR
- Pre-computed AFR LD scores are not bundled by the upstream LDSC project; consider [PolyFun's ancestry-specific scores](https://github.com/omerwe/polyfun) or compute your own with `ldsc.py --l2`.

## Partitioned (stratified) LD scores
- `1000G_Phase3_baselineLD_v2.2_ldscores.tgz` — https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
  - For stratified h2 partitioning by functional annotation (`--h2-cts`, `--ref-ld-chr-cts`)

## Citation

When you use LDSC outputs in a manuscript, cite:

- **Heritability / intercept**: Bulik-Sullivan et al. 2015 *Nat Genet* (`paper:ldsc_2015`). DOI: `10.1038/ng.3211`.
- **Genetic correlation**: Bulik-Sullivan et al. 2015 *Nat Genet* (rg paper). DOI: `10.1038/ng.3406`.
- **Stratified h2**: Finucane et al. 2015 *Nat Genet*. DOI: `10.1038/ng.3404`.

The skill emits these in the `okg_node_ids` block of the sidecar manifest when `$OKG_REPO` is set, so you can pull them into a BibTeX via the OKG.

## Fork notes

- **`CBIIT/ldsc`** (default for this skill) — Python 3 + macOS-arm64 compatible. Maintained by NCI's Center for Biomedical Informatics.
- **`bulik/ldsc`** (canonical, OKG default) — Python 2.7. Useful for reproducing pre-2020 results that ran on the upstream code; otherwise prefer CBIIT.
- **`belowlab/ldsc`** — alternative Python 3 fork; works the same.

Pass `--repo-url <url>` to override the default.
