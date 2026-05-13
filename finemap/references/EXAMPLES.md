# Examples and bundled test data

## The 3-ancestry tutorial (run by `--verify-install`)

When the skill is first installed (or when the user passes `--verify-install`), it runs the canonical tutorial that ships in the sushie repository at `~/.cache/sushie/repo/data/`. This is the example from <https://mancusolab.github.io/sushie/manual.html>.

### Inputs

| Ancestry | VCF                  | Pheno         | Covar         | N samples | N SNPs |
|----------|----------------------|---------------|---------------|-----------|--------|
| EUR      | `vcf/EUR.vcf`        | `EUR.pheno`   | `EUR.covar`   | 489       | 123    |
| AFR      | `vcf/AFR.vcf`        | `AFR.pheno`   | `AFR.covar`   | 639       | 129    |
| EAS      | `vcf/EAS.vcf`        | `EAS.pheno`   | `EAS.covar`   | 481       | 113    |

Variants are HapMap SNPs on a single small locus; samples are from 1000 Genomes Project (real genotypes, simulated phenotype with a planted causal).

### Command

```bash
cd ~/.cache/sushie/repo/data/
sushie finemap \
    --pheno EUR.pheno AFR.pheno EAS.pheno \
    --vcf vcf/EUR.vcf vcf/AFR.vcf vcf/EAS.vcf \
    --covar EUR.covar AFR.covar EAS.covar \
    --output ~/.cache/sushie/verify_install/test_result
```

### Expected output files

After a successful run the skill checks:

| File                          | Meaning                                                    |
|-------------------------------|------------------------------------------------------------|
| `test_result.cs.tsv`          | Credible sets (one row per CS × variant)                   |
| `test_result.weight.tsv`      | Per-variant posterior weights / PIPs                       |
| `test_result.cv.tsv`          | Cross-validated prediction accuracy                        |
| `test_result.corr.tsv`        | Estimated cross-ancestry effect-size correlations          |
| `test_result.log`             | Log file with iteration trace and final ELBO               |

The verifier passes if `test_result.cs.tsv` exists and contains at least one credible set, and there are no Python tracebacks in stderr.

### Pheno / covar file format

sushie expects whitespace-separated text with no header. The pheno file has two columns (sample_id, phenotype), the covar file has 2+columns (sample_id, covariate1, covariate2, ...). Sample IDs must match the VCF.

```text
# EUR.pheno (first 3 lines)
HG00096   -0.241
HG00097    1.183
HG00099   -0.045

# EUR.covar (first 3 lines, 2 covariates)
HG00096    0.512  -1.044
HG00097   -0.876   0.213
HG00099    0.018   0.901
```

## Single-ancestry use

For `python3 scripts/finemap.py susie ...`, supply exactly one VCF and one pheno file. Internally this is sushie's K=1 mode (degenerate SuShiE), which the sushie authors note is numerically equivalent to SuSiE-RSS at typical tolerances.

```bash
python3 scripts/finemap.py susie \
    --vcf GWAS/genotypes/MSH3.vcf \
    --pheno GWAS/phenos/T01.pheno \
    --out results/finemap/MSH3_T01
```

Outputs: `results/finemap/MSH3_T01.cs.tsv`, `.weight.tsv`, `.cv.tsv`, `.log`, `.finemap.json`.

## Multi-ancestry use

For `sushie` subcommand, supply K ≥ 2 VCFs and K pheno files (and optionally K covar files). Order must match across the three flags — the i-th VCF, pheno, and covar must all be the same ancestry.

```bash
python3 scripts/finemap.py sushie \
    --vcf  EUR_MSH3.vcf  AFR_MSH3.vcf  EAS_MSH3.vcf \
    --pheno EUR_T01.pheno AFR_T01.pheno EAS_T01.pheno \
    --covar EUR_T01.covar AFR_T01.covar EAS_T01.covar \
    --out  results/finemap/MSH3_T01_3anc
```

The extra output relative to the single-ancestry case is `*.corr.tsv` — the estimated effect-size correlation matrix between ancestries. Values near 1 mean the causal effect transfers cleanly; values near 0 mean ancestry-specific architecture.

## Pinning a sushie version

To reproduce a prior run exactly, pass `--sushie-commit <sha>`:

```bash
python3 scripts/finemap.py sushie --sushie-commit a1b2c3d --vcf ... --pheno ... --out ...
```

The skill checks out that SHA in `~/.cache/sushie/repo` before running. The SHA is captured in the `.finemap.json` sidecar so downstream consumers know which version produced the result.
