# SBayesRC notes

Operational notes for the SBayesRC backend of the `prs` skill. Canonical
reference: [Zheng et al. 2024, *Nat Genet*](https://doi.org/10.1038/s41588-024-01704-y)
(PMID 38689000). Code: <https://github.com/zhilizheng/SBayesRC>.

## What the orchestrator does

`scripts/prs.py --method sbayesrc` runs through:

1. **Resolve LD panel** — via OKG: `--ancestry eur|eas|afr` → matching
   `ld_panel:sbayesrc_hm3_<anc>` node; reads its `source_url` +
   `local_path_hint`. Override with `--okg-ld-panel-id`.
2. **Cache + extract LD reference** — downloads the zip to
   `~/.cache/sbayesrc/` on first use, unzips to a sibling directory. Cache
   keyed by the panel's `local_path_hint` so subsequent runs are offline.
3. **Convert sumstats → COJO** — auto-detects harmonised TSV / LDSC munged
   / GWAS-SSF, emits the 8-column COJO TSV. See [COJO_FORMAT.md](COJO_FORMAT.md).
4. **Hand off to R** — calls `scripts/sbayesrc_run.R` via `Rscript`, which
   loads the SBayesRC R package and runs the 3-step `tidy → impute →
   sbayesrc` pipeline.
5. **Parse + manifest** — reads the R-side log for fit stats (`hsq`, `Pi`,
   MCMC iter count) and writes `<out>.prs.json` with full OKG provenance.

## The R-side pipeline

SBayesRC's R API exposes three sequential entry points:

```r
SBayesRC::tidy(mafile, LDdir, output)        # align sumstats to LD scaffold
SBayesRC::impute(mafile, LDdir, output)       # impute missing SNPs
SBayesRC::sbayesrc(mafile, LDdir, annot, output)  # MCMC sampler
```

The final call writes `<output>.snpRes` — a per-SNP TSV with posterior
mean effect sizes (the PRS weights). The skill's manifest cites this file
as `weights_path`.

## R install

On first invocation the R script auto-installs:

- `optparse` (for argument parsing) from CRAN
- `SBayesRC` from `github.com/zhilizheng/SBayesRC` via `remotes::install_github`

R 4.0+ is required. The package builds native C++ code so a working
toolchain is needed (Xcode CLT on macOS; gfortran for the underlying
matrix routines).

### macOS arm64 build fix (R 4.5 + Apple clang 14+)

Two compile errors block a stock install:

1. **`-Wc++11-narrowing` promoted to error** by clang 14 — fails on
   `SBayesRC.cpp:157  uint32_t size[2] = {thinIter, m};` with
   *non-constant-expression cannot be narrowed from int to uint32_t*.

2. **BH (Boost) headers require C++17** — `std::is_null_pointer`,
   `std::is_final`, `std::remove_cv_t` are unavailable under
   `-std=gnu++11` which the package declares.

Both are fixed by appending to `~/.R/Makevars`:

```make
CXX11FLAGS += -Wno-c++11-narrowing
CXX14FLAGS += -Wno-c++11-narrowing
CXX17FLAGS += -Wno-c++11-narrowing
CXX11 = $(CXX17)
CXX11STD = -std=gnu++17
```

This demotes the narrowing rule back to a warning and reroutes any
`CXX11` request through the C++17 toolchain (safe because C++17 is a
superset of C++11). After this patch the package installs cleanly:

```r
remotes::install_github("zhilizheng/SBayesRC", upgrade = "never")
# * DONE (SBayesRC)
library(SBayesRC)
# version 0.2.6
```

The upstream fix would be a one-line `static_cast<uint32_t>(...)` in
`SBayesRC.cpp` plus declaring `SystemRequirements: C++17` in `DESCRIPTION`.

## BaselineLD annotation file

SBayesRC requires a functional-genomic annotation file (the BaselineLD
2.2 model, ~50 MB after unzip). The R script downloads it from
`https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotations/baselineLD_2.2.annot.txt.gz`
to `~/.cache/sbayesrc/annot_baseline2.2.txt` on first use. Subsequent runs
reuse the cache.

This is a separate artifact from the per-ancestry LD eigendecomposition
references. The annotation is ancestry-agnostic.

## LD references in the OKG

| OKG node | Source URL | Build | Size |
|---|---|---|---|
| `ld_panel:sbayesrc_hm3_eur` | `gctbhub.cloud.edu.au/.../ukbEUR_HM3.zip` | GRCh37 | 3.1 GB |
| `ld_panel:sbayesrc_hm3_eas` | `gctbhub.cloud.edu.au/.../ukbEAS_HM3.zip` | GRCh37 | 2.4 GB |
| `ld_panel:sbayesrc_hm3_afr` | `gctbhub.cloud.edu.au/.../ukbAFR_HM3.zip` | GRCh37 | 5.0 GB |

All three derive from UK Biobank genotypes intersected with HapMap3 SNPs.
Each `.zip` contains per-LD-block eigen-decomposition matrices
(`block*.info` + `block*.eigen.bin`), pre-rotated low-rank factors that
SBayesRC consumes directly.

## Why the LD reference matters

SBayesRC's MCMC sampler operates on the eigendecomposition of the LD
matrix per block — not the raw correlation matrix. The low-rank truncation
(top eigenvalues) is what makes it tractable on 1.1M+ HapMap3 SNPs.
**Using an LD reference whose ancestry doesn't match the GWAS will
produce biased weights** — a non-trivial form of cross-ancestry portability
failure that LDSC papers over via the regression intercept but SBayesRC
does not. The OKG's `ld_panel.ancestry_scope` attr is the discriminator
the skill checks against.

## Container alternative (Apptainer)

`docker://zhiliz/sbayesrc` is an x86_64-only Apptainer/Docker image that
bundles SBayesRC + all R deps. The orchestrator currently defaults to the
**pure-R backend** for portability (works on macOS arm64). To use the
container on a Linux x86_64 host (e.g. an HPC cluster), pass
`--backend apptainer` *(not implemented yet — TODO)*.

## OKG manifest

Each run writes `<out>.prs.json` citing five OKG nodes:

```
method:sbayesrc         (the algorithm)
software:sbayesrc_r     (the R+Apptainer distribution actually used)
paper:zheng_2024_sbayesrc
ld_panel:sbayesrc_hm3_<anc>
dataset:<gcst-of-input>  (only if --okg-dataset-id was passed)
```

The companion `software:sbayesrc` node (GCTB C++ upstream) is the canonical
implementation citation; we record both, mirroring the
`software:ldsc` + `software:ldsc_cbiit` dual-software pattern.
