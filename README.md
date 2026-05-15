# statgen-skills

A collection of Claude Code / agent skills for statistical-genetics workflows, built to the [Agent Skills](https://agentskills.io) standard. Each skill is a self-contained directory with a `SKILL.md` + helper scripts.

## Available skills

| Skill | Purpose |
|---|---|
| [`liftover/`](liftover/) | KG-aware genomic liftover for GWAS sumstats. Resolves source/target build via the OKG when available; falls back to the GWAS Catalog REST API. Refuses unknown-provenance lifts and emits a coverage-gap proposal stub. Supports `--gcst <accession>` for direct GWAS-Catalog fetch+lift. |
| [`gwas-fetch/`](gwas-fetch/) | Fetch harmonised GWAS Catalog sumstats by accession (`GCST*`). OKG-first metadata resolution with REST-API fallback. Caches downloads + writes a `.fetch.json` provenance manifest. Compose with `liftover` for build alignment. |
| [`ldsc/`](ldsc/) | Run LD Score Regression: `munge` raw sumstats, estimate SNP heritability (`h2`), compute genetic correlation (`rg`). Auto-installs the CBIIT LDSC fork; caches the canonical EUR LD-score reference; writes a `.ldsc.json` sidecar with OKG provenance and parsed key results. |
| [`finemap/`](finemap/) | SuSiE / SuShiE fine-mapping via the `mancusolab/sushie` package. Subcommands: `susie` (single-ancestry, individual genotypes), `sushie` (multi-ancestry K≥2), `sumstats` (GWAS Z + LD via `infer_sushie_ss`, works for K=1 SuSiE-RSS or K≥2 SuShiE). Auto-installs sushie + runs the bundled 3-ancestry tutorial as a smoke test; writes a `.finemap.json` sidecar with OKG provenance and parsed CS/PIP summary. |
| [`prs/`](prs/) | Polygenic risk scoring via a pluggable `--method` flag. Currently implements **SBayesRC** (Zheng 2024) through the `zhilizheng/SBayesRC` R package. OKG-resolves the matching HapMap3 LD eigendecomposition by `--ancestry eur\|eas\|afr` (or explicit `--okg-ld-panel-id`), auto-pulls `N` from `--okg-dataset-id`, converts harmonised / LDSC / GWAS-SSF sumstats to COJO format, runs SBayesRC's 3-step `tidy → impute → sbayesrc` pipeline, and writes per-SNP PRS weights + `.prs.json` sidecar citing 5 OKG nodes (method, software, paper, ld_panel, dataset). LDpred2 / PRS-CS easily wired in behind the same flag. |
| [`locuszoom/`](locuszoom/) | LocusZoom-style regional association plots via the [`locuszoomr`](https://github.com/myles-lewis/locuszoomr) R package. Selects a window by gene symbol, region (`chr:start-end`), or lead-SNP rsid; fetches pairwise r² from the **LDlink** REST API (1000G, configurable population); renders a multi-track PDF + PNG with the association scatter on top and Ensembl gene tracks below. Requires a free LDlink API token — the skill prompts on first use and caches it under `~/.cache/locuszoom/ldlink_token`. Writes a `.locuszoom.json` sidecar with locus selection, LD coverage stats, and OKG provenance. |
| [`twas/`](twas/) | Transcriptome-wide association study via **TWAS-FUSION** (Gusev et al. 2016, *Nat Genet*) against the GTEx v8 multi-tissue weights. Resolves the (`--tissue`, `--ancestry`) pair to an OKG `dataset:fusion_gtex_v8_<eur\|all>:<tissue>` node, downloads + caches the panel + the FUSION 1000G EUR LDREF on first use, runs `FUSION.assoc_test.R` per chromosome, concatenates per-gene results, and writes a `.twas.json` sidecar citing the method / software / paper / ld_panel / tissue / cohort / dataset OKG nodes. Auto-installs FUSION's R deps (`plink2R`, `glmnet`, `here`, `coloc`) on first use. With `--okg-trait-id`, queries the OKG's `tissue → trait relevant_to` mapping and iterates over every tissue relevant to that trait. |

## Installation

The skills follow the [Agent Skills](https://agentskills.io/specification) specification, so they should work with any compatible agent (Claude Code, Cursor, Goose, OpenCode, etc.).

### Claude Code (user-level)

Symlink each skill into `~/.claude/skills/`:

```bash
git clone https://github.com/svdorn/statgen-skills.git
cd statgen-skills
./install.sh                  # symlinks all skills into ~/.claude/skills/
```

Or symlink a single skill:

```bash
ln -s ~/statgen-skills/liftover ~/.claude/skills/liftover
```

After symlinking, restart Claude Code (or the agent client) so it picks up the new skills at startup.

### Other agentskills.io-compatible clients

Point your client at `statgen-skills/<skill>` or copy individual skills into the client's skills directory. Each skill's `SKILL.md` documents its `compatibility` field for any environment requirements.

## Skill development

Each skill is structured per the agentskills.io spec:

```
<skill-name>/
├── SKILL.md          # required: YAML frontmatter (name, description, ...) + body
├── scripts/          # executable helpers (Python, bash, etc.)
├── references/       # detailed docs loaded on demand
└── assets/           # templates, schemas, etc.
```

Validate a skill before committing:

```bash
# optional: uses agentskills' reference validator
skills-ref validate ./liftover
```

## OKG integration

Several skills here are designed around an **OKG** (Operational Knowledge Graph) — the upstream reference is [github.com/mitdbg/okg](https://github.com/mitdbg/okg/tree/main/deployments/statgen-analysis), but any fork or compatible deployment works. To use OKG-aware features, set the `OKG_REPO` environment variable to your local clone:

```bash
git clone https://github.com/mitdbg/okg ~/Lab/KG/okg   # or your fork
export OKG_REPO=~/Lab/KG/okg
```

The skills then spawn the OKG MCP server at `$OKG_REPO/deployments/<deployment>/server.py` (default deployment: `statgen-analysis`) and query for provenance metadata before acting. When the OKG can't resolve a needed entity, skills emit a coverage-gap proposal stub instead of guessing — see each skill's `references/COVERAGE_GAPS.md` for the workflow.

If `$OKG_REPO` is not set and a skill is invoked with an `--okg-*` flag, the skill refuses with a clear error pointing you here. Skills can still run **without** the OKG by passing all required metadata explicitly (e.g. `--source hg38 --target hg19` for liftover).

## License

MIT — see [LICENSE](LICENSE).
