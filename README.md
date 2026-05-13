# statgen-skills

A collection of Claude Code / agent skills for statistical-genetics workflows, built to the [Agent Skills](https://agentskills.io) standard. Each skill is a self-contained directory with a `SKILL.md` + helper scripts.

## Available skills

| Skill | Purpose |
|---|---|
| [`liftover/`](liftover/) | KG-aware genomic liftover for GWAS sumstats. Resolves source/target build via the OKG when available; falls back to the GWAS Catalog REST API. Refuses unknown-provenance lifts and emits a coverage-gap proposal stub. Supports `--gcst <accession>` for direct GWAS-Catalog fetch+lift. |
| [`gwas-fetch/`](gwas-fetch/) | Fetch harmonised GWAS Catalog sumstats by accession (`GCST*`). OKG-first metadata resolution with REST-API fallback. Caches downloads + writes a `.fetch.json` provenance manifest. Compose with `liftover` for build alignment. |
| [`ldsc/`](ldsc/) | Run LD Score Regression: `munge` raw sumstats, estimate SNP heritability (`h2`), compute genetic correlation (`rg`). Auto-installs the CBIIT Python 3 / Mac-compatible LDSC fork; caches the canonical EUR LD-score reference; writes a `.ldsc.json` sidecar with OKG provenance and parsed key results. |

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
