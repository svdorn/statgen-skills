# statgen-skills

A collection of Claude Code / agent skills for statistical-genetics workflows, built to the [agentskills.io](https://agentskills.io) standard. Each skill is a self-contained directory with a `SKILL.md` + helper scripts.

## Available skills

| Skill | Purpose |
|---|---|
| [`liftover/`](liftover/) | KG-aware genomic liftover for GWAS sumstats. Resolves source/target build via the statgen-analysis OKG when available; falls back to the GWAS Catalog REST API. Refuses unknown-provenance lifts and emits a coverage-gap proposal stub. Supports `--gcst <accession>` for direct GWAS-Catalog fetch+lift. |

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
ln -s ~/Lab/KG/statgen-skills/liftover ~/.claude/skills/liftover
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

Several skills here are designed around the **statgen-analysis OKG** at [github.com/mitdbg/okg](https://github.com/mitdbg/okg/tree/main/deployments/statgen-analysis). They spawn the OKG MCP server (default path: `/Users/stephen/Lab/KG/okg`; override via `OKG_REPO` env var) and query for provenance metadata before acting. When the OKG can't resolve a needed entity, skills emit a coverage-gap proposal stub instead of guessing — see each skill's `references/COVERAGE_GAPS.md` for the workflow.

## License

MIT — see [LICENSE](LICENSE).
