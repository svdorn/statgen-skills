# Coverage-gap workflow

When the liftover skill cannot resolve a dataset's `genome_build` from the OKG and the user didn't supply `--source` explicitly, it refuses the lift and writes a paired-OpenSpec proposal stub at:

```
<input-dir>/okg-coverage-stubs/add-dataset-<slug>/proposal.md
```

This mirrors the pattern used by the statgen-analysis benchmark harness's `okg_provenance.py`: rather than guess at provenance, surface the gap as a stub the operator can copy into the okg repo.

## Operator workflow

1. **Read the stub** — it describes the missing OKG node and what attrs it should carry (`genome_build`, `provider`, `source_url`, `release`, `access_posture`).

2. **Decide where the node belongs**. For most GWAS Catalog studies, the natural home is `okg/deployments/statgen-analysis/methods/methods_literature.yaml` (under `papers:` if it's a published study; under a future `datasets:` block if it's a standalone sumstats record). For LD panels, fixtures/proof_slice.json is the established home (see `add-ukb-ldm13m-panel`).

3. **Scaffold a paired OpenSpec change** in the okg repo:

   ```bash
   cd /Users/stephen/Lab/KG/okg
   openspec new change add-dataset-<slug>
   ```

4. **Author the proposal** — copy the stub body into the new change's `proposal.md` and expand. Include:
   - **Why**: a sentence on why the dataset is now in scope (e.g. "the liftover skill refused trait T01 because the OKG had no node for GCST90704615").
   - **What Changes**: the exact YAML rows being added.
   - **Capabilities** (modified): typically `finemap-method-literature-coverage` extended to cover the new dataset.

5. **Edit the OKG content**:
   - Add the row to `methods_literature.yaml` (or appropriate fixture).
   - Confirm no schema/narrowing/invariant files touched (content-only edits).

6. **Publish a new generation**:

   ```bash
   uv run okg catalog load --apply --deployment deployments/statgen-analysis --dsn "$DSN"
   uv run okg ingest --deployment statgen-analysis --dsn "$DSN" --no-publish --progress
   uv run okg run --once --apply --deployment statgen-analysis --dsn "$DSN"
   ```

7. **Verify**:

   ```bash
   uv run okg doctor --check integrity --generation <new-gen> --json
   ```

   Plus MCP `get_node` on the new id to confirm `attrs.genome_build` is populated.

8. **Re-run the liftover skill** with the same `--okg-dataset-id <new-id>` (or `--gcst <accession>`, which will now resolve via the OKG path). The lift should proceed cleanly.

## When NOT to use this workflow

- **Ad-hoc lifts where provenance doesn't matter** — pass `--source` and `--target` explicitly; the skill won't try OKG resolution and won't refuse.
- **Datasets that are already in the OKG but the alias just doesn't match** — try `--okg-dataset-id <node-id>` directly. If that resolves, the issue is just alias coverage, not a real OKG gap.

## Where this pattern came from

The `add-finemap-method-coverage` (gen 4) and `add-ukb-ldm13m-panel` (gen 5) changes in the okg repo are the canonical worked examples of this paired-OpenSpec pattern. Use them as templates when authoring new coverage-gap fixes.
