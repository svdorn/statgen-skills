# GWAS Catalog REST API endpoints

## Study record
```
GET https://www.ebi.ac.uk/gwas/rest/api/studies/<accession>
```
Returns the study JSON including `publicationInfo.pubmedId`, `publicationInfo.author.fullname`, and (newer studies) `summaryStatisticsAssembly`.

## Summary-statistics record
```
GET https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/<accession>
```
Returns sumstats-specific metadata including `summary_statistics_url` / `harmonised_summary_statistics_url` when available.

## FTP layout
Harmonised sumstats live under:
```
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<author>_<pubmed>/<accession>/harmonised/*.h.tsv.gz
```
where `<author>` is the first author's lowercased surname.

## File schema (harmonised)
Standard columns the GWAS Catalog harmonised TSV emits:
- `chromosome` (numeric, no `chr` prefix)
- `base_pair_location` (1-based)
- `effect_allele`, `other_allele`
- `beta`, `standard_error`
- `effect_allele_frequency`
- `p_value`
- variable extras: `odds_ratio`, `ci_lower`, `ci_upper`, etc.

## Build conventions
- Most studies from ~2020 onward use **GRCh38**.
- Pre-2020 studies are typically **GRCh37**.
- Always verify via the API's `summaryStatisticsAssembly` field rather than assuming.

## Citation

When using GWAS Catalog data:
> Sollis E, Mosaku A, Abid A et al. *The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource.* Nucleic Acids Res. 2023.

When `$OKG_REPO` is set the sidecar manifest cites the OKG node that resolved the metadata, so a downstream BibTeX generator (e.g. a future `cite-methods` skill) can pull the canonical citations from the OKG.
