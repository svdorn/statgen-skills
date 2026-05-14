# LD options for `sumstats`-driven fine-mapping

When you have GWAS summary statistics (not individual-level genotypes) and
want to fine-map with SuSiE-RSS or SuShiE, sushie needs **LD information**
to compute credible sets. The skill exposes two paths:

## Path A — `region` subcommand (recommended)

Pass a **reference genotype panel** (`--ref-vcf`, `--ref-plink`, or
`--ref-bgen`) plus a window (`--chrom --start --end`). Sushie computes the
LD matrix internally from the reference and runs SuSiE-RSS:

```bash
python3 scripts/finemap.py region \
    --gwas-sumstats GCST004132.h.tsv.gz \
    --chrom 16 --start 50222970 --end 51222970 \
    --N 40266 \
    --ref-vcf 1000G/ALL.chr16.phase3.vcf.gz \
    --okg-dataset-id dataset:gcst004132_cd \
    --out finemap/CD_NOD2
```

You don't think about LD format — sushie reads the genotypes, computes
the SNP × SNP correlation matrix on the intersection of the GWAS SNPs
and the reference, and runs the regression.

## Path B — `sumstats` subcommand

Pre-compute LD yourself and pass it as a matrix file:

```bash
python3 scripts/finemap.py sumstats \
    --z   NOD2_eur.z.tsv \
    --ld  NOD2_eur.ld.npy \
    --n   40266 \
    --out finemap/CD_NOD2_precomp
```

Useful when:
- You already have an LD matrix from elsewhere (PRS-CS-style precomputed
  blocks, GCTB output, etc.)
- You're running multi-ancestry SuShiE (`sumstats` accepts K files per
  flag; `region` doesn't multi-yet)
- You want full control over LD-block boundaries

## Reference panels — where to get them

The skill is **build-agnostic**: it doesn't ship a reference panel,
because the right one depends on your GWAS's ancestry and genome build.

| Need | Source | Notes |
|---|---|---|
| **1000G phase3 EUR, hg19** | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/` | Per-chromosome VCFs (~1-2 GB each). 503 EUR samples (CEU/GBR/FIN/IBS/TSI). Use `tabix` to slice a locus before passing to sushie. |
| **1000G phase3 EUR, GRCh38** | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/` | Same samples, lifted to GRCh38. |
| **1000G phase3 EAS / AFR / SAS / AMR** | same FTP, filter `--keep <population.txt>` via plink | Subset by super-population. |
| **UK Biobank** | Local plink files (DNANexus / approved access) | Larger N (~350-500k) but access-restricted. |
| **HRC / TOPMed** | Same — access-restricted | |
| **GCTB UKB LDm13M (`ld_panel:gctb_ukb_ldm13m`)** | per-block binary LD files on hg37 | Already in the OKG. Blockwise binary format — would need a per-block adapter (not currently wired into this skill). |
| **SBayesRC HM3 LD eigendecomp (`ld_panel:sbayesrc_hm3_eur/eas/afr`)** | low-rank eigendecomposition per block | Already in the OKG. Used by the `prs` skill; not directly consumable by sushie (would need r-square reconstruction). |

## Slicing a 1000G chromosome VCF to a locus

The 1000G per-chromosome VCFs are huge. For a single-locus run, slice
with `tabix` first:

```bash
URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
chr=16; from=50222970; to=51222970
tabix -h "${URL}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
      "${chr}:${from}-${to}" | bgzip > NOD2.eur.vcf.gz
tabix -p vcf NOD2.eur.vcf.gz
```

Subset to EUR samples with `bcftools view -S eur.samples.txt` (sample IDs
from the 1000G integrated panel file at the same FTP host).

## Choosing the LD reference

A few rules of thumb:

- **Build must match** — hg19 GWAS sumstats + hg19 reference, or both
  GRCh38. Use the [`liftover`](../../liftover/SKILL.md) skill to align.
- **Ancestry must match** — EUR sumstats + EUR reference. Multi-ancestry
  meta-analyses (mixed EUR + EAS) are tricky; ideally fine-map per
  ancestry sub-sumstats with the matching reference, or use SuShiE in
  the `sumstats` subcommand.
- **Sample size matters less than ancestry match** — 503 EUR 1000G
  samples is enough for stable LD on common variants in HapMap3, even
  though the GWAS itself might have 40k+ cases + controls.
- **For rare variants** (MAF < 1%), 1000G LD is unstable — UKB / TOPMed
  is preferred. SuSiE-RSS is generally not the right tool for rare
  variants anyway (it assumes well-estimated LD).

## OKG provenance

When `$OKG_REPO` is set, the `region` subcommand records the LD-panel
node ID (if `--okg-ld-panel-id` is passed) plus the GWAS dataset node ID
(`--okg-dataset-id`) in `<out>.finemap.json`. Future work: add an OKG
node `ld_panel:1000g_phase3_eur_hg19_vcf` so the skill can auto-resolve a
canonical EUR reference from the graph.
