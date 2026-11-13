# LPS_typing_Illumina
Bioinformatics pipeline for *Pasteurella multocida* LPS typing using Illumina sequencing data

- [Quick start](#quick-start)
- [Overall pipeline](#overall-pipeline)
- [Step-by-step user guide](#step-by-step-user-guide)
- [Database setup](#database-setup)
- [Optional parameters](#optional-parameters)
- [Output files](#structure-of-the-output-folders)
- [Publication figures](#publication-figures)
- [Advanced use](#advanced-use)
- [Acknowledgements / citations / credits](#acknowledgements--citations--credits)

---

## Quick start

**Requirements:** [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 25.04.6, and either [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) or [Apptainer](https://apptainer.org/docs/user/latest/).

```bash
# 1. Clone the repository
git clone https://github.com/julianzaugg/LPS_typing_Illumina.git
cd LPS_typing_Illumina

# 2. Obtain large databases (see Database setup below)
#    The LPS and Kaptive databases are already included in databases/

# 3. Create your samplesheet (see Samplesheet section below)

# 4. Run (local execution with Apptainer)
nextflow run main.nf -profile apptainer \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  -resume

# On a SLURM cluster, add the slurm profile and your account name:
nextflow run main.nf -profile apptainer,slurm \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  --slurm_account YOUR_ACCOUNT \
  -resume
```

---

## Overall pipeline

![Pipeline overview](pipeline_flowchart.svg)

### 1. Read trimming

Raw Illumina reads are trimmed using [fastp](https://github.com/OpenGene/fastp) v0.24.0.

### 2. Illumina reads quality metrics

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.12.1 computes quality metrics for each sample on the trimmed reads. [MultiQC](https://github.com/MultiQC/MultiQC) v1.28 produces a combined report across all samples.

### 3. Genome assembly using Shovill

Paired-end reads are assembled using [Shovill](https://github.com/tseemann/shovill) v1.4.2, which uses SPAdes at its core.

### 4. Assembly quality assessment with QUAST

[QUAST](https://quast.sourceforge.net/quast.html) v5.2.0 computes assembly metrics on the polished assemblies.

### 5. Assembly quality assessment with CheckM

[CheckM](https://github.com/Ecogenomics/CheckM) v1.2.2 (command `lineage_wf`) estimates genome completeness and contamination from marker genes.

### 6. Sylph taxonomy classification

Trimmed reads are classified with the containment-based taxonomy profiler [Sylph](https://sylph-docs.github.io), and taxonomic labels are assigned with [sylph-tax](https://github.com/bluenote-1577/sylph-tax). Both the GTDB R232 and RefSeq Fungi databases are used, downloaded from the [Sylph pre-built databases](https://sylph-docs.github.io/pre%E2%80%90built-databases/).

### 7. LPS typing using Kaptive

The LPS type is determined from the assembled genome using [Kaptive](https://kaptive.readthedocs.io/en/latest/) v3 with the 9-LPS database (included in `databases/kaptive3_LPS_db_v1/`).

### 8. Variant calling using Snippy

- Trimmed reads are mapped to the LPS reference sequence identified by Kaptive using BWA-mem as part of [Snippy](https://github.com/tseemann/snippy) v4.6.0.
- Snippy calls substitutions and indels relative to the reference LPS locus.
- [SnpEff](https://pcingola.github.io/SnpEff/) annotates variant effects on genes and proteins.
- High-impact variants (frameshift, stop_gained) are extracted to a separate file.

#### Majority-allele 0/1 calls

Snippy runs freebayes in diploid mode and keeps only homozygous-alt (`GT=1/1`) calls, which is its intended design ([tseemann/snippy#343](https://github.com/tseemann/snippy/issues/343)).
A real allele with a minority of reference-supporting reads is genotyped `0/1` and never reaches `snps.tab`.
This typically affects indels in tandem repeats, where a fraction of reads support a different repeat length, and can hide a typing allele.
No snippy option recovers such calls: `--minfrac` already defaults to 0, and haploid freebayes (`-p 1`) writes `GT=1`, which snippy's filter then drops entirely.

The pipeline therefore rebuilds `snps.tab` from `snps.raw.vcf` with `bin/snippy_typing_calls.sh`:

- Records passing snippy's own rule are kept unchanged (`GT=1/1`, QUAL >= `--snippy_minqual`, depth >= `--snippy_mincov`, AO/DP >= `--snippy_minfrac`).
- Biallelic `GT=0/1` records with the same QUAL and depth and AO/DP >= `--snippy_het_min_af` (default 0.5) are added.
- The selected records go through the same `vt normalize`, snpEff and `snippy-vcf_to_tab` steps snippy uses, so the table schema and subtype matching are unchanged.

The 0.5 default keeps the majority allele, which is the genotype of a haploid isolate.
Per-position sequencing noise at typical depths sits well below this, so the value mainly decides how minority alleles are treated rather than separating signal from noise.

A rescued call means the majority allele is present.
It can also reflect a minor subpopulation, library-prep stutter or a repeat-length difference between preparations.
Each one is listed in `sample_id_snps.rescued.tsv` with its allele fraction and a `heterogeneous` flag (AF < 0.9).

Notes:

- `snps.raw.vcf`, `snps.filt.vcf` and the BAM are snippy's own and untouched; snippy's original table is kept as `sample_id_snps.snippy.tab`.
- Snippy's consensus and alignment FASTA files (not published) do not include rescued alleles.
- Set `--snippy_het_min_af false` to reproduce plain snippy output exactly.
- The snippy process ignores errors, so a sample whose snippy or typing step fails is missing from the subtype report rather than stopping the run; check `sample_id_snippy.log`.
- `tests/test_snippy_typing_calls.sh` checks the script against a fixture snippy output (see the header of that file for how to run it in the snippy container).

### 9. MLST typing

[mlst](https://github.com/tseemann/mlst) scans assemblies against the PubMLST scheme `pmultocida_2` (RIRDC) by default. The scheme can be changed with `--mlst_scheme`.

The database bundled inside the mlst container is frozen at the time of the container build. PubMLST now requires a login to download updated alleles and profiles, so the tool's built-in update mechanism no longer works. To use the latest scheme, download files manually and supply them via `--mlst_datadir` (see below).

#### Updating the MLST database

1. Log in to [PubMLST](https://pubmlst.org/organisms/pasteurella-multocida) and navigate to the **pmultocida_2 (RIRDC)** scheme.
2. Download allele sequences for each of the 7 loci (`RIRDC_adk`, `RIRDC_est`, `RIRDC_gdh`, `RIRDC_mdh`, `RIRDC_pgi`, `RIRDC_pmi`, `RIRDC_zwf`). PubMLST provides these as `<locus>.fas` FASTA files.
3. Download the allele profiles as `pmultocida_2.txt`.
4. Place all files in `databases/mlst_db/pmultocida_2/`:
   ```
   databases/mlst_db/pmultocida_2/
     RIRDC_adk.fas
     RIRDC_est.fas
     ...
     pmultocida_2.txt
   ```
5. Run the pipeline with `--mlst_datadir databases/mlst_db`. On first run, the pipeline automatically renames allele files to the `.tfa` extension that mlst expects, and builds the BLAST index. The prepared database is saved to `databases/mlst_db_prepared/`. Use `-resume` on subsequent runs to skip rebuilding.

A pre-downloaded database (2026-08-18) is included at `databases/mlst_db/`.

### 10. petG detection

[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) v2.17.0 searches the assembly with the X73 *petG* CDS (`petG_X73_NZ_CM001580.fasta`, 1581 bp, in `databases/LPS/`).
Only X73 is used as the reference, so petG-like paralogs are excluded by the identity and length thresholds.
A hit counts as *petG* when it has at least 95% identity, a genomic span greater than 1570 bp, and covers the X73 gene to within 3 bp of both ends.
The best such hit is then checked for function from its BLAST alignment, without translating the gene:

* `nonfunctional`: the net gap length is not a multiple of 3 (frameshift), a substituted codon becomes a stop codon before the X73 stop, or the start codon is lost
* `uncertain`: the reading frame is restored by compensating frameshifts, or there is an in-frame indel
* `functional`: none of the above
* `absent`: no hit passes the thresholds (treated like `nonfunctional` by the phenotype rules)

Hit sequences and the raw BLAST table (with alignment strings) are written to `<sample>/13_petG/`.
The status and the reason (for example `frameshift: net -1 bp (first -1 bp at ref 988)`) are reported in `10_report/10_Illumina_petG_status.tsv` and the phenotype report.

### 11. LPS subtype and phenotype reports

The subtype report (`10_Illumina_subtype_report.tsv`) lists every Snippy variant that matches the [LPS subtype database](databases/LPS/LPS_subtype_database_v3.txt) at the same reference sequence, position, reference allele and alternate allele.

The phenotype report (`10_Illumina_phenotype_report.tsv`) gives one predicted phenotype per sample.
It is produced by `bin/lps_phenotype.py`, which is shared with the ONT pipeline, from three inputs:

* `databases/LPS/phenotype_rules.tsv`: the phenotype definitions (the only place phenotypes are defined)
* `databases/LPS/LPS_subtype_database_v3.txt`: the curated subtypes, with an optional curated `EFFECT` per variant
* the sample's variants (all Snippy calls in `snps.tab` (including the majority-allele 0/1 calls)), the per-gene read coverage and the petG status

#### How a phenotype is assigned

Before these steps the sample gets a subtype status (see "Subtype status"), which can withhold the phenotype altogether.
After that each sample goes through four steps.

1. Every gene named in `phenotype_rules.tsv` gets one **gene state** from the sample's variants and read coverage (see "Gene states").
2. *petG* gets a **petG state** from the BLAST result (see "petG states").
3. The **lowest step** that contains an `inactivating` gene decides the step. If no gene is inactivating, the `full` row is used.
4. Among the rule rows for that step, the first one whose `PETG` and `REQUIRES` conditions are met gives the phenotype.

Worked examples with the L1 rules:

* *pcgB* has a frameshift and *petG* is functional. Step 3 is the lowest inactivated step, so the phenotype is `L1_P3`.
* *hptE* (step 1) and *pcgB* (step 3) are both inactivating. Step 1 is lower, so the phenotype is `L1_P6`, whatever *petG* is.
* *pcgB* has an in-frame indel (state `uncertain`) and *petG* is nonfunctional. If *pcgB* were functional the `full` row gives `L1_P2`. If it were inactivating, step 3 gives `L1_P4`. Both are reported as `L1_P2 | L1_P4`.

#### Curated variants

A variant is **curated** when it is listed in the [LPS subtype database](databases/LPS/LPS_subtype_database_v3.txt).
A sample variant matches a database row when the reference sequence, position, reference allele and alternate allele are all the same.
Curators can set an optional `EFFECT` on a database row:

* `inactivating`: the variant is known to abolish the gene's function
* `modifying`: the variant is known to change the structure the gene makes (for example an amino acid change that makes the enzyme act like a different one)
* `none`: the variant is known to have no effect
* blank: no curated effect, so the gene state is derived from the variant itself (stop gained, frameshift and so on)

A subtype row with `VARTYPE` = `reference` and no position is the type's reference sequence, for example L3-1 (Pm70).
It is never matched by a variant; it is reported when a sample has no variants in the locus (see "Subtype status").
Other rows without a position (for example L3-3 and L3-4) are subtypes whose causative variants are not known yet, so they cannot be matched either.

A variant that is not in the database is **uncurated**.
Uncurated variants can still make a gene `inactivating` or `uncertain` through the rules below, and the phenotype basis flags them as candidates for curation.

#### Subtype status

Every sample also gets a `SUBTYPE_STATUS` that says how far the subtype assignment could be taken.
It is decided before the phenotype and can stop a phenotype from being given.

| Status | When | Subtype | Phenotype |
|---|---|---|---|
| `assigned` | At least one variant matches the subtype database | The matched subtype(s) | From the rules (for L3 see "Types that need a database subtype") |
| `reference_identical` | No variants anywhere in the LPS locus, and every CDS of the reference is covered (see "Coverage") | The reference subtype (`VARTYPE = reference` in the database, for example L3-1), or `reference` if the database has none | From the `full` row |
| `uncurated_variants` | There are variants in the locus and none of them is in the database | `Not assigned` | L1: from the rules. L3: `Not predicted` |
| `not_assessed` | Variant calling produced no output for the sample (Snippy did not run for it (no `_depth.tsv` was produced)), or there are no variants but the locus is not fully covered or has no coverage data | `Not assigned` | `Not predicted` |
| `not_applicable` | The sample is untypeable | `Not assigned` | `Not predicted` |

`reference_identical` and `not_assessed` are separated by coverage.
A missing gene leaves no variant, so it looks the same as an unchanged one; without the coverage check a failed or partial run would be reported as identical to the reference.

`Not predicted` is written to `PREDICTED_PHENOTYPE` (and its description) instead of a blank whenever no phenotype is given: untypeable, no rules for the type, not assessed, a type that needs a database subtype and has none, or no matching rule row.
`Not assigned` is written to `SUBTYPES` and `SUBTYPE` instead of `NA` or a blank.
Both are filterable text, and `SUBTYPE_STATUS` says why.

#### Types that need a database subtype

L3 structures are not explained by inactivating mutations alone.
For example, curators report that P1059 (L3_P2) and P1662 (L3_P4) have chemically determined structures without a known inactivating mutation.
So a type with `CURATED_ONLY = yes` in `phenotype_rules.tsv` (currently L3) gives a phenotype only when the database supports it:

* `reference_identical`: the `full` row is used (L3_P1), because the sequence is the reference, whose structure was determined.
* `assigned`: the rules are run twice, once with all variants and once with only the curated ones.
  If both give the same phenotype it is reported, and any uncurated variants are flagged.
  If they differ, the phenotype is `Not predicted`.
* `uncurated_variants`: `Not predicted`, however the rules would score the variants.
* In every `Not predicted` case the basis says what the rules alone would give and lists the uncurated variants, so the isolate can be curated.

Examples (hypothetical isolates):

* L3-6 (*natC*, step 6) plus a new *hptE* stop (step 1). With all variants the rules give L3_P7, with only the curated ones L3_P2. The two differ, so the phenotype is `Not predicted`.
* L3-12 (*hptE*, step 1) plus a new *natC* stop (step 6). *hptE* wins either way, so L3_P7 is reported and the *natC* variant is flagged as uncurated.

L1 does not need a database subtype: X73 (L1_P1) and VP161 (L1_P2) are both structure-determined, so the L1 rules apply in every status except `not_assessed` and `not_applicable`.

#### Coverage

Both pipelines measure coverage the same way: `samtools depth -a` on the reads mapped to the type's LPS reference (Snippy's `snps.bam`).
A position is covered at a depth of at least `--phenotype_min_depth` (default 10), and a CDS is adequately covered when at least `--phenotype_min_gene_coverage` (default 90%) of its positions are covered.
The depth table is kept as `<sample>/8_snippy/<sample>_depth.tsv`.

It is used in two places, both because a deleted or unread gene leaves no variant and would otherwise read as intact:

* a rule gene below the threshold is `uncertain`
* a sample with no variants is `reference_identical` only if every CDS of the reference passes (not only the rule genes); otherwise it is `not_assessed`, and the basis names the genes that failed

#### Gene states

Each rule gene gets exactly one state, computed from all of the sample's variants inside its CDS (coordinates from the reference GenBank file).
If several states apply, the first one in this table wins.

| State | When | Effect on the phenotype |
|---|---|---|
| `inactivating` | A variant is curated `EFFECT=inactivating`, or snpEff reports `stop_gained` or `start_lost`, or an indel touches the start or stop codon, or the net indel length in the gene is not a multiple of 3 (net frameshift) | Selects the step of the gene |
| `uncertain` | The reading frame is intact but altered (an in-frame indel, or frameshifts that cancel out, such as L1-3 + L1-4 in *pcgB*, net +9 bp), or snpEff reports `stop_lost`, or less than `--phenotype_min_gene_coverage` (default 90%) of the CDS is covered by reads | Evaluated both ways, see "Uncertain calls" |
| `modifying` | A variant is curated `EFFECT=modifying`. One such variant is enough, and every `modifying` variant counts the same | None on its own. Only a `REQUIRES` condition can use it |
| `functional` | Anything else, including plain substitutions | None |

Variants curated `EFFECT=none` are ignored, which is how an `uncertain` call is resolved as functional once there is evidence.

#### Confidence gate for uncurated calls

Curated variants always count.
An uncurated variant call that is not confident cannot decide a gene state on its own.
A call is low-confidence when its allele fraction is below `--phenotype_min_confident_af` (default 0.9).
Confidence uses the allele fraction from Snippy's rescued 0/1 calls (`snps.rescued.tsv`); calls Snippy keeps as `1/1` are treated as confident.
If leaving the low-confidence calls out changes the gene state, the gene is `uncertain` and the basis names the call, for example `uncertain: hptE depends on a low-confidence call (1bp ins at 5155: AF=0.63, extends a GGG run)`.
The gate can only widen a call to `uncertain`; it never makes a call more definite.
No translation is done, so a stop codon between two compensating frameshifts cannot be seen; this is why such genes are `uncertain` rather than `functional`.

#### petG states

| petG status | Meaning | Value used by the rules |
|---|---|---|
| `functional` | Hit found, no frameshift, premature stop or lost start | `functional` |
| `nonfunctional` | Net frameshift, premature stop codon or lost start codon | `nonfunctional` |
| `absent` | No hit passes the thresholds | `nonfunctional` |
| `uncertain` | Reading frame restored by compensating frameshifts, or an in-frame indel | both `functional` and `nonfunctional` |
| `not assessed` | petG detection was skipped | both `functional` and `nonfunctional` |

*petG* only matters for types that have a rule row with `PETG` set to `functional` or `nonfunctional` (currently L1).
For types whose rows are all `any` (currently L3), the *petG* status is ignored.

#### Phenotype rules

`phenotype_rules.tsv` has one row per phenotype:

| Column | Meaning |
|---|---|
| `TYPE` | LPS outer core type (L1, L3, ...) |
| `STEP` | Position of the gene(s) in biosynthetic order, 1 = the residue closest to the inner core. `full` is the row used when no rule gene is inactivated |
| `GENES` | The interchangeable genes of that step, for example `pcgA,pcgB,pcgC,pcgD`. `-` for `full` rows |
| `PETG` | `functional`, `nonfunctional` or `any` (see "petG states") |
| `REQUIRES` | Optional extra condition, `-` if none (see below) |
| `PHENOTYPE` | The phenotype reported |
| `CURATED_ONLY` | `yes` or `no` (default `no`), the same on every row of a type. `yes` means a phenotype is only given when a database subtype supports it (see "Types that need a database subtype") |

The current file:

```
TYPE  STEP  GENES                PETG           REQUIRES  PHENOTYPE  CURATED_ONLY
L1    1     hptE                 any            -         L1_P6      no
L1    2     gatA                 any            -         L1_P5      no
L1    3     pcgA,pcgB,pcgC,pcgD  functional     -         L1_P3      no
L1    3     pcgA,pcgB,pcgC,pcgD  nonfunctional  -         L1_P4      no
L1    full  -                    functional     -         L1_P1      no
L1    full  -                    nonfunctional  -         L1_P2      no
L3    1     hptE                 any            -         L3_P7      yes
L3    2     gctC                 any            -         L3_P6      yes
L3    3     gatF                 any            -         L3_P5      yes
L3    4     gatG                 any            -         L3_P4      yes
L3    5     natB                 any            -         L3_P3      yes
L3    6     natC                 any            -         L3_P2      yes
L3    full  -                    any            -         L3_P1      yes
```

Types without rules (currently L2 and L4 to L9) report `Not predicted`.
Adding a type, a step or a petG split is a change to `phenotype_rules.tsv` only.
The file is validated at run time: gene names must be CDS in the type's reference, every type needs a `full` row, a gene can appear at only one step, `CURATED_ONLY` must be the same on every row of a type, and duplicate rows are rejected.

##### REQUIRES

`REQUIRES` makes a row more specific than the plain row for the same step.

* Format: `gene:state`, with several conditions joined by `;` (all must hold), for example `hptE:modifying;gatA:functional`.
* The state must be `inactivating`, `modifying` or `functional`.
* The gene must be a CDS in the type's reference. It does not need to be a step gene.
* A condition is met when the gene's computed state equals the stated state (see "Gene states"). `REQUIRES` cannot name individual variants. For `gene:modifying`, a single curated `modifying` variant in the gene is enough, and several such variants have the same effect as one. Any `inactivating` or `uncertain` signal in the same gene overrides `modifying`, so the condition is then not met.
* Rows with a `REQUIRES` are tried first, then rows without. The first row whose `PETG` and `REQUIRES` both match is used. If no `REQUIRES` row matches, the plain row applies.
* `REQUIRES` also works on `full` rows.
* If several rows match, the one earlier in the file wins, so avoid overlapping rows.

Example (hypothetical):

```
L1  3  pcgA,pcgB,pcgC,pcgD  functional  hptE:modifying  L1_Px
```

This applies only when the lowest inactivated step is 3, *petG* is functional and *hptE* carries a curated modifying variant.
Otherwise the plain step-3 row applies.

#### Uncertain calls

An `uncertain` gene (and an uncertain or unassessed *petG*) is evaluated both ways, and every distinct phenotype is reported, for example `L1_P1 | L1_P3`.
For a gene, the rules are run twice, once treating it as `functional` and once as `inactivating`, and `REQUIRES` is checked against those assumed states.
So a condition like `hptE:modifying` is never met while *hptE* is `uncertain`, and `hptE:uncertain` can never be met at all.

#### Phenotype basis

`PHENOTYPE_BASIS` states how much to trust each call:

* `curated`: every variant that changed a gene state is in the subtype database
* `uncertain: ...`: the call depends on an uncertain gene or petG; the text names the gene, the reason and the phenotype for each alternative
* `uncurated: ...`: a variant that is not in the subtype database changed a gene state (the variants are listed, and `SUBTYPES` includes `novel`); these are candidates for curation
* `mixed call: ...` (suffix): a curated variant that changed a gene state has an allele fraction below 0.9 (rescued 0/1 calls with alt fraction < 0.9 (`snps.rescued.tsv`))
* `uncurated variants: ...`: the sample has variants in the locus that are not in the subtype database and did not change a gene state (listed with their gene, up to 8)
* `not predicted: ...; the rules alone would give ...`: the type needs a database subtype and none supports the call (see "Types that need a database subtype")
* `not assessed: ...`: the subtype status is `not_assessed`; the text gives the reason
* `no variants in the locus (identical to the reference)`: the subtype status is `reference_identical`
* `no phenotype rules for type`, `untypeable`

`LPS_GENE_STATES` lists only the genes that are not plainly functional, as `gene:state` separated by `;`, for example `gatG:inactivating;natB:inactivating`.

### 12. Genome annotation using Bakta

[Bakta](https://github.com/oschwengers/bakta) v1.12.0 annotates the genome assemblies. The default database is v6.0 (2025-02-24), available from [Zenodo record 10.5281/zenodo.14916843](https://zenodo.org/records/14916843).

Bakta is run with `--proteins databases/trusted_bakta_proteins.fasta`, a set of trusted proteins used to prioritise annotations. This file is built with [`bin/bakta_proteins_format.py`](bin/bakta_proteins_format.py), which combines a type-strain GenBank file and/or additional protein FASTA files into the combined FASTA format Bakta expects for `--proteins`:

```bash
python bin/bakta_proteins_format.py type_strain.gbk additional_proteins.fasta -o databases/trusted_bakta_proteins.fasta
```

### 13. Antimicrobial resistance genes

[AMRFinderPlus](https://github.com/ncbi/amr) v4.0.23 identifies AMR genes in the assemblies. The tested database version is 2025-03-25.1.

---

## Step-by-step user guide

### 1. Clone the repository

```bash
git clone https://github.com/julianzaugg/LPS_typing_Illumina.git
cd LPS_typing_Illumina
```

The repository includes three key files:

**a) `nextflow.config`** — default parameters and container profiles. Nextflow reads this automatically from the working directory. Select a container engine with `-profile singularity` or `-profile apptainer`. Add `-profile slurm` on SLURM clusters. Container images are pulled automatically and cached in `./singularity/` by default.

**b) `main.nf`** — the pipeline code. Not user-modifiable.

**c) `nextflow.sh`** — an optional SLURM submission script template. Edit it to set your account, partition, samplesheet path, and output directory, then submit with `sbatch nextflow.sh`.

### 2. Obtain databases

See the [Database setup](#database-setup) section below.

### 3. Prepare a samplesheet

The samplesheet is a comma-separated file listing each sample and its paired FASTQ files. Paths are resolved relative to the directory where you run `nextflow run` (the launch directory).

```bash
mkdir samplesheet
# Create/edit samplesheet/samples.csv:
```

```
sample_id,short_fastq_1,short_fastq_2
PM3034,fastq/PM3034_R1.fastq.gz,fastq/PM3034_R2.fastq.gz
PM3065,fastq/PM3065_R1.fastq.gz,fastq/PM3065_R2.fastq.gz
```

### 4. Run the pipeline

```bash
# Local execution with Apptainer:
nextflow run main.nf -profile apptainer \
  --samplesheet samplesheet/samples.csv \
  --outdir results

# SLURM cluster with Apptainer:
nextflow run main.nf -profile apptainer,slurm \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  --slurm_account YOUR_ACCOUNT

# Use -resume to restart from the last successful step:
nextflow run main.nf -profile apptainer,slurm \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  --slurm_account YOUR_ACCOUNT \
  -resume
```

---

## Database setup

The LPS reference data and Kaptive database are **included in the repository** under `databases/` and require no additional setup. Large third-party databases must be obtained separately; you can either download them automatically using pipeline flags, or provide pre-downloaded copies.

| Database | Required for | Default path | Tested version | Source | Approx. size | How to obtain |
|----------|-------------|--------------|---------------|--------|-------------|---------------|
| LPS references | Steps 7, 8, 10, 11 | `databases/LPS/` | Included in repo | This repository | < 1 MB | Included — no action needed |
| Kaptive3 LPS DB | Step 7 | `databases/kaptive3_LPS_db_v1/9lps.gbk` | Included in repo | This repository | < 200 KB | Included — no action needed |
| Sylph GTDB + Fungi | Step 6 | `databases/sylph/` | GTDB R232, Fungi RefSeq 2025-10-11 | [Sylph pre-built databases](https://sylph-docs.github.io/pre%E2%80%90built-databases/) | ~ 25 GB | Add `--skip_download_sylph_db false`, or download manually (see below) |
| CheckM | Step 5 | `databases/checkm_data_2015_01_16` | 2015_01_16 | [CheckM installation docs](https://github.com/Ecogenomics/CheckM/wiki/Installation) | ~ 1.4 GB | Add `--download_checkm_db`, or download manually (see below) |
| Bakta | Step 12 | `databases/bakta/db` | v6.0 (2025-02-24) | [Zenodo 10.5281/zenodo.14916843](https://zenodo.org/records/14916843) | ~ 65 GB | Add `--download_bakta_db`, or download manually (see below) |
| AMRFinderPlus | Step 13 | `databases/amrfinderplus/2025-03-25.1` | 2025-03-25.1 | [NCBI AMRFinderPlus wiki](https://github.com/ncbi/amr/wiki/AMRFinderPlus-database) | ~ 300 MB | Add `--download_amrfinder_db`, or download manually (see below) |

> **Reproducibility note:** The `--download_*` flags (and `--skip_download_sylph_db false`) retrieve the latest available database version, which may differ from the tested versions listed above and could affect results. For reproducible analyses, use pinned database versions and specify their paths explicitly with `--sylph_db`, `--sylph_metadata`, `--checkm_db`, `--bakta_db`, and `--amrfinder_db`.

### Downloading databases via pipeline flags

Add the relevant flag(s) on the first run. The database will be downloaded into `databases/` and reused automatically on subsequent runs (omit the flag after the first download).

```bash
# Download Sylph databases (GTDB R232 + Fungi RefSeq, ~25 GB):
nextflow run main.nf -profile apptainer --samplesheet samplesheet/samples.csv \
  --outdir results --skip_download_sylph_db false

# Download CheckM database (~1.4 GB):
nextflow run main.nf -profile apptainer --samplesheet samplesheet/samples.csv \
  --outdir results --download_checkm_db

# Download Bakta database (~65 GB, slow):
nextflow run main.nf -profile apptainer --samplesheet samplesheet/samples.csv \
  --outdir results --download_bakta_db

# Download AMRFinderPlus database (~300 MB):
nextflow run main.nf -profile apptainer --samplesheet samplesheet/samples.csv \
  --outdir results --download_amrfinder_db
```

> **Zenodo connection errors:** Downloads from Zenodo (Sylph-tax metadata, Bakta) can occasionally fail on HPC clusters with `Connection reset by peer`. The pipeline will automatically retry up to 5 times with 30-second gaps. If downloads still fail, download the files manually (see below) and use `--sylph_metadata` / `--bakta_db` to point the pipeline to the local copies.

### Manual database download

If you prefer to download databases yourself (e.g. for speed or reproducibility), download them to any location and point the pipeline to them:

```bash
# Sylph databases (GTDB R232 + Fungi RefSeq) and matching sylph-tax metadata:
mkdir -p /path/to/databases/sylph
cd /path/to/databases/sylph
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r232-c200-dbv1.syldb
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/fungi-refseq-2025-10-11-c200-dbv1.syldb
wget https://zenodo.org/records/19646381/files/gtdb_r232_metadata.tsv.gz
wget https://zenodo.org/records/17330476/files/fungi_refseq_2025-10-11_metadata.tsv.gz
# Then run with: --sylph_db "/path/to/databases/sylph/*.syldb" --sylph_metadata "/path/to/databases/sylph/*.tsv.gz"

# CheckM (tested version checkm_data_2015_01_16):
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
mkdir -p /path/to/databases/checkm_data_2015_01_16
tar -xvzf checkm_data_2015_01_16.tar.gz -C /path/to/databases/checkm_data_2015_01_16/
# Then run with: --checkm_db /path/to/databases/checkm_data_2015_01_16

# Bakta (tested version v6.0 from Zenodo):
wget https://zenodo.org/records/14916843/files/db.tar.gz
tar -xvzf db.tar.gz -C /path/to/databases/bakta/
# Then run with: --bakta_db /path/to/databases/bakta/db

# AMRFinderPlus (tested version 2025-03-25.1):
# Use amrfinder_update within the container, or specify an existing database directory
# Then run with: --amrfinder_db /path/to/databases/amrfinderplus/2025-03-25.1
```

---

## Optional parameters

### 1. Read trimming

* `--skip_fastp`: skip the read trimming step (default=false). Not recommended.

### 2. FastQC reads quality metrics

* `--skip_fastqc`: skip the FastQC step (default=false)
* `--skip_summary_fastqc`: skip the MultiQC summary step (default=false)

### 3. Genome assembly

* `--shovill_threads`: number of threads for Shovill (default=4)
* `--shovill_args`: additional Shovill parameters (default="", example: `"--minlen 200 --mincov 10"`)
* `--genome_size`: estimated genome size (default="2.3M")

### 4. Assembly quality assessment with QUAST

* `--skip_quast`: skip the QUAST step (default=false)
* `--quast_threads`: number of threads for QUAST (default=2)

### 5. Assembly quality assessment with CheckM

* `--skip_checkm`: skip the CheckM step (default=false)
* `--checkm_db`: path to the CheckM database folder (default=`databases/checkm_data_2015_01_16`)
* `--download_checkm_db`: download the CheckM database automatically (default=false)

### 6. Sylph taxonomy classification

* `--skip_sylph`: skip the Sylph classification step (default=false)
* `--skip_download_sylph_db`: skip downloading the Sylph databases (default=true — assumes databases are already present at `--sylph_db`)
* `--sylph_db_gtdb_file` and `--sylph_db_fungal_file`: URLs for Sylph GTDB and Fungi RefSeq database files to download
* `--sylph_tax_gtdb_metadata` and `--sylph_tax_fungal_metadata`: URLs for Sylph-tax metadata files to download (must match the database files)
* `--sylph_db`: glob pattern for pre-downloaded Sylph database files (default=`databases/sylph/*.syldb`)
* `--sylph_metadata`: glob pattern for pre-downloaded Sylph-tax metadata files (default=`databases/sylph/*.tsv.gz`)
* `--sylph_threads`: number of threads for the Sylph classification step (default=6)

### 7. LPS typing using Kaptive

* `--skip_kaptive3`: skip the Kaptive typing step (default=false). Note: skipping this also skips variant calling (Snippy) and petG detection.
* `--kaptive_db_9lps`: path to the Kaptive database file (default=`databases/kaptive3_LPS_db_v1/9lps.gbk`)

### 8. Variant calling using Snippy

* `--skip_snippy`: skip the Snippy variant calling step (default=false)
* `--snippy_threads`: number of threads for Snippy (default=6)
* `--snippy_args`: additional Snippy parameters (default=""). Do not pass `--minqual`, `--mincov` or `--minfrac` here; use the options below.
* `--snippy_minqual`: minimum QUAL for a call (default=100, as Snippy)
* `--snippy_mincov`: minimum read depth for a call (default=10, as Snippy)
* `--snippy_minfrac`: minimum alt fraction for `1/1` calls (default=0, as Snippy)
* `--snippy_het_min_af`: minimum alt fraction for freebayes `0/1` calls to be kept in `snps.tab`, or `false` to keep Snippy's `1/1`-only filter (default=0.5). See [Majority-allele 0/1 calls](#majority-allele-01-calls).
* `--reference_LPS_directory`: path to the LPS reference directory containing `reference_LPS.txt`, `LPS_subtype_database_v3.txt`, `phenotype_rules.tsv`, `petG_X73_NZ_CM001580.fasta`, and FASTA/GB files for all LPS types (default=`databases/LPS`)

### 9. MLST typing

* `--skip_mlst`: skip the MLST typing step (default=false)
* `--mlst_scheme`: MLST typing scheme (default="pmultocida_2")
* `--mlst_datadir`: path to a directory of custom MLST scheme files (default="" = use bundled database). The directory must contain a subdirectory named after the scheme (matching `--mlst_scheme`) holding allele FASTA files (`*.fas` from PubMLST, or `*.tfa`) and a `<scheme>.txt` profile. File renaming and BLAST indexing are handled automatically by the pipeline.

### 10. petG detection

* `--skip_petg`: skip petG detection with BLAST (default=false)
* `--petg_threads`: number of threads for the BLAST step (default=2)
* `--petg_min_length`: minimum genomic hit span in bp for a petG hit; hits must exceed this value (default=1570)
* `--petg_min_identity`: minimum percent identity for a petG hit (default=95)

### 11. Report

The subtype and phenotype reports use `LPS_subtype_database_v3.txt`, `phenotype_rules.tsv` and `phenotype_lookup.tsv` from `--reference_LPS_directory` (see [LPS subtype and phenotype reports](#11-lps-subtype-and-phenotype-reports)).

* `--subtype_database`: subtype database (default="" = `LPS_subtype_database_v3.txt` in `--reference_LPS_directory`)
* `--phenotype_rules`: phenotype rules file (default="" = `phenotype_rules.tsv` in `--reference_LPS_directory`)
* `--phenotype_min_gene_coverage`: minimum covered fraction of a rule gene's CDS; below it the gene is `uncertain` (default=0.9)
* `--phenotype_min_depth`: minimum read depth for a position to count as covered, measured with `samtools depth` on the Snippy BAM (default=10)
* `--phenotype_min_confident_af`: uncurated variant calls below this allele fraction are low-confidence (default=0.9)
* `--phenotype_homopolymer_min_run`: uncurated indels of a base in a reference run of at least this length are low-confidence (default=0, i.e. off)

The pipeline also generates a self-contained HTML summary report from the combined `10_report` outputs when Kaptive and Snippy are enabled.

* `--skip_html_report`: skip generation of the combined HTML report (default=false)

### 12. Genome annotation using Bakta

* `--skip_bakta`: skip genome annotation (default=false)
* `--bakta_threads`: number of threads for Bakta (default=8)
* `--bakta_db`: path to the Bakta database folder (default=`databases/bakta/db`)
* `--download_bakta_db`: download the latest Bakta database automatically (default=false)
* `--bakta_args`: additional Bakta parameters (default="")

### 13. AMR gene identification using AMRFinderPlus

* `--skip_amrfinder`: skip AMR gene identification (default=false)
* `--amrfinder_db`: path to the AMRFinderPlus database folder (default=`databases/amrfinderplus/2025-03-25.1`)
* `--download_amrfinder_db`: download the latest AMRFinderPlus database automatically (default=false)
* `--amrfinder_args`: additional AMRFinderPlus parameters (default="")

---

## Structure of the output folders

The pipeline creates a folder per sample (named by `sample_id`) inside `--outdir`, plus a combined `10_report` folder.

Each sample folder contains:

* **1_trimming:** Trimmed paired-end FASTQ files (`sample_id_R1_trimmed.fastq.gz`, `sample_id_R2_trimmed.fastq.gz`).
* **2_fastqc:** FastQC results — HTML report and ZIP archive for each read set.
* **3_assembly:** Shovill assembly output. See [Shovill output docs](https://github.com/tseemann/shovill?tab=readme-ov-file#output-files).
  * Final assembly in FASTA format (`sample_id_contigs.fa`)
  * SPAdes assembly graph (`sample_id_contigs.gfa`)
* **4_quast:** QUAST report (`sample_id_report.tsv`).
* **5_checkm:** CheckM lineage workflow results (`sample_id_checkm_lineage_wf_results.tsv`).
* **6_sylph:** Sylph taxonomy classification results, see [Sylph output format](https://sylph-docs.github.io/Output-format/) and [sylph-tax output format](https://sylph-docs.github.io/sylph-tax-output-format/).
  * Sylph list of genomes detected (`sample_id_sylph_profile.tsv`)
  * Sylph-tax combined taxonomic abundances (`sample_id_merged_taxonomic_abundance.tsv`)
  * Sylph-tax combined sequence abundances (`sample_id_merged_sequence_abundance.tsv`)
* **7_kaptive_v3:** Kaptive output files. See [Kaptive output docs](https://kaptive.readthedocs.io/en/latest/Outputs.html).
  * LPS type results (`sample_id_kaptive_results.tsv`)
  * LPS sequence in FASTA format (`sample_id_kaptive_results.fna`)
* **8_snippy:** Snippy mapping and variant calling results. See [Snippy output docs](https://github.com/tseemann/snippy?tab=readme-ov-file#output-files).
  * BAM alignment file (`sample_id_snps.bam` and `.bai` index)
  * Unfiltered variants VCF (`sample_id_snps.raw.vcf`)
  * Filtered variants VCF, Snippy's `1/1` calls (`sample_id_snps.filt.vcf`)
  * Variant summary in tabular format, including majority-allele `0/1` calls (`sample_id_snps.tab`)
  * High-impact variants (frameshift, stop_gained) (`sample_id_snps.high_impact.tab`)
  * Annotated VCF behind `sample_id_snps.tab` (`sample_id_snps.typing.vcf`)
  * Snippy's own `1/1`-only table (`sample_id_snps.snippy.tab`)
  * Calls added over Snippy's filter, with allele fraction (`sample_id_snps.rescued.tsv`)
  * Snippy log (`sample_id_snippy.log`)
* **9_mlst:** MLST typing result (`sample_id_mlst_pmultocida_rirdc.csv`).
* **13_petG:** petG BLAST output files.
  * BLAST tabular output for all hits (`sample_id_petG_blast.tsv`)
  * BLAST tabular output for accepted hits (`sample_id_petG_blast.filtered.tsv`)
  * FASTA sequences for accepted hits (`sample_id_petG_hits.fasta`)
  * petG presence summary (`sample_id_petG_summary.tsv`)
* **11_bakta:** Bakta annotation output. See [Bakta output docs](https://github.com/oschwengers/bakta?tab=readme-ov-file#output).
  * Annotations in GenBank format (`sample_id_bakta.gbff`)
  * Inference metrics (`sample_id_bakta.inference.tsv`)
  * Annotation summary (`sample_id_bakta.txt`)
* **12_amrfinder:** AMRFinderPlus results (`sample_id_amrfinder.tsv`). See [output format](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#output-format).

The `10_report` folder contains combined results across all samples:

* MultiQC report (`2_Illumina_multiqc_report.html`) and general statistics (`2_Illumina_multiqc_general_stats.txt`)
* Shovill assembly statistics: coverage, contig count, assembly size (`3_Illumina_shovill_stats.tsv`)
* Combined QUAST report (`4_Illumina_quast_report.tsv`)
* Combined CheckM results (`5_Illumina_checkm_lineage_wf_results.tsv`)
* Sylph taxonomy results:
  * Abundance of *P. multocida* reads, and information on the most abundant species (if not *P. multocida*): `6_Illumina_sylph_summary.tsv`
* Kaptive results (`7_Illumina_kaptive_results.tsv`)
* Snippy variant results:
  * All variants (`8_Illumina_snippy_snps.tsv`)
  * High-impact variants only (`8_Illumina_snippy_snps.high_impact.tsv`)
* MLST results (`9_Illumina_mlst.csv`)
* Self-contained HTML summary report (`LPS_typing_report.html`) with sample-level LPS calls, QC summaries, phenotype information, and variant/locus visualisations. This report is generated when `--skip_html_report false`, `--skip_kaptive3 false`, and `--skip_snippy false`. The per-LPS-type lollipop plots show the mutations observed across the run (lollipop height = number of genomes carrying each mutation); the gene track below each plot is drawn as directional arrows, so gene orientation on the locus is visible at a glance. Hovering a lollipop or gene shows its details, a gene-colour legend accompanies each plot, and the full mutation list is available in a collapsible table beneath it. For publication-quality PDF/SVG versions of these plots see [Publication figures](#publication-figures).
* Subtype report (`10_Illumina_subtype_report.tsv`). Column descriptions:
  * **SAMPLE**: sample identifier
  * **MLST**: MLST sequence type
  * **TYPE**: LPS type assigned by Kaptive (or `untypeable`)
  * **SUBTYPE**: LPS subtype from the subtype database; the reference subtype for a `reference_identical` sample, or `Not assigned`
  * **SUBTYPE_STATUS**: `assigned`, `reference_identical`, `uncurated_variants`, `not_assessed` or `not_applicable` (see "Subtype status")
  * **VARTYPE**: description of the variant
  * **ISOLATE_DATABASE**: reference isolate from the subtype database
  * **CHROM**: reference sequence name for the LPS locus
  * **POS**: variant position in the reference
  * **REF**: reference allele
  * **ALT**: alternate allele identified in the sample
  * **GENE**: gene containing the variant (curated label from the subtype database)
  * **EFFECT**: curated effect of the variant (`none`, `modifying`, `inactivating`, or blank = derived by the gene-state rules)
  * **PETG_STATUS**: petG status (`functional`, `nonfunctional`, `uncertain`, `absent`, or `not assessed`)
  * **NOTE**: subtype database note for the matched variant, when available
* Phenotype report (`10_Illumina_phenotype_report.tsv`), one row per sample. Column descriptions:
  * **SAMPLE**, **MLST**, **TYPE**: as in the subtype report
  * **SUBTYPES**: matched subtypes, separated by `;` (`novel` is added when an uncurated variant changed a gene state); the reference subtype for a `reference_identical` sample; `Not assigned` otherwise
  * **SUBTYPE_STATUS**: as in the subtype report
  * **PETG_STATUS**, **PETG_DETAIL**: petG status and the reason or hit location
  * **LPS_GENE_STATES**: rule genes that are not plainly functional, e.g. `gatG:inactivating;natB:inactivating`
  * **DETERMINING_GENE**: the inactivated gene of the lowest step (`none` when no rule gene is inactivated)
  * **PREDICTED_PHENOTYPE**: predicted phenotype; alternatives are separated by ` | ` when the call is uncertain; `Not predicted` when no phenotype is given
  * **PREDICTED_PHENOTYPE_DESCRIPTION**: description(s) from `phenotype_lookup.tsv`, or `Not predicted`
  * **PHENOTYPE_BASIS**: `curated`, `uncertain: ...`, `uncurated: ...`, `uncurated variants: ...`, `not predicted: ...`, `not assessed: ...`, `no phenotype rules for type` or `untypeable`, plus an optional `mixed call: ...` suffix
* petG status per sample (`10_Illumina_petG_status.tsv`)
* AMRFinderPlus combined results (`12_Illumina_amrfinder.tsv`)

---

## Publication figures

`bin/plot_lps_mutations.R` draws the per-LPS-type mutation lollipops as standalone
vector files for use in papers.
Each figure shows the same data as the HTML report (lollipop height = number of genomes carrying the mutation)
above a directional gene-arrow map of the reference locus, with gene direction taken from the `complement(...)`
annotation in the LPS database GenBank files.

This script is run by hand against a finished results directory.
It is not part of the Nextflow workflow and does not use a container.

Requirements (R, once):

```bash
Rscript -e 'install.packages(c("ggplot2", "dplyr", "scales", "svglite"))'
```

Alternatively pass `--install-deps` and the script installs whatever is missing.

Basic use, writing one PDF and one SVG per LPS type into `figures/`:

```bash
Rscript bin/plot_lps_mutations.R \
  --report-dir results/10_report \
  --lps-db-dir databases/LPS \
  --outdir figures
```

Common options:

```bash
# a single LPS type, colour-blind-safe palette, also write the plotted variants as TSV
Rscript bin/plot_lps_mutations.R --report-dir results/10_report --lps-db-dir databases/LPS \
  --outdir figures --types L3 --palette okabe-ito --write-table

# one genome only (every lollipop is then height 1: a per-genome mutation map)
Rscript bin/plot_lps_mutations.R --report-dir results/10_report --lps-db-dir databases/LPS \
  --outdir figures --samples PM1205

# combine an Illumina and an ONT run in one figure (genomes present in both are counted once)
Rscript bin/plot_lps_mutations.R --lps-db-dir databases/LPS --outdir figures \
  --subtype-report illumina_run/10_report/10_Illumina_subtype_report.tsv,ont_run/10_report/10_ONT_subtype_report.tsv
```

`--help` lists every option, including `--min-genomes`, `--exclude-samples`, `--legend`,
`--gene-colors`, `--font`, `--width`/`--height`, `--prefix` and `--formats`.
Gene colours default to `gene_colors.tsv` in the LPS database, so the figures match the HTML report.

---

## Advanced use

### Running in assembly-only mode for other organisms

The default parameters are optimised for *Pasteurella multocida*. LPS typing, variant calling, and petG detection are species-specific. To use the pipeline for genome assembly, QC, taxonomy, and MLST for another species, skip the *P. multocida*-specific steps:

```bash
nextflow run main.nf -profile apptainer \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  --genome_size 2.5M \
  --mlst_scheme your_scheme \
  --skip_kaptive3 \
  --skip_bakta \
  --skip_amrfinder
```

Note: `--skip_kaptive3` automatically skips Snippy variant calling and petG detection.

### Using a custom cluster configuration

For clusters with non-standard SLURM partitions, memory limits, or other requirements, create a custom config file and pass it with `-c`:

```bash
# my_cluster.config
process {
    executor = 'slurm'
    clusterOptions = '--account=my_account --partition=high_mem'
    time = '12h'
    withLabel: high_memory { memory = 512.GB }
}
```

```bash
nextflow run main.nf -profile apptainer \
  --samplesheet samplesheet/samples.csv \
  --outdir results \
  -c my_cluster.config
```

---

## Acknowledgements / citations / credits

Please cite the following tools when using this pipeline:

- [fastp](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://github.com/MultiQC/MultiQC)
- [Shovill](https://github.com/tseemann/shovill) / [SPAdes](https://github.com/ablab/spades)
- [QUAST](https://quast.sourceforge.net/quast.html)
- [CheckM](https://github.com/Ecogenomics/CheckM)
- [Sylph](https://github.com/bluenote-1577/sylph)
- [sylph-tax](https://github.com/bluenote-1577/sylph-tax)
- [Kaptive](https://kaptive.readthedocs.io/en/latest/)
- [Snippy](https://github.com/tseemann/snippy) / [SnpEff](https://pcingola.github.io/SnpEff/)
- [mlst](https://github.com/tseemann/mlst)
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [Bakta](https://github.com/oschwengers/bakta)
- [AMRFinderPlus](https://github.com/ncbi/amr)

Pipeline developed by Valentine Murigneux and Julian Zaugg.
