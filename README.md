# longread-cohort

Long-read sequencing (Oxford Nanopore, PacBio) captures things short-read
cancer sequencing cannot: large and complex structural variants, which parental
copy of a chromosome a variant sits on (phasing), and DNA methylation read
directly from the same molecules. Open-source pipelines now call all of these
per sample, but nothing turns the per-sample output into cohort-level,
annotated, clinically joined results. This repository is the start of that layer.

A cohort-level interpretation and annotation layer for long-read somatic cancer
sequencing. Per-sample variant calling is left to existing pipelines
(nf-core/pacsomatic, LRSomatic); this project aggregates long-read structural
variants (SVs), phasing and methylation across a cohort, reconciles them against
short-read calls on the same tumours, annotates against cancer knowledge bases
(COSMIC now; CIViC, OncoKB planned), and joins to clinical data. Research use only.

Status (2026-09-30): first dataset loaded and inventoried; SV reconciliation
reproduces the source paper; first annotation pass done on 43 tumours.

## Data

Long-Read POG (BC Cancer; O'Neill et al., *Cell Genomics* 2024,
[PMC11605692](https://pmc.ncbi.nlm.nih.gov/articles/PMC11605692/)): 189 tumours
from 181 patients, Oxford Nanopore PromethION, GRCh38. Processed data from
https://www.bcgsc.ca/downloads/nanopore_pog/ (599 of 603 files, 1.3 GB; the
short-read MAVIS and Mutect2 archives, 1.2 GB; and, from the 17.3 GB
copy-number archive, the Ploidetect segment files for the 43 nanopore-normal
tumours only, stream-extracted without keeping the archive). 43 tumours have a nanopore-sequenced blood normal and
therefore somatic long-read SV calls (nanomonsv, SAVANA); all 189 have Illumina
tumour/normal calls. Gene model: GENCODE v50. Cancer genes: COSMIC Cancer Gene
Census as shipped with the dataset.

Full file-by-file description, sample-ID mapping and chromosome-naming table:
`results/data_inventory.md` (generated locally; `results/` is not committed).

## Results

### 1. Long-read vs short-read SV reconciliation reproduces the paper

The dataset ships a MAVIS table merging short-read (delly, manta) and long-read
(nanomonsv, SAVANA) somatic SVs for the 43 tumours. Applying the authors' own
filters (no DGV overlap, >50 bp, autosomes, Illumina high-quality flag) gives
the "consistent between platforms" counts the paper reports for Figure 2B:

| SV type | Paper | This repo |
|---|---|---|
| deletions | 3,358 (37.6%) | 3,355 (37.6%) |
| duplications | 1,919 (54.1%) | 1,919 (54.1%) |
| inversions | 1,943 (57.1%) | 1,943 (57.2%) |
| insertions | 7 (<1%) | 7 (0.1%) |

### 2. A `nanopore_rescued` category

The paper's filter discards any merged call whose *Illumina* quality flag is
low, even when nanopore independently supports it. We keep those as
`nanopore_rescued`: 2,724 calls (753 deletions, 769 duplications, 569
inversions, 587 translocations, 46 insertions), 23% as many as the 9,137
concordant calls. The "both" counts above are unchanged by this.

### 3. Nanopore-supported SVs in cancer genes

12,449 nanopore-only + rescued SVs across 43 tumours give 904 focal hits
(events ≤ 5 Mb) in 355 COSMIC Census genes: 437 breakpoints in introns, 409
genes wholly inside an event, 58 breakpoints in exons. Every tumour has at
least one. A further 7,865 hits lie inside > 5 Mb events and are held back
for copy-number analysis.

- Summary and five hits read one by one: [docs/cancer_gene_sv_summary.md](docs/cancer_gene_sv_summary.md)
- Per-sample counts: [docs/cancer_gene_sv_per_sample.tsv](docs/cancer_gene_sv_per_sample.tsv)
- Recurrent genes (top 50): [docs/cancer_gene_sv_recurrent_genes.tsv](docs/cancer_gene_sv_recurrent_genes.tsv)
- Chart: [docs/cancer_gene_sv_top_genes.png](docs/cancer_gene_sv_top_genes.png)

![top genes](docs/cancer_gene_sv_top_genes.png)

**Length bias.** Ranked by number of samples, the list is led by very long
genes — MUC4 (tandem-repeat exon insertions), LRP1B (1.9 Mb), CSMD3, PTPRD,
GPC5, ROBO2, CTNND2, FHIT (FRA3B fragile site), CNTNAP2 (2.3 Mb) — against a
median protein-coding gene of 29 kb. Normalised by gene length the entries
that stand out are ERBB2, IKZF3 and RARA (the 17q12 amplicon in breast
cancers), MYC, NF1 and SND1. A background-corrected recurrence test is the
next step; the current table is descriptive.

### 4. MSI tumours drive the biggest platform disagreement

Five tumours have < 12% platform agreement. They share no sequencing-QC
cause (nanopore depth 22–75×, normal N50 and error rate). Two — POG884
(gastro-oesophageal, TMB 172) and POG986 (lung squamous, TMB 125) — are the
two MSI tumours among the 43; 93–95% of their nanopore-only calls are small
insertions from SAVANA, consistent with microsatellite indels that short-read
SV callers do not report. The other three (POG530, POG044, POG1000) are
SV-quiet genomes (TMB 1.3–2.3) where a small number of nanopore-only calls
dominates a small denominator. All five are flagged `low_agreement` in every
per-sample table.

### 5. ERBB2 in POG137: nanopore adds structure, not detection

POG137 (breast IDC) carries an ERBB2 amplicon at 187–217 copies across the
gene body in short-read Ploidetect, with the highest ERBB2 expression in the
cohort (2,508 TPM; median 46). Short reads called the amplification and two
high-quality rearrangement breakpoints inside ERBB2 plus seven more within
0.5 Mb, all also seen by nanopore. Nanopore (SAVANA) adds five further
breakpoints inside the gene — three nested inversions of 21–125 kb, a
translocation to chr5, a 9.3 Mb inversion — and these coincide with
copy-number step boundaries in the short-read depth signal. The short-read
data did not miss the event; nanopore resolved sub-structure of an amplicon
short reads had already detected. Read-level review remains the proper
confirmation. Detail in
[docs/cancer_gene_sv_summary.md §5b](docs/cancer_gene_sv_summary.md).

## Caveats

- **Copy number not yet integrated.** Ploidetect segments are local for the
  43 nanopore-normal tumours (not the other 146), but only the ERBB2 locus
  has been checked by hand; the 7,865 Census hits inside > 5 Mb events are
  still unclassified.
- **Recurrence is length-biased** (section 3); no background model yet.
- **`whole_gene` conflates** loss (deletion), gain (duplication) and intact
  relocation (inversion).
- **Insertion sizes** are not in the main MAVIS table (breakpoints coincide);
  insertions are exempt from the > 50 bp filter here, so Illumina-only and
  nanopore-only insertion counts differ from the paper's, which used the
  separate `_ins` table.
- **Intronic nanopore-only insertions** are mostly repeat/mobile-element
  insertions; they inflate hits in long genes.
- **Tumour-type "driver" matching** is a keyword lookup against the Census
  free-text field; "no" means not matched, not not-a-driver (NF1 in sarcoma
  reads "no").
- **Platform category comes from the MAVIS `tools` string**; breakpoints were
  not re-matched, and the MAVIS merge tolerance is not stated in the paper
  text retrieved (STAR Methods not accessible via the APIs used).
- **Sample-ID mapping** between VCF names (`POG117-1`) and clinical names
  (`POG117-OCT-2`) was derived from library IDs in the SAVANA VCF headers;
  the paper and repository do not document it.
- **Sex chromosomes** are excluded by the paper's filter; this implementation
  drops events with either end on X/Y, the authors' R code only events with
  both ends on the same sex chromosome (likely the 3-deletion difference).
- **Low-agreement samples** are included and flagged, not removed.
- **Not yet validated** against a truth set; HCC1395 / SEQC2 is planned.

## Next

1. **Copy-number integration.** Classify the 7,865 Census hits inside > 5 Mb
   events against the Ploidetect segments (gain / loss / copy-neutral), and
   test systematically whether nanopore-only breakpoints coincide with
   copy-number step boundaries, as they do at ERBB2 in POG137.
2. **Methylation and phasing on the same genes.** Join promoter methylation
   (`AveragedMethylation_Promoters`), allelic DMRs and phase-block coverage to
   the genes hit by nanopore SVs, per sample.
3. **Background-corrected recurrence.** Replace the sample-count ranking with
   a length- and fragile-site-aware expectation so long genes stop leading.
4. **HCC1395 / SEQC2 validation.** Run the same reconciliation on the
   HCC1395 cell-line pair against the SEQC2 somatic truth set to measure
   caller precision and recall per SV type.
5. **Raw-read access via EGA** (study EGAS00001001159) for read-level review
   of nanopore-only calls, starting with the ERBB2 breakpoints and the
   CREBBP deletion in POG044.

## Layout and reproduction

- `src/pogcohort/` — `chrom.py` (chromosome-name normalization),
  `sv_reconciliation.py`, `genes.py` (GENCODE parser), `cosmic.py`,
  `sv_annotation.py`
- `src/run_sv_reconciliation.py`, `src/run_cancer_gene_annotation.py` — drivers
  (`--samples A B C` for a subset)
- `src/data_acquisition/` — site crawl and download
- `tests/` — pytest (33 tests)
- `data/` — inputs, not committed; `results/` — large outputs, not committed;
  `docs/` — small tables and figures, committed

```bash
python3.12 -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/python -m pytest
.venv/bin/python src/run_sv_reconciliation.py
.venv/bin/python src/run_cancer_gene_annotation.py
```
