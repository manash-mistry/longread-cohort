# CLAUDE.md

## Project

A cohort-level interpretation and annotation engine for long-read somatic cancer sequencing data.

Open-source pipelines (nf-core/pacsomatic, LRSomatic) handle per-sample variant calling. We build the layer above:

- aggregating long-read structural variants, phasing, and methylation across a cohort
- reconciling against short-read calls on the same tumors
- annotating against cancer knowledge bases (CIViC, OncoKB, COSMIC)
- joining to clinical data

**Research use only.**

### Datasets

- **First dataset:** the public Long-Read POG cohort (189 tumors, Oxford Nanopore), processed data from bcgsc.ca/downloads/nanopore_pog/.
- **Later:** per-sample validation on HCC1395 against the SEQC2 truth set.

## Rules

- Inventory and describe data files before writing any analysis code.
- Build and test on a small subset (3 samples) before running on the full cohort.
- For any count or summary statistic, show five raw examples behind it and explain how each was classified.
- Never invent file paths, URLs, accessions, or gene names. If unsure, say so and ask.
- Explain every genomics file format and concept the first time it appears, assuming the reader has a biology degree but no bioinformatics background.
- Python 3.12, pandas, pytest. Keep code in `src/`, tests in `tests/`, data in `data/` (never committed to git), outputs in `results/`.
- Commit after each working step with a clear message.
