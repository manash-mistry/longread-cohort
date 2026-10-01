# longread-cohort

A cohort-level interpretation and annotation engine for long-read somatic cancer sequencing data. Open-source pipelines (nf-core/pacsomatic, LRSomatic) handle per-sample variant calling; this project builds the layer above: aggregating long-read structural variants, phasing, and methylation across a cohort, reconciling against short-read calls on the same tumors, annotating against cancer knowledge bases (CIViC, OncoKB, COSMIC), and joining to clinical data. The first dataset is the public Long-Read POG cohort (189 tumors, Oxford Nanopore), with per-sample validation later on HCC1395 against the SEQC2 truth set. Research use only.

## Layout

- `src/` code (`pogcohort/` package, `run_*.py` drivers, `data_acquisition/`)
- `tests/` pytest suite
- `data/` raw inputs, not committed (`data/pog/` = Long-Read POG download, `data/reference/` = GENCODE)
- `results/` large outputs, not committed
- `docs/` small summary tables and figures (<1 MB each), committed

## Summaries

- [docs/cancer_gene_sv_summary.md](docs/cancer_gene_sv_summary.md) — nanopore-only / rescued SVs in COSMIC Census genes
