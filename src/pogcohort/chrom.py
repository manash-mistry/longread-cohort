"""Chromosome-name normalization shared by every loader.

The POG files use two conventions for the same GRCh38 chromosomes:

UCSC style ("chr1", "chrX", "chrM"):
    SVs/somatic_vcfs/NanomonSV/*.vcf, SVs/somatic_vcfs/SAVANA/*.vcf,
    short_read/small_mutations/vcfs/*/*.vcf, Methylation/* (all),
    phasing/ensembl100-genes+promoters-blocks-summary-tumours.csv,
    ecDNA/POG816_NRG1_haplotype_methylation_source_table.txt

Ensembl style ("1", "X", "MT"):
    SVs/SV_nanopore_illumina_combined*.tsv.gz (break1/break2_chromosome),
    short_read/mavis/**/mavis_summary_somatic_gd-*.tab,
    ecDNA/POG816_NRG1_CN_ploidetect_source_table.txt,
    phasing/cosmic-gene-census.csv (Genome Location column)

Every loader should call `to_ucsc` (the project default) so joins never
fail on a "chr" prefix.
"""
from __future__ import annotations

import pandas as pd

# Mitochondrial naming differs beyond the prefix: Ensembl "MT", UCSC "chrM".
_MITO_ALIASES = {"M", "MT", "CHRM", "CHRMT"}


def to_ucsc(name) -> str | None:
    """Return the UCSC form of a chromosome name ("chr1", "chrX", "chrM").

    Names that are already UCSC-style are returned unchanged; unplaced and
    random contigs (e.g. "chrUn_KI270302v1", "GL000009.2") are passed through
    with only the prefix fixed. None/NaN/empty input returns None.
    """
    if name is None or (isinstance(name, float) and pd.isna(name)):
        return None
    s = str(name).strip()
    if not s or s.lower() in {"none", "nan"}:
        return None
    if s.upper() in _MITO_ALIASES:
        return "chrM"
    if s.lower().startswith("chr"):
        return "chr" + s[3:]
    return "chr" + s


def to_ensembl(name) -> str | None:
    """Return the Ensembl form ("1", "X", "MT"). Inverse of `to_ucsc`."""
    u = to_ucsc(name)
    if u is None:
        return None
    if u == "chrM":
        return "MT"
    return u[3:]


def normalize_series(s: pd.Series, style: str = "ucsc") -> pd.Series:
    """Normalize a pandas Series of chromosome names to one style."""
    fn = {"ucsc": to_ucsc, "ensembl": to_ensembl}[style]
    return s.map(fn)


def is_autosome(name) -> bool:
    """True for chr1..chr22 in either naming style."""
    e = to_ensembl(name)
    return e is not None and e.isdigit() and 1 <= int(e) <= 22
