"""Gene and exon coordinates from a GENCODE GTF.

GTF (Gene Transfer Format): tab-separated, one feature per line — chromosome,
source, feature type (gene/transcript/exon/...), 1-based start, end, score,
strand, frame, and a free-text attribute column (`key "value"; ...`).
GENCODE GRCh38 GTFs use UCSC chromosome names ("chr1").

We keep two tables:
  genes — one row per gene: chrom, start, end, strand, gene_id (version
          stripped), gene_name, gene_type
  exons — one row per merged exonic interval per gene (union of all
          transcripts' exons), so "is this base exonic?" has one answer.
Both are cached as TSV next to the GTF because parsing takes ~1 min.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

from pogcohort.chrom import normalize_series

GTF = Path("data/reference/gencode.v50.annotation.gtf.gz")
GENCODE_VERSION = "GENCODE v50 (GRCh38, Ensembl 116; file dated 2026-04-08)"

_ATTR = re.compile(r'(\w+) "([^"]*)"')


def _parse(gtf: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes, exons = [], []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feat, start, end, _, strand, _, attr = line.rstrip("\n").split("\t")
            if feat == "gene":
                a = dict(_ATTR.findall(attr))
                genes.append((chrom, int(start), int(end), strand,
                              a["gene_id"].split(".")[0], a.get("gene_name", ""), a.get("gene_type", "")))
            elif feat == "exon":
                gid = attr.split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]
                exons.append((gid, chrom, int(start), int(end)))
    g = pd.DataFrame(genes, columns=["chrom", "start", "end", "strand", "gene_id", "gene_name", "gene_type"])
    e = pd.DataFrame(exons, columns=["gene_id", "chrom", "start", "end"])
    return g, _merge_intervals(e)


def _merge_intervals(e: pd.DataFrame) -> pd.DataFrame:
    """Union overlapping exon intervals within each gene."""
    e = e.sort_values(["gene_id", "start", "end"])
    out = []
    for gid, grp in e.groupby("gene_id", sort=False):
        chrom = grp["chrom"].iloc[0]
        cur_s = cur_e = None
        for s, t in zip(grp["start"].values, grp["end"].values):
            if cur_e is None or s > cur_e + 1:
                if cur_e is not None:
                    out.append((gid, chrom, cur_s, cur_e))
                cur_s, cur_e = s, t
            else:
                cur_e = max(cur_e, t)
        out.append((gid, chrom, cur_s, cur_e))
    return pd.DataFrame(out, columns=["gene_id", "chrom", "start", "end"])


def load_genes(gtf: Path = GTF) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (genes, exons) with UCSC chromosome names, using a TSV cache."""
    gc, ec = gtf.with_suffix(".genes.tsv"), gtf.with_suffix(".exons.tsv")
    if gc.exists() and ec.exists():
        g, e = pd.read_csv(gc, sep="\t", dtype={"chrom": str}), pd.read_csv(ec, sep="\t", dtype={"chrom": str})
    else:
        g, e = _parse(gtf)
        g.to_csv(gc, sep="\t", index=False)
        e.to_csv(ec, sep="\t", index=False)
    g["chrom"] = normalize_series(g["chrom"])
    e["chrom"] = normalize_series(e["chrom"])
    return g, e
