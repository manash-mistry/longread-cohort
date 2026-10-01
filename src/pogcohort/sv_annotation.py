"""Annotate nanopore-supported SVs against genes and the COSMIC Census.

For each SV in a nanopore category (nanopore_only, nanopore_rescued) and each
gene it touches, one row with `effect`:

  whole_gene        the gene lies entirely inside the SV span (del/dup/inv)
  breakpoint_exon   a breakpoint falls inside the gene, within a merged exon
  breakpoint_intron a breakpoint falls inside the gene, between exons

Translocations and insertions have no span, so only their breakpoints are
tested. For del/dup/inv the span is [break1_position_start, break2_position_end].
Breakpoints are tested with their start..end uncertainty window.

Spans above LARGE_EVENT_BP are flagged `large_event`; a 70 Mb "duplication"
contains hundreds of genes and says nothing specific about any one of them.
"""
from __future__ import annotations

import pandas as pd

from pogcohort.cosmic import load_census, tumour_type_match
from pogcohort.genes import load_genes
from pogcohort.sv_reconciliation import NANOPORE_CATEGORIES, sample_type_lookup

LARGE_EVENT_BP = 5_000_000
SPAN_TYPES = {"deletion", "duplication", "inversion"}


def _overlaps(genes_c: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    return genes_c[(genes_c["start"] <= end) & (genes_c["end"] >= start)]


def _is_exonic(exons_g: pd.DataFrame, start: int, end: int) -> bool:
    return bool(((exons_g["start"] <= end) & (exons_g["end"] >= start)).any())


def annotate(svs: pd.DataFrame, genes: pd.DataFrame, exons: pd.DataFrame) -> pd.DataFrame:
    """One row per (SV, gene) hit for the nanopore categories in `svs`."""
    svs = svs[svs["category"].isin(NANOPORE_CATEGORIES)]
    genes_by_chrom = {c: g.sort_values("start") for c, g in genes.groupby("chrom")}
    exons_by_gene = {g: e for g, e in exons.groupby("gene_id")}
    rows = []
    for sv in svs.itertuples(index=False):
        b1 = (sv.break1_chromosome, int(sv.break1_position_start), int(sv.break1_position_end))
        b2 = (sv.break2_chromosome, int(sv.break2_position_start), int(sv.break2_position_end))
        hits: dict[str, str] = {}
        if sv.event_type in SPAN_TYPES and b1[0] == b2[0]:
            span_s, span_e = b1[1], b2[2]
            for g in _overlaps(genes_by_chrom.get(b1[0], genes.iloc[0:0]), span_s, span_e).itertuples():
                if g.start >= span_s and g.end <= span_e:
                    hits[g.gene_id] = "whole_gene"
        for chrom, s, e in (b1, b2):
            for g in _overlaps(genes_by_chrom.get(chrom, genes.iloc[0:0]), s, e).itertuples():
                if g.gene_id in hits:
                    continue
                ex = exons_by_gene.get(g.gene_id)
                hits[g.gene_id] = "breakpoint_exon" if ex is not None and _is_exonic(ex, s, e) else "breakpoint_intron"
        for gid, effect in hits.items():
            rows.append((sv.library, sv.tracking_id, sv.category, sv.event_type,
                         b1[0], b1[1], b2[0], b2[1], sv.event_size, gid, effect))
    out = pd.DataFrame(rows, columns=[
        "library", "tracking_id", "category", "event_type", "break1_chrom", "break1_pos",
        "break2_chrom", "break2_pos", "event_size", "gene_id", "effect"])
    out = out.merge(genes[["gene_id", "gene_name", "gene_type", "start", "end"]], on="gene_id", how="left")
    out["large_event"] = out["event_size"].fillna(0) > LARGE_EVENT_BP
    return out


def add_cosmic_and_tumour_type(hits: pd.DataFrame) -> pd.DataFrame:
    """Join Census columns and the sample's tumour type; keep Census genes only."""
    census = load_census()
    samples = sample_type_lookup()[["tumour_type_cohort", "oncotree_tumour_type", "tumour_content"]]
    out = hits.merge(census, on="gene_id", how="inner")
    out = out.merge(samples, left_on="library", right_index=True, how="left")
    out["driver_in_tumour_type"] = [
        tumour_type_match(t, c) for t, c in zip(out["cosmic_tumour_types"], out["tumour_type_cohort"])]
    return out


def per_sample_counts(cosmic_hits: pd.DataFrame, all_samples: list[str], low: list[str]) -> pd.DataFrame:
    """Per sample: distinct Census genes hit, split by category and by
    breakpoint-level vs whole-gene-in-large-event."""
    h = cosmic_hits.copy()
    h["kind"] = h["large_event"].map({True: "whole_gene_large", False: "focal"})
    tab = (h.groupby(["library", "category", "kind"], observed=True)["gene_id"].nunique()
             .unstack(["category", "kind"]).reindex(all_samples).fillna(0).astype(int))
    tab.columns = ["_".join(c) for c in tab.columns]
    tab["cosmic_genes_focal"] = h[~h["large_event"]].groupby("library")["gene_id"].nunique().reindex(all_samples).fillna(0).astype(int)
    tab["cosmic_genes_driver_match_focal"] = (
        h[(~h["large_event"]) & (h["driver_in_tumour_type"] == "yes")]
        .groupby("library")["gene_id"].nunique().reindex(all_samples).fillna(0).astype(int))
    tab["low_agreement"] = tab.index.isin(low)
    return tab


def recurrent_genes(cosmic_hits: pd.DataFrame) -> pd.DataFrame:
    """Genes ranked by number of distinct samples with a focal hit."""
    h = cosmic_hits[~cosmic_hits["large_event"]]
    g = (h.groupby(["gene_id", "gene_name", "cosmic_tier", "cosmic_role"], observed=True)
           .agg(n_samples=("library", "nunique"), n_svs=("tracking_id", "nunique"),
                samples=("library", lambda s: ",".join(sorted(set(s)))),
                effects=("effect", lambda s: ",".join(f"{k}:{v}" for k, v in s.value_counts().items())),
                sv_types=("event_type", lambda s: ",".join(f"{k}:{v}" for k, v in s.value_counts().items())),
                driver_match=("driver_in_tumour_type", lambda s: int((s == "yes").sum())))
           .reset_index().sort_values(["n_samples", "n_svs"], ascending=False))
    g["n_samples_large_event"] = g["gene_id"].map(
        cosmic_hits[cosmic_hits["large_event"]].groupby("gene_id")["library"].nunique()).fillna(0).astype(int)
    return g
