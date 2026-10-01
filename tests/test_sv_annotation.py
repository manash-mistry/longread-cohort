import pandas as pd

from pogcohort.cosmic import tumour_type_match
from pogcohort.genes import _merge_intervals
from pogcohort.sv_annotation import annotate

GENES = pd.DataFrame([
    # gene A chr1 1000-2000 (exons 1000-1100, 1500-1600); gene B chr1 5000-6000; gene C chr2 100-900
    ("chr1", 1000, 2000, "+", "GA", "A", "protein_coding"),
    ("chr1", 5000, 6000, "+", "GB", "B", "protein_coding"),
    ("chr2", 100, 900, "-", "GC", "C", "lncRNA"),
], columns=["chrom", "start", "end", "strand", "gene_id", "gene_name", "gene_type"])
EXONS = pd.DataFrame([
    ("GA", "chr1", 1000, 1100), ("GA", "chr1", 1500, 1600),
    ("GB", "chr1", 5000, 5100), ("GC", "chr2", 100, 900),
], columns=["gene_id", "chrom", "start", "end"])


def sv(tid, etype, c1, s1, c2, s2, cat="nanopore_only", size=None):
    return dict(library="POG1", tracking_id=tid, category=cat, event_type=etype,
                break1_chromosome=c1, break1_position_start=s1, break1_position_end=s1,
                break2_chromosome=c2, break2_position_start=s2, break2_position_end=s2,
                event_size=(s2 - s1 if c1 == c2 else None) if size is None else size)


def test_effects():
    svs = pd.DataFrame([
        sv("del_whole", "deletion", "chr1", 900, "chr1", 2100),      # A whole, breakpoints outside
        sv("del_exon", "deletion", "chr1", 1050, "chr1", 1200),      # break1 in exon of A
        sv("del_intron", "deletion", "chr1", 1300, "chr1", 3000),    # break1 in intron of A
        sv("tra", "translocation", "chr1", 5050, "chr2", 500),       # B exon, C exon
        sv("ins", "insertion", "chr1", 1550, "chr1", 1550, size=0),  # A exon
        sv("ill", "deletion", "chr1", 900, "chr1", 2100, cat="illumina_only"),  # ignored
    ])
    out = annotate(svs, GENES, EXONS).set_index(["tracking_id", "gene_name"])["effect"]
    assert out[("del_whole", "A")] == "whole_gene"
    assert out[("del_exon", "A")] == "breakpoint_exon"
    assert out[("del_intron", "A")] == "breakpoint_intron"
    assert out[("tra", "B")] == "breakpoint_exon" and out[("tra", "C")] == "breakpoint_exon"
    assert out[("ins", "A")] == "breakpoint_exon"
    assert "ill" not in out.index.get_level_values(0)


def test_large_event_flag():
    svs = pd.DataFrame([sv("big", "duplication", "chr1", 1, "chr1", 10_000_000)])
    out = annotate(svs, GENES, EXONS)
    assert set(out["gene_name"]) == {"A", "B"} and out["large_event"].all()


def test_merge_intervals():
    e = pd.DataFrame([("G", "chr1", 10, 20), ("G", "chr1", 15, 30), ("G", "chr1", 40, 50)],
                     columns=["gene_id", "chrom", "start", "end"])
    m = _merge_intervals(e)
    assert m[["start", "end"]].values.tolist() == [[10, 30], [40, 50]]


def test_tumour_type_match():
    assert tumour_type_match("breast, ovarian", "BRCA") == "yes"
    assert tumour_type_match("AML", "BRCA") == "no"
    assert tumour_type_match("AML", "MISC") == "unknown"
    assert tumour_type_match("", "BRCA") == "unknown"
