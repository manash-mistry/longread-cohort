import gzip

import pandas as pd
import pytest

from pogcohort.sv_reconciliation import (
    USECOLS, apply_paper_filters, count_cohort, load_combined)

HEADER = USECOLS
ROWS = [
    # tracking_id, library, dgv, event_type, tools, b1chr, b1s, b1e, b2chr, b2s, b2e, gene1, gene2
    ["POG1_False_delly-DEL1", "POG1", "", "deletion", "Illumina", "1", 100, 101, "1", 500, 501, "None", "None"],
    ["POG1_True_manta-DEL2", "POG1", "", "deletion", "Illumina;Nanopore_0", "2", 100, 101, "2", 900, 901, "A", "A"],
    ["POG1_nanomonSV_d_1", "POG1", "", "deletion", "Nanopore_0", "3", 100, 101, "3", 2000, 2001, "None", "None"],
    ["POG1_ID_5_1", "POG1", "", "translocation", "Nanopore_savana", "4", 100, 101, "7", 100, 101, "None", "None"],
    ["POG1_True_delly-INV1", "POG1", "", "inverted translocation", "Illumina", "5", 100, 101, "6", 100, 101, "None", "None"],
    ["POG1_True_delly-DEL9", "POG1", "nsv123", "deletion", "Illumina", "8", 100, 101, "8", 10000, 10001, "None", "None"],
    ["POG1_True_delly-DEL10", "POG1", "", "deletion", "Illumina", "X", 100, 101, "X", 10000, 10001, "None", "None"],
    ["POG1_True_delly-DEL11", "POG1", "", "deletion", "Illumina", "9", 100, 101, "9", 120, 121, "None", "None"],
    ["POG1_True_delly-DEL11", "POG1", "", "deletion", "Illumina", "9", 100, 101, "9", 120, 121, "None", "None"],  # dup id
    ["POG2_nanomonSV_i_1", "POG2", "", "insertion", "Nanopore_0", "1", 100, 100, "1", 100, 100, "None", "None"],
]


@pytest.fixture
def table(tmp_path):
    p = tmp_path / "combined.tsv.gz"
    with gzip.open(p, "wt") as f:
        f.write("\t".join(HEADER) + "\n")
        for r in ROWS:
            f.write("\t".join(map(str, r)) + "\n")
    return p


def test_load_categories_and_normalization(table):
    df = load_combined(table).set_index("tracking_id")
    assert len(df) == 9  # duplicate tracking_id dropped
    assert df.loc["POG1_False_delly-DEL1", "category"] == "illumina_only"
    assert df.loc["POG1_True_manta-DEL2", "category"] == "both"
    assert df.loc["POG1_nanomonSV_d_1", "category"] == "nanopore_only"
    assert df.loc["POG1_ID_5_1", "category"] == "nanopore_only"
    assert df.loc["POG1_ID_5_1", "savana"] and not df.loc["POG1_ID_5_1", "nanomonsv"]
    assert df.loc["POG1_True_delly-INV1", "event_type"] == "translocation"
    assert df.loc["POG1_False_delly-DEL1", "illumina_hq"] is False
    assert df.loc["POG1_nanomonSV_d_1", "illumina_hq"] is None
    assert df.loc["POG1_True_delly-DEL10", "break1_chromosome"] == "chrX"
    assert df.loc["POG1_False_delly-DEL1", "event_size"] == 401
    assert pd.isna(df.loc["POG1_ID_5_1", "event_size"])  # inter-chromosomal


def test_subset_by_sample(table):
    assert load_combined(table, samples=["POG2"])["library"].unique().tolist() == ["POG2"]


def test_paper_filters(table):
    df = load_combined(table)
    kept = set(apply_paper_filters(df)["tracking_id"])
    assert "POG1_False_delly-DEL1" not in kept   # Illumina-only, low quality
    assert "POG1_True_delly-DEL9" not in kept    # in DGV
    assert "POG1_True_delly-DEL10" not in kept   # chrX
    assert "POG1_True_delly-DEL11" not in kept   # 21 bp
    assert {"POG1_True_manta-DEL2", "POG1_nanomonSV_d_1", "POG1_ID_5_1",
            "POG1_True_delly-INV1"} <= kept
    assert "POG2_nanomonSV_i_1" in kept  # insertions exempt from the size filter


def test_merged_tracking_id_quality_flag(table):
    df = load_combined(table)
    # a merged row carrying both False and True Illumina IDs counts as high quality
    merged = pd.DataFrame([{
        "tracking_id": "POG1_False_delly-INV1;POG1_True_delly-INV2;POG1_nanomonSV_r_1",
        "library": "POG1", "dgv": "", "event_type": "inversion",
        "tools": "Illumina;Illumina;Nanopore_0", "break1_chromosome": "1",
        "break1_position_start": 100, "break1_position_end": 101,
        "break2_chromosome": "1", "break2_position_start": 5000,
        "break2_position_end": 5001, "gene1": "None", "gene2": "None"}])
    p = table.parent / "merged.tsv.gz"
    pd.concat([pd.read_csv(table, sep="\t", dtype=str), merged]).to_csv(p, sep="\t", index=False)
    row = load_combined(p).set_index("tracking_id").iloc[-1]
    assert row["illumina_hq"] is True and row["category"] == "both"
    assert row["nanomonsv"] and not row["savana"]
    # ...and a both-row with only a False flag is dropped by the paper filters
    both_false = df[df["tracking_id"] == "POG1_True_manta-DEL2"].copy()
    both_false["tracking_id"] = "POG1_False_manta-DEL3"
    both_false["illumina_hq"] = False
    assert apply_paper_filters(both_false).empty


def test_count_cohort(table):
    tab = count_cohort(load_combined(table))
    assert tab.loc["deletion", "illumina_only"] == 4
    assert tab.loc["deletion", "both"] == 1
    assert tab.loc["deletion", "nanopore_only"] == 1
    assert tab.loc["translocation", "total"] == 2
    assert tab.loc["insertion", "nanopore_only"] == 1
