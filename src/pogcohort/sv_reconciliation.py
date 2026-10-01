"""Long-read vs short-read SV reconciliation from the MAVIS combined table.

Input: SVs/SV_nanopore_illumina_combined.tsv.gz — one row per MAVIS-merged SV
for the 43 tumours with a nanopore-sequenced normal.

How a row's platform support is read (mirrors the authors' Figures/SV/SV_longPOG.Rmd
in github.com/bcgsc/long_read_pog):

* `tools` lists the inputs that supported the merged call, ";"-separated:
  "Illumina" / "Illumina_<N>" = short-read callers (delly, manta) — N is an
  insertion-size estimate; "Nanopore_0" = nanomonsv; "Nanopore_savana" = SAVANA.
* `tracking_id` is "<sample>_<high_quality>_<caller-id>" for Illumina-origin
  rows, where <high_quality> is True/False copied from the per-sample MAVIS
  `high_quality` column; nanopore-origin rows are "<sample>_nanomonSV_<id>" or
  "<sample>_ID_<n>_<m>" (SAVANA).

Platform category:
    both              "Illumina" and "Nanopore" both appear in `tools`
    illumina_only     only "Illumina"
    nanopore_only     only "Nanopore"
    nanopore_rescued  (filtered mode only) nanopore-supported call whose
                      Illumina partner carries the low-quality flag; the
                      authors drop these, we keep them as their own category

Low-agreement samples (<12% "both" after filters; see
results/sv_reconciliation_summary.md §4) are flagged in `low_agreement` so
every downstream table can carry the flag until the cause is explained.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from pogcohort.chrom import is_autosome, normalize_series

COMBINED_TSV = Path("data/pog/SVs/SV_nanopore_illumina_combined.tsv.gz")

# The authors fold "inverted translocation" into "translocation".
EVENT_TYPE_MAP = {"inverted translocation": "translocation"}
EVENT_ORDER = ["deletion", "duplication", "insertion", "inversion", "translocation"]
CATEGORY_ORDER = ["illumina_only", "nanopore_only", "both", "nanopore_rescued"]
NANOPORE_CATEGORIES = ["nanopore_only", "nanopore_rescued"]

# Samples with <12% platform agreement after the paper filters. Two are
# hypermutated/MSI (POG884, POG986), three are SV-quiet genomes; no shared
# sequencing-QC cause found. Flagged, not excluded.
LOW_AGREEMENT_SAMPLES = ["POG884", "POG986", "POG530", "POG044", "POG1000"]

USECOLS = [
    "tracking_id", "library", "dgv", "event_type", "tools",
    "break1_chromosome", "break1_position_start", "break1_position_end",
    "break2_chromosome", "break2_position_start", "break2_position_end",
    "gene1", "gene2",
]


def load_combined(path: Path = COMBINED_TSV, samples: list[str] | None = None,
                  filtered: bool = True) -> pd.DataFrame:
    """Load the MAVIS combined table and add platform/quality columns.

    Chromosomes are normalized to UCSC style. Rows are de-duplicated on
    `tracking_id` (the authors do the same: `distinct(x, tracking_id)`).
    With `filtered=True` (default) the paper filters are applied and
    low-quality-Illumina rows with nanopore support become `nanopore_rescued`.
    """
    df = pd.read_csv(path, sep="\t", usecols=USECOLS, dtype=str, keep_default_na=False)
    if samples is not None:
        df = df[df["library"].isin(samples)]
    df = df.drop_duplicates("tracking_id").copy()

    for c in ("break1_chromosome", "break2_chromosome"):
        df[c] = normalize_series(df[c])
    for c in ("break1_position_start", "break1_position_end",
              "break2_position_start", "break2_position_end"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["event_type"] = df["event_type"].replace(EVENT_TYPE_MAP)
    df["illumina"] = df["tools"].str.contains("Illumina")
    df["nanopore"] = df["tools"].str.contains("Nanopore")
    df["nanomonsv"] = df["tools"].str.contains("Nanopore_0")
    df["savana"] = df["tools"].str.contains("Nanopore_savana")
    df["category"] = pd.Categorical(
        df.apply(_category, axis=1), categories=CATEGORY_ORDER)
    df["illumina_hq"] = df["tracking_id"].map(_illumina_hq)
    # Size as the authors compute it (break2 end − break1 start); NaN for
    # inter-chromosomal events.
    same_chrom = df["break1_chromosome"] == df["break2_chromosome"]
    df["event_size"] = (df["break2_position_end"] - df["break1_position_start"]).where(same_chrom)
    df["low_agreement"] = df["library"].isin(LOW_AGREEMENT_SAMPLES)
    if filtered:
        df = apply_paper_filters(df)
    return df


def _category(row) -> str:
    if row["illumina"] and row["nanopore"]:
        return "both"
    if row["nanopore"]:
        return "nanopore_only"
    return "illumina_only"


def _illumina_hq(tracking_id: str):
    """Illumina high_quality flag read the way SV_longPOG.Rmd reads it.

    A merged row's tracking_id can concatenate several caller IDs with ";", so
    both "True" and "False" may appear; the authors treat that as True.
    None when no Illumina ID is present (nanopore-only rows).
    """
    has_true = "_True_" in tracking_id or tracking_id.startswith("True_")
    has_false = "_False_" in tracking_id
    if has_true:
        return True
    if has_false:
        return False
    return None


def apply_paper_filters(df: pd.DataFrame) -> pd.DataFrame:
    """The filters SV_longPOG.Rmd applies before its concordance figure (Fig 2B).

    * `dgv` blank — drop events overlapping the Database of Genomic Variants
      (catalogue of common germline SVs; overlap suggests not somatic).
    * event_size > 50 bp. Insertions are exempt here: their break1/break2 are
      the same base so event_size is 0, and the authors size them from the
      separate `_ins` table instead (not reproduced in this pass).
    * autosomes only (no X/Y).
    * Illumina high_quality flag False: the Rmd drops every such row
      (`filter(Illumina_High_Quality == 'True' | ... == 'Nanopore_only')`).
      We drop them only when Illumina is the sole support; when nanopore
      also supports the call it is kept and re-labelled `nanopore_rescued`.
    """
    base = (
        (df["dgv"] == "")
        & ((df["event_type"] == "insertion") | df["event_size"].isna() | (df["event_size"] > 50))
        & df["break1_chromosome"].map(is_autosome)
        & df["break2_chromosome"].map(is_autosome)
    )
    low_q = df["illumina_hq"] == False  # noqa: E712
    out = df[base & ~(low_q & ~df["nanopore"])].copy()
    rescued = (out["illumina_hq"] == False) & out["nanopore"]  # noqa: E712
    out["category"] = out["category"].astype(str).where(~rescued, "nanopore_rescued")
    out["category"] = pd.Categorical(out["category"], categories=CATEGORY_ORDER)
    return out


def count_by_sample(df: pd.DataFrame) -> pd.DataFrame:
    """Rows = (library, event_type); columns = category counts + total."""
    tab = (df.groupby(["library", "event_type", "category"], observed=False)
             .size().unstack("category").fillna(0).astype(int))
    tab = tab.reindex(columns=CATEGORY_ORDER)
    tab["total"] = tab.sum(axis=1)
    tab["low_agreement"] = tab.index.get_level_values("library").isin(LOW_AGREEMENT_SAMPLES)
    return tab


def count_cohort(df: pd.DataFrame) -> pd.DataFrame:
    """Rows = event_type; columns = category counts, total, percent in both."""
    tab = (df.groupby(["event_type", "category"], observed=False)
             .size().unstack("category").fillna(0).astype(int))
    tab = tab.reindex(index=[e for e in EVENT_ORDER if e in tab.index], columns=CATEGORY_ORDER)
    tab["total"] = tab.sum(axis=1)
    tab["pct_both"] = (100 * tab["both"] / tab["total"]).round(1)
    return tab


def sample_type_lookup() -> pd.DataFrame:
    """Map MAVIS/VCF sample names to clinical rows (tumour type, content...).

    The MAVIS name is POG_ID, or POG_ID-<biopsy> for the two patients with
    two nanopore-sequenced biopsies (POG117, POG415) — verified via SAVANA
    library IDs, see results/data_inventory.md §4.1.
    """
    s = pd.read_csv("data/pog/supplementary_tables/Supplementary_Table_1_samples.tsv",
                    sep="\t", dtype=str, keep_default_na=False)
    s = s[s["has_normal"] == "Y"].copy()
    multi = s["POG_ID"].duplicated(keep=False)
    s["library"] = s["POG_ID"].where(~multi, s["POG_ID"] + "-" + s["biopsy"])
    return s.set_index("library")


def examples(df: pd.DataFrame, n: int = 5) -> pd.DataFrame:
    """First `n` rows of each category, with the columns needed to read them."""
    cols = ["category", "library", "tracking_id", "event_type", "tools", "illumina_hq",
            "break1_chromosome", "break1_position_start", "break2_chromosome",
            "break2_position_start", "event_size", "dgv", "gene1", "gene2"]
    return (df.sort_values(["library", "tracking_id"])
              .groupby("category", observed=False).head(n)[cols]
              .sort_values(["category", "library"]))
