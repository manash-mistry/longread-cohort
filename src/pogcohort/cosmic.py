"""COSMIC Cancer Gene Census (the copy shipped in data/pog/phasing/).

The Census is COSMIC's curated list of genes with causal evidence in cancer.
Columns used: "Gene Symbol", "Ensembl ID" (joined to GENCODE on the
version-less Ensembl gene ID), "Tier" (1 = strong evidence, 2 = emerging),
"Role in Cancer" (oncogene / TSG / fusion), "Tumour Types(Somatic)" — a
free-text, comma-separated list of tumour types in which somatic mutations
are documented. The POG table gives tumour type as a short cohort code
(BRCA, COLO, ...), so we match with a hand-written keyword list per code.
This is crude: a "no" means no keyword matched, not that the gene is
irrelevant to that tumour type.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

CENSUS = Path("data/pog/phasing/cosmic-gene-census.csv")

# POG `tumour_type_cohort` code -> substrings searched (case-insensitive) in
# the Census "Tumour Types(Somatic)" text.
COHORT_KEYWORDS: dict[str, list[str]] = {
    "BRCA": ["breast"],
    "SARC": ["sarcoma", "gist", "rhabdo", "leiomyo", "liposarcoma", "ewing"],
    "COLO": ["colorectal", "colon", "rectal", "colorectum"],
    "LUNG": ["lung", "nsclc", "sclc"],
    "PANC": ["pancrea"],
    "OV": ["ovar"],
    "CNS-PNS": ["glioma", "glioblastoma", "astrocytoma", "oligodendroglioma", "medulloblastoma", "cns", "brain"],
    "CHOL": ["cholangio", "biliary", "bile duct", "gallbladder"],
    "SKCM": ["melanoma"],
    "LYMP": ["lymphoma", "dlbcl", "b-cell", "t-cell", "hodgkin", "all", "cll"],
    "HNSC": ["head and neck", "hnscc", "oral", "laryn", "pharyn", "salivary"],
    "UVM": ["uveal", "melanoma"],
    "UCEC": ["endometri", "uterine"],
    "THYM": ["thymoma", "thymic"],
    "STAD": ["gastric", "stomach", "oesophag", "esophag", "gej"],
    "ESCA": ["oesophag", "esophag", "gastric", "gej"],
    "KDNY": ["renal", "kidney", "rcc", "wilms"],
    "CERV": ["cervi"],
    "ACC": ["adrenocortical", "adrenal"],
    "THCA": ["thyroid"],
    "PRAD": ["prostate"],
    "HCC": ["hepatocellular", "liver", "hcc"],
    "BLCA": ["bladder", "urothelial"],
    "BCC": ["basal cell", "skin"],
    "MISC": [],
    "SECR": [],
}


def load_census(path: Path = CENSUS) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    df = df.rename(columns={
        "Gene Symbol": "cosmic_symbol", "Ensembl ID": "gene_id", "Tier": "cosmic_tier",
        "Role in Cancer": "cosmic_role", "Tumour Types(Somatic)": "cosmic_tumour_types",
        "Hallmark": "cosmic_hallmark"})
    df["gene_id"] = df["gene_id"].str.split(".").str[0]
    return df[["gene_id", "cosmic_symbol", "cosmic_tier", "cosmic_role",
               "cosmic_hallmark", "cosmic_tumour_types"]]


def tumour_type_match(cosmic_tumour_types: str, cohort_code: str) -> str:
    """'yes' if any keyword for the cohort code occurs in the Census text,
    'no' if keywords exist but none match, 'unknown' if the code has none."""
    kws = COHORT_KEYWORDS.get(cohort_code)
    if not kws or not cosmic_tumour_types:
        return "unknown"
    text = cosmic_tumour_types.lower()
    return "yes" if any(k in text for k in kws) else "no"
