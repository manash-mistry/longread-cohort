"""Annotate nanopore-only / nanopore-rescued SVs against GENCODE genes and the
COSMIC Cancer Gene Census; summarize per sample and per gene; draw a bar chart.

Usage:
    python src/run_cancer_gene_annotation.py --samples POG044 POG049 POG068
    python src/run_cancer_gene_annotation.py            # all 43

Writes results/cancer_gene_sv/<tag>/ (full hit tables) and, for the all-sample
run, small summary tables + chart to docs/.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pogcohort.cosmic import CENSUS  # noqa: E402
from pogcohort.genes import GENCODE_VERSION, load_genes  # noqa: E402
from pogcohort.sv_annotation import (  # noqa: E402
    add_cosmic_and_tumour_type, annotate, per_sample_counts, recurrent_genes)
from pogcohort.sv_reconciliation import LOW_AGREEMENT_SAMPLES, load_combined  # noqa: E402

# Categorical slots 1 and 2 of the reference palette (docs: dataviz skill).
COLORS = {"nanopore_only": "#2a78d6", "nanopore_rescued": "#eb6834"}


def bar_chart(rec: pd.DataFrame, hits: pd.DataFrame, path: Path, top: int = 20) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    genes = rec.head(top)["gene_name"].tolist()
    focal = hits[~hits["large_event"] & hits["gene_name"].isin(genes)]
    counts = (focal.groupby(["gene_name", "category"], observed=True)["library"].nunique()
                   .unstack("category").reindex(genes).fillna(0))
    for c in COLORS:
        if c not in counts:
            counts[c] = 0

    fig, ax = plt.subplots(figsize=(7, 0.32 * len(genes) + 1.4), dpi=150)
    y = range(len(genes))[::-1]
    left = pd.Series(0.0, index=genes)
    for cat, color in COLORS.items():
        ax.barh(list(y), counts[cat].values, left=left.values, color=color, height=0.62,
                edgecolor="white", linewidth=1.5, label=cat.replace("_", " "))
        left = left + counts[cat]
    ax.set_yticks(list(y)); ax.set_yticklabels(genes, fontsize=8, style="italic")
    ax.set_xlabel("samples with a focal nanopore-supported SV in gene", fontsize=8)
    ax.tick_params(axis="x", labelsize=8)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="x", color="#e5e4e0", linewidth=0.6); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("COSMIC Census genes most often hit by nanopore-only / rescued SVs\n"
                 f"(43 tumours, focal events ≤5 Mb; {GENCODE_VERSION.split(' (')[0]})", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="*")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    tag = args.tag or ("subset_" + "_".join(args.samples) if args.samples else "all")
    out = Path("results/cancer_gene_sv") / tag
    out.mkdir(parents=True, exist_ok=True)

    svs = load_combined(samples=args.samples)          # paper filters + rescued
    genes, exons = load_genes()
    print(f"{GENCODE_VERSION}: {len(genes):,} genes, {len(exons):,} merged exons; census {CENSUS}")

    hits = annotate(svs, genes, exons)
    cosmic = add_cosmic_and_tumour_type(hits)
    all_samples = sorted(svs["library"].unique())
    per_sample = per_sample_counts(cosmic, all_samples, LOW_AGREEMENT_SAMPLES)
    rec = recurrent_genes(cosmic)

    hits.to_csv(out / "gene_hits_all.tsv", sep="\t", index=False)
    cosmic.to_csv(out / "cosmic_hits.tsv", sep="\t", index=False)
    per_sample.to_csv(out / "per_sample_counts.tsv", sep="\t")
    rec.to_csv(out / "recurrent_genes.tsv", sep="\t", index=False)

    n_sv = svs["category"].isin(["nanopore_only", "nanopore_rescued"]).sum()
    print(f"nanopore SVs: {n_sv:,}; gene hits: {len(hits):,}; census hits: {len(cosmic):,} "
          f"(focal {(~cosmic.large_event).sum():,}, in large events {cosmic.large_event.sum():,})")
    print("\n== per sample ==");  print(per_sample.to_string())
    print("\n== top recurrent census genes (focal) ==")
    print(rec.head(25)[["gene_name", "cosmic_tier", "cosmic_role", "n_samples", "n_svs",
                        "driver_match", "n_samples_large_event", "effects", "sv_types"]].to_string(index=False))

    if not args.samples:
        docs = Path("docs"); docs.mkdir(exist_ok=True)
        per_sample.to_csv(docs / "cancer_gene_sv_per_sample.tsv", sep="\t")
        rec.head(50).to_csv(docs / "cancer_gene_sv_recurrent_genes.tsv", sep="\t", index=False)
        bar_chart(rec, cosmic, docs / "cancer_gene_sv_top_genes.png")
        print("wrote docs/ tables and chart")


if __name__ == "__main__":
    main()
