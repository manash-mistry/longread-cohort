"""Count MAVIS-merged SVs by platform support and type, per sample and cohort-wide.

Usage:
    python src/run_sv_reconciliation.py --samples POG044 POG049 POG068   # subset
    python src/run_sv_reconciliation.py                                  # all 43

Writes to results/sv_reconciliation/<tag>/ :
    counts_per_sample_raw.tsv, counts_per_sample_filtered.tsv,
    counts_cohort_raw.tsv, counts_cohort_filtered.tsv,
    examples_raw.tsv (five rows per category)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pogcohort.sv_reconciliation import (  # noqa: E402
    apply_paper_filters, count_by_sample, count_cohort, examples, load_combined)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="*", help="library IDs to keep (default: all)")
    ap.add_argument("--tag", default=None, help="output subfolder name")
    ap.add_argument("--unfiltered", action="store_true",
                    help="also write the no-filter view (default: paper filters + rescued)")
    args = ap.parse_args()

    tag = args.tag or ("subset_" + "_".join(args.samples) if args.samples else "all")
    out = Path("results/sv_reconciliation") / tag
    out.mkdir(parents=True, exist_ok=True)

    df = load_combined(samples=args.samples, filtered=False)
    print(f"loaded {len(df):,} distinct SV rows for {df['library'].nunique()} samples")

    flt = apply_paper_filters(df)
    print(f"after paper filters (+ nanopore_rescued): {len(flt):,} rows")
    if not args.unfiltered:
        df = flt  # "raw" outputs below then equal the filtered ones

    count_by_sample(df).to_csv(out / "counts_per_sample_raw.tsv", sep="\t")
    count_by_sample(flt).to_csv(out / "counts_per_sample_filtered.tsv", sep="\t")
    count_cohort(df).to_csv(out / "counts_cohort_raw.tsv", sep="\t")
    count_cohort(flt).to_csv(out / "counts_cohort_filtered.tsv", sep="\t")
    examples(df).to_csv(out / "examples_raw.tsv", sep="\t", index=False)

    print("\n== cohort, raw ==");       print(count_cohort(df).to_string())
    print("\n== cohort, filtered ==");  print(count_cohort(flt).to_string())
    print("\n== per sample, raw ==");   print(count_by_sample(df).to_string())
    print(f"\nwrote {out}/")


if __name__ == "__main__":
    main()
