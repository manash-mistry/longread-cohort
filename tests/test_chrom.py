import math

import pandas as pd
import pytest

from pogcohort.chrom import is_autosome, normalize_series, to_ensembl, to_ucsc


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("1", "chr1"),
        ("chr1", "chr1"),
        ("X", "chrX"),
        ("chrY", "chrY"),
        ("MT", "chrM"),
        ("M", "chrM"),
        ("chrM", "chrM"),
        ("chrUn_KI270302v1", "chrUn_KI270302v1"),
        (" 7 ", "chr7"),
        (12, "chr12"),  # pandas may read Ensembl-style chromosomes as ints
    ],
)
def test_to_ucsc(raw, expected):
    assert to_ucsc(raw) == expected


@pytest.mark.parametrize("raw", [None, math.nan, "", "None", "nan"])
def test_to_ucsc_missing(raw):
    assert to_ucsc(raw) is None


@pytest.mark.parametrize(
    "raw, expected",
    [("chr1", "1"), ("1", "1"), ("chrM", "MT"), ("MT", "MT"), ("chrX", "X")],
)
def test_to_ensembl(raw, expected):
    assert to_ensembl(raw) == expected


def test_round_trip():
    for n in ["1", "22", "X", "Y", "MT"]:
        assert to_ensembl(to_ucsc(n)) == n


def test_normalize_series_mixed_input():
    s = pd.Series(["1", "chr2", "X", "MT", None])
    out = normalize_series(s)
    assert out[:4].tolist() == ["chr1", "chr2", "chrX", "chrM"] and pd.isna(out[4])
    out = normalize_series(s, "ensembl")
    assert out[:4].tolist() == ["1", "2", "X", "MT"] and pd.isna(out[4])


def test_is_autosome():
    assert is_autosome("chr1") and is_autosome("22")
    assert not is_autosome("X") and not is_autosome("chrM") and not is_autosome(None)
