"""Label integrity: decoy purity + complex-disjoint split (SPEC test plan)."""
from __future__ import annotations
from pathlib import Path

import pandas as pd
import pytest

from emperor import config as C

LAB = C.INTERIM / "labels.parquet"
pytestmark = pytest.mark.skipif(not LAB.exists(), reason="run `make labels` first")


def _labels():
    return pd.read_parquet(LAB)


def test_positive_not_in_negatives():
    lab = _labels()
    pos = {frozenset((r.uniprot_a, r.uniprot_b)) for r in lab[lab.label == 1].itertuples()}
    neg = {frozenset((r.uniprot_a, r.uniprot_b)) for r in lab[lab.label == 0].itertuples()}
    assert not (pos & neg)


def test_complex_disjoint_split():
    lab = _labels()
    pos = lab[lab.label == 1]
    span = pos.groupby("complex_id")["split"].nunique()
    assert (span <= 1).all(), "a CORUM complex spans cal and test"


def test_splits_present():
    lab = _labels()
    for split in ("cal", "test"):
        s = lab[lab.split == split]
        assert (s.label == 1).any() and (s.label == 0).any()
