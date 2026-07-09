"""Synthetic-null validity of conformal p-values + BH FDR control (SPEC test plan).

We generate calibration negatives and a test pool that mixes true positives
(lower nonconformity) with nulls drawn from the SAME distribution as the
negatives. Under exchangeability the conformal p-values of nulls are ~Uniform,
so BH must control FDR <= q (up to finite-sample slack) across many seeds.
"""
from __future__ import annotations
import numpy as np
import pytest

from emperor.conformal import benjamini_hochberg, conformal_pvalues


def test_null_pvalues_uniform():
    """With test == calibration distribution, conformal p-values are ~Uniform."""
    rng = np.random.default_rng(0)
    n_cal, n_test = 2000, 5000
    cal_neg = rng.normal(size=n_cal)
    test_null = rng.normal(size=n_test)          # same law as negatives
    p = conformal_pvalues(test_null, cal_neg)
    # mean ~ 0.5, and P(p <= t) ~ t (super-uniform / valid)
    assert abs(p.mean() - 0.5) < 0.02
    for t in (0.05, 0.1, 0.2):
        assert (p <= t).mean() <= t + 0.03


@pytest.mark.parametrize("q", [0.05, 0.10, 0.20])
def test_bh_fdr_control(q):
    """Empirical FDR <= q (within slack) averaged over 30 seeds."""
    fdrs = []
    for seed in range(30):
        rng = np.random.default_rng(seed)
        n_cal = 1500
        cal_neg = rng.normal(0, 1, n_cal)
        # test pool: 200 true positives (shifted lower = more interaction-like)
        # + 800 nulls from the negative law
        pos = rng.normal(-2.5, 1, 200)
        null = rng.normal(0, 1, 800)
        s_test = np.concatenate([pos, null])
        is_null = np.concatenate([np.zeros(200, bool), np.ones(800, bool)])
        p = conformal_pvalues(s_test, cal_neg)
        rej = benjamini_hochberg(p, q)
        if rej.sum():
            fdrs.append(is_null[rej].mean())
        else:
            fdrs.append(0.0)
    mean_fdr = float(np.mean(fdrs))
    # BH controls FDR in expectation; allow small Monte-Carlo slack
    assert mean_fdr <= q + 0.03, f"mean empirical FDR {mean_fdr:.3f} > q={q}"


def test_pvalue_monotone_in_confidence():
    """Lower nonconformity (more confident) -> smaller p-value."""
    rng = np.random.default_rng(1)
    cal_neg = rng.normal(size=1000)
    s = np.array([-3.0, 0.0, 3.0])
    p = conformal_pvalues(s, cal_neg)
    assert p[0] < p[1] < p[2]


def test_pvalue_range():
    cal_neg = np.zeros(99)
    p = conformal_pvalues(np.array([-1.0, 1.0]), cal_neg)
    assert (p >= 1 / 100).all() and (p <= 1.0).all()
