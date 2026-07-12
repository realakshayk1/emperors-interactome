"""Adversarial + validity tests for dependence-robust FDR (PLAN_V3 T0.4)."""
import numpy as np
import pytest

from emperor.conformal import benjamini_hochberg, conformal_pvalues
from emperor.dependence import benjamini_yekutieli, conformal_evalues, ebh, _harmonic


def test_by_subset_of_bh():
    """BY (arbitrary dependence) must reject a subset of BH on any p-vector/q."""
    rng = np.random.default_rng(0)
    for _ in range(50):
        p = rng.random(200)
        for q in (0.05, 0.1, 0.2):
            r_bh = benjamini_hochberg(p, q)
            r_by = benjamini_yekutieli(p, q)
            assert set(np.nonzero(r_by)[0]) <= set(np.nonzero(r_bh)[0])


def test_evalue_null_mean_one():
    """Harmonic-calibrated conformal e-values have null mean ~1 under exchangeability."""
    rng = np.random.default_rng(1)
    n = 2000
    # test points and calibration negatives drawn from the SAME law -> exchangeable null
    means = []
    for _ in range(200):
        cal = rng.normal(size=n)
        test = rng.normal(size=500)          # all true nulls
        means.append(conformal_evalues(test, cal).mean())
    assert abs(np.mean(means) - 1.0) < 0.05


def test_ebh_fdr_control_arbitrary():
    """e-BH controls FDR at q on a synthetic mix with a known null (exact-control check)."""
    rng = np.random.default_rng(2)
    n = 3000
    fdps = []
    q = 0.2
    for _ in range(300):
        cal = rng.normal(0, 1, n)                       # calibration negatives
        s_null = rng.normal(0, 1, 400)                  # true nulls (lower=more signal-like)
        s_sig = rng.normal(-3, 1, 100)                  # true signals
        s = np.concatenate([s_null, s_sig])
        y = np.array([0] * 400 + [1] * 100)
        e = conformal_evalues(s, cal)
        r = ebh(e, q)
        if r.sum():
            fdps.append((1 - y[r]).mean())
    assert np.mean(fdps) <= q + 0.02


def test_harmonic():
    assert _harmonic(1) == 1.0
    assert abs(_harmonic(4) - (1 + 0.5 + 1/3 + 0.25)) < 1e-12
