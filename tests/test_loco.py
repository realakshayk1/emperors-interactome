"""test_loco.py — leave-one-complex-member-out recovery (G-VALID) sanity + purity."""
import os
import tempfile

from emperor import loco_validation as L
from emperor import config as C

_TMP = os.path.join(tempfile.gettempdir(), "loco_test.json")  # never clobber the shipped artifact


def test_loco_gate_met_and_fair_null():
    o = L.run(n_perm=2000, seed=C.SEED, out_path=_TMP)
    # members recovered above eligible impostors, clearing pre-registered thresholds
    assert o["odds_ratio_member_vs_impostor"] >= 2.0
    assert o["permutation_p"] < 0.05
    assert o["gate_met"] is True
    # the impostor null must be non-trivial (not the tautological "random protein has no edge")
    assert o["observed_impostor_recovery"] > 0.05, "impostor null collapsed -> test not discriminating"
    assert o["observed_member_recovery"] > o["observed_impostor_recovery"]


def test_loco_no_depmap_dependency():
    # LOCO must not READ any DepMap artifact; certification is ipTM-vs-decoy only.
    # Scan executable code, not the module docstring (which legitimately says "no DepMap").
    import ast, inspect
    tree = ast.parse(inspect.getsource(L))
    tree.body = [n for n in tree.body
                 if not (isinstance(n, ast.Expr) and isinstance(getattr(n, "value", None), ast.Constant)
                         and isinstance(n.value.value, str))]
    code = ast.unparse(tree).lower()
    for tok in ("depmap", "coess", "co_ess", "gls", "essential"):
        assert tok not in code, f"purity firewall: '{tok}' present in loco_validation code"


def test_loco_recovery_is_membership_blind():
    # Certification depends only on nonconformity (1-ipTM) vs calibration decoys, so a
    # trial's recovery can't depend on its own membership label. Recompute recovery with
    # the certified-edge set and confirm it's a pure graph lookup (deterministic).
    o1 = L.run(n_perm=500, seed=1, out_path=_TMP)
    o2 = L.run(n_perm=500, seed=2, out_path=_TMP)
    # observed rates are permutation-independent (only the null uses the RNG)
    assert o1["observed_member_recovery"] == o2["observed_member_recovery"]
    assert o1["observed_impostor_recovery"] == o2["observed_impostor_recovery"]
