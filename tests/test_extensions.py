"""Tests for the Phase 2/4 extension modules: pDockQ formula + FDR-controlled nomination sets."""
import math
import json
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_pdockq_formula_monotone():
    """pDockQ must be monotone increasing in x = mean_if_plddt * ln(n_contacts), bounded (0.018, 0.742)."""
    def pd(x):
        return 0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018
    lo, mid, hi = pd(50), pd(152), pd(400)
    assert lo < mid < hi
    assert 0.018 <= lo and hi <= 0.724 + 0.018 + 1e-9
    # midpoint of the sigmoid is at x=152.611 -> ~0.5*0.724 + 0.018
    assert abs(pd(152.611) - (0.362 + 0.018)) < 1e-3


@pytest.mark.skipif(not (ROOT / "data/processed/nomination_sets.json").exists(),
                    reason="nomination_sets.json not built")
def test_nomination_sets_selection_is_conformal_only():
    """Nomination sets must be selected purely conformally; DepMap only ranks (purity firewall)."""
    r = json.loads((ROOT / "data/processed/nomination_sets.json").read_text())
    assert r["q"] == 0.1
    assert "certified" in r["selection"].lower()
    # every nomination must carry a conf_pvalue (the selection basis); coess may be None
    for s in r["complexes"]:
        for n in s["nominations"]:
            assert n["conf_pvalue"] is not None
    # total nominations equals sum over complexes (dedup consistency)
    assert r["n_total_nominations"] == sum(s["n_nominations"] for s in r["complexes"])


@pytest.mark.skipif(not (ROOT / "data/mimicry/auditability_boundary.json").exists(),
                    reason="auditability_boundary.json not built")
def test_second_map_not_auditable():
    """The host-pathogen deposit must be recorded as positives-only / not auditable."""
    d = json.loads((ROOT / "data/mimicry/auditability_boundary.json").read_text())
    assert d["map2_hostpathogen"]["ships_native_null"] is False
    assert d["map2_hostpathogen"]["auditable"] is False
    assert d["map1_cm4ai"]["auditable"] is True
