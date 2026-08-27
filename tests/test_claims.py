"""Every numeric claim in the manuscript still matches its result artifact.

The claim table lives in scripts/verify_claims.py, which re-reads each artifact at
check time rather than comparing against a transcribed value. This wrapper puts it
under `make test` and CI so a number cannot silently drift from its artifact.

Kept as one aggregate test rather than one test per claim: the claims are rows in a
data table, not distinct behaviours, and a per-claim parametrisation would swamp the
suite's test count without telling you anything the aggregate failure message does not.

Unlike the pipeline tests this needs no interim data -- the artifacts it reads are
committed, so it runs on a bare clone.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERIFIER = ROOT / "scripts" / "verify_claims.py"


def _verifier():
    spec = importlib.util.spec_from_file_location("verify_claims", VERIFIER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_verifier_is_present_and_nonempty():
    assert VERIFIER.exists(), "scripts/verify_claims.py is the machine record of every claim"
    assert len(_verifier().CLAIMS) >= 50


def test_every_claim_matches_its_artifact():
    problems = []
    for claim in _verifier().CLAIMS:
        try:
            ok, actual = claim.check()
        except (FileNotFoundError, KeyError, IndexError, TypeError) as e:
            problems.append(f"{claim.text}: cannot resolve {claim.file}:{claim.path} ({e!r})")
            continue
        if not ok:
            problems.append(
                f"{claim.text}: claimed {claim.expect}, {claim.file}:{claim.path} = {actual}"
            )
    assert not problems, "\n".join(problems)
