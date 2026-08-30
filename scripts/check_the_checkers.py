"""Negative control for the verification harness: perturb the inputs and confirm it complains.

Six review rounds on this repo found ten assertions that could not fail -- several of them
inside scripts written specifically to catch that. The pattern is always the same: the check
compares a quantity against itself, or against something the constructor forces, and passes
forever. A green harness is evidence of nothing unless a broken input turns it red.

So this runs the claim checker against deliberately corrupted copies of the artifacts and
requires it to fail. It restores every file afterwards, and it refuses to run against a dirty
tree so an interrupted run cannot leave a mutated artifact behind.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"

# (artifact, dotted path into the JSON, how to break it)
MUTATIONS = [
    ("audit_summary.json", ["n_candidates"], lambda v: v + 1),
    ("dependence_robustness.json", ["n_cal_neg"], lambda v: v - 1),
    ("calibrator_comparison.json", ["bh_cutoff", "0.05", "largest_rejected_rank"], lambda v: v + 1),
    # A float field, so the tolerance branch is exercised and not only the integer equality
    # branch. verify_claims derives its tolerance from repr(), which is weaker for values
    # ending in zero, and that arm went untested until this line existed.
    ("shift_attribution.json", ["per_covariate_ks", "degree", "ks"], lambda v: v + 0.05),
]


def _dig(d, path):
    for k in path[:-1]:
        d = d[k]
    return d, path[-1]


def _run_checker() -> int:
    return subprocess.run([sys.executable, str(ROOT / "scripts" / "verify_claims.py")],
                          capture_output=True, text=True).returncode


def main() -> int:
    dirty = subprocess.run(["git", "status", "--porcelain", str(PROC)],
                           cwd=ROOT, capture_output=True, text=True).stdout.strip()
    if dirty:
        print("  refusing to run: data/processed has uncommitted changes")
        print("  commit or stash them first, so a failed run cannot lose work")
        return 1

    if _run_checker() != 0:
        print("  verify_claims already fails on clean data; fix that first")
        return 1
    print("  baseline: verify_claims passes on clean artifacts")

    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        for fname, path, breaker in MUTATIONS:
            src = PROC / fname
            backup = Path(tmp) / fname
            shutil.copy2(src, backup)
            try:
                obj = json.loads(src.read_text())
                parent, key = _dig(obj, path)
                original = parent[key]
                parent[key] = breaker(original)
                src.write_text(json.dumps(obj, indent=2))

                rc = _run_checker()
                where = "/".join(path)
                if rc == 0:
                    print(f"  {fname}:{where}  {original} -> {parent[key]}   NOT CAUGHT")
                    failures.append(f"{fname}:{where}")
                else:
                    print(f"  {fname}:{where}  {original} -> {parent[key]}   caught")
            finally:
                shutil.copy2(backup, src)

    if _run_checker() != 0:
        print("  ERROR: artifacts not restored cleanly; run `git checkout data/processed`")
        return 1

    print()
    if failures:
        print("FAIL: the harness did not notice " + "; ".join(failures))
        return 1
    print(f"the harness caught all {len(MUTATIONS)} injected errors and the artifacts are restored")
    return 0


if __name__ == "__main__":
    sys.exit(main())
