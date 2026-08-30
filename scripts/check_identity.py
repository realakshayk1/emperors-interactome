"""Check Proposition 3 by brute force: the largest set e-BH can certify over the monotone
budgeted class equals the number BH rejects, on every design and every rank configuration.

The paper proves this, so the script is a guard against the proof and the code drifting apart,
not evidence for the claim. It sweeps designs whose calibration-to-test ratio spans four orders
of magnitude, since the ratio is the quantity one might wrongly suppose the result depends on.

Proposition 2(iv) is what makes the sweep finite: a threshold calibrator attains the budget
bound at its own index, so it certifies whenever any admissible calibrator with the same
down-set does, and no t beyond q(n_cal+1) can certify anything.
"""
from __future__ import annotations

import math
import random
import sys

import numpy as np


def bh_count(ranks: np.ndarray, m: int, n1: int, q: float) -> int:
    p = np.sort(ranks) / n1
    ok = np.nonzero(p <= np.arange(1, m + 1) * q / m)[0]
    return int(ok[-1] + 1) if ok.size else 0


def bh_set(ranks: np.ndarray, m: int, n1: int, q: float) -> set:
    """The rejected set itself. Proposition 3 claims set equality, so counts are too weak."""
    p = np.sort(ranks) / n1
    ok = np.nonzero(p <= np.arange(1, m + 1) * q / m)[0]
    if not ok.size:
        return set()
    return set(np.nonzero(ranks / n1 <= p[int(ok[-1])])[0].tolist())


def class_max_set(ranks: np.ndarray, m: int, n1: int, q: float) -> set:
    best = set()
    for t in range(1, min(n1, int(math.floor(q * n1)) + 1) + 1):
        n_t = int((ranks <= t).sum())
        if n_t and n_t >= math.ceil(m * t / (q * n1)) and n_t > len(best):
            best = set(np.nonzero(ranks <= t)[0].tolist())
    return best


def class_max(ranks: np.ndarray, m: int, n1: int, q: float) -> int:
    best = 0
    for t in range(1, min(n1, int(math.floor(q * n1)) + 1) + 1):
        n_t = int((ranks <= t).sum())
        if n_t >= math.ceil(m * t / (q * n1)):
            best = max(best, n_t)
    return best


def main() -> int:
    rng = random.Random(0)
    bad = 0
    trials = 0

    for _ in range(2000):
        n1 = rng.randint(5, 300)
        m = rng.randint(5, 600)
        q = rng.choice([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5])
        ranks = np.array([rng.randint(1, n1) for _ in range(m)])
        trials += 1
        if class_max(ranks, m, n1, q) != bh_count(ranks, m, n1, q):
            bad += 1
    print(f"  random designs                 {trials - bad}/{trials} agree")

    # the proposition claims the SETS coincide, which equal counts would not establish
    set_bad = 0
    for _ in range(800):
        n1 = rng.randint(5, 120)
        m = rng.randint(5, 200)
        q = rng.choice([0.05, 0.1, 0.2, 0.3])
        ranks = np.array([rng.randint(1, n1) for _ in range(m)])
        if bh_set(ranks, m, n1, q) != class_max_set(ranks, m, n1, q):
            set_bad += 1
    print(f"  set equality (not just counts)  {800 - set_bad}/800 agree")
    bad += set_bad

    for ratio in (0.002, 0.01, 0.03, 0.543, 10, 200):
        sub = 0
        for _ in range(200):
            m = rng.randint(50, 400)
            n1 = max(3, int(m * ratio))
            q = rng.choice([0.05, 0.1, 0.2])
            ranks = np.array([rng.randint(1, n1) for _ in range(m)])
            trials += 1
            if class_max(ranks, m, n1, q) != bh_count(ranks, m, n1, q):
                sub += 1
                bad += 1
        print(f"  (n_cal+1)/m = {ratio:<7}          {200 - sub}/200 agree")

    # Dropping monotonicity only helps where the rank distribution is lumpy. Under the exact
    # uniformity that a.s.-distinct scores would give -- the premise under which the monotone
    # restriction stops being free -- the unconstrained class certifies nothing.
    # All three levels the manuscript quotes, across independent seeds. Running only two
    # levels left the third claim with no code path, and the number that reached the paper was
    # reproducible only by continuing this script's RNG stream past the earlier loops.
    M, mm = 905, 1666
    for q, ceiling in ((0.05, 0.005), (0.10, 0.005), (0.20, 0.040)):
        esc, draws = 0, 3000
        for seed in range(3):
            rng2 = np.random.default_rng(seed)
            for _ in range(draws // 3):
                occ = np.bincount(rng2.integers(1, M + 1, size=mm), minlength=M + 1)[1:]
                order = np.argsort(-occ)
                best = 0
                for a_ in range(1, M + 1):
                    n_a = int(occ[order[:a_]].sum())
                    if n_a and (M / a_) >= mm / (n_a * q):
                        best = max(best, n_a)
                        break
                esc += best > 0
        rate = esc / draws
        print(f"  exactly-uniform ranks, q={q}: class escapes in {esc}/{draws} ({100 * rate:.1f}%)")
        if rate > ceiling:
            bad += 1
            print(f"    escape rate {rate:.3f} exceeds the {ceiling:.3f} the manuscript implies")

    print()
    if bad:
        print(f"FAIL: {bad} of {trials} designs contradict Proposition 3")
        return 1
    print(f"Proposition 3 holds on all {trials} designs tested")
    return 0


if __name__ == "__main__":
    sys.exit(main())
