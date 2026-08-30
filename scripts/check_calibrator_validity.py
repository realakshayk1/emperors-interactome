"""Guard against claiming a rank calibrator is an e-value when it is not.

A numeric verifier cannot catch this: a calibrator that overspends its budget still produces
arithmetically correct rejection counts, and those counts reproduce exactly. The error is in
whether the object is an e-value at all.

Exchangeability constrains the rank only through P(R <= t) <= t/(n_cal+1). Under that
constraint the least favourable null puts its mass wherever the calibrator is largest from
that rank onward, so

    sup E[e(R)] = (1/(n_cal+1)) * sum_r g(r),   g(r) = max_{s >= r} e(s).

For non-increasing e we have g == e and the budget sum_r e(r) <= n_cal+1 is exactly validity.
Off the monotone class the budget is not sufficient: the calibrator supported on
{1,...,6} u {15} spends exactly n_cal+1 and still reaches E[e] = 15/7.

The same identity is why monotonicity costs nothing. g is non-increasing, is valid whenever e
is, and dominates e pointwise; e-BH's rejected set only grows when e-values grow, so g
certifies whatever e does. Every deterministic rank calibrator is therefore dominated by a
monotone one, which is what lets Proposition 3 range over all of them.
"""
from __future__ import annotations

import math
import sys

import numpy as np


def envelope(e: np.ndarray) -> np.ndarray:
    """Right-running maximum: the smallest non-increasing function dominating e."""
    return np.maximum.accumulate(e[::-1])[::-1]


def sup_expectation(e: np.ndarray, n1: int) -> float:
    """sup E[e(R)] over every null obeying P(R <= t) <= t/n1."""
    return float(envelope(e).sum() / n1)


def is_valid(e: np.ndarray, n1: int, tol: float = 1e-9) -> bool:
    return sup_expectation(e, n1) <= 1.0 + tol


def worst_case_null(e: np.ndarray, n1: int) -> np.ndarray:
    """A null attaining the supremum, for exhibiting a violation concretely."""
    g = envelope(e)
    p = np.zeros(e.size)
    # spend each unit of rank mass at the first rank whose envelope it attains
    for r in range(e.size):
        if g[r] == e[r]:
            p[r] += 1.0 / n1
        else:
            p[np.argmax(e[r:] == g[r]) + r] += 1.0 / n1
    return p


def threshold(n1: int, t: int) -> np.ndarray:
    e = np.zeros(n1)
    e[:t] = n1 / t
    return e


def harmonic(n1: int) -> np.ndarray:
    H = sum(1.0 / i for i in range(1, n1 + 1))
    return np.array([n1 / (r * H) for r in range(1, n1 + 1)])


def main() -> int:
    n1 = 905
    failures = []

    # every calibrator the paper reports must be a genuine e-value
    reported = {"harmonic": harmonic(n1)}
    reported.update({f"threshold t={t}": threshold(n1, t) for t in (1, 2, 3, 7, 8, 19, 20, 25)})
    for name, e in reported.items():
        s = sup_expectation(e, n1)
        ok = s <= 1 + 1e-9
        print(f"  {name:18s} sup E[e] = {s:.6f}  {'ok' if ok else 'INVALID'}")
        if not ok:
            failures.append(name)

    # the budget alone does not imply validity off the monotone class
    S = [1, 2, 3, 4, 5, 6, 15]
    e = np.zeros(n1)
    e[[r - 1 for r in S]] = n1 / len(S)
    budget, sup = e.sum() / n1, sup_expectation(e, n1)
    print(f"\n  non-monotone support {S}:")
    print(f"    budget sum e(r)/(n+1) = {budget:.6f}  (looks admissible)")
    print(f"    sup E[e]              = {sup:.6f}  = {15}/{7}")
    if not (abs(budget - 1.0) < 1e-9 and abs(sup - 15 / 7) < 1e-9):
        failures.append("budget-is-not-validity demonstration")
    q = worst_case_null(e, n1)
    cdf = np.cumsum(q)
    if not np.all(cdf <= np.arange(1, n1 + 1) / n1 + 1e-12):
        failures.append("worst-case null violates the dominance constraint")

    # monotonicity is free: the envelope is valid whenever the calibrator is, and dominates it
    rng = np.random.default_rng(0)
    for _ in range(2000):
        m = int(rng.integers(4, 60))
        a = rng.random(m) * rng.integers(1, 40)
        g = envelope(a)
        if not np.all(g >= a - 1e-12):
            failures.append("envelope fails to dominate")
            break
        if abs(sup_expectation(a, m) - sup_expectation(g, m)) > 1e-9:
            failures.append("envelope changes the supremum")
            break

    print()
    if failures:
        print("FAIL: " + "; ".join(sorted(set(failures))))
        return 1
    print(f"all {len(reported)} reported calibrators are valid e-values; "
          f"budget-vs-validity and envelope properties hold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
