"""§A identifying experiment: distribution-free control is present IFF the null
is exchangeable. Asserts the ON->OFF toggle as delta (null non-exchangeability)
increases, everything else held fixed."""
from emperor.identifying import run


def test_control_on_when_exchangeable_off_when_not():
    o = run(n_splits=150, deltas=[0.0, 1.0, 2.0])
    qp = o["q_primary"]
    by_delta = {r["delta"]: r for r in o["rows"]}
    # ON: exchangeable null (delta=0) controls FDR at/below q (+ small MC slack)
    assert by_delta[0.0][f"fdr@{qp}"] <= qp + 0.03, by_delta[0.0]
    assert by_delta[0.0]["control_holds@primary"]
    # OFF: non-exchangeable null blows the guarantee, monotonically
    assert by_delta[1.0][f"fdr@{qp}"] > qp + 0.1
    assert by_delta[2.0][f"fdr@{qp}"] > by_delta[1.0][f"fdr@{qp}"] > by_delta[0.0][f"fdr@{qp}"]
    assert not by_delta[2.0]["control_holds@primary"]
