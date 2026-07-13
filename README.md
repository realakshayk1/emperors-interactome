# The Emperor's Interactome

**When can you audit an AI protein interactome? A distribution-free reliability audit of a published AlphaFold-Multimer map, refereed by held-out biology.**

---

## Summary

The Krogan/Ideker CM4AI multimodal cell map (Gladstone Institutes, *Nature* 2025) is among the few published AI interactomes that can be audited externally, because it ships a raw confidence axis (AlphaFold-Multimer ipTM) alongside a native random-pair null — the exchangeable calibration set that honest, distribution-free error control requires. That matters because a "high-confidence" cutoff calibrated on a class-balanced benchmark hides a much larger false-discovery rate in the real regime of a genome-scale screen, where true interactions are rare — and we measure exactly that gap.

Applying **distribution-free conformal FDR control** to the map's published high-confidence tier, **35 of 161 edges (22%) fail honest error control at q=0.10**. The audit's direction is confirmed by an independent referee the structural model never saw: **certified edges are 41% co-essential (DepMap) versus 17% for dropped edges** (permutation p=0.016), under a strict purity firewall that keeps the referee out of every score, label, and calibration set. We then interrogate the guarantee itself — measuring the calibration-to-candidate exchangeability violation and converting it into a conditional robustness statement — and validate the certified set against orthogonal wet-lab physical interactions and by held-out complex-member recovery.

**Full paper:** [`reports/MANUSCRIPT.md`](reports/MANUSCRIPT.md). **Short summary:** [`reports/summary.md`](reports/summary.md).

## Reproduce

```bash
make reproduce    # raw → idmap → labels → conformal audit → validate → nominate → figures
make audit-self   # the full hardened analysis suite (dependence, shift, sensitivity, second map, recovery)
make test         # unit + adversarial firewall + held-out-recovery tests
```

The audit is reproducible from a bare clone (interim tables are regenerated from raw data; the second-map score table is auto-fetched from its public source).

## Data and configuration

- **Primary map**: Krogan/Ideker CM4AI, Schaffer et al., *Nature* 2025 (doi:10.1038/s41586-025-08878-3); AlphaFold-Multimer ipTM as the confidence axis, with 1,788 native random-pair decoys as calibration negatives.
- **Positive labels**: CORUM human complexes, complex-disjoint calibration/test split.
- **Held-out referee**: DepMap co-essentiality (Wainberg et al. 2021) — strict purity firewall.
- **Target FDR**: q = 0.10 (reported across a q ∈ {0.05, 0.10, 0.20} sweep).

## Repository

- `reports/MANUSCRIPT.md` — the full paper.
- `reports/summary.md` — short summary.
- `src/emperor/` — the pipeline (`conformal.py`, `validate.py`, `nominate.py`, and the hardening modules).
- `data/processed/` — every result JSON; each numeric claim in the paper traces to one.
- `results/figures/` — generated figures.
- `tests/` — firewall, dependence, and held-out-recovery tests.

Licensed MIT.