# Prospective Nomination Register (pre-registered 2026-07-10)

**Status:** frozen, timestamped. These are conformal-BH *certified* nominations (q=0.10, FDR ≤ 0.10 upper bound) ranked by ORTHOGONAL DepMap co-essentiality support (ranking only — DepMap never entered selection; purity firewall intact).

**Positive control (NOT a discovery claim):** KANSL3 → MLL1-WDR5 / NSL / MSL complexes (CORUM 5386/610/7221). The engine recovers KANSL3's known NSL-family membership — a sanity check that the nomination procedure finds true biology, demoted from 'flagship discovery' to labeled control per PLAN_V3 T3.2. KANSL3→MSL(610) carries the highest co-essentiality of any nomination (signed −log10p = 19.6).

## Top 15 prospective NOVEL nominations (KANSL3 excluded)

| rank | nominee | via member | complex | conf p | AF-M | co-ess (signed −log10p) |
|---|---|---|---|---|---|---|
| 1 | YPEL5 | WDR26 | CUL4B-DDB1-WDR26 complex (6490) | 0.0011 | 0.844 | 68.25 |
| 2 | MAEA | WDR26 | CUL4B-DDB1-WDR26 complex (6490) | 0.0055 | 0.566 | 49.67 |
| 3 | CHUK | IKBKB | IKBKB-MAP3K14 complex (2104) | 0.0022 | 0.626 | 31.80 |
| 4 | IKBKB | CHUK | TNF-alpha/NF-kappa B signaling complex 9 (5285) | 0.0022 | 0.626 | 31.80 |
| 5 | CHUK | IKBKB | BCL10-CARD11-DPP4-IKBKB complex (8868) | 0.0022 | 0.626 | 31.80 |
| 6 | YPEL5 | WDR26 | CTLH complex (8273) | 0.0011 | 0.844 | 25.59 |
| 7 | WDR77 | PRMT5 | Brm-associated complex (711) | 0.0011 | 0.898 | 20.53 |
| 8 | WDR77 | PRMT5 | BRG1-SIN3A complex (713) | 0.0011 | 0.898 | 20.53 |
| 9 | WDR77 | PRMT5 | BRM-SIN3A complex (714) | 0.0011 | 0.898 | 20.53 |
| 10 | WDR77 | PRMT5 | BRG1-SIN3A-HDAC containing SWI/SNF remodeling complex I (803) | 0.0011 | 0.898 | 20.53 |
| 11 | WDR77 | PRMT5 | BRM-SIN3A-HDAC complex (806) | 0.0011 | 0.898 | 20.53 |
| 12 | WDR77 | PRMT5 | BRG1-associated complex (807) | 0.0011 | 0.898 | 20.53 |
| 13 | WDR77 | PRMT5 | BRM-associated complex (808) | 0.0011 | 0.898 | 20.53 |
| 14 | WDR77 | PRMT5 | PRMT5 complex (2764) | 0.0011 | 0.898 | 20.53 |
| 15 | WDR77 | PRMT5 | Brg1-associated complex II (3063) | 0.0011 | 0.898 | 20.53 |

## What would confirm / refute each

- **Confirm:** (a) direct co-IP or proximity labeling of nominee↔via-member in a relevant cell line; (b) high-fidelity AF3/Boltz-2 re-fold with ipTM ≥ 0.6 and a consistent interface; (c) appearance as a co-complex member in a future CORUM/complex-portal release.
- **Refute:** (a) co-IP negative under conditions where the via-member's known partners are recovered; (b) AF3/Boltz-2 re-fold ipTM < 0.3 with no consistent interface; (c) nominee localizes to a compartment incompatible with the complex.

**Freeze basis:** CM4AI map (Schaffer et al. 2025), CORUM 5.3, DepMap co-essentiality (Wainberg 2021). Register frozen 2026-07-10; any nomination later found in a database predating this date is disqualified from the prospective set.

**Reproducibility:** ranking = `data/processed/nomination_sets.json` sorted by (coess desc, conf_pvalue asc), KANSL3 rows removed.