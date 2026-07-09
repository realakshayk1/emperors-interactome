# The Emperor's Interactome — the 22% that the audit cannot certify

## Headline

Of the 161 high-confidence edges in the CM4AI multimodal cell map (Schaffer et al., *Nature* 2025; doi:10.1038/s41586-025-08878-3), **35 (22%) are not supported at q=0.10 under distribution-free conformal calibration** audited on raw AlphaFold ipTM/Score alone. The number sweeps with the target FDR: **83 not supported at q=0.05, 35 at q=0.10, 12 at q=0.20.** These are edges whose conformal p-value is too large for Benjamini–Hochberg to certify at the stated error rate — the audit withholds a certificate, it does not label the interaction false.

## Notable named casualties (dropped at q=0.10)

The seven dropped edges whose two proteins co-occur in the *same* annotated CORUM complex are the most conspicuous, because these are edges a database-driven pipeline would treat as "known-good":

- **PIH1D1–RPAP3** (R2TP complex, CORUM 8735; conf p=0.036) — the core HSP90 co-chaperone scaffold that assembles box C/D snoRNPs, PIKK kinases and the U5 snRNP.
- **DIS3–EXOSC5** (RNA Exosome, CORUM 788/789; conf p=0.017) — the catalytic ribonuclease subunit paired with a core barrel subunit of the 3′→5′ RNA-degradation exosome.
- **DDX23–TXNL4A** (Spliceosome pre-B complex, CORUM 8370; conf p=0.017) — two spliceosomal components brought together during tri-snRNP/pre-B assembly.
- **KANSL3–MCRS1** (NSL histone-acetyltransferase complex, CORUM 7221; conf p=0.017) — both NSL subunits; notable because KANSL3 is the *certified* nominee elsewhere in this audit (KANSL3→MLL1-WDR5), so the method is discriminating within the same protein's edge set rather than blanket-passing it.
- **FBRSL1–PCGF5** and **AUTS2–PCGF5** (PRC1.5 Polycomb complex, CORUM 8387/8388; conf p=0.013–0.012) — two edges into the non-canonical PRC1.5 transcriptional-repression module.
- **STX16–VAMP4** (Golgi/endosomal SNARE complexes, CORUM 875/877; conf p=0.013) — cognate SNAREs mediating endosome-to-TGN fusion.

Beyond these, the least-supported dropped edges include several signalling/adaptor pairs where only one partner has a CORUM annotation (e.g. **GIT2–PAK3**, conf p=0.041, with GIT2 in the ARHGEF7–GIT2–PAK1 module; **TAF4–TAF8**, conf p=0.034, both TFIID-associated but not co-annotated).

## Caveat

Held-out FDR here is an **upper bound**: CORUM is incomplete, so a genuine interaction that CORUM simply hasn't catalogued will also read as unsupported. "Not supported at q under distribution-free calibration" means the raw-metric evidence is insufficient to certify the edge at that error rate — it is a call for orthogonal validation (e.g. the held-out DepMap co-essentiality referee), not a determination that any of these complexes is wrong.
