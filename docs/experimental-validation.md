# Experimental validation plan

## Design

The proposed experiment tests whether literature-trained molecular scores help select spacer candidates under one fabrication protocol. The primary outcome is reverse-scan PCE for each spacer, summarized across fabrication batches. A same-batch reference is used to account for run-to-run variation.

The candidate list and scores will be fixed before fabrication. Twelve new spacer precursors will be selected from the screened, synthetically feasible pool:

| Group | Compounds | Selection |
| --- | ---: | --- |
| High-scoring | 4 | Top 20% of the model score among eligible candidates |
| Matched middle/low-scoring | 4 | Similar molecular size, charge, accessibility and training-set similarity |
| Score-blind controls | 4 | Randomly selected from the same eligible pool, stratified by size and scaffold |

Two established spacer references, provisionally benzylammonium and phenethylammonium, will be included after confirming compatibility with the selected perovskite device architecture.

## Fabrication and measurements

One PSC architecture, absorber composition, solvent system, spacer loading and post-treatment procedure will be used across the study. Candidate-specific recipe optimization will be reserved for follow-up work.

The initial design includes three fabrication batches on separate days and six devices per compound per batch, giving 14 compounds × 3 batches × 6 devices = 252 devices. Positions and fabrication order will be randomized within batches. Failed devices and compounds that cannot be processed will be recorded.

PCE will be measured from calibrated J–V scans, with scan direction, rate, hysteresis and stabilized operation documented. XRD or GIWAXS, steady-state PL, TRPL and microscopy will be used on representative films to examine phase formation, recombination and morphology. Stability measurements can follow the primary ranking test.

A smaller feasibility run may first establish device reproducibility and material availability. Its results will be kept separate from the prespecified ranking comparison.

## Analysis

The principal analyses use the twelve chemical identities as the ranking units:

- Spearman correlation between model score and batch-adjusted median PCE, with permutation-based uncertainty.
- Difference in batch-adjusted PCE between the four high-scoring candidates and four matched middle/low-scoring controls.
- Recovery of compounds exceeding a prespecified reference-relative PCE threshold when selecting four molecules by model score versus four from the eligible pool at random.

Batch and substrate effects will be examined with device-level mixed-effects models. All candidate identities and fabrication failures will remain in the analysis record.

The reference-relative performance threshold and minimum replication criteria will be set before device fabrication. The measurements described here are planned; no prospective device results are included in this repository.
