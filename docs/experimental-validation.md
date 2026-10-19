# Experimental validation of molecular prioritization

## Hypothesis and endpoint

The retrospective molecular ranker is intended to allocate a fixed experimental budget among candidate spacer chemistries. Its practical value is the ability to prioritize higher-performing materials under a controlled fabrication protocol, not to reproduce literature processing differences.

Primary endpoint: device-level reverse-scan PCE (%) summarized by *chemical identity* and fabrication batch, with a prespecified comparison against a known-spacer reference made in the same batch. Experimental hit and relative-performance criteria are registered before making the first device.

## Sample selection (prior to synthesis)

Freeze a versioned CSV containing candidate CID, precursor/cation SMILES, ranker score, model source commit and applicability-domain similarity. Hold chemical identities and model predictions fixed after preparation.

Choose **12 distinct, feasible new spacer precursors** from the same structurally and synthetically admissible PubChem pool:

- **Four high-ranked** candidates from the top 20% of the model score within the common eligibility conditions.
- **Four low-ranked or median-ranked, chemically matched controls** chosen by similar heavy-atom count, ionic valence, synthetic-accessibility band and nearest-spacer similarity where possible. These probe rank-order discrimination rather than trivial availability differences.
- **Four score-blind random candidates** sampled reproducibly from the *same eligibility pool*. Stratify by molecule size and scaffold family so random does not systematically select unmakeable or obviously different precursors.

Add **two known spacer references**, provisionally benzylammonium (BZA) and phenethylammonium (PEA), subject to confirmed compatibility with the lab's perovskite device architecture. The references are processed in every batch and are not counted among the 12 new chemical identities.

Choose a single PSC architecture, A-site/cation ratio, absorber composition, spacer concentration, solvent system and post-treatment schedule *before* study initiation. New spacer chemical identities must not be tuned separately to maximize their own PCE in the primary comparison. If a chemistry fails to form a usable device, record it as a failure rather than silently removing it from the denominator.

## Measurement blocks and replication

Prepare at least **three independent fabrication batches** on separate days. For each chemistry and reference, produce at least **six devices per batch**, distributed across two substrates or substrate positions where practical.

With 12 candidates + two reference spacers, the nominal primary design is **14 × 3 × 6 = 252 devices**. Randomize substrate layout within each batch and distribute candidate groups throughout the fabrication order. Where synthesis or equipment capacity is limited, a staged four-chemistry pilot should establish a reproducible reference, measurement noise and material feasibility **without unblinding model rankings**. The planned validation denominator remains prespecified.

Measure J–V PCE with an established scan rate, direction, light calibration and stabilized maximum power point where possible. Record hysteresis and device failures. Characterize representative films with XRD/GIWAXS (phase/orientation), steady PL and TRPL (recombination), AFM/SEM (morphology), and controlled stability testing when available.

## Statistical analysis

Use the *molecule* as the unit of validation, not each device as an independent molecular replicate.

- Primary statistic: Spearman correlation between frozen molecular scores and the **batch-adjusted median PCE** of the 12 new chemistries. Report the rank correlation with permutation-based uncertainty.
- Comparison: the difference in median batch-adjusted PCE between four high-ranked and four score-matched low/median candidates, evaluated by an exact small-sample/permutation procedure at the chemistry level. Random candidates provide a score-blind decision baseline.
- Screening-efficiency endpoint: for a fixed budget of four new compounds, compare hit recovery by the ranker against draws of four from the same eligible pool. A hit is defined before experiment, using a specified improvement margin over same-batch BZA/PEA reference and reproducibility across at least two batches.
- Sensitivity: mixed-effects device-level model with batch and substrate random effects, separated from the chemistry-level primary analysis; do not infer 252 independent observations for the molecule ranking.
- Missingness: report devices that fail or are unmeasurable and chemistry availability/synthesis failures, with all candidate identities kept in an intention-to-test table.

If practical, an external lab that has not participated in model fitting should execute fabrication and measurement from sealed instructions. The scorer and experimental analyst should not see outcome-driven candidate reorderings before the endpoint file is locked.

## What would establish a useful result?

A useful molecular prioritization method should show consistent improvement over equal-budget random selection and chemically matched controls under a single standardized process. If the model performs comparably to simple logP/carbon baselines, that would indicate that costly multifeature models do not add practical value in this dataset. If the effects disappear on previously unseen families, limit deployment to the supported chemical domain.

Structural observations are measured alongside PCE to distinguish chemistry-associated film formation and defect/transport changes; SHAP attributions alone do not identify those mechanisms.

This study should be repeated in a second material/process regime if a claim of broad perovskite-spacer transfer is intended. The PSC experiments in this document are planned rather than completed.
