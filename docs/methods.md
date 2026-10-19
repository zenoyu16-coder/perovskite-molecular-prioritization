# Methods and interpretation

## Cohort and prediction unit

The source contains 792 reported device entries. After restricting to a single verified spacer per record, 749 devices from 134 publications and 42 molecular spacer identities remain. Of these identities, 31 are monovalent and form the main prospective candidate-screening cohort.

A device record is not an independent molecular-chemistry replicate. Targets, model splits and experimental decisions operate at the level of **spacer identity**; repeated devices and common DOIs create correlated evidence.

## Context-adjusted chemistry labels

A standardized ridge model (alpha = 10) uses a predefined set of composition, device architecture and processing features. DOI-grouped out-of-fold predictions provide PCE residuals, which are aggregated as a median by spacer. The upper quintile defines an exploratory high-residual label.

Context fitting and label thresholds are reconstructed inside training-chemistry partitions. The resulting residual associations are **not causal effects** of spacer substitution; publication-level chemistry and processing effects cannot be fully separated in retrospective data.

## Molecular representations and validation

All models use fixed ordered RDKit descriptors derived from charged ionic SMILES. The evaluated rankers are logistic regression with two or 35 descriptors and XGBoost with 35 descriptors, against constant-score and seeded random-order baselines.

Model evaluation comprises ten repeated five-fold held-out chemistry partitions, plus leave-one-family-out evaluation over eight monovalent families. Average precision is reported for pooled test predictions. Pooled family AP mixes scores from separately fitted models and cannot be read as eight independent within-family experiments.

The two-descriptor logistic model was selected after inspecting exploratory validation performance. Its reported AP is therefore **not independent of model selection**, and requires prospective validation.

## Literature-source sensitivity

Five of seven positive monovalent molecular identities occur in only one DOI. Removing that DOI removes direct evidence for the chemical identity. A separate source-disjoint test removes publications shared by held-out chemistries; this also changes training composition and label support. Only two monovalent held-out families contain both classes, making general claims about within-family transfer premature.

## Candidate selection

13,595 PubChem records were retrieved; 2,087 eligible monovalent ionic candidates remain after chemical filtering. A fixed logistic ranker is bootstrapped 30 times. A combined priority incorporates model-score ranking, model variability, synthetic accessibility, empirical domain proximity and novelty; a diversity selection produces ten candidates.

The score is an **uncalibrated ranking output**, not a probability of high-PCE fabrication. Final shortlist ordering is diversity-selection order rather than descending score. All ten candidates fall in the current empirical moderate-extrapolation category, whose thresholds have not been prospectively calibrated.

## Next experiment

Prespecify identities and ranking scores. Test higher-ranked compounds against matched lower-ranked and score-blind controls under one fixed fabrication recipe, randomized by device position and manufacturing batch. Analyze the molecule, not individual devices, as the independent ranking unit.

See [experimental validation](experimental-validation.md) and [reproducibility](reproducibility.md).
