# Methods

## Data cohorts and prediction unit

The source dataset contains 792 photovoltaic device records. Selecting records with one identified, structure-verified spacer retains **749 devices, 134 DOI and 42 spacer identities**. The main `mono31` cohort contains 31 monovalent identities represented by 647 devices. The supplementary `all42` cohort also includes 11 divalent identities. Comparisons in the README all use `mono31`.

Spacers are the ranking units. The outcomes are literature-derived molecular labels, not absolute PCE predictions for new devices. The original reported reverse-scan PCE values remain unchanged. Registry identities and formal charges define the cohorts; names alone do not establish new molecular identities.

## Publication-aware context adjustment

The context estimator is `StandardScaler` followed by `Ridge(alpha=10)`. Its 78 inputs are listed in [device_cohort.json](../results/clean_study1/device_cohort.json): composition coefficients and element/ion categories, nominal large-organic/metal ratio, device area and architecture, ETL/HTL stacks and their deposition methods. Spacer-identity columns and molecular descriptors are excluded. Detailed perovskite solvent, annealing and additive conditions are absent from this context matrix.

For a specified chemistry cohort, five-fold `GroupKFold` by DOI gives out-of-publication context predictions. The device residual is `reported PCE − predicted context PCE`; its median across a spacer's devices defines that spacer's performance representation. These associations can retain unrecorded recipe and source effects.

For an outer training/test split by spacer:

1. Keep only devices belonging to the training spacer identities. Fit context scalers and ridge models within their DOI cross-fitting partitions, then aggregate the training out-of-DOI residuals by spacer.
2. Set the positive threshold to the training residuals' empirical 80th percentile, using NumPy `method="higher"`. Values **at or above** the threshold are positive.
3. Fit a context model on all outer-training devices and apply it to held-out spacer devices. Their observed PCE residual medians are compared with the training threshold to construct evaluation labels.
4. Fit molecular scaling and the ranker only on training spacers and their training labels. Evaluate the held-out molecular scores against the fold-specific labels.

Held-out PCE observations define evaluation labels, not training features or the threshold. Full-cohort labels are used to stratify the initial repeated splits; they are not reused as the ranker's fold-specific training targets. Final screening labels are separately reconstructed from the full 31-spacer cohort. These two uses are recorded in [training_labels.csv](../results/clean_study1/training_labels.csv) and [outer_folds.csv](../results/clean_study1/outer_folds.csv).

## Descriptors, rankers and metrics

The charged ionic structures have 35 ordered RDKit descriptors, with exact order in the [manifest](../results/clean_study1/manifest.json). Logistic2 uses MolLogP and CarbonHeavyAtomFraction. Logistic35 uses all 35. Both use training-fold standardization and logistic regression with `C=1`, `class_weight="balanced"`, and `max_iter=1000`.

XGBoost35 uses 140 trees, maximum depth 2, learning rate 0.05, subsample 0.85, column subsample 0.8 and L2 regularization 2. These ranker parameters are fixed. The stored protocol name `nested` refers to target reconstruction inside the outer splits; there is no nested hyperparameter search.

The evaluation comprises ten repetitions of five-fold stratified spacer holdout and leave-one-family-out holdout across eight monovalent families. For repeated holdout, scores are pooled across five test folds to make one complete 31-spacer panel per repetition; the benchmark averages the ten panel metrics. Family metrics use one pooled 31-spacer panel from the eight excluded-family fits. They are not the mean of eight within-family AP values.

Reported metrics are AP, AUROC and top-20% enrichment. For enrichment, select `ceil(0.2 × n_spacers)` highest scores and divide their positive fraction by panel prevalence. Seeded random rankings supply 300 repeated-holdout panels (10 repetitions × 30 draws) and 30 family panels. Their finite-list AP exceeds prevalence on average; the random AP is not a theoretical prevalence line. The constant-score baseline is also retained in [benchmark.csv](../results/clean_study1/benchmark.csv).

The 35-descriptor logistic ranker has the highest monovalent repeated-holdout AP. Logistic2 was retained for its compact representation and exploratory performance across both settings after the outer results had been inspected. [model_selection.json](../results/clean_study1/model_selection.json) records this selection exposure. Frozen AP values are model comparisons, not an independent post-selection performance estimate.

## Publication and family sensitivity

DOI-aware context cross-fitting does not make the primary spacer holdout source-disjoint: one paper can report several spacers. The separate source-disjoint diagnostic removes every publication containing a held-out spacer from context/ranker training, recalculates labels and refits. It covers one five-fold spacer repetition and eight family folds. Its changed training population and labels prevent a direct interpretation as a same-target performance gain. The [source-disjoint summary](../results/clean_study1/source_disjoint_summary.csv) preserves both label sets and the resulting AP.

Five of the seven positive monovalent identities in the full-cohort labels have a single supporting DOI. Removing that source leaves no target observation; missing support is not imputed as a negative label. The two multiply reported positives are examined by DOI deletion with context and threshold refitting. [Publication-deletion results](../results/clean_study1/positive_source_stability.csv) record both label changes and unavailable cases.

Only two of eight excluded families contain both classes under the fold-specific labels. Pooled family AP therefore includes between-family ordering. [Family-level](../results/clean_study1/family_level_diagnostics.csv) and [spacer-level](../results/clean_study1/family_spacer_diagnostics.csv) diagnostics make this distinction inspectable.

## Virtual screening and candidate order

The PubChem snapshot contains 13,595 retrieval records. Chemical eligibility, precursor deduplication and the `amine_only_v2` protonation policy yield 3,296 candidate representations; 2,087 are monovalent and enter the main screen. See [chemical representation](chemical-representation.md) for the protonation rules.

Thirty class-stratified bootstrap fits of Logistic2 yield the mean classifier score and its variation. Candidate similarity uses radius-2, 256-bit Morgan fingerprints. Descriptor distances use training-standardized 35-descriptor space; a fingerprint PCA representation supplies a second distance. These inputs define an empirical applicability score and interpolation/moderate/strong-extrapolation categories. The moderate category requires maximum Tanimoto ≥ 0.30 and descriptor distance ≤ 2.5 times the reference distance; interpolation requires similarity ≥ 0.60 and a tighter distance criterion. They are empirical coverage rules, not calibrated prediction-error bounds.

| Gate or selection | Molecules |
| --- | ---: |
| Monovalent pool | 2,087 |
| Top score quintile | 418 |
| Non-strong-extrapolation domain | 715 |
| Score/domain/accessibility intersection, SA ≤ 5 | 121 |
| Pareto-optimal | 74 |
| Preliminary diversity selection | 20 |
| Final diversity shortlist | 10 |

The 418 and 715 counts are parallel gates on the same pool, not successive exclusions. Pareto objectives are model score, bootstrap variation, accessibility, applicability and novelty. Composite candidate priority uses weights 0.40/0.15/0.15/0.15/0.15; fingerprint clustering uses similarity 0.55. The final ten balance priority with distance from already selected molecules.

`score` is the bootstrap mean classifier output; `logistic2_score` is a single full-cohort fit. `candidate_priority` adds the other criteria. `final_rank` is diversity-selection order and need not match descending score. These distinct quantities are preserved in the [candidate table](../results/clean_study1/candidate_scores.csv) and [shortlist](../results/clean_study1/shortlist.csv).

## Calculation and experimental scope

The shortlisted rank-tied isomers motivate the [conformer comparison](quantum_chemistry.md) and [iodide association study](iodide_association.md). Those calculations concern isolated cations and one-counterion pairs. A separate [prospective device protocol](experimental-validation.md) describes a proposed ranking test; no fabricated-candidate PCE results are included. Source-data requirements and frozen-build checks are in [reproducibility.md](reproducibility.md).
