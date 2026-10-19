# Methods

## Literature cohort

The source dataset contains 792 photovoltaic device records. We retained 749 records with one identified organic spacer, representing 42 spacer cations across 134 publications. The molecular-level dataset includes 31 monovalent and 11 divalent spacer identities.

Published device performance depends on absorber composition, architecture and processing. We therefore use a context-adjusted residual for the literature analysis, with the spacer molecule as the unit of ranking.

## Context-adjusted labels

A standardized ridge regression model (alpha = 10) is fitted to composition, device-architecture and processing features. DOI-grouped out-of-fold predictions are subtracted from measured PCE to obtain device residuals. Their median is calculated for each spacer; molecules above the training-cohort 80th percentile are labeled positive.

Within each held-out chemistry evaluation, the context model and residual threshold are refitted using the training subset. The residual captures associations with reported device performance after adjustment for available covariates. Processing differences that coincide with spacer identity may still contribute to the observed labels.

## Molecular descriptors and evaluation

The training structures retain their reported ionic charge states. RDKit provides 35 ordered molecular descriptors. We compared logistic regression with two and 35 descriptors and an XGBoost model with 35 descriptors, alongside a constant-score baseline and seeded random rankings.

Two validation settings were used:

- Ten repetitions of five-fold splits by spacer identity.
- Leave-one-molecular-family-out evaluation across eight monovalent families.

Average precision (AP), AUROC and top-20% enrichment are recorded in [benchmark.csv](../results/clean_study1/benchmark.csv), with individual test predictions in [heldout_predictions.csv](../results/clean_study1/heldout_predictions.csv).

The two-descriptor logistic model uses MolLogP and carbon-heavy-atom fraction. It was selected after exploratory comparisons of the held-out results. Its measured AP is consequently an exploratory estimate. Prospective screening will test the selected model on new compounds.

## Publication and family sensitivity

Five of the seven positive monovalent spacer identities are supported by a single publication. For these molecules, deleting the only DOI leaves no independent observation with which to recompute the label. For multiply reported spacers, publication deletion recalculates both the context model and the threshold.

A separate source-disjoint analysis removes publications that contain held-out spacer identities from each training set. This procedure also reduces the available training population. Family-holdout diagnostics are included because only two of the eight monovalent families contain both positive and negative labels under the fold-specific targets.

[Publication-deletion results](../results/clean_study1/positive_source_stability.csv) · [Source-disjoint summary](../results/clean_study1/source_disjoint_summary.csv) · [Family diagnostics](../results/clean_study1/family_spacer_diagnostics.csv)

## Molecular screening

A PubChem retrieval produced 13,595 records. The candidate-ionization procedure applies an amine-site protonation rule to neutral precursors and retains compatible molecular species. Of 3,296 eligible ions, 2,087 are monovalent and enter the main screen.

Thirty bootstrap fits of the selected two-descriptor logistic model provide mean candidate scores and score variation. The selection procedure combines the score with estimated synthetic accessibility, descriptor/fingerprint similarity to measured spacers and structural novelty:

| Screen | Molecules |
| --- | ---: |
| Monovalent candidates | 2,087 |
| Top score quintile | 418 |
| Non-strong-extrapolation domain | 715 |
| Score, domain and accessibility intersection | 121 |
| Pareto selection | 74 |
| Diversity shortlist | 10 |

The score and domain filters are evaluated on the same 2,087-compound pool; their individual counts should not be read as successive exclusions. Model scores have not been calibrated to device success, and the domain thresholds have not been tested against prospective prediction errors. Final ranks indicate diversity-selection order.

[Candidate score table](../results/clean_study1/candidate_scores.csv) · [Shortlist](../results/clean_study1/shortlist.csv)

## Experimental test

The planned experiment evaluates twelve new compounds: four high-ranked, four matched middle/low-ranked and four score-blind selections. Two known spacer references are included in every fabrication batch. Batch-adjusted PCE, rank correlation and enrichment at a fixed selection budget are the primary comparisons.

The complete experimental design is in [experimental-validation.md](experimental-validation.md). See [reproducibility.md](reproducibility.md) for the source files and analysis environment.
