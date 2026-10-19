# Literature-guided prioritization of organic spacer cations for layered perovskite solar cells

Research code, retrospective evaluation and a prospective molecular shortlist.

Organic spacer selection in layered halide perovskites is usually studied one chemistry at a time. This project asks whether published photovoltaic device results contain enough information to **prioritize candidate spacer cations**, after accounting for variation in device composition, processing and publication source.

The analysis is exploratory. It is intended to guide a controlled experiment, not to predict a new molecule's power conversion efficiency (PCE) from its structure alone.

## At a glance

| Dataset | Evaluation | Screening |
| --- | --- | --- |
| 749 qualifying devices from 134 publications | Spacer-held-out and leave-one-family-out tests | 13,595 PubChem retrieval records |
| 42 verified cations; 31 monovalent | Fixed logistic and XGBoost baselines | 2,087 eligible monovalent candidates |
| 35 RDKit molecular descriptors | Source and family sensitivity analyses | 10 diversity-selected candidates |

### Retrospective evaluation

![Average precision under held-out chemistry and family evaluations](results/clean_study1/evaluation.png)

For the 31 monovalent spacers, the two-descriptor logistic model obtained **AP 0.420** in repeated chemistry holdouts and **AP 0.462** in leave-one-family-out evaluation, versus **0.285** and **0.292** for random ranking, respectively.

These values are **not unbiased prospective performance estimates**: model selection was informed by exploratory examination of the evaluation results. In addition, five of seven positive monovalent identities are supported by only one publication, and family-level performance is heterogeneous. The present evidence does not establish transferability to new molecular families or device processes.

[Evaluation table](results/clean_study1/benchmark.csv) · [Held-out predictions](results/clean_study1/heldout_predictions.csv) · [Source sensitivity](results/clean_study1/positive_source_stability.csv)

### Molecular screening

![Ten selected monovalent spacer cations](results/clean_study1/candidate_structures.png)

A fixed candidate-ionization and screening workflow selects ten structurally diverse monovalent cations. Their scores are **uncalibrated ranking scores, not device-success probabilities**. All ten candidates occupy an empirical *moderate-extrapolation* region rather than a validated applicability domain.

[Shortlist with PubChem IDs](results/clean_study1/shortlist.csv) · [Scored candidate pool](results/clean_study1/candidate_scores.csv) · [Selection counts](results/clean_study1/screening_counts.json)

## Approach

1. Curate ionic spacer identities and corresponding published device records.
2. Fit a DOI-aware device-context model; aggregate context-adjusted PCE residuals by molecular identity.
3. Derive training-fold-specific spacer labels; evaluate molecular descriptors with held-out chemistry and family partitions.
4. Apply the selected exploratory ranker to a filtered PubChem-derived monovalent ion pool.
5. Select a chemically diverse experimental shortlist subject to synthesis and applicability constraints.

The code preserves ionic charge states for measured spacer species and uses an explicitly defined amine-protonation rule for neutral candidate precursors. Model performance is assessed at the **molecule** level rather than treating devices fabricated with the same spacer as independent chemical examples.

See [methods and interpretation](docs/methods.md), [chemical representation](docs/chemical-representation.md) and the [prospective test plan](docs/experimental-validation.md).

## Reproducing the analysis

Validate the published artifact hashes without source-data access:

```bash
python scripts/verify_results.py
```

The published directory [`results/clean_study1/`](results/clean_study1/) contains model comparisons, source diagnostics, candidate scores, figures and an input/output fingerprint manifest.

```bash
# Python 3.13; use a clean virtual environment
pip install -e . pytest

# Requires the original, locally prepared source matrices and PubChem pool
python scripts/run_clean_study1.py --check
python scripts/run_clean_study1.py --rebuild
```

**Data availability:** the original device data are described in [Choi et al., *Solar Energy* 303 (2026) 114112](https://doi.org/10.1016/j.solener.2025.114112) and the associated [Duke repository](https://doi.org/10.7924/r40g3xh1d). This GitHub repository provides derived results and software, **not** the full source matrices, verified ionic registry or PubChem retrieval snapshot. Rebuilding requires reconstructing the exact source inputs listed in [data requirements](data/README.md); it cannot be run from a fresh clone alone.

The input/output SHA-256 values and dependency versions are recorded in [manifest.json](results/clean_study1/manifest.json). Current result files have not been independently validated by a new device-fabrication experiment.

## Repository guide

| Path | Contents |
| --- | --- |
| [`src/perovskite_prioritization/`](src/perovskite_prioritization/) | Ionic chemistry, feature construction, model evaluation, diagnostics |
| [`scripts/`](scripts/) | Analysis entry point and PubChem pool retrieval |
| [`results/clean_study1/`](results/clean_study1/) | Figures, validation tables, shortlist and provenance manifest |
| [`docs/`](docs/) | Methods, caveats and prospective validation |
| [`tests/`](tests/) | Unit and input-dependent integration checks |

## Status

**Retrospective computational study; experimental screening proposed, not completed.** The immediate next scientific test is a controlled, batch-replicated comparison of high-ranked candidates, matched low-ranked candidates and score-blind selections.

Released as research software under the [MIT License](LICENSE). Citation metadata are in [CITATION.cff](CITATION.cff).
