# Reproducibility and data access

This release includes analysis code, generated result tables and figures, candidate screening scores, and a manifest recording source SHA-256 hashes, model configuration, dependency versions and generated-output hashes.

It **does not** include the original device matrices, curated ionic registry or exact PubChem retrieval snapshot. A clean checkout can inspect the derived evidence but cannot rebuild the complete analysis without reconstructing those inputs.

## Install

Use Python 3.13 and the pinned dependencies in `pyproject.toml`:

```bash
python -m venv .venv
python -m pip install -e . pytest
```

Once the source files listed in [data/README.md](../data/README.md) are present:

```bash
python scripts/run_clean_study1.py --check
python scripts/run_clean_study1.py --rebuild
pytest -q tests/test_clean_study1.py tests/test_study1_diagnostics.py tests/test_ionization.py
```

`--check` validates that exact source/code/dependency fingerprint and checks the generated files against the manifest. `--rebuild` regenerates the entire results directory. Fresh PubChem retrievals may differ from the recorded pool and therefore are not a byte-for-byte reproduction.

## Inspection without the raw data

- [Evaluation](../results/clean_study1/benchmark.csv) and [held-out predictions](../results/clean_study1/heldout_predictions.csv).
- [DOI sensitivity](../results/clean_study1/positive_source_stability.csv), [family diagnostics](../results/clean_study1/family_level_diagnostics.csv) and [source-disjoint diagnostic](../results/clean_study1/source_disjoint_summary.csv).
- [Shortlist](../results/clean_study1/shortlist.csv), [candidate pool](../results/clean_study1/candidate_scores.csv) and [provenance manifest](../results/clean_study1/manifest.json).

Prospective experimental verification has not yet been performed.
