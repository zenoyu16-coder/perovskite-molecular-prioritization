# Reproducing the analysis

The source code, model evaluation tables, candidate scores and figures are included in this repository. The [results manifest](../results/clean_study1/manifest.json) records the input SHA-256 hashes, feature order, model settings, library versions and hashes of the generated outputs.

## Inspecting released results

With Python 3.13 and the packages listed in `pyproject.toml`:

```bash
pip install -e . pytest
python scripts/verify_results.py
pytest -q tests/test_public_release.py tests/test_ionization.py tests/test_study1_diagnostics.py
```

The verification script checks the released figures and tables against the recorded hashes. The input-independent tests check ion representation, model-selection behavior and stored diagnostics.

The main results are available without rerunning the analysis:

- [Model comparison](../results/clean_study1/benchmark.csv) and [held-out predictions](../results/clean_study1/heldout_predictions.csv)
- [Publication sensitivity](../results/clean_study1/positive_source_stability.csv) and [family diagnostics](../results/clean_study1/family_level_diagnostics.csv)
- [Ten candidates](../results/clean_study1/shortlist.csv) and [full candidate scores](../results/clean_study1/candidate_scores.csv)

## Quantum-chemical calculations

The [molecular conformer study](quantum_chemistry.md) and [iodide association study](iodide_association.md) include optimized coordinate files, input identifiers, calculation settings and electronic-energy tables.

The iodide reference can be checked without xTB:

```bash
python scripts/run_iodide_reference.py --check-only
pytest -q tests/test_iodide_pilot.py
```

The tests recompute ionic association energies and DFT double differences from published total energies. The full reference-state calculation requires xTB 6.7.1. See the ion-pair methods for the steps to regenerate it and the raw PySCF energy records.

## Rebuilding from source data

The full build also requires the device matrices, ionic structure registry and PubChem retrieval snapshot listed in [data/README.md](../data/README.md). These files are stored separately from the public code.

With the recorded inputs in place:

```bash
python scripts/run_clean_study1.py --check
python scripts/run_clean_study1.py --rebuild
pytest -q tests/test_clean_study1.py
```

`--check` compares the installed source files, code, environment and saved outputs with the manifest. `--rebuild` regenerates the results directory. A new PubChem retrieval or a different device-data release should be treated as a new input dataset; it may change the compound pool and model outputs.
