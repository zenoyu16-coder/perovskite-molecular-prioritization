# Reproducing and checking the research

The public repository contains analysis code, frozen predictions and screening tables, molecular coordinates, energy tables and figures. A public checkout supports numerical and provenance checks without fitting ranking models or running xTB/DFT. Rebuilding the literature analysis additionally needs the exact source package in [data/README.md](../data/README.md).

## Environment

Use Python 3.13. Scientific package versions are pinned in [pyproject.toml](../pyproject.toml) and [requirements.txt](../requirements.txt); the test extra pins pytest. Installation from the repository root:

```bash
python -m pip install -e ".[test]"
```

An optional Conda environment is defined in [environment.yml](../environment.yml). CREST 3.0.2, GFN2-xTB 6.7.1 and PySCF 2.11.0 are separate calculation dependencies, not required for inspecting published energies. A complete Study1 rebuild also calls `pdflatex` twice to produce the bundled internal research report; install the LaTeX packages used by [manuscript.tex](../results/clean_study1/manuscript.tex).

## Public-checkout checks: no model training or quantum jobs

```bash
python scripts/verify_results.py
python scripts/verify_repository.py
python scripts/run_iodide_reference.py --check-only
python -m pytest -q
```

`verify_results.py` checks the frozen Study1 output hashes. `verify_repository.py` checks local Markdown/image links, README numbers against the CSV/JSON evidence, figure readability and package/citation metadata. It also recalculates AP from stored molecular predictions and electronic-energy differences from stored totals. Neither verifier fits a model.

The tests cover ionic representation, frozen selection behavior, molecular/source diagnostics, conformer records and ion-pair reference arithmetic. Source-dependent Study1 tests skip when the nonredistributed matrices are absent; pytest reports that distinction. A successful public-checkout suite does not imply that the full source build was run.

## Source-locked literature analysis

The [manifest](../results/clean_study1/manifest.json) records nine input hashes, six core analysis-code hashes, ordered descriptors, model parameters and package versions. With all nine matching files in their documented locations:

```bash
python scripts/run_clean_study1.py --check
python -m pytest -q tests/test_clean_study1.py
```

`--check` compares source/configuration/software and output fingerprints; it does not train models. Changes to inputs, core code or the recorded environment produce a stale-output error. README and documentation edits do not alter the frozen scientific fingerprint. `.gitattributes` fixes Python source files to LF line endings so a Windows checkout retains the recorded code SHA-256 values.

The following command launches the full analysis and replaces `results/clean_study1/`. Use a separate checkout when deliberately rebuilding:

```bash
python scripts/run_clean_study1.py --rebuild
```

It regenerates the complete dataset-derived target construction, holdout fits, screening outputs, diagnostics and internal report. A later PubChem query cannot substitute for the exact13,595-record snapshot. No rebuild is needed to inspect the released scientific results.

## Quantum records and replay

The [conformer study](quantum_chemistry.md) and [iodide association study](iodide_association.md) specify charges, continuum models, optimized geometries and single-point settings. Their coordinate and total-energy files are public. Raw external-program logs and intermediate scratch files are stored separately. Tests verify the published arithmetic and matching coordinates; they do not regenerate those logs or establish exhaustive conformer sampling.

All quantum scripts use `STUDY1_QM_WORKDIR`, defaulting to the ignored `data/external/study1_qm/`. Use a writable absolute path in Linux/WSL or Windows and keep it consistent across preparation, calculation and analysis. `STUDY1_XTB_BIN` selects the executable; otherwise xTB is located on PATH. The calculation commands in the two quantum documents launch new work and can overwrite corresponding derived outputs. They were not run during repository curation.

## Repository map and provenance boundaries

| Location | Contents |
| --- | --- |
| `src/perovskite_prioritization/` | Ionic representation, descriptors, context labels, ranking, screening and diagnostics |
| `scripts/` | Analysis entry point, artifact checks, PubChem acquisition and quantum replay |
| `tests/` | Public-artifact and source-dependent verification |
| `results/clean_study1/` | Frozen molecular predictions, benchmark, screening and internal generated report |
| `results/quantum_pilot/` | Isolated-cation conformers and fixed-geometry DFT energies |
| `results/iodide_pilot/` | Counterion-pair coordinates, references and energy comparisons |
| `docs/` | Methods, representations, calculation settings and proposed experimental protocol |

The frozen internal report in `results/clean_study1/` belongs to this project; it is not the published Choi paper. Current documentation specifies the exact scope of its general references to processing context: transport-layer deposition methods are available; detailed perovskite recipes are absent. The [MIT license](../LICENSE) covers source code and scripts. Third-party data retain their own provenance and permissions. [CITATION.cff](../CITATION.cff) identifies the software; cite the original device paper and dataset when using the literature-derived outputs.
