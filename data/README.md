# Source data and access

The device measurements originate from Choi et al., *Interpretable machine learning on a curated dataset identifies chemical descriptors governing 2D perovskite solar cell performance*, *Solar Energy* 303 (2026), 114112 ([paper DOI](https://doi.org/10.1016/j.solener.2025.114112)). The associated [Duke data deposit](https://doi.org/10.7924/r40g3xh1d) and [authors' code/data repository](https://github.com/migonchoi/ML-for-2D-Perovskite-Solar-Cell-Performance) provide the upstream device data and original device-level modeling context.

This project's molecular analysis starts from the released 792-row matrices and retains 749 single-identified-spacer devices, covering 134 DOI and 42 spacer identities. The primary 31-monovalent-spacer cohort contains 647 of those devices. These cohort counts have different denominators from the original paper's device-level PCE task.

## Exact build inputs

The following files are intentionally outside the tracked public release:

| Path | Role |
| --- | --- |
| `data/source_v2/paper1_X_author_faithful.csv` | Author-faithful 137-column device input matrix |
| `data/source_v2/paper1_y.csv` | Four unchanged reverse-scan device targets |
| `data/source_v2/paper1_metadata.csv` | Row-aligned record IDs, DOI, spacer and device metadata |
| `data/source_v2/spacer_master.csv` | Curated, structure-verified ionic identity registry |
| `data/source_v2/spacer_families.csv` | Curated molecular-family assignments |
| `data/source_v2/feature_groups.yaml` | Frozen device-context column groups |
| `data/source_v2/paper1_matrix_manifest.json` | Matrix provenance and author-encoding checks |
| `data/source_v2/molecular_representation_manifest.json` | Ionic-structure provenance |
| `data/external/pubchem_amine_pool.csv` | Frozen 13,595-record PubChem retrieval snapshot |

All nine SHA-256 values are recorded in [manifest.json](../results/clean_study1/manifest.json), under `provenance.files`. The matrices, registries and snapshot must match these hashes for a source-locked check. The upstream archive alone does not provide every curated derivative under these filenames. Obtain the matching curated input package from the project maintainer; this repository currently has no public download for that complete package and no automated reconstruction of every curated file.

Source data remain governed by their upstream terms. The Duke deposit is identified as CC BY 4.0 in the source audit; the repository's MIT code license does not replace third-party data terms. No new redistribution permission is asserted here.

## Available without the source package

[results/clean_study1](../results/clean_study1/) includes modeled ionic SMILES and descriptors, molecular training labels, fold-specific held-out predictions, source/family diagnostics, candidate scores and the ten-candidate shortlist. Coordinates and electronic-energy tables for the quantum case studies are also public. These support inspection and arithmetic checks even when the original build inputs are unavailable.

[download_pubchem_library.py](../scripts/download_pubchem_library.py) documents a repeatable retrieval procedure. A new download can change as PubChem changes; it is a new snapshot, not a replacement for the frozen input. See [reproducibility](../docs/reproducibility.md) for the separate public-verification and full-rebuild commands.
