# Source data

The photovoltaic device records were taken from Choi et al., *Solar Energy* 303 (2026) 114112 ([DOI](https://doi.org/10.1016/j.solener.2025.114112)) and the accompanying [Duke data deposit](https://doi.org/10.7924/r40g3xh1d).

The analysis uses the following input files:

```text
data/source_v2/paper1_X_author_faithful.csv
data/source_v2/paper1_y.csv
data/source_v2/paper1_metadata.csv
data/source_v2/spacer_master.csv
data/source_v2/spacer_families.csv
data/source_v2/feature_groups.yaml
data/source_v2/paper1_matrix_manifest.json
data/source_v2/molecular_representation_manifest.json
data/external/pubchem_amine_pool.csv
```

The original files include device variables and PCE outcomes, DOI metadata, ionic spacer structures, molecular-family assignments and a snapshot of the PubChem retrieval. Several are curated or reformatted from their upstream sources, and the complete set is stored outside this public repository pending redistribution permissions.

The [build manifest](../results/clean_study1/manifest.json) lists the SHA-256 hashes for all nine inputs. The raw PubChem search can be repeated with [download_pubchem_library.py](../scripts/download_pubchem_library.py), although a later retrieval may contain different compounds.
