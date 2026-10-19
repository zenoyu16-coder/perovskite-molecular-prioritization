# Input data and availability

The device literature dataset originates from Choi et al., *Solar Energy* 303 (2026) 114112, DOI [10.1016/j.solener.2025.114112](https://doi.org/10.1016/j.solener.2025.114112), and the [associated Duke research-data deposit](https://doi.org/10.7924/r40g3xh1d).

The full reconstructed source matrices, curated ionic registry and original PubChem retrieval snapshot are not redistributed in this repository. The rebuild script expects these relative paths:

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

Several files contain curated or reformatted information rather than a byte-for-byte copy of the upstream publication deposit. The input SHA-256 hashes and exact build configuration are in [results/clean_study1/manifest.json](../results/clean_study1/manifest.json).

The `scripts/download_pubchem_library.py` utility can retrieve a new pool, but PubChem search results may change over time. A new snapshot cannot be substituted for the recorded source without treating it as a new analysis.
