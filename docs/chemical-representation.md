# Chemical representation

The 42 literature spacers are represented using their reported ionic structures. The dataset includes monovalent and divalent cations, with charged nitrogen environments spanning alkylammonium, aromatic and other ionic groups. RDKit descriptors and Morgan fingerprints are calculated from the charged structures.

For screening, PubChem supplies neutral precursor SMILES. The `amine_only_v2` rule identifies eligible amines and generates corresponding cationic structures. Aromatic ring nitrogens and resonance-deactivated sites, such as sulfonamides, are excluded from automatic amine protonation.

The main candidate pool is restricted to singly charged cations. Their predicted ranking, synthetic accessibility and similarity to previously measured spacers are calculated on the same ionic representation.

The generated ionic SMILES specify the structures used for modelling. Precursor salts, solution speciation and incorporation into a layered perovskite remain subjects for experimental characterization.

[Training descriptors](../results/clean_study1/ionic_training_descriptors.csv) · [Model configuration](../results/clean_study1/manifest.json)
