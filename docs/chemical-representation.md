# Ionic chemical representation

The molecular rankers operate on descriptors of **explicitly charged spacer species**. The training registry retains the reported ionic identities, including monovalent and divalent cations; uniformly protonating those structures would change their chemistry.

Neutral PubChem precursors are converted to candidate ions through the fixed `amine_only_v2` rule. It selects eligible amine sites while avoiding indiscriminate protonation of aromatic or resonance-deactivated nitrogen environments. The main virtual screening is restricted to singly charged candidates.

This transformation is a reproducible cheminformatics convention, **not** experimental evidence of solution speciation, precursor salt chemistry or the crystallized ionic state of any proposed molecule.

[Training ion descriptors](../results/clean_study1/ionic_training_descriptors.csv) · [Build manifest](../results/clean_study1/manifest.json)
