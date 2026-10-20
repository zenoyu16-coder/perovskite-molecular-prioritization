# Quantum-chemical comparison of rank-tied spacer candidates

## Molecular pairs

The two-descriptor logistic ranker uses MolLogP and carbon-heavy-atom fraction. Both descriptors are insensitive to the position of a substituent on the aromatic ring. Two isomer pairs appear in the [ten-molecule shortlist](../results/clean_study1/shortlist.csv):

| Compounds (PubChem CIDs) | Substituent pattern | Logistic2 score |
| --- | --- | ---: |
| 2760937 / 16218335 | meta / ortho methylsulfonylbenzylammonium | 0.4875 / 0.4875 |
| 40786958 / 24721203 | Fluoro/methoxybenzylammonium regioisomers | 0.2809 / 0.2809 |

The shortlist's combined `score` is an ensemble average (0.4713 for the first pair and 0.2533 for the second pair); it differs slightly from the fitted Logistic2 column. We used the registered monovalent cation SMILES for every calculation.

## Conformer calculations

ETKDGv3 conformers were minimized with MMFF94s and subsequently optimized using GFN2-xTB 6.7.1. All 13 selected starting geometries converged in both the gas phase and ALPB(DMF), producing 26 optimization records. The screening used 90 embedding attempts per molecule and a 0.65 Å heavy-atom RMS filter on initial geometries.

For the two methylsulfonylbenzylammonium isomers, CREST 3.0.2 sampled additional conformations using quick iMTD-GC at GFN2-xTB/ALPB(DMF). The resulting low-energy ensembles show different intramolecular N–H···O contacts.

![CREST conformational energy and intramolecular N···O distance](../results/quantum_pilot/crest_isomer_conformations.png)

| CREST / ALPB(DMF) | meta | ortho |
| --- | ---: | ---: |
| Conformers retained by CREST | 11 | 5 |
| Conformers within 3 kcal/mol of the minimum | 11 | 4 |
| Of these, meeting the N–H···O criterion | 0 | 4 |
| Minimum-energy N···O separation | 5.35 Å | 2.61 Å |
| Minimum-energy H···O separation | 4.58 Å | 1.63 Å |

The contact criterion is N···O < 3.3 Å, H···O < 2.5 Å, and N–H···O angle ≥ 130°. The ortho cation favors a geometry with a local ammonium–sulfonyl oxygen contact; the meta cation does not access an analogous contact in the sampled low-energy ensemble.

The second isomer pair shows the same kind of geometric contrast in the initial xTB optimizations (minimum-energy N···O separations of 2.72 and 5.41 Å in ALPB(DMF)). CREST sampling has not yet been performed for that pair.

[CREST conformer data](../results/quantum_pilot/crest_conformers.csv) · [GFN2-xTB records](../results/quantum_pilot/xtb_conformers.csv) · [Optimized CREST ensembles](../results/quantum_pilot/crest_ensembles/)

## DFT single-point comparison

The lowest CREST structure of each methylsulfonyl isomer was evaluated using PySCF 2.11.0, density-fitted PBE0/def2-SVP, and ddCOSMO with dielectric constant 36.7 to approximate DMF. Both +1 cation calculations converged at an energy tolerance of 10⁻⁸ Hartree. The geometries were fixed at the CREST/GFN2-xTB minima.

| Method | Ortho minus meta electronic energy |
| --- | ---: |
| GFN2-xTB / ALPB(DMF) | −5.56 kcal/mol |
| PBE0/def2-SVP / ddCOSMO(DMF), single-point | −3.46 kcal/mol |

Both calculations favor the ortho isomer in the specified continuum environment. These are total electronic-energy differences between molecules with the same composition. They include all substituent-position effects; they do not isolate an intramolecular hydrogen-bond energy.

[DFT energies](../results/quantum_pilot/dft_singlepoints.csv) · [Geometry input files](../results/quantum_pilot/dft_inputs/) · [Method comparison](../results/quantum_pilot/dft_comparison.json)

## Interpretation and next calculation

The machine-learning tie reflects the two-descriptor representation. The conformational analysis identifies an additional distinction: the position of the polar substituent changes whether the ammonium group can contact an oxygen atom within the same cation.

The next atomistic question is whether intramolecular N–H···O contacts alter the geometry or strength of ammonium–iodide interactions. A consistent ion-pair or Pb–I interface model, evaluated across several starting arrangements, would test this. Isolated-cation calculations cannot establish spacer incorporation, defect passivation, or photovoltaic performance.

The DFT calculation is currently a fixed-geometry single-point comparison. Periodic models, DFT geometry relaxation, dispersion/basis convergence and explicit counterions have not been calculated.

## Reproduction

Python 3.13 with RDKit generates starting conformers; CREST and xTB run under Linux. The external working directory is set by `STUDY1_QM_WORKDIR`; `STUDY1_XTB_BIN` identifies the xTB executable when it is absent from PATH.

```bash
# Prepare with RDKit (Windows or Linux environment with RDKit)
python scripts/prepare_quantum_conformers.py --nembed 90 --limit-primary 12 --limit-secondary 7

# GFN2-xTB gas and ALPB(DMF) optimization (Linux)
python scripts/run_quantum_xtb.py --workers 3
```

CREST runs for the matched sulfonyl pair used the lowest ALPB(DMF) xTB geometry:

```bash
crest input.xyz --gfn2 --chrg 1 --alpb dmf --quick -T 4
```

Analyze results using `scripts/analyze_quantum_xtb.py`, `scripts/analyze_quantum_crest.py` and `scripts/plot_quantum_crest.py`.

PBE0 single-point calculations can be rerun from the published [DFT input coordinates](../results/quantum_pilot/dft_inputs/), using PySCF 2.11.0 and the documented ddCOSMO settings:

```bash
python scripts/run_quantum_dft.py --group sulfonyl_meta
python scripts/run_quantum_dft.py --group sulfonyl_ortho
```

The xTB and DFT SCF logs and intermediate files are stored separately from Git. The published coordinates, energies, source IDs and analysis scripts are included in this branch.
