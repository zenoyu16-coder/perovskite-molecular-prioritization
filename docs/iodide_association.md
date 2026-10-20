# Iodide association of rank-tied spacer isomers

The molecular ranker assigns identical Logistic2 scores to two pairs of aromatic spacer-cation regioisomers. GFN2-xTB conformer calculations revealed different ammonium–oxygen contact geometries within each pair. We examined whether these differences are accompanied by changes in the interaction with iodide.

## Ion-pair calculations

We used the CREST/GFN2-xTB conformer ensembles from [the molecular study](quantum_chemistry.md). An iodide ion was placed along each available ammonium N–H vector, initially 2.70 Å from H. Starting poses with excessive nonbonded overlap were excluded using identical distance cutoffs, followed by full ion-pair geometry optimization with GFN2-xTB 6.7.1 and ALPB(DMF).

The calculated energy for each reaction RNH₃⁺ + I⁻ → [RNH₃···I] was referenced to the lowest isolated-cation conformer and isolated iodide, using the same electronic-structure and solvent method. Reported quantities are electronic association energies, in kcal/mol.

| Positional isomer | Optimized starts | GFN2-xTB / ALPB(DMF) association energy |
| --- | ---: | ---: |
| Methylsulfonylbenzylammonium, meta (CID 2760937) | 33 | −20.50 |
| Methylsulfonylbenzylammonium, ortho (CID 16218335) | 9 | −16.92 |
| Fluoro/methoxybenzylammonium A (CID 40786958) | 7 | −17.44 |
| Fluoro/methoxybenzylammonium B (CID 24721203) | 33 | −18.98 |

The meta sulfonyl isomer binds iodide 3.58 kcal/mol more favorably than the ortho isomer in this model. For the fluoro/methoxy pair, isomer B binds iodide 1.54 kcal/mol more favorably than isomer A.

The number of accepted starting poses differs across isomers because the CREST ensembles contain different numbers of conformers and some initial iodide placements fail the overlap criteria. Each lowest-energy geometry was reached from multiple starting configurations or closely related conformers.

## DFT check on the sulfonyl pair

The lowest GFN2-xTB/ALPB(DMF) ion pair of each sulfonyl isomer was evaluated at PBE0/def2-SVP using PySCF 2.11.0, an iodine def2-SVP ECP28 potential and ddCOSMO with a dielectric constant of 36.7. Both SCF calculations converged. The cation reference energies were taken from the PBE0/def2-SVP/ddCOSMO calculations described in the [conformer analysis](quantum_chemistry.md).

For two isomers, the difference in iodide association energies can be written directly as:

`ΔΔE_assoc = [E(pair, ortho) − E(pair, meta)] − [E(cation, ortho) − E(cation, meta)]`

The isolated iodide energy cancels in this comparison.

| Method | ΔΔE_assoc: ortho minus meta |
| --- | ---: |
| GFN2-xTB / ALPB(DMF) | +3.58 kcal/mol |
| PBE0/def2-SVP / ddCOSMO(DMF), single point | +3.64 kcal/mol |

The two methods agree on the ordering. All DFT geometries remain at their GFN2-xTB optimized positions.

## Coordination geometry

In the lowest-energy sulfonyl ion pairs, N···I separations are 3.29 Å (meta) and 3.27 Å (ortho). The ortho cation retains an intramolecular N–H···O contact while also coordinating external iodide through an ammonium hydrogen. The intramolecular contact is absent in the meta isomer.

The fluoro/methoxy pair follows the same association-energy ordering predicted from its free-cation conformers, although the energy difference is smaller. Isomer A retains an N···O separation of 2.78 Å in the optimized ion pair, but its N–H···O angle decreases to about 127°, below the 130° contact threshold used in the conformer analysis.

The two pairs illustrate how substituent position changes molecular conformations and ion-pair interactions in ways not captured by the selected two-descriptor ranker. The calculations do not establish that a single hydrogen-bond competition mechanism explains both pairs.

## Scope

The calculations describe one iodide counterion and an isolated organic cation in a continuum DMF model. Gas-phase GFN2-xTB calculations were also attempted, but their optimized structures showed substantially different iodide partial charges and short H···I contacts; they are excluded from the relative-solution analysis.

Association free energies, explicit solvent coordination, basis-set superposition correction, DFT geometry relaxation and binding in a Pb–I perovskite framework have not been evaluated. The second isomer pair has only GFN2-xTB results. These values can guide a future solution-phase interaction study but do not measure device performance or prove a change in crystal incorporation.

## Files

The reproducible numerical record includes [all optimized sulfonyl ion pairs](../results/iodide_pilot/optimized_ion_pairs.csv), [the corresponding energy comparisons](../results/iodide_pilot/association_differences.json), [DFT differences](../results/iodide_pilot/dft_association_check.json), and [the fluoro/methoxy transfer test](../results/iodide_pilot/fluoro_transfer.json). Starting-pose metadata and selected optimized XYZ structures are in [results/iodide_pilot](../results/iodide_pilot/).

The coordinate-generation, xTB optimization and analysis scripts are `scripts/prepare_iodide_pairs.py`, `scripts/run_iodide_xtb.py`, `scripts/analyze_iodide_pairs.py`, and `scripts/analyze_iodide_transfer.py`. DFT input structures are the lowest-energy `*_dmf_lowest.xyz` files. Use `scripts/run_iodide_dft.py` for the two sulfonyl single points and `scripts/analyze_iodide_dft.py` to compute the double difference.

The raw GFN2-xTB and PySCF logs remain in the local calculation scratch directory.

### Rebuilding the pilot

The following commands use the published CREST ensembles. Set `STUDY1_QM_WORKDIR` to a writable scratch directory and `STUDY1_XTB_BIN` to the xTB 6.7.1 executable, or use the recorded defaults on the originating computer.

```bash
# Prepare the two molecule pairs (RDKit/Python)
python scripts/prepare_iodide_pairs.py --pair sulfonyl
python scripts/prepare_iodide_pairs.py --pair fluoromethoxy

# Optimize ion-pair starts in Linux with GFN2-xTB
python scripts/run_iodide_xtb.py --media dmf gas --workers 4
STUDY1_IONPAIR_SUBDIR=iodide_association_fluoromethoxy python scripts/run_iodide_xtb.py --media dmf --workers 4

# Analyze on a Python environment with RDKit and pandas
python scripts/analyze_iodide_pairs.py
python scripts/analyze_iodide_transfer.py
python scripts/run_iodide_dft.py --group sulfonyl_meta
python scripts/run_iodide_dft.py --group sulfonyl_ortho
python scripts/analyze_iodide_dft.py
```

The iodide single-point reference in GFN2-xTB uses one iodine atom with charge -1 and ALPB(DMF), as specified in the energy-reference fields of the released result table. It is required for absolute GFN2-xTB association energies; it cancels from the DFT double difference.
