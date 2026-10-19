"""Versioned, conservative amine protonation for virtual spacer screening.

This is a molecular representation convention, not a pH/speciation prediction.
The old study's frozen molecular registry uses a separate legacy convention.
"""
from __future__ import annotations
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

PROTONATION_POLICY_V2 = "amine_only_v2"


@dataclass(frozen=True)
class ProtonationV2:
    source_smiles: str
    modeled_smiles: str
    protonated_atom_indices: tuple[int, ...]
    formal_charge_before: int
    formal_charge_after: int
    warnings: tuple[str, ...] = ()


def _is_resonance_deactivated_nitrogen(atom: Chem.Atom) -> bool:
    """Skip N coupled to C=O, C=S, C=N, or S(=O) groups."""
    for neighbor in atom.GetNeighbors():
        element = neighbor.GetSymbol()
        if element not in {"C", "S", "P"}:
            continue
        for bond in neighbor.GetBonds():
            if bond.GetOtherAtom(neighbor).GetIdx() == atom.GetIdx():
                continue
            if bond.GetBondTypeAsDouble() >= 2.0 and bond.GetOtherAtom(neighbor).GetSymbol() in {"O", "N", "S"}:
                return True
    return False


def protonate_candidate_v2(smiles: str, *, require_amine: bool = True) -> ProtonationV2:
    """Protonate neutral nonaromatic primary/secondary/tertiary amine sites only.

    Excludes aromatic [n], sulfonamides, amides, imine-like N and existing cations.
    May protonate >1 independent amine (putative diammonium). Real protonation
    fractions, counterions, and solid-state valence require separate validation.
    """
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    before = int(Chem.GetFormalCharge(mol))
    editable = Chem.RWMol(mol)
    protonated = []
    for atom in editable.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 0:
            continue
        if atom.GetIsAromatic() or atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        if atom.GetTotalValence() > 3 or _is_resonance_deactivated_nitrogen(atom):
            continue
        atom.SetFormalCharge(+1)
        protonated.append(atom.GetIdx())
    if require_amine and not protonated:
        raise ValueError("No unprotonated, non-resonance-deactivated amine identified")
    result = editable.GetMol()
    Chem.SanitizeMol(result)
    warnings = []
    if len(protonated) > 1:
        warnings.append("multiple_amine_sites: assess mono-/di-cation speciation and salt stoichiometry")
    if Chem.GetFormalCharge(result) > 2:
        warnings.append("charge_gt_2: requires explicit chemical review")
    return ProtonationV2(
        source_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        modeled_smiles=Chem.MolToSmiles(result, canonical=True, isomericSmiles=True),
        protonated_atom_indices=tuple(protonated),
        formal_charge_before=before,
        formal_charge_after=int(Chem.GetFormalCharge(result)),
        warnings=tuple(warnings),
    )
