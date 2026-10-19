"""Deterministic chemistry utilities used by the analysis and screening."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    Crippen,
    Descriptors,
    Lipinski,
    rdFingerprintGenerator,
    rdMolDescriptors,
)
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from .protonation import PROTONATION_POLICY_V2, protonate_candidate_v2

DESCRIPTOR_COLUMNS = (
    "MolWt", "HeavyAtomMolWt", "HeavyAtomCount", "TPSA", "HBD", "HBA",
    "FormalCharge", "NitrogenCount", "PositiveNitrogenCount", "OxygenCount",
    "SulfurCount", "FluorineCount", "ChlorineCount", "BromineCount",
    "IodineCount", "HeteroatomCount", "RotatableBondCount", "FractionCSP3",
    "RingCount", "AromaticRingCount", "AliphaticRingCount", "BertzCT", "Chi0v",
    "Chi1v", "Chi2v", "Chi3v", "Kappa1", "Kappa2", "Kappa3", "MolLogP",
    "MolMR", "CarbonCount", "CarbonHeavyAtomFraction",
    "NitrogenHeavyAtomFraction", "HeteroatomHeavyAtomFraction",
)

ALLOWED_ELEMENTS = {"C", "H", "N", "O", "F", "Cl", "Br", "S"}
PROTONATABLE = Chem.MolFromSmarts(
    "[$([N;H1,H2,H3;+0;!$(N-C=O);!$(N-S(=O)=O);!$(N=*)]),$([nH0;+0])]"
)
EXTREME_REACTIVE = {
    "acid_halide": Chem.MolFromSmarts("[CX3](=[OX1])[F,Cl,Br,I]"),
    "isocyanate": Chem.MolFromSmarts("N=C=O"),
    "isothiocyanate": Chem.MolFromSmarts("N=C=S"),
    "peroxide": Chem.MolFromSmarts("[OX2]-[OX2]"),
    "azide": Chem.MolFromSmarts("[N-]=[N+]=N"),
    "diazonium": Chem.MolFromSmarts("[N+]#N"),
}


def mol_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def protonate_basic_sites(source_smiles: str) -> str:
    """Apply the basic-site convention used in the frozen molecular registry."""
    molecule = mol_from_smiles(source_smiles)
    editable = Chem.RWMol(molecule)
    for atom in editable.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetFormalCharge() != 0:
            continue
        if atom.GetIsAromatic():
            if atom.GetTotalNumHs() == 0 and atom.GetDegree() == 2:
                atom.SetFormalCharge(1)
                atom.SetNumExplicitHs(1)
                atom.SetNoImplicit(True)
            continue
        if atom.GetHybridization() == Chem.HybridizationType.SP2:
            neighbors = list(atom.GetNeighbors())
            if any(
                n.GetSymbol() in {"C", "S"}
                and any(
                    editable.GetBondBetweenAtoms(n.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() >= 2
                    and o.GetSymbol() in {"O", "N", "S"}
                    for o in n.GetNeighbors() if o.GetIdx() != atom.GetIdx()
                )
                for n in neighbors
            ):
                continue
        atom.SetFormalCharge(1)
    result = editable.GetMol()
    Chem.SanitizeMol(result)
    return Chem.MolToSmiles(result, canonical=True, isomericSmiles=False)


def descriptors(smiles: str) -> dict[str, float]:
    mol = mol_from_smiles(smiles)
    atoms = list(mol.GetAtoms())
    heavy = mol.GetNumHeavyAtoms()
    count = lambda symbol: sum(a.GetSymbol() == symbol for a in atoms)
    carbon = count("C")
    nitrogen = count("N")
    hetero = int(Lipinski.NumHeteroatoms(mol))
    values = {
        "MolWt": Descriptors.MolWt(mol),
        "HeavyAtomMolWt": Descriptors.HeavyAtomMolWt(mol),
        "HeavyAtomCount": heavy,
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "FormalCharge": Chem.GetFormalCharge(mol),
        "NitrogenCount": nitrogen,
        "PositiveNitrogenCount": sum(a.GetSymbol() == "N" and a.GetFormalCharge() > 0 for a in atoms),
        "OxygenCount": count("O"), "SulfurCount": count("S"),
        "FluorineCount": count("F"), "ChlorineCount": count("Cl"),
        "BromineCount": count("Br"), "IodineCount": count("I"),
        "HeteroatomCount": hetero,
        "RotatableBondCount": Lipinski.NumRotatableBonds(mol),
        "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
        "RingCount": Lipinski.RingCount(mol),
        "AromaticRingCount": Lipinski.NumAromaticRings(mol),
        "AliphaticRingCount": Lipinski.NumAliphaticRings(mol),
        "BertzCT": Descriptors.BertzCT(mol),
        "Chi0v": Descriptors.Chi0v(mol), "Chi1v": Descriptors.Chi1v(mol),
        "Chi2v": Descriptors.Chi2v(mol), "Chi3v": Descriptors.Chi3v(mol),
        "Kappa1": Descriptors.Kappa1(mol), "Kappa2": Descriptors.Kappa2(mol),
        "Kappa3": Descriptors.Kappa3(mol), "MolLogP": Crippen.MolLogP(mol),
        "MolMR": Crippen.MolMR(mol), "CarbonCount": carbon,
        "CarbonHeavyAtomFraction": carbon / heavy,
        "NitrogenHeavyAtomFraction": nitrogen / heavy,
        "HeteroatomHeavyAtomFraction": hetero / heavy,
    }
    if not np.isfinite([values[c] for c in DESCRIPTOR_COLUMNS]).all():
        raise ValueError("Non-finite molecular descriptors")
    return {c: float(values[c]) for c in DESCRIPTOR_COLUMNS}


def fingerprint(smiles: str, n_bits: int = 256) -> np.ndarray:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    fp = generator.GetFingerprint(mol_from_smiles(smiles))
    arr = np.zeros(n_bits, dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def tanimoto_matrix(query: np.ndarray, reference: np.ndarray) -> np.ndarray:
    q = query.astype(bool)
    r = reference.astype(bool)
    intersection = q.astype(np.uint16) @ r.astype(np.uint16).T
    union = q.sum(axis=1)[:, None] + r.sum(axis=1)[None, :] - intersection
    return np.divide(intersection, union, out=np.zeros_like(intersection, dtype=float), where=union > 0)


def synthetic_accessibility(mol: Chem.Mol) -> float:
    """Return Ertl SA score (1 easy, 10 difficult) from RDKit's contributed scorer."""
    try:
        from rdkit.Contrib.SA_Score import sascorer
        return float(sascorer.calculateScore(mol))
    except (ImportError, AttributeError):
        # Auditable descriptor fallback if a minimal RDKit build omits Contrib.
        rings = Lipinski.RingCount(mol)
        stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        score = 1.0 + 0.012 * Descriptors.MolWt(mol) + 0.18 * Lipinski.NumRotatableBonds(mol)
        score += 0.22 * rings + 0.35 * stereo + 0.002 * Descriptors.BertzCT(mol)
        return float(np.clip(score, 1.0, 10.0))


@dataclass(frozen=True)
class FilterResult:
    accepted: bool
    reason: str
    canonical_smiles: str | None = None
    modeled_cation_smiles: str | None = None
    molecular_formula: str | None = None
    molecular_weight: float | None = None
    sa_score: float | None = None


def _catalog() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


def filter_candidate(
    smiles: str, catalog: FilterCatalog | None = None, *, protonation_policy: str = PROTONATION_POLICY_V2
) -> FilterResult:
    """Screen a precursor with the amine-only charge-state convention.

    Use ``protonation_policy="legacy"`` only when reproducing the frozen
    model and its original feature/score tables.
    """
    if protonation_policy not in {"legacy", PROTONATION_POLICY_V2}:
        raise ValueError(f"Unknown protonation policy: {protonation_policy}")
    try:
        mol = mol_from_smiles(smiles)
    except ValueError:
        return FilterResult(False, "invalid_smiles")
    if len(Chem.GetMolFrags(mol)) != 1:
        return FilterResult(False, "salt_or_mixture")
    elements = {a.GetSymbol() for a in mol.GetAtoms()}
    if not elements <= ALLOWED_ELEMENTS:
        return FilterResult(False, "disallowed_element")
    mw = float(Descriptors.MolWt(mol))
    if not 70 <= mw <= 350:
        return FilterResult(False, "mw_outside_70_350")
    if protonation_policy == "legacy":
        if PROTONATABLE is None or not mol.HasSubstructMatch(PROTONATABLE):
            return FilterResult(False, "no_plausible_protonatable_amine")
    else:
        # The historic SMARTS excludes tertiary aliphatic amines and admits
        # aromatic nitrogens. Use the same amine criterion as the current
        # representation, not the original model's prefilter.
        try:
            protonate_candidate_v2(Chem.MolToSmiles(mol, isomericSmiles=True))
        except ValueError:
            return FilterResult(False, "no_plausible_amine_v2")
    for name, pattern in EXTREME_REACTIVE.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return FilterResult(False, f"reactive_{name}")
    active_catalog = catalog or _catalog()
    match = active_catalog.GetFirstMatch(mol)
    if match is not None:
        return FilterResult(False, f"rdkit_alert_{match.GetDescription()}")
    canonical = Chem.MolToSmiles(
        mol, canonical=True, isomericSmiles=(protonation_policy == PROTONATION_POLICY_V2)
    )
    try:
        if protonation_policy == PROTONATION_POLICY_V2:
            cation = protonate_candidate_v2(canonical).modeled_smiles
        else:
            cation = protonate_basic_sites(canonical)
    except (ValueError, RuntimeError):
        return FilterResult(False, "no_plausible_amine_v2" if protonation_policy == PROTONATION_POLICY_V2 else "protonation_failed")
    return FilterResult(
        True, "accepted", canonical, cation,
        rdMolDescriptors.CalcMolFormula(mol), mw, synthetic_accessibility(mol),
    )
