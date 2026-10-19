"""Unit checks for the ionic representation used by the screening workflow."""

import pytest
from rdkit import Chem

from perovskite_prioritization.chemistry import filter_candidate
from perovskite_prioritization.protonation import protonate_candidate_v2


@pytest.mark.parametrize("precursor,expected_charge", [
    ("NCc1ccccc1", 1),
    ("NCCN", 2),
    ("NCc1ccccc1S(N)(=O)=O", 1),
    ("Cc1nc(-c2ccc(CN)cc2)no1", 1),
    ("CC(=O)N", 0),
    ("c1ccncc1", 0),
])
def test_amine_protonation_charge(precursor, expected_charge):
    result = protonate_candidate_v2(precursor, require_amine=False)
    assert result.formal_charge_after == expected_charge
    assert Chem.GetFormalCharge(Chem.MolFromSmiles(result.modeled_smiles)) == expected_charge


def test_monoion_candidate_filter():
    primary = filter_candidate("NCc1ccccc1")
    tertiary = filter_candidate("CN(C)Cc1ccccc1")
    assert primary.accepted and tertiary.accepted
    assert Chem.GetFormalCharge(Chem.MolFromSmiles(tertiary.modeled_cation_smiles)) == 1


def test_aromatic_nitrogen_not_protonated_by_default():
    result = filter_candidate("Cc1nc(-c2ccc(CN)cc2)no1")
    assert result.accepted
    assert Chem.GetFormalCharge(Chem.MolFromSmiles(result.modeled_cation_smiles)) == 1


def test_stereochemistry_is_retained():
    result = protonate_candidate_v2("N[C@@H](C)CO")
    assert "@" in result.modeled_smiles


def test_unsupported_policy_rejected():
    with pytest.raises(ValueError, match="Unknown protonation policy"):
        filter_candidate("NCc1ccccc1", protonation_policy="unsupported")
