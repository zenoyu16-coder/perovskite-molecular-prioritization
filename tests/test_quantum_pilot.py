"""Checks for frozen-candidate provenance and the published quantum pilot."""
from pathlib import Path
import json

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/quantum_pilot"


def test_ranker_cannot_distinguish_regioisomers():
    table=pd.read_csv(ROOT/"results/clean_study1/shortlist.csv").set_index("pubchem_cid")
    for a,b in ((2760937,16218335),(40786958,24721203)):
        assert table.loc[a,"charge"] == table.loc[b,"charge"] == 1
        assert abs(table.loc[a,"logistic2_score"]-table.loc[b,"logistic2_score"]) < 1e-8
        ma=Chem.MolFromSmiles(table.loc[a,"modeled_cation_smiles"])
        mb=Chem.MolFromSmiles(table.loc[b,"modeled_cation_smiles"])
        assert rdMolDescriptors.CalcMolFormula(ma) == rdMolDescriptors.CalcMolFormula(mb)


def test_xtb_and_crest_results():
    xtb=pd.read_csv(OUT/"xtb_conformers.csv")
    assert len(xtb)==26
    assert set(xtb.medium)=={"gas","dmf"}
    assert xtb.success.all()
    crest=pd.read_csv(OUT/"crest_conformers.csv")
    assert len(crest)==16
    for group,n,contacts in (("sulfonyl_meta",11,0),("sulfonyl_ortho",4,4)):
        low=crest[(crest.group==group)&(crest.relative_kcal_mol <= 3)]
        assert len(low)==n
        assert int(low.hb_strict.sum())==contacts
    assert (OUT/"crest_isomer_conformations.png").stat().st_size>10000


def test_dft_source_and_energy():
    from pytest import approx

    report=json.loads((OUT/"dft_comparison.json").read_text(encoding="utf-8"))
    dft=pd.read_csv(OUT/"dft_singlepoints.csv").set_index("group")
    assert len(dft)==2
    assert dft["converged"].all()
    assert report["n_converged_scf"] == int(dft["converged"].sum())

    group_a, group_b = "sulfonyl_meta", "sulfonyl_ortho"
    delta_dft=(dft.loc[group_b,"energy_hartree"] -
               dft.loc[group_a,"energy_hartree"])*627.509474
    assert report["dft_ortho_minus_meta_kcal_mol"] == approx(delta_dft, abs=2e-6)

    crest=pd.read_csv(OUT/"crest_summary.csv").set_index("group")
    delta_xtb=(crest.loc[group_b,"E_min_hartree"] -
               crest.loc[group_a,"E_min_hartree"])*627.509474
    assert report["xtb_ortho_minus_meta_kcal_mol"] == approx(delta_xtb, abs=2e-5)

    for group in (group_a,group_b):
        geom=OUT/"dft_inputs"/f"{group}.xyz"
        lines=geom.read_text(encoding="utf-8").splitlines()
        assert int(lines[0])==int(dft.loc[group,"n_atoms"])
        assert len(lines)>=int(dft.loc[group,"n_atoms"])+2
