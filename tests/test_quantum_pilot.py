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
    s=json.loads((OUT/"dft_comparison.json").read_text(encoding="utf-8"))
    assert s["n_converged_scf"]==2
    assert s["dft_ortho_minus_meta_kcal_mol"]<0
    assert s["xtb_ortho_minus_meta_kcal_mol"]<0
    dft=pd.read_csv(OUT/"dft_singlepoints.csv")
    assert len(dft)==2 and dft.converged.all()
    for group in ("sulfonyl_meta","sulfonyl_ortho"):
        geom=OUT/"dft_inputs"/f"{group}.xyz"
        lines=geom.read_text(encoding="utf-8").splitlines()
        assert int(lines[0])==24 and len(lines)>=26
