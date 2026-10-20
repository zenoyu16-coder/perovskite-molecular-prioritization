"""Check numerical provenance and scope of the iodide association pilot."""
import json
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/iodide_pilot"

def test_iodide_pair_convergence_and_charge():
    s=pd.read_csv(OUT/"optimized_ion_pairs.csv")
    f=pd.read_csv(OUT/"fluoro_ion_pairs.csv")
    assert len(s)==84 and len(f)==40
    dmf=s[s.medium.eq("dmf")]
    assert len(dmf)==42
    assert set(dmf.group)=={"sulfonyl_meta","sulfonyl_ortho"}
    assert set(f.group)=={"fluoromethoxy_A","fluoromethoxy_B"}
    assert dmf.iodide_xtb_partial_charge.between(-1.01,-.60).all()
    assert f.I_partial_charge.between(-1.01,-.60).all()
    assert dmf.external_Hbond_any.all()
    assert f.external_Hbond_any.all()

def test_energy_differences_and_dft_gate():
    p=json.loads((OUT/"association_differences.json").read_text(encoding="utf-8"))
    dmf=next(q for q in p["energy_differences"] if q["medium"]=="dmf")
    assert 3.0<dmf["ortho_minus_meta_association_kcal_mol"]<4.1
    dft=json.loads((OUT/"dft_association_check.json").read_text(encoding="utf-8"))
    assert dft["pbe0_def2svp_ddcosmo"]["n_ion_pair_scf_converged"]==2
    assert dft["agreement_on_sign"]
    assert 3.0<dft["pbe0_def2svp_ddcosmo"]["delta_association_kcal_mol"]<4.1
    fluoro=json.loads((OUT/"fluoro_transfer.json").read_text(encoding="utf-8"))
    assert fluoro["observed_positive"]
    assert 1.0<fluoro["a_minus_b_electronic_association_kcal_mol"]<2.1

def test_coordinate_files_and_input_counts():
    x=pd.read_csv(OUT/"sulfonyl_initial_poses.csv")
    y=pd.read_csv(OUT/"fluoro_initial_poses.csv")
    assert len(x)==42 and len(y)==40
    for group,n in (("sulfonyl_meta",25),("sulfonyl_ortho",25),
                    ("fluoromethoxy_A",23),("fluoromethoxy_B",23)):
        p=OUT/f"{group}_dmf_lowest.xyz"
        lines=p.read_text(encoding="utf-8").splitlines()
        assert int(lines[0])==n
        assert len(lines)>=n+2
        assert lines[n+1].split()[0]=="I"
