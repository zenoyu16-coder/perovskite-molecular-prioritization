"""Evaluate the one-pose-per-isomer DFT countercheck of xTB ion-pair ranking.

Double difference: (E_pair_ortho-E_pair_meta) -
                   (E_cation_ortho-E_cation_meta).
The iodine reference cancels identically. DFT uses the same functional,
basis, ECP and dielectric for the two ion pairs, and the same settings
except for I ECP for the isolated cations.

The DFT calculations are fixed-geometry single points on GFN2-xTB
structures and are not an independent conformational/geometry refinement.
"""
from __future__ import annotations
import json,os
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/iodide_pilot"
WORK=Path(os.environ.get("STUDY1_QM_WORKDIR","D:/AI4S/study1_quantum_pilot"))
HARTREE_TO_KCAL=627.509474

def main():
    isomers=("sulfonyl_meta","sulfonyl_ortho")
    input_energies={}
    for isomer in isomers:
        pc=WORK/"dft"/isomer/"pbe0_def2svp_dmf.json"
        pi=WORK/"iodide_association"/"dft"/isomer/"pbe0_pair.json"
        c=json.loads(pc.read_text(encoding="utf-8"))
        d=json.loads(pi.read_text(encoding="utf-8"))
        if not c["converged"] or not d["converged"]:raise RuntimeError("SCF unconverged")
        if c["basis"].lower()!=d["basis"].lower() or c["xc"].lower()!=d["xc"].lower():
            raise RuntimeError("Incompatible DFT settings")
        input_energies[isomer]={"cation_eh":c["energy_hartree"],
                                "pair_eh":d["energy_hartree"],
                                "ion_pair_scf_converged":d["converged"]}
    meta,ortho=input_energies["sulfonyl_meta"],input_energies["sulfonyl_ortho"]
    de_pair=(ortho["pair_eh"]-meta["pair_eh"])*HARTREE_TO_KCAL
    de_cat=(ortho["cation_eh"]-meta["cation_eh"])*HARTREE_TO_KCAL
    dd=de_pair-de_cat
    xtb=json.loads((OUT/"association_differences.json").read_text(encoding="utf-8"))
    x=[r for r in xtb["energy_differences"] if r["medium"]=="dmf"][0]
    print("PBE0: ortho - meta cation energy:",de_cat,"kcal/mol")
    print("PBE0: ortho - meta ion-pair total energy:",de_pair,"kcal/mol")
    print("PBE0: ortho - meta association energy:",dd,"kcal/mol")
    print("xTB: ortho - meta association energy:",
          x["ortho_minus_meta_association_kcal_mol"],"kcal/mol")
    report=dict(
        model="(+1 organic cation and -1 iodide) ion pair in implicit DMF",
        relative_association_energy_definition=
            "(ortho pair-meta pair)-(ortho cation-meta cation)",
        pbe0_def2svp_ddcosmo=dict(delta_pair_kcal_mol=de_pair,
                                   delta_cation_kcal_mol=de_cat,
                                   delta_association_kcal_mol=dd,
                                   n_ion_pair_scf_converged=2,
                                   iodine_ECP="def2-svp ECP28"),
        gfn2_xtb_alpb_dmf=dict(delta_association_kcal_mol=
                    x["ortho_minus_meta_association_kcal_mol"]),
        agreement_on_sign=(dd*x["ortho_minus_meta_association_kcal_mol"]>0),
        interpretation="Positive sign: stronger electronic iodide association for meta; negative: ortho",
        fixed_geometry_limitation="DFT only single points at GFN2-xTB optimized geometries",
        uncertainty="No free energies, solvent molecules, method convergence, BSSE correction or periodic perovskite"
    )
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"dft_association_check.json").write_text(json.dumps(report,indent=2)+"\n",
                                                  encoding="utf-8")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
