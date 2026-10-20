"""Tabulate the converged PySCF PBE0 single-point results.

Different positional isomers have the same formula, but their total-energy
difference is an isomerization electronic-energy difference. It is not an
adsorption energy, crystal formation energy, or photovoltaic prediction.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/quantum_pilot"
RAW=Path(os.environ.get("STUDY1_QM_WORKDIR",str(ROOT/"data/external/study1_qm")))
TO_KCAL=627.509474


def run():
    cases=[]
    for name in ("sulfonyl_meta","sulfonyl_ortho"):
        p=RAW/"dft"/name/"pbe0_def2svp_dmf.json"
        result=json.loads(p.read_text(encoding="utf-8"))
        if not result["converged"]:raise RuntimeError(f"DFT not converged: {name}")
        if result["xc"]!="pbe0" or result["basis"]!="def2-svp":
            raise ValueError("inconsistent method")
        cases.append(result)
    frame=pd.DataFrame(cases)
    frame["relative_energy_kcal_mol"]=(frame.energy_hartree-frame.energy_hartree.min())*TO_KCAL
    frame.to_csv(OUT/"dft_singlepoints.csv",index=False,float_format="%.10f")
    c=pd.read_csv(OUT/"crest_summary.csv")
    cmin=c.set_index("group")["E_min_hartree"]
    dmin=frame.set_index("group")["energy_hartree"]
    record=dict(
        model="two isolated isomeric cations in implicit DMF",
        dft_method="PySCF 2.11.0 PBE0/def2-SVP, density fitting, ddCOSMO eps=36.7",
        geometry="lowest CREST GFN2-xTB/ALPB(DMF) conformer; fixed",
        dft_ortho_minus_meta_kcal_mol=float(
            (dmin["sulfonyl_ortho"]-dmin["sulfonyl_meta"])*TO_KCAL),
        xtb_ortho_minus_meta_kcal_mol=float(
            (cmin["sulfonyl_ortho"]-cmin["sulfonyl_meta"])*TO_KCAL),
        n_converged_scf=int(frame.converged.sum()),
        interpretation="Relative electronic energies of positional isomers at two electronic-structure levels",
        not_assessed=["periodic incorporation","iodide binding","formation free energy",
                      "device power conversion efficiency","DFT optimized geometries",
                      "dispersion correction and basis set convergence"])
    (OUT/"dft_comparison.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    print(frame[["group","energy_hartree","relative_energy_kcal_mol",
                 "wall_seconds"]].to_string(index=False))
    print(json.dumps(record,indent=2))


if __name__=="__main__":run()
