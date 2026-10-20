"""Read the CREST GFN2-xTB/ALPB(DMF) conformer ensembles.

Input structures are generated from the frozen molecular shortlist. CREST
sampling is quick mode and cannot establish an exhaustive conformer landscape.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Geometry import Point3D

from prepare_quantum_conformers import ROOT, PRIVATE, RESULTS, CANDIDATES, hbond_features

KCAL_PER_HARTREE = 627.509474


def load_ensemble(path: Path, smiles: str):
    lines = path.read_text(encoding="utf-8").splitlines()
    mol0 = Chem.AddHs(Chem.MolFromSmiles(smiles))
    n = mol0.GetNumAtoms()
    frames = []
    i=0
    while i<len(lines):
        nn=int(lines[i].strip())
        if nn != n:
            raise ValueError(f"Atom count mismatch in {path}")
        energy=float(lines[i+1].strip().split()[0])
        mol=Chem.Mol(mol0)
        conf=Chem.Conformer(n)
        for a,line in enumerate(lines[i+2:i+2+n]):
            fields=line.split()
            symbol=fields[0]
            if symbol.upper()!=mol.GetAtomWithIdx(a).GetSymbol().upper():
                raise ValueError(f"Atom order mismatch in {path} atom {a}")
            xyz=tuple(map(float,fields[1:4]))
            conf.SetAtomPosition(a,Point3D(*xyz))
        mol.AddConformer(conf)
        frames.append((energy,mol))
        i+=n+2
        while i<len(lines) and not lines[i].strip():
            i+=1
    return frames


def main():
    inputs = {CANDIDATES[int(row["pubchem_cid"])][0]:row for row in
              csv.DictReader((ROOT/"results/clean_study1/shortlist.csv").open(encoding="utf-8"))
              if int(row["pubchem_cid"]) in CANDIDATES}
    rows=[]
    summary=[]
    for group in ("sulfonyl_meta","sulfonyl_ortho"):
        p=RESULTS/"crest_ensembles"/f"{group}.xyz"
        log=PRIVATE/"crest"/group/"crest.stdout"
        if not p.is_file():
            raise FileNotFoundError(p)
        if log.exists() and "CREST terminated normally." not in log.read_text(encoding="utf-8",errors="replace"):
            raise RuntimeError("Incomplete CREST run "+group)
        frames=load_ensemble(p,inputs[group]["modeled_cation_smiles"])
        low=min(e for e,_ in frames)
        for rank,(energy,mol) in enumerate(frames,1):
            rows.append(dict(group=group,conformer=rank,method="CREST-3.0.2-GFN2-xTB",
                environment="ALPB-dmf",energy_hartree=energy,
                relative_kcal_mol=(energy-low)*KCAL_PER_HARTREE,
                **hbond_features(mol,0)))
        sub=pd.DataFrame([row for row in rows if row["group"]==group])
        near=sub[sub.relative_kcal_mol<=3.]
        best=sub.loc[sub.energy_hartree.idxmin()]
        summary.append(dict(group=group,crest_conformers=len(frames),
                            E_min_hartree=best.energy_hartree,
                            low_3kcal_conformers=len(near),
                            low_3kcal_hb_count=int(near.hb_strict.sum()),
                            best_N_O_A=best.min_NO_A,
                            best_H_O_A=best.min_HO_A,
                            best_N_H_O_deg=best.max_NHO_deg,
                            best_hbond=bool(best.hb_strict)))
    RESULTS.mkdir(exist_ok=True,parents=True)
    pd.DataFrame(rows).to_csv(RESULTS/"crest_conformers.csv",index=False,float_format="%.8f")
    pd.DataFrame(summary).to_csv(RESULTS/"crest_summary.csv",index=False,float_format="%.8f")
    print(pd.DataFrame(summary).to_string(index=False,float_format=lambda x:f"{x:.4f}"))
    print("Near-minimum contact counts:",{x["group"]:(x["low_3kcal_hb_count"],
                                                    x["low_3kcal_conformers"]) for x in summary})


if __name__=="__main__":
    main()
