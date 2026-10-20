"""Prepare matched (+1 organic cation, -1 iodide) ion-pair starting geometries.

Uses the complete published CREST/ALPB(DMF) conformer ensembles. For every
cation conformer, place iodide along each of the three N-H vectors at
an initial H..I distance of 2.7 Angstrom. Retain geometries meeting the
same, predeclared distance criteria across both regioisomers.

The site labels are *initial-placement identifiers*; resulting binding
poses are allowed to relax completely in GFN2-xTB.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import numpy as np
from rdkit import Chem
from analyze_quantum_crest import load_ensemble
from prepare_quantum_conformers import ROOT, RESULTS, CANDIDATES

SCRATCH=Path(os.environ.get("STUDY1_QM_WORKDIR","D:/AI4S/study1_quantum_pilot"))
OUT=SCRATCH/"iodide_association"
PAIR_LIST=OUT/"pair_inputs.csv"
H_I_INITIAL=2.70
MIN_ION_OTHER_HEAVY=2.90
MIN_ION_OTHER_H=1.95

def xyz(mol, iodide_coord, label):
    c=mol.GetConformer()
    atoms=[(a.GetSymbol(),np.array(c.GetAtomPosition(a.GetIdx()))) for a in mol.GetAtoms()]
    atoms.append(("I", np.array(iodide_coord)))
    data=[str(len(atoms)),label]
    for sym,pos in atoms:
        data.append(f"{sym:<2s} {pos[0]: .9f} {pos[1]: .9f} {pos[2]: .9f}")
    return "\n".join(data)+"\n"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pair", choices=("sulfonyl", "fluoromethoxy"), default="sulfonyl")
    p.add_argument("--max-conformers",type=int,default=0,
                   help="0: all CREST conformers, or max per molecule")
    p.add_argument("--h-i",type=float,default=H_I_INITIAL)
    args=p.parse_args()
    if args.h_i < 2.4 or args.h_i > 3.0:
        raise ValueError("Initial H..I must be between 2.4 and 3.0 Å")
    global OUT, PAIR_LIST
    pair_ids = (2760937,16218335) if args.pair == "sulfonyl" else (40786958,24721203)
    OUT = SCRATCH / ("iodide_association" if args.pair == "sulfonyl" else "iodide_association_fluoromethoxy")
    PAIR_LIST = OUT / "pair_inputs.csv"
    OUT.mkdir(parents=True, exist_ok=True)
    source={int(row["pubchem_cid"]):row for row in csv.DictReader(
        (ROOT/"results/clean_study1/shortlist.csv").open(encoding="utf-8"))}
    rows=[]; attempted=0; rejects=[]
    groups=[CANDIDATES[cid][0] for cid in pair_ids]
    for cid in pair_ids:
        group=CANDIDATES[cid][0]
        smi=source[cid]["modeled_cation_smiles"]
        ens=load_ensemble(RESULTS/"crest_ensembles"/f"{group}.xyz",smi)
        if args.max_conformers: ens=ens[:args.max_conformers]
        for confidx,(energy,mol) in enumerate(ens,1):
            nids=[a.GetIdx() for a in mol.GetAtoms() if
                  a.GetSymbol()=="N" and a.GetFormalCharge()>0
                  and sum(x.GetSymbol()=="H" for x in a.GetNeighbors())>=2]
            if len(nids)!=1:raise RuntimeError((group,confidx,nids))
            nidx=nids[0]
            hids=sorted(a.GetIdx() for a in mol.GetAtomWithIdx(nidx).GetNeighbors()
                        if a.GetSymbol()=="H")
            conf=mol.GetConformer()
            n=np.array(conf.GetAtomPosition(nidx))
            for pose,hidx in enumerate(hids,1):
                attempted+=1
                h=np.array(conf.GetAtomPosition(hidx))
                unit=(h-n)/np.linalg.norm(h-n)
                i=h+unit*args.h_i
                dother=[]
                for atom in mol.GetAtoms():
                    idx=atom.GetIdx()
                    if idx in (nidx,hidx):continue
                    pt=np.array(conf.GetAtomPosition(idx))
                    distance=float(np.linalg.norm(i-pt))
                    limit=MIN_ION_OTHER_H if atom.GetSymbol()=="H" else MIN_ION_OTHER_HEAVY
                    if distance<limit:dother.append((idx,atom.GetSymbol(),round(distance,3)))
                if dother:
                    rejects.append(dict(group=group,conf=confidx,pose=pose,
                                        why="clash",bad=dother))
                    continue
                work=OUT/group/f"conf_{confidx:03d}"/f"site_{pose}"
                work.mkdir(parents=True,exist_ok=True)
                rel=work.relative_to(OUT)
                file=work/"input.xyz"
                file.write_text(xyz(mol,i,
                    f"CID {cid}; conformer {confidx}; NH_site {pose}; charge 0; GFN2-xTB/ALPB(DMF)"),encoding="utf-8")
                rows.append(dict(group=group,pubchem_cid=cid,conformer=confidx,
                                 site=pose,organic_charge=1,iodide_charge=-1,
                                 total_charge=0,
                                 crest_energy_hartree=energy,
                                 initial_HI_A=args.h_i,
                                 initial_NI_A=float(np.linalg.norm(i-n)),
                                 input_file=str(rel/"input.xyz").replace("\\","/")))
    if not rows:raise RuntimeError("No valid initial ion-pair geometries")
    with PAIR_LIST.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (OUT/"setup.json").write_text(json.dumps(dict(
        method="CREST GFN2-xTB/ALPB(DMF) conformers + NH-axis iodide placement",
        group_counts={k:sum(r["group"]==k for r in rows) for k in groups},
        n_initial_attempts=attempted,retained=len(rows),rejected=len(rejects),
        rejected_geometries=rejects,initial_HI_A=args.h_i,
        initial_clash_filters=dict(other_heavy_A=MIN_ION_OTHER_HEAVY,other_H_A=MIN_ION_OTHER_H),
        note="Not a physical ion-pair energy result; initial geometry enumeration"
    ),indent=2)+"\n",encoding="utf-8")
    print("attempted",attempted,"retained",len(rows),"rejected",len(rejects))
    print("counts", {g:sum(r["group"]==g for r in rows) for g in groups})

if __name__=="__main__":
    main()
