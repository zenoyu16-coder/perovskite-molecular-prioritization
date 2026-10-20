"""Compare optimized (+1 cation, -1 iodide) GFN2-xTB ion pairs.

For each medium:
 association electronic energy = E(pair) - E(cation) - E(I-).
Only compare isomers at the same level/solvation model. The difference
in association energy between isomers cancels the isolated iodide
reference. These quantities omit dissociation entropy and explicit solvent.
"""
from __future__ import annotations
import csv,json,os,re,math,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Geometry import Point3D

from prepare_quantum_conformers import ROOT, RESULTS, hbond_features

WORK=Path(os.environ.get("STUDY1_QM_WORKDIR","D:/AI4S/study1_quantum_pilot"))
RAW=WORK/"iodide_association"
OUT=ROOT/"results/iodide_pilot"
HARTREE_TO_KCAL=627.509474
E_PATTERN=re.compile(r"TOTAL ENERGY\s+([-+\d.]+)\s+Eh")

def dot_angle(a,b):
    cos=float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))
    return math.degrees(math.acos(max(-1,min(1,cos))))

def bond_summary(cation_smiles, xyzfile):
    mol=Chem.AddHs(Chem.MolFromSmiles(cation_smiles))
    ncat=mol.GetNumAtoms()
    xyz=xyzfile.read_text(encoding="utf-8").splitlines()
    nat=int(xyz[0].strip())
    if nat!=ncat+1: raise ValueError((str(xyzfile),nat,ncat))
    entries=[line.split() for line in xyz[2:nat+2]]
    if entries[-1][0].upper()!="I":
        raise ValueError("Last atom is not iodide")
    xyz=np.array([[float(v) for v in e[1:4]] for e in entries])
    for i,a in enumerate(mol.GetAtoms()):
        if a.GetSymbol().upper()!=entries[i][0].upper():
            raise ValueError("Cation atom order mismatch")
    conf=Chem.Conformer(ncat)
    for i,coord in enumerate(xyz[:-1]):
        conf.SetAtomPosition(i,Point3D(*coord))
    mol.AddConformer(conf)
    nids=[a.GetIdx() for a in mol.GetAtoms() if
          a.GetSymbol()=="N" and a.GetFormalCharge()>0]
    if len(nids)!=1:raise ValueError("Expected one protonated N")
    n=nids[0]
    ncoord=xyz[n]
    hs=[a.GetIdx() for a in mol.GetAtomWithIdx(n).GetNeighbors()
        if a.GetSymbol()=="H"]
    iodide=xyz[-1]
    nhis=[]
    for hid in hs:
        h=xyz[hid]
        angle=dot_angle(ncoord-h,iodide-h)
        nhis.append((float(np.linalg.norm(iodide-h)),angle,hid))
    contact=sorted(nhis)[0]
    hbond=hbond_features(mol,0)
    return dict(N_I_A=float(np.linalg.norm(ncoord-iodide)),
                min_H_I_A=contact[0],
                N_H_I_angle_deg=contact[1],
                external_Hbond_any=any(d<3.1 and a>=125 for d,a,_ in nhis),
                internal_Hbond=bool(hbond["hb_strict"]),
                internal_N_O_A=hbond["min_NO_A"],
                internal_H_O_A=hbond["min_HO_A"],
                internal_N_H_O_max_deg=hbond["max_NHO_deg"])

def read_atomic_energy(medium):
    """Read the independently computed -1 iodide reference."""
    path = OUT / "iodide_reference.json"
    ref = json.loads(path.read_text(encoding="utf-8"))
    record = ref["results"][medium]
    if not record["converged"] or record["charge"] != -1 or record["gfn"] != 2:
        raise ValueError("Incompatible isolated-iodide reference for " + medium)
    expected = "ALPB(DMF)" if medium == "dmf" else "gas"
    if record["solvation"] != expected:
        raise ValueError("Solvent mismatch in iodide reference")
    return float(record["energy_hartree"])


def cation_energy(group,medium):
    if medium=="dmf":
        df=pd.read_csv(RESULTS/"crest_conformers.csv")
        sub=df[df.group.eq(group)]
        return float(sub.energy_hartree.min())
    df=pd.read_csv(RESULTS/"xtb_conformers.csv")
    return float(df[(df.group.eq(group))&(df.medium.eq(medium))].energy_hartree.min())

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    smiles={("sulfonyl_meta"):"CS(=O)(=O)c1cccc(C[NH3+])c1",
            ("sulfonyl_ortho"):"CS(=O)(=O)c1ccccc1C[NH3+]"}
    ion={m:read_atomic_energy(m) for m in ("dmf","gas")}
    rows=[]
    for rec in csv.DictReader((RAW/"pair_inputs.csv").open(encoding="utf-8")):
        for medium in ("dmf","gas"):
            runpath=RAW/Path(rec["input_file"]).parent/medium/"run.json"
            run=json.loads(runpath.read_text(encoding="utf-8"))
            if not run["success"]:raise RuntimeError("Unsuccessful run "+str(runpath))
            xyz=WORK/run["optimized_xyz"]
            e=float(run["energy_hartree"])
            group=rec["group"]
            mono=cation_energy(group,medium)
            charge_path=xyz.parent/"charges"
            icharge=None
            if charge_path.exists():
                vals=[float(x.strip()) for x in charge_path.read_text().splitlines() if x.strip()]
                icharge=vals[-1]
            d=dict(group=group,conformer=int(rec["conformer"]),
                   site=int(rec["site"]),medium=medium,
                   pair_energy_hartree=e,cation_reference_Eh=mono,
                   iodide_reference_Eh=ion[medium],
                   electronic_association_kcal_mol=(e-mono-ion[medium])*HARTREE_TO_KCAL,
                   iodide_xtb_partial_charge=icharge,
                   raw_optimized_xyz=run["optimized_xyz"],
                   **bond_summary(smiles[group],xyz))
            rows.append(d)
    df=pd.DataFrame(rows)
    df.to_csv(OUT/"optimized_ion_pairs.csv",index=False,float_format="%.9f")
    winners=[]; records={}
    for (group,medium),sub in df.groupby(["group","medium"]):
        best=sub.loc[sub.pair_energy_hartree.idxmin()]
        close=sub[(sub.pair_energy_hartree-best.pair_energy_hartree)*HARTREE_TO_KCAL<=2.]
        valid=sub[sub.external_Hbond_any]
        win=dict(group=group,medium=medium,n_optimized=len(sub),
                 n_external_Hbond=len(valid),
                 lowest_pair_Eh=float(best.pair_energy_hartree),
                 cation_reference_Eh=float(best.cation_reference_Eh),
                 iodide_reference_Eh=float(best.iodide_reference_Eh),
                 association_kcal_mol=float(best.electronic_association_kcal_mol),
                 best_N_I_A=float(best.N_I_A),
                 best_min_H_I_A=float(best.min_H_I_A),
                 best_N_H_I_angle_deg=float(best.N_H_I_angle_deg),
                 best_internal_Hbond=bool(best.internal_Hbond),
                 n_within_2kcal=len(close),
                 n_low_ion_pairs_retain_internalHB=int(close.internal_Hbond.sum()),
                 lowest_iodide_partial_charge=float(best.iodide_xtb_partial_charge),
                 best_source=str(best.raw_optimized_xyz))
        winners.append(win)
        records[(group,medium)]=best
        f=OUT/f"{group}_{medium}_lowest.xyz"
        source=WORK/best.raw_optimized_xyz
        raw=source.read_text(encoding="utf-8").splitlines()
        raw[1]=f"GFN2-xTB 6.7.1; group {group}; {medium}; charge 0; E={best.pair_energy_hartree} Eh"
        f.write_text("\n".join(raw)+"\n",encoding="utf-8")
    summary=pd.DataFrame(winners)
    summary.to_csv(OUT/"ion_pair_summary.csv",index=False,float_format="%.9f")
    deltas=[]
    for medium in ("dmf","gas"):
        a=records[("sulfonyl_ortho",medium)]
        b=records[("sulfonyl_meta",medium)]
        deltas.append(dict(medium=medium,
                    ortho_minus_meta_association_kcal_mol=float(
                        a.electronic_association_kcal_mol - b.electronic_association_kcal_mol),
                    ortho_minus_meta_pair_total_kcal_mol=float(
                        (a.pair_energy_hartree-b.pair_energy_hartree)*HARTREE_TO_KCAL),
                    ortho_minus_meta_monomer_total_kcal_mol=float(
                        (a.cation_reference_Eh-b.cation_reference_Eh)*HARTREE_TO_KCAL)))
    (OUT/"association_differences.json").write_text(json.dumps(dict(
        method="GFN2-xTB/ALPB(DMF) vs gas",
        model="monocation (+1) and isolated iodide (-1) to neutral contact ion pair",
        coefficients="all species at GFN2-xTB 6.7.1; electronic energies only; no thermal corrections",
        charge_check="xTB charge on iodine; inspect ion-pair covalency/charge transfer",
        energy_differences=deltas,
        important_limits=["implicit solvent pair and separated-ion energies lack association entropy",
            "unequal number of accepted initial poses: meta 33, ortho 9",
            "starting geometries along three NH axes only",
            "free cations sampled with CREST in DMF and ETKDG+xTB in gas",
            "only one explicit I- and no periodic PbI4 environment"],
    ),indent=2)+"\n",encoding="utf-8")
    print(summary.to_string(index=False,float_format=lambda x:f"{x:.4f}"))
    print("DELTA association",json.dumps(deltas,indent=2))
if __name__=="__main__":main()
