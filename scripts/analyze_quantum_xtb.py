"""Analyze actual GFN2-xTB optimized conformers without inflating repeated minima.

All energies are electronic energies of isolated +1 cations, in gas or
ALPB(DMF). Relative energies are taken within the same isomer and solvent.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Geometry import Point3D

from prepare_quantum_conformers import ROOT, PRIVATE, RESULTS, hbond_features

HARTREE_TO_KCAL = 627.509474
HARTREE_TO_EV = 27.2113862459


def rdkit_from_xyz(smiles: str, file: Path):
    lines = file.read_text(encoding="utf-8").splitlines()
    nat = int(lines[0].strip())
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if mol.GetNumAtoms() != nat:
        raise ValueError(f"Atom count mismatch: {file}")
    if len(lines) < nat+2:
        raise ValueError(f"Truncated XYZ: {file}")
    conf = Chem.Conformer(nat)
    for i, line in enumerate(lines[2:nat+2]):
        sym, *xyz = line.split()
        if mol.GetAtomWithIdx(i).GetSymbol().upper() != sym.upper():
            raise ValueError(f"Atom index/element mismatch: {file}:{i}")
        x, y, z = map(float, xyz[:3])
        conf.SetAtomPosition(i, Point3D(x, y, z))
    mol.AddConformer(conf)
    return mol


def read_data():
    sources = {(row["group"],int(row["conformer"])):row
               for row in csv.DictReader(
                   (PRIVATE/"conformers.csv").open(encoding="utf-8"))}
    rows = []
    mols = {}
    for (group, conf), origin in sources.items():
        for medium in ("gas", "dmf"):
            p = PRIVATE/"qm"/group/f"conf_{conf:03d}"/medium/"run.json"
            run = json.loads(p.read_text(encoding="utf-8"))
            if not run["success"]:
                rows.append(dict(group=group, conformer=conf, medium=medium,
                                 success=False))
                continue
            geom = rdkit_from_xyz(origin["smiles"], PRIVATE/run["final_xyz"])
            feature = hbond_features(geom, 0)
            key=(group,conf,medium)
            mols[key]=geom
            rows.append(dict(group=group, pubchem_cid=int(origin["pubchem_cid"]),
                             conformer=conf, medium=medium, success=True,
                             energy_hartree=run["energy_hartree"],
                             energy_eV=run["energy_hartree"]*HARTREE_TO_EV,
                             **feature))
    df = pd.DataFrame(rows)
    if not df.success.all():
        raise RuntimeError(f"Some xTB runs failed: {df[~df.success].to_dict('records')}")
    df["relative_kcal_mol"] = (
        (df.energy_hartree-df.groupby(["group","medium"]).energy_hartree.transform("min"))
         * HARTREE_TO_KCAL)
    return df,mols


def save_representatives(df, mols):
    out = RESULTS/"representative_xyz"
    out.mkdir(parents=True,exist_ok=True)
    for (group,medium), sub in df.groupby(["group","medium"]):
        row=sub.sort_values("energy_hartree").iloc[0]
        source=(PRIVATE/"qm"/group/f"conf_{int(row.conformer):03d}"/medium/"xtbopt.xyz")
        dest=out/f"{group}_{medium}_lowest.xyz"
        lines=source.read_text(encoding="utf-8").splitlines()
        lines[1]=f"GFN2-xTB optimized; group={group}; medium={medium}; charge=+1; E={row.energy_hartree:.10f} Eh"
        dest.write_text("\n".join(lines)+"\n",encoding="utf-8")


def save_figure(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    primary=df[df.group.str.startswith("sulfonyl")]
    fig,axes=plt.subplots(1,2,figsize=(9,4.0),sharey=True)
    for ax,medium in zip(axes,("gas","dmf")):
        for group,marker in (("sulfonyl_meta","o"),("sulfonyl_ortho","s")):
            sub=primary[(primary.group==group)&(primary.medium==medium)&
                        primary.relative_kcal_mol.ge(0)]
            ax.scatter(sub["min_NO_A"],sub["relative_kcal_mol"],s=65,
                       marker=marker,label=group.replace("_"," "),alpha=.9)
            for _,r in sub.iterrows():
                if r.hb_strict:
                    ax.annotate("N–H···O",(r.min_NO_A,r.relative_kcal_mol),
                                xytext=(6,6),textcoords="offset points",fontsize=7)
        ax.set(xlabel="Closest N···O distance (Å)",
               title=("Gas phase" if medium=="gas" else "ALPB(DMF)"))
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Relative GFN2-xTB energy (kcal/mol)")
    axes[0].legend(loc="upper left",fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS/"sulfonyl_conformer_landscape.png",dpi=180)
    plt.close(fig)


def main():
    RESULTS.mkdir(parents=True,exist_ok=True)
    df,mols=read_data()
    df.to_csv(RESULTS/"xtb_conformers.csv",index=False,float_format="%.9f")
    mini=[]
    for (group,medium),sub in df.groupby(["group","medium"]):
        mi=sub.loc[sub.energy_hartree.idxmin()]
        low=sub[sub.relative_kcal_mol<=3.]
        mini.append(dict(group=group,medium=medium,n_runs=len(sub),
                         best_conformer=int(mi.conformer),
                         best_energy_Eh=float(mi.energy_hartree),
                         best_N_O_A=float(mi.min_NO_A),
                         best_H_O_A=float(mi.min_HO_A),
                         best_N_H_O_deg=float(mi.max_NHO_deg),
                         best_hbond=bool(mi.hb_strict),
                         low_3kcal_starts=len(low),
                         low_3kcal_hbond_starts=int(low.hb_strict.sum())))
    pd.DataFrame(mini).to_csv(RESULTS/"conformer_summary.csv",index=False,float_format="%.7f")
    save_representatives(df,mols)
    save_figure(df)
    manifest=dict(
        source="Study 1 frozen shortlist, CID groups",
        method="RDKit ETKDGv3/MMFF94s -> GFN2-xTB geometry optimization",
        xtb_version="6.7.1", total_runs=len(df),
        n_successful=int(df.success.sum()),
        environments=["gas","ALPB(dmf)"],
        charge=1, unpaired_electrons=0,
        hartree_to_kcal=HARTREE_TO_KCAL,
        intramolecular_contact="same-cation N-H..O, N-O<3.3 A, H-O<2.5 A, angle N-H-O>=130deg",
        sampling_note="Number of optimized starting structures; do not count them as unique energy minima or Boltzmann populations",
        raw_output_path="external scratch directory configured via STUDY1_QM_WORKDIR",
        calculation_level="semiempirical quantum chemistry, not DFT",
        scope="gas-phase isolated cations and implicitly solvated cations; no counterions, explicit solvent or Pb-I environment",
        limitation="ETKDG/MMFF starting conformers; conformational sampling is provisional, no CREST exhaustive search")
    (RESULTS/"xtb_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(pd.DataFrame(mini).to_string(index=False, float_format=lambda x:f"{x:.4f}"))
    print("Generated",len(df),"successful xTB optimizations from initial MMFF starts")


if __name__=="__main__":
    main()
