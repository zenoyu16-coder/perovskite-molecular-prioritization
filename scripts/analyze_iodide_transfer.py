"""Held-out fluoro/methoxy positional-isomer ion-pair comparison.

Uses the same GFN2-xTB ALPB(DMF) calculation and initial-site protocol
as the sulfonyl pair. The chemical pair and qualitative prediction were
chosen before inspecting these ion-pair energies.
"""
from __future__ import annotations
import csv,json,os,shutil
from pathlib import Path
import pandas as pd
from analyze_iodide_pairs import bond_summary

ROOT=Path(__file__).resolve().parents[1]
RAW=Path(os.environ.get("STUDY1_QM_WORKDIR","D:/AI4S/study1_quantum_pilot"))
PAIR=RAW/"iodide_association_fluoromethoxy"
OUT=ROOT/"results/iodide_pilot"
KCAL=627.509474
groups=(("fluoromethoxy_A",40786958),("fluoromethoxy_B",24721203))

def crest_min(group):
    text=(ROOT/"results/quantum_pilot/crest_ensembles"/f"{group}.xyz").read_text().splitlines()
    n=int(text[0])
    return min(float(text[i+1]) for i in range(0,len(text),n+2) if i+1<len(text))

def main():
    sr={int(r["pubchem_cid"]):r for r in csv.DictReader(
        (ROOT/"results/clean_study1/shortlist.csv").open(encoding="utf-8"))}
    reference = json.loads((OUT/"iodide_reference.json").read_text(encoding="utf-8"))
    iodide = reference["results"]["dmf"]
    if not iodide["converged"] or iodide["charge"] != -1 or iodide["solvation"] != "ALPB(DMF)":
        raise ValueError("Incompatible isolated-iodide reference")
    ion = float(iodide["energy_hartree"])
    runs=[]
    stats=[]
    for group,cid in groups:
        source=sr[cid]["modeled_cation_smiles"]
        mono=crest_min(group)
        subset=[]
        for item in (PAIR/group).rglob("run.json"):
            rec=json.loads(item.read_text(encoding="utf-8"))
            if not rec["success"] or rec["medium"]!="dmf":raise RuntimeError(item)
            original_xyz=RAW/rec["optimized_xyz"]
            e=rec["energy_hartree"]
            charges=[float(x) for x in (original_xyz.parent/"charges").read_text().split()]
            data=dict(group=group,cid=cid,conformer=rec["conformer"],
                      site=rec["site"],pair_Eh=e,cation_Eh=mono,
                      association_kcal_mol=(e-mono-ion)*KCAL,
                      I_partial_charge=charges[-1],
                      raw_xyz=rec["optimized_xyz"],
                      **bond_summary(source,original_xyz))
            subset.append(data)
            runs.append(data)
        best=min(subset,key=lambda r:r["pair_Eh"])
        near=[r for r in subset if (r["pair_Eh"]-best["pair_Eh"])*KCAL<=0.5]
        summary=dict(group=group,cid=cid,n_success=len(subset),
                n_distinct_CREST_conformers=len(set(r["conformer"] for r in subset)),
                n_within_0_5kcal=len(near),
                cation_min_Eh=mono,pair_min_Eh=best["pair_Eh"],
                association_kcal_mol=best["association_kcal_mol"],
                best_I_partial_charge=best["I_partial_charge"],
                best_N_I_A=best["N_I_A"],
                best_internal_N_O_A=best["internal_N_O_A"],
                best_internal_N_H_O_deg=best["internal_N_H_O_max_deg"],
                best_internal_Hbond=best["internal_Hbond"])
        stats.append(summary)
        xyz=RAW/best["raw_xyz"]
        lines=xyz.read_text(encoding="utf-8").splitlines()
        lines[1]=f"GFN2-xTB/ALPB(DMF) {group}, charge 0, E={best['pair_Eh']:.10f} Eh"
        (OUT/f"{group}_dmf_lowest.xyz").write_text("\n".join(lines)+"\n",encoding="utf-8")
    OUT.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(runs).to_csv(OUT/"fluoro_ion_pairs.csv",index=False,float_format="%.9f")
    pd.DataFrame(stats).to_csv(OUT/"fluoro_ion_pair_summary.csv",index=False,float_format="%.9f")
    a,b=stats
    delta=a["association_kcal_mol"]-b["association_kcal_mol"]
    report=dict(
        test="Fluoromethoxybenzylammonium A vs B: external preselected rank-tied regioisomer challenge",
        environment="GFN2-xTB 6.7.1 ALPB(DMF), monocation+I- to neutral ion pair",
        prediction="Isomer A forms a nearby intramolecular ammonium-O contact in the free cation and should have less favorable iodide association than B",
        a_minus_b_electronic_association_kcal_mol=delta,
        expected_positive=True,observed_positive=delta>0,
        n_optimized_ion_pairs=len(runs),
        pose_sampling_counts={g["group"]:g["n_success"] for g in stats},
        limitations=["One continuum solvent model and only GFN2-xTB for transfer pair",
          "Unequal accepted starting-pose counts (A=7, B=33), so lowest-state search may favor B",
          "Ammonium-O contact in A weakens upon external iodide coordination and may fail angular contact criterion",
          "Absolute association enthalpies/free energies not obtained; no Pb-I framework"],
        method="Lowest optimized ion-pair electronic energy minus lowest CREST free-cation energy and isolated I- electronic energy, same xTB/solvent model"
    )
    (OUT/"fluoro_transfer.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(pd.DataFrame(stats).to_string(index=False,float_format=lambda x:f"{x:.4f}"))
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
