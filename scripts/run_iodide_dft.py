"""Minimal PBE0/ddCOSMO single-point cross-check of iodide association.

Cation-ion-pair geometries are lowest GFN2-xTB/ALPB(DMF) minima selected
from 42 ion-pair initial geometries. Input coordinates are published and
do not change during these single-point DFT calculations.
"""
from __future__ import annotations
import argparse,contextlib,json,os,time
from pathlib import Path
from pyscf import dft,gto,lib,solvent
import pyscf

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/iodide_pilot"
WORK=Path(os.environ.get("STUDY1_QM_WORKDIR","/mnt/d/AI4S/study1_quantum_pilot"))
XC="pbe0"
BASIS="def2-svp"
EPS=36.7

def read_xyz(p):
    t=p.read_text(encoding="utf-8").splitlines()
    n=int(t[0]);atoms=[]
    for x in t[2:2+n]:
        f=x.split()
        atoms.append([f[0],[float(v) for v in f[1:4]]])
    if len(atoms)!=n or atoms[-1][0].upper()!="I":
        raise ValueError("Invalid cation+iodide geometry")
    return atoms

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--group",required=True,choices=["sulfonyl_meta","sulfonyl_ortho"])
    arg.add_argument("--force",action="store_true")
    a=arg.parse_args()
    coord=OUT/f"{a.group}_dmf_lowest.xyz"
    run=WORK/"iodide_association"/"dft"/a.group
    run.mkdir(parents=True,exist_ok=True)
    result=run/"pbe0_pair.json"
    if result.is_file() and not a.force:
        print(result.read_text(encoding="utf-8"))
        return
    atoms=read_xyz(coord)
    lib.num_threads(4)
    mol=gto.M(atom=atoms,basis=BASIS,ecp={"I":"def2-svp"},
        charge=0,spin=0,unit="Angstrom",verbose=3,max_memory=3400)
    # I: ECP28, 25 electrons explicitly for neutral iodine; anion adds one
    if mol.nelectron%2: raise ValueError("Open-shell system unexpected")
    calc=dft.RKS(mol).density_fit()
    calc.xc=XC
    calc.grids.level=3
    calc.conv_tol=1e-8
    calc.max_cycle=100
    calc.chkfile=str(run/"pair.chk")
    calc=solvent.ddCOSMO(calc)
    calc.with_solvent.eps=EPS
    start=time.monotonic()
    with (run/"pyscf.stdout").open("w",encoding="utf-8") as f:
        with contextlib.redirect_stdout(f):
            e=calc.kernel()
    meta=dict(group=a.group,energy_hartree=float(e),
        converged=bool(calc.converged),method="PBE0/def2-SVP ddCOSMO",
        iodine_ecp="def2-svp ECP28",xc=XC,basis=BASIS,
        eps=EPS,n_atoms=len(atoms),n_effective_electrons=mol.nelectron,
        charge=0,spin=0,solvent="ddCOSMO eps36.7",
        geometry="GFN2-xTB ALPB(DMF) optimized contact ion pair (not DFT optimized)",
        input_xyz=str(coord.relative_to(ROOT)),
        wall_seconds=round(time.monotonic()-start,2))
    result.write_text(json.dumps(meta,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(meta,indent=2))
    if not calc.converged:raise RuntimeError("DFT SCF failed to converge")

if __name__=="__main__":main()
