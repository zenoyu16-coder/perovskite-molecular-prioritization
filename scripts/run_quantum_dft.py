"""PBE0/def2-SVP single-point cross-check of two CREST-selected cations.

Runs under WSL with PySCF 2.11. The +1 ionic species are evaluated using
ddCOSMO with epsilon=36.7 (approximate DMF dielectric). Geometry is fixed
at the GFN2-xTB/ALPB(DMF) optimized structure; no DFT geometry optimization.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import time

from pyscf import dft, gto, lib, solvent
import pyscf

ROOT = Path(__file__).resolve().parents[1]
BASE = Path(os.environ.get("STUDY1_QM_WORKDIR", str(ROOT/"data/external/study1_qm")))
BASIS = "def2-svp"
XC = "pbe0"
EPS_DMF = 36.7


def first_xyz(path):
    lines=path.read_text(encoding="utf-8").splitlines()
    n=int(lines[0].strip())
    atoms=[]
    for ln in lines[2:2+n]:
        parts=ln.split()
        atoms.append([parts[0], [float(x) for x in parts[1:4]]])
    if len(atoms)!=n: raise ValueError("truncated XYZ")
    return atoms


def run(group, force=False):
    inp=ROOT/"results/quantum_pilot/dft_inputs"/f"{group}.xyz"
    work=BASE/"dft"/group
    work.mkdir(exist_ok=True,parents=True)
    summary=work/"pbe0_def2svp_dmf.json"
    if summary.exists() and not force:
        print(summary.read_text(encoding="utf-8"))
        return
    coords=first_xyz(inp)
    lib.num_threads(4)
    mol=gto.M(atom=coords, basis=BASIS, charge=1, spin=0,
              unit="Angstrom", verbose=4, max_memory=3400)
    if mol.nelectron%2 !=0: raise ValueError("Unexpected electron parity")
    mf=dft.RKS(mol).density_fit()
    mf.xc=XC
    mf.grids.level=3
    mf.conv_tol=1e-8
    mf.max_cycle=100
    mf.chkfile=str(work/"scf.chk")
    mf=solvent.ddCOSMO(mf)
    mf.with_solvent.eps=EPS_DMF
    t=time.perf_counter()
    with (work/"pyscf.stdout").open("w",encoding="utf-8") as out:
        import contextlib
        with contextlib.redirect_stdout(out):
            e=mf.kernel()
    elapsed=time.perf_counter()-t
    rec=dict(group=group,method="PBE0/def2-SVP density-fitted DFT single-point",
             pyscf_version=pyscf.__version__,charge=1,spin=0,
             n_atoms=len(coords),n_electrons=mol.nelectron,
             geometry="lowest CREST GFN2-xTB/ALPB(DMF) conformer; fixed",
             solvent="ddCOSMO dielectric=36.7 (DMF approximation)",
             energy_hartree=float(e),converged=bool(mf.converged),
             wall_seconds=round(elapsed,2),memory_limit_mb=3400,
             basis=BASIS,xc=XC,grid_level=3)
    summary.write_text(json.dumps(rec,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(rec,indent=2))
    if not mf.converged: raise RuntimeError("DFT SCF failed to converge")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--group",choices=["sulfonyl_meta","sulfonyl_ortho"],required=True)
    p.add_argument("--force",action="store_true")
    a=p.parse_args()
    run(a.group,a.force)
