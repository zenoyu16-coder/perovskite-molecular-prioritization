"""Run GFN2-xTB optimizations for Study 1 positional-isomer conformers.

Intended for the connected WSL Ubuntu environment, with the vetted local xTB
binary. Source files and raw outputs are kept on D:, outside the Git repository.
This program performs actual electronic-structure calculations.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT = Path(os.environ.get("STUDY1_QM_WORKDIR", str(Path(__file__).resolve().parents[1]/"data/external/study1_qm")))
XTB = Path(os.environ.get("STUDY1_XTB_BIN", shutil.which("xtb") or "xtb"))
XTB_VER = "6.7.1"
ENERGY_RE = re.compile(r"TOTAL ENERGY\s+([-+\d.]+)\s+Eh")


def run_one(source: dict, solvent: str, force: bool):
    name = source["group"]
    rank = int(source["conformer"])
    work = ROOT / "qm" / name / f"conf_{rank:03d}" / solvent
    work.mkdir(parents=True, exist_ok=True)
    summary_path = work / "run.json"
    if summary_path.exists() and not force:
        existing = json.loads(summary_path.read_text(encoding="utf-8"))
        if existing.get("success"):
            return existing
    initial = ROOT / source["xyz_path"]
    if not initial.is_file():
        raise FileNotFoundError(initial)
    shutil.copy2(initial, work / "input.xyz")
    cmd = [str(XTB), "input.xyz", "--gfn", "2", "--chrg", "1",
           "--opt", "tight", "--cycles", "300"]
    if solvent != "gas":
        cmd.extend(["--alpb", solvent])
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", OMP_STACKSIZE="128M")
    start = time.monotonic()
    try:
        result = subprocess.run(cmd, cwd=work, env=env,
                                capture_output=True, text=True, timeout=240)
        out = result.stdout + "\n" + result.stderr
        (work / "xtb.stdout").write_text(out, encoding="utf-8")
        matches = ENERGY_RE.findall(out)
        e = float(matches[-1]) if matches else None
        ok = (result.returncode == 0 and
              "normal termination of xtb" in out and
              (work / "xtbopt.xyz").exists() and e is not None)
        status = dict(group=name, pubchem_cid=int(source["pubchem_cid"]),
                      conformer=rank, solvent=solvent, xtb_version=XTB_VER,
                      formal_charge=1, electron_model="GFN2-xTB",
                      solvation=("ALPB:" + solvent if solvent != "gas" else "none"),
                      energy_hartree=e, success=ok,
                      exit_code=result.returncode,
                      runtime_seconds=round(time.monotonic() - start, 3),
                      final_xyz=str((work / "xtbopt.xyz").relative_to(ROOT)) if ok else None,
                      stdout_file=str((work / "xtb.stdout").relative_to(ROOT)),
                      failure_reason="" if ok else "nonzero exit, missing optimized xyz or energy")
    except subprocess.TimeoutExpired:
        status = dict(group=name, pubchem_cid=int(source["pubchem_cid"]),
                      conformer=rank, solvent=solvent, xtb_version=XTB_VER,
                      formal_charge=1, electron_model="GFN2-xTB",
                      solvation=("ALPB:" + solvent if solvent != "gas" else "none"),
                      energy_hartree=None, success=False, exit_code=124,
                      runtime_seconds=round(time.monotonic()-start, 3),
                      failure_reason="timeout")
    summary_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    return status


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--solvents", nargs="+", default=["gas", "dmf"])
    args = parser.parse_args()
    if not XTB.is_file():
        raise FileNotFoundError(XTB)
    source = list(csv.DictReader((ROOT / "conformers.csv").open(encoding="utf-8")))
    jobs = [(row, s) for row in source for s in args.solvents]
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, row, solvent, args.force):
                   (row["group"], row["conformer"], solvent)
                   for row, solvent in jobs}
        for future in as_completed(futures):
            try:
                r = future.result()
                results.append(r)
                print(f"{r['group']:20s} {r['conformer']:>3} "
                      f"{r['solvent']:4s} success={r['success']} "
                      f"E={r.get('energy_hartree')}", flush=True)
            except Exception as exc:
                print("ERROR", futures[future], str(exc), flush=True)
                raise
    result_dir = ROOT / "qm"
    (result_dir / "run_summary.json").write_text(json.dumps({
        "model": "GFN2-xTB", "xtb_version": XTB_VER,
        "solvation": args.solvents, "n_jobs": len(results),
        "n_success": sum(r["success"] for r in results),
        "n_failed": sum(not r["success"] for r in results),
        "total_runtime_seconds": round(sum(r["runtime_seconds"] for r in results), 2)
    }, indent=2) + "\n", encoding="utf-8")
    print("ALL JOBS", len(results), "SUCCESS", sum(r["success"] for r in results))


if __name__ == "__main__":
    main()
