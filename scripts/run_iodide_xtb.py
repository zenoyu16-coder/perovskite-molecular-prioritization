"""GFN2-xTB ion-pair pilot: varied cation conformers and iodide binding sites.

Use with WSL Python 3.12+, on the local filesystem/scratch. All raw
calculations stay outside Git. Iodide and cation charges sum to neutral.
The computations assess ion-pair rearrangement, not PbI4 adsorption.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse,csv,json,os,re,shutil,subprocess,time
from pathlib import Path

ROOT=Path(os.environ.get("STUDY1_QM_WORKDIR","/mnt/d/AI4S/study1_quantum_pilot"))
BASE=ROOT / os.environ.get("STUDY1_IONPAIR_SUBDIR", "iodide_association")
XTB=Path(os.environ.get("STUDY1_XTB_BIN","/mnt/d/AI4S/tools/xtb-6.7.1/xtb-dist/bin/xtb"))
RE_E=re.compile(r"TOTAL ENERGY\s+([-+\d.]+)\s+Eh")

def run(job,medium):
    p=BASE/job["input_file"]
    work=p.parent/medium
    work.mkdir(parents=True,exist_ok=True)
    summary=work/"run.json"
    if summary.exists():
        r=json.loads(summary.read_text(encoding="utf-8"))
        if r.get("success") and (work/"xtbopt.xyz").is_file():
            return r
    shutil.copy2(p,work/"input.xyz")
    cmd=[str(XTB),"input.xyz","--gfn","2","--chrg","0",
         "--opt","tight","--cycles","350"]
    if medium!="gas":cmd.extend(["--alpb",medium])
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1",OMP_STACKSIZE="128M")
    start=time.monotonic()
    try:
        process=subprocess.run(cmd,cwd=work,env=env,text=True,
                               capture_output=True,timeout=180)
        output=process.stdout+"\n"+process.stderr
        (work/"xtb.stdout").write_text(output,encoding="utf-8")
        energies=RE_E.findall(output)
        e=float(energies[-1]) if energies else None
        ok=bool(process.returncode==0 and e is not None and
                "normal termination of xtb" in output and
                (work/"xtbopt.xyz").is_file())
        r=dict(group=job["group"],conformer=int(job["conformer"]),
               site=int(job["site"]),medium=medium,
               energy_hartree=e,success=ok,exit_code=process.returncode,
               wall_seconds=round(time.monotonic()-start,3),
               optimized_xyz=str((work/"xtbopt.xyz").relative_to(ROOT)) if ok else None,
               raw_stdout=str((work/"xtb.stdout").relative_to(ROOT)),
               method="GFN2-xTB 6.7.1",
               model="(+1 cation, -1 iodide); total charge 0",
               environment=("ALPB("+medium+")" if medium!="gas" else "gas"))
    except subprocess.TimeoutExpired:
        r=dict(group=job["group"],conformer=int(job["conformer"]),
               site=int(job["site"]),medium=medium,
               energy_hartree=None,success=False,exit_code=124,
               wall_seconds=round(time.monotonic()-start,3),
               failure_reason="timeout")
    summary.write_text(json.dumps(r,indent=2)+"\n",encoding="utf-8")
    return r

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--media",nargs="+",default=["dmf","gas"])
    p.add_argument("--workers",type=int,default=4)
    a=p.parse_args()
    jobs=list(csv.DictReader((BASE/"pair_inputs.csv").open(encoding="utf-8")))
    tasks=[(j,m) for j in jobs for m in a.media]
    out=[]
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futures={ex.submit(run,j,m):(j,m) for j,m in tasks}
        for f in as_completed(futures):
            r=f.result()
            out.append(r)
            print(f"{r['group']:18s} c{r['conformer']:02d} "
                  f"s{r['site']:d} {r['medium']:3s}: {r['success']} "
                  f"{r.get('energy_hartree')}",flush=True)
    (BASE/"run_summary.json").write_text(json.dumps(dict(
        n_jobs=len(out),n_success=sum(v["success"] for v in out),
        n_failed=sum(not v["success"] for v in out),
        media=a.media,method="GFN2-xTB 6.7.1",
        status="optimize gas and ALPB(dmf) ion-pair structures"),indent=2)+"\n",
        encoding="utf-8")
    print("ALL",len(out),"SUCCESS",sum(v["success"] for v in out))

if __name__=="__main__":
    main()
