"""Export raw converged ion-pair PBE0 energies from local PySCF JSON records.

This explicit export step is used after rerunning run_iodide_dft.py on
the published optimized ion-pair XYZs. The resulting CSV is publicly
readable and is the input to analyze_iodide_dft.py and unit tests.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "results" / "iodide_pilot" / "dft_pair_singlepoints.csv"


def collect(scratch: Path):
    rows = []
    for group in ("sulfonyl_meta", "sulfonyl_ortho"):
        source = scratch / "iodide_association" / "dft" / group / "pbe0_pair.json"
        data = json.loads(source.read_text(encoding="utf-8"))
        expected_input = f"results/iodide_pilot/{group}_dmf_lowest.xyz"
        if not data["converged"] or data["input_xyz"] != expected_input:
            raise ValueError(f"SCF or geometry provenance failed for {group}")
        if (data["xc"].lower() != "pbe0" or
            data["basis"].lower() != "def2-svp" or
            data["iodine_ecp"] != "def2-svp ECP28" or
            float(data["eps"]) != 36.7 or
            data["charge"] != 0 or data["spin"] != 0 or data["n_atoms"] != 25):
            raise ValueError(f"Unexpected DFT parameters for {group}")
        rows.append({
            "group": group,
            "energy_hartree": repr(float(data["energy_hartree"])),
            "converged": True,
            "xc": data["xc"].lower(),
            "basis": data["basis"].lower(),
            "solvation": "ddCOSMO(DMF;eps=36.7)",
            "iodine_ecp": data["iodine_ecp"],
            "charge": 0,
            "n_atoms": 25,
            "input_xyz": expected_input,
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scratch", type=Path,
                        default=Path(os.environ.get("STUDY1_QM_WORKDIR",
                                                    str(ROOT / "data" / "external" / "study1_qm"))))
    args = parser.parse_args()
    rows = collect(args.scratch)
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    with TARGET.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for row in rows:
        print(row["group"], row["energy_hartree"], "Eh")


if __name__ == "__main__":
    main()
