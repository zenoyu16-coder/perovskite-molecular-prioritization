"""Compute isolated iodide (-1) reference energies used by the ion-pair pilot.

Run GFN2-xTB 6.7.1 single points on one iodine atom with charge -1,
in the gas phase and ALPB(DMF). Separate-ion and ion-pair calculations
must use the same GFN2-xTB/solvation setup.

The generated input and structured reference energies are kept in the
repository; raw xTB logs and temporary data remain in scratch storage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PUBLIC = PROJECT_ROOT / "results" / "iodide_pilot"
DEFAULT_SCRATCH = PROJECT_ROOT / "data" / "external" / "study1_qm"
IODIDE_XYZ = "1\nIodide anion; charge -1\nI    0.000000000  0.000000000  0.000000000\n"
HARTREE_PATTERN = re.compile(r"TOTAL ENERGY\s+([-+]?\d+(?:\.\d+)?)\s+Eh", re.IGNORECASE)


def extract_energy(output: str) -> float:
    """Require a normally terminated calculation with a final total energy."""
    if "normal termination of xtb" not in output.lower():
        raise RuntimeError("xTB did not report normal termination")
    values = HARTREE_PATTERN.findall(output)
    if not values:
        raise RuntimeError("No total energy in xTB output")
    return float(values[-1])


def calculate_reference(executable: str, scratch: Path, medium: str) -> dict:
    if medium not in ("gas", "dmf"):
        raise ValueError("Only the two original reference environments are supported")
    work = scratch / "iodide_association" / "isolated_iodide"
    work.mkdir(parents=True, exist_ok=True)
    structure = work / "iodide.xyz"
    structure.write_text(IODIDE_XYZ, encoding="utf-8")
    command = [executable, "iodide.xyz", "--gfn", "2", "--chrg", "-1", "--sp"]
    if medium == "dmf":
        command.extend(["--alpb", "dmf"])
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1")
    result = subprocess.run(command, cwd=work, env=env, text=True,
                            capture_output=True, timeout=120, check=False)
    output = result.stdout + "\n" + result.stderr
    (work / f"iodide_{medium}.stdout").write_text(output, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"xTB {medium} reference failed, exit={result.returncode}")
    energy = extract_energy(output)
    return dict(energy_hartree=energy, charge=-1, gfn=2,
                solvation="ALPB(DMF)" if medium == "dmf" else "gas",
                optimized=False, converged=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch", type=Path,
                        default=Path(os.environ.get("STUDY1_QM_WORKDIR", str(DEFAULT_SCRATCH))))
    parser.add_argument("--xtb", default=os.environ.get("STUDY1_XTB_BIN", shutil.which("xtb") or "xtb"))
    parser.add_argument("--check-only", action="store_true",
                        help="Verify published reference records without launching xTB")
    args = parser.parse_args()
    output_file = PUBLIC / "iodide_reference.json"
    if args.check_only:
        data = json.loads(output_file.read_text(encoding="utf-8"))
        assert data["input_sha256"] == hashlib.sha256(IODIDE_XYZ.encode()).hexdigest()
        assert set(data["results"]) == {"gas", "dmf"}
        assert all(data["results"][medium]["converged"] for medium in ("gas", "dmf"))
        print("Verified published iodide reference metadata")
        return
    result = {
        "method": "GFN2-xTB 6.7.1 single-point, charge -1",
        "geometry": "single iodine atom at Cartesian origin, Angstrom",
        "input_sha256": hashlib.sha256(IODIDE_XYZ.encode()).hexdigest(),
        "electronically_isolated_atom": True,
        "results": {
            medium: calculate_reference(args.xtb, args.scratch, medium)
            for medium in ("dmf", "gas")
        },
        "note": "Reference energies are used with the matching xTB/solvent model; gas results are not used for the main solution comparison."
    }
    PUBLIC.mkdir(parents=True, exist_ok=True)
    (PUBLIC / "iodide.xyz").write_text(IODIDE_XYZ, encoding="utf-8")
    output_file.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
