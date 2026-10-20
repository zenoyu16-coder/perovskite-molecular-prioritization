"""Recompute DFT iodide association double differences from released raw energies.

All four total energies needed for a double difference are committed as
CSV: two isolated cations and two contact ion pairs. No private scratch
directory is required to run this analysis.
"""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "iodide_pilot"
HARTREE_TO_KCAL = 627.509474


def compute_energy_comparison():
    cations = pd.read_csv(ROOT / "results" / "quantum_pilot" /
                          "dft_singlepoints.csv").set_index("group")
    pairs = pd.read_csv(OUT / "dft_pair_singlepoints.csv").set_index("group")
    keys = ("sulfonyl_meta", "sulfonyl_ortho")
    for group in keys:
        cat, pair = cations.loc[group], pairs.loc[group]
        if not bool(cat["converged"]) or not bool(pair["converged"]):
            raise RuntimeError(f"DFT convergence missing: {group}")
        if (str(cat["xc"]).lower() != str(pair["xc"]).lower() or
            str(cat["basis"]).lower() != str(pair["basis"]).lower()):
            raise RuntimeError(f"Mismatch in DFT functional or basis: {group}")
        if (int(cat["charge"]) != 1 or int(pair["charge"]) != 0 or
            int(pair["n_atoms"]) != int(cat["n_atoms"]) + 1):
            raise RuntimeError(f"Unexpected ion-pair composition: {group}")
        if (pair["solvation"] != "ddCOSMO(DMF;eps=36.7)" or
            pair["iodine_ecp"] != "def2-svp ECP28" or
            int(pair["n_atoms"]) != 25 or
            pair["input_xyz"] != f"results/iodide_pilot/{group}_dmf_lowest.xyz"):
            raise RuntimeError(f"Iodine ECP / solvent / structure mismatch: {group}")
        if not (ROOT / pair["input_xyz"]).is_file():
            raise FileNotFoundError(pair["input_xyz"])
    a,b = (pairs.loc[key] for key in keys)
    c,d = (cations.loc[key] for key in keys)
    delta_pair = (float(b["energy_hartree"]) -
                  float(a["energy_hartree"])) * HARTREE_TO_KCAL
    delta_cation = (float(d["energy_hartree"]) -
                    float(c["energy_hartree"])) * HARTREE_TO_KCAL
    return delta_pair, delta_cation, delta_pair - delta_cation


def main():
    delta_pair, delta_cation, delta_assoc = compute_energy_comparison()
    xtb = json.loads((OUT / "association_differences.json").read_text(encoding="utf-8"))
    xtb_dmf = next(item for item in xtb["energy_differences"] if item["medium"] == "dmf")
    xtb_delta = float(xtb_dmf["ortho_minus_meta_association_kcal_mol"])
    report = {
        "model": "(+1 organic cation and -1 iodide) ion pair in implicit DMF",
        "relative_association_energy_definition":
            "(ortho pair-meta pair)-(ortho cation-meta cation)",
        "pbe0_def2svp_ddcosmo": {
            "delta_pair_kcal_mol": delta_pair,
            "delta_cation_kcal_mol": delta_cation,
            "delta_association_kcal_mol": delta_assoc,
            "n_ion_pair_scf_converged": 2,
            "iodine_ECP": "def2-svp ECP28",
        },
        "gfn2_xtb_alpb_dmf": {
            "delta_association_kcal_mol": xtb_delta
        },
        "agreement_on_sign": bool(delta_assoc * xtb_delta > 0),
        "interpretation":
            "Positive sign: stronger electronic iodide association for meta; negative: ortho",
        "fixed_geometry_limitation":
            "DFT single-point energies at GFN2-xTB optimized geometries",
        "uncertainty":
            "No association free energies, explicit solvent, basis convergence, BSSE correction or periodic Pb-I framework",
        "raw_energy_sources": [
            "results/quantum_pilot/dft_singlepoints.csv",
            "results/iodide_pilot/dft_pair_singlepoints.csv",
        ],
    }
    (OUT / "dft_association_check.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"DFT pair delta: {delta_pair:.7f} kcal/mol")
    print(f"DFT cation delta: {delta_cation:.7f} kcal/mol")
    print(f"DFT association double difference: {delta_assoc:.7f} kcal/mol")


if __name__ == "__main__":
    main()
