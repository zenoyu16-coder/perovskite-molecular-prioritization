"""Generate charge-consistent conformers of two matched positional-isomer pairs.

Preparation uses ETKDGv3 and MMFF94s only. No QM energies are inferred
at this stage. Raw coordinates remain on the local D: drive.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign

ROOT = Path(__file__).resolve().parents[1]
PRIVATE = Path(os.environ.get("STUDY1_QM_WORKDIR", str(ROOT/"data/external/study1_qm")))
RESULTS = ROOT / "results" / "quantum_pilot"
CANDIDATES = {
    2760937: ("sulfonyl_meta", 16),
    16218335: ("sulfonyl_ortho", 16),
    40786958: ("fluoromethoxy_A", 8),
    24721203: ("fluoromethoxy_B", 8),
}
SEED = 20261020
N_EMBED = 100
RMS_CUTOFF = 0.65
ENERGY_WINDOW_KCAL = 18.0


def hbond_features(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    nids = [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == "N" and a.GetFormalCharge() > 0
            and any(x.GetSymbol() == "H" for x in a.GetNeighbors())]
    aids = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O"]
    if not nids or not aids:
        return dict(min_NO_A=None, min_HO_A=None, max_NHO_deg=None, hb_strict=False)
    rec = []
    for n in nids:
        hs = [h.GetIdx() for h in mol.GetAtomWithIdx(n).GetNeighbors()
              if h.GetSymbol() == "H"]
        for o in aids:
            no = conf.GetAtomPosition(n).Distance(conf.GetAtomPosition(o))
            for h in hs:
                nh = conf.GetAtomPosition(n) - conf.GetAtomPosition(h)
                oh = conf.GetAtomPosition(o) - conf.GetAtomPosition(h)
                denom = nh.Length() * oh.Length()
                cosine = max(-1, min(1, nh.DotProduct(oh) / denom))
                ang = math.degrees(math.acos(cosine))
                rec.append((no, oh.Length(), ang))
    return dict(
        min_NO_A=min(v[0] for v in rec),
        min_HO_A=min(v[1] for v in rec),
        max_NHO_deg=max(v[2] for v in rec),
        hb_strict=any(no < 3.3 and ho < 2.5 and angle >= 130
                      for no, ho, angle in rec)
    )


def xyz_text(mol, cid, title):
    cf = mol.GetConformer(cid)
    rows = [str(mol.GetNumAtoms()), title]
    for atom in mol.GetAtoms():
        p = cf.GetAtomPosition(atom.GetIdx())
        rows.append(f"{atom.GetSymbol():2s}  {p.x: .9f} {p.y: .9f} {p.z: .9f}")
    return "\n".join(rows) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nembed", type=int, default=N_EMBED)
    parser.add_argument("--limit-primary", type=int, default=16)
    parser.add_argument("--limit-secondary", type=int, default=8)
    args = parser.parse_args()
    input_path = ROOT / "results/clean_study1/shortlist.csv"
    chem_rows = {int(row["pubchem_cid"]): row
                 for row in csv.DictReader(input_path.open(encoding="utf-8"))}
    RESULTS.mkdir(parents=True, exist_ok=True)
    PRIVATE.mkdir(parents=True, exist_ok=True)
    rows = []
    groups = []
    for cid, (name, _) in CANDIDATES.items():
        if cid not in chem_rows:
            raise ValueError(f"CID {cid} absent from frozen shortlist")
        smiles = chem_rows[cid]["modeled_cation_smiles"]
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if mol is None or Chem.GetFormalCharge(mol) != 1:
            raise ValueError(f"Incorrect cation charge for {cid}")
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            raise ValueError(f"MMFF94s parameters missing for {cid}")
        ps = AllChem.ETKDGv3()
        ps.randomSeed = SEED + (cid % 100000)
        ps.pruneRmsThresh = 0.4
        ps.numThreads = 4
        confids = list(AllChem.EmbedMultipleConfs(mol, numConfs=args.nembed, params=ps))
        if not confids:
            raise RuntimeError(f"ETKDG failed for {cid}")
        results = list(AllChem.MMFFOptimizeMoleculeConfs(
            mol, numThreads=4, maxIters=500, mmffVariant="MMFF94s"))
        viable = [(cid0, float(e)) for cid0, (status, e) in zip(confids, results)
                  if status == 0 and math.isfinite(e)]
        if not viable:
            raise RuntimeError(f"MMFF minimization failed for {cid}")
        viable.sort(key=lambda v: v[1])
        base = viable[0][1]
        heavy = [(a.GetIdx(), a.GetIdx()) for a in mol.GetAtoms()
                 if a.GetAtomicNum() > 1]
        limit = args.limit_primary if name.startswith("sulfonyl") else args.limit_secondary
        selected = []
        for m_idx, energy in viable:
            if energy > base + ENERGY_WINDOW_KCAL:
                continue
            if any(rdMolAlign.GetAlignmentTransform(
                    mol, mol, prbCid=m_idx, refCid=k, atomMap=heavy)[0]
                   < RMS_CUTOFF for k, _ in selected):
                continue
            selected.append((m_idx, energy))
            if len(selected) >= limit:
                break
        group_dir = PRIVATE / name
        group_dir.mkdir(parents=True, exist_ok=True)
        for rank, (idx, energy) in enumerate(selected, 1):
            title = f"CID={cid} conformer={rank} formal_charge=+1 mmff94s={energy:.5f}"
            geom = group_dir / f"conf_{rank:03d}.xyz"
            geom.write_text(xyz_text(mol, idx, title), encoding="utf-8")
            rows.append(dict(group=name, pubchem_cid=cid, conformer=rank,
                             smiles=smiles, charge=1, n_atoms=mol.GetNumAtoms(),
                             mmff_kcal_mol=energy,
                             mmff_relative_kcal_mol=energy-base,
                             **hbond_features(mol, idx),
                             xyz_path=str(geom.relative_to(PRIVATE)).replace("\\", "/")))
        groups.append(dict(group=name, cid=cid, smiles=smiles,
                           n_embedded=len(confids), n_optimized=len(viable),
                           n_selected=len(selected),
                           best_MMFF94s_kcal_mol=base))
    with (PRIVATE / "conformers.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    meta = dict(title="Position-isomer pilot from frozen Study 1 shortlist",
                initial_method="RDKit ETKDGv3 + MMFF94s",
                rdkit_version=rdBase.rdkitVersion, seed=SEED,
                n_requested_embeddings=args.nembed,
                heavy_atom_rms_cutoff_A=RMS_CUTOFF,
                energy_window_kcal_mol=ENERGY_WINDOW_KCAL,
                charge=1, total_selected=len(rows),
                shortlist_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
                groups=groups,
                status="MMFF_CONFORMERS_ONLY; no xTB/DFT yet")
    (RESULTS / "preparation.json").write_text(json.dumps(meta, indent=2) + "\n",
                                             encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
