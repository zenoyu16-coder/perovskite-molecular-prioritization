"""Recompute iodide association results from released total energies.

Tests check arithmetic, matching sources, references, calculation parameters
and structural file provenance. They do not assume a desired energy ordering.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "iodide_pilot"
KCAL = 627.509474
sys.path.insert(0, str(ROOT / "scripts"))
from run_iodide_reference import (IODIDE_XYZ, calculate_reference, extract_energy)


def load_json(name):
    return json.loads((OUT / name).read_text(encoding="utf-8"))


def energy_minima(df, energy_key, group_key="group", medium_key="medium"):
    # Actual minimized electronic energy is the selection criterion.
    return df.loc[df.groupby([group_key, medium_key])[energy_key].idxmin()].set_index(
        [group_key, medium_key]
    )


def test_iodide_reference_matches_input_and_raw_record():
    ref = load_json("iodide_reference.json")
    assert hashlib.sha256((OUT / "iodide.xyz").read_bytes()).hexdigest() == ref["input_sha256"]
    assert (OUT / "iodide.xyz").read_text(encoding="utf-8") == IODIDE_XYZ
    data = pd.read_csv(OUT / "optimized_ion_pairs.csv")
    for medium, solvation in (("dmf", "ALPB(DMF)"), ("gas", "gas")):
        ion = ref["results"][medium]
        assert ion["converged"] and ion["charge"] == -1 and ion["gfn"] == 2
        assert ion["solvation"] == solvation
        np.testing.assert_allclose(
            data.loc[data.medium == medium, "iodide_reference_Eh"].to_numpy(),
            ion["energy_hartree"], rtol=0, atol=1e-8
        )


@pytest.mark.parametrize("medium", ["gas", "dmf"])
def test_reference_runner_constructs_correct_xtb_command(tmp_path, monkeypatch, medium):
    import run_iodide_reference as runner
    seen = {}
    def fake_run(command, **kwargs):
        seen.update(command=command, **kwargs)
        return SimpleNamespace(returncode=0, stdout=
            " | TOTAL ENERGY  -4.142471526080 Eh |\nnormal termination of xtb\n",
            stderr="")
    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    result = calculate_reference("xtb-test", tmp_path, medium)
    assert result["converged"] is True and result["charge"] == -1
    assert seen["command"][:7] == ["xtb-test", "iodide.xyz", "--gfn", "2",
                                    "--chrg", "-1", "--sp"]
    assert ("--alpb" in seen["command"]) == (medium == "dmf")
    assert seen["cwd"] == tmp_path / "iodide_association" / "isolated_iodide"
    assert (seen["cwd"] / "iodide.xyz").read_text(encoding="utf-8") == IODIDE_XYZ
    assert result["energy_hartree"] == pytest.approx(-4.14247152608, abs=1e-12)


def test_reference_parser_rejects_incomplete_xtb_run():
    with pytest.raises(RuntimeError, match="normal termination"):
        extract_energy("TOTAL ENERGY -4.14 Eh\nFAILURE")
    with pytest.raises(RuntimeError, match="No total energy"):
        extract_energy("normal termination of xtb")


def test_xtb_row_energies_and_summaries():
    raw = pd.read_csv(OUT / "optimized_ion_pairs.csv")
    assert len(raw) == 84
    assert set(raw.medium) == {"gas", "dmf"}
    assert set(raw.group) == {"sulfonyl_meta", "sulfonyl_ortho"}
    calculated = (raw["pair_energy_hartree"] - raw["cation_reference_Eh"] -
                  raw["iodide_reference_Eh"]) * KCAL
    np.testing.assert_allclose(calculated, raw["electronic_association_kcal_mol"],
                               rtol=0, atol=2e-6)
    minima = energy_minima(raw, "pair_energy_hartree")
    reported = pd.read_csv(OUT / "ion_pair_summary.csv").set_index(["group", "medium"])
    for key, row in minima.iterrows():
        record = reported.loc[key]
        group, medium = key
        assert record["n_optimized"] == len(raw[(raw.group == group) & (raw.medium == medium)])
        assert record["n_external_Hbond"] == raw[
            (raw.group == group) & (raw.medium == medium)
        ]["external_Hbond_any"].sum()
        assert record["lowest_pair_Eh"] == pytest.approx(row["pair_energy_hartree"], abs=2e-8)
        assert record["association_kcal_mol"] == pytest.approx(
            row["electronic_association_kcal_mol"], abs=2e-6
        )
        assert (OUT / f"{group}_{medium}_lowest.xyz").is_file()

    comparisons = load_json("association_differences.json")
    for rec in comparisons["energy_differences"]:
        medium = rec["medium"]
        a = minima.loc[("sulfonyl_ortho", medium)]
        b = minima.loc[("sulfonyl_meta", medium)]
        delta_pair = (a["pair_energy_hartree"] - b["pair_energy_hartree"]) * KCAL
        delta_cat = (a["cation_reference_Eh"] - b["cation_reference_Eh"]) * KCAL
        delta_assoc = delta_pair - delta_cat
        assert rec["ortho_minus_meta_pair_total_kcal_mol"] == pytest.approx(
            delta_pair, abs=2e-6)
        assert rec["ortho_minus_meta_monomer_total_kcal_mol"] == pytest.approx(
            delta_cat, abs=2e-6)
        assert rec["ortho_minus_meta_association_kcal_mol"] == pytest.approx(
            delta_assoc, abs=3e-6)


def test_fluoromethoxy_energy_recomputation():
    raw = pd.read_csv(OUT / "fluoro_ion_pairs.csv")
    assert len(raw) == 40
    assert set(raw.group) == {"fluoromethoxy_A", "fluoromethoxy_B"}
    ref = load_json("iodide_reference.json")["results"]["dmf"]["energy_hartree"]
    expected = (raw["pair_Eh"] - raw["cation_Eh"] - ref) * KCAL
    np.testing.assert_allclose(expected, raw["association_kcal_mol"], rtol=0, atol=2e-6)
    minima = raw.loc[raw.groupby("group")["pair_Eh"].idxmin()].set_index("group")
    reported = pd.read_csv(OUT / "fluoro_ion_pair_summary.csv").set_index("group")
    for group, row in minima.iterrows():
        assert reported.loc[group, "pair_min_Eh"] == pytest.approx(row["pair_Eh"], abs=2e-8)
        assert reported.loc[group, "association_kcal_mol"] == pytest.approx(
            row["association_kcal_mol"], abs=2e-6)
    double_difference = (minima.loc["fluoromethoxy_A", "pair_Eh"] -
                         minima.loc["fluoromethoxy_B", "pair_Eh"] -
                         minima.loc["fluoromethoxy_A", "cation_Eh"] +
                         minima.loc["fluoromethoxy_B", "cation_Eh"]) * KCAL
    declared = load_json("fluoro_transfer.json")
    assert declared["a_minus_b_electronic_association_kcal_mol"] == pytest.approx(
        double_difference, abs=3e-6)
    assert declared["observed_positive"] == bool(double_difference > 0)


def test_dft_raw_energies_reproduce_double_difference():
    # Both raw ion-pair energies must be public and provenance checked.
    pairs = pd.read_csv(OUT / "dft_pair_singlepoints.csv").set_index("group")
    cations = pd.read_csv(ROOT / "results/quantum_pilot/dft_singlepoints.csv").set_index("group")
    keys = {"sulfonyl_meta", "sulfonyl_ortho"}
    assert set(pairs.index) == keys and keys.issubset(set(cations.index))
    for group in keys:
        pair, cat = pairs.loc[group], cations.loc[group]
        assert bool(pair["converged"]) and bool(cat["converged"])
        assert pair["xc"].lower() == cat["xc"].lower() == "pbe0"
        assert pair["basis"].lower() == cat["basis"].lower() == "def2-svp"
        assert pair["iodine_ecp"] == "def2-svp ECP28"
        assert pair["solvation"] == "ddCOSMO(DMF;eps=36.7)"
        assert int(pair["charge"]) == 0 and int(cat["charge"]) == 1
        assert int(pair["n_atoms"]) == int(cat["n_atoms"]) + 1
        assert (ROOT / pair["input_xyz"]).is_file()
    pair_delta = (pairs.loc["sulfonyl_ortho", "energy_hartree"] -
                  pairs.loc["sulfonyl_meta", "energy_hartree"]) * KCAL
    cat_delta = (cations.loc["sulfonyl_ortho", "energy_hartree"] -
                 cations.loc["sulfonyl_meta", "energy_hartree"]) * KCAL
    double_difference = pair_delta - cat_delta
    report = load_json("dft_association_check.json")["pbe0_def2svp_ddcosmo"]
    assert report["delta_pair_kcal_mol"] == pytest.approx(pair_delta, abs=2e-6)
    assert report["delta_cation_kcal_mol"] == pytest.approx(cat_delta, abs=2e-6)
    assert report["delta_association_kcal_mol"] == pytest.approx(double_difference, abs=3e-6)


def test_coordinate_metadata_and_counts():
    one = pd.read_csv(OUT / "sulfonyl_initial_poses.csv")
    two = pd.read_csv(OUT / "fluoro_initial_poses.csv")
    assert len(one) == 42 and len(two) == 40
    for group, n in (("sulfonyl_meta", 25), ("sulfonyl_ortho", 25),
                     ("fluoromethoxy_A", 23), ("fluoromethoxy_B", 23)):
        xyz = OUT / f"{group}_dmf_lowest.xyz"
        lines = xyz.read_text(encoding="utf-8").splitlines()
        assert int(lines[0]) == n and len(lines) >= n + 2
        assert lines[n + 1].split()[0] == "I"
