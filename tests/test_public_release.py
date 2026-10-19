"""Public-release checks that do not depend on nonredistributed input data."""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from verify_results import verify


def test_output_hashes_match_manifest():
    count, fingerprint = verify()
    assert count >= 20
    assert len(fingerprint) == 64


def test_reported_cohort_and_ranker_results():
    cohort = pd.read_json(ROOT / "results/clean_study1/device_cohort.json", typ="series")
    assert int(cohort["eligible_device_rows"]) == 749
    assert int(cohort["all42_n"]) == 42
    assert int(cohort["mono31_n"]) == 31
    results = pd.read_csv(ROOT / "results/clean_study1/benchmark.csv")
    assert {"logistic2", "random_order"}.issubset(results.model.unique())
    assert {"nested", "family_ood"}.issubset(results.protocol.unique())
