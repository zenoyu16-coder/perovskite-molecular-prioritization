"""Clean Study 1 invariants and stale-result rejection."""
from __future__ import annotations
import copy
import json
from pathlib import Path

import pandas as pd
import pytest

from perovskite_prioritization import clean_study1 as clean
from perovskite_prioritization.chemistry import DESCRIPTOR_COLUMNS

ROOT=Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(not (ROOT / "data/source_v2/paper1_X_author_faithful.csv").exists(), reason="Requires nonredistributed source matrices")
OUTPUT=ROOT/"results/clean_study1"


def test_raw_matrix_validation_and_ion_distribution():
    dev,spacers,cols=clean.load_sources(ROOT)
    assert len(dev)==749 and len(spacers)==42
    assert spacers.ionic_charge.value_counts().to_dict()=={1:31,2:11}
    assert dev.Ref_DOI_number.nunique()==134
    assert not any(x.startswith("large_organic_") for x in cols)
    assert (spacers[[*DESCRIPTOR_COLUMNS]].notna()).all().all()


def test_clean_build_does_not_depend_on_any_historical_result():
    code=(ROOT/"src/perovskite_prioritization/clean_study1.py").read_text(encoding="utf-8")
    for forbidden in ("results/final/","results/study1_current/","data/processed/",
                      "model_oof_predictions.csv","results/chemical_space_atlas.csv"):
        assert forbidden not in code


def test_complete_clean_output_and_one_reproducible_entry():
    record=clean.run(ROOT,check=True)
    assert record["provenance"]["config"]["chemistry_rule"]=="amine_only_v2"
    assert record["provenance"]["config"]["descriptor_order"]==list(DESCRIPTOR_COLUMNS)
    cohort=json.loads((OUTPUT/"device_cohort.json").read_text())
    assert cohort["eligible_device_rows"]==749
    assert cohort["all42_n"]==42 and cohort["mono31_n"]==31
    table=pd.read_csv(OUTPUT/"benchmark.csv")
    assert set(table.cohort)=={"all42","mono31"}
    assert set(table.protocol)=={"nested","family_ood"}
    assert set(table.model)=={"xgboost35","logistic35","logistic2","constant","random_order"}
    assert (OUTPUT/"manuscript.pdf").is_file()
    assert (OUTPUT/"manuscript.pdf").stat().st_size>10000


def test_full_and_mono_heldout_coverage():
    predictions=pd.read_csv(OUTPUT/"heldout_predictions.csv")
    for cohort,number in (("all42",42),("mono31",31)):
        sub=predictions[predictions.cohort.eq(cohort)]
        for protocol,rep_count in (("nested",10),("family_ood",1)):
            for rep in range(rep_count):
                for name in ("xgboost35","logistic35","logistic2","constant"):
                    p=sub[(sub.protocol==protocol)&(sub.repeat==rep)&(sub.model==name)]
                    assert len(p)==p.spacer_id.nunique()==number
    shortlist=pd.read_csv(OUTPUT/"shortlist.csv")
    assert len(shortlist)==10 and shortlist.charge.eq(1).all()
    assert shortlist.pubchem_cid.is_unique


@pytest.mark.parametrize("changed",["chemistry_rule","descriptor_order","models","max_sa",
                                         "candidate_pool_scope","outer_nested_repeats"])
def test_fingerprint_detects_change(monkeypatch,changed):
    state=clean.source_state(ROOT)
    original=copy.deepcopy(clean.CONFIG[changed])
    change=original+[999] if isinstance(original,list) else str(original)+"-changed"
    monkeypatch.setitem(clean.CONFIG,changed,change)
    mutated=clean.source_state(ROOT)
    assert state["fingerprint"]!=mutated["fingerprint"]
    with pytest.raises(ValueError,match="STALE_OUTPUT"):
        clean.run(ROOT,check=True)


def test_tampered_file_fails_manifest(tmp_path,monkeypatch):
    real=json.loads((OUTPUT/"manifest.json").read_text())
    dst=tmp_path/"results/clean_study1"
    dst.mkdir(parents=True)
    (dst/"manifest.json").write_text(json.dumps(real))
    (dst/"benchmark.csv").write_text("tampered")
    monkeypatch.setattr(clean,"source_state",lambda path:real["provenance"])
    with pytest.raises(ValueError,match="STALE_OUTPUT"):
        clean.run(tmp_path,check=True)
