"""Acceptance tests for publication dependence and post-selection accounting."""
from pathlib import Path
import json

import pandas as pd

from perovskite_prioritization.clean_study1 import CONFIG,choose_model
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/clean_study1"


def test_frozen_choice_does_not_reuse_outer_test_for_model_selection():
    # Reversing the AP comparison must not silently change the frozen ranker.
    fabricated=pd.DataFrame([
       {"cohort":"mono31","protocol":protocol,"model":model,"AP_mean":score}
       for protocol in ("nested","family_ood")
       for model,score in (
           ("logistic2",0.01),("logistic35",0.99),
           ("xgboost35",0.99),("random_order",0.25))
    ])
    assert choose_model(fabricated)[0]=="logistic2"
    assert CONFIG["selection_bias_status"]=="outer_scores_previously_consulted_not_selection_unbiased"
    decision=json.loads((OUT/"model_selection.json").read_text())
    assert decision["selected_model"]=="logistic2"
    assert decision["selection_test_independent"] is False
    assert decision["unbiased_generalization_estimate_available"] is False
    assert decision["prospective_validation_required"] is True


def test_single_source_positive_status_is_not_imputed():
    stability=pd.read_csv(OUT/"positive_source_stability.csv")
    assert stability.spacer_id.nunique()==7 and len(stability)==13
    singleton=stability[stability.n_original_doi.eq(1)]
    assert len(singleton)==5
    assert singleton.status.eq("no_remaining_source").all()
    assert singleton.label_after_drop.isna().all()
    assert stability.query("spacer_id == 'GU'").label_after_drop.sum()==1
    assert stability.query("spacer_id == 'ThMA'").label_after_drop.sum()==1


def test_family_table_reveals_imbalanced_holds_and_score_errors():
    detail=pd.read_csv(OUT/"family_spacer_diagnostics.csv")
    family=pd.read_csv(OUT/"family_level_diagnostics.csv")
    assert detail.spacer_id.nunique()==31
    assert len(family)==8
    assert int(family.ranking_identifiable.sum())==2
    functional=family[family.molecular_family.eq("functionalized_aliphatic_ammonium")].iloc[0]
    assert functional.n_spacers==7 and functional.heldout_positive==5
    examples=detail[detail.molecular_family.eq("functionalized_aliphatic_ammonium")].sort_values(
        "within_family_score_rank")
    assert examples.iloc[0].spacer_id=="HEA" and examples.iloc[0].y_true==0


def test_doi_disjoint_sensitivity_stores_all_folds_and_targets():
    f=pd.read_csv(OUT/"source_disjoint_folds.csv")
    s=pd.read_csv(OUT/"source_disjoint_predictions.csv")
    m=pd.read_csv(OUT/"source_disjoint_summary.csv")
    assert len(f)==13 and f.status.eq("evaluated").all()
    assert f.overlapping_DOIs.ge(0).all()
    assert f.source_train_rows_after.le(f.source_train_rows_before).all()
    assert set(m.protocol)=={"nested","family_ood"}
    assert m.complete_panel.all() and m.n_heldout_spacers.eq(31).all()
    for _,group in s.groupby("protocol"):
        assert group.spacer_id.nunique()==31
    assert f[f.protocol.eq("family_ood")].shape[0]==8
    assert f[f.protocol.eq("nested")].shape[0]==5


def test_score_priority_and_final_order_are_distinct():
    final=pd.read_csv(OUT/"shortlist.csv")
    assert len(final)==10
    assert final.applicability_class.eq("moderate extrapolation").all()
    assert final.max_tanimoto_similarity.between(.30,.42).all()
    assert not final.score.is_monotonic_decreasing
    assert final.model_score_is_calibrated_probability.eq(False).all()
    assert final.rank_semantics.eq("diversity_selection_order").all()
    assert set(final.model_score_descending_rank)!=set(final.final_rank)
    assert final.score_semantics.str.contains("uncalibrated").all()
    assert final.set_index("pubchem_cid").loc[16460979,"score"] > final.set_index(
        "pubchem_cid").loc[2760937,"score"]
