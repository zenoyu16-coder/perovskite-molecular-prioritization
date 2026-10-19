"""Source-disjoint, label-stability and family diagnostics for Clean Study 1."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

from .chemistry import DESCRIPTOR_COLUMNS


def positive_publication_stability(device,mono,cols,context_targets,labels):
    """Recompute the *entire* mono-cohort context when one positive's DOI is omitted."""
    names=mono.spacer_id.tolist()
    original=context_targets(device,cols,names)
    before,cut=labels(original)
    positive=[name for name,y in zip(names,before) if y]
    if len(positive)!=7:
        raise ValueError("Expected seven full-cohort monocationic positive identities")
    original_map=dict(zip(names,before))
    baseline=dict(zip(names,original))
    rows=[]
    doi=device.Ref_DOI_number.astype(str)
    for name in positive:
        publications=sorted(device.loc[device.large_organic.eq(name),"Ref_DOI_number"]
                            .astype(str).unique())
        for pub in publications:
            available=device[doi.ne(pub)]
            supported=[x for x in names if x in set(available.large_organic)]
            if name not in supported:
                rows.append({"spacer_id":name,"removed_DOI":pub,"original_label":1,
                   "n_original_doi":len(publications),"remaining_n_doi":0,
                   "status":"no_remaining_source","label_after_drop":pd.NA,
                   "median_uplift_original":baseline[name],
                   "median_uplift_after_drop":np.nan,"label_cutoff_after_drop":np.nan,
                   "n_spacers_with_support":len(supported)})
                continue
            try:
                uplift=context_targets(available,cols,supported)
                y,threshold=labels(uplift)
                mapping=dict(zip(supported,range(len(supported))))
                present=np.asarray(available.loc[available.large_organic.eq(name),
                                                  "Ref_DOI_number"].astype(str).unique())
                rows.append({"spacer_id":name,"removed_DOI":pub,"original_label":1,
                   "n_original_doi":len(publications),"remaining_n_doi":len(present),
                   "status":"recomputed","label_after_drop":int(y[mapping[name]]),
                   "median_uplift_original":baseline[name],
                   "median_uplift_after_drop":float(uplift[mapping[name]]),
                   "label_cutoff_after_drop":float(threshold),
                   "n_spacers_with_support":len(supported)})
            except (ValueError,IndexError) as exc:
                rows.append({"spacer_id":name,"removed_DOI":pub,"original_label":1,
                   "n_original_doi":len(publications),"remaining_n_doi":np.nan,
                   "status":"recompute_failed:"+type(exc).__name__,
                   "label_after_drop":pd.NA,"median_uplift_original":baseline[name],
                   "median_uplift_after_drop":np.nan,"label_cutoff_after_drop":np.nan,
                   "n_spacers_with_support":len(supported)})
    return pd.DataFrame(rows)


def family_diagnostics(predictions,labels_df,cohort="mono31"):
    """Use actual outer-fold test labels and scores, not full-data positive labels."""
    m=predictions[(predictions.cohort==cohort)&(predictions.protocol=="family_ood")
                   &(predictions.model=="logistic2")].copy()
    lookup=labels_df[labels_df.cohort==cohort][
        ["spacer_id","molecular_family","n_devices","n_doi","label_full",
         "median_context_uplift"]]
    m=m.merge(lookup,on="spacer_id",validate="one_to_one")
    if not m.spacer_id.is_unique:raise AssertionError("Family OOD requires unique heldout spacer")
    m["within_family_score_rank"]=m.groupby("molecular_family")["score"].rank(
        ascending=False,method="min").astype(int)
    family=[]
    for name,group in m.groupby("molecular_family",sort=True):
        pos=int(group.y_true.sum());neg=len(group)-pos
        mixed=pos>0 and neg>0
        family.append({"molecular_family":name,"n_spacers":len(group),
             "heldout_positive":pos,"heldout_negative":neg,
             "within_family_AP":float(average_precision_score(group.y_true,group.score))
                  if mixed else np.nan,
             "ranking_identifiable":bool(mixed),
             "mean_score":float(group.score.mean()),
             "min_score":float(group.score.min()),
             "max_score":float(group.score.max()),
             "n_single_DOI":int(group.n_doi.eq(1).sum()),
             "n_label_changed_vs_full":int(group.y_true.ne(group.label_full).sum())})
    return m.sort_values(["molecular_family","within_family_score_rank"]),pd.DataFrame(family)


def doi_support_strata(predictions,labels_df,cohort="mono31"):
    """Small-group audits; undefined single-class AP is recorded as missing."""
    lookup=labels_df[labels_df.cohort==cohort][["spacer_id","n_doi"]]
    m=predictions[(predictions.cohort==cohort)&(predictions.draw.eq(-1))].merge(
        lookup,on="spacer_id",validate="many_to_one")
    m["source_support"]=np.where(m.n_doi.eq(1),"one_DOI","two_plus_DOI")
    rows=[]
    for (protocol,model,source),sub in m.groupby(
            ["protocol","model","source_support"]):
        reps=range(10) if protocol=="nested" else range(1)
        for repeat in reps:
            group=sub[sub.repeat.eq(repeat)]
            mixed=group.y_true.nunique()==2
            rows.append({"protocol":protocol,"model":model,"source_support":source,
                "repeat":repeat,"n_spacers":len(group),"n_positive":int(group.y_true.sum()),
                "AP_when_both_classes":float(average_precision_score(
                    group.y_true,group.score)) if mixed else np.nan,
                "class_mixed":bool(mixed)})
    return pd.DataFrame(rows)


def source_disjoint_sensitivity(device,mono,cols,predictions,
                                context_targets,labels,context_model,model,cols_for):
    """Remove every source DOI belonging to outer-test spacers from outer training.

    This is a first-repetition and family-wise sensitivity analysis, not an
    independent prospective estimate. If source exclusion eliminates a training
    species, that species is dropped from that fold's training and counted.
    """
    mono=mono.reset_index(drop=True)
    ids=mono.spacer_id.to_numpy(str)
    X=mono[list(DESCRIPTOR_COLUMNS)].to_numpy(float)
    fixed=predictions[(predictions.cohort=="mono31")&
                      (predictions.model=="logistic2")&
                      (predictions.draw.eq(-1))&
                      (predictions.repeat.eq(0))]
    records=[];folds=[]
    all_doi=device.Ref_DOI_number.astype(str)
    for (protocol,fold),original in fixed.groupby(["protocol","fold"]):
        held=original.spacer_id.tolist()
        train=[s for s in ids if s not in set(held)]
        test_doi=set(device.loc[device.large_organic.isin(held),
                               "Ref_DOI_number"].astype(str))
        eligible=device[device.large_organic.isin(train)]
        source_doi=set(eligible.Ref_DOI_number.astype(str))
        clean_train=eligible.loc[~eligible.Ref_DOI_number.astype(str).isin(test_doi)]
        if set(clean_train.Ref_DOI_number.astype(str)) & test_doi:
            raise AssertionError("Training/test publication source overlap")
        kept=sorted(set(clean_train.large_organic))
        statistic={"protocol":protocol,"fold":int(fold),
                   "test_spacers":len(held),"test_DOIs":len(test_doi),
                   "overlapping_DOIs":len(test_doi&source_doi),
                   "source_train_rows_before":len(eligible),
                   "source_train_rows_after":len(clean_train),
                   "chemistry_train_before":len(train),
                   "chemistry_train_after":len(kept),
                   "status":"evaluated"}
        if len(kept)<8 or clean_train.Ref_DOI_number.nunique()<5:
            statistic["status"]="insufficient_training_support"
            folds.append(statistic);continue
        try:
            residuals=context_targets(clean_train,cols,kept)
            ytr,threshold=labels(residuals)
            if np.unique(ytr).size<2:
                statistic["status"]="training_single_class"
                folds.append(statistic);continue
            mapping={x:i for i,x in enumerate(ids)}
            ix=np.asarray([mapping[x] for x in kept],int)
            it=np.asarray([mapping[x] for x in held],int)
            clf=model("logistic2").fit(X[ix][:,cols_for("logistic2")],ytr)
            pred=clf.predict_proba(X[it][:,cols_for("logistic2")])[:,1]
            held_dev=device[device.large_organic.isin(held)]
            fitted=context_model().fit(
                clean_train[cols].to_numpy(float),
                clean_train.JV_reverse_scan_PCE.to_numpy(float))
            lifted=held_dev.assign(_uplift=held_dev.JV_reverse_scan_PCE.to_numpy(float)-
                fitted.predict(held_dev[cols].to_numpy(float)))
            target=lifted.groupby("large_organic")["_uplift"].median().reindex(held)
            if target.isna().any():raise AssertionError("Unscored test chemistry")
            yt=(target.to_numpy(float)>=threshold).astype(int)
            for name,t,p in zip(held,yt,pred,strict=True):
                records.append({"protocol":protocol,"fold":int(fold),
                 "spacer_id":name,"y_true_source_disjoint":int(t),
                 "score_source_disjoint":float(p),
                 "y_true_original_protocol":int(original.set_index("spacer_id").loc[name,"y_true"]),
                 "score_original_protocol":float(original.set_index("spacer_id").loc[name,"score"])})
            statistic["target_cutoff"]=float(threshold)
            statistic["heldout_positives"]=int(yt.sum())
            statistic["test_label_changes"]=int(
                (yt!=original.set_index("spacer_id").reindex(held).y_true.to_numpy()).sum())
        except (ValueError,AssertionError,IndexError) as exc:
            statistic["status"]="failed:"+type(exc).__name__
        folds.append(statistic)
    score=pd.DataFrame(records)
    report=pd.DataFrame(folds)
    return score,report


def disjoint_summary(scores,folds):
    "Two complete diagnostic panels (one nested repeat, eight families)."
    output=[]
    for protocol in ("nested","family_ood"):
        panel=scores[scores.protocol.eq(protocol)].copy()
        base=panel.y_true_original_protocol.astype(int).to_numpy()
        test=panel.y_true_source_disjoint.astype(int).to_numpy()
        original=panel.score_original_protocol.to_numpy(float)
        source=panel.score_source_disjoint.to_numpy(float)
        n_expected=31
        n_scored=panel.spacer_id.nunique()
        if len(panel)!=n_scored:raise AssertionError("Repeated heldout chemistry")
        ap=lambda y,p:float(average_precision_score(y,p)) if (len(np.unique(y))==2) else np.nan
        output.append({"protocol":protocol,"n_heldout_spacers":n_scored,
         "n_expected_spacers":n_expected,"complete_panel":bool(n_scored==n_expected),
         "n_source_disjoint_folds":int(folds[folds.protocol.eq(protocol)].status.eq("evaluated").sum()),
         "n_positive_original_labels":int(base.sum()),
         "n_positive_source_disjoint_labels":int(test.sum()),
         "n_test_label_changed":int((base!=test).sum()),
         "original_protocol_AP":ap(base,original),
         "source_disjoint_AP":ap(test,source),
         "interpretation":"sensitivity_relabel_and_retrain_not_prospective_performance"})
    return pd.DataFrame(output)


def annotate_candidate_priority(monocation,shortlist):
    mono=monocation.copy();ranking=shortlist.copy()
    mono["model_score_percentile_rank"]=mono.score.rank(pct=True,method="average")
    mono["model_score_descending_rank"]=mono.score.rank(ascending=False,method="min").astype(int)
    mono["candidate_priority_descending_rank"]=mono.candidate_priority.rank(
        ascending=False,method="min").astype(int)
    ranking=ranking.merge(
        mono[["pubchem_cid","model_score_percentile_rank","model_score_descending_rank",
              "candidate_priority_descending_rank"]],
        on="pubchem_cid",validate="one_to_one")
    ranking["model_score_is_calibrated_probability"]=False
    ranking["rank_semantics"]="diversity_selection_order"
    ranking["score_semantics"]="uncalibrated_bootstrap_mean_classifier_output"
    ranking["applicability_semantics"]="heuristic_no_error_rate_calibration"
    return mono,ranking
