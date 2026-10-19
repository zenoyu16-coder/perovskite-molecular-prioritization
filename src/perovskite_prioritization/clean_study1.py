"""Isolated, source-locked Study 1 analysis; never reads historical results.

All public tables, figures, candidate scores and manuscript are built in one
temporary output directory, then atomically installed under results/clean_study1.
Only explicit --rebuild may replace that complete directory. --check verifies
source, chemistry, feature-order, config, software and output fingerprints.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdkit
import sklearn
import xgboost
import yaml
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from sklearn.base import clone
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from .chemistry import DESCRIPTOR_COLUMNS, descriptors, filter_candidate, fingerprint, tanimoto_matrix
from .metrics import pareto_mask, percentile_scale
from .protonation import PROTONATION_POLICY_V2

RAW_FILES = (
    "data/source_v2/paper1_X_author_faithful.csv",
    "data/source_v2/paper1_y.csv",
    "data/source_v2/paper1_metadata.csv",
    "data/source_v2/spacer_master.csv",
    "data/source_v2/spacer_families.csv",
    "data/source_v2/feature_groups.yaml",
    "data/source_v2/paper1_matrix_manifest.json",
    "data/source_v2/molecular_representation_manifest.json",
    "data/external/pubchem_amine_pool.csv",
)
CODE_FILES = (
    "src/perovskite_prioritization/clean_study1.py",
    "src/perovskite_prioritization/chemistry.py",
    "src/perovskite_prioritization/protonation.py",
    "src/perovskite_prioritization/metrics.py",
    "scripts/run_clean_study1.py",
    "src/perovskite_prioritization/study1_diagnostics.py",
)
CONFIG = {
    "schema": "clean-study1-2",
    "chemistry_rule": PROTONATION_POLICY_V2,
    "training_ion_policy": "source_verified_charged_species",
    "candidate_pool_scope": "monocations_only",
    "device_filter": "one_verified_spacer_no_comma",
    "context_estimator": "StandardScaler+Ridge(alpha=10)",
    "doi_context_folds": 5,
    "outer_nested_folds": 5,
    "outer_nested_repeats": 10,
    "family_ood": "leave_one_family_out",
    "target": "training_fold_p80_median_out_of_DOI_context_residual",
    "candidate_target": "cohort_full_DOI_OOF_median_residual_p80",
    "random_repeats": 30,
    "ensemble_bootstraps": 30,
    "models": ["xgboost35", "logistic35", "logistic2", "constant", "random_order"],
    "random_seed": 20260822,
    "selection_mode": "exploratory_frozen_logistic2_after_inspecting_outer_scores",
    "frozen_screening_model": "logistic2",
    "selection_bias_status": "outer_scores_previously_consulted_not_selection_unbiased",
    "source_diagnostic_protocol": "mono_positive_one_DOI_at_a_time_and_source_disjoint_one_nested_repeat_plus_family",
    "score_gate_percentile": .80,
    "similarity_interpolation": .60,
    "similarity_moderate": .30,
    "domain_descriptor_factor": 2.5,
    "max_sa": 5.,
    "pareto_objectives": ["score_rank","uncertainty_rank","sa_rank","ad_rank","novelty_rank"],
    "diversity_similarity": .55,
    "n_preliminary": 20,
    "n_final": 10,
    "morgan_radius": 2,
    "morgan_bits": 256,
    "descriptor_order": list(DESCRIPTOR_COLUMNS),
}

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source_state(root: Path) -> dict:
    missing = [p for p in (*RAW_FILES, *CODE_FILES) if not (root/p).is_file()]
    if missing: raise FileNotFoundError("Clean Study 1 sources missing: " + repr(missing))
    files = {p:digest(root/p) for p in (*RAW_FILES,*CODE_FILES)}
    versions = {"rdkit":rdkit.__version__,"sklearn":sklearn.__version__,
                "xgboost":xgboost.__version__,
                "pandas":pd.__version__,"numpy":np.__version__}
    data={"files":files,"config":CONFIG,"software":versions}
    data["fingerprint"]=hashlib.sha256(json.dumps(data,sort_keys=True,
                     separators=(",",":")).encode()).hexdigest()
    return data

def load_sources(root:Path):
    src=root/"data/source_v2"
    checks=json.loads((src/"paper1_matrix_manifest.json").read_text())["outputs"]
    for filename,key in (("paper1_X_author_faithful.csv","author_faithful_X"),
                         ("paper1_y.csv","y"),("paper1_metadata.csv","metadata")):
        if digest(src/filename)!=checks[key]["sha256"]:
            raise ValueError("Published matrix provenance digest mismatch: "+filename)
    meta=pd.read_csv(src/"paper1_metadata.csv",keep_default_na=False)
    X=pd.read_csv(src/"paper1_X_author_faithful.csv")
    y=pd.read_csv(src/"paper1_y.csv")
    master=pd.read_csv(src/"spacer_master.csv")
    fam=pd.read_csv(src/"spacer_families.csv")
    reg=yaml.safe_load((src/"feature_groups.yaml").read_text())
    groups=reg["feature_groups"]
    cols=[*groups["composition"]["columns"],*groups["device_architecture"]["columns"],
          *groups["processing"]["columns"]]
    if len(cols)!=len(set(cols)) or any(c not in X for c in cols):
        raise ValueError("Invalid context feature registry")
    if len(meta)!=len(X) or len(X)!=len(y) or len(X)!=792:raise ValueError("Source row alignment")
    if not master.spacer_id.is_unique or not fam.spacer_id.is_unique:
        raise ValueError("Source molecular identifier duplication")
    master=master.loc[master.structure_verified.astype(bool)].copy()
    mask=(~meta.large_organic.astype(str).str.contains(",",regex=False) &
          meta.large_organic.isin(master.spacer_id))
    if mask.sum()!=749:raise ValueError("Expected 749 source-selected devices")
    indices=np.flatnonzero(mask.to_numpy())
    device=pd.concat([meta.iloc[indices].reset_index(drop=True),
        y.iloc[indices][["JV_reverse_scan_PCE"]].reset_index(drop=True),
        X.iloc[indices][cols].reset_index(drop=True)],axis=1)
    if len(device)!=749 or device[cols+["JV_reverse_scan_PCE"]].isna().any().any():
        raise ValueError("Missing source device data")
    used=sorted(device.large_organic.unique())
    chem=master[master.spacer_id.isin(used)][["spacer_id","SMILES"]].merge(
        fam[["spacer_id","molecular_family"]],on="spacer_id",validate="one_to_one"
    ).sort_values("spacer_id").reset_index(drop=True)
    if len(chem)!=42:raise ValueError("Expected 42 source-verified chemical identities")
    smiles=[];charges=[]
    for row in chem.itertuples(index=False):
        mol=Chem.MolFromSmiles(row.SMILES)
        if mol is None:raise ValueError("Unparseable original cation "+row.spacer_id)
        n=int(Chem.GetFormalCharge(mol))
        if n not in (1,2):raise ValueError("Unexpected cation charge: "+row.spacer_id)
        smiles.append(Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True))
        charges.append(n)
    chem["modeled_cation_smiles"]=smiles
    chem["ionic_charge"]=charges
    for index,k in enumerate(DESCRIPTOR_COLUMNS):
        chem[k]=[descriptors(s)[k] for s in smiles]
    count=device.groupby("large_organic").agg(n_devices=("JV_reverse_scan_PCE","size"),
                   n_doi=("Ref_DOI_number","nunique"))
    chem=chem.merge(count,left_on="spacer_id",right_index=True,validate="one_to_one")
    if len(chem.query("ionic_charge == 1"))!=31:raise ValueError("Expected 31 monocations")
    return device,chem,cols

def context():
    return make_pipeline(StandardScaler(),Ridge(alpha=10.0))

def labels(values):
    cutoff=float(np.quantile(np.asarray(values,float),.80,method="higher"))
    return (np.asarray(values,float)>=cutoff).astype(int),cutoff

def uplift_by_doi(device:pd.DataFrame,cols:list[str],ids:list[str],*,heldout=None):
    selected=device[device.large_organic.isin(ids)].copy()
    X=selected[cols].to_numpy(float); y=selected.JV_reverse_scan_PCE.to_numpy(float)
    doi=selected.Ref_DOI_number.astype(str).to_numpy()
    if len(np.unique(doi))<5:raise ValueError("Insufficient DOI for context cross-fitting")
    pred=np.full(len(y),np.nan)
    for tr,te in GroupKFold(n_splits=5).split(X,y,doi):
        pred[te]=clone(context()).fit(X[tr],y[tr]).predict(X[te])
    if not np.isfinite(pred).all():raise AssertionError("Context coverage")
    uplift=selected.assign(_uplift=y-pred).groupby("large_organic")["_uplift"].median().reindex(ids)
    if uplift.isna().any():raise ValueError("Missing chemistry residual")
    if heldout is None:return uplift.to_numpy(float)
    check=device[device.large_organic.isin(heldout)].copy()
    expected=context().fit(X,y).predict(check[cols].to_numpy(float))
    held=check.assign(_uplift=check.JV_reverse_scan_PCE.to_numpy(float)-expected).groupby(
        "large_organic")["_uplift"].median().reindex(heldout)
    if held.isna().any():raise ValueError("Missing heldout residual")
    return uplift.to_numpy(float),held.to_numpy(float)

def model(name):
    if name=="xgboost35":
        return XGBClassifier(n_estimators=140,max_depth=2,learning_rate=.05,
            subsample=.85,colsample_bytree=.8,reg_lambda=2.0,
            eval_metric="logloss",random_state=CONFIG["random_seed"],n_jobs=1)
    if name in ("logistic35","logistic2"):
        return make_pipeline(StandardScaler(),LogisticRegression(C=1.,
          class_weight="balanced",random_state=CONFIG["random_seed"],max_iter=1000))
    if name=="constant":return DummyClassifier(strategy="prior")
    raise KeyError(name)

def cols_for(name):
    if name=="logistic2":
        return [list(DESCRIPTOR_COLUMNS).index("MolLogP"),
                list(DESCRIPTOR_COLUMNS).index("CarbonHeavyAtomFraction")]
    return list(range(len(DESCRIPTOR_COLUMNS)))

def metric(g):
    y=g.y_true.to_numpy(int); s=g.score.to_numpy(float)
    if len(y)!=len(set(g.spacer_id)):raise AssertionError("Chemistry-level duplication")
    prevalence=float(y.mean())
    k=max(1,math.ceil(len(y)*.20))
    ordered=np.argsort(-s,kind="stable")
    precision=float(y[ordered[:k]].mean())
    return {"AP":float(average_precision_score(y,s)) if y.sum() else 0.,
      "AUROC":float(roc_auc_score(y,s)) if len(np.unique(y))==2 else np.nan,
      "EF20":precision/prevalence if prevalence else 0.,"prevalence":prevalence}

def holdout_scores(cohort_name, chem,device,ctx_cols):
    ids=chem.spacer_id.to_numpy(str)
    X=chem[list(DESCRIPTOR_COLUMNS)].to_numpy(float)
    full_uplift=uplift_by_doi(device,ctx_cols,list(ids))
    original_y,_=labels(full_uplift)
    groups=chem.molecular_family.astype(str).to_numpy()
    splits=list(RepeatedStratifiedKFold(n_splits=5,n_repeats=10,random_state=20260842).split(X,original_y))
    splits.extend((np.flatnonzero(groups!=fam),np.flatnonzero(groups==fam))
                  for fam in sorted(set(groups)))
    records=[];folds=[]
    for ix,(tr,te) in enumerate(splits):
        protocol="nested" if ix<50 else "family_ood"
        repeat=ix//5 if ix<50 else 0
        fold=ix%5 if ix<50 else ix-50
        train_u,test_u=uplift_by_doi(device,ctx_cols,list(ids[tr]),heldout=list(ids[te]))
        ytr,cut=labels(train_u); yte=(test_u>=cut).astype(int)
        if len(np.unique(ytr))<2:
            raise ValueError(f"One-class outer training cohort {cohort_name} {protocol} {fold}")
        folds.append({"cohort":cohort_name,"protocol":protocol,"repeat":repeat,"fold":fold,
          "n_train":len(tr),"n_test":len(te),"train_pos":int(ytr.sum()),
          "test_pos":int(yte.sum()),"cutoff":cut,
          "changed_test_labels":int((original_y[te]!=yte).sum()),
          "excluded_family":str(groups[te[0]]) if protocol=="family_ood" else ""})
        for name in ("xgboost35","logistic35","logistic2","constant"):
            cols_idx=cols_for(name)
            fitted=model(name).fit(X[tr][:,cols_idx],ytr)
            pred=fitted.predict_proba(X[te][:,cols_idx])[:,1]
            for a,val in zip(range(len(te)),pred,strict=True):
                records.append({"cohort":cohort_name,"protocol":protocol,
                  "repeat":repeat,"fold":fold,"model":name,"draw":-1,
                  "spacer_id":ids[te[a]],"y_true":int(yte[a]),"score":float(val)})
        for draw in range(CONFIG["random_repeats"]):
            rnd=np.random.default_rng(20260822+ix*1009+draw)
            for a,val in enumerate(rnd.random(len(te))):
                records.append({"cohort":cohort_name,"protocol":protocol,
                  "repeat":repeat,"fold":fold,"model":"random_order","draw":draw,
                  "spacer_id":ids[te[a]],"y_true":int(yte[a]),"score":float(val)})
    return pd.DataFrame(records),pd.DataFrame(folds),pd.DataFrame({
        "spacer_id":ids,"cohort":cohort_name,"modeled_cation_smiles":chem.modeled_cation_smiles,
        "molecular_family":groups,"n_devices":chem.n_devices,
        "n_doi":chem.n_doi,"median_context_uplift":full_uplift,
        "label_full":original_y,
    })

def summarize(predictions):
    rows=[]
    for (cohort,protocol,name),part in predictions.groupby(["cohort","protocol","model"]):
        metrics=[]
        repeats=range(10) if protocol=="nested" else range(1)
        draws=range(CONFIG["random_repeats"]) if name=="random_order" else [-1]
        for rep in repeats:
            for draw in draws:
                g=part[(part.repeat==rep)&(part.draw==draw)]
                expected=42 if cohort=="all42" else 31
                if len(g)!=expected or g.spacer_id.nunique()!=expected:
                    raise AssertionError(f"Incomplete heldout coverage {cohort} {protocol} {name}")
                metrics.append(metric(g.rename(columns={"y_true":"y_true"})))
        rows.append({"cohort":cohort,"protocol":protocol,"model":name,
          "n_scoring_runs":len(metrics),"n_spacers":expected,
          **{f"{k}_mean":float(np.nanmean([m[k] for m in metrics]))
             for k in ("AP","AUROC","EF20","prevalence")},
          "AP_p10":float(np.quantile([m["AP"] for m in metrics],.10)),
          "AP_p90":float(np.quantile([m["AP"] for m in metrics],.90))})
    return pd.DataFrame(rows).sort_values(["cohort","protocol","model"])

def choose_model(scores):
    """Fixed exploratory choice, explicitly conditional on previously inspected OOD results.

    The AP comparisons remain model diagnostics, not an independent post-choice
    generalization estimate. Prospective devices are the next independent test.
    """
    a=scores[scores.cohort.eq("mono31")]
    valid={}
    for name in ("xgboost35","logistic35","logistic2"):
        z=a[a.model.eq(name)].set_index("protocol").AP_mean
        rand=a[a.model.eq("random_order")].set_index("protocol").AP_mean
        valid[name]=min(float(z.loc[x]-rand.loc[x]) for x in ("nested","family_ood"))
    selected=CONFIG["frozen_screening_model"]
    if selected not in valid:raise ValueError("Frozen model not in documented comparison")
    return selected,valid

def candidate_pool(root):
    raw=pd.read_csv(root/"data/external/pubchem_amine_pool.csv")
    params=FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    cat=FilterCatalog(params)
    RDLogger.DisableLog("rdApp.warning")
    RDLogger.DisableLog("rdApp.error")
    records=[]; seen=set(); count={}
    for row in raw.itertuples(index=False):
        res=filter_candidate(row.ConnectivitySMILES,cat,
                             protonation_policy=PROTONATION_POLICY_V2)
        count[res.reason]=count.get(res.reason,0)+1
        if not res.accepted:continue
        if res.canonical_smiles in seen:
            count["duplicate_precursor"]=count.get("duplicate_precursor",0)+1
            continue
        seen.add(res.canonical_smiles)
        mol=Chem.MolFromSmiles(res.modeled_cation_smiles)
        charge=int(Chem.GetFormalCharge(mol))
        records.append({"pubchem_cid":int(row.CID),"precursor_smiles":res.canonical_smiles,
         "modeled_cation_smiles":res.modeled_cation_smiles,"charge":charge,
         "sa_score":res.sa_score})
    p=pd.DataFrame(records)
    if len(p)<100 or p.pubchem_cid.duplicated().any():
        raise ValueError("Invalid candidate pool")
    return p,{"source_pubchem_rows":len(raw),"eligible_ions":len(p),
              "monocations":int((p.charge==1).sum()),
              "multications":int((p.charge!=1).sum()),
              "filter_rejection_counts":count}

def fit_candidates(cohort,selection,pool):
    train=cohort.copy()
    X=train[list(DESCRIPTOR_COLUMNS)].to_numpy(float)
    y,_=labels(train.median_context_uplift.to_numpy(float))
    mono=pool.loc[pool.charge==1].copy().reset_index(drop=True)
    xnew=np.asarray([[descriptors(s)[c] for c in DESCRIPTOR_COLUMNS]
                   for s in mono.modeled_cation_smiles],float)
    boot=[];rng=np.random.default_rng(20260822)
    ipos=np.flatnonzero(y==1); ineg=np.flatnonzero(y==0)
    for _ in range(CONFIG["ensemble_bootstraps"]):
        indices=np.r_[rng.choice(ipos,len(ipos),replace=True),
                      rng.choice(ineg,len(ineg),replace=True)]
        ix=cols_for(selection)
        fitted=model(selection).fit(X[indices][:,ix],y[indices])
        boot.append(fitted.predict_proba(xnew[:,ix])[:,1])
    arr=np.array(boot);mono["score"]=arr.mean(axis=0)
    mono["ensemble_sd"]=arr.std(axis=0,ddof=1)
    names=["xgboost35","logistic35","logistic2"]
    for name in names:
        m=model(name).fit(X[:,cols_for(name)],y)
        mono[name+"_score"]=m.predict_proba(xnew[:,cols_for(name)])[:,1]
    trainfp=np.vstack([fingerprint(s) for s in train.modeled_cation_smiles])
    candfp=np.vstack([fingerprint(s) for s in mono.modeled_cation_smiles])
    similarity=tanimoto_matrix(candfp,trainfp)
    near=similarity.argmax(axis=1)
    mono["nearest_training_spacer"]=train.iloc[near].spacer_id.to_numpy()
    mono["max_tanimoto_similarity"]=similarity[np.arange(len(mono)),near]
    from scipy.spatial.distance import cdist
    from sklearn.decomposition import PCA
    scaler=StandardScaler().fit(X)
    xt=scaler.transform(X);xc=scaler.transform(xnew)
    dt=cdist(xt,xt);np.fill_diagonal(dt,np.inf)
    dc=cdist(xc,xt).min(axis=1)
    pca=PCA(n_components=min(16,len(train)-1),random_state=20260822).fit(trainfp)
    at=pca.transform(trainfp);ac=pca.transform(candfp)
    em=cdist(at,at);np.fill_diagonal(em,np.inf)
    edm=cdist(ac,at).min(axis=1)
    dscale=max(float(np.quantile(dt.min(axis=1),.9)),1e-6)
    escale=max(float(np.quantile(em.min(axis=1),.9)),1e-6)
    sim=mono.max_tanimoto_similarity.to_numpy()
    mono["ad_score"]=.45*sim+.35*np.exp(-dc/dscale)+.20*np.exp(-edm/escale)
    mono["applicability_class"]=np.where(
      (sim>=.6)&(dc<=dscale),"interpolation-like",
      np.where((sim>=.3)&(dc<=2.5*dscale),"moderate extrapolation","strong extrapolation"))
    mono["score_rank"]=percentile_scale(mono.score.to_numpy())
    mono["uncertainty_rank"]=percentile_scale(mono.ensemble_sd.to_numpy(),False)
    mono["sa_rank"]=percentile_scale(np.clip((10-mono.sa_score.to_numpy())/9,0,1))
    mono["ad_rank"]=percentile_scale(mono.ad_score.to_numpy())
    mono["novelty_rank"]=percentile_scale(1-sim)
    mono["candidate_priority"]=(.40*mono.score_rank+.15*mono.uncertainty_rank+
       .15*mono.sa_rank+.15*mono.ad_rank+.15*mono.novelty_rank)
    top=mono[mono.score_rank>=.8]
    in_domain=mono[mono.applicability_class!="strong extrapolation"]
    eligible=mono[(mono.score_rank>=.8)&(mono.applicability_class!="strong extrapolation")&
                  (mono.sa_score<=5)].copy()
    objectives=eligible[list(CONFIG["pareto_objectives"])].to_numpy(float)
    eligible["pareto_optimal"]=pareto_mask(objectives)
    kept=eligible[eligible.pareto_optimal].copy()
    if len(kept)<20:
        kept=pd.concat([kept,eligible.loc[~eligible.pareto_optimal].nlargest(20-len(kept),"candidate_priority")])
    top=kept.nlargest(min(1000,len(kept)),"candidate_priority")
    subfp=np.vstack([fingerprint(s) for s in top.modeled_cation_smiles])
    sims=tanimoto_matrix(subfp,subfp)
    assigned=set(range(len(top)));clusters=np.full(len(top),-1,int)
    for center in range(len(top)):
        if center not in assigned:continue
        members={j for j in assigned if sims[center,j]>=CONFIG["diversity_similarity"]}
        for j in members:clusters[j]=center
        assigned-=members
    top=top.copy();top["diversity_cluster"]=clusters
    prelim=top.sort_values("candidate_priority",ascending=False).drop_duplicates(
        "diversity_cluster").head(20).copy()
    if len(prelim)<10:raise ValueError("Diversity pool cannot supply ten")
    ix=list(prelim.index)
    pref=np.vstack([fingerprint(s) for s in prelim.modeled_cation_smiles])
    pref_sim=tanimoto_matrix(pref,pref)
    chosen=[0]
    while len(chosen)<10:
        others=[i for i in range(len(prelim)) if i not in chosen]
        dist=np.asarray([1-pref_sim[i,chosen].max() for i in others])
        pri=prelim.iloc[others].candidate_priority.to_numpy()
        chosen.append(others[int(np.argmax(.65*pri+.35*dist))])
    shortlist=prelim.iloc[chosen].copy()
    shortlist.insert(0,"final_rank",np.arange(1,len(shortlist)+1))
    columns=["final_rank","pubchem_cid","precursor_smiles","modeled_cation_smiles",
             "charge","score","ensemble_sd","max_tanimoto_similarity",
             "nearest_training_spacer","ad_score","applicability_class","sa_score",
             "candidate_priority","xgboost35_score","logistic35_score","logistic2_score"]
    sizes={"source_raw_query":None,"eligible_ions":len(pool),"monocation_pool":len(mono),
           "top_score_quintile":int((mono.score_rank>=.8).sum()),
           "non_strong_domain":len(in_domain),"score_domain_sa_intersection":len(eligible),
           "pareto":int(eligible.pareto_optimal.sum()),"diversity_preliminary":len(prelim),
           "final":len(shortlist)}
    return mono,shortlist[columns],sizes

def make_figures(folder,summary,shortlist,selected):
    fig,ax=plt.subplots(figsize=(8,4.5))
    chart=summary[summary.cohort=="mono31"]
    order=["random_order","constant","logistic2","logistic35","xgboost35"]
    nm={"random_order":"Random ordering","constant":"Constant prevalence",
        "logistic2":"Logistic (2)","logistic35":"Logistic (35)",
        "xgboost35":"XGBoost (35)"}
    yy=np.arange(len(order))
    for k,protocol in enumerate(("nested","family_ood")):
        s=chart[chart.protocol==protocol].set_index("model").reindex(order)
        ax.scatter(s.AP_mean,yy+(.11 if k==0 else -.11),
                   label=("Nested" if k==0 else "Leave-one-family-out"),s=50)
    ax.set_yticks(yy,[nm[i] for i in order]);ax.invert_yaxis()
    ax.set_xlim(0,1.);ax.grid(axis="x",alpha=.2)
    ax.set_xlabel("Average precision, monovalent spacer holdouts")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(folder/"evaluation.png",dpi=160)
    plt.close(fig)
    from rdkit.Chem import Draw
    ms=[Chem.MolFromSmiles(s) for s in shortlist.modeled_cation_smiles]
    legend=[f"CID {int(row.pubchem_cid)} | {row.score:.3f}"
            for row in shortlist.itertuples(index=False)]
    img=Draw.MolsToGridImage(ms,molsPerRow=5,subImgSize=(260,230),legends=legend)
    img.save(str(folder/"candidate_structures.png"))

def build(root:Path,target:Path):
    device,chemical,ctx=load_sources(root)
    chemical.to_csv(target/"ionic_training_descriptors.csv",index=False)
    device_summary={"source_device_rows":792,"eligible_device_rows":len(device),
                    "n_DOI":int(device.Ref_DOI_number.nunique()),
                    "all42_n":len(chemical),
                    "mono31_n":int((chemical.ionic_charge==1).sum()),
                    "di_n":int((chemical.ionic_charge==2).sum()),
                    "context_columns":ctx}
    (target/"device_cohort.json").write_text(json.dumps(device_summary,indent=2))
    allpred=[];allfolds=[];allsumm=[]
    cohorts={"all42":chemical,"mono31":chemical[chemical.ionic_charge==1].copy().reset_index(drop=True)}
    for name,sub in cohorts.items():
        pred,folds,summary=holdout_scores(name,sub,device,ctx)
        allpred.append(pred);allfolds.append(folds);allsumm.append(summary)
    predictions=pd.concat(allpred,ignore_index=True)
    folds=pd.concat(allfolds,ignore_index=True)
    labels_df=pd.concat(allsumm,ignore_index=True)
    scores=summarize(predictions)
    selected,margins=choose_model(scores)
    from .study1_diagnostics import (
        positive_publication_stability,family_diagnostics,doi_support_strata,
        source_disjoint_sensitivity,disjoint_summary,annotate_candidate_priority)
    positive_doi=positive_publication_stability(
        device,cohorts["mono31"],ctx,uplift_by_doi,labels)
    family_species,family_summary=family_diagnostics(predictions,labels_df)
    doi_strata=doi_support_strata(predictions,labels_df)
    disjoint_predictions,disjoint_folds=source_disjoint_sensitivity(
        device,cohorts["mono31"],ctx,predictions,uplift_by_doi,labels,
        context,model,cols_for)
    disjoint_metrics=disjoint_summary(disjoint_predictions,disjoint_folds)
    pool,counts=candidate_pool(root)
    model_training=cohorts["mono31"].copy()
    full_u=uplift_by_doi(device,ctx,list(model_training.spacer_id))
    model_training["median_context_uplift"]=full_u
    monocation,shortlist,funnel=fit_candidates(model_training,selected,pool)
    monocation,shortlist=annotate_candidate_priority(monocation,shortlist)
    funnel["source_raw_query"]=counts["source_pubchem_rows"]
    predictions.to_csv(target/"heldout_predictions.csv",index=False)
    folds.to_csv(target/"outer_folds.csv",index=False)
    labels_df.to_csv(target/"training_labels.csv",index=False)
    scores.to_csv(target/"benchmark.csv",index=False)
    monocation.to_csv(target/"candidate_scores.csv",index=False)
    shortlist.to_csv(target/"shortlist.csv",index=False)
    positive_doi.to_csv(target/"positive_source_stability.csv",index=False)
    family_species.to_csv(target/"family_spacer_diagnostics.csv",index=False)
    family_summary.to_csv(target/"family_level_diagnostics.csv",index=False)
    doi_strata.to_csv(target/"doi_support_strata.csv",index=False)
    disjoint_predictions.to_csv(target/"source_disjoint_predictions.csv",index=False)
    disjoint_folds.to_csv(target/"source_disjoint_folds.csv",index=False)
    disjoint_metrics.to_csv(target/"source_disjoint_summary.csv",index=False)
    (target/"screening_counts.json").write_text(json.dumps(funnel,indent=2))
    (target/"chemical_filter_counts.json").write_text(json.dumps(counts,indent=2))
    (target/"model_selection.json").write_text(json.dumps({
      "selected_model":selected,
      "selection_rule":CONFIG["selection_mode"],
      "selection_evidence":"exploratory_post_selection_comparison",
      "selection_test_independent":False,
      "unbiased_generalization_estimate_available":False,
      "prospective_validation_required":True,
      "exploratory_AP_minimum_gain_over_same_cohort_random":margins,
      "training_cohort":"mono31",
      "training_n":len(model_training),"train_positive":int(labels(full_u)[0].sum()),
      "source_training_labels":"DOI OOF context uplift within mono31 only",
      "probability_note":"Uncalibrated classifier output; not a probability of experimental success",
      "candidate_rank_note":"final_rank is the order of diversity selection, not a sorted model score",
      "applicability_note":"Empirical descriptor and fingerprint rules, not calibrated prediction error"},indent=2))
    make_figures(target,scores,shortlist,selected)
    create_manuscript(target,scores,shortlist,selected,device_summary,counts,funnel,
                     positive_doi,family_summary,disjoint_folds)

def create_manuscript(output,bench,shortlist,selected,dev,counts,funnel,
                      positive_doi,family_summary,disjoint_folds):
    def ap(cohort,protocol,name):
        row=bench[(bench.cohort==cohort)&(bench.protocol==protocol)&(bench.model==name)]
        return f"{float(row.AP_mean.iloc[0]):.3f}"
    n_positive=positive_doi.spacer_id.nunique()
    n_single=(positive_doi.drop_duplicates("spacer_id").n_original_doi==1).sum()
    n_fold_disjoint=(disjoint_folds.status=="evaluated").sum()
    n_mixed=(family_summary.ranking_identifiable).sum()
    lines=[
        r"\\documentclass[11pt,a4paper]{article}",
        r"\\usepackage[margin=25mm]{geometry}",
        r"\\usepackage{graphicx,booktabs,amsmath,hyperref}",
        r"\\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}",
        r"\\title{Molecular prioritization of layered-perovskite spacers from device literature}",
        r"\\author{Research working manuscript}",
        r"\\date{19 October 2026}",
        r"\\begin{document}\\maketitle",
        r"\\begin{abstract}",
        f"We reconstructed {dev['eligible_device_rows']} reported solar-cell devices across {dev['n_DOI']} publications and 42 source-verified ionic spacer structures (31 monovalent).",
        "A DOI-cross-fitted context model generates chemistry-level PCE residual labels in each outer training fold. ",
        "Fixed models and chance baselines are evaluated with chemistry holdouts and leave-one-family-out testing.",
        f"For monovalent chemistry, the exploratory nested AP of {selected} is {ap('mono31','nested',selected)} versus {ap('mono31','nested','random_order')} for random ranking. Model choice was informed by these same external comparisons, so this is not a selection-independent estimate.",
        f"Out of {counts['source_pubchem_rows']} PubChem records, {counts['monocations']} distinct screened monovalent precursor ions were considered; ten molecules were selected for controlled fabrication.",
        r"\\end{abstract}",
        r"\\section{Introduction}",
        "Organic spacer cations affect the processing and properties of layered halide perovskites. ",
        "Existing literature offers many device records but few chemically independent spacer identities. ",
        "Our aim is to determine whether those records support reliable decisions on which untested cations to examine. ",
        "The analysis treats spacers as the prediction units and records published processing conditions as context.",
        r"\\section{Data and chemical representation}",
        f"The source contains 792 device records; {dev['eligible_device_rows']} satisfy single verified spacer requirements.",
        "Original charged species, including aromatic and diammonium ions, are retained as defined by the source chemical registry.",
        "The 35 ordered RDKit descriptors are generated anew for each ion. ",
        "A separate amine-site transformation supplies monovalent candidate structures from neutral PubChem precursors.",
        r"\\section{Validation design}",
        "The device-context estimator is a standardized ridge regression with alpha 10. DOI-grouped out-of-fold residuals define spacer median uplift.",
        "Within each outer chemistry fold, both the context model and the upper-quintile label threshold are regenerated using training chemistry only.",
        "We report ten repetitions of five-fold stratified chemistry holdout and leave-one-molecular-family-out validation.",
        "The candidate ranker is frozen to two-descriptor logistic regression after exploratory inspection of outer validation scores. The AP values in this paper are exploratory model-comparison statistics; an unbiased post-selection assessment requires genuinely external prospective data.",
        "The candidate ranking models are fixed logistic regression with 2 or 35 descriptors and a fixed XGBoost ranker, alongside constant-score and random-order baselines.",
        r"\\section{Results}",
        r"\\begin{table}[h]\\centering",
        r"\\caption{Pooled average precision (AP) for held-out spacer identities.}",
        r"\\begin{tabular}{lrrrr}\\toprule",
        r"Model & All nested & All family & Mono nested & Mono family\\\\ \\midrule",
    ]
    names=[("Random ordering","random_order"),("Constant score","constant"),
           ("Logistic (2)","logistic2"),("Logistic (35)","logistic35"),
           ("XGBoost (35)","xgboost35")]
    for name,key in names:
        lines.append(f"{name} & {ap('all42','nested',key)} & {ap('all42','family_ood',key)} & {ap('mono31','nested',key)} & {ap('mono31','family_ood',key)}"+chr(92)*4)
    lines += [
       r"\\bottomrule\\end{tabular}\\end{table}",
       r"\\section{Literature-source and family diagnostics}",
       f"Of {n_positive} positive monovalent cations in full-data labels, {n_single} have only one supporting DOI. Removing their only publication removes the target observation, rather than providing an independent validation of positive status.",
       "For supported multisource cations, the full-context model and upper-quintile threshold are recalculated after deleting individual publications. The accompanying source-stability table records label changes and cases in which no evidence remains.",
       f"Only {n_mixed} of {len(family_summary)} monovalent held-out molecular families contain both positive and negative labels under fold-specific context refitting. Scores originate from different family-exclusion fits; pooled AP includes between-family ordering and is not equivalent to within-family ranking quality.",
       f"A publication-disjoint diagnostic removes every publication containing a held-out spacer from the corresponding training devices before context fitting. {n_fold_disjoint} source-disjoint folds were evaluable. This stress test changes available training chemistry and source support; it is reported alongside, rather than substituted for, the original model comparison.",
       r"\\begin{figure}[h]\\centering\\includegraphics[width=.89\\linewidth]{evaluation.png}",
       r"\\caption{Fixed models under monovalent chemistry and family holdout protocols.}\\end{figure}",
       r"\\section{Molecular screening}",
       f"The chemical filter retains {counts['eligible_ions']} unique candidate ions, including {counts['monocations']} monovalent and {counts['multications']} non-monovalent structures.",
       f"The frozen screening model is {selected}. The shared score and domain gates yield {funnel['score_domain_sa_intersection']} candidates, followed by {funnel['pareto']} Pareto-optimal entries and ten diversity-selected molecules.",
       f"All ten shortlist molecules are in the moderate extrapolation category under the current empirical fingerprint and descriptor criteria. This class has not been calibrated to an expected prediction error.",
       "The model score is an uncalibrated mean classifier output, not the probability of successful device fabrication. Candidate priority combines rank, uncertainty, synthetic accessibility, applicability and novelty. Final rank is the order of diversity selection and need not follow descending model score.",
       r"\\begin{figure}[h]\\centering\\includegraphics[width=\\linewidth]{candidate_structures.png}",
       r"\\caption{Ten current monovalent candidate cations. Scores are retrospective model rankings.}\\end{figure}",
       r"\\section{Experimental test}",
       "We propose 12 prespecified new spacer compounds: four high-ranked, four chemistry-matched lower-ranked, and four random eligible choices.",
       "Two known spacers serve as batch references. The primary PCE endpoint is evaluated across independent fabrication batches with the spacer, not device, as the independent ranking unit.",
       "Prospective devices supply the first model-selection-independent comparison of the frozen two-descriptor ranker with random ordering and prespecified chemical baselines.",
       "Diffraction, PL/TRPL and morphology measurements provide secondary mechanistic evidence.",
       r"\\section{Data and computation}",
       "The complete clean analysis is regenerated from source-verified device matrices, chemical registries and the PubChem query pool. ",
       "Code and data-file digests, feature order, package versions and model configuration are recorded in the accompanying manifest.",
       r"\\end{document}",
    ]
    tex="\n".join(lines)+"\n"
    tex=tex.replace(chr(92)*2,chr(92))
    (output/"manuscript.tex").write_text(tex,encoding="utf-8")
    compiler=shutil.which("pdflatex")
    if not compiler:
        raise FileNotFoundError("pdflatex required for complete Clean Study 1 manuscript")
    for _ in range(2):
        p=subprocess.run([compiler,"-interaction=nonstopmode","-halt-on-error","manuscript.tex"],
            cwd=output,capture_output=True,text=True,timeout=90)
        if p.returncode:
            raise RuntimeError("Manuscript LaTeX build failed: "+p.stdout[-3200:])
    for ext in (".log",".aux",".out"):
        (output/("manuscript"+ext)).unlink(missing_ok=True)

def run(root:Path,*,rebuild:bool=False,check:bool=False):
    target=root/"results/clean_study1"
    current=source_state(root)
    if check:
        if not (target/"manifest.json").is_file():raise FileNotFoundError("No clean-study manifest")
        doc=json.loads((target/"manifest.json").read_text())
        if doc["provenance"]["fingerprint"]!=current["fingerprint"]:
            raise ValueError("STALE_OUTPUT: source, chemistry, features, model, config or software changed. Run --rebuild.")
        expected=doc["output_hashes"]
        found={p.relative_to(target).as_posix():digest(p)
               for p in target.rglob("*") if p.is_file() and p.name!="manifest.json"}
        if found!=expected:
            raise ValueError("STALE_OUTPUT: output missing, modified or added. Run --rebuild.")
        return doc
    if target.exists() and not rebuild:
        raise FileExistsError("Clean output exists. Run --check, or --rebuild to regenerate it completely.")
    with tempfile.TemporaryDirectory(prefix=".clean_study1_build_",dir=root/"results") as name:
        staging=Path(name)
        build(root,staging)
        source_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=root,text=True).strip()
        records={p.relative_to(staging).as_posix():digest(p)
                 for p in staging.rglob("*") if p.is_file()}
        record={"provenance":current,"source_git_commit":source_commit,
                "historical_tag":"study1-legacy-2026-10-19",
                "output_hashes":records}
        (staging/"manifest.json").write_text(json.dumps(record,indent=2),encoding="utf-8")
        if target.exists():
            backup=root/"results/.clean_study1_previous"
            if backup.exists():shutil.rmtree(backup)
            target.rename(backup)
            try:shutil.copytree(staging,target)
            except BaseException:
                if target.exists():shutil.rmtree(target)
                backup.rename(target)
                raise
            shutil.rmtree(backup)
        else:
            shutil.copytree(staging,target)
    return run(root,check=True)
