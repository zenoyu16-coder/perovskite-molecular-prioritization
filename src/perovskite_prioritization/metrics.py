"""Ranking metrics and simple multi-objective selection helpers."""

from __future__ import annotations

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


def ranking_metrics(y_true: np.ndarray, score: np.ndarray) -> dict[str, float]:
    y = np.asarray(y_true, dtype=int)
    s = np.asarray(score, dtype=float)
    prevalence = float(y.mean())
    order = np.argsort(-s)
    result: dict[str, float] = {
        "auroc": float(roc_auc_score(y, s)) if np.unique(y).size == 2 else np.nan,
        "average_precision": float(average_precision_score(y, s)) if y.sum() else np.nan,
        "prevalence": prevalence,
    }
    for fraction in (0.10, 0.20):
        k = max(1, int(np.ceil(len(y) * fraction)))
        hits = int(y[order[:k]].sum())
        precision = hits / k
        recall = hits / max(1, int(y.sum()))
        result[f"precision_at_{int(fraction * 100)}pct"] = float(precision)
        result[f"recall_at_{int(fraction * 100)}pct"] = float(recall)
        result[f"ef_at_{int(fraction * 100)}pct"] = float(precision / prevalence) if prevalence else np.nan
    result["top_decile_recovery"] = result["recall_at_10pct"]
    return result


def percentile_scale(values: np.ndarray, higher_is_better: bool = True) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    order = np.argsort(np.argsort(arr, kind="mergesort"), kind="mergesort")
    scaled = order / max(1, len(arr) - 1)
    return scaled if higher_is_better else 1.0 - scaled


def pareto_mask(objectives: np.ndarray) -> np.ndarray:
    """Return nondominated rows for objectives already oriented to maximize."""
    values = np.asarray(objectives, dtype=float)
    keep = np.ones(len(values), dtype=bool)
    for i, row in enumerate(values):
        if not keep[i]:
            continue
        dominated = np.all(values >= row, axis=1) & np.any(values > row, axis=1)
        if dominated.any():
            keep[i] = False
    return keep
