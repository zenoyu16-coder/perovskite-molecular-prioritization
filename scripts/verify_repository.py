"""Check public documentation against frozen artifacts; never fit or calculate energies."""
from __future__ import annotations

import json
from pathlib import Path
import re
import tomllib
from urllib.parse import unquote, urlsplit

import numpy as np
import pandas as pd
from PIL import Image
from sklearn.metrics import average_precision_score
import yaml

ROOT = Path(__file__).resolve().parents[1]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def verify() -> dict:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    study = ROOT / "results/clean_study1"
    cohort = json.loads((study / "device_cohort.json").read_text())
    for key in ("source_device_rows", "eligible_device_rows", "n_DOI", "all42_n", "mono31_n"):
        require(re.search(rf"\b{cohort[key]}\b", readme) is not None, f"README cohort missing: {key}")

    benchmark = pd.read_csv(study / "benchmark.csv", float_precision="round_trip")
    rows = [line.split("|")[1:-1] for line in readme.splitlines() if line.startswith("| ")]
    for label, model in (("Random ranking", "random_order"),
                         ("Logistic regression, 2 descriptors", "logistic2"),
                         ("Logistic regression, 35 descriptors", "logistic35"),
                         ("XGBoost, 35 descriptors", "xgboost35")):
        row = next((v for v in rows if v[0].strip() == label), None)
        require(row is not None, f"Missing README model: {label}")
        for column, protocol in enumerate(("nested", "family_ood"), 1):
            value = benchmark.loc[(benchmark.cohort == "mono31") & (benchmark.protocol == protocol)
                                  & (benchmark.model == model), "AP_mean"].item()
            require(row[column].strip() == f"{value:.3f}", f"README AP mismatch: {model}/{protocol}")

    # Recover complete held-out panels; do not average metrics across tiny folds.
    predictions = pd.read_csv(study / "heldout_predictions.csv", float_precision="round_trip")
    panels = []
    for (cohort_name, protocol, model, repeat, draw), z in predictions.groupby(
            ["cohort", "protocol", "model", "repeat", "draw"]):
        expected = 31 if cohort_name == "mono31" else 42
        require(len(z) == z.spacer_id.nunique() == expected, "Incomplete molecular panel")
        score = average_precision_score(z.y_true, z.score) if z.y_true.sum() else 0.0
        panels.append((cohort_name, protocol, model, float(score)))
    ap = pd.DataFrame(panels, columns=["cohort", "protocol", "model", "AP"])
    recomputed = ap.groupby(["cohort", "protocol", "model"]).AP.mean()
    for row in benchmark.itertuples(index=False):
        require(abs(recomputed.loc[(row.cohort, row.protocol, row.model)] - row.AP_mean) < 1e-12,
                f"Stored-prediction AP mismatch: {row.cohort}/{row.protocol}/{row.model}")

    funnel = json.loads((study / "screening_counts.json").read_text())
    for key in ("source_raw_query", "eligible_ions", "monocation_pool", "score_domain_sa_intersection", "pareto"):
        require(f"{funnel[key]:,}" in readme, f"README screening count missing: {key}")
    shortlist = pd.read_csv(study / "shortlist.csv")
    require(len(shortlist) == funnel["final"] == 10 and shortlist.pubchem_cid.is_unique, "Shortlist identity/count")
    require(shortlist.charge.eq(1).all() and shortlist.applicability_class.eq("moderate extrapolation").all(), "Shortlist cohort/domain")
    require(shortlist.final_rank.tolist() == list(range(1, 11)), "Shortlist order")
    require(not shortlist.model_score_is_calibrated_probability.any(), "Score semantics")
    bounds = shortlist.max_tanimoto_similarity.agg(["min", "max"])
    require(f"{bounds['min']:.3f}–{bounds['max']:.3f}" in readme, "README similarity interval")

    q = ROOT / "results/quantum_pilot"
    ion = ROOT / "results/iodide_pilot"
    cat = pd.read_csv(q / "dft_singlepoints.csv").set_index("group")
    pair = pd.read_csv(ion / "dft_pair_singlepoints.csv").set_index("group")
    conversion = 627.509474
    cation_delta = (cat.loc["sulfonyl_ortho", "energy_hartree"] - cat.loc["sulfonyl_meta", "energy_hartree"]) * conversion
    association_delta = ((pair.loc["sulfonyl_ortho", "energy_hartree"] - pair.loc["sulfonyl_meta", "energy_hartree"]) * conversion) - cation_delta
    xtb = json.loads((ion / "association_differences.json").read_text())
    xtb_delta = next(r["ortho_minus_meta_association_kcal_mol"] for r in xtb["energy_differences"] if r["medium"] == "dmf")
    fluoro = json.loads((ion / "fluoro_transfer.json").read_text())["a_minus_b_electronic_association_kcal_mol"]
    require(cation_delta < 0 and association_delta > 0 and xtb_delta > 0 and fluoro > 0, "Energy-direction wording needs review")
    for value in (-cation_delta, association_delta, xtb_delta, fluoro):
        require(f"{value:.2f}" in readme, "README electronic-energy number")
    crest = pd.read_csv(q / "crest_conformers.csv")
    for group, count, contacts in (("sulfonyl_meta", 11, 0), ("sulfonyl_ortho", 4, 4)):
        z = crest[(crest.group == group) & (crest.relative_kcal_mol <= 3)]
        require(len(z) == count and int(z.hb_strict.sum()) == contacts, "CREST contact/count mismatch")

    docs = [ROOT / "README.md", ROOT / "data/README.md", *sorted((ROOT / "docs").glob("*.md"))]
    links = 0
    images = set()
    for doc in docs:
        text = doc.read_text(encoding="utf-8")
        for match in re.finditer(r"(!?)\[[^\]]*\]\(([^\s)]+)\)", text):
            is_image, target = match.groups()
            target = target.strip("<>")
            parts = urlsplit(target)
            if parts.scheme or parts.netloc:
                continue  # Network availability is separate from public offline checks.
            path = (doc.parent / unquote(parts.path)).resolve() if parts.path else doc
            require(path.exists(), f"Broken link: {doc.relative_to(ROOT)} -> {target}")
            if parts.fragment:
                heading_text = path.read_text(encoding="utf-8")
                anchors = {re.sub(r"[^\w\- ]", "", h.lower()).replace(" ", "-")
                           for h in re.findall(r"^#+\s+(.+)$", heading_text, flags=re.M)}
                require(parts.fragment in anchors, f"Unknown Markdown anchor: {target}")
            links += 1
            if is_image:
                images.add(path)
                with Image.open(path) as im:
                    require(im.width > 100 and im.height > 100, f"Image too small: {path.name}")
                    im.verify()

    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    cite = yaml.safe_load((ROOT / "CITATION.cff").read_text())
    require(cite["cff-version"] == "1.2.0" and cite["type"] == "software", "CFF metadata")
    require(cite["title"] == readme.splitlines()[0].removeprefix("# "), "Citation title")
    require(cite["version"] == cfg["project"]["version"], "Citation/package version")
    require(cite["repository-code"] == cfg["project"]["urls"]["Repository"], "Repository URL metadata")
    require(cfg["project"]["optional-dependencies"]["test"] == ["pytest==9.0.2"], "Test extra")
    require("MIT License" in (ROOT / "LICENSE").read_text(), "Code license")
    return dict(status="verified", local_links=links, readable_linked_images=len(images),
                AP_panels_recomputed=len(panels), source_inputs_required=False,
                model_fits=0, quantum_jobs=0)


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
