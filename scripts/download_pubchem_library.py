"""Retrieve a reproducible, query-defined PubChem amine candidate pool."""

from __future__ import annotations

import argparse
import csv
import io
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/CanonicalSMILES,MolecularFormula,MolecularWeight,IUPACName/CSV"
QUERIES = (
    "amine", "primary amine", "benzylamine", "phenethylamine", "aniline",
    "pyridine", "diamine", "amino alcohol", "fluoro amine", "bromo amine",
)


def get(url: str, attempts: int = 5) -> bytes:
    last: OSError | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "perovskite-molecular-prioritization/1.0"},
            )
            with urllib.request.urlopen(request, timeout=90) as response:
                return response.read()
        except OSError as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Download failed after {attempts} attempts: {url}") from last


def search_ids(query: str, retmax: int) -> list[str]:
    term = f"{query}[All Fields] AND 70:350[MW]"
    url = ESEARCH + "?" + urllib.parse.urlencode(
        {"db": "pccompound", "term": term, "retmode": "json", "retmax": retmax}
    )
    payload = json.loads(get(url))
    return list(payload["esearchresult"]["idlist"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("data/external/pubchem_amine_pool.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("data/external/pubchem_amine_pool.manifest.json"))
    parser.add_argument("--target", type=int, default=15000)
    parser.add_argument("--per-query", type=int, default=2500)
    args = parser.parse_args()
    ids: list[str] = []
    seen: set[str] = set()
    query_counts: dict[str, int] = {}
    for query in QUERIES:
        found = search_ids(query, args.per_query)
        query_counts[query] = len(found)
        for cid in found:
            if cid not in seen:
                ids.append(cid)
                seen.add(cid)
            if len(ids) >= args.target:
                break
        if len(ids) >= args.target:
            break
        time.sleep(0.34)
    rows: list[dict[str, str]] = []
    for start in range(0, len(ids), 200):
        batch = ids[start:start + 200]
        body = get(PUG.format(cids=",".join(batch))).decode("utf-8", errors="replace")
        rows.extend(csv.DictReader(io.StringIO(body)))
        time.sleep(0.34)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["CID", "ConnectivitySMILES", "MolecularFormula", "MolecularWeight", "IUPACName"])
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "source": "PubChem Compound via NCBI E-utilities and PUG REST",
        "retrieval_date": "2026-08-22",
        "queries": list(QUERIES),
        "query_counts": query_counts,
        "requested_unique_cids": len(ids),
        "retrieved_rows": len(rows),
        "selection_note": "Query-defined pool, not an exhaustive or random sample of PubChem.",
        "usage_policy": "Requests throttled below five per second.",
    }
    args.manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
