"""Verify the published result files without requiring restricted source matrices."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "clean_study1"


def verify():
    manifest = json.loads((OUTPUT / "manifest.json").read_text(encoding="utf-8"))
    errors = []
    for filename, expected in manifest["output_hashes"].items():
        path = OUTPUT / filename
        if not path.is_file():
            errors.append(f"Missing: {filename}")
            continue
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed != expected:
            errors.append(f"SHA-256 mismatch: {filename}")
    if errors:
        raise RuntimeError("Published result verification failed:\n" + "\n".join(errors))
    return len(manifest["output_hashes"]), manifest["provenance"]["fingerprint"]


if __name__ == "__main__":
    count, fingerprint = verify()
    print(f"Verified {count} published result files.")
    print(f"Recorded source/config fingerprint: {fingerprint}")
    print("Note: This checks published artifacts, not the unavailable raw input files.")
