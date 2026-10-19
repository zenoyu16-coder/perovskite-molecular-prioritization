"""Run or validate a complete Clean Study 1 research build."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from perovskite_prioritization.clean_study1 import run

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--rebuild",action="store_true")
    group.add_argument("--check",action="store_true")
    options=parser.parse_args()
    m=run(ROOT,rebuild=options.rebuild,check=options.check)
    print(json.dumps({"status":"verified","fingerprint":m["provenance"]["fingerprint"],
                      "source_commit":m["source_git_commit"],
                      "files":len(m["output_hashes"])},indent=2))
