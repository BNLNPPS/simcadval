#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Grade a submitted diagnosis against the known answer of an eval case.

    grade.py CASE_DIR --answer shift10um

Reads CASE_DIR/cmp/diagnosis.json (written by svtool submit_diagnosis) and prints the checks
and the score as JSON. Works for any harness: Claude Code, ollama_agent.py, or a person.
Score = passed checks minus one per false finding (another category with share > 0.01).
See docs/agent-eval.md for how eval cases are made and what the scores so far were.
"""
import argparse
import json
import os
import re

# Known answers, keyed by the fault (the eval copy itself gets a neutral name).
ANSWERS = {
    # scripts/selftest.sh shift10um on examples/MiniSiD.xml: the G4OCCT navigation artifact
    # (eic/G4OCCT#405, thin ECal slices, ~100 % of the mismatch) plus the inner vertex modules
    # component0@4 and component1@5 moved by D = (0, 0, +0.010) mm.
    "shift10um": {
        "nav_parts": r"slice\d@",
        "moved": {"component0@4", "component1@5", "component0@*", "component1@*"},
        "moved_names": {"component0", "component1"},
        "D_mm": (0.0, 0.0, 0.010),
        "tol_mm": (0.003, 0.003, 0.002),
    },
}
EXPECTED = ("navigation_artifact", "placement")


def grade(case, answer="shift10um"):
    key = ANSWERS[answer]
    p = os.path.join(case, "cmp", "diagnosis.json")
    if not os.path.isfile(p):
        return {"submitted": False, "score": None}
    with open(p) as f:
        fs = json.load(f)["findings"]
    nav = [f for f in fs if f["category"] == "navigation_artifact"]
    pl = [f for f in fs if f["category"] == "placement"]
    parts_pl = [q for f in pl for q in f["parts"]]
    g = {"submitted": True,
         "nav_artifact_found": bool(nav) and sum(f["share_of_mismatch"] for f in nav) >= 0.9,
         "nav_parts_right": bool(nav) and all(re.match(key["nav_parts"], q) for f in nav for q in f["parts"]),
         "placement_found": bool(pl),
         "placement_parts_right": bool(pl) and all(q in key["moved"] for q in parts_pl),
         "placement_parts_complete": bool(pl) and {q.split("@")[0] for q in parts_pl} >= key["moved_names"]}
    D = next((f["magnitude"]["displacement_mm"] for f in pl
              if (f.get("magnitude") or {}).get("displacement_mm")), None)
    g["displacement_mm"] = D
    g["displacement_right"] = bool(D) and all(abs(d - e) <= t for d, e, t in zip(D, key["D_mm"], key["tol_mm"]))
    g["false_findings"] = [f["category"] for f in fs
                           if f["category"] not in EXPECTED and f["share_of_mismatch"] > 0.01]
    checks = ("nav_artifact_found", "nav_parts_right", "placement_found", "placement_parts_right",
              "placement_parts_complete", "displacement_right")
    g["score"] = sum(bool(g[k]) for k in checks) - len(g["false_findings"])
    g["max_score"] = len(checks)
    return g


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case")
    ap.add_argument("--answer", required=True, choices=sorted(ANSWERS))
    a = ap.parse_args()
    print(json.dumps(grade(a.case, a.answer)))


if __name__ == "__main__":
    main()
