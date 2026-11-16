#!/usr/bin/env python3
"""Markdown table of headline metrics across cases.  summarize.py out/CASE1/cmp out/CASE2/cmp ..."""
import json
import os
import sys


def f(x, fmt):
    return "–" if x is None else format(x, fmt)


rows = []
for d in sys.argv[1:]:
    m = json.load(open(os.path.join(d, "metrics.json")))
    h, r, p = m["headline"], m["rays"], m["points"]
    exc_b = sum(m["b"].get("g4_exceptions", {}).values())
    rows.append([
        os.path.basename(os.path.dirname(d.rstrip("/"))),
        f"{r['n']} {m['a']['ray_mode']}",
        f(1 - h["material_agreement_nonworld_length"], ".2e"),
        f(r["volume_mismatch_fraction_nonworld"], ".2e"),
        f(1 - p["material_agreement"], ".2e") if p["n"] else "–",
        f"{p['n_disagree']}/{p['n']}" if p["n"] else "–",
        f(h["boundary_f1_at_1nm"], ".5f"),
        f(h["boundary_f1_at_1um"], ".5f"),
        f(h["max_boundary_shift_mm"], ".2e"),
        f(r["x_over_x0"]["rel_diff"]["p99"], ".1e"),
        f(h["max_rel_dX0"], ".1e"),
        f(h["sequence_identical_fraction"], ".4f"),
        f(h["cpu_ratio_ray_b_over_a"], ".0f"),
        str(r.get("rays_length_differ", "–")),
        str(len(m["b"].get("killed_events", []))),
        str(len(m.get("hung", {}).get("b", []))),
        str(exc_b),
    ])
hdr = ["case", "rays", "1−mat.agree (len)", "1−vol.agree (len)", "1−point agree", "pts differ", "bnd F1 @1nm",
       "bnd F1 @1µm", "max bnd shift [mm]", "|ΔX/X₀| p99", "|ΔX/X₀| max", "same seq.", "CPU ×", "rays len≠", "killed (B)", "hung (B)", "G4 exc. (B)"]
print("| " + " | ".join(hdr) + " |")
print("|" + "---|" * len(hdr))
for r in rows:
    print("| " + " | ".join(r) + " |")
