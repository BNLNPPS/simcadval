#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""svtool: model-neutral investigation tools for one simcadval case.

Each tool takes named arguments and prints one JSON object on stdout. The same
tools serve Claude Code (skill investigate-mismatch), an SDK agent or a local
model with OpenAI-style tool calls; `svtool.py spec` prints their definitions.

    svtool.py TOOL --arg value ...         e.g. svtool.py ray --case out/x --evt 12 --smin 100 --smax 110
    svtool.py call TOOL '{"case": ...}'    JSON arguments (or @file.json, or - for stdin)
    svtool.py spec [--format anthropic|openai]
    svtool.py list                         tool names and first description lines

Outside the eic container the command re-runs itself through ../eic.sh.
Every call on a case is logged with an id (c001, c002, ...) to
<case>/investigation/calls.jsonl; a diagnosis cites these ids as evidence.
Library output (DD4hep, Geant4) goes to <case>/investigation/svtool.log, so
stdout carries only the JSON result.
"""
import json
import math
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SCHEMA_PATH = os.path.join(REPO, "schemas", "diagnosis.schema.json")
MAX_RESCAN_RAYS = 20000
SAME_MM = 1e-7  # two step edges closer than this are the same boundary
sys.path.insert(0, HERE)

TOOLS = {}


class ToolError(Exception):
    pass


def tool(name, description, properties, required=("case",), log=True):
    def deco(fn):
        TOOLS[name] = {"fn": fn, "description": description, "log": log,
                       "schema": {"type": "object", "properties": properties, "required": list(required),
                                  "additionalProperties": False}}
        return fn
    return deco


def P(type_, description, **kw):
    return {"type": type_, "description": description, **kw}


CASE = P("string", "case directory made by run_case.sh, e.g. out/minisid_ip, or a rescan under "
                   "out/NAME/investigation/")
NUMS = lambda n, d: {"type": "array", "items": {"type": "number"}, "minItems": n, "maxItems": n, "description": d}  # noqa: E731


# ── output helpers ────────────────────────────────────────────────────────────

EXACT = {"s0_mm", "s_mm", "smin", "smax", "s_range_mm", "first_mismatch_s_mm", "last_mismatch_s_mm"}


def compact(x, key=None):
    """JSON-safe copy: numpy -> python, NaN/inf -> None, floats to 7 significant digits
    (positions along a ray to 1e-9 mm)."""
    if type(x).__module__ == "numpy":
        x = x.tolist() if hasattr(x, "shape") and x.shape else x.item()
    if isinstance(x, bool) or x is None or isinstance(x, (int, str)):
        return x
    if isinstance(x, float):
        if not math.isfinite(x) or abs(x) >= 1e90:
            return None
        return round(x, 9) if key in EXACT else float(f"{x:.7g}")
    if isinstance(x, dict):
        return {k: compact(v, k) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [compact(v, key) for v in x]
    return str(x)


def rel(p):
    if not p:
        return p
    r = os.path.relpath(p)
    return r if not r.startswith("../..") else p


# ── case ──────────────────────────────────────────────────────────────────────

class Case:
    def __init__(self, path):
        d = os.path.normpath(path)
        for w in ("dd4hep", "step"):
            if not os.path.isfile(os.path.join(d, w + ".meta.json")):
                raise ToolError(f"{path}: not a case directory (no {w}.meta.json); make one with scripts/run_case.sh")
        self.dir = d
        self.meta_a = json.load(open(os.path.join(d, "dd4hep.meta.json")))
        self.meta_b = json.load(open(os.path.join(d, "step.meta.json")))
        parts = d.split(os.sep + "investigation" + os.sep)
        self.root = parts[0] if len(parts) > 1 else d
        self.step = self.meta_b.get("step", "")
        mm = self.meta_b.get("matmap")
        if mm is None:  # scans made before the scanner recorded the map
            mm = os.path.join(d, "materials.xml") if os.path.isfile(os.path.join(d, "materials.xml")) \
                else self.step + ".materials.xml"
        self.matmap = mm
        rep = os.path.join(d, "export.json")
        if not os.path.isfile(rep):
            rep = os.path.join(os.path.dirname(self.step), "export.json")
        self.report = rep if os.path.isfile(rep) else None
        self.volumes = json.load(open(rep))["volumes"] if self.report else []
        self.vol = {v["id"]: v for v in self.volumes}
        self.world_material = self.meta_a["world_material"]

    def vname(self, vid, mat=None):
        """volume name without id; -1 is the world (or a volume the scan does not map)."""
        if vid < 0:
            return "<outside>" if vid < -1 else ("<world>" if mat in (None, self.world_material) else "<unmapped>")
        v = self.vol.get(vid)
        return v["name"] if v else "?"

    def pname(self, vid, mat=None):
        n = self.vname(vid, mat)
        return f"{n}@{vid}" if vid >= 0 else n

    def resolve(self, q):
        """'name@id', 'name@*', 'name', an id, or a name substring -> list of volume ids."""
        q = str(q).strip()
        m = re.fullmatch(r"(.+)@(\d+|\*)", q)
        if m:
            name, i = m.groups()
            if i == "*":
                return [v["id"] for v in self.volumes if v["name"] == name]
            v = self.vol.get(int(i))
            return [v["id"]] if v and v["name"] == name else []
        if q.isdigit():
            return [int(q)] if int(q) in self.vol else []
        ids = [v["id"] for v in self.volumes if v["name"] == q]
        return ids or [v["id"] for v in self.volumes if q.lower() in v["name"].lower()]


def load_scans(case):
    """Both scans restricted to the rays present in both (hung rays dropped); original ray numbers kept.
    Returns A, B, material names (index -1 -> '<outside>')."""
    import numpy as np
    from compare import Scan
    mats = {}
    A, B = Scan(os.path.join(case.dir, "dd4hep"), mats), Scan(os.path.join(case.dir, "step"), mats)
    common = np.intersect1d(A.rays["evt"], B.rays["evt"])
    for S in (A, B):
        S.rays = S.rays[np.isin(S.rays["evt"], common)]
        keep = np.isin(S.steps["evt"], common)
        S.steps, S.smat = S.steps[keep], S.smat[keep]
    names = [None] * len(mats)
    for n, m in mats.items():
        names[m["id"]] = n
    return A, B, names + ["<outside>"]


def mismatch_intervals(case, A, B):
    """Mismatching intervals of the exact overlay (evt, s, len, matA, matB, volA, volB), total mismatch length.
    Cached per case in <case>/.svtool_mismatch.npz (keyed by the scans' size and mtime)."""
    import numpy as np
    from compare import overlay
    key = repr([(os.path.getsize(f), os.path.getmtime(f)) for f in
                (os.path.join(case.dir, w + ".steps.bin") for w in ("dd4hep", "step"))])
    cache = os.path.join(case.dir, ".svtool_mismatch.npz")
    if os.path.isfile(cache):
        try:
            z = np.load(cache)
            if str(z["key"]) == key:
                return tuple(z[k] for k in ("evt", "s", "ln", "ma", "mb", "va", "vb")), float(z["total"])
        except Exception:
            pass
    evt, s, ln, ma, mb, va, vb = overlay(A, B, 0, 0)
    mis = ma != mb
    out = tuple(x[mis] for x in (evt, s, ln, ma, mb, va, vb))
    total = float(ln[mis].sum())
    try:
        np.savez(cache, key=key, total=total, **dict(zip(("evt", "s", "ln", "ma", "mb", "va", "vb"), out)))
    except OSError:
        pass
    return out, total


def step_edges(S):
    """(evt, s) of every step start and ray end, sorted."""
    import numpy as np
    st = S.steps
    last = np.r_[st["evt"][1:] != st["evt"][:-1], True] if len(st) else np.zeros(0, bool)
    e = np.r_[st["evt"], st["evt"][last]].astype(np.int64)
    s = np.r_[st["s0"], (st["s0"] + st["len"])[last]]
    o = np.lexsort((s, e))
    return e[o], s[o]


def nearest_in_group(ge, gs, qe, qs):
    """Nearest gs to each qs within the same group (ge, gs sorted by (ge, gs)): index (-1: none), |distance|."""
    import numpy as np
    n = len(qs)
    if len(gs) == 0:
        return np.full(n, -1), np.full(n, np.inf)
    lo, hi = np.searchsorted(ge, qe, "left"), np.searchsorted(ge, qe, "right")
    a, b = lo.copy(), hi.copy()
    while True:  # vectorised bisection inside each group: a = first index with gs >= qs
        act = a < b
        if not act.any():
            break
        mid = (a + b) // 2
        right = act & (gs[np.minimum(mid, len(gs) - 1)] < qs)
        a = np.where(right, mid + 1, a)
        b = np.where(act & ~right, mid, b)
    best, bd = np.full(n, -1), np.full(n, np.inf)
    for c in (a - 1, a):
        ok = (c >= lo) & (c < hi)
        cc = np.clip(c, 0, len(gs) - 1)
        d = np.where(ok, np.abs(gs[cc] - qs), np.inf)
        better = d < bd
        bd, best = np.where(better, d, bd), np.where(better, cc, best)
    return best, bd


def top_parts(case, vids, mats, weights, names, k=5):
    out = {}
    for v, m, w in zip(vids.tolist(), mats.tolist(), weights.tolist()):
        key = case.pname(v, names[m])
        out[key] = out.get(key, 0.0) + w
    return [{"part": p, "mm": w} for p, w in sorted(out.items(), key=lambda x: -x[1])[:k]] if len(vids) else []


# ── geometry (TGeo, reference world) ──────────────────────────────────────────

GEO_CPP = r"""
#include <TGeoManager.h>
#include <TGeoNavigator.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
namespace svgeo {
std::vector<TGeoVolume*> gVols;
std::unordered_map<const TGeoVolume*, int> gVid;
std::vector<std::string> gPath;
std::vector<double> gPl;
// same numbering as src/geo_ids.hh (DFS pre-order over TGeoVolumes, first visit wins)
int Init() {
  gVols.clear(); gVid.clear();
  std::unordered_set<const TGeoVolume*> seen;
  std::vector<TGeoVolume*> st{gGeoManager->GetTopVolume()};
  while (!st.empty()) {
    TGeoVolume* v = st.back(); st.pop_back();
    if (!seen.insert(v).second) continue;
    gVid[v] = (int)gVols.size(); gVols.push_back(v);
    for (int i = v->GetNdaughters() - 1; i >= 0; --i) st.push_back(v->GetNode(i)->GetVolume());
  }
  return (int)gVols.size();
}
int Vid(const TGeoVolume* v) { auto it = gVid.find(v); return it == gVid.end() ? -2 : it->second; }
// Boundary normals: for each point p on a boundary (mm) with ray direction u, find the node just
// after (side 1) or before (side 0) the boundary, and take the normal of that node's own shape at p
// (global frame). dist = the shape's safety at p (mm): about 0 when p is on that shape's surface.
void Normals(int n, const double* p, const double* u, const int* side, double* nrm, double* dist, int* vid) {
  TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
  if (!nav) nav = gGeoManager->AddNavigator();
  gPath.assign(n, "");
  const double eps = 1e-4;  // mm
  for (int i = 0; i < n; ++i) {
    const double* pi = p + 3 * i; const double* ui = u + 3 * i;
    double b[3], gp[3], lp[3], ld[3], ln[3];
    for (int k = 0; k < 3; ++k) { b[k] = (pi[k] + (side[i] ? eps : -eps) * ui[k]) * 0.1; gp[k] = pi[k] * 0.1; }
    TGeoNode* node = nav->FindNode(b[0], b[1], b[2]);
    vid[i] = node ? Vid(node->GetVolume()) : -2;
    gPath[i] = nav->GetPath();
    for (int k = 0; k < 3; ++k) nrm[3 * i + k] = 0.0;
    dist[i] = -1;
    if (!node) continue;
    const TGeoHMatrix* m = nav->GetCurrentMatrix();
    m->MasterToLocal(gp, lp);
    m->MasterToLocalVect(ui, ld);
    TGeoShape* sh = node->GetVolume()->GetShape();
    sh->ComputeNormal(lp, ld, ln);
    m->LocalToMasterVect(ln, nrm + 3 * i);
    dist[i] = sh->Safety(lp, sh->Contains(lp)) * 10.0;
  }
}
// Placements of volume vid: per placement 18 doubles (global AABB centre[3], half[3] in mm,
// translation[3] in mm, rotation[9]) and the path; returns the total number of placements.
int Placements(int vid, int maxn) {
  gPl.clear(); gPath.clear();
  if (vid <= 0 || vid >= (int)gVols.size()) return 0;
  TGeoVolume* target = gVols[vid];
  TGeoIterator it(gGeoManager->GetTopVolume());
  TGeoNode* nd; int total = 0;
  while ((nd = it.Next())) {
    if (nd->GetVolume() != target) continue;
    ++total;
    if ((int)gPath.size() >= maxn) continue;
    TString path; it.GetPath(path); gPath.push_back(path.Data());
    const TGeoMatrix* m = it.GetCurrentMatrix();
    TGeoBBox* bb = (TGeoBBox*)target->GetShape();
    const double* o = bb->GetOrigin(); double h[3] = {bb->GetDX(), bb->GetDY(), bb->GetDZ()};
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (int c = 0; c < 8; ++c) {
      double l[3] = {o[0] + ((c & 1) ? h[0] : -h[0]), o[1] + ((c & 2) ? h[1] : -h[1]), o[2] + ((c & 4) ? h[2] : -h[2])}, g[3];
      m->LocalToMaster(l, g);
      for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], g[k]); hi[k] = std::max(hi[k], g[k]); }
    }
    for (int k = 0; k < 3; ++k) gPl.push_back(5.0 * (lo[k] + hi[k]));
    for (int k = 0; k < 3; ++k) gPl.push_back(5.0 * (hi[k] - lo[k]));
    const double* t = m->GetTranslation(); const double* r = m->GetRotationMatrix();
    for (int k = 0; k < 3; ++k) gPl.push_back(10.0 * t[k]);
    for (int k = 0; k < 9; ++k) gPl.push_back(r[k]);
  }
  return total;
}
std::vector<double> Half(int vid) {
  TGeoBBox* bb = (TGeoBBox*)gVols[vid]->GetShape();
  return {10 * bb->GetDX(), 10 * bb->GetDY(), 10 * bb->GetDZ()};
}
}
"""
_GEO = {}


def geometry(compact):
    if compact not in _GEO:
        import ROOT
        import dd4hep  # noqa: F401
        ROOT.gErrorIgnoreLevel = ROOT.kWarning
        dd4hep.Detector.getInstance().fromXML(compact)
        ROOT.gInterpreter.Declare(GEO_CPP)
        ROOT.svgeo.Init()
        _GEO[compact] = ROOT.svgeo
    return _GEO[compact]


# ── tools ─────────────────────────────────────────────────────────────────────

@tool("case_summary",
      "Headline of a case: ray/point agreement, mismatch length, boundary F1, worst parts, materials whose "
      "IoU < 1, robustness counters (hung/killed rays, Geant4 exceptions), export status (failed parts; a "
      "'perturb' entry means the STEP was moved on purpose by dd4hep2step --perturb) and which further outputs "
      "exist. no_difference_test.passes = true means the worlds agree at floating-point level. Start here.",
      {"case": CASE})
def case_summary(case):
    c = Case(case)
    mp = os.path.join(c.dir, "cmp", "metrics.json")
    if not os.path.isfile(mp):
        raise ToolError(f"{mp} missing: run scripts/compare.py (run_case.sh does)")
    m = json.load(open(mp))
    h, r, p = m["headline"], m["rays"], m["points"]
    ma, mb = m["a"], m["b"]
    agree = r["material_agreement_nonworld"]
    nd = (agree is not None and 1 - agree < 1e-8) and not p["n_disagree"] and h["boundary_f1_at_1nm"] == 1
    iou = {k: v for k, v in r["material_iou"].items() if v["iou"] < 1 - 1e-12}
    worst = [{"part": f"{v['name']}@{v['id']}", "material": c.vol.get(v["id"], {}).get("material"),
              "iou": v["iou"], "mismatch_mm": v["mismatch_mm"]} for v in m["worst_volumes"] if v["mismatch_mm"] > 0][:8]
    exp = None
    if c.report:
        R = json.load(open(c.report))
        exp = {"report": rel(c.report), "volumes": len(R["volumes"]),
               "parts_exported": sum(1 for v in R["volumes"] if v["exported"]),
               "parts_failed": [f"{v['name']}@{v['id']}" for v in R["volumes"]
                                if v["status"] not in ("ok", "world", "world-material", "assembly")
                                or (v["exported"] and not v["readback_ok"])][:20]}
        if R.get("perturb"):
            exp["perturb"] = R["perturb"]
    inv = os.path.join(c.root, "investigation")
    out = {
        "case": rel(c.dir), "rays": r["n"], "ray_mode": ma["ray_mode"], "points": p["n"], "seed": ma["seed"],
        "box_mm": ma["box_mm"], **({"target_mm": ma["target_mm"]} if ma["ray_mode"] == "aim" else {}),
        "world_material": c.world_material,
        "material_agreement_nonworld": agree,
        "mismatch_fraction_nonworld": r["material_mismatch_fraction_nonworld"],
        "mismatch_length_mm": r["material_mismatch_length_mm"], "nonworld_length_mm": r["nonworld_length_mm"],
        "rays_with_mismatch_over_1nm": r["rays_with_material_mismatch"].get("> 1e-06 mm"),
        "volume_id_mismatch_fraction_nonworld": r["volume_mismatch_fraction_nonworld"],
        "point_disagreements": p["n_disagree"],
        "point_pairs": dict(sorted(p["disagreements"].items(), key=lambda x: -x[1])[:8]),
        "boundary_f1_at_1nm": h["boundary_f1_at_1nm"], "boundary_f1_at_1um": h["boundary_f1_at_1um"],
        "max_boundary_shift_mm": h["max_boundary_shift_mm"],
        "sequence_identical_fraction": h["sequence_identical_fraction"],
        "materials_iou_below_1": {k: {"one_minus_iou": 1 - v["iou"], "len_a_mm": v["len_a_mm"],
                                      "len_b_mm": v["len_b_mm"]}
                                  for k, v in sorted(iou.items(), key=lambda x: x[1]["iou"])},
        "worst_parts": worst,
        "robustness": {"rays_missing": {k: len(v) for k, v in m["rays_missing"].items()},
                       "hung": {k: len(v) for k, v in m["hung"].items()},
                       "killed_events": {"a": len(ma.get("killed_events", [])), "b": len(mb.get("killed_events", []))},
                       "g4_exceptions": {"a": ma.get("g4_exceptions", {}), "b": mb.get("g4_exceptions", {})},
                       "zero_length_steps": {"a": ma.get("zero_length_steps"), "b": mb.get("zero_length_steps")}},
        "export": exp,
        "step_file": rel(c.step), "material_map": rel(c.matmap),
        "available": {"triage": os.path.isfile(os.path.join(c.dir, "cmp", "triage", "triage.json")),
                      "rescans": sorted(d for d in os.listdir(inv) if os.path.isdir(os.path.join(inv, d)))
                      if os.path.isdir(inv) else [],
                      "other_comparisons": sorted(d for d in os.listdir(c.dir) if d.startswith("cmp_"))},
        "no_difference_test": {"passes": bool(nd),
                               "criteria": "1 - material agreement < 1e-8, no point disagreements, boundary F1 @1nm = 1"},
    }
    return out


@tool("mismatch_families",
      "Groups every mismatching ray interval (longer than min_mm) by (material in A -> material in B, volume "
      "in A, volume in B) and describes each family: length and share of the case's mismatch, rays, parts "
      "involved, where (r, z, phi of interval midpoints), example intervals to replay, and boundary_pattern: "
      "how the interval ends line up with the step boundaries of the two worlds. same_boundaries = both ends "
      "are boundaries in A and in B (only the material differs: wrong material, or navigation that walks a "
      "part as the world). only_in_a = A has boundaries B lacks (part absent in B). only_in_b = B has "
      "boundaries A lacks (extra part in B). one_end_each = one end is an A boundary, the other a B boundary "
      "(a boundary moved: placement or size). Covers small families that triage's modes skip.",
      {"case": CASE, "min_mm": P("number", "ignore intervals shorter than this (mm); default 1e-6"),
       "limit": P("integer", "families to list, largest first; default 12"),
       "examples": P("integer", "example intervals per family; default 3")})
def mismatch_families(case, min_mm=1e-6, limit=12, examples=3):
    import numpy as np
    c = Case(case)
    A, B, names = load_scans(c)
    (evt, s, ln, ma, mb, va, vb), total = mismatch_intervals(c, A, B)
    sel = ln > min_mm
    below = total - float(ln[sel].sum())
    evt, s, ln, ma, mb, va, vb = (x[sel] for x in (evt, s, ln, ma, mb, va, vb))
    res = {"case": rel(c.dir), "mismatch_length_mm": total, "length_below_min_mm": below, "intervals": int(len(ln))}
    if not len(ln):
        res["families"] = []
        return res
    ea, sa = step_edges(A)
    eb, sb = step_edges(B)

    def edge(ge, gs, x):
        return nearest_in_group(ge, gs, evt, x)[1] <= SAME_MM

    ends_a = edge(ea, sa, s) & edge(ea, sa, s + ln)
    ends_b = edge(eb, sb, s) & edge(eb, sb, s + ln)
    pattern = np.where(ends_a & ends_b, 0, np.where(ends_a, 1, np.where(ends_b, 2, 3)))
    pat_names = ["same_boundaries", "only_in_a", "only_in_b", "one_end_each"]
    ri = np.searchsorted(A.rays["evt"], evt)
    o = np.stack([A.rays["ox"], A.rays["oy"], A.rays["oz"]], 1)[ri]
    u = np.stack([A.rays["dx"], A.rays["dy"], A.rays["dz"]], 1)[ri]
    mid = o + u * (s + 0.5 * ln)[:, None]
    r, z, phi = np.hypot(mid[:, 0], mid[:, 1]), mid[:, 2], np.degrees(np.arctan2(mid[:, 1], mid[:, 0]))
    # family label: materials and volume names (without ids)
    nm = len(names)
    lab = {}

    def vlabel(vids, mats):
        key = (vids.astype(np.int64) + 10) * (nm + 2) + (mats + 1)
        uk, inv = np.unique(key, return_inverse=True)
        idx = []
        for k in uk.tolist():
            v, mm = k // (nm + 2) - 10, k % (nm + 2) - 1
            n = c.vname(v, names[mm])
            idx.append(lab.setdefault(n, len(lab)))
        return np.array(idx, np.int64)[inv]

    la, lb = vlabel(va, ma), vlabel(vb, mb)
    L = len(lab) + 1
    fkey = ((ma + 1) * (nm + 2) + (mb + 1)) * L * L + la * L + lb
    uk, inv = np.unique(fkey, return_inverse=True)
    flen = np.bincount(inv, weights=ln)
    order = np.argsort(-flen)
    labn = {i: n for n, i in lab.items()}
    rng = np.random.default_rng(1)
    fams = []
    for f in order[:limit]:
        k = np.nonzero(inv == f)[0]
        i0 = k[0]
        w = ln[k]
        pl = np.bincount(pattern[k], weights=w, minlength=4) / w.sum()
        best = int(np.argmax(pl))
        # examples: longest interval, median, random; distinct rays
        cand = [k[np.argmax(w)], k[np.argsort(w)[len(k) // 2]]] + list(rng.permutation(k)[:20])
        ex, seen = [], set()
        for j in cand:
            if evt[j] in seen:
                continue
            seen.add(evt[j])
            ex.append({"evt": int(evt[j]), "s0_mm": float(s[j]), "len_mm": float(ln[j]),
                       "part_a": c.pname(int(va[j]), names[ma[j]]), "part_b": c.pname(int(vb[j]), names[mb[j]]),
                       "s_range_mm": [max(0.0, float(s[j]) - 1.0), float(s[j] + ln[j]) + 1.0]})
            if len(ex) >= examples:
                break
        fams.append({
            "family": f"{names[ma[i0]]} -> {names[mb[i0]]} in {labn[la[i0]]} (A) / {labn[lb[i0]]} (B)",
            "material_a": names[ma[i0]], "material_b": names[mb[i0]],
            "length_mm": float(w.sum()), "share": float(w.sum() / total), "intervals": int(len(k)),
            "rays": int(len(np.unique(evt[k]))),
            "interval_mm": {"median": float(np.median(w)), "max": float(w.max())},
            "boundary_pattern": pat_names[best] if pl[best] >= 0.8 else "mixed",
            "boundary_pattern_shares": {pat_names[i]: float(pl[i]) for i in range(4) if pl[i] > 0},
            "parts_a": top_parts(c, va[k], ma[k], w, names), "parts_b": top_parts(c, vb[k], mb[k], w, names),
            "region": {"r_mm": [float(r[k].min()), float(r[k].max())], "z_mm": [float(z[k].min()), float(z[k].max())],
                       "phi_deg": [float(phi[k].min()), float(phi[k].max())]},
            "examples": ex,
        })
    res["families_total"] = int(len(uk))
    res["families"] = fams
    return res


@tool("point_mismatches",
      "Points where the two worlds report different materials (G4Navigator point location), grouped by "
      "(material A -> material B, volume A, volume B) with counts and example coordinates. Points rarely fall "
      "on a surface, so a navigation artifact usually shows no point disagreements while a wrong material or "
      "a moved or missing part does.",
      {"case": CASE, "limit": P("integer", "groups to list; default 10")})
def point_mismatches(case, limit=10):
    import numpy as np
    c = Case(case)
    A, B, names = load_scans(c)
    pa, pb = A.pmat, B.pmat
    bad = pa != pb
    res = {"case": rel(c.dir), "points": int(len(pa)), "disagree": int(bad.sum()),
           "same_material_other_volume": int(((pa == pb) & (A.points["vol"] != B.points["vol"])).sum())}
    groups = {}
    for i in np.nonzero(bad)[0][:200000]:
        key = (names[pa[i]], names[pb[i]], c.vname(int(A.points["vol"][i]), names[pa[i]]),
               c.vname(int(B.points["vol"][i]), names[pb[i]]))
        g = groups.setdefault(key, {"n": 0, "parts_a": {}, "parts_b": {}, "examples": []})
        g["n"] += 1
        for side, S, m in (("parts_a", A, pa), ("parts_b", B, pb)):
            pn = c.pname(int(S.points["vol"][i]), names[m[i]])
            g[side][pn] = g[side].get(pn, 0) + 1
        if len(g["examples"]) < 3:
            g["examples"].append([float(A.points["x"][i]), float(A.points["y"][i]), float(A.points["z"][i])])
    res["groups"] = [{"pair": f"{k[0]} -> {k[1]}", "volume_a": k[2], "volume_b": k[3], "points": g["n"],
                      "parts_a": dict(sorted(g["parts_a"].items(), key=lambda x: -x[1])[:5]),
                      "parts_b": dict(sorted(g["parts_b"].items(), key=lambda x: -x[1])[:5]),
                      "examples_xyz_mm": g["examples"]}
                     for k, g in sorted(groups.items(), key=lambda x: -x[1]["n"])[:limit]]
    return res


@tool("ray",
      "Side-by-side listing of one ray in both worlds as the exact overlay of their step boundaries: each row "
      "is an interval [s0, s0+len] along the ray with A's and B's material and part, whether they differ, and "
      "which world has a step boundary at the row start (edge A, B or A+B). Without smin/smax, shows the rows "
      "within context_mm of each mismatch. Shows whether boundaries match and only the material differs, a "
      "boundary is offset (a short row whose start is an edge of one world only) or a part is absent.",
      {"case": CASE, "evt": P("integer", "ray number"), "smin": P("number", "start of the s range (mm)"),
       "smax": P("number", "end of the s range (mm)"), "max_rows": P("integer", "default 40"),
       "context_mm": P("number", "rows kept around each mismatch when no range is given; default 1")},
      required=("case", "evt"))
def ray(case, evt, smin=None, smax=None, max_rows=40, context_mm=1.0):
    import numpy as np
    from types import SimpleNamespace
    from compare import overlay
    c = Case(case)
    A, B, names = load_scans(c)
    sub, hdr = [], {}
    for k, S in (("a", A), ("b", B)):
        hit = S.rays[S.rays["evt"] == evt]
        if not len(hit):
            raise ToolError(f"ray {evt} is not in both scans (hung, or out of range 0..{int(S.rays['evt'].max())})")
        keep = S.steps["evt"] == evt
        sub.append(SimpleNamespace(steps=S.steps[keep], smat=S.smat[keep]))
        hdr[k] = {"steps": int(keep.sum()), "length_mm": float(hit[0]["total_len"])}
    r0 = A.rays[A.rays["evt"] == evt][0]
    e, s, ln, ma, mb, va, vb = overlay(sub[0], sub[1], 0, 0)
    ea, sa = step_edges(sub[0])
    eb, sb = step_edges(sub[1])
    big = ln > SAME_MM  # shorter intervals are floating-point slivers between equal boundaries
    e, s, ln, ma, mb, va, vb = (x[big] for x in (e, s, ln, ma, mb, va, vb))
    at_a = nearest_in_group(ea, sa, e, s)[1] <= SAME_MM
    at_b = nearest_in_group(eb, sb, e, s)[1] <= SAME_MM
    mis = ma != mb
    res = {"case": rel(c.dir), "evt": int(evt), "origin_mm": [float(r0["ox"]), float(r0["oy"]), float(r0["oz"])],
           "direction": [float(r0["dx"]), float(r0["dy"]), float(r0["dz"])], "a": hdr["a"], "b": hdr["b"],
           "mismatch_mm": float(ln[mis].sum()), "mismatch_intervals": int(mis.sum()),
           "note": f"intervals shorter than {SAME_MM:g} mm (rounding) are left out"}
    if mis.any():
        res["first_mismatch_s_mm"] = float(s[mis][0])
        res["last_mismatch_s_mm"] = float((s + ln)[mis][-1])
    if smin is not None or smax is not None:
        lo, hi = (smin if smin is not None else -1.0), (smax if smax is not None else 1e99)
        keep = (s + ln >= lo) & (s <= hi)
    else:
        keep = np.zeros(len(s), bool)
        for a0, a1 in zip(s[mis], (s + ln)[mis]):
            keep |= (s + ln >= a0 - context_mm) & (s <= a1 + context_mm)
    idx = np.nonzero(keep)[0]
    if len(idx) > max_rows:  # keep mismatching rows first, then their neighbours
        pri = np.where(mis[idx], 0, 1)
        idx = np.sort(idx[np.argsort(pri, kind="stable")[:max_rows]])
        res["rows_omitted"] = int(keep.sum() - max_rows)

    def side(m, v):
        return "<outside>" if m < 0 else f"{names[m]} {c.pname(int(v), names[m])}"

    res["rows"] = [{"s0_mm": float(s[i]), "len_mm": float(ln[i]), "a": side(ma[i], va[i]), "b": side(mb[i], vb[i]),
                    "mismatch": bool(mis[i]), "edge": ("A+B" if at_a[i] and at_b[i] else "A" if at_a[i] else "B")}
                   for i in idx]
    return res


PV_RE = re.compile(r"^av_\d+_impr_\d+_(.+)_pv_\d+$")
DBG_PT = re.compile(r"^s=(\S+) step=(\S+) in (\S+)\s+id=(-?\d+) mat=(\S+) p=\((\S+),(\S+),(\S+)\)")
DBG_LOC = re.compile(r"^\s+local p=\((.*?)\)\s+bbox=\((.*?)\)")
DBG_PART = re.compile(r"^\s+(\S+)\s+id=(-?\d+)\s+(\S+)\s+safety=(\S+)\s+DistIn\(p,v\)=(\S+)\s+DistOut\(p,v\)=(\S+)")


def scan_ray_args(c):
    m = c.meta_a
    a = ["--compact", m["compact"], "--ray-mode", m["ray_mode"], "--seed", str(m["seed"]),
         "--box", ",".join(repr(float(v)) for v in m["box_mm"])]
    if m["ray_mode"] == "aim":
        a += ["--target", ",".join(repr(float(v)) for v in m["target_mm"])]
    return a


def part_of(pv):
    if pv in ("simcadval_world", "world_volume"):
        return "<world>"
    mm = PV_RE.match(pv)
    return mm.group(1) if mm else pv


@tool("debug_ray",
      "Replays one ray in one world (default the STEP world B) with a bare G4Navigator and, at every step "
      "start in [smin, smax], asks each part whose solid is within 1 mm what it answers: inside/surface/"
      "outside, safety, distance in and out along the ray, the point in the part's local frame and the "
      "part's local bounding box. Shows which parts are at a mismatch and whether a part answers wrongly. A "
      "part answering outside with safety 0 on its own face also happens where B enters correctly, so compare "
      "a failing ray with one that enters the same kind of part correctly. Keep smax - smin <= 20 mm.",
      {"case": CASE, "evt": P("integer", "ray number"), "smin": P("number", "mm"), "smax": P("number", "mm"),
       "world": P("string", "step (B, default) or dd4hep (A)", enum=["step", "dd4hep"]),
       "max_points": P("integer", "step points to return; default 12"),
       "max_parts": P("integer", "nearby parts per point; default 4")},
      required=("case", "evt", "smin", "smax"))
def debug_ray(case, evt, smin, smax, world="step", max_points=12, max_parts=4):
    c = Case(case)
    if smax - smin > 20 or smax <= smin:
        raise ToolError("need 0 < smax - smin <= 20 mm")
    scanner = os.environ.get("SCANNER") or os.path.join(REPO, "build", "simcadval_scan")
    if not os.access(scanner, os.X_OK):
        raise ToolError(f"{scanner} missing: build first (README, Quick start)")
    with tempfile.TemporaryDirectory() as td:
        cmd = [scanner, "--world", world, "--out", os.path.join(td, "x")] + scan_ray_args(c)
        if world == "step":
            cmd += ["--step", c.step, "--matmap", c.matmap or "none"]
        cmd += ["--debug-ray", f"{int(evt)},{smin!r},{smax!r}"]
        pr = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    if pr.returncode != 0:
        raise ToolError(f"scanner failed (rc {pr.returncode}): {pr.stderr[-800:]}")
    pts, loc = [], None
    for line in pr.stdout.splitlines():
        if (mm := DBG_PT.match(line)):
            s, st, pv, vid, mat, x, y, z = mm.groups()
            pts.append({"s_mm": float(s), "step_mm": float(st), "in": part_of(pv), "id": int(vid), "material": mat,
                        "p_mm": [float(x), float(y), float(z)], "near": []})
        elif (mm := DBG_LOC.match(line)):
            lp = [float(v) for v in mm.group(1).split(",")]
            bb = [[float(a) for a in seg.strip().split("..")] for seg in mm.group(2).split(", ")]
            loc = (lp, bb)
        elif (mm := DBG_PART.match(line)) and pts:
            pv, vid, ins, saf, din, dout = mm.groups()
            pts[-1]["near"].append({"part": part_of(pv), "id": int(vid), "inside": ins.lower(),
                                    "safety_mm": float(saf), "dist_in_mm": float(din),
                                    "dist_out_mm": None if float(dout) < 0 else float(dout),
                                    **({"local_p_mm": loc[0], "local_bbox_mm": loc[1]} if loc else {})})
            loc = None
    for p in pts:
        p["near"] = sorted(p["near"], key=lambda q: q["safety_mm"])[:max_parts]
    return {"case": rel(c.dir), "world": world, "evt": int(evt), "points_in_range": len(pts),
            "points": pts[:max_points], "truncated": len(pts) > max_points,
            "note": "s/step: navigator step start and length; in: volume the navigator is in; near: parts with "
                    "safety <= 1 mm, nearest first; dist_in/out along the ray, null = infinite / not inside"}


@tool("material_map",
      "Compares the material map the STEP world uses (stepName -> geant4Name) with the DD4hep material of "
      "each part (export report). differences lists every part whose map material differs from its DD4hep "
      "material; parts gives the entries for the requested parts. Use for any finding where a part has the "
      "wrong material in B.",
      {"case": CASE, "parts": {"type": "array", "items": {"type": "string"},
                               "description": "parts to look up: name@id, name@* or name"},
       "limit": P("integer", "differences to list; default 30")})
def material_map(case, parts=None, limit=30):
    c = Case(case)
    if not c.matmap or not os.path.isfile(c.matmap):
        raise ToolError(f"no material map file ({c.matmap!r}); the STEP world then uses the XDE materials")
    mp = dict(re.findall(r'stepName="([^"]+)"\s+geant4Name="([^"]+)"', open(c.matmap).read()))
    ref = {f"{v['name']}@{v['id']}": v["material"] for v in c.volumes if v["exported"]}
    diff = [{"part": k, "material_a": ref[k], "material_b": v} for k, v in mp.items() if k in ref and ref[k] != v]
    res = {"case": rel(c.dir), "map": rel(c.matmap), "entries": len(mp), "parts_exported": len(ref),
           "differences_total": len(diff), "differences": diff[:limit],
           "exported_parts_missing_from_map": [k for k in ref if k not in mp][:limit],
           "map_entries_unknown_to_export": [k for k in mp if k not in ref][:limit]}
    if parts:
        res["parts"], seen = [], set()
        for q in parts:
            ids = c.resolve(q)
            if not ids:
                res["parts"].append({"part": q, "error": "not in the export report"})
            for i in ids[:limit]:
                k = c.pname(i)
                if k not in seen:
                    seen.add(k)
                    res["parts"].append({"part": k, "material_a": c.vol[i]["material"], "material_b": mp.get(k)})
    return res


@tool("part_info",
      "Looks up parts in the export report by name@id, name@*, name, id or name substring: material, shape, "
      "export status, volumes (TGeo expected, OCCT, STEP read-back). With placements = true (default, at most "
      "4 matched volumes) also lists each placement in the reference geometry: path, global bounding-box centre "
      "and half-sizes (mm), translation and rotation. Use the centre and half-sizes as a rescan target.",
      {"case": CASE, "query": P("string", "name@id, name@*, name, id or substring"),
       "placements": P("boolean", "default true"), "max_placements": P("integer", "per volume; default 6"),
       "limit": P("integer", "volumes to list; default 20")},
      required=("case", "query"))
def part_info(case, query, placements=True, max_placements=6, limit=20):
    c = Case(case)
    ids = c.resolve(query)
    if not ids:
        raise ToolError(f"no part matches {query!r}")
    vols = []
    for i in ids[:limit]:
        v = c.vol[i]
        vols.append({"part": c.pname(i), "material": v["material"], "shape": v["shape"], "status": v["status"],
                     "exported": v["exported"], "readback_ok": v["readback_ok"],
                     "volume_mm3": {"expected": v["expect_mm3"], "occt": v["occt_mm3"], "readback": v["readback_mm3"]},
                     "daughters": v["ndaughters"]})
    res = {"case": rel(c.dir), "matches": len(ids), "parts": vols}
    if placements and len(ids) <= 4:
        g = geometry(c.meta_a["compact"])
        for v, i in zip(vols, ids):
            n = g.Placements(i, max_placements)
            pl = list(g.gPl)
            v["local_half_mm"] = list(g.Half(i))
            v["placements_total"] = n
            v["placements"] = []
            for k, path in enumerate(list(g.gPath)):
                q = pl[18 * k:18 * (k + 1)]
                rot = [q[9:12], q[12:15], q[15:18]]
                ident = all(abs(rot[a][b] - (a == b)) < 1e-12 for a in range(3) for b in range(3))
                v["placements"].append({"path": str(path), "bbox_centre_mm": q[0:3], "bbox_half_mm": q[3:6],
                                        "translation_mm": q[6:9], "rotation": "identity" if ident else rot})
    return res


@tool("boundary_shifts",
      "Measures how far boundaries moved between A and B. Each material-change boundary of A is matched to the "
      "nearest B boundary with the same transition on the same ray: exact (|ds| <= min_mm), shifted "
      "(min_mm < |ds| <= max_mm) or unmatched (no such B boundary within max_mm: lost face). Filter by part "
      "(a boundary of the part in A) and/or transition. With fit = true, takes the face normal n at each "
      "shifted boundary from the reference geometry; a face moved by D shifts the boundary by ds = (D.n)/(u.n) "
      "along a ray with direction u, so D.n = ds (u.n) for every ray through that face. Returns per-face "
      "groups (normal, rays, mean and spread of D.n), the least-squares translation D, and a rigid fit "
      "(translation + small rotation about the boundaries' centroid). A small translation residual means a "
      "rigid shift; opposite faces moving apart means a size change. per_placement fits each placement "
      "separately (copies rotated differently).",
      {"case": CASE, "part": P("string", "name@id, name@* or name"),
       "transition": P("string", "e.g. 'Carbon->Air' (material before -> after along the ray, in A)"),
       "min_mm": P("number", "default 1e-7"), "max_mm": P("number", "default 1.0"),
       "fit": P("boolean", "default true"), "per_placement": P("boolean", "default false"),
       "max_boundaries": P("integer", "shifted boundaries used for the fit; default 4000")})
def boundary_shifts(case, part=None, transition=None, min_mm=1e-7, max_mm=1.0, fit=True, per_placement=False,
                    max_boundaries=4000):
    import numpy as np
    c = Case(case)
    A, B, names = load_scans(c)

    def bnds(S):
        st, gm = S.steps, S.smat
        k = np.nonzero(np.r_[False, (st["evt"][1:] == st["evt"][:-1]) & (gm[1:] != gm[:-1])])[0]
        return {"evt": st["evt"][k].astype(np.int64), "s": st["s0"][k], "fm": gm[k - 1], "to": gm[k],
                "vbefore": st["vol"][k - 1], "vafter": st["vol"][k]}

    ba, bb = bnds(A), bnds(B)
    K = len(names) + 2

    def gid(b):
        return (b["evt"] * K + (b["fm"] + 1)) * K + (b["to"] + 1)

    def match(q, t):
        gt = gid(t)
        o = np.lexsort((t["s"], gt))
        j, d = nearest_in_group(gt[o], t["s"][o], gid(q), q["s"])
        ds = np.where(j >= 0, t["s"][o][np.maximum(j, 0)] - q["s"], np.nan)
        return ds

    sel = np.ones(len(ba["s"]), bool)
    selb = np.ones(len(bb["s"]), bool)
    ids = None
    if part:
        ids = c.resolve(part)
        if not ids:
            raise ToolError(f"no part matches {part!r}")
        sel &= np.isin(ba["vbefore"], ids) | np.isin(ba["vafter"], ids)
        selb &= np.isin(bb["vbefore"], ids) | np.isin(bb["vafter"], ids)
    if transition:
        mm = re.split(r"\s*(?:->|→)\s*", transition.strip())
        if len(mm) != 2 or mm[0] not in names or mm[1] not in names:
            raise ToolError(f"transition {transition!r}: use 'MatA->MatB' with materials {names[:-1]}")
        fa, ta = names.index(mm[0]), names.index(mm[1])
        sel &= (ba["fm"] == fa) & (ba["to"] == ta)
        selb &= (bb["fm"] == fa) & (bb["to"] == ta)
    ds = match(ba, bb)
    dsb = match(bb, ba)
    ad = np.abs(ds)
    exact = sel & (ad <= min_mm)
    shifted = sel & (ad > min_mm) & (ad <= max_mm)
    unmatched = sel & ~(ad <= max_mm)
    res = {"case": rel(c.dir), "part": part, "transition": transition, "part_ids": ids[:20] if ids else None,
           "boundaries_a": int(sel.sum()), "exact": int(exact.sum()), "shifted": int(shifted.sum()),
           "unmatched_a": int(unmatched.sum()), "boundaries_b": int(selb.sum()),
           "unmatched_b": int((selb & ~(np.abs(dsb) <= max_mm)).sum())}
    if shifted.any():
        d = ds[shifted]
        res["ds_mm"] = {"p5": float(np.percentile(d, 5)), "median": float(np.median(d)),
                        "p95": float(np.percentile(d, 95)), "max_abs": float(np.abs(d).max())}
        tr = {}
        for f, t, x in zip(ba["fm"][shifted].tolist(), ba["to"][shifted].tolist(), d.tolist()):
            tr.setdefault(f"{names[f]}->{names[t]}", []).append(x)
        res["by_transition"] = [{"transition": k, "n": len(v), "median_ds_mm": float(np.median(v))}
                                for k, v in sorted(tr.items(), key=lambda x: -len(x[1]))[:8]]
    if not (fit and shifted.sum() >= 3):
        if fit:
            res["fit"] = "need at least 3 shifted boundaries"
        return res
    k = np.nonzero(shifted)[0]
    if len(k) > max_boundaries:
        k = np.sort(np.random.default_rng(1).choice(k, max_boundaries, replace=False))
    ri = np.searchsorted(A.rays["evt"], ba["evt"][k])
    o = np.stack([A.rays["ox"], A.rays["oy"], A.rays["oz"]], 1)[ri]
    u = np.stack([A.rays["dx"], A.rays["dy"], A.rays["dz"]], 1)[ri]
    p = o + u * ba["s"][k][:, None]
    side = (np.isin(ba["vafter"][k], ids) if ids else ba["vafter"][k] >= 0).astype(np.int32)
    g = geometry(c.meta_a["compact"])
    n = len(k)
    nrm, dist, vid = np.zeros(3 * n), np.zeros(n), np.zeros(n, np.int32)
    g.Normals(n, np.ascontiguousarray(p).ravel(), np.ascontiguousarray(u).ravel(), side, nrm, dist, vid)
    nrm = nrm.reshape(n, 3)
    paths = [str(x) for x in g.gPath]
    ok = (np.abs(np.linalg.norm(nrm, axis=1) - 1) < 1e-6) & (dist >= 0) & (dist < 1e-6)
    if ids:
        ok &= np.isin(vid, ids)
    cdot = ds[k] * np.einsum("ij,ij->i", u, nrm)
    nrm, cdot, p, kk = nrm[ok], cdot[ok], p[ok], np.nonzero(ok)[0]
    paths = [paths[i] for i in kk]
    res_fit = {"boundaries_used": int(ok.sum()), "normal_lookup_failed": int((~ok).sum())}
    if ok.sum() < 3:
        res["fit"] = res_fit
        return res

    def solve(N, cvec, P=None):
        """least squares with outlier rejection (residual > 5 MAD, twice): x, rms, rank, outliers"""
        M = N if P is None else np.hstack([N, np.cross(P - P.mean(0), N)])
        inl = np.ones(len(cvec), bool)
        for _ in range(3):
            x, *_ = np.linalg.lstsq(M[inl], cvec[inl], rcond=None)
            r = np.abs(M @ x - cvec)
            mad = np.median(r[inl])
            new = r <= max(5 * 1.4826 * mad, 1e-9)
            if new.sum() < 3 or np.array_equal(new, inl):
                break
            inl = new
        _, sv, vt = np.linalg.svd(M[inl], full_matrices=False)
        rms = float(np.sqrt(np.mean((M[inl] @ x - cvec[inl]) ** 2)))
        rank = int((sv > 1e-6 * sv[0]).sum())
        return x, rms, rank, int((~inl).sum()), vt[:rank]

    D, rms_t, rank_t, out_t, dirs_t = solve(nrm, cdot)
    x6, rms_r, rank_r, out_r, _ = solve(nrm, cdot, p)
    scale = float(np.median(np.abs(cdot)))
    # faces: sign-normalised normals rounded to 0.01
    big = np.argmax(np.abs(nrm), axis=1)
    sgn = np.sign(nrm[np.arange(len(nrm)), big])
    ns, cs = nrm * sgn[:, None], cdot * sgn
    key = np.round(ns, 2)
    faces = {}
    for kk_, cc in zip(map(tuple, key.tolist()), cs.tolist()):
        faces.setdefault(kk_, []).append(cc)
    res_fit.update({
        "faces": [{"normal": list(f), "rays": len(v), "median_D_dot_n_mm": float(np.median(v)),
                   "mad_mm": float(np.median(np.abs(np.array(v) - np.median(v))))}
                  for f, v in sorted(faces.items(), key=lambda x: -len(x[1]))[:8]],
        "translation": {"D_mm": D.tolist(), "rms_residual_mm": rms_t, "rank": rank_t, "outliers": out_t,
                        "determined_along": np.round(dirs_t, 3).tolist()},
        "rigid": {"D_mm": x6[:3].tolist(), "rotation_mrad": (1e3 * x6[3:]).tolist(), "about_mm": p.mean(0).tolist(),
                  "rms_residual_mm": rms_r, "rank": rank_r, "outliers": out_r},
        "median_abs_D_dot_n_mm": scale,
        "note": "D is measured only along determined_along (the normals of the faces seen); components "
                "across those directions are unmeasured and reported as 0",
    })
    if rms_t > 0.1 * scale:
        res_fit["hint"] = ("translation fit leaves a large residual: faces move by different amounts "
                           "(size change, rotation, or copies moved differently: try per_placement=true)")
    if per_placement:
        byp = {}
        for i, pa in enumerate(paths):
            byp.setdefault(pa, []).append(i)
        res_fit["placements"] = []
        for pa, ii in sorted(byp.items(), key=lambda x: -len(x[1]))[:6]:
            if len(ii) < 3:
                continue
            Dp, rp, kp, op, _ = solve(nrm[ii], cdot[ii])
            res_fit["placements"].append({"path": pa, "rays": len(ii), "D_mm": Dp.tolist(), "rms_residual_mm": rp,
                                          "rank": kp, "outliers": op})
    res["fit"] = res_fit
    return res


@tool("triage",
      "Statistical triage of the mismatch (scripts/triage.py; runs it first if needed, which takes a few "
      "minutes): strata (every material -> material and solid-class family, including tiny ones), a readable "
      "decision rule with held-out precision/recall, the most predictive single features, and failure modes "
      "(clusters of mismatching steps) with a representative step each. Rules and modes are correlations: "
      "confirm them with ray / debug_ray.",
      {"case": CASE, "run": P("boolean", "run triage.py if its output is missing; default true")})
def triage(case, run=True):
    c = Case(case)
    tj = os.path.join(c.dir, "cmp", "triage", "triage.json")
    if not os.path.isfile(tj):
        if not run:
            raise ToolError("triage has not been run on this case")
        pr = subprocess.run([sys.executable, os.path.join(HERE, "triage.py"), c.dir, "--no-plots"],
                            capture_output=True, text=True, timeout=3600)
        if pr.returncode != 0 or not os.path.isfile(tj):
            raise ToolError(f"triage.py failed (rc {pr.returncode}): {(pr.stderr or pr.stdout)[-800:]}")
    d = json.load(open(tj))
    out = {"case": rel(c.dir), "mismatch_length_mm": d.get("mismatch_length_mm"),
           "positive_steps": d.get("positive_steps"), "strata": d.get("strata", [])[:20],
           "relocation_agreement": d.get("relocation_agreement")}
    if d.get("single_feature"):
        out["top_features"] = [{"feature": f["label"], "auc": f["auc"], "higher_means": f.get("higher_means")}
                               for f in d["single_feature"][:6]]
    tr = d.get("trees", {}).get("mechanism")
    if tr:
        out["rule_tree"] = {"test_auc": tr.get("test_auc"), "rules": [
            {"rule": r["rule"], "test_precision": r["test"]["precision"], "test_recall": r["test"]["recall"],
             "test_mismatch_share": r["test"]["mismatch_length_share"]} for r in tr.get("rules", [])[:4]]}
    out["modes"] = [{"id": m["id"], "share": m["mismatch_length_share"], "steps": m["steps"],
                     "signature": {k: v[:2] for k, v in m.get("signature", {}).items()},
                     "distinguishing": m.get("distinguishing", [])[:4],
                     "representative": {k: m["medoid"][k] for k in ("evt", "s0_mm", "len_mm", "volume", "material")
                                        if k in m.get("medoid", {})}}
                    for m in d.get("clusters", [])]
    if "listed" in d:
        out["listed_steps"] = d["listed"][:30]
    return out


@tool("rescan",
      "Runs a targeted rescan (scripts/run_case.sh) of this case's two worlds with the same STEP file and "
      "material map: rays start uniformly in box and aim at a uniform point of target, no points. Use it to "
      "measure a shift over many rays through one part, or to check that a pattern holds across a part; take "
      "target from part_info placements or from mismatch positions. Results go to "
      "<case>/investigation/<name>, which is itself a case for every other tool (e.g. boundary_shifts). At "
      "most 20000 rays; 5000 take about a minute.",
      {"case": CASE, "name": P("string", "new rescan name, letters/digits/_/-"),
       "target": NUMS(6, "cx,cy,cz,hx,hy,hz (mm): centre and half-sizes of the box rays aim at"),
       "box": NUMS(3, "hx,hy,hz (mm): half-sizes of the origin box (default: the case's box)"),
       "rays": P("integer", "default 5000, at most 20000"), "seed": P("integer", "default: the case's seed")},
      required=("case", "name", "target"))
def rescan(case, name, target, box=None, rays=5000, seed=None):
    c = Case(case)
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,40}", name):
        raise ToolError("name: 1-40 letters, digits, _ or -")
    if not 1 <= rays <= MAX_RESCAN_RAYS:
        raise ToolError(f"rays must be 1..{MAX_RESCAN_RAYS}")
    inv = os.path.join(c.root, "investigation")
    dest = os.path.join(inv, name)
    if os.path.exists(dest):
        raise ToolError(f"{rel(dest)} exists; choose another name")
    box = box or c.meta_a["box_mm"]
    env = dict(os.environ, OUT=inv, STEP=c.step, MATMAP=c.matmap or "")
    cmd = [os.path.join(HERE, "run_case.sh"), name, c.meta_a["compact"], "--rays", str(int(rays)),
           "--ray-mode", "aim", "--box", ",".join(repr(float(v)) for v in box),
           "--target", ",".join(repr(float(v)) for v in target), "--points", "0",
           "--seed", str(int(seed if seed is not None else c.meta_a["seed"]))]
    t0 = time.time()
    pr = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=3600)
    os.makedirs(dest, exist_ok=True)
    open(os.path.join(dest, "run.log"), "w").write(pr.stdout + pr.stderr)
    if pr.returncode != 0:
        raise ToolError(f"run_case.sh failed (rc {pr.returncode}): {(pr.stdout + pr.stderr)[-800:]}")
    s = case_summary(dest)
    keep = ("case", "rays", "material_agreement_nonworld", "mismatch_length_mm", "nonworld_length_mm",
            "rays_with_mismatch_over_1nm", "boundary_f1_at_1nm", "max_boundary_shift_mm", "worst_parts",
            "materials_iou_below_1", "robustness")
    fam = mismatch_families(dest, limit=6, examples=1)
    return {"rescan": rel(dest), "seconds": round(time.time() - t0, 1), **{k: s[k] for k in keep},
            "families": [{k: f[k] for k in ("family", "length_mm", "share", "rays", "boundary_pattern", "parts_a",
                                            "parts_b")} for f in fam["families"]]}


# ── diagnosis ─────────────────────────────────────────────────────────────────

def read_calls(root):
    p = os.path.join(root, "investigation", "calls.jsonl")
    calls = {}
    if os.path.isfile(p):
        for ln in open(p):
            try:
                r = json.loads(ln)
                calls[r["id"]] = r
            except (ValueError, KeyError):
                pass
    return calls


def validate_diagnosis(c, d):
    errors, warnings = [], []
    try:
        import jsonschema
        v = jsonschema.Draft202012Validator(json.load(open(SCHEMA_PATH)))
        for e in sorted(v.iter_errors(d), key=lambda e: list(e.path)):
            errors.append(f"{'/'.join(map(str, e.path)) or '(top)'}: {e.message}")
    except ImportError:
        warnings.append("jsonschema not installed: schema not checked")
    if not isinstance(d, dict) or not isinstance(d.get("findings"), list):
        return errors or ["(top): findings must be a list"], warnings
    calls = read_calls(c.root)
    fs = [f for f in d["findings"] if isinstance(f, dict)]
    cats = [f.get("category") for f in fs]
    if "no_difference" in cats and len(cats) > 1:
        errors.append("findings: a no_difference finding must stand alone")
    tot = 0.0
    for i, f in enumerate(fs):
        w = f"findings/{i}"
        if isinstance(f.get("share_of_mismatch"), (int, float)):
            tot += f["share_of_mismatch"]
        for e in f.get("evidence") or []:
            cid = e.get("call") if isinstance(e, dict) else None
            if cid not in calls:
                errors.append(f"{w}: evidence call {cid!r} is not in investigation/calls.jsonl")
            elif not calls[cid].get("ok"):
                errors.append(f"{w}: evidence call {cid} failed; cite a successful call")
        for pn in f.get("parts") or []:
            if isinstance(pn, str) and not c.resolve(pn):
                errors.append(f"{w}: part {pn!r} is not in the export report")
        if not f.get("parts") and f.get("category") not in ("no_difference", None):
            errors.append(f"{w}: parts is empty (allowed only for no_difference)")
        mg = f.get("magnitude") or {}
        if f.get("category") == "placement" and not any(k in mg for k in ("displacement_mm", "boundary_shift_mm",
                                                                           "rotation_mrad")):
            warnings.append(f"{w}: placement without displacement_mm / boundary_shift_mm / rotation_mrad")
        if f.get("category") == "material" and "material_b" not in mg and "materials" not in f:
            warnings.append(f"{w}: material finding without materials.b or magnitude.material_b")
    if tot > 1.02:
        errors.append(f"findings: share_of_mismatch adds up to {tot:.3f} > 1")
    dc = d.get("case")
    if isinstance(dc, str) and os.path.abspath(dc) != os.path.abspath(c.root):
        warnings.append(f"case {dc!r} differs from {rel(c.root)!r}; using the latter")
    return errors, warnings


def md_text(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def call_cmd(r):
    a = " ".join(f"--{k} {shlex.quote(json.dumps(v) if isinstance(v, (list, dict)) else str(v))}"
                 for k, v in r.get("args", {}).items())
    return f"scripts/svtool.py {r['tool']} {a}".strip()


def render_md(c, d):
    calls = read_calls(c.root)
    L = [f"# Diagnosis: {rel(c.root)}", "", f"**Verdict:** {md_text(d['verdict'])}", "",
         "| # | category | parts | A → B | magnitude | share of mismatch | confidence |",
         "|---|---|---|---|---|---|---|"]
    for i, f in enumerate(d["findings"], 1):
        parts = ", ".join(f"`{p}`" for p in f["parts"][:6]) + (f" (+{len(f['parts']) - 6} more)"
                                                               if len(f["parts"]) > 6 else "")
        mat = f"{f['materials']['a']} → {f['materials']['b']}" if f.get("materials") else ""
        mg = f.get("magnitude", {})
        mag = md_text(mg.get("summary") or ", ".join(f"{k} {v}" for k, v in mg.items()))
        L.append(f"| {i} | {f['category'].replace('_', ' ')} | {parts or '-'} | {mat} | {mag} | "
                 f"{100 * f['share_of_mismatch']:.4g} % | {f['confidence']} |")
    for i, f in enumerate(d["findings"], 1):
        L += ["", f"## Finding {i}: {f['category'].replace('_', ' ')}", ""]
        if f.get("explanation"):
            L += [md_text(f["explanation"]), ""]
        mg = f.get("magnitude", {})
        extra = {k: v for k, v in mg.items() if k != "summary"}
        if extra:
            L += ["Magnitude: " + "; ".join(f"{k} = {v}" for k, v in extra.items()), ""]
        L.append("Evidence:")
        for e in f["evidence"]:
            r = calls.get(e["call"], {})
            L.append(f"- `{e['call']}` `{call_cmd(r) if r else '?'}`: {md_text(e['shows'])}")
    if d.get("unexplained"):
        u = d["unexplained"]
        L += ["", "## Not explained", ""]
        if "share_of_mismatch" in u:
            L.append(f"Share of mismatch: {100 * u['share_of_mismatch']:.4g} %.")
        if u.get("notes"):
            L.append(md_text(u["notes"]))
    if d.get("next_steps"):
        L += ["", "## Next steps", ""] + [f"- {md_text(s)}" for s in d["next_steps"]]
    L += ["", "---", f"Rendered by `svtool submit_diagnosis` from `cmp/diagnosis.json` (schema {d['schema']}); "
                     "call ids refer to `investigation/calls.jsonl`."
          + (f" Agent: {d['agent'].get('name', '')} {d['agent'].get('model', '')}." if d.get("agent") else "")]
    return "\n".join(L) + "\n"


@tool("submit_diagnosis",
      "Checks a diagnosis against schemas/diagnosis.schema.json and the case (evidence call ids exist in this "
      "investigation's call log, parts exist, shares add up to at most 1), then writes "
      "<case>/cmp/diagnosis.json and renders <case>/cmp/diagnosis.md. On errors nothing final is written "
      "(the draft goes to cmp/diagnosis.draft.json): fix them and submit again. Diagnosis fields: schema = "
      "'simcadval-diagnosis/1', case, verdict (one sentence), findings[] each with category (placement | "
      "material | missing_part | extra_part | navigation_artifact | no_difference | unclear), parts "
      "(['name@id' or 'name@*']), share_of_mismatch (0..1), confidence (high|medium|low), evidence "
      "([{call: 'c007', shows: 'quoted values'}]), optional materials {a, b}, magnitude {summary, "
      "displacement_mm [3], rotation_mrad [3], boundary_shift_mm, material_b, mismatch_length_mm}, explanation; "
      "optional unexplained {share_of_mismatch, notes} and next_steps [].",
      {"case": CASE, "diagnosis": {"type": "object", "description": "the diagnosis (see schemas/diagnosis.schema.json)"}},
      required=("case", "diagnosis"))
def submit_diagnosis(case, diagnosis):
    c = Case(case)
    if c.root != c.dir:
        c = Case(c.root)
    d = json.loads(diagnosis) if isinstance(diagnosis, str) else diagnosis
    errors, warnings = validate_diagnosis(c, d)
    cmp_dir = os.path.join(c.root, "cmp")
    if errors:
        json.dump(d, open(os.path.join(cmp_dir, "diagnosis.draft.json"), "w"), indent=1)
        return {"ok": False, "errors": errors[:25], "warnings": warnings}
    d["case"] = rel(c.root)
    json.dump(d, open(os.path.join(cmp_dir, "diagnosis.json"), "w"), indent=1, ensure_ascii=False)
    open(os.path.join(cmp_dir, "diagnosis.md"), "w").write(render_md(c, d))
    dr = os.path.join(cmp_dir, "diagnosis.draft.json")
    if os.path.exists(dr):
        os.remove(dr)
    return {"ok": True, "warnings": warnings, "written": [rel(os.path.join(cmp_dir, "diagnosis.json")),
                                                          rel(os.path.join(cmp_dir, "diagnosis.md"))]}


# ── driver ────────────────────────────────────────────────────────────────────

def spec(fmt="anthropic"):
    out = []
    for n, t in TOOLS.items():
        if fmt == "openai":
            out.append({"type": "function", "function": {"name": n, "description": t["description"],
                                                         "parameters": t["schema"]}})
        else:
            out.append({"name": n, "description": t["description"], "input_schema": t["schema"]})
    return out


def coerce(name, args):
    """Check and coerce arguments against the tool's schema (lenient with types, strict with names)."""
    sch = TOOLS[name]["schema"]
    props = sch["properties"]
    unknown = [k for k in args if k not in props]
    if unknown:
        raise ToolError(f"unknown argument(s) {unknown} for {name}; allowed: {list(props)}")
    miss = [k for k in sch["required"] if k not in args]
    if miss:
        raise ToolError(f"missing argument(s) {miss} for {name}")
    out = {}
    for k, v in args.items():
        t = props[k].get("type")
        try:
            if v is None:
                out[k] = None
            elif t == "integer":
                out[k] = int(float(v))
            elif t == "number":
                out[k] = float(v)
            elif t == "boolean":
                out[k] = v if isinstance(v, bool) else str(v).lower() in ("1", "true", "yes", "y")
            elif t == "array":
                if isinstance(v, str):
                    v = json.loads(v) if v.strip().startswith("[") else [x for x in v.split(",") if x.strip()]
                it = props[k].get("items", {}).get("type")
                out[k] = [float(x) for x in v] if it == "number" else [str(x) for x in v]
                if "minItems" in props[k] and not props[k]["minItems"] <= len(out[k]) <= props[k]["maxItems"]:
                    raise ValueError(f"needs {props[k]['minItems']} values")
            elif t == "object":
                out[k] = json.loads(v) if isinstance(v, str) else v
            else:
                out[k] = str(v)
                if "enum" in props[k] and out[k] not in props[k]["enum"]:
                    alias = ENUM_ALIASES.get(k, {}).get(out[k].strip().lower())
                    if alias is None:
                        raise ValueError(f"{out[k]!r} is not one of {props[k]['enum']}")
                    out[k] = alias
        except (ValueError, TypeError) as e:
            raise ToolError(f"argument {k}: {e}") from None
    return out


# Other names agents use for enum values; anything else outside the enum is an error. Without this check
# gpt-oss sent world="B" to debug_ray, the scanner silently found nothing, and 13 calls were lost
# (docs/agent-eval.md).
ENUM_ALIASES = {"world": {"a": "dd4hep", "b": "step", "dd4hep": "dd4hep", "step": "step"}}


def parse_flags(name, argv):
    props = TOOLS[name]["schema"]["properties"]
    args, i = {}, 0
    while i < len(argv):
        a = argv[i]
        if not a.startswith("--"):
            raise ToolError(f"unexpected {a!r}; use --name value")
        k = a[2:].replace("-", "_")
        if props.get(k, {}).get("type") == "boolean" and (i + 1 == len(argv) or argv[i + 1].startswith("--")):
            args[k], i = True, i + 1
            continue
        if i + 1 == len(argv):
            raise ToolError(f"{a} needs a value")
        v = argv[i + 1]
        if v.startswith("@") and props.get(k, {}).get("type") == "object":
            v = open(v[1:]).read()
        args[k], i = v, i + 2
    return args


def log_call(root, name, args, result, ok, seconds):
    import fcntl
    inv = os.path.join(root, "investigation")
    os.makedirs(inv, exist_ok=True)
    with open(os.path.join(inv, "calls.jsonl"), "a+") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.seek(0)
        n = sum(1 for _ in f)
        cid = f"c{n + 1:03d}"
        f.write(json.dumps({"id": cid, "time": time.strftime("%Y-%m-%dT%H:%M:%S"), "tool": name, "args": args,
                            "ok": ok, "seconds": round(seconds, 2), "result": result}, ensure_ascii=False) + "\n")
        fcntl.flock(f, fcntl.LOCK_UN)
    return cid


def run(name, args):
    """Run one tool; returns the JSON-safe result (with call_id when logged)."""
    if name not in TOOLS:
        return {"error": f"unknown tool {name!r}; tools: {list(TOOLS)}"}
    t0 = time.time()
    root = None
    try:
        args = coerce(name, args)
        if TOOLS[name]["log"] and "case" in args:
            root = Case(args["case"]).root
        res, ok = compact(TOOLS[name]["fn"](**args)), True
    except ToolError as e:
        res, ok = {"error": str(e)}, False
    except Exception as e:  # report, do not crash the agent loop
        res, ok = {"error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-1500:]}, False
    if root:
        logged = {k: v for k, v in args.items() if k != "diagnosis"} if name == "submit_diagnosis" else args
        res = {"call_id": log_call(root, name, logged, res, ok, time.time() - t0), **res}
    return res


def main(argv):
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(__doc__)
        print("tools: " + ", ".join(TOOLS))
        return 0
    if argv[0] == "spec":
        fmt = argv[2] if len(argv) > 2 and argv[1] == "--format" else "anthropic"
        print(json.dumps(spec(fmt), indent=1))
        return 0
    if argv[0] == "list":
        for n, t in TOOLS.items():
            print(f"{n:18s} {t['description'].split('. ')[0]}")
        return 0
    if "DD4hepINSTALL" not in os.environ:  # not in the container: re-run through eic.sh
        eic = os.path.join(REPO, "eic.sh")
        os.execv(eic, [eic, shlex.join(["python3", os.path.abspath(__file__)] + argv)])
    try:
        if argv[0] == "call":
            if len(argv) < 2:
                raise ToolError("usage: svtool.py call TOOL JSON")
            name = argv[1]
            raw = argv[2] if len(argv) > 2 else "{}"
            raw = sys.stdin.read() if raw == "-" else open(raw[1:]).read() if raw.startswith("@") else raw
            args = json.loads(raw)
        else:
            name = argv[0]
            if name not in TOOLS:
                raise ToolError(f"unknown tool {name!r}; tools: {list(TOOLS)}")
            args = parse_flags(name, argv[1:])
    except (ToolError, ValueError, OSError) as e:
        print(json.dumps({"error": str(e)}))
        return 2
    # keep stdout for the JSON result: library chatter goes to a log file
    out = os.fdopen(os.dup(1), "w")
    logdir = None
    try:
        logdir = os.path.join(Case(args["case"]).root, "investigation") if "case" in args else None
    except Exception:
        pass
    logf = os.path.join(logdir, "svtool.log") if logdir else os.devnull
    if logdir:
        os.makedirs(logdir, exist_ok=True)
    fd = os.open(logf, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    sys.stdout.flush()
    sys.stderr.flush()
    os.dup2(fd, 1)
    os.dup2(fd, 2)
    res = run(name, args)
    out.write(json.dumps(res, ensure_ascii=False) + "\n")
    out.flush()
    return 0 if "error" not in res else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
