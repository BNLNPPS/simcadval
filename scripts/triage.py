#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Mismatch triage: which geometric conditions predict a DD4hep-vs-STEP disagreement?

    triage.py CASE_DIR [--out DIR] [--min-mm 1e-6] [--neg-ratio 20] [--depth 4] [--seed 1]

Run inside the container from the repo root, after run_case.sh / compare.py, e.g.
    ./eic.sh 'python3 scripts/triage.py out/minisid_ip'

Unit of analysis = one step of the reference (A, DD4hep) world: ray `evt`, step `k`,
inside one placed volume. A step is *positive* if the exact overlay (compare.overlay)
shows more than --min-mm of material mismatch inside it. Every positive step and a
random sample of negative steps (weights restore the population) get generic features
computed in the reference TGeo geometry (no detector-specific names):

  what  - material, solid class, which surface the step entered through (own face or a
          daughter's), incidence angle, entry offset from that face, distance to the
          face edge, face size, solid size, placement rotation (angle, alignment to the
          axes), translation, local coordinate magnitude, exit surface, step length.
  where - global position (r, z, phi), ray eta, step index, path length, volume name.

Outputs in CASE_DIR/cmp/triage/: triage.md (report), triage.json, features.csv.gz
(the per-step table, one row per sampled step), feature_auc.png, clusters.png.

Analyses (numpy/scipy only):
  1. single-feature separation: weighted ROC AUC of each feature, mismatch vs. match;
  2. shallow decision tree (CART, weighted Gini, balanced classes), trained on even rays,
     scored on odd rays -> human-readable rules with held-out precision/recall;
     once on "what" features (mechanism) and once on all features (mechanism+location);
  3. failure-mode clusters of positive steps (Ward linkage on standardised "what" features
     plus mismatch shape; k by silhouette, k=1 if no structure), each with its signature
     and a medoid ray to replay with show_ray.py / simcadval_scan --debug-ray.

These are correlations in the reference frame. They show where to look; they do not prove
a cause. Entry offsets are computed from the reconstructed point origin + s0*dir, which has
its own ~1e-13 mm rounding, so they are a noisy proxy for what the candidate navigator saw.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare import Scan, overlay  # noqa: E402

# ── reference-geometry feature extractor (TGeo, compiled by cling) ─────────────

CPP = r"""
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TGeoShape.h>
#include <TGeoMaterial.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace rdtri {
std::unordered_map<const TGeoVolume*, int> gVid;
std::vector<TGeoVolume*> gVols;

// same numbering as src/geo_ids.hh (DFS pre-order over TGeoVolumes, first visit wins)
int Init() {
  gVid.clear(); gVols.clear();
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
std::string VolName(int i) { return gVols[i]->GetName(); }
std::string ShapeClass(int i) { return gVols[i]->GetShape()->ClassName(); }
std::string MatName(int i) { auto* m = gVols[i]->GetMaterial(); return m ? m->GetName() : "<none>"; }

struct Frame { TGeoVolume* vol = nullptr; TGeoHMatrix m; int level = 0; };

bool Locate(const double* pmm, Frame& f) {
  TGeoNode* n = gGeoManager->FindNode(pmm[0] * 0.1, pmm[1] * 0.1, pmm[2] * 0.1);
  if (!n) return false;
  f.vol = n->GetVolume(); f.m = *gGeoManager->GetCurrentMatrix(); f.level = gGeoManager->GetLevel();
  return true;
}
void ToLocal(const Frame& f, const double* pmm, double* lcm) {
  double p[3] = {pmm[0] * 0.1, pmm[1] * 0.1, pmm[2] * 0.1};
  f.m.MasterToLocal(p, lcm);
}
// rotation angle [deg], phi_z of the local x axis [deg], its distance to a multiple of 90 deg,
// max matrix-entry distance from {0,+-1}
void RotFeat(const TGeoHMatrix& m, double& angle, double& phimod, double& alignErr, double& phiz) {
  const double* r = m.GetRotationMatrix();
  const double c = std::max(-1.0, std::min(1.0, 0.5 * (r[0] + r[4] + r[8] - 1)));
  angle = std::acos(c) * 180 / M_PI;
  phiz = std::atan2(r[3], r[0]) * 180 / M_PI;
  const double md = std::fmod(std::fabs(phiz), 90.0);
  phimod = std::min(md, 90.0 - md);
  alignErr = 0;
  for (int i = 0; i < 9; ++i) { const double a = std::fabs(r[i]); alignErr = std::max(alignErr, std::min(a, std::fabs(1 - a))); }
}
double SignedSafetyMM(TGeoShape* s, const double* lcm) {
  const bool in = s->Contains(lcm);
  const double d = s->Safety(lcm, in) * 10;
  return in ? d : -d;
}

const int NF = 29;
std::vector<std::string> Names() {
  return {"level", "own_entry", "cos_in", "normal_axis", "normal_align", "entry_offset_mm",
          "u_x", "u_y", "u_z", "edge_dist_mm", "face_area_mm2", "own_hmin_mm", "own_hmax_mm",
          "own_rot_deg", "own_phimod90_deg", "own_align_err", "own_trans_mm", "own_local_mag_mm",
          "cur_hmin_mm", "cur_hmax_mm", "cur_rot_deg", "cur_phimod90_deg", "cur_align_err",
          "exit_own", "cos_out", "mid_safety_mm", "cur_entry_safety_mm", "own_phiz_deg", "cur_phiz_deg"};
}

// p0/p1: step entry/exit (mm), pprev: a point inside the previous step's volume (NaN if none),
// dir: unit ray direction. out[n*NF]; ivol/iown: volume id of the step / of the entry-face owner.
void Extract(long n, const double* p0, const double* p1, const double* pprev, const double* dir,
             double* out, int* ivol, int* iown) {
  const double tol = 1e-6;  // mm: "this point is on that solid's surface"
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (long i = 0; i < n; ++i) {
    double* o = out + i * NF;
    for (int k = 0; k < NF; ++k) o[k] = nan;
    ivol[i] = iown[i] = -1;
    const double *a = p0 + 3 * i, *b = p1 + 3 * i, *d = dir + 3 * i, *q = pprev + 3 * i;
    const double mid[3] = {0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
    Frame cur;
    if (!Locate(mid, cur)) continue;
    auto it = gVid.find(cur.vol);
    ivol[i] = it == gVid.end() ? -1 : it->second;
    TGeoShape* cs = cur.vol->GetShape();
    double l0[3], l1[3], lm[3], ld[3];
    ToLocal(cur, a, l0); ToLocal(cur, b, l1); ToLocal(cur, mid, lm);
    cur.m.MasterToLocalVect(d, ld);
    const double safe0 = std::fabs(SignedSafetyMM(cs, l0));
    o[26] = safe0;

    // Which solid owns the entry surface: the step's own volume (entering it through its
    // own face) or the previous one (leaving a daughter back into this mother).
    Frame prev;
    const bool hasPrev = !std::isnan(q[0]) && Locate(q, prev);
    const Frame* own = &cur;
    double ol[3] = {l0[0], l0[1], l0[2]}, od[3] = {ld[0], ld[1], ld[2]};
    int ownEntry = safe0 <= tol ? 1 : -1;
    if (safe0 > tol && hasPrev) {
      double pl[3], pd[3];
      ToLocal(prev, a, pl); prev.m.MasterToLocalVect(d, pd);
      const double s = std::fabs(SignedSafetyMM(prev.vol->GetShape(), pl));
      if (s < safe0) {
        own = &prev; ownEntry = 0;
        for (int k = 0; k < 3; ++k) { ol[k] = pl[k]; od[k] = pd[k]; }
      }
    }
    auto jt = gVid.find(own->vol);
    iown[i] = jt == gVid.end() ? -1 : jt->second;
    TGeoShape* os = own->vol->GetShape();
    o[0] = cur.level;
    o[1] = ownEntry;
    double nr[3] = {0, 0, 0};
    os->ComputeNormal(ol, od, nr);
    const double nn = std::sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
    if (nn > 0) {
      for (double& x : nr) x /= nn;
      o[2] = std::fabs(nr[0] * od[0] + nr[1] * od[1] + nr[2] * od[2]);
      int ax = 0;
      for (int k = 1; k < 3; ++k) if (std::fabs(nr[k]) > std::fabs(nr[ax])) ax = k;
      o[3] = ax;
      o[4] = std::fabs(nr[ax]);
      auto* bb = static_cast<TGeoBBox*>(os);
      const double h[3] = {bb->GetDX() * 10, bb->GetDY() * 10, bb->GetDZ() * 10};
      const double* org = bb->GetOrigin();
      double u[3], edge = 1e30, area = 4;
      for (int k = 0; k < 3; ++k) {
        const double c = std::fabs(ol[k] * 10 - org[k] * 10);
        u[k] = h[k] > 0 ? c / h[k] : nan;
        if (k != ax) { edge = std::min(edge, h[k] - c); area *= h[k]; }
      }
      o[5] = SignedSafetyMM(os, ol);
      o[6] = u[0]; o[7] = u[1]; o[8] = u[2];
      o[9] = edge;
      o[10] = area;
      o[11] = std::min({h[0], h[1], h[2]});
      o[12] = std::max({h[0], h[1], h[2]});
    }
    RotFeat(own->m, o[13], o[14], o[15], o[27]);
    const double* t = own->m.GetTranslation();
    o[16] = 10 * std::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
    o[17] = 10 * std::sqrt(ol[0] * ol[0] + ol[1] * ol[1] + ol[2] * ol[2]);
    auto* cb = static_cast<TGeoBBox*>(cs);
    o[18] = 10 * std::min({cb->GetDX(), cb->GetDY(), cb->GetDZ()});
    o[19] = 10 * std::max({cb->GetDX(), cb->GetDY(), cb->GetDZ()});
    RotFeat(cur.m, o[20], o[21], o[22], o[28]);
    const double safe1 = std::fabs(SignedSafetyMM(cs, l1));
    o[23] = safe1 <= tol ? 1 : 0;
    if (safe1 <= tol) {
      double n1[3] = {0, 0, 0};
      cs->ComputeNormal(l1, ld, n1);
      const double m1 = std::sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
      if (m1 > 0) o[24] = std::fabs(n1[0] * ld[0] + n1[1] * ld[1] + n1[2] * ld[2]) / m1;
    }
    o[25] = cs->Safety(lm, true) * 10;
  }
}
}  // namespace rdtri
"""


class Geometry:
    def __init__(self, compact):
        import ROOT
        import dd4hep  # noqa: F401
        ROOT.gErrorIgnoreLevel = ROOT.kWarning
        dd4hep.Detector.getInstance().fromXML(compact)
        ROOT.gInterpreter.Declare(CPP)
        self.R = ROOT.rdtri
        n = self.R.Init()
        self.vol_name = [str(self.R.VolName(i)) for i in range(n)]
        self.vol_shape = [str(self.R.ShapeClass(i)) for i in range(n)]
        self.vol_mat = [str(self.R.MatName(i)) for i in range(n)]
        self.names = [str(x) for x in self.R.Names()]

    def extract(self, p0, p1, pprev, d):
        n = len(p0)
        c = lambda a: np.ascontiguousarray(a, dtype=np.float64).ravel()  # noqa: E731
        p0, p1, pprev, d = c(p0), c(p1), c(pprev), c(d)
        out = np.zeros(n * len(self.names))
        iv = np.zeros(n, dtype=np.int32)
        io = np.zeros(n, dtype=np.int32)
        self.R.Extract(n, p0, p1, pprev, d, out, iv, io)
        return out.reshape(n, -1), iv, io


# ── data: scans, overlay, per-step labels ──────────────────────────────────────

def load_pair(case):
    mats = {}
    A, B = Scan(os.path.join(case, "dd4hep"), mats), Scan(os.path.join(case, "step"), mats)
    common = np.intersect1d(A.rays["evt"], B.rays["evt"])  # as compare.py: drop hung rays
    for S in (A, B):
        S.rays = S.rays[np.isin(S.rays["evt"], common)]
        keep = np.isin(S.steps["evt"], common)
        S.steps, S.smat = S.steps[keep], S.smat[keep]
        S.rays["evt"] = np.searchsorted(common, S.rays["evt"])
        S.steps["evt"] = np.searchsorted(common, S.steps["evt"])
    L = 4.0 * math.sqrt(sum(h * h for h in A.meta["world_half_mm"])) + 1.0
    names = [None] * len(mats)
    for n, m in mats.items():
        names[m["id"]] = n
    return A, B, names, L, common


def assign_steps(evtA, s0A, evtI, sI):
    """Index of the A step containing each interval start (A steps sorted by evt, s0)."""
    n = len(evtA)
    evt = np.r_[evtA, evtI]
    s = np.r_[s0A, sI]
    typ = np.r_[np.zeros(n, np.int8), np.ones(len(evtI), np.int8)]
    order = np.lexsort((typ, s, evt))
    last = np.maximum.accumulate(np.where(order < n, order, -1))
    res = np.empty(len(evtI), np.int64)
    isI = order >= n
    res[order[isI] - n] = last[isI]
    return res


def step_labels(A, B, L, min_mm):
    S = A.steps
    key_ok = np.all((S["evt"][1:] > S["evt"][:-1]) | ((S["evt"][1:] == S["evt"][:-1]) & (S["s0"][1:] >= S["s0"][:-1])))
    if not key_ok:
        sys.exit("reference steps are not ordered by (evt, s0)")
    evt, s, ln, ma, mb, va, vb = overlay(A, B, len(A.rays), L)
    k = assign_steps(S["evt"].astype(np.int64), S["s0"], evt, s)
    bad = (k < 0) | (S["evt"][np.maximum(k, 0)] != evt)
    mis = (ma != mb) & (ma >= 0) & (mb >= 0) & ~bad
    ns = len(S)
    mis_len = np.bincount(k[mis], weights=ln[mis], minlength=ns)
    pos = mis_len > min_mm
    # per positive step: dominant candidate material, where the mismatch starts
    ki, si, li, mbi = k[mis], s[mis], ln[mis], mb[mis]
    o = np.lexsort((-li, ki))
    first = np.r_[True, ki[o][1:] != ki[o][:-1]]
    dom_mb = np.full(ns, -1, np.int64)
    dom_mb[ki[o][first]] = mbi[o][first]
    o2 = np.lexsort((si, ki))
    f2 = np.r_[True, ki[o2][1:] != ki[o2][:-1]]
    start = np.full(ns, np.nan)
    start[ki[o2][f2]] = si[o2][f2]
    return mis_len, pos, dom_mb, start, int(bad.sum())


# ── statistics: weighted AUC, CART, clustering ─────────────────────────────────

def wauc(score, y, w):
    """Weighted Mann-Whitney AUC (ties count 1/2); NaN scores ignored."""
    ok = np.isfinite(score)
    score, y, w = score[ok], y[ok], w[ok]
    wp, wn = w[y].sum(), w[~y].sum()
    if wp == 0 or wn == 0:
        return float("nan")
    u, inv = np.unique(score, return_inverse=True)
    gp = np.bincount(inv, weights=w * y, minlength=len(u))
    gn = np.bincount(inv, weights=w * ~y, minlength=len(u))
    below = np.cumsum(gn) - gn
    return float((gp * (below + 0.5 * gn)).sum() / (wp * wn))


def waterfill(groups, total, rng):
    """Stratified sample: min(n_s, q) rows from each group, q set so about `total` rows are taken.
    Small groups are kept whole, so a rare failure mode is not drowned by a common one.
    Returns sorted indices and weights n_s / taken_s (restore the population)."""
    sizes = np.array([len(g) for g in groups])
    q = int(sizes.max()) if sizes.sum() <= total else 0
    if not q:
        lo, hi = 1, int(sizes.max())
        while lo < hi:
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if np.minimum(sizes, mid).sum() >= total else (mid + 1, hi)
        q = lo
    idx, wt = [], []
    for g in groups:
        k = min(len(g), q)
        idx.append(g if k == len(g) else rng.choice(g, k, replace=False))
        wt.append(np.full(k, len(g) / k))
    idx, wt = np.concatenate(idx), np.concatenate(wt)
    o = np.argsort(idx)
    return idx[o], wt[o]


def wmedian(x, w):
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return float(x[o][np.searchsorted(c, 0.5 * c[-1])])


def gini_imp(wp, wn):
    t = wp + wn
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(t > 0, wp / t, 0.0)
    return t * 2 * p * (1 - p)  # weight * Gini


def fit_tree(X, y, w, depth, min_w, min_n, min_gain, _d=0):
    wp, wn = float(w[y].sum()), float(w[~y].sum())
    node = {"n": int(len(y)), "wp": wp, "wn": wn}
    if _d >= depth or wp == 0 or wn == 0 or len(y) < 2 * min_n:
        return node
    parent = float(gini_imp(wp, wn))
    best = (0.0, None, None, None, None)
    for j in range(X.shape[1]):
        x = X[:, j]
        o = np.argsort(x, kind="stable")
        xs = x[o]
        if xs[0] == xs[-1]:
            continue
        cp = np.cumsum((w * y)[o])
        cn = np.cumsum((w * ~y)[o])
        cut = np.nonzero(xs[:-1] < xs[1:])[0]  # split after position cut
        nl = cut + 1
        lp, ln_ = cp[cut], cn[cut]
        rp, rn = wp - lp, wn - ln_
        ok = (nl >= min_n) & (len(y) - nl >= min_n) & (lp + ln_ >= min_w) & (rp + rn >= min_w)
        if not ok.any():
            continue
        gain = parent - gini_imp(lp, ln_) - gini_imp(rp, rn)
        gain = np.where(ok, gain, -np.inf)
        i = int(np.argmax(gain))
        if gain[i] > best[0]:
            best = (float(gain[i]), j, 0.5 * (xs[cut[i]] + xs[cut[i] + 1]), float(xs[cut[i]]), float(xs[cut[i] + 1]))
    if best[1] is None or best[0] < min_gain:
        return node
    j, thr = best[1], best[2]
    m = X[:, j] <= thr
    node.update(feat=j, thr=thr, lo=best[3], hi=best[4], gain=best[0],
                left=fit_tree(X[m], y[m], w[m], depth, min_w, min_n, min_gain, _d + 1),
                right=fit_tree(X[~m], y[~m], w[~m], depth, min_w, min_n, min_gain, _d + 1))
    return node


def leaves(node, path=()):
    if "feat" not in node:
        yield path, node
        return
    s = (node["thr"], node["lo"], node["hi"])  # threshold and the data values either side of it
    yield from leaves(node["left"], path + ((node["feat"], "<=", s),))
    yield from leaves(node["right"], path + ((node["feat"], ">", s),))


def apply_tree(node, X):
    """Leaf id (index in leaves() order) per row."""
    out = np.zeros(len(X), np.int64)
    counter = [0]

    def rec(nd, idx):
        if "feat" not in nd:
            out[idx] = counter[0]
            counter[0] += 1
            return
        m = X[idx, nd["feat"]] <= nd["thr"]
        rec(nd["left"], idx[m])
        rec(nd["right"], idx[~m])

    rec(node, np.arange(len(X)))
    return out


def silhouette(D, lab):
    ks = np.unique(lab)
    if len(ks) < 2:
        return float("nan")
    a = np.zeros(len(lab))
    b = np.full(len(lab), np.inf)
    for k in ks:
        m = lab == k
        dk = D[:, m].sum(1)
        own = np.where(m, (dk) / np.maximum(m.sum() - 1, 1), np.inf)
        a = np.where(m, own, a)
        b = np.where(m, b, np.minimum(b, dk / m.sum()))
    s = np.where(np.isfinite(b), (b - a) / np.maximum(np.maximum(a, b), 1e-12), 0)
    return float(s.mean())


# ── feature table ──────────────────────────────────────────────────────────────

class Features:
    """Column store with metadata: kind num|log|bool|onehot, group what|where."""

    def __init__(self):
        self.cols, self.meta = {}, {}

    def add(self, name, values, kind="num", group="what", label=None, src=None):
        self.cols[name] = np.asarray(values, dtype=np.float64)
        self.meta[name] = {"kind": kind, "group": group, "label": label or name, "src": src}

    def onehot(self, base, cats, group, top=12, prefix=None):
        cats = np.asarray(cats, dtype=object)
        vals, cnt = np.unique(cats, return_counts=True)
        for v in vals[np.argsort(-cnt)][:top]:
            self.add(f"{prefix or base}={v}", cats == v, "onehot", group, label=f"{base} = {v}", src=(base, v))

    def matrix(self, groups):
        names = [n for n in self.cols if self.meta[n]["group"] in groups]
        X = np.column_stack([self.cols[n] for n in names])
        return np.nan_to_num(X, nan=-1e3, posinf=1e30, neginf=-1e30), names


def log10f(x, floor=1e-18):
    return np.log10(np.maximum(np.abs(x), floor))


def fmt_val(v):
    if v == 0 or (1e-3 <= abs(v) < 1e5):
        return f"{v:.4g}"
    return f"{v:.2e}"


def fmt_between(lo, hi):
    """Shortest decimal d with lo < d < hi (a split threshold that keeps the data on the right side)."""
    for n in range(2, 16):
        for spec in (f".{n}g", f".{n - 1}e"):
            t = format(0.5 * (lo + hi), spec)
            if lo < float(t) < hi:
                v = float(t)
                return fmt_val(v) if float(fmt_val(v)) == v else t
    return format(0.5 * (lo + hi), ".17g")


def cond_text(F, name, lo, hi):
    m = F.meta[name]
    if m["kind"] in ("onehot", "bool"):
        neg = hi is not None and hi[0] < 1  # "<= 0.5" means false
        if m["kind"] == "onehot":
            base, v = m["src"]
            return f"{base} {'≠' if neg else '='} {v}"
        return f"not {m['label']}" if neg else m["label"]
    tr = (lambda x: 10 ** x) if m["kind"] == "log" else (lambda x: x)
    lab = m["label"]
    if lo is not None and lo[0] <= -500:  # NaN sentinel (-1e3) excluded: "defined"
        lo = None
        if hi is None:
            return f"{lab} defined"
    if hi is not None and hi[0] <= -500:
        return f"{lab} undefined"
    fv = lambda s: fmt_between(tr(s[1]), tr(s[2]))  # noqa: E731
    if lo is not None and hi is not None:
        return f"{fv(lo)} < {lab} ≤ {fv(hi)}"
    if lo is not None:
        return f"{lab} > {fv(lo)}"
    return f"{lab} ≤ {fv(hi)}"


def rule_text(F, names, path):
    bounds = {}
    for j, op, sp in path:
        lo, hi = bounds.get(j, (None, None))
        if op == ">":
            lo = sp if lo is None or sp[0] > lo[0] else lo
        else:
            hi = sp if hi is None or sp[0] < hi[0] else hi
        bounds[j] = (lo, hi)
    return " AND ".join(cond_text(F, names[j], *bounds[j]) for j in bounds) or "(all steps)"


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", help="case directory, e.g. out/minisid_ip")
    ap.add_argument("--out", default=None, help="default CASE/cmp/triage")
    ap.add_argument("--compact", default=None, help="default: compact file recorded in the scan metadata")
    ap.add_argument("--min-mm", type=float, default=1e-6, help="a step is positive above this mismatch length")
    ap.add_argument("--neg-ratio", type=float, default=20, help="negatives sampled per positive")
    ap.add_argument("--max-pos", type=int, default=50000)
    ap.add_argument("--max-neg", type=int, default=400000)
    ap.add_argument("--depth", type=int, default=4)
    ap.add_argument("--max-k", type=int, default=8)
    ap.add_argument("--min-pos", type=int, default=30, help="below this many positive steps: list them, no rules/clusters")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    case = args.case.rstrip("/")
    out = args.out or os.path.join(case, "cmp", "triage")
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    A, B, matn, L, common = load_pair(case)
    compact = args.compact or A.meta["compact"]
    S, R = A.steps, A.rays
    mis_len, pos, dom_mb, mstart, nbad = step_labels(A, B, L, args.min_mm)
    ns, npos_all = len(S), int(pos.sum())
    report = {"case": case, "compact": compact, "min_mm": args.min_mm, "steps": ns, "positive_steps": npos_all,
              "mismatch_length_mm": float(mis_len.sum()),
              "mismatch_length_in_positive_steps_mm": float(mis_len[pos].sum()),
              "intervals_not_assigned": nbad}
    print(f"[triage] {case}: {ns} reference steps, {npos_all} with mismatch > {args.min_mm:g} mm "
          f"({mis_len[pos].sum():.4g} of {mis_len.sum():.4g} mm)")
    if npos_all == 0:
        write_empty(out, report)
        return

    print(f"[triage] loading {compact} ...")
    G = Geometry(compact)
    nm = lambda i: matn[i] if i >= 0 else "<none>"  # noqa: E731

    # strata of positive steps: (reference material -> candidate material, solid class)
    pall = np.nonzero(pos)[0]
    vs = S["vol"][pall]
    skey = np.array([f"{nm(a)} → {nm(b)}, {G.vol_shape[v] if v >= 0 else 'unmapped volume'}"
                     for a, b, v in zip(A.smat[pall], dom_mb[pall], vs)], dtype=object)
    keys, inv = np.unique(skey, return_inverse=True)
    tot_len = mis_len[pall].sum()
    strata = []
    for k, key in enumerate(keys):
        m = inv == k
        strata.append({"stratum": str(key), "steps": int(m.sum()), "mismatch_length_mm": float(mis_len[pall][m].sum()),
                       "mismatch_length_share": float(mis_len[pall][m].sum() / tot_len)})
    strata.sort(key=lambda r: -r["mismatch_length_share"])
    report["strata"] = strata

    # sample: positives stratified (rare strata kept whole), negatives uniform
    pidx, wpos = waterfill([pall[inv == k] for k in range(len(keys))], args.max_pos, rng)
    pkey = dict(zip(pall, skey))
    nidx = np.nonzero(~pos)[0]
    nneg = int(min(len(nidx), args.max_neg, max(args.neg_ratio * len(pidx), 2000)))
    nsel = np.sort(rng.choice(nidx, nneg, replace=False))
    idx = np.r_[pidx, nsel]
    y = np.r_[np.ones(len(pidx), bool), np.zeros(len(nsel), bool)]
    w = np.r_[wpos, np.full(len(nsel), len(nidx) / max(len(nsel), 1))]
    stratum = np.array([pkey[i] for i in pidx] + [""] * len(nsel), dtype=object)
    report["sampled"] = {"positive": int(len(pidx)), "negative": int(len(nsel)), "strata": int(len(keys)),
                         "weight_positive_range": [float(wpos.min()), float(wpos.max())],
                         "weight_negative": float(w[~y][0]) if len(nsel) else 0}

    # geometry of each sampled step
    st = S[idx]
    r = R[np.searchsorted(R["evt"], st["evt"])]
    o = np.column_stack([r["ox"], r["oy"], r["oz"]])
    d = np.column_stack([r["dx"], r["dy"], r["dz"]])
    d = d / np.linalg.norm(d, axis=1, keepdims=True)
    p0 = o + d * st["s0"][:, None]
    p1 = o + d * (st["s0"] + st["len"])[:, None]
    has_prev = (idx > 0) & (S["evt"][np.maximum(idx - 1, 0)] == st["evt"])
    sp = S[np.maximum(idx - 1, 0)]
    pprev = o + d * (sp["s0"] + 0.5 * sp["len"])[:, None]
    pprev[~has_prev] = np.nan
    prev_mat = np.where(has_prev, A.smat[np.maximum(idx - 1, 0)], -1)
    first_in_ray = np.searchsorted(S["evt"], st["evt"])
    kstep = idx - first_in_ray

    print(f"[triage] extracting features for {len(idx)} steps ({len(keys)} strata of positives) ...")
    Fx, ivol, iown = G.extract(p0, p1, pprev, d)
    fx = {n: Fx[:, i] for i, n in enumerate(G.names)}
    mapped = st["vol"] >= 0  # the scanner writes -1 for volumes it does not map (e.g. some Air volumes)
    report["scan_volume_unmapped_fraction"] = float(np.average(~mapped, weights=w))
    report["relocation_agreement"] = float(np.average(ivol[mapped] == st["vol"][mapped], weights=w[mapped]))
    print(f"[triage] TGeo relocation agrees with the scan's volume id for {100 * report['relocation_agreement']:.4f} % "
          f"of mapped steps ({100 * report['scan_volume_unmapped_fraction']:.1f} % unmapped in the scan, located by TGeo)")

    mat_a = np.array([nm(i) for i in A.smat[idx]], dtype=object)
    vname = np.array([G.vol_name[v] if v >= 0 else "<?>" for v in ivol], dtype=object)
    shape = np.array([G.vol_shape[v] if v >= 0 else "<?>" for v in ivol], dtype=object)
    oshape = np.array([G.vol_shape[v] if v >= 0 else "<?>" for v in iown], dtype=object)
    oname = np.array([G.vol_name[v] if v >= 0 else "<?>" for v in iown], dtype=object)
    pmat = np.array([nm(i) for i in prev_mat], dtype=object)
    axis = np.array([("x", "y", "z")[int(a)] if np.isfinite(a) else "?" for a in fx["normal_axis"]], dtype=object)
    own_kind = np.array([{1: "own face", 0: "daughter face"}.get(int(v), "unknown") if np.isfinite(v) else "unknown"
                         for v in fx["own_entry"]], dtype=object)

    F = Features()
    # what: mechanism
    F.onehot("material", mat_a, "what")
    F.onehot("previous material", pmat, "what", prefix="prev_material")
    F.onehot("solid", shape, "what")
    F.onehot("entry face solid", oshape, "what", prefix="face_solid")
    F.onehot("entry surface", own_kind, "what", prefix="entry_surface")
    F.onehot("entry normal axis (local)", axis, "what", prefix="normal_axis")
    F.add("cos_in", fx["cos_in"], label="|cos incidence|")
    F.add("normal_align", fx["normal_align"], label="max |local normal component|")
    F.add("log_entry_offset", log10f(fx["entry_offset_mm"]), "log", label="|entry offset from face| [mm]")
    F.add("entry_inside", fx["entry_offset_mm"] > 0, "bool", label="entry point inside face owner")
    for k in "xyz":
        F.add(f"u_{k}", fx[f"u_{k}"], label=f"|local {k}| / half-length")
    F.add("log_edge_dist", log10f(fx["edge_dist_mm"]), "log", label="entry distance to face edge [mm]")
    F.add("log_face_area", log10f(fx["face_area_mm2"]), "log", label="entry face area (bbox) [mm²]")
    F.add("log_own_hmin", log10f(fx["own_hmin_mm"]), "log", label="face owner min half-length [mm]")
    F.add("log_own_hmax", log10f(fx["own_hmax_mm"]), "log", label="face owner max half-length [mm]")
    F.add("own_rot_deg", fx["own_rot_deg"], label="face owner rotation angle [deg]")
    F.add("own_phimod90", np.round(fx["own_phimod90_deg"], 6), label="face owner φz off 90° grid [deg]")
    F.add("own_aligned", fx["own_align_err"] < 1e-9, "bool", label="face owner rotation axis-aligned")
    F.add("own_phiz", np.round(fx["own_phiz_deg"], 6), label="face owner φz [deg]")
    F.add("log_own_trans", log10f(fx["own_trans_mm"], 1e-6), "log", label="face owner |translation| [mm]")
    F.add("log_own_local", log10f(fx["own_local_mag_mm"], 1e-6), "log", label="|entry point, owner frame| [mm]")
    F.add("log_cur_hmin", log10f(fx["cur_hmin_mm"]), "log", label="solid min half-length [mm]")
    F.add("log_cur_hmax", log10f(fx["cur_hmax_mm"]), "log", label="solid max half-length [mm]")
    F.add("cur_rot_deg", fx["cur_rot_deg"], label="placement rotation angle [deg]")
    F.add("cur_phimod90", np.round(fx["cur_phimod90_deg"], 6), label="placement φz off 90° grid [deg]")
    F.add("cur_aligned", fx["cur_align_err"] < 1e-9, "bool", label="placement rotation axis-aligned")
    F.add("cur_phiz", np.round(fx["cur_phiz_deg"], 6), label="placement φz [deg]")
    F.add("level", fx["level"], label="geometry depth")
    F.add("exit_own", fx["exit_own"], "bool", label="exit through own face")
    F.add("cos_out", fx["cos_out"], label="|cos exit|")
    F.add("log_mid_safety", log10f(fx["mid_safety_mm"]), "log", label="midpoint safety [mm]")
    F.add("log_len", log10f(st["len"]), "log", label="step length [mm]")
    F.add("log_thick", log10f(st["len"] * fx["cos_in"]), "log", label="step length × |cos in| [mm]")
    # where: location
    rr = np.hypot(p0[:, 0], p0[:, 1])
    F.add("r", rr, group="where", label="entry r [mm]")
    F.add("z", p0[:, 2], group="where", label="entry z [mm]")
    F.add("abs_z", np.abs(p0[:, 2]), group="where", label="entry |z| [mm]")
    F.add("phi", np.degrees(np.arctan2(p0[:, 1], p0[:, 0])), group="where", label="entry φ [deg]")
    th = np.arccos(np.clip(d[:, 2], -1, 1))
    F.add("eta", -np.log(np.tan(np.clip(th, 1e-9, np.pi - 1e-9) / 2)), group="where", label="ray η")
    F.add("s0", st["s0"], group="where", label="path length at entry [mm]")
    F.add("step_index", kstep, group="where", label="step index on ray")
    F.add("log_gmag", log10f(np.linalg.norm(p0, axis=1), 1e-6), "log", "where", label="|entry point, global| [mm]")
    F.onehot("volume", vname, "where", top=25, prefix="volume")

    # per-step table
    write_table(out, F, idx, st, common, y, w, mis_len[idx], np.where(y, [nm(i) for i in dom_mb[idx]], "-"),
                mat_a, vname, oname, shape, axis, own_kind)

    # 1. single-feature AUC
    auc_rows = []
    for n in F.cols:
        a = wauc(F.cols[n], y, w)
        if np.isfinite(a):
            auc_rows.append({"feature": n, "label": F.meta[n]["label"], "group": F.meta[n]["group"], "auc": a,
                             "separation": abs(a - 0.5) * 2, "higher_means": "mismatch" if a > 0.5 else "match"})
    auc_rows.sort(key=lambda x: -x["separation"])
    report["single_feature"] = auc_rows[:25]

    small = len(pidx) < args.min_pos
    if small:
        report["listed_steps"] = [
            {"evt": int(common[st["evt"][i]]), "s0_mm": float(st["s0"][i]), "len_mm": float(st["len"][i]),
             "volume": vname[i], "material": mat_a[i], "candidate_material": nm(dom_mb[idx[i]]),
             "mismatch_mm": float(mis_len[idx[i]]), "entry_surface": own_kind[i], "solid": shape[i],
             "cos_in": float(fx["cos_in"][i]), "edge_dist_mm": float(fx["edge_dist_mm"][i]),
             "entry": [float(v) for v in p0[i]],
             "show_ray": f"python3 scripts/show_ray.py {case}/dd4hep {case}/step {int(common[st['evt'][i]])} "
                         f"{st['s0'][i] - 1:.3f} {st['s0'][i] + st['len'][i] + 1:.3f} --report {case}/export.json"}
            for i in np.nonzero(y)[0]]
        report["trees"], report["clusters"] = {}, []
        json.dump(report, open(os.path.join(out, "triage.json"), "w"), indent=1, default=float)
        write_md(out, report)
        print(f"[triage] only {len(pidx)} positive steps (< --min-pos {args.min_pos}): listed, no rules or clusters; "
              f"wrote {out}/triage.md")
        return

    # 2. trees
    train = (common[st["evt"]] % 2) == 0
    report["trees"] = {}
    Xall = F.matrix(("what", "where"))
    for tag, groups in (("mechanism", ("what",)), ("mechanism+location", ("what", "where"))):
        X, names = F.matrix(groups)
        report["trees"][tag] = tree_report(F, X, names, y, w, mis_len[idx], train, args.depth,
                                           drill=Xall if tag == "mechanism" else None)

    # 3. clusters of positive steps
    report["clusters"], lab, med = clusters(F, y, w, mis_len[idx], st, fx, args.max_k, rng,
                                            np.where(y, [nm(i) for i in dom_mb[idx]], "-"), mat_a, mstart[idx], stratum)
    for c in report["clusters"]:
        i = c.pop("_medoid")
        e, s0, ln = int(common[st["evt"][i]]), float(st["s0"][i]), float(st["len"][i])
        c["medoid"] = {"evt": e, "s0_mm": s0, "len_mm": ln, "volume": vname[i], "material": mat_a[i],
                       "show_ray": f"python3 scripts/show_ray.py {case}/dd4hep {case}/step {e} {s0 - 1:.3f} {s0 + ln + 1:.3f} "
                                   f"--report {case}/export.json",
                       "debug_ray": debug_cmd(A, B, e, s0 - 0.01, s0 + ln + 0.01)}

    json.dump(report, open(os.path.join(out, "triage.json"), "w"), indent=1, default=float)
    write_md(out, report)
    if not args.no_plots:
        plots(out, report, F, y, lab, p0)
    print(f"[triage] wrote {out}/triage.md")


def debug_cmd(A, B, evt, smin, smax):
    m, mb = A.meta, B.meta
    args = [f"--compact {m['compact']}", "--world step", f"--step {mb.get('step', '')}", f"--ray-mode {m['ray_mode']}",
            f"--seed {m['seed']}", "--box " + ",".join(f"{v:g}" for v in m["box_mm"])]
    if m["ray_mode"] == "aim":
        args.append("--target " + ",".join(f"{v:g}" for v in m["target_mm"]))
    return "${SCANNER:-build/simcadval_scan} " + " ".join(args) + f" --debug-ray {evt},{max(smin, 0):.3f},{smax:.3f}"


def tree_report(F, X, names, y, w, mlen, train, depth, drill=None):
    def bal(m):
        wp, wn = w[m & y].sum(), w[m & ~y].sum()
        return np.where(y, 0.5 / max(wp, 1e-300), 0.5 / max(wn, 1e-300)) * w

    wb = bal(train)
    tr = train
    tree = fit_tree(X[tr], y[tr], wb[tr], depth, min_w=0.01, min_n=10, min_gain=0.002)
    leaf = apply_tree(tree, X)
    out = {"depth": depth, "rules": []}
    # held-out ranking quality: leaf positive rate as score
    lv = list(leaves(tree))
    score_tr = np.array([nd["wp"] / max(nd["wp"] + nd["wn"], 1e-300) for _, nd in lv])
    te = ~train
    out["test_auc"] = wauc(score_tr[leaf[te]].astype(float), y[te], w[te]) if te.any() else float("nan")
    tot = {}
    for part, m in (("train", tr), ("test", te)):
        tot[part] = (w[m & y].sum(), (w * mlen)[m & y].sum())
    for li, (path, nd) in enumerate(lv):
        if score_tr[li] <= 0.5:
            continue  # leaf predicts "match"
        row = {"rule": rule_text(F, names, path)}
        for part, m in (("train", tr), ("test", te)):
            mm = m & (leaf == li)
            wp, wn = w[mm & y].sum(), w[mm & ~y].sum()
            row[part] = {"steps_pos": int((mm & y).sum()), "steps_neg_est": float(wn),
                         "precision": float(wp / (wp + wn)) if wp + wn > 0 else float("nan"),
                         "recall": float(wp / tot[part][0]) if tot[part][0] > 0 else float("nan"),
                         "mismatch_length_share": float((w * mlen)[mm & y].sum() / tot[part][1]) if tot[part][1] > 0 else float("nan")}
        out["rules"].append(row)
    if drill is not None:
        Xd, nd_ = drill
        for li, row in zip([li for li in range(len(lv)) if score_tr[li] > 0.5], out["rules"]):
            if np.isfinite(row["test"]["precision"]) and row["test"]["precision"] < 0.8:
                row["drilldown"] = drilldown(F, Xd, nd_, y, w, mlen, train, leaf == li)
    out["rules"].sort(key=lambda r: -r["train"]["recall"])
    out["covered_recall_test"] = float(sum(r["test"]["recall"] for r in out["rules"] if np.isfinite(r["test"]["recall"])))
    # the first split, as the single most informative condition
    if "feat" in tree:
        out["root_split"] = rule_text(F, names, ((tree["feat"], ">", (tree["thr"], tree["lo"], tree["hi"])),))
    return out


def drilldown(F, X, names, y, w, mlen, train, inrule, depth=3):
    """Inside one rule's population: which conditions separate failing from passing steps?"""
    m = inrule
    if (m & y).sum() < 100 or (m & ~y).sum() < 100:
        return None
    auc = []
    for j, n in enumerate(names):
        a = wauc(X[m, j], y[m], w[m])
        if np.isfinite(a) and abs(a - 0.5) > 0.02:
            auc.append({"label": F.meta[n]["label"], "group": F.meta[n]["group"], "auc": a})
    auc.sort(key=lambda r: -abs(r["auc"] - 0.5))
    tr = m & train
    wp, wn = w[tr & y].sum(), w[tr & ~y].sum()
    wb = np.where(y, 0.5 / wp, 0.5 / wn) * w
    sub = fit_tree(X[tr], y[tr], wb[tr], depth, min_w=0.02, min_n=50, min_gain=0.002)
    leaf = apply_tree(sub, X)
    te = m & ~train
    tot = (w * mlen)[te & y].sum()
    rows = []
    for li, (path, _) in enumerate(leaves(sub)):
        mm = te & (leaf == li)
        p, q = w[mm & y].sum(), w[mm & ~y].sum()
        if p + q == 0:
            continue
        rows.append({"rule": rule_text(F, names, path), "test_fail_rate": float(p / (p + q)),
                     "test_share_of_rule_population": float((p + q) / (w[te].sum())),
                     "test_mismatch_length_share": float((w * mlen)[mm & y].sum() / tot) if tot > 0 else float("nan")})
    rows.sort(key=lambda r: -r["test_fail_rate"])
    base = float(w[te & y].sum() / w[te].sum()) if w[te].sum() > 0 else float("nan")
    return {"fail_rate": base, "single_feature": auc[:8], "rules": rows}


def clusters(F, y, w, mlen, st, fx, max_k, rng, mat_b, mat_a, mstart, stratum):
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist, squareform

    P = np.nonzero(y)[0]
    # mechanism features + mismatch shape (only defined for positives). The entry offset and its sign
    # are left out: at ~1e-12 mm they are dominated by the rounding of origin + s0*dir and would split
    # one population into "inside"/"outside" halves.
    X, names = F.matrix(("what",))
    noise = [j for j, n in enumerate(names) if n in ("entry_inside", "log_entry_offset")]
    X = np.delete(X, noise, axis=1)
    names = [n for j, n in enumerate(names) if j not in noise]
    frac = np.clip(mlen / np.maximum(st["len"], 1e-300), 0, 1)
    sfrac = np.clip((mstart - st["s0"]) / np.maximum(st["len"], 1e-300), 0, 1)
    extra = {"mis_frac": frac, "mis_start_frac": np.nan_to_num(sfrac)}
    mb_cats = np.asarray(mat_b, dtype=object)
    vals, cnt = np.unique(mb_cats[P], return_counts=True)
    for v in vals[np.argsort(-cnt)][:8]:
        extra[f"mat_b={v}"] = (mb_cats == v).astype(float)
    Xc = np.column_stack([X[P]] + [extra[k][P] for k in extra])
    cn = names + list(extra)
    kinds = [F.meta[n]["kind"] for n in names] + ["num", "num"] + ["onehot"] * (len(extra) - 2)
    keep = np.array([np.ptp(Xc[:, j]) > 0 for j in range(Xc.shape[1])])
    Xc, cn, kinds = Xc[:, keep], [n for n, k in zip(cn, keep) if k], [k for k, kk in zip(kinds, keep) if kk]
    Z = np.empty_like(Xc)
    for j in range(Xc.shape[1]):
        x = Xc[:, j]
        if kinds[j] in ("onehot", "bool"):
            Z[:, j] = 1.5 * x
        else:
            med = np.median(x)
            q = np.subtract(*np.percentile(x, [75, 25]))
            sc = q if q > 0 else (x.std() or 1.0)
            Z[:, j] = np.clip((x - med) / sc, -3, 3)
    # stratified subsample for Ward (needs all pairwise distances): rare strata kept whole
    sk, sinv = np.unique(stratum[P], return_inverse=True)
    si, _ = waterfill([np.nonzero(sinv == k)[0] for k in range(len(sk))], 2500, rng)
    D = squareform(pdist(Z[si]))
    lk = linkage(Z[si], "ward")
    best_k, best_s, sil = 1, float("nan"), {}
    for k in range(2, max_k + 1):
        lab = fcluster(lk, k, "maxclust")
        if len(np.unique(lab)) < k:
            continue
        s = silhouette(D, lab)
        sil[k] = s
        if not np.isfinite(best_s) or s > best_s:
            best_k, best_s = k, s
    if not np.isfinite(best_s) or best_s < 0.25:
        best_k = 1
    lab_s = fcluster(lk, best_k, "maxclust") - 1 if best_k > 1 else np.zeros(len(si), int)
    cent = np.array([Z[si][lab_s == c].mean(0) for c in range(best_k)])
    labP = np.argmin(((Z[:, None, :] - cent[None]) ** 2).sum(-1), axis=1)
    lab_all = np.full(len(y), -1)
    lab_all[P] = labP
    wl = (w * mlen)[P]
    wP = w[P]
    rows = []
    for c in range(best_k):
        m = labP == c
        row = {"id": c, "steps": int(round(wP[m].sum())), "sampled_steps": int(m.sum()),
               "mismatch_length_share": float(wl[m].sum() / wl.sum()),
               "silhouette_k": best_k, "silhouette": best_s if best_k > 1 else None}
        # categorical signature
        sig = {}
        for base in ("material", "previous material", "solid", "entry face solid", "entry surface",
                     "entry normal axis (local)"):
            cols = [n for n in F.cols if F.meta[n]["kind"] == "onehot" and F.meta[n]["src"][0] == base]
            if not cols:
                continue
            shares = {F.meta[n]["src"][1]: float(np.average(F.cols[n][P][m], weights=wP[m])) for n in cols}
            top = sorted(shares.items(), key=lambda kv: -kv[1])[:2]
            sig[base] = [f"{k} {100 * v:.0f}%" for k, v in top if v > 0.05]
        pairs, pinv = np.unique(np.char.add(np.char.add(mat_a[P][m].astype(str), " → "), mat_b[P][m].astype(str)),
                                return_inverse=True)
        pc = np.bincount(pinv, weights=wP[m])
        sig["mismatch"] = [f"{p} {100 * n / pc.sum():.0f}%" for p, n in sorted(zip(pairs, pc), key=lambda t: -t[1])[:3]]
        row["signature"] = sig
        # numeric differences vs other positives
        diffs = []
        for j, n in enumerate(cn):
            if kinds[j] in ("onehot",) or n.startswith("mat_b="):
                continue
            x = Xc[:, j]
            q = np.subtract(*np.percentile(x, [75, 25])) or (x.std() or 1.0)
            a = wmedian(x[m], wP[m])
            b = wmedian(x[~m], wP[~m]) if (~m).any() else a
            diffs.append((abs(a - b) / q, n, a, b))
        diffs.sort(reverse=True)
        row["medians"] = {}
        show = [t for t in diffs[:6]] if best_k > 1 else []
        for key in ("cos_in", "log_edge_dist", "log_face_area", "own_phimod90", "mis_frac",
                    "mis_start_frac"):
            if key in cn:
                show.append((0, key, wmedian(Xc[m, cn.index(key)], wP[m]), None))
        seen = set()
        for _, n, a, b in show:
            if n in seen:
                continue
            seen.add(n)
            meta = F.meta.get(n, {"kind": "num", "label": {"mis_frac": "fraction of step mismatched",
                                                            "mis_start_frac": "mismatch start, fraction of step"}.get(n, n)})
            tr = (lambda v: 10 ** v) if meta["kind"] == "log" else (lambda v: v)
            row["medians"][meta["label"]] = {"cluster": float(tr(a)), "other_positives": None if b is None else float(tr(b))}
        if best_k > 1:
            row["distinguishing"] = [F.meta.get(n, {"label": n})["label"] for _, n, _, _ in diffs[:4]]
        mi = np.nonzero(m)[0]
        dd = ((Z[mi] - cent[c]) ** 2).sum(1)
        row["_medoid"] = int(P[mi[np.argmin(dd)]])
        rows.append(row)
    rows.sort(key=lambda r: -r["mismatch_length_share"])
    remap = {r["id"]: i for i, r in enumerate(rows)}
    for r in rows:
        r["id"] = remap[r["id"]]
    lab_all = np.where(lab_all >= 0, np.vectorize(lambda v: remap.get(v, -1))(lab_all), -1)
    meta = {"k": best_k, "silhouette_by_k": sil}
    if rows:
        rows[0]["_k_selection"] = meta
    return rows, lab_all, None


def write_table(out, F, idx, st, common, y, w, mlen, mat_b, mat_a, vname, oname, shape, axis, own_kind):
    cols = {"evt": common[st["evt"]], "step_global": idx, "s0_mm": st["s0"], "len_mm": st["len"], "label": y.astype(int),
            "weight": w, "mismatch_mm": mlen, "material": mat_a, "material_b": mat_b, "volume": vname,
            "face_owner": oname, "solid": shape, "normal_axis": axis, "entry_surface": own_kind}
    for n, v in F.cols.items():
        if F.meta[n]["kind"] != "onehot":
            cols[n] = v
    try:
        import pandas as pd
        pd.DataFrame(cols).to_csv(os.path.join(out, "features.csv.gz"), index=False, float_format="%.9g")
    except ImportError:
        import gzip
        keys = list(cols)
        with gzip.open(os.path.join(out, "features.csv.gz"), "wt") as f:
            f.write(",".join(keys) + "\n")
            for i in range(len(y)):
                f.write(",".join(str(cols[k][i]) for k in keys) + "\n")


def write_empty(out, report):
    report["note"] = "no positive steps: nothing to triage"
    json.dump(report, open(os.path.join(out, "triage.json"), "w"), indent=1)
    with open(os.path.join(out, "triage.md"), "w") as f:
        f.write(f"# Mismatch triage: {report['case']}\n\nNo reference step has more than {report['min_mm']:g} mm of "
                f"material mismatch ({report['steps']} steps; total mismatch {report['mismatch_length_mm']:.3g} mm). "
                "Nothing to triage.\n")
    print(f"[triage] no positive steps; wrote {out}/triage.md")


def pct(x):
    return "–" if x is None or not np.isfinite(x) else f"{100 * x:.1f} %"


def write_md(out, R):
    L = []
    s = R["sampled"]
    L.append(f"# Mismatch triage: {R['case']}\n")
    L.append(f"Reference steps: {R['steps']:,}. Steps with more than {R['min_mm']:g} mm of material mismatch: "
             f"**{R['positive_steps']:,}**, holding {R['mismatch_length_in_positive_steps_mm']:.4g} mm of the "
             f"{R['mismatch_length_mm']:.4g} mm total mismatch. Sampled {s['positive']:,} positive steps "
             f"(stratified over {s['strata']} material/solid strata, rare strata kept whole) and "
             f"{s['negative']:,} negative steps; weights restore the population. TGeo relocation agrees with the "
             f"scan's volume id for {100 * R['relocation_agreement']:.4f} % of the steps the scan maps "
             f"({100 * R['scan_volume_unmapped_fraction']:.1f} % are unmapped there and located by TGeo only).\n")
    L.append("These are correlations in the reference (DD4hep/TGeo) frame. Confirm them by replaying "
             "rays with `show_ray.py` / `--debug-ray`.\n")
    L.append("## Mismatch by stratum\n")
    L.append("All positive steps, grouped by reference material → candidate material and reference solid. "
             "Small strata are listed too: a rare failure mode shows up here even when it is too small to "
             "move the rules.\n")
    L.append("| stratum | steps | mismatch length | share of length |")
    L.append("|---|---|---|---|")
    for r in R["strata"][:20]:
        L.append(f"| {r['stratum']} | {r['steps']:,} | {r['mismatch_length_mm']:.4g} mm | {pct(r['mismatch_length_share'])} |")
    if len(R["strata"]) > 20:
        rest = R["strata"][20:]
        L.append(f"| ({len(rest)} more) | {sum(r['steps'] for r in rest):,} | "
                 f"{sum(r['mismatch_length_mm'] for r in rest):.4g} mm | "
                 f"{pct(sum(r['mismatch_length_share'] for r in rest))} |")
    L.append("")
    if R.get("listed_steps"):
        L.append(f"\n## All {len(R['listed_steps'])} positive steps (too few for rules or clusters)\n")
        L.append("| ray | s0 [mm] | length [mm] | volume | material → candidate | mismatch [mm] | entry surface | "
                 "|cos in| | edge dist [mm] | entry (x, y, z) [mm] |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for r in R["listed_steps"]:
            L.append(f"| {r['evt']} | {r['s0_mm']:.3f} | {r['len_mm']:.4g} | {r['volume']} | {r['material']} → "
                     f"{r['candidate_material']} | {r['mismatch_mm']:.3g} | {r['entry_surface']} | {r['cos_in']:.3f} | "
                     f"{fmt_val(r['edge_dist_mm'])} | ({', '.join(f'{v:.1f}' for v in r['entry'])}) |")
        L.append(f"\nReplay the first one: `{R['listed_steps'][0]['show_ray']}`\n")
    L.append("## Strongest single features\n")
    L.append("AUC of the feature as a score for \"this step mismatches\" (0.5 = no information, 1 or 0 = perfect).\n")
    L.append("| feature | group | AUC | higher value means |")
    L.append("|---|---|---|---|")
    for r in R["single_feature"][:15]:
        L.append(f"| {r['label']} | {r['group']} | {r['auc']:.3f} | {r['higher_means']} |")
    for tag, T in R["trees"].items():
        L.append(f"\n## Rules: {tag} features (depth {T['depth']}, trained on even rays, scored on odd rays)\n")
        L.append(f"Held-out AUC of the tree: {T['test_auc']:.3f}. Held-out recall covered by the rules below: "
                 f"{pct(T['covered_recall_test'])}." + (f" Root split: `{T['root_split']}`." if 'root_split' in T else "") + "\n")
        if not T["rules"]:
            L.append("No leaf predicts a mismatch.\n")
            continue
        L.append("| rule | test precision | test recall | test mismatch-length share | train recall |")
        L.append("|---|---|---|---|---|")
        for r in T["rules"]:
            L.append(f"| {r['rule']} | {pct(r['test']['precision'])} | {pct(r['test']['recall'])} | "
                     f"{pct(r['test']['mismatch_length_share'])} | {pct(r['train']['recall'])} |")
        L.append("\nPrecision is population-weighted: the share of all steps matching the rule that mismatch.")
        for i, r in enumerate(T["rules"]):
            D = r.get("drilldown")
            if not D:
                continue
            L.append(f"\n### Inside rule {i + 1}: what separates the {pct(D['fail_rate'])} that fail from the rest?\n")
            L.append(f"Rule: `{r['rule']}`. Depth-3 tree on all features, restricted to steps matching the rule "
                     "(trained on even rays, scored on odd rays).\n")
            L.append("| sub-rule | test fail rate | share of rule's steps | share of rule's mismatch length |")
            L.append("|---|---|---|---|")
            for s in D["rules"]:
                L.append(f"| {s['rule']} | {pct(s['test_fail_rate'])} | {pct(s['test_share_of_rule_population'])} | "
                         f"{pct(s['test_mismatch_length_share'])} |")
            if D["single_feature"]:
                L.append("\nStrongest single features inside the rule: " + "; ".join(
                    f"{a['label']} (AUC {a['auc']:.2f})" for a in D["single_feature"][:6]) + ".")
    C = R["clusters"]
    ksel = C[0].get("_k_selection", {}) if C else {}
    if C:
        L.append(f"\n## Failure modes (clusters of positive steps): k = {ksel.get('k', len(C))}\n")
    sil = ksel.get("silhouette_by_k", {})
    if C and sil:
        L.append("Silhouette by k: " + ", ".join(f"{k}: {v:.2f}" for k, v in sil.items()) +
                 ". k = 1 is kept when no k reaches 0.25.\n")
    for c in C:
        L.append(f"### Mode {c['id']}: ~{c['steps']:,} steps ({c['sampled_steps']:,} sampled), "
                 f"{pct(c['mismatch_length_share'])} of mismatch length\n")
        for k, v in c["signature"].items():
            if v:
                L.append(f"- {k}: {', '.join(v)}")
        if c.get("distinguishing"):
            L.append(f"- differs most from other modes in: {', '.join(c['distinguishing'])}")
        L.append("\n| median of | this mode | other positives |")
        L.append("|---|---|---|")
        for k, v in c["medians"].items():
            L.append(f"| {k} | {fmt_val(v['cluster'])} | {'–' if v['other_positives'] is None else fmt_val(v['other_positives'])} |")
        m = c["medoid"]
        L.append(f"\nRepresentative step: ray {m['evt']}, s0 = {m['s0_mm']:.3f} mm, length {m['len_mm']:.4g} mm, "
                 f"in {m['volume']} ({m['material']}).\n")
        L.append(f"```bash\n{m['show_ray']}\n{m['debug_ray']}\n```\n")
    L.append("## Files\n")
    L.append("- `features.csv.gz`: one row per sampled step (label, weight, mismatch length, all features)")
    L.append("- `triage.json`: everything above, machine-readable")
    if not R.get("listed_steps"):
        L.append("- `feature_auc.png`, `clusters.png`")
    open(os.path.join(out, "triage.md"), "w").write("\n".join(L) + "\n")


def plots(out, R, F, y, lab, p0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, ink2, grid, blue, ctx = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#c9c8c2"
    plt.rcParams.update({"axes.edgecolor": ink2, "axes.labelcolor": ink, "xtick.color": ink2, "ytick.color": ink2,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9})
    rows = R["single_feature"][:15][::-1]
    fig, ax = plt.subplots(figsize=(7.5, 0.32 * len(rows) + 1.2))
    v = [abs(r["auc"] - 0.5) * 2 for r in rows]
    ax.barh(range(len(rows)), v, height=0.6, color=blue)
    ax.set_yticks(range(len(rows)), [f"{r['label']}  ({r['group']})" for r in rows])
    for i, (r, x) in enumerate(zip(rows, v)):
        ax.text(x + 0.01, i, f"AUC {r['auc']:.2f}", va="center", color=ink2, fontsize=8)
    ax.set_xlim(0, 1.15)
    ax.set_xlabel("separation |2·AUC − 1|  (0 = no information, 1 = perfect)")
    ax.grid(axis="x", color=grid, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("Single-feature separation: mismatching vs matching steps", loc="left", color=ink)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "feature_auc.png"), dpi=140)
    plt.close(fig)

    k = len(R["clusters"])
    P = y
    rr = np.hypot(p0[:, 0], p0[:, 1])
    fig, axs = plt.subplots(2, k, figsize=(3.4 * k + 0.6, 6.2), squeeze=False)
    for c in range(k):
        m = lab == c
        info = R["clusters"][c]
        for row, (h, v, xl, yl) in enumerate(((p0[:, 2], rr, "z [mm]", "r [mm]"), (p0[:, 0], p0[:, 1], "x [mm]", "y [mm]"))):
            a = axs[row, c]
            a.scatter(h[P], v[P], s=1, color=ctx, lw=0, rasterized=True)
            a.scatter(h[m], v[m], s=2 if m.sum() > 500 else 14, color=blue, lw=0, rasterized=True)
            a.set_xlabel(xl)
            if c == 0:
                a.set_ylabel(yl)
            if row == 1:
                a.set_aspect("equal", adjustable="datalim")
        sh = 100 * info["mismatch_length_share"]
        sh = f"{sh:.0f}" if sh >= 1 else f"{sh:.1f}"
        axs[0, c].set_title(f"mode {c}: ~{info['steps']:,} steps, {sh} % of length",
                            loc="left", color=ink, fontsize=9)
    fig.suptitle("Where each failure mode sits (blue) among all mismatching steps (grey)", x=0.01, ha="left", color=ink)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "clusters.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
