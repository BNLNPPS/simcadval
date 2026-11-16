#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Diff two simcadval_scan outputs (reference A, candidate B) and compute similarity metrics.

    compare.py A_PREFIX B_PREFIX OUTDIR [--label-a dd4hep] [--label-b step] [--export-report X.step.json]

Both scans must have been run with identical --rays/--ray-mode/--points/--box/--seed so that
ray i and point j are the same query in both worlds (checked).

Ray metrics use an exact overlay of the two piecewise-constant material functions along each
ray (union of all step boundaries), so no along-ray sampling resolution enters.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

STEP_DT = np.dtype([("evt", "<i4"), ("mat", "<i4"), ("vol", "<i4"), ("pad", "<i4"), ("s0", "<f8"), ("len", "<f8")])
RAY_DT = np.dtype([("evt", "<i4"), ("nsteps", "<i4"), ("ox", "<f8"), ("oy", "<f8"), ("oz", "<f8"),
                   ("dx", "<f8"), ("dy", "<f8"), ("dz", "<f8"), ("total_len", "<f8"), ("cpu_us", "<f8")])
PT_DT = np.dtype([("x", "<f8"), ("y", "<f8"), ("z", "<f8"), ("mat", "<i4"), ("vol", "<i4"), ("cpu_us", "<f8")])
TOLS_MM = [1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]  # 1 nm .. 1 mm


class Scan:
    def __init__(self, prefix, matnames):
        self.prefix = prefix
        self.meta = json.load(open(prefix + ".meta.json"))
        self.rays = np.fromfile(prefix + ".rays.bin", dtype=RAY_DT)
        self.steps = np.fromfile(prefix + ".steps.bin", dtype=STEP_DT)
        self.points = np.fromfile(prefix + ".points.bin", dtype=PT_DT)
        self.hung = []
        if os.path.exists(prefix + ".hung"):
            self.hung = [ln.split()[:2] for ln in open(prefix + ".hung") if ln.strip()]
        # local material index -> global id (shared name table)
        lut = []
        assert [m["idx"] for m in self.meta["materials"]] == list(range(len(self.meta["materials"])))
        for m in self.meta["materials"]:
            if m["name"] not in matnames:
                matnames[m["name"]] = dict(id=len(matnames), X0=m["X0_mm"], lam=m["lambdaI_mm"], rho=m["density_g_cm3"])
            lut.append(matnames[m["name"]]["id"])
        self.lut = np.array(lut + [-1], dtype=np.int64)  # index -1 (no material) -> -1
        self.smat = self.lut[self.steps["mat"]]
        self.pmat = self.lut[self.points["mat"]]


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def pct(a, qs=(50, 90, 99, 100)):
    if len(a) == 0:
        return {f"p{q}": None for q in qs}
    v = np.percentile(a, qs)
    return {f"p{q}": float(x) for q, x in zip(qs, v)}


def overlay(A, B, nrays, L):
    """Elementary intervals of the union of A and B step boundaries.
    Returns evt, s_start, length, matA, matB, volA, volB (global material ids; -1 = outside world)."""
    recs = []
    for src, S, gm in ((0, A.steps, A.smat), (1, B.steps, B.smat)):
        end = S["s0"] + S["len"]
        # last step of each event -> an END record (material -1 beyond the world)
        last = np.r_[S["evt"][1:] != S["evt"][:-1], True]
        evt = np.r_[S["evt"], S["evt"][last]]
        s = np.r_[S["s0"], end[last]]
        mat = np.r_[gm, np.full(last.sum(), -1)]
        vol = np.r_[S["vol"], np.full(last.sum(), -9)]
        recs.append((np.full(len(evt), src), evt, s, mat, vol))
    src = np.concatenate([r[0] for r in recs])
    evt = np.concatenate([r[1] for r in recs]).astype(np.int64)
    s = np.concatenate([r[2] for r in recs])
    mat = np.concatenate([r[3] for r in recs])
    vol = np.concatenate([r[4] for r in recs])
    order = np.lexsort((src, s, evt))  # by evt, then s, A before B on ties
    src, evt, s, mat, vol = src[order], evt[order], s[order], mat[order], vol[order]
    idx = np.arange(len(s))
    lastA = np.maximum.accumulate(np.where(src == 0, idx, -1))
    lastB = np.maximum.accumulate(np.where(src == 1, idx, -1))
    # records before the first A (or B) record of an event would inherit the previous event:
    first = np.r_[True, evt[1:] != evt[:-1]]
    evstart = np.maximum.accumulate(np.where(first, idx, 0))
    okA = lastA >= evstart
    okB = lastB >= evstart
    matA = np.where(okA, mat[np.maximum(lastA, 0)], -1)
    matB = np.where(okB, mat[np.maximum(lastB, 0)], -1)
    volA = np.where(okA, vol[np.maximum(lastA, 0)], -9)
    volB = np.where(okB, vol[np.maximum(lastB, 0)], -9)
    same = np.r_[evt[1:] == evt[:-1], False]
    length = np.where(same, np.r_[s[1:] - s[:-1], 0.0], 0.0)
    keep = length > 0
    return evt[keep], s[keep], length[keep], matA[keep], matB[keep], volA[keep], volB[keep]


def boundaries(S, gm):
    """Material-change boundaries along each ray (consecutive same-material steps merged)."""
    chg = np.r_[False, (S["evt"][1:] == S["evt"][:-1]) & (gm[1:] != gm[:-1])]
    k = np.nonzero(chg)[0]
    return S["evt"][k].astype(np.int64), S["s0"][k], gm[k - 1], gm[k]


def match_boundaries(qa, qb, L):
    """For each boundary in qa, nearest boundary in qb on the same ray: distance and whether the
    material transition (from, to) is the same. Ordering key evt*L+s (L > ray length)."""
    ea, sa, fa, ta = qa
    eb, sb, fb, tb = qb
    if len(sa) == 0:
        return np.zeros(0), np.zeros(0, bool)
    if len(sb) == 0:
        return np.full(len(sa), np.inf), np.zeros(len(sa), bool)
    kb = eb * L + sb
    ob = np.argsort(kb, kind="stable")
    kb, eb, sb, fb, tb = kb[ob], eb[ob], sb[ob], fb[ob], tb[ob]
    pos = np.searchsorted(kb, ea * L + sa)
    best_d = np.full(len(sa), np.inf)
    best_same = np.zeros(len(sa), bool)
    for c in (pos - 1, pos):
        valid = (c >= 0) & (c < len(sb))
        cc = np.clip(c, 0, len(sb) - 1)
        valid &= eb[cc] == ea
        d = np.where(valid, np.abs(sb[cc] - sa), np.inf)
        same = valid & (fb[cc] == fa) & (tb[cc] == ta)
        better = d < best_d
        best_d = np.where(better, d, best_d)
        best_same = np.where(better, same, best_same)
    return best_d, best_same


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("outdir")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--export-report", default=None, help="dd4hep2step JSON report (volume names)")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--sig-mm", type=float, default=1e-6,
                    help="mismatch intervals shorter than this are floating-point noise (map/CSV only)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    mats = {}
    A, B = Scan(args.a, mats), Scan(args.b, mats)
    # Rays lost to a hung query (watchdog) exist in only one scan: compare the
    # common set and renumber events 0..n-1.
    common = np.intersect1d(A.rays["evt"], B.rays["evt"])
    missing = {"a": np.setdiff1d(B.rays["evt"], common).tolist(), "b": np.setdiff1d(A.rays["evt"], common).tolist()}
    for S in (A, B):
        S.rays = S.rays[np.isin(S.rays["evt"], common)]
        keep = np.isin(S.steps["evt"], common)
        S.steps = S.steps[keep]
        S.smat = S.smat[keep]
        S.rays["evt"] = np.searchsorted(common, S.rays["evt"])
        S.steps["evt"] = np.searchsorted(common, S.steps["evt"])
    for key in ("rays", "ray_mode", "points", "seed", "box_mm", "target_mm"):
        if A.meta.get(key) != B.meta.get(key):
            sys.exit(f"query sets differ: {key} {A.meta[key]} != {B.meta[key]}")
    if not (np.allclose(A.rays["ox"], B.rays["ox"]) and np.allclose(A.rays["dz"], B.rays["dz"])):
        sys.exit("ray origins/directions differ between scans")
    nm = len(mats)
    names = [None] * nm
    X0 = np.zeros(nm)
    LAM = np.zeros(nm)
    for n, m in mats.items():
        names[m["id"]] = n
        X0[m["id"]] = m["X0"]
        LAM[m["id"]] = m["lam"]
    world_mat = mats[A.meta["world_material"]]["id"]
    nr = len(A.rays)
    L = 4.0 * math.sqrt(sum(h * h for h in A.meta["world_half_mm"])) + 1.0

    vnames = {}
    if args.export_report and os.path.exists(args.export_report):
        for v in json.load(open(args.export_report))["volumes"]:
            vnames[v["id"]] = v["name"]

    R = {"label_a": args.label_a, "label_b": args.label_b, "a": A.meta, "b": B.meta,
         "rays_missing": missing, "hung": {"a": A.hung, "b": B.hung}}
    for k in ("a", "b"):
        R[k] = {kk: vv for kk, vv in R[k].items() if kk != "materials"}
    R["materials"] = {n: {"X0_mm": m["X0"], "lambdaI_mm": m["lam"], "density_g_cm3": m["rho"]} for n, m in mats.items()}

    # ── Rays: exact overlay ───────────────────────────────────────────────
    evt, s, ln, ma, mb, va, vb = overlay(A, B, nr, L)
    total = np.bincount(evt, weights=ln, minlength=nr)
    mism = ma != mb
    mis_len = np.bincount(evt, weights=ln * mism, minlength=nr)
    nonworld = (ma != world_mat) | (mb != world_mat)
    nonworld = nonworld & (ma >= 0) & (mb >= 0)
    nw_len = np.bincount(evt, weights=ln * nonworld, minlength=nr)
    vmis = va != vb
    vmis_len = np.bincount(evt, weights=ln * vmis, minlength=nr)

    # confusion (length, mm); index nm == "outside world"
    C = np.zeros((nm + 1, nm + 1))
    np.add.at(C, (np.where(ma < 0, nm, ma), np.where(mb < 0, nm, mb)), ln)
    lab = names + ["<outside>"]
    iou = {}
    for i in range(nm):
        inter = C[i, i]
        union = C[i, :].sum() + C[:, i].sum() - inter
        if union > 0:
            iou[names[i]] = {"iou": inter / union, "len_a_mm": C[i, :].sum(), "len_b_mm": C[:, i].sum()}

    tl_a = A.rays["total_len"]
    tl_b = B.rays["total_len"]

    def budget(S, gm, table):
        w = np.where(gm >= 0, 1.0 / table[np.maximum(gm, 0)], 0.0)
        return np.bincount(S["evt"], weights=S["len"] * w, minlength=nr)

    xa, xb = budget(A.steps, A.smat, X0), budget(B.steps, B.smat, X0)
    la, lb = budget(A.steps, A.smat, LAM), budget(B.steps, B.smat, LAM)
    pos = xa > 1e-12
    relx = np.where(pos, (xb - xa) / np.where(pos, xa, 1), 0.0)
    rell = np.where(la > 1e-12, (lb - la) / np.where(la > 1e-12, la, 1), 0.0)

    ba, bb = boundaries(A.steps, A.smat), boundaries(B.steps, B.smat)
    dab, sab = match_boundaries(ba, bb, L)
    dba, sba = match_boundaries(bb, ba, L)
    bmetrics = {}
    for t in TOLS_MM:
        rec = float(np.mean((dab <= t) & sab)) if len(dab) else float("nan")
        prec = float(np.mean((dba <= t) & sba)) if len(dba) else float("nan")
        f1 = 2 * rec * prec / (rec + prec) if rec + prec > 0 else 0.0
        bmetrics[f"{t:g}mm"] = {"recall": rec, "precision": prec, "f1": f1}

    # collapsed material sequences per ray
    def seqs(S, gm):
        keep = np.r_[True, (S["evt"][1:] != S["evt"][:-1]) | (gm[1:] != gm[:-1])]
        e, g = S["evt"][keep], gm[keep]
        cuts = np.searchsorted(e, np.arange(nr + 1))
        return [g[cuts[i]:cuts[i + 1]] for i in range(nr)]

    qa, qb = seqs(A.steps, A.smat), seqs(B.steps, B.smat)
    seq_equal = np.array([len(x) == len(y) and np.array_equal(x, y) for x, y in zip(qa, qb)])

    tot = total.sum()
    nwt = nw_len.sum()
    R["rays"] = {
        "n": int(nr),
        "steps_a": int(len(A.steps)), "steps_b": int(len(B.steps)),
        "total_length_mm": float(tot),
        "nonworld_length_mm": float(nwt),
        "material_mismatch_length_mm": float(mis_len.sum()),
        "material_mismatch_fraction_total": float(mis_len.sum() / tot),
        "material_mismatch_fraction_nonworld": float(mis_len.sum() / nwt) if nwt > 0 else None,
        "material_agreement_nonworld": float(1 - mis_len.sum() / nwt) if nwt > 0 else None,
        "volume_mismatch_fraction_nonworld": float(vmis_len.sum() / nwt) if nwt > 0 else None,
        "rays_with_material_mismatch": {f"> {t:g} mm": int((mis_len > t).sum()) for t in TOLS_MM},
        "mismatch_intervals_longer_than": {f"{t:g} mm": int((mism & (ln > t)).sum()) for t in TOLS_MM},
        "per_ray_mismatch_mm": pct(mis_len),
        "ray_total_length_absdiff_mm": pct(np.abs(tl_a - tl_b)),
        "rays_length_differ": int((np.abs(tl_a - tl_b) > 1e-6).sum()),  # e.g. events aborted on stuck tracks
        "x_over_x0": {"mean_a": float(xa.mean()), "mean_b": float(xb.mean()),
                      "rel_diff": pct(np.abs(relx[pos])),
                      "rays_rel_diff_gt": {f"{t:g}": int((np.abs(relx[pos]) > t).sum()) for t in (1e-9, 1e-6, 1e-4, 1e-2)},
                      "abs_diff": pct(np.abs(xb - xa))},
        "lambda_over_lambdaI": {"mean_a": float(la.mean()), "mean_b": float(lb.mean()),
                                "rel_diff": pct(np.abs(rell[la > 1e-12]))},
        "boundaries": {"n_a": int(len(dab)), "n_b": int(len(dba)),
                       "nearest_dist_a_to_b_mm": pct(dab[np.isfinite(dab)]),
                       "nearest_dist_b_to_a_mm": pct(dba[np.isfinite(dba)]),
                       "unmatched_a": int((~np.isfinite(dab)).sum()), "unmatched_b": int((~np.isfinite(dba)).sum()),
                       "by_tolerance": bmetrics},
        "sequence_identical_fraction": float(seq_equal.mean()),
        "cpu_us_per_ray": {"a": pct(A.rays["cpu_us"], (50, 99)), "b": pct(B.rays["cpu_us"], (50, 99)),
                           "mean_a": float(A.rays["cpu_us"].mean()), "mean_b": float(B.rays["cpu_us"].mean())},
        "material_iou": iou,
    }

    # per-volume diagnostics (length-weighted recall of A volumes in B)
    vol_rows = []
    if np.any(va >= 0):
        vmax = int(max(va.max(), vb.max())) + 1
        la_v = np.bincount(np.maximum(va, 0), weights=ln * (va >= 0), minlength=vmax)
        lb_v = np.bincount(np.maximum(vb, 0), weights=ln * (vb >= 0), minlength=vmax)
        both = np.bincount(np.maximum(va, 0), weights=ln * (va >= 0) * (va == vb), minlength=vmax)
        for v in range(vmax):
            if la_v[v] + lb_v[v] <= 0:
                continue
            miss = la_v[v] + lb_v[v] - 2 * both[v]
            vol_rows.append({"id": v, "name": vnames.get(v, ""), "len_a_mm": float(la_v[v]), "len_b_mm": float(lb_v[v]),
                             "iou": float(both[v] / (la_v[v] + lb_v[v] - both[v])), "mismatch_mm": float(miss)})
        vol_rows.sort(key=lambda r: -r["mismatch_mm"])
    R["worst_volumes"] = vol_rows[:15]

    # where are the mismatches? midpoint of each mismatching elementary interval
    rmap = np.searchsorted(A.rays["evt"], evt)
    mid = s + 0.5 * ln
    sel = mism & (ln > args.sig_mm)
    org = np.stack([A.rays["ox"], A.rays["oy"], A.rays["oz"]], 1)[rmap[sel]]
    dirs = np.stack([A.rays["dx"], A.rays["dy"], A.rays["dz"]], 1)[rmap[sel]]
    mpos = org + dirs * mid[sel, None]
    mism_tab = np.rec.fromarrays([evt[sel], mpos[:, 0], mpos[:, 1], mpos[:, 2], ln[sel], ma[sel], mb[sel]],
                                 names="evt,x,y,z,len,matA,matB")
    order = np.argsort(-mism_tab["len"])
    with open(os.path.join(args.outdir, "ray_mismatches.csv"), "w") as f:
        f.write("evt,x_mm,y_mm,z_mm,len_mm,mat_a,mat_b\n")
        for r in mism_tab[order][:5000]:
            f.write(f"{r.evt},{r.x:.6f},{r.y:.6f},{r.z:.6f},{r.len:.9g},{lab[r.matA]},{lab[r.matB]}\n")

    # ── Points ────────────────────────────────────────────────────────────
    npts = len(A.points)
    pa, pb = A.pmat, B.pmat
    agree = pa == pb
    k = int(agree.sum())
    p, lo, hi = wilson(k, npts)
    vagree = A.points["vol"] == B.points["vol"]
    pv, vlo, vhi = wilson(int(vagree.sum()), npts)
    pconf = {}
    for i, j in zip(*np.nonzero(np.ones((nm, nm)))):
        c = int(((pa == i) & (pb == j)).sum())
        if c and i != j:
            pconf[f"{names[i]} -> {names[j]}"] = c
    vf = {names[i]: {"a": float((pa == i).mean()), "b": float((pb == i).mean())} for i in range(nm)}
    R["points"] = {"n": int(npts), "material_agreement": p, "material_agreement_ci95": [lo, hi],
                   "n_disagree": int(npts - k), "volume_agreement": pv, "volume_agreement_ci95": [vlo, vhi],
                   "disagreements": pconf, "volume_fraction_by_material": vf,
                   "cpu_us_per_point": {"mean_a": float(A.points["cpu_us"].mean()), "mean_b": float(B.points["cpu_us"].mean())}}
    bad = ~agree
    with open(os.path.join(args.outdir, "point_mismatches.csv"), "w") as f:
        f.write("x_mm,y_mm,z_mm,mat_a,mat_b,vol_a,vol_b\n")
        for x, y, z, i, j, u, w in zip(A.points["x"][bad], A.points["y"][bad], A.points["z"][bad], pa[bad], pb[bad],
                                       A.points["vol"][bad], B.points["vol"][bad]):
            f.write(f"{x:.6f},{y:.6f},{z:.6f},{lab[i]},{lab[j]},{u},{w}\n")

    # headline
    R["headline"] = {
        "material_agreement_nonworld_length": R["rays"]["material_agreement_nonworld"],
        "point_material_agreement": p,
        "boundary_f1_at_1um": bmetrics["0.001mm"]["f1"],
        "boundary_f1_at_1nm": bmetrics["1e-06mm"]["f1"],
        "max_rel_dX0": R["rays"]["x_over_x0"]["rel_diff"]["p100"],
        "max_boundary_shift_mm": max(R["rays"]["boundaries"]["nearest_dist_a_to_b_mm"]["p100"] or 0,
                                     R["rays"]["boundaries"]["nearest_dist_b_to_a_mm"]["p100"] or 0),
        "sequence_identical_fraction": R["rays"]["sequence_identical_fraction"],
        "cpu_ratio_ray_b_over_a": R["rays"]["cpu_us_per_ray"]["mean_b"] / R["rays"]["cpu_us_per_ray"]["mean_a"],
    }
    json.dump(R, open(os.path.join(args.outdir, "metrics.json"), "w"), indent=1, default=float)
    np.savez_compressed(os.path.join(args.outdir, "per_ray.npz"), evt=A.rays["evt"], dx=A.rays["dx"], dy=A.rays["dy"],
                        dz=A.rays["dz"], ox=A.rays["ox"], oy=A.rays["oy"], oz=A.rays["oz"], xa=xa, xb=xb, la=la, lb=lb,
                        mis_len=mis_len, nw_len=nw_len, cpu_a=A.rays["cpu_us"], cpu_b=B.rays["cpu_us"],
                        seq_equal=seq_equal)

    print_summary(R)
    if not args.no_plots:
        plots(args, A, B, R, xa, xb, relx, pos, dab, dba, mism_tab, C, lab, bad, pa, pb)


def print_summary(R):
    h, r, p = R["headline"], R["rays"], R["points"]
    a, b = R["label_a"], R["label_b"]
    print(f"== {a} vs {b}: {r['n']} rays ({R['a']['ray_mode']}), {p['n']} points")
    if R["hung"]["a"] or R["hung"]["b"]:
        print(f"  HUNG queries (excluded): {a}: {R['hung']['a']}  {b}: {R['hung']['b']}")
    print(f"  material agreement (non-world ray length): {h['material_agreement_nonworld_length']:.9f}"
          f"   mismatch {r['material_mismatch_length_mm']:.4g} mm of {r['nonworld_length_mm']:.4g} mm")
    print(f"  point material agreement: {p['material_agreement']:.6f}  CI95 {p['material_agreement_ci95'][0]:.6f}-"
          f"{p['material_agreement_ci95'][1]:.6f}  ({p['n_disagree']} disagree)")
    print(f"  boundary F1 @1nm {h['boundary_f1_at_1nm']:.6f}  @1um {h['boundary_f1_at_1um']:.6f}"
          f"   max shift {h['max_boundary_shift_mm']:.3g} mm   (n_a={r['boundaries']['n_a']}, n_b={r['boundaries']['n_b']})")
    print(f"  X/X0 per ray: mean {r['x_over_x0']['mean_a']:.6g} vs {r['x_over_x0']['mean_b']:.6g};"
          f"  |rel diff| p50 {r['x_over_x0']['rel_diff']['p50']:.3g} p99 {r['x_over_x0']['rel_diff']['p99']:.3g}"
          f" max {r['x_over_x0']['rel_diff']['p100']:.3g}")
    mi = r["mismatch_intervals_longer_than"]
    print(f"  mismatching intervals longer than 1 nm / 1 um / 1 mm: {mi['1e-06 mm']} / {mi['0.001 mm']} / {mi['1 mm']}"
          f"   rays with >1 um mismatch: {r['rays_with_material_mismatch']['> 0.001 mm']}")
    print(f"  identical material sequence: {h['sequence_identical_fraction']:.5f} of rays;"
          f"  volume-id mismatch (non-world length) {r['volume_mismatch_fraction_nonworld']:.3g}")
    print(f"  cpu/ray mean {r['cpu_us_per_ray']['mean_a']:.1f} us vs {r['cpu_us_per_ray']['mean_b']:.1f} us"
          f" (x{h['cpu_ratio_ray_b_over_a']:.1f});  cpu/point {p['cpu_us_per_point']['mean_a']:.2f} vs"
          f" {p['cpu_us_per_point']['mean_b']:.2f} us")
    if R["worst_volumes"] and R["worst_volumes"][0]["mismatch_mm"] > 0:
        print("  worst volumes (mismatch length):")
        for v in R["worst_volumes"][:5]:
            print(f"    {v['name']}@{v['id']:<5d} IoU {v['iou']:.6f}  mismatch {v['mismatch_mm']:.4g} mm")


def plots(args, A, B, R, xa, xb, relx, pos, dab, dba, mism_tab, C, lab, bad, pa, pb):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    la, lb = args.label_a, args.label_b
    out = lambda n: os.path.join(args.outdir, n)  # noqa: E731
    ca, cb = "#3b6fb6", "#d1495b"

    # 1. material budget per ray
    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    if A.meta["ray_mode"] == "ip":
        th = np.arccos(np.clip(A.rays["dz"], -1, 1))
        xv = -np.log(np.tan(np.clip(th, 1e-9, np.pi - 1e-9) / 2))
        xl = "pseudorapidity $\\eta$ of ray"
    else:
        xv = np.arange(len(xa))
        xl = "ray index"
    o = np.argsort(xv)
    ax[0].plot(xv[o], xa[o], ".", ms=2, color=ca, label=la)
    ax[0].plot(xv[o], xb[o], ".", ms=1, color=cb, label=lb, alpha=0.7)
    ax[0].set_ylabel("X / X$_0$ along ray")
    ax[0].legend(markerscale=6)
    ax[1].plot(xv[pos], np.abs(relx[pos]) + 1e-17, ".", ms=2, color="#444")
    ax[1].set_yscale("log")
    ax[1].set_ylabel("|ΔX/X$_0$| / X/X$_0$")
    ax[1].set_xlabel(xl)
    fig.tight_layout()
    fig.savefig(out("material_budget.png"), dpi=130)
    plt.close(fig)

    # 2. boundary displacement
    fig, ax = plt.subplots(figsize=(7, 4))
    for d, c, l in ((dab, ca, f"{la} → nearest {lb}"), (dba, cb, f"{lb} → nearest {la}")):
        d = d[np.isfinite(d)]
        d = np.maximum(d, 1e-12)
        ax.hist(np.log10(d), bins=120, range=(-12, 3), histtype="step", color=c, label=l)
    ax.set_xlabel("log10 |Δs| of material-change boundary along ray [mm]  (floor 1e-12)")
    ax.set_ylabel("boundaries")
    ax.set_yscale("log")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out("boundary_shift.png"), dpi=130)
    plt.close(fig)

    # 3. where do they disagree
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
    pts_x, pts_y, pts_z = A.points["x"][bad], A.points["y"][bad], A.points["z"][bad]
    for a_, (h1, h2, n1, n2) in zip(ax, ((mism_tab["z"], np.hypot(mism_tab["x"], mism_tab["y"]), "z [mm]", "r [mm]"),
                                          (mism_tab["x"], mism_tab["y"], "x [mm]", "y [mm]"))):
        if len(mism_tab):
            sc = a_.scatter(h1, h2, s=4, c=np.log10(np.maximum(mism_tab["len"], 1e-9)), cmap="viridis")
        a_.set_xlabel(n1)
        a_.set_ylabel(n2)
    if len(pts_x):
        ax[0].plot(pts_z, np.hypot(pts_x, pts_y), "x", color=cb, ms=4, label="point disagreement")
        ax[1].plot(pts_x, pts_y, "x", color=cb, ms=4)
        ax[0].legend(loc="upper right")
    if len(mism_tab):
        fig.colorbar(sc, ax=ax, label="log10 mismatch interval length [mm]")
    ax[0].set_title(f"{len(mism_tab)} mismatching ray intervals > {args.sig_mm:g} mm, {len(pts_x)} points")
    fig.savefig(out("mismatch_map.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 4. confusion matrix (length)
    used = np.nonzero(C.sum(0) + C.sum(1) > 0)[0]
    Cu = C[np.ix_(used, used)]
    fig, ax = plt.subplots(figsize=(1.2 + 0.55 * len(used), 1 + 0.5 * len(used)))
    im = ax.imshow(np.log10(Cu + 1e-12), cmap="Blues", vmin=-6)
    ax.set_xticks(range(len(used)), [lab[i] for i in used], rotation=60, ha="right")
    ax.set_yticks(range(len(used)), [lab[i] for i in used])
    ax.set_xlabel(f"material in {lb}")
    ax.set_ylabel(f"material in {la}")
    for i in range(len(used)):
        for j in range(len(used)):
            if i != j and Cu[i, j] > 0:
                ax.text(j, i, f"{Cu[i, j]:.2g}", ha="center", va="center", fontsize=7, color="#b00")
    fig.colorbar(im, label="log10 ray length [mm]")
    fig.tight_layout()
    fig.savefig(out("confusion.png"), dpi=130)
    plt.close(fig)

    # 5. timing
    fig, ax = plt.subplots(figsize=(7, 4))
    for r, c, l in ((A.rays, ca, la), (B.rays, cb, lb)):
        ax.hist(np.log10(np.maximum(r["cpu_us"], 1e-3)), bins=100, histtype="step", color=c, label=l)
    ax.set_xlabel("log10 CPU time per geantino [µs]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out("timing.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
