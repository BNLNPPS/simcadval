#!/usr/bin/env python3
"""Side-by-side step listing of one ray in two scans.
    show_ray.py A_PREFIX B_PREFIX EVT [s_min s_max] [--report X.step.json]"""
import argparse
import json

import numpy as np

from compare import RAY_DT, STEP_DT

ap = argparse.ArgumentParser()
ap.add_argument("a")
ap.add_argument("b")
ap.add_argument("evt", type=int)
ap.add_argument("smin", type=float, nargs="?", default=-1)
ap.add_argument("smax", type=float, nargs="?", default=1e99)
ap.add_argument("--report")
args = ap.parse_args()

vn = {}
if args.report:
    vn = {v["id"]: v["name"] for v in json.load(open(args.report))["volumes"]}
for pre in (args.a, args.b):
    meta = json.load(open(pre + ".meta.json"))
    mats = [m["name"] for m in meta["materials"]]
    r = np.fromfile(pre + ".rays.bin", dtype=RAY_DT)[args.evt]
    st = np.fromfile(pre + ".steps.bin", dtype=STEP_DT)
    st = st[st["evt"] == args.evt]
    print(f"== {pre} [{meta['world']}] evt {args.evt} origin ({r['ox']:.3f},{r['oy']:.3f},{r['oz']:.3f}) "
          f"dir ({r['dx']:.6f},{r['dy']:.6f},{r['dz']:.6f}) steps {len(st)} len {r['total_len']:.6f}")
    o = np.array([r["ox"], r["oy"], r["oz"]])
    d = np.array([r["dx"], r["dy"], r["dz"]])
    for s in st:
        if s["s0"] + s["len"] < args.smin or s["s0"] > args.smax:
            continue
        p = o + d * s["s0"]
        print(f"  s0 {s['s0']:14.9f}  len {s['len']:14.9f}  {mats[s['mat']]:14s} vol {s['vol']:4d} "
              f"{vn.get(s['vol'], ''):22s} at ({p[0]:9.3f},{p[1]:9.3f},{p[2]:9.3f}) r={np.hypot(p[0], p[1]):8.3f}")
