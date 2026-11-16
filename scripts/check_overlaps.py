#!/usr/bin/env python3
"""TGeo overlap/extrusion check of a DD4hep compact. Usage: check_overlaps.py compact.xml [tol_cm]"""
import sys
import ROOT
import dd4hep  # noqa: F401
d = dd4hep.Detector.getInstance()
d.fromXML(sys.argv[1])
tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4  # cm
mgr = ROOT.gGeoManager
mgr.CheckOverlaps(tol)
ov = mgr.GetListOfOverlaps()
print(f"overlaps/extrusions > {tol*1e4:g} um: {ov.GetEntries()}")
for o in ov:
    print(f"  {'EXTRUSION' if o.IsExtrusion() else 'OVERLAP  '} {o.GetOverlap()*10:.4g} mm  "
          f"{o.GetFirstVolume().GetName()} / {o.GetSecondVolume().GetName()}  :: {o.GetTitle()}")
