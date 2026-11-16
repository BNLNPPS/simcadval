// SPDX-License-Identifier: LGPL-3.0-or-later
// OCCT side of dd4hep2step: neutral GeometryDesc -> XCAF document -> STEP.
//
// Export semantics (chosen so the STEP file encodes the same material map as
// the Geant4 geometry DD4hep builds):
//   * A volume becomes a part ("body") iff its material differs from the world
//     material. World-material regions stay empty and resolve to the world
//     volume on import, i.e. to the same material.
//   * A body is the volume's own shape minus all daughter shapes (assemblies
//     looked through), so exported solids never overlap and G4OCCT's flat
//     import is a valid placement.
//   * One prototype per DD4hep volume, placed via XCAF references (STEP
//     instancing; shared G4 logical volumes on import).
//   * Each body carries an XDE material attribute (the DD4hep material name)
//     and the part name "<volume>@<id>".
#include "geo_desc.hh"

#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <GProp_GProps.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <TDF_ChildIterator.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TCollection_HAsciiString.hxx>
#include <TDataStd_Name.hxx>
#include <TDocStd_Application.hxx>
#include <TDocStd_Document.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_MaterialTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <gp_Ax2.hxx>
#include <gp_Elips.hxx>
#include <gp_GTrsf.hxx>
#include <gp_Trsf.hxx>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

namespace simcadval {
namespace {

constexpr double kDeg = M_PI / 180.0;

gp_Trsf ToGp(const Trsf& t) {
  gp_Trsf g;
  const double* m = t.m;
  g.SetValues(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10], m[11]);
  return g;
}

TopoDS_Shape Moved(const TopoDS_Shape& s, const gp_Trsf& t) {
  return BRepBuilderAPI_Transform(s, t, /*copy=*/true).Shape();
}

// ── Primitive builders ──────────────────────────────────────────────────────

TopoDS_Shape FuseAll(const std::vector<TopoDS_Shape>& parts);

// Body of revolution about Z from an (r, z) profile (Tube, TubeSeg, Cone,
// ConeSeg, Pcon), built from analytic cone/cylinder sections. A single
// BRepPrimAPI_MakeRevol of the profile is exact in memory but its
// surfaces-of-revolution do not survive the OCCT STEP write/read round trip
// (parts come back as shells with wrong volume), so avoid it.
TopoDS_Shape RevolveProfile(const ShapeDesc& d) {
  const double phi1 = d.p[0], dphi = d.p[1];
  const bool full = dphi >= 360.0 - 1e-9;
  const double ang = full ? 2 * M_PI : dphi * kDeg;
  auto frustum = [&](double r1, double r2, double z1, double z2) -> TopoDS_Shape {
    gp_Ax2 ax(gp_Pnt(0, 0, z1), gp_Dir(0, 0, 1), gp_Dir(std::cos(phi1 * kDeg), std::sin(phi1 * kDeg), 0));
    if (std::abs(r1 - r2) < 1e-12) return BRepPrimAPI_MakeCylinder(ax, r1, z2 - z1, ang).Shape();
    if (r1 > r2) {
      // Keep the cone semi-angle positive (STEP conical_surface requires it; a
      // negative one does not survive the OCCT write/read round trip): build
      // from z2 along -Z. Sweeping about -Z runs clockwise seen from +Z, so
      // start at phi1+dphi to cover the same [phi1, phi1+dphi] segment.
      const double a = (phi1 + (full ? 0.0 : dphi)) * kDeg;
      gp_Ax2 flipped(gp_Pnt(0, 0, z2), gp_Dir(0, 0, -1), gp_Dir(std::cos(a), std::sin(a), 0));
      return BRepPrimAPI_MakeCone(flipped, r2, r1, z2 - z1, ang).Shape();
    }
    return BRepPrimAPI_MakeCone(ax, r1, r2, z2 - z1, ang).Shape();
  };
  std::vector<TopoDS_Shape> sections;
  for (std::size_t i = 0; i + 1 < d.z.size(); ++i) {
    const double z1 = d.z[i], z2 = d.z[i + 1];
    if (z2 - z1 < 1e-12) continue; // radial step: no volume
    if (d.rmax[i] < 1e-12 && d.rmax[i + 1] < 1e-12) continue;
    TopoDS_Shape sec = frustum(d.rmax[i], d.rmax[i + 1], z1, z2);
    if (d.rmin[i] > 1e-12 || d.rmin[i + 1] > 1e-12)
      sec = BRepAlgoAPI_Cut(sec, frustum(d.rmin[i], d.rmin[i + 1], z1, z2)).Shape();
    sections.push_back(sec);
  }
  if (sections.empty()) throw std::runtime_error("degenerate revolution profile");
  return FuseAll(sections);
}

TopoDS_Wire PolygonWire(const std::vector<gp_Pnt>& pts) {
  BRepBuilderAPI_MakePolygon poly;
  for (const auto& p : pts) poly.Add(p);
  poly.Close();
  return poly.Wire();
}

// Ruled loft through planar polygon sections with equal vertex counts (vertex
// i joined to vertex i of the next section). Exact for Arb8/Trd/Para/Pgon/Xtru:
// TGeo side faces are bilinear between corresponding edges.
TopoDS_Shape RuledLoft(const std::vector<std::vector<gp_Pnt>>& sections) {
  BRepOffsetAPI_ThruSections loft(/*isSolid=*/true, /*ruled=*/true, 1e-6);
  loft.CheckCompatibility(false);
  for (const auto& s : sections) loft.AddWire(PolygonWire(s));
  loft.Build();
  if (!loft.IsDone()) throw std::runtime_error("ThruSections failed");
  return loft.Shape();
}

TopoDS_Shape EightPoint(const ShapeDesc& d) {
  std::vector<gp_Pnt> lo, hi;
  for (int i = 0; i < 4; ++i) lo.emplace_back(d.xy[3 * i], d.xy[3 * i + 1], d.xy[3 * i + 2]);
  for (int i = 4; i < 8; ++i) hi.emplace_back(d.xy[3 * i], d.xy[3 * i + 1], d.xy[3 * i + 2]);
  for (int i = 0; i < 4; ++i)
    if (lo[i].Distance(lo[(i + 1) % 4]) < 1e-9 || hi[i].Distance(hi[(i + 1) % 4]) < 1e-9)
      throw std::runtime_error("degenerate 8-point solid (coincident vertices) not supported");
  return RuledLoft({lo, hi});
}

// Runs of z-planes without a repeated z (a repeated z is a radial step).
std::vector<std::pair<int, int>> ZRuns(const std::vector<double>& z) {
  std::vector<std::pair<int, int>> runs;
  int start = 0;
  for (int i = 1; i < static_cast<int>(z.size()); ++i) {
    if (std::abs(z[i] - z[i - 1]) < 1e-12) {
      if (i - 1 > start) runs.emplace_back(start, i - 1);
      start = i;
    }
  }
  if (static_cast<int>(z.size()) - 1 > start) runs.emplace_back(start, static_cast<int>(z.size()) - 1);
  return runs;
}

TopoDS_Shape FuseAll(const std::vector<TopoDS_Shape>& parts) {
  if (parts.empty()) throw std::runtime_error("nothing to fuse");
  if (parts.size() == 1) return parts.front();
  TopTools_ListOfShape args, tools;
  args.Append(parts.front());
  for (std::size_t i = 1; i < parts.size(); ++i) tools.Append(parts[i]);
  BRepAlgoAPI_Fuse fuse;
  fuse.SetArguments(args);
  fuse.SetTools(tools);
  fuse.Build();
  if (!fuse.IsDone()) throw std::runtime_error("fuse failed");
  ShapeUpgrade_UnifySameDomain unify(fuse.Shape(), true, true, false);
  unify.Build();
  return unify.Shape();
}

TopoDS_Shape Pgon(const ShapeDesc& d) {
  const double phi1 = d.p[0], dphi = d.p[1];
  const int ne = static_cast<int>(d.p[2]);
  const bool full = dphi >= 360.0 - 1e-9;
  const double step = dphi / ne * kDeg;
  const double fac = 1.0 / std::cos(0.5 * step); // TGeoPgon radii are apothems
  auto ring = [&](double r, double z) {
    std::vector<gp_Pnt> v;
    for (int j = 0; j < (full ? ne : ne + 1); ++j) {
      const double a = phi1 * kDeg + j * step;
      v.emplace_back(r * fac * std::cos(a), r * fac * std::sin(a), z);
    }
    return v;
  };
  std::vector<TopoDS_Shape> pieces;
  for (auto [a, b] : ZRuns(d.z)) {
    bool anyInner = false, allInner = true;
    for (int i = a; i <= b; ++i) {
      if (d.rmin[i] > 1e-12) anyInner = true;
      else allInner = false;
    }
    if (anyInner && !allInner) throw std::runtime_error("Pgon with mixed zero/non-zero rmin in one z-run");
    if (full) {
      std::vector<std::vector<gp_Pnt>> outer, inner;
      for (int i = a; i <= b; ++i) {
        outer.push_back(ring(d.rmax[i], d.z[i]));
        inner.push_back(ring(d.rmin[i], d.z[i]));
      }
      TopoDS_Shape body = RuledLoft(outer);
      if (anyInner) body = BRepAlgoAPI_Cut(body, RuledLoft(inner)).Shape();
      pieces.push_back(body);
    } else {
      // Partial phi: the cross-section is one simple polygon (outer arc, then
      // the inner arc reversed or the centre point), so one loft is exact.
      std::vector<std::vector<gp_Pnt>> secs;
      for (int i = a; i <= b; ++i) {
        auto o = ring(d.rmax[i], d.z[i]);
        if (anyInner) {
          auto in = ring(d.rmin[i], d.z[i]);
          o.insert(o.end(), in.rbegin(), in.rend());
        } else {
          o.emplace_back(0.0, 0.0, d.z[i]);
        }
        secs.push_back(o);
      }
      pieces.push_back(RuledLoft(secs));
    }
  }
  return FuseAll(pieces);
}

TopoDS_Shape Xtru(const ShapeDesc& d) {
  std::vector<TopoDS_Shape> pieces;
  for (auto [a, b] : ZRuns(d.z)) {
    std::vector<std::vector<gp_Pnt>> secs;
    for (int i = a; i <= b; ++i) {
      std::vector<gp_Pnt> s;
      for (std::size_t j = 0; j + 1 < d.xy.size(); j += 2)
        s.emplace_back(d.x0[i] + d.scale[i] * d.xy[j], d.y0[i] + d.scale[i] * d.xy[j + 1], d.z[i]);
      secs.push_back(s);
    }
    pieces.push_back(RuledLoft(secs));
  }
  return FuseAll(pieces);
}

TopoDS_Shape Sphere(const ShapeDesc& d) {
  const double rmin = d.p[0], rmax = d.p[1], t1 = d.p[2], t2 = d.p[3], p1 = d.p[4], dphi = d.p[5] - d.p[4];
  gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(std::cos(p1 * kDeg), std::sin(p1 * kDeg), 0));
  // OCCT takes latitudes (from the XY plane); TGeo polar angles (from +Z).
  auto make = [&](double r) {
    return BRepPrimAPI_MakeSphere(ax, r, (90.0 - t2) * kDeg, (90.0 - t1) * kDeg,
                                  dphi >= 360.0 - 1e-9 ? 2 * M_PI : dphi * kDeg)
        .Shape();
  };
  TopoDS_Shape body = make(rmax);
  if (rmin > 0) body = BRepAlgoAPI_Cut(body, make(rmin)).Shape();
  return body;
}

TopoDS_Shape Torus(const ShapeDesc& d) {
  const double R = d.p[0], rmin = d.p[1], rmax = d.p[2], phi1 = d.p[3], dphi = d.p[4];
  gp_Ax2 ax(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(std::cos(phi1 * kDeg), std::sin(phi1 * kDeg), 0));
  const double ang = dphi >= 360.0 - 1e-9 ? 2 * M_PI : dphi * kDeg;
  TopoDS_Shape body = BRepPrimAPI_MakeTorus(ax, R, rmax, ang).Shape();
  if (rmin > 0) body = BRepAlgoAPI_Cut(body, BRepPrimAPI_MakeTorus(ax, R, rmin, ang).Shape()).Shape();
  return body;
}

TopoDS_Shape Eltu(const ShapeDesc& d) {
  double a = d.p[0], b = d.p[1];
  const double dz = d.p[2];
  gp_Ax2 ax(gp_Pnt(0, 0, -dz), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
  if (a < b) { // gp_Elips needs major >= minor; put the major axis along Y
    ax.SetXDirection(gp_Dir(0, 1, 0));
    std::swap(a, b);
  }
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(gp_Elips(ax, a, b));
  TopoDS_Face face = BRepBuilderAPI_MakeFace(BRepBuilderAPI_MakeWire(edge).Wire());
  return BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, 2 * dz)).Shape();
}

class ShapeBuilder {
public:
  explicit ShapeBuilder(const std::vector<ShapeDesc>& s) : fDesc(s) {}

  const TopoDS_Shape& Get(int idx) {
    auto it = fCache.find(idx);
    if (it == fCache.end()) it = fCache.emplace(idx, Build(fDesc.at(idx))).first;
    return it->second;
  }

private:
  TopoDS_Shape Build(const ShapeDesc& d) {
    switch (d.kind) {
    case ShapeKind::Box:
      return BRepPrimAPI_MakeBox(gp_Pnt(d.p[3] - d.p[0], d.p[4] - d.p[1], d.p[5] - d.p[2]),
                                 gp_Pnt(d.p[3] + d.p[0], d.p[4] + d.p[1], d.p[5] + d.p[2]))
          .Shape();
    case ShapeKind::Revolve: return RevolveProfile(d);
    case ShapeKind::Pgon: return Pgon(d);
    case ShapeKind::EightPoint: return EightPoint(d);
    case ShapeKind::Xtru: return Xtru(d);
    case ShapeKind::Sphere: return Sphere(d);
    case ShapeKind::Torus: return Torus(d);
    case ShapeKind::Eltu: return Eltu(d);
    case ShapeKind::Boolean: {
      const TopoDS_Shape l = Moved(Get(d.a), ToGp(d.ma));
      const TopoDS_Shape r = Moved(Get(d.b), ToGp(d.mb));
      if (d.p[0] == 0) return BRepAlgoAPI_Fuse(l, r).Shape();
      if (d.p[0] == 1) return BRepAlgoAPI_Cut(l, r).Shape();
      return BRepAlgoAPI_Common(l, r).Shape();
    }
    case ShapeKind::Scaled: {
      gp_GTrsf g;
      g.SetValue(1, 1, d.p[0]);
      g.SetValue(2, 2, d.p[1]);
      g.SetValue(3, 3, d.p[2]);
      return BRepBuilderAPI_GTransform(Get(d.a), g, true).Shape();
    }
    case ShapeKind::Unsupported: break;
    }
    throw std::runtime_error("unsupported TGeo shape class " + d.tgeoClass);
  }

  const std::vector<ShapeDesc>& fDesc;
  std::map<int, TopoDS_Shape> fCache;
};

// ── Exporter ────────────────────────────────────────────────────────────────

struct VolumeReport {
  std::string status = "skipped";
  bool exported = false, valid = true;
  int ncutters = 0;
  double expect_mm3 = 0, occt_mm3 = 0, build_ms = 0;
  // STEP read-back check: the part must come back as exactly one solid with the same volume.
  int rt_solids = -1;
  double rt_mm3 = 0;
  bool rt_ok = false;
};

class Exporter {
public:
  Exporter(const GeometryDesc& g, Handle(XCAFDoc_ShapeTool) st, Handle(XCAFDoc_MaterialTool) mt,
           Perturbation* perturb)
      : fG(g), fShapes(g.shapes), fST(st), fMT(mt), fPerturb(perturb), fRep(g.volumes.size()),
        fProto(g.volumes.size()), fDone(g.volumes.size(), false) {}

  TDF_Label Build(int id) {
    if (fDone[id]) return fProto[id];
    fDone[id] = true;
    const VolumeDesc& v = fG.volumes[id];
    VolumeReport& rep = fRep[id];

    TDF_Label body;
    if (id != 0 && !v.assembly && v.material != fG.worldMaterial) {
      const auto t0 = std::chrono::steady_clock::now();
      try {
        TopoDS_Shape shape = fShapes.Get(v.shape);
        rep.expect_mm3 = fG.shapes[v.shape].capacity_mm3;
        TopTools_ListOfShape cutters;
        CollectCutters(v, gp_Trsf(), cutters, rep.expect_mm3);
        rep.ncutters = cutters.Extent();
        if (!cutters.IsEmpty()) {
          TopTools_ListOfShape args;
          args.Append(shape);
          BRepAlgoAPI_Cut cut;
          cut.SetArguments(args);
          cut.SetTools(cutters);
          cut.Build();
          if (!cut.IsDone()) throw std::runtime_error("mother-minus-daughters cut failed");
          shape = cut.Shape();
        }
        GProp_GProps props;
        BRepGProp::VolumeProperties(shape, props);
        rep.occt_mm3 = props.Mass();
        rep.valid = BRepCheck_Analyzer(shape).IsValid();
        body = fST->AddShape(shape, /*makeAssembly=*/false, /*makePrepare=*/false);
        TDataStd_Name::Set(body, PartName(v.name, id).c_str());
        fMT->SetMaterial(body, MaterialLabel(v));
        rep.exported = true;
        rep.status = "ok";
      } catch (const std::exception& e) {
        rep.status = std::string("FAILED: ") + e.what();
        std::cerr << "dd4hep2step: volume " << PartName(v.name, id) << ": " << e.what() << "\n";
      }
      rep.build_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    } else {
      rep.status = v.assembly ? "assembly" : (id == 0 ? "world" : "world-material");
    }

    std::vector<std::pair<TDF_Label, TopLoc_Location>> comps;
    for (const DaughterDesc& dd : v.daughters) {
      TDF_Label child = Build(dd.vol);
      if (child.IsNull()) continue;
      gp_Trsf loc = ToGp(dd.t);
      if (fPerturb && fG.volumes[dd.vol].name.find(fPerturb->match) != std::string::npos) {
        gp_Trsf sh;
        sh.SetTranslation(gp_Vec(fPerturb->shift[0], fPerturb->shift[1], fPerturb->shift[2]));
        loc = sh * loc;
        ++fPerturb->applied;
      }
      comps.emplace_back(child, TopLoc_Location(loc));
    }
    if (comps.empty()) return fProto[id] = body;

    TDF_Label asmL = fST->NewShape();
    TDataStd_Name::Set(asmL, (PartName(v.name, id) + "/asm").c_str());
    if (!body.IsNull()) fST->AddComponent(asmL, body, TopLoc_Location());
    for (auto& [lbl, loc] : comps) fST->AddComponent(asmL, lbl, loc);
    return fProto[id] = asmL;
  }

  const std::vector<VolumeReport>& Reports() const { return fRep; }
  std::vector<VolumeReport>& Reports() { return fRep; }

private:
  // Daughter shapes in the mother frame; assemblies are transparent.
  void CollectCutters(const VolumeDesc& v, const gp_Trsf& frame, TopTools_ListOfShape& out, double& expect) {
    for (const DaughterDesc& dd : v.daughters) {
      const VolumeDesc& d = fG.volumes[dd.vol];
      const gp_Trsf t = frame * ToGp(dd.t);
      if (d.assembly) {
        CollectCutters(d, t, out, expect);
        continue;
      }
      out.Append(Moved(fShapes.Get(d.shape), t));
      expect -= fG.shapes[d.shape].capacity_mm3;
    }
  }

  TDF_Label MaterialLabel(const VolumeDesc& v) {
    if (auto it = fMats.find(v.material); it != fMats.end()) return it->second;
    TDF_Label l = fMT->AddMaterial(new TCollection_HAsciiString(v.material.c_str()),
                                   new TCollection_HAsciiString("DD4hep material"), v.density_g_cm3,
                                   new TCollection_HAsciiString("density"),
                                   new TCollection_HAsciiString("POSITIVE_RATIO_MEASURE"));
    return fMats[v.material] = l;
  }

  const GeometryDesc& fG;
  ShapeBuilder fShapes;
  Handle(XCAFDoc_ShapeTool) fST;
  Handle(XCAFDoc_MaterialTool) fMT;
  Perturbation* fPerturb;
  std::vector<VolumeReport> fRep;
  std::vector<TDF_Label> fProto;
  std::vector<bool> fDone;
  std::map<std::string, TDF_Label> fMats;
};

std::string Esc(const std::string& s) {
  std::string o;
  for (char c : s) {
    if (c == '"' || c == '\\') o += '\\';
    o += c;
  }
  return o;
}

} // namespace

int WriteStep(const GeometryDesc& g, const std::string& stepPath, const std::string& reportPath,
              Perturbation* perturb) {
  Handle(TDocStd_Application) app = new TDocStd_Application;
  Handle(TDocStd_Document) doc;
  app->NewDocument("MDTV-XCAF", doc);
  Handle(XCAFDoc_ShapeTool) st = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
  Handle(XCAFDoc_MaterialTool) mt = XCAFDoc_DocumentTool::MaterialTool(doc->Main());

  const auto t0 = std::chrono::steady_clock::now();
  Exporter ex(g, st, mt, perturb);
  if (ex.Build(0).IsNull()) {
    std::cerr << "dd4hep2step: nothing to export\n";
    return 1;
  }
  st->UpdateAssemblies();

  STEPCAFControl_Writer writer;
  writer.SetNameMode(true);
  writer.SetMaterialMode(true);
  writer.SetColorMode(false);
  if (!writer.Transfer(doc, STEPControl_AsIs) || writer.Write(stepPath.c_str()) != IFSelect_RetDone) {
    std::cerr << "dd4hep2step: STEP write failed\n";
    return 1;
  }
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  // Read the file back the way G4OCCT does and check every part survived.
  auto& reps = ex.Reports();
  {
    Handle(TDocStd_Document) rdoc;
    app->NewDocument("MDTV-XCAF", rdoc);
    STEPCAFControl_Reader reader;
    reader.SetNameMode(true);
    if (reader.ReadFile(stepPath.c_str()) != IFSelect_RetDone || !reader.Transfer(rdoc)) {
      std::cerr << "dd4hep2step: cannot read back " << stepPath << "\n";
      return 1;
    }
    Handle(XCAFDoc_ShapeTool) rst = XCAFDoc_DocumentTool::ShapeTool(rdoc->Main());
    for (TDF_ChildIterator it(rst->Label(), false); it.More(); it.Next()) {
      const TDF_Label l = it.Value();
      if (!rst->IsSimpleShape(l) || rst->IsAssembly(l)) continue;
      Handle(TDataStd_Name) nm;
      if (!l.FindAttribute(TDataStd_Name::GetID(), nm)) continue;
      const int id = IdFromPartName(TCollection_AsciiString(nm->Get()).ToCString());
      if (id < 0 || id >= static_cast<int>(reps.size())) continue;
      const TopoDS_Shape sh = rst->GetShape(l);
      TopTools_IndexedMapOfShape solids;
      TopExp::MapShapes(sh, TopAbs_SOLID, solids);
      GProp_GProps props;
      BRepGProp::VolumeProperties(sh, props);
      VolumeReport& r = reps[id];
      r.rt_solids = solids.Extent();
      r.rt_mm3 = props.Mass();
      r.rt_ok = r.rt_solids == 1 && std::abs(r.rt_mm3 - r.occt_mm3) <= 1e-6 * std::max(1.0, std::abs(r.occt_mm3));
    }
  }
  int nExp = 0, nFail = 0, nInvalid = 0, nReadback = 0;
  for (const auto& r : reps) nReadback += r.exported && !r.rt_ok;
  std::ofstream js(reportPath);
  js << "{\n  \"step\": \"" << Esc(stepPath) << "\",\n  \"world_material\": \"" << Esc(g.worldMaterial)
     << "\",\n  \"seconds\": " << secs;
  if (perturb)
    js << ",\n  \"perturb\": {\"match\": \"" << Esc(perturb->match) << "\", \"shift_mm\": [" << perturb->shift[0]
       << ", " << perturb->shift[1] << ", " << perturb->shift[2] << "], \"placements\": " << perturb->applied
       << "}";
  js << ",\n  \"volumes\": [\n";
  for (std::size_t i = 0; i < reps.size(); ++i) {
    const auto& r = reps[i];
    const auto& v = g.volumes[i];
    nExp += r.exported;
    nFail += r.status.rfind("FAILED", 0) == 0;
    nInvalid += r.exported && !r.valid;
    js << "    {\"id\": " << i << ", \"name\": \"" << Esc(v.name) << "\", \"material\": \"" << Esc(v.material)
       << "\", \"shape\": \"" << (v.shape >= 0 ? g.shapes[v.shape].tgeoClass : "assembly")
       << "\", \"status\": \"" << Esc(r.status) << "\", \"exported\": " << (r.exported ? "true" : "false")
       << ", \"ndaughters\": " << v.daughters.size() << ", \"ncutters\": " << r.ncutters
       << ", \"expect_mm3\": " << r.expect_mm3 << ", \"occt_mm3\": " << r.occt_mm3
       << ", \"brep_valid\": " << (r.valid ? "true" : "false") << ", \"build_ms\": " << r.build_ms
       << ", \"readback_solids\": " << r.rt_solids << ", \"readback_mm3\": " << r.rt_mm3
       << ", \"readback_ok\": " << (r.rt_ok ? "true" : "false") << "}"
       << (i + 1 < reps.size() ? "," : "") << "\n";
  }
  js << "  ]\n}\n";

  // OCCT 8.0 writes XDE material attributes but its reader does not re-attach
  // them to shape labels, so G4OCCT falls back to part names. Emit the
  // part-name -> material map in G4OCCT's materials.xml format.
  std::ofstream mx(stepPath + ".materials.xml");
  mx << "<materials>\n";
  for (std::size_t i = 0; i < reps.size(); ++i)
    if (reps[i].exported)
      mx << "  <material stepName=\"" << PartName(g.volumes[i].name, static_cast<int>(i)) << "\" geant4Name=\""
         << g.volumes[i].material << "\"/>\n";
  mx << "</materials>\n";

  std::printf("dd4hep2step: %zu volumes, %d exported, %d failed, %d invalid BRep, %d failed STEP read-back, "
              "%.2f s -> %s\n", reps.size(), nExp, nFail, nInvalid, nReadback, secs, stepPath.c_str());
  for (std::size_t i = 0; i < reps.size(); ++i)
    if (reps[i].exported && !reps[i].rt_ok)
      std::printf("  read-back FAILED %-28s %-14s solids=%d V_read/V_written-1 = %+.3e\n",
                  PartName(g.volumes[i].name, static_cast<int>(i)).c_str(), g.shapes[g.volumes[i].shape].tgeoClass.c_str(),
                  reps[i].rt_solids, reps[i].occt_mm3 > 0 ? reps[i].rt_mm3 / reps[i].occt_mm3 - 1 : 0.0);
  if (perturb) std::printf("dd4hep2step: perturbation applied to %d placements\n", perturb->applied);
  for (std::size_t i = 0; i < reps.size(); ++i) {
    const auto& r = reps[i];
    if (!r.exported || r.expect_mm3 <= 0) continue;
    const double rel = (r.occt_mm3 - r.expect_mm3) / r.expect_mm3;
    if (std::abs(rel) > 1e-6)
      std::printf("  volume check %-28s %-14s V_occt/V_tgeo-1 = %+.3e%s\n", PartName(g.volumes[i].name, i).c_str(),
                  g.shapes[g.volumes[i].shape].tgeoClass.c_str(), rel, r.valid ? "" : "  [invalid BRep]");
  }
  return nFail ? 1 : (nReadback ? 4 : 0);
}

} // namespace simcadval
