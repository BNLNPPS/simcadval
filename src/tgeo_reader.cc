// SPDX-License-Identifier: LGPL-3.0-or-later
// ROOT/DD4hep side of dd4hep2step: compact file -> neutral GeometryDesc.
#include "geo_desc.hh"
#include "geo_ids.hh"

#include <DD4hep/Detector.h>

#include <TGeoArb8.h>
#include <TGeoBBox.h>
#include <TGeoBoolNode.h>
#include <TGeoCompositeShape.h>
#include <TGeoCone.h>
#include <TGeoEltu.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoPcon.h>
#include <TGeoPgon.h>
#include <TGeoScaledShape.h>
#include <TGeoSphere.h>
#include <TGeoTorus.h>
#include <TGeoTube.h>
#include <TGeoXtru.h>

#include <map>

namespace simcadval {
namespace {

constexpr double kCm = 10.0; // TGeo cm -> mm

Trsf ToTrsf(const TGeoMatrix* m) {
  Trsf t;
  if (m == nullptr || m->IsIdentity()) return t;
  const double* r = m->GetRotationMatrix();
  const double* x = m->GetTranslation();
  const double v[12] = {r[0], r[1], r[2], x[0] * kCm, r[3], r[4], r[5], x[1] * kCm,
                        r[6], r[7], r[8], x[2] * kCm};
  for (int i = 0; i < 12; ++i) t.m[i] = v[i];
  return t;
}

class ShapeTable {
public:
  explicit ShapeTable(std::vector<ShapeDesc>& out) : fOut(out) {}

  int Get(const TGeoShape* s) {
    if (auto it = fIdx.find(s); it != fIdx.end()) return it->second;
    ShapeDesc d = Describe(s);
    d.tgeoClass = s->ClassName();
    d.capacity_mm3 = s->Capacity() * 1000.0;
    fOut.push_back(std::move(d));
    return fIdx[s] = static_cast<int>(fOut.size()) - 1;
  }

private:
  static void Revolve(ShapeDesc& d, std::vector<double> z, std::vector<double> rmin, std::vector<double> rmax,
                      double phi1, double dphi) {
    d.kind = ShapeKind::Revolve;
    for (auto* v : {&z, &rmin, &rmax})
      for (double& x : *v) x *= kCm;
    d.z = std::move(z);
    d.rmin = std::move(rmin);
    d.rmax = std::move(rmax);
    d.p = {phi1, dphi};
  }

  ShapeDesc Describe(const TGeoShape* s) {
    ShapeDesc d;
    const std::string cls = s->ClassName();
    if (cls == "TGeoBBox") {
      auto* b = static_cast<const TGeoBBox*>(s);
      const double* o = b->GetOrigin();
      d.kind = ShapeKind::Box;
      d.p = {b->GetDX() * kCm, b->GetDY() * kCm, b->GetDZ() * kCm, o[0] * kCm, o[1] * kCm, o[2] * kCm};
    } else if (cls == "TGeoTube" || cls == "TGeoTubeSeg") {
      auto* t = static_cast<const TGeoTube*>(s);
      double p1 = 0, dp = 360;
      if (cls == "TGeoTubeSeg") {
        auto* ts = static_cast<const TGeoTubeSeg*>(s);
        p1 = ts->GetPhi1();
        dp = ts->GetPhi2() - ts->GetPhi1();
      }
      Revolve(d, {-t->GetDz(), t->GetDz()}, {t->GetRmin(), t->GetRmin()}, {t->GetRmax(), t->GetRmax()}, p1, dp);
    } else if (cls == "TGeoCone" || cls == "TGeoConeSeg") {
      auto* c = static_cast<const TGeoCone*>(s);
      double p1 = 0, dp = 360;
      if (cls == "TGeoConeSeg") {
        auto* cs = static_cast<const TGeoConeSeg*>(s);
        p1 = cs->GetPhi1();
        dp = cs->GetPhi2() - cs->GetPhi1();
      }
      Revolve(d, {-c->GetDz(), c->GetDz()}, {c->GetRmin1(), c->GetRmin2()}, {c->GetRmax1(), c->GetRmax2()}, p1,
              dp);
    } else if (cls == "TGeoPcon" || cls == "TGeoPgon") {
      auto* p = static_cast<const TGeoPcon*>(s);
      std::vector<double> z, rmin, rmax;
      for (int i = 0; i < p->GetNz(); ++i) {
        z.push_back(p->GetZ(i));
        rmin.push_back(p->GetRmin(i));
        rmax.push_back(p->GetRmax(i));
      }
      Revolve(d, z, rmin, rmax, p->GetPhi1(), p->GetDphi());
      if (cls == "TGeoPgon") {
        d.kind = ShapeKind::Pgon;
        d.p.push_back(static_cast<const TGeoPgon*>(s)->GetNedges());
      }
    } else if (cls == "TGeoTrd1" || cls == "TGeoTrd2" || cls == "TGeoArb8" || cls == "TGeoTrap" ||
               cls == "TGeoGtra" || cls == "TGeoPara") {
      // TGeo's own 8 mesh vertices: 4 at -dz then 4 at +dz, corresponding.
      d.kind = ShapeKind::EightPoint;
      d.xy.resize(24);
      s->SetPoints(d.xy.data());
      for (double& x : d.xy) x *= kCm;
    } else if (cls == "TGeoXtru") {
      auto* x = static_cast<const TGeoXtru*>(s);
      d.kind = ShapeKind::Xtru;
      for (int j = 0; j < x->GetNvert(); ++j) {
        d.xy.push_back(x->GetX(j) * kCm);
        d.xy.push_back(x->GetY(j) * kCm);
      }
      for (int i = 0; i < x->GetNz(); ++i) {
        d.z.push_back(x->GetZ(i) * kCm);
        d.scale.push_back(x->GetScale(i));
        d.x0.push_back(x->GetXOffset(i) * kCm);
        d.y0.push_back(x->GetYOffset(i) * kCm);
      }
    } else if (cls == "TGeoSphere") {
      auto* sp = static_cast<const TGeoSphere*>(s);
      d.kind = ShapeKind::Sphere;
      d.p = {sp->GetRmin() * kCm, sp->GetRmax() * kCm, sp->GetTheta1(), sp->GetTheta2(), sp->GetPhi1(),
             sp->GetPhi2()};
    } else if (cls == "TGeoTorus") {
      auto* t = static_cast<const TGeoTorus*>(s);
      d.kind = ShapeKind::Torus;
      d.p = {t->GetR() * kCm, t->GetRmin() * kCm, t->GetRmax() * kCm, t->GetPhi1(), t->GetDphi()};
    } else if (cls == "TGeoEltu") {
      auto* e = static_cast<const TGeoEltu*>(s);
      d.kind = ShapeKind::Eltu;
      d.p = {e->GetA() * kCm, e->GetB() * kCm, e->GetDz() * kCm};
    } else if (cls == "TGeoCompositeShape") {
      const TGeoBoolNode* n = static_cast<const TGeoCompositeShape*>(s)->GetBoolNode();
      d.kind = ShapeKind::Boolean;
      d.p = {n->GetBooleanOperator() == TGeoBoolNode::kGeoUnion         ? 0.0
             : n->GetBooleanOperator() == TGeoBoolNode::kGeoSubtraction ? 1.0
                                                                         : 2.0};
      d.a = Get(n->GetLeftShape());
      d.b = Get(n->GetRightShape());
      d.ma = ToTrsf(n->GetLeftMatrix());
      d.mb = ToTrsf(n->GetRightMatrix());
    } else if (cls == "TGeoScaledShape") {
      auto* sc = static_cast<const TGeoScaledShape*>(s);
      const double* f = sc->GetScale()->GetScale();
      d.kind = ShapeKind::Scaled;
      d.p = {f[0], f[1], f[2]};
      d.a = Get(sc->GetShape());
    }
    return d;
  }

  std::vector<ShapeDesc>& fOut;
  std::map<const TGeoShape*, int> fIdx;
};

} // namespace

GeometryDesc ReadDD4hep(const std::string& compact) {
  dd4hep::Detector& desc = dd4hep::Detector::getInstance();
  desc.fromCompact(compact);
  const TGeoVolume* top = desc.manager().GetTopVolume();

  GeometryDesc g;
  g.worldMaterial = top->GetMaterial()->GetName();
  const auto order = EnumerateVolumes(top);
  std::map<const TGeoVolume*, int> ids;
  for (std::size_t i = 0; i < order.size(); ++i) ids[order[i]] = static_cast<int>(i);

  ShapeTable shapes(g.shapes);
  for (std::size_t i = 0; i < order.size(); ++i) {
    const TGeoVolume* v = order[i];
    VolumeDesc vd;
    vd.id = static_cast<int>(i);
    vd.name = v->GetName();
    vd.material = v->GetMaterial() ? v->GetMaterial()->GetName() : "";
    vd.density_g_cm3 = v->GetMaterial() ? v->GetMaterial()->GetDensity() : 0;
    vd.assembly = v->IsAssembly();
    if (!vd.assembly) vd.shape = shapes.Get(v->GetShape());
    for (int k = 0; k < v->GetNdaughters(); ++k) {
      const TGeoNode* n = v->GetNode(k);
      vd.daughters.push_back({ids.at(n->GetVolume()), ToTrsf(n->GetMatrix()), n->GetName()});
    }
    g.volumes.push_back(std::move(vd));
  }
  return g;
}

} // namespace simcadval
