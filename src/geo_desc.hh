// SPDX-License-Identifier: LGPL-3.0-or-later
// Neutral (std-only) geometry description passed between the ROOT/DD4hep
// translation unit and the OCCT translation unit. ROOT's global Printf() and
// OCCT's extern "C" Printf() cannot be declared in the same TU, and OCCT's
// Handle() macro breaks DD4hep's Handle<> template, so the two never meet.
//
// Units: mm and degrees. Matrices: row-major 3x4 [R | t], p_mother = R p + t.
#pragma once

#include <cctype>
#include <string>
#include <vector>

namespace simcadval {

struct Trsf {
  double m[12] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
};

enum class ShapeKind { Box, Revolve, Pgon, EightPoint, Xtru, Sphere, Torus, Eltu, Boolean, Scaled, Unsupported };

struct ShapeDesc {
  ShapeKind kind = ShapeKind::Unsupported;
  std::string tgeoClass;
  // Box {dx,dy,dz,ox,oy,oz}; Revolve {phi1,dphi}; Pgon {phi1,dphi,nedges};
  // Sphere {rmin,rmax,theta1,theta2,phi1,phi2}; Torus {R,rmin,rmax,phi1,dphi};
  // Eltu {a,b,dz}; Scaled {sx,sy,sz}; Boolean {op: 0 union,1 subtraction,2 intersection}
  std::vector<double> p;
  std::vector<double> z, rmin, rmax;   // Revolve, Pgon; Xtru uses z
  std::vector<double> xy;              // EightPoint: 24 coords; Xtru: polygon (x,y)*
  std::vector<double> scale, x0, y0;   // Xtru per-section
  int a = -1, b = -1;                  // Boolean operands / Scaled operand (shape indices)
  Trsf ma, mb;                         // Boolean operand placements
  double capacity_mm3 = 0;             // TGeoShape::Capacity()
};

struct DaughterDesc {
  int vol = -1;
  Trsf t;
  std::string node;
};

struct VolumeDesc {
  int id = -1;
  std::string name, material;
  double density_g_cm3 = 0;
  bool assembly = false;
  int shape = -1;
  std::vector<DaughterDesc> daughters;
};

struct GeometryDesc {
  std::vector<ShapeDesc> shapes;
  std::vector<VolumeDesc> volumes; // index == volume id (geo_ids.hh numbering); [0] is the top
  std::string worldMaterial;
};

inline std::string PartName(const std::string& volName, int id) { return volName + "@" + std::to_string(id); }

// Inverse of PartName; -1 when there is no "@<id>" suffix. G4OCCT may append
// "_<n>" to duplicate names, so parsing stops at the first non-digit.
inline int IdFromPartName(const std::string& name) {
  const auto at = name.rfind('@');
  if (at == std::string::npos) return -1;
  int id = -1;
  for (std::size_t i = at + 1; i < name.size() && std::isdigit(static_cast<unsigned char>(name[i])); ++i)
    id = (id < 0 ? 0 : id * 10) + (name[i] - '0');
  return id;
}

// ROOT side (tgeo_reader.cc): load a compact file and describe the TGeo tree.
GeometryDesc ReadDD4hep(const std::string& compact);

// OCCT side (occt_writer.cc): build the XCAF document and write STEP + report.
struct Perturbation {
  std::string match; // substring of the daughter volume name
  double shift[3] = {0, 0, 0};
  int applied = 0;
};
int WriteStep(const GeometryDesc& g, const std::string& stepPath, const std::string& reportPath,
              Perturbation* perturb);

} // namespace simcadval
