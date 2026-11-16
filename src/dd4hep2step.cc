// dd4hep2step: export a DD4hep compact geometry to a STEP (AP214) assembly
// with XDE material attributes. Export semantics are documented in
// occt_writer.cc; the TGeo walk is in tgeo_reader.cc.
//
// Usage: dd4hep2step <compact.xml> <out.step> [--report out.json]
//                    [--perturb <volume-name-substring>:<dx>,<dy>,<dz>]   (mm)
//                    [--export-world-material 1]  (also export world-material volumes;
//                     for shape-fidelity tests, not for ray comparison)
//
// --perturb shifts every placement of matching volumes. It exists only to
// produce deliberately wrong CAD for validating the comparison metrics.
#include "geo_desc.hh"

#include <cstdio>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: dd4hep2step <compact.xml> <out.step> [--report out.json] [--perturb name:dx,dy,dz]\n";
    return 2;
  }
  const std::string compact = argv[1], out = argv[2];
  std::string report = out + ".json";
  simcadval::Perturbation perturb;
  bool doPerturb = false, exportWorldMat = false;
  for (int i = 3; i + 1 < argc; i += 2) {
    const std::string k = argv[i], v = argv[i + 1];
    if (k == "--report") {
      report = v;
    } else if (k == "--perturb") {
      const auto colon = v.rfind(':');
      if (colon == std::string::npos ||
          std::sscanf(v.c_str() + colon + 1, "%lf,%lf,%lf", &perturb.shift[0], &perturb.shift[1],
                      &perturb.shift[2]) != 3) {
        std::cerr << "bad --perturb spec " << v << "\n";
        return 2;
      }
      perturb.match = v.substr(0, colon);
      doPerturb = true;
    } else if (k == "--export-world-material") {
      exportWorldMat = v != "0";
    } else {
      std::cerr << "unknown option " << k << "\n";
      return 2;
    }
  }
  simcadval::GeometryDesc g = simcadval::ReadDD4hep(compact);
  if (exportWorldMat) g.worldMaterial = "<none>";
  return simcadval::WriteStep(g, out, report, doPerturb ? &perturb : nullptr);
}
