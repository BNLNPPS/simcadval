// STEP world construction (G4OCCT/OCCT translation unit), kept apart from the
// DD4hep/ROOT translation unit; see geo_desc.hh for why.
#pragma once

#include <string>
#include <utility>
#include <vector>

class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;
class G4VSolid;

namespace simcadval {

// Imports the STEP assembly with G4OCCTAssemblyVolume into a world box of the
// given half-lengths and material. The material map holds every G4Material
// under its own name plus, if matMapXml is non-empty, the part-name entries of
// that G4OCCT materials.xml file. Returns the world PV; fills (logical volume, volume id)
// pairs for each imported part (id parsed from "<name>@<id>"), plus the world
// LV with id -1.
G4VPhysicalVolume* BuildStepWorld(const std::string& stepPath, const std::string& matMapXml, const double half[3],
                                  G4Material* worldMat,
                                  std::vector<std::pair<const G4LogicalVolume*, int>>& ids);

// Debug aid: Inside() of a freshly constructed G4OCCTSolid with the same OCCT
// shape as `solid` (no per-solid query history / caches). Returns -1 if
// `solid` is not a G4OCCTSolid, else the EInside value.
int FreshOCCTInside(const G4VSolid* solid, double x, double y, double z);

} // namespace simcadval
