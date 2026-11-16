#include "step_world.hh"

#include "geo_desc.hh"

#include "G4OCCT/G4OCCTAssemblyVolume.hh"
#include "G4OCCT/G4OCCTMaterialMap.hh"
#include "G4OCCT/G4OCCTMaterialMapReader.hh"
#include "G4OCCT/G4OCCTSolid.hh"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4PVPlacement.hh>

namespace simcadval {

G4VPhysicalVolume* BuildStepWorld(const std::string& stepPath, const std::string& matMapXml, const double half[3],
                                  G4Material* worldMat,
                                  std::vector<std::pair<const G4LogicalVolume*, int>>& ids) {
  G4OCCTMaterialMap matMap;
  for (G4Material* m : *G4Material::GetMaterialTable()) matMap.Add(m->GetName(), m);
  if (!matMapXml.empty()) {
    G4OCCTMaterialMapReader reader;
    matMap.Merge(reader.ReadFile(matMapXml));
  }
  auto* assembly = G4OCCTAssemblyVolume::FromSTEP(stepPath, matMap);

  auto* worldSolid = new G4Box("simcadval_world", half[0], half[1], half[2]);
  auto* worldLV = new G4LogicalVolume(worldSolid, worldMat, "simcadval_world");
  G4ThreeVector pos;
  G4RotationMatrix rot;
  assembly->MakeImprint(worldLV, pos, &rot);

  ids.emplace_back(worldLV, -1);
  for (const auto& [name, lv] : assembly->GetLogicalVolumes()) ids.emplace_back(lv, IdFromPartName(name));
  return new G4PVPlacement(nullptr, {}, worldLV, "simcadval_world", nullptr, false, 0);
}

int FreshOCCTInside(const G4VSolid* solid, double x, double y, double z) {
  auto* occt = dynamic_cast<const G4OCCTSolid*>(solid);
  if (occt == nullptr) return -1;
  G4OCCTSolid fresh("fresh_" + occt->GetName(), occt->GetOCCTShape());
  return static_cast<int>(fresh.Inside(G4ThreeVector(x, y, z)));
}

} // namespace simcadval
