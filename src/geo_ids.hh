// Shared volume-numbering convention between dd4hep2step and simcadval_scan
// (ROOT-side translation units only).
//
// Every distinct TGeoVolume reachable from the top volume gets an integer id in
// depth-first pre-order (daughters in node index order, first visit wins). The
// exporter names STEP parts "<volume name>@<id>", and the scanner maps Geant4
// logical volumes back to the same id in both worlds, so "which DD4hep volume
// am I in?" is comparable across representations.
#pragma once

#include <TGeoNode.h>
#include <TGeoVolume.h>

#include <unordered_set>
#include <vector>

namespace simcadval {

inline std::vector<const TGeoVolume*> EnumerateVolumes(const TGeoVolume* top) {
  std::vector<const TGeoVolume*> order;
  std::unordered_set<const TGeoVolume*> seen;
  std::vector<const TGeoVolume*> stack{top};
  // Explicit pre-order DFS; push daughters in reverse so index order is kept.
  while (!stack.empty()) {
    const TGeoVolume* v = stack.back();
    stack.pop_back();
    if (!seen.insert(v).second) continue;
    order.push_back(v);
    for (int i = v->GetNdaughters() - 1; i >= 0; --i) stack.push_back(v->GetNode(i)->GetVolume());
  }
  return order;
}

} // namespace simcadval
