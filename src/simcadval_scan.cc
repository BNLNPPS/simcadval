// SPDX-License-Identifier: LGPL-3.0-or-later
// simcadval_scan: geantino ray scan + navigator point sampling of ONE world,
// built either from a DD4hep compact file (DDG4 Geant4Converter, exactly what
// ddsim does) or from a STEP assembly (G4OCCTAssemblyVolume). Run it once per
// world with identical arguments; compare.py diffs the two outputs.
//
// The two Geant4 questions being sampled:
//   "where is the next boundary along this ray?"  -> geantino tracking steps
//   "what material am I in?"                     -> G4Navigator::LocateGlobalPointAndSetup
//
// Rays and points are a pure function of (seed, index), independent of the
// world, so both runs ask exactly the same queries.
//
// Usage:
//   simcadval_scan --compact C.xml --world dd4hep|step [--step F.step]
//                [--matmap F.step.materials.xml (default) | none]
//                --out PREFIX [--rays N] [--ray-mode ip|box|aim] [--points M]
//                [--box hx,hy,hz (mm; ray origins / points, default 0.98*world)]
//                [--target cx,cy,cz,hx,hy,hz (mm; aim mode: rays point at a uniform
//                 point of this box, from origins uniform in --box)]
//                [--seed S] [--max-steps N] [--max-event-sec T] [--hang-sec T]
//                [--first-ray K (resume; appends)] [--debug-ray EVT,SMIN,SMAX]
// A watchdog exits with status 3 and appends "ray|point <id>" to PREFIX.hung if one
// query runs longer than --hang-sec; run_case.sh resumes after the hung ray.
// Outputs: PREFIX.rays.bin, PREFIX.steps.bin, PREFIX.points.bin, PREFIX.meta.json

#include "geo_ids.hh"
#include "step_world.hh"

#include <DD4hep/Detector.h>
#include <DDG4/Geant4Converter.h>
#include <DDG4/Geant4GeometryInfo.h>

#include <G4Box.hh>
#include <G4Event.hh>
#include <G4Geantino.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Navigator.hh>
#include <G4ParticleGun.hh>
#include <G4RunManager.hh>
#include <G4StateManager.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4UImanager.hh>
#include <G4UserEventAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4VExceptionHandler.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ios.hh>

#include <TGeoManager.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <atomic>
#include <thread>
#include <unordered_map>
#include <vector>

namespace {

struct Options {
  std::string compact, world = "dd4hep", step, matmap, out = "scan";
  long rays = 10000, points = 100000;
  long maxSteps = 200000;     // per event; a geantino looping on tiny steps is killed
  double maxEventSec = 20.0;  // per event wall-clock guard (checked between steps)
  double hangSec = 60.0;      // watchdog: a single query running this long = hang
  long firstRay = 0;          // resume after a hang: rays [firstRay, rays) appended
  std::string rayMode = "ip";
  double box[3] = {-1, -1, -1};
  double target[6] = {0, 0, 0, 0, 0, 0}; // aim mode: centre + half-lengths (mm)
  std::uint64_t seed = 12345;
};

// ── Deterministic query generation ──────────────────────────────────────────

std::uint64_t Mix(std::uint64_t x) { // splitmix64 finaliser
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

struct Ray {
  G4ThreeVector origin, dir;
};

G4ThreeVector Isotropic(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> u(0.0, 1.0);
  const double c = 2 * u(rng) - 1, phi = 2 * M_PI * u(rng), s = std::sqrt(1 - c * c);
  return {s * std::cos(phi), s * std::sin(phi), c};
}

G4ThreeVector InBox(std::mt19937_64& rng, const double* h) {
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  return {u(rng) * h[0], u(rng) * h[1], u(rng) * h[2]};
}

Ray MakeRay(const Options& o, long i) {
  std::mt19937_64 rng(Mix(o.seed ^ Mix(static_cast<std::uint64_t>(i) * 2 + 1)));
  Ray r;
  r.origin = o.rayMode == "ip" ? G4ThreeVector() : InBox(rng, o.box);
  if (o.rayMode == "aim") {
    const G4ThreeVector tgt = G4ThreeVector(o.target[0], o.target[1], o.target[2]) + InBox(rng, o.target + 3);
    r.dir = (tgt - r.origin).unit();
  } else {
    r.dir = Isotropic(rng);
  }
  return r;
}

G4ThreeVector MakePoint(const Options& o, long i) {
  std::mt19937_64 rng(Mix(o.seed ^ Mix(static_cast<std::uint64_t>(i) * 2 + 2)));
  return InBox(rng, o.box);
}

// ── Output records (little-endian, read with numpy structured dtypes) ───────

#pragma pack(push, 1)
struct StepRec {
  std::int32_t evt, mat, vol, pad;
  double s0, len;
};
struct RayRec {
  std::int32_t evt, nsteps;
  double ox, oy, oz, dx, dy, dz, total_len, cpu_us;
};
struct PointRec {
  double x, y, z;
  std::int32_t mat, vol;
  double cpu_us;
};
#pragma pack(pop)

// ── Shared state between user actions ───────────────────────────────────────

struct State {
  Options opt;
  std::unordered_map<const G4LogicalVolume*, int> volId; // -1 = world-equivalent
  std::unordered_map<const G4Material*, int> matIdx;
  std::vector<const G4Material*> mats;
  std::vector<StepRec> evSteps;
  std::ofstream raysOut, stepsOut;
  Ray current;
  std::chrono::steady_clock::time_point evStart;
  // Watchdog state: what the main thread is doing and since when (ns since epoch).
  std::atomic<long> wdItem{-1};      // ray id, or -(point index)-2 while sampling points
  std::atomic<long long> wdSince{0};
  long totalSteps = 0, stuckSteps = 0, evStepCount = 0;
  std::vector<int> killedEvents;
  G4VPhysicalVolume* world = nullptr;
  double worldHalf[3] = {0, 0, 0};
  std::string worldMaterial;
  int nVolumesMapped = 0;

  // G4Material::GetIndex(): deterministic across processes (same construction
  // order), so resumed scan segments share one material numbering.
  int MatIndex(const G4Material* m) {
    if (matIdx.emplace(m, static_cast<int>(m->GetIndex())).second) mats.push_back(m);
    return static_cast<int>(m->GetIndex());
  }
  int VolIndex(const G4LogicalVolume* lv) {
    auto it = volId.find(lv);
    return it == volId.end() ? -2 : it->second; // -2 = unmapped (should not happen)
  }
};
State* gState = nullptr;

// Count Geant4 warnings by code instead of aborting; fatal still aborts.
class CountingExceptionHandler : public G4VExceptionHandler {
public:
  G4bool Notify(const char* origin, const char* code, G4ExceptionSeverity sev, const char* desc) override {
    const int n = ++counts[code];
    if (n <= 3)
      G4cerr << "G4Exception [" << code << "] from " << origin << ": " << desc << G4endl;
    return sev == FatalException || sev == FatalErrorInArgument;
  }
  std::map<std::string, int> counts;
};

// ── Geometry ────────────────────────────────────────────────────────────────

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  explicit DetectorConstruction(State& s) : fS(s) {}

  G4VPhysicalVolume* Construct() override {
    dd4hep::Detector& desc = dd4hep::Detector::getInstance();
    // Both modes convert the DD4hep geometry: the dd4hep world IS the result,
    // the step world borrows its G4Materials and world box (materials are not
    // what is under test, and identical G4Material objects keep X0 identical).
    dd4hep::sim::Geant4Converter conv(desc, dd4hep::WARNING);
    dd4hep::sim::Geant4GeometryInfo* info = conv.create(desc.world()).detach();
    G4VPhysicalVolume* ddWorld = info->world();
    auto* box = dynamic_cast<G4Box*>(ddWorld->GetLogicalVolume()->GetSolid());
    if (box == nullptr) G4Exception("simcadval", "SCV001", FatalException, "DD4hep world is not a G4Box");
    fS.worldHalf[0] = box->GetXHalfLength();
    fS.worldHalf[1] = box->GetYHalfLength();
    fS.worldHalf[2] = box->GetZHalfLength();
    G4Material* worldMat = ddWorld->GetLogicalVolume()->GetMaterial();
    fS.worldMaterial = worldMat->GetName();

    const TGeoVolume* top = desc.manager().GetTopVolume();
    const auto order = simcadval::EnumerateVolumes(top);
    const std::string tgeoWorldMat = top->GetMaterial()->GetName();

    if (fS.opt.world == "dd4hep") {
      for (std::size_t id = 0; id < order.size(); ++id) {
        const TGeoVolume* v = order[id];
        auto it = info->g4Volumes.find(dd4hep::Volume(const_cast<TGeoVolume*>(v)));
        if (it == info->g4Volumes.end() || it->second == nullptr) continue; // assemblies
        const bool worldLike = v == top || v->GetMaterial()->GetName() == tgeoWorldMat;
        fS.volId[it->second] = worldLike ? -1 : static_cast<int>(id);
        ++fS.nVolumesMapped;
      }
      fS.world = ddWorld;
      return ddWorld;
    }

    std::vector<std::pair<const G4LogicalVolume*, int>> ids;
    fS.world = simcadval::BuildStepWorld(fS.opt.step, fS.opt.matmap, fS.worldHalf, worldMat, ids);
    for (const auto& [lv, id] : ids) {
      fS.volId[lv] = id;
      fS.nVolumesMapped += id >= 0;
    }
    return fS.world;
  }

private:
  State& fS;
};

class GeantinoPhysics : public G4VUserPhysicsList {
public:
  void ConstructParticle() override { G4Geantino::GeantinoDefinition(); }
  void ConstructProcess() override { AddTransportation(); }
  void SetCuts() override {}
};

// ── Actions ─────────────────────────────────────────────────────────────────

class PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGenerator() : fGun(1) {
    fGun.SetParticleDefinition(G4Geantino::Definition());
    fGun.SetParticleEnergy(1 * GeV);
  }
  void GeneratePrimaries(G4Event* ev) override {
    gState->current = MakeRay(gState->opt, ev->GetEventID() + gState->opt.firstRay);
    fGun.SetParticlePosition(gState->current.origin);
    fGun.SetParticleMomentumDirection(gState->current.dir);
    fGun.GeneratePrimaryVertex(ev);
  }

private:
  G4ParticleGun fGun;
};

class EventAction : public G4UserEventAction {
public:
  void BeginOfEventAction(const G4Event* ev) override {
    gState->evSteps.clear();
    gState->evStepCount = 0;
    gState->evStart = std::chrono::steady_clock::now();
    gState->wdItem = ev->GetEventID() + gState->opt.firstRay;
    gState->wdSince = gState->evStart.time_since_epoch().count();
  }
  void EndOfEventAction(const G4Event* ev) override {
    State& s = *gState;
    const double us =
        std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - s.evStart).count();
    RayRec r{};
    r.evt = ev->GetEventID() + s.opt.firstRay;
    r.nsteps = static_cast<std::int32_t>(s.evSteps.size());
    r.ox = s.current.origin.x();
    r.oy = s.current.origin.y();
    r.oz = s.current.origin.z();
    r.dx = s.current.dir.x();
    r.dy = s.current.dir.y();
    r.dz = s.current.dir.z();
    r.total_len = s.evSteps.empty() ? 0 : s.evSteps.back().s0 + s.evSteps.back().len;
    r.cpu_us = us;
    s.raysOut.write(reinterpret_cast<const char*>(&r), sizeof r);
    s.stepsOut.write(reinterpret_cast<const char*>(s.evSteps.data()),
                     static_cast<std::streamsize>(s.evSteps.size() * sizeof(StepRec)));
    s.totalSteps += r.nsteps;
  }
};

class SteppingAction : public G4UserSteppingAction {
public:
  void UserSteppingAction(const G4Step* step) override {
    State& s = *gState;
    const G4StepPoint* pre = step->GetPreStepPoint();
    const double len = step->GetStepLength();
    if (++s.evStepCount % 1000 == 0 &&
        (s.evStepCount > s.opt.maxSteps ||
         std::chrono::duration<double>(std::chrono::steady_clock::now() - s.evStart).count() > s.opt.maxEventSec)) {
      const int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() + s.opt.firstRay;
      step->GetTrack()->SetTrackStatus(fStopAndKill);
      s.killedEvents.push_back(evt);
      if (s.killedEvents.size() <= 10)
        G4cerr << "simcadval: killing event " << evt << " after " << s.evStepCount << " steps at "
               << step->GetPostStepPoint()->GetPosition() << " in "
               << pre->GetPhysicalVolume()->GetName() << G4endl;
    }
    if (len <= 0) {
      ++s.stuckSteps;
      return;
    }
    StepRec rec{};
    rec.evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() + s.opt.firstRay;
    rec.mat = s.MatIndex(pre->GetMaterial());
    rec.vol = s.VolIndex(pre->GetPhysicalVolume()->GetLogicalVolume());
    rec.s0 = step->GetTrack()->GetTrackLength() - len;
    rec.len = len;
    s.evSteps.push_back(rec);
  }
};

void ParseBox(const std::string& v, double* b) {
  if (std::sscanf(v.c_str(), "%lf,%lf,%lf", &b[0], &b[1], &b[2]) != 3) {
    std::cerr << "bad --box " << v << "\n";
    std::exit(2);
  }
}

void ParseTarget(const std::string& v, double* t) {
  if (std::sscanf(v.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf", &t[0], &t[1], &t[2], &t[3], &t[4], &t[5]) != 6) {
    std::cerr << "bad --target " << v << "\n";
    std::exit(2);
  }
}

std::string Esc(const std::string& s) {
  std::string o;
  for (char c : s) {
    if (c == '"' || c == '\\') o += '\\';
    o += c;
  }
  return o;
}

// Replay ray `evt` with a bare G4Navigator (same call pattern as G4Transportation
// for a neutral particle) and, at every step start in [smin, smax], ask every
// world daughter whose solid is near the point for Inside()/DistanceToIn().
// Shows which parts are at each step point and what their solids answer.
void DebugRay(const State& s, long evt, double smin, double smax) {
  const Ray ray = MakeRay(s.opt, evt);
  G4Navigator nav;
  nav.SetWorldVolume(s.world);
  G4ThreeVector p = ray.origin;
  const G4ThreeVector& v = ray.dir;
  G4VPhysicalVolume* pv = nav.LocateGlobalPointAndSetup(p, &v, false, false);
  double sTot = 0, safety = 0;
  const G4LogicalVolume* worldLV = s.world->GetLogicalVolume();
  for (int n = 0; pv != nullptr && n < 100000; ++n) {
    const double step = nav.ComputeStep(p, v, kInfinity, safety);
    if (sTot + step >= smin && sTot <= smax) {
      auto it = s.volId.find(pv->GetLogicalVolume());
      std::printf("s=%.9f step=%.9f in %-28s id=%d mat=%s p=(%.6f,%.6f,%.6f)\n", sTot, step, pv->GetName().c_str(),
                  it == s.volId.end() ? -2 : it->second, pv->GetLogicalVolume()->GetMaterial()->GetName().c_str(),
                  p.x(), p.y(), p.z());
      std::fflush(stdout);
      for (std::size_t i = 0; i < worldLV->GetNoDaughters(); ++i) {
        const G4VPhysicalVolume* d = worldLV->GetDaughter(i);
        G4AffineTransform toLocal(d->GetRotation(), d->GetTranslation());
        toLocal.Invert();
        const G4ThreeVector lp = toLocal.TransformPoint(p), lv = toLocal.TransformAxis(v);
        const G4VSolid* so = d->GetLogicalVolume()->GetSolid();
        const double safe = so->DistanceToIn(lp);
        if (safe > 1.0) continue;
        const EInside in = so->Inside(lp);
        auto jt = s.volId.find(d->GetLogicalVolume());
        G4ThreeVector bmin, bmax;
        so->BoundingLimits(bmin, bmax);
        std::printf("    local p=(%.15g, %.15g, %.15g)  bbox=(%.15g..%.15g, %.15g..%.15g, %.15g..%.15g)\n", lp.x(),
                    lp.y(), lp.z(), bmin.x(), bmax.x(), bmin.y(), bmax.y(), bmin.z(), bmax.z());
        auto inName = [](int e) { return e == kInside ? "INSIDE" : e == kSurface ? "SURFACE" : e == kOutside ? "outside" : "n/a"; };
        std::printf("    %-32s id=%-4d %-8s safety=%.3e  DistIn(p,v)=%.9g  DistOut(p,v)=%.9g\n",
                    d->GetName().c_str(), jt == s.volId.end() ? -2 : jt->second, inName(in), safe,
                    so->DistanceToIn(lp, lv), in == kOutside ? -1.0 : so->DistanceToOut(lp, lv));
        std::fflush(stdout);
      }
    }
    if (step >= kInfinity || sTot > smax) break;
    sTot += step;
    p += step * v;
    nav.SetGeometricallyLimitedStep();
    pv = nav.LocateGlobalPointAndSetup(p, &v, true, false);
  }
}

} // namespace

int main(int argc, char** argv) {
  State state;
  gState = &state;
  Options& o = state.opt;
  long debugEvt = -1;
  double debugRange[2] = {0, 1e99};
  for (int i = 1; i + 1 < argc; i += 2) {
    const std::string k = argv[i], v = argv[i + 1];
    if (k == "--compact") o.compact = v;
    else if (k == "--world") o.world = v;
    else if (k == "--step") o.step = v;
    else if (k == "--matmap") o.matmap = v;
    else if (k == "--out") o.out = v;
    else if (k == "--rays") o.rays = std::stol(v);
    else if (k == "--points") o.points = std::stol(v);
    else if (k == "--ray-mode") o.rayMode = v;
    else if (k == "--box") ParseBox(v, o.box);
    else if (k == "--target") ParseTarget(v, o.target);
    else if (k == "--seed") o.seed = std::stoull(v);
    else if (k == "--max-steps") o.maxSteps = std::stol(v);
    else if (k == "--hang-sec") o.hangSec = std::stod(v);
    else if (k == "--first-ray") o.firstRay = std::stol(v);
    else if (k == "--max-event-sec") o.maxEventSec = std::stod(v);
    else if (k == "--debug-ray") std::sscanf(v.c_str(), "%ld,%lf,%lf", &debugEvt, &debugRange[0], &debugRange[1]);
    else {
      std::cerr << "unknown option " << k << "\n";
      return 2;
    }
  }
  if (o.compact.empty() || (o.world == "step" && o.step.empty()) ||
      (o.world != "dd4hep" && o.world != "step") || (o.rayMode != "ip" && o.rayMode != "box" && o.rayMode != "aim")) {
    std::cerr << "usage: simcadval_scan --compact C.xml --world dd4hep|step [--step F.step] --out PREFIX\n"
                 "       [--rays N] [--ray-mode ip|box|aim] [--target cx,cy,cz,hx,hy,hz] [--points M]\n"
                 "       [--box hx,hy,hz] [--seed S]\n";
    return 2;
  }

  if (o.world == "step" && o.matmap.empty()) o.matmap = o.step + ".materials.xml";
  if (o.matmap == "none") o.matmap.clear();
  dd4hep::Detector::getInstance().fromCompact(o.compact);

  auto* handler = new CountingExceptionHandler;
  G4StateManager::GetStateManager()->SetExceptionHandler(handler);

  auto* rm = new G4RunManager;
  rm->SetUserInitialization(new DetectorConstruction(state));
  rm->SetUserInitialization(new GeantinoPhysics);
  rm->SetUserAction(new PrimaryGenerator);
  rm->SetUserAction(new EventAction);
  rm->SetUserAction(new SteppingAction);
  G4UImanager::GetUIpointer()->ApplyCommand("/run/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");

  const auto tInit = std::chrono::steady_clock::now();
  rm->Initialize();
  const double initSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - tInit).count();
  for (int k = 0; k < 3; ++k)
    if (o.box[k] <= 0) o.box[k] = 0.98 * state.worldHalf[k];
  if (debugEvt >= 0) {
    DebugRay(state, debugEvt, debugRange[0], debugRange[1]);
    return 0;
  }

  const auto mode = o.firstRay > 0 ? std::ios::binary | std::ios::app : std::ios::binary;
  state.raysOut.open(o.out + ".rays.bin", mode);
  state.stepsOut.open(o.out + ".steps.bin", mode);
  // Watchdog: some G4OCCT/OCCT queries never return. Record which ray/point
  // hung, flush completed events and exit(3) so the driver can resume after it.
  std::thread([&state, &o] {
    for (;;) {
      std::this_thread::sleep_for(std::chrono::seconds(1));
      const long long since = state.wdSince.load();
      const long item = state.wdItem.load();
      if (since == 0) continue;
      const double age = 1e-9 * static_cast<double>(std::chrono::steady_clock::now().time_since_epoch().count() - since);
      if (age < o.hangSec) continue;
      std::ofstream(o.out + ".hung", std::ios::app) << (item >= 0 ? "ray " : "point ") << (item >= 0 ? item : -item - 2)
                                                    << " " << age << "\n";
      state.raysOut.flush();
      state.stepsOut.flush();
      std::fprintf(stderr, "simcadval watchdog: %s %ld hung for %.0f s, exiting\n", item >= 0 ? "ray" : "point",
                   item >= 0 ? item : -item - 2, age);
      std::_Exit(3);
    }
  }).detach();
  const auto tRays = std::chrono::steady_clock::now();
  if (o.rays > o.firstRay) rm->BeamOn(static_cast<G4int>(o.rays - o.firstRay));
  state.wdSince = 0;
  const double raySec = std::chrono::duration<double>(std::chrono::steady_clock::now() - tRays).count();
  state.raysOut.close();
  state.stepsOut.close();

  // Point sampling with a private navigator on the closed geometry.
  std::ofstream ptsOut(o.out + ".points.bin", std::ios::binary);
  G4Navigator nav;
  nav.SetWorldVolume(state.world);
  const auto tPts = std::chrono::steady_clock::now();
  for (long i = 0; i < o.points; ++i) {
    const G4ThreeVector p = MakePoint(o, i);
    state.wdItem = -i - 2;
    state.wdSince = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto t0 = std::chrono::steady_clock::now();
    G4VPhysicalVolume* pv = nav.LocateGlobalPointAndSetup(p, nullptr, false, true);
    const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
    PointRec r{p.x(), p.y(), p.z(), -1, -3, us};
    if (pv != nullptr) {
      r.mat = state.MatIndex(pv->GetLogicalVolume()->GetMaterial());
      r.vol = state.VolIndex(pv->GetLogicalVolume());
    }
    ptsOut.write(reinterpret_cast<const char*>(&r), sizeof r);
  }
  const double ptSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - tPts).count();

  std::ofstream js(o.out + ".meta.json");
  js << "{\n  \"world\": \"" << o.world << "\", \"compact\": \"" << Esc(o.compact) << "\", \"step\": \""
     << Esc(o.step) << "\",\n  \"rays\": " << o.rays << ", \"ray_mode\": \"" << o.rayMode
     << "\", \"target_mm\": [" << o.target[0] << ", " << o.target[1] << ", " << o.target[2] << ", "
     << o.target[3] << ", " << o.target[4] << ", " << o.target[5] << "], \"points\": " << o.points << ", \"seed\": " << o.seed << ",\n  \"box_mm\": [" << o.box[0]
     << ", " << o.box[1] << ", " << o.box[2] << "], \"world_half_mm\": [" << state.worldHalf[0] << ", "
     << state.worldHalf[1] << ", " << state.worldHalf[2] << "], \"world_material\": \""
     << Esc(state.worldMaterial) << "\",\n  \"init_s\": " << initSec << ", \"rays_s\": " << raySec
     << ", \"points_s\": " << ptSec << ", \"total_steps\": " << state.totalSteps
     << ", \"zero_length_steps\": " << state.stuckSteps << ", \"killed_events\": [";
  for (std::size_t i = 0; i < state.killedEvents.size(); ++i) js << (i ? ", " : "") << state.killedEvents[i];
  js << "], \"volumes_mapped\": " << state.nVolumesMapped
     << ",\n  \"g4_exceptions\": {";
  bool first = true;
  for (const auto& [code, n] : handler->counts) {
    js << (first ? "" : ", ") << "\"" << Esc(code) << "\": " << n;
    first = false;
  }
  js << "},\n  \"materials\": [\n";
  const G4MaterialTable& table = *G4Material::GetMaterialTable(); // position == GetIndex()
  for (std::size_t i = 0; i < table.size(); ++i) {
    const G4Material* m = table[i];
    js << "    {\"idx\": " << i << ", \"name\": \"" << Esc(m->GetName())
       << "\", \"density_g_cm3\": " << m->GetDensity() / (g / cm3) << ", \"X0_mm\": " << m->GetRadlen() / mm
       << ", \"lambdaI_mm\": " << m->GetNuclearInterLength() / mm << ", \"used\": "
       << (state.matIdx.count(m) ? "true" : "false") << "}" << (i + 1 < table.size() ? "," : "") << "\n";
  }
  js << "  ]\n}\n";

  std::printf("simcadval_scan[%s]: init %.1f s | %ld rays, %ld steps in %.2f s (%.1f us/ray) | "
              "%ld points in %.2f s (%.2f us/pt) | %zu materials, %d volumes mapped\n",
              o.world.c_str(), initSec, o.rays, state.totalSteps, raySec, 1e6 * raySec / std::max(1L, o.rays),
              o.points, ptSec, 1e6 * ptSec / std::max(1L, o.points), state.mats.size(), state.nVolumesMapped);
  for (const auto& [code, n] : handler->counts) std::printf("  G4Exception %s x%d\n", code.c_str(), n);
  if (!state.killedEvents.empty()) std::printf("  killed (looping) events: %zu\n", state.killedEvents.size());
  delete rm;
  return 0;
}
