// WrapperWorld.hh — a small navigated world for driving
// BBSimOpBoundaryProcess::PostStepDoIt on a hand-built step, with no run
// manager, physics list or event loop.
//
// The world is a 50 mm G4_Galactic box holding one vacuum_wg G4Box named
// "crack" (the HFSS dataset id: the test writes <root>/waveguides/crack_<f>GHz_
// Ephi={0,1} and calls BBRConfigManager::SetDataDir(<root>) before the first
// Shoot). Crack-local axes: x = the crack axis (entry and exit faces), y = the
// long dimension, z = the gap. With slabMat the crack is a daughter of a slab of
// that material, half (crackHalf.x, 20 mm, 20 mm), centred at the origin, so the
// crack's x faces lie in the slab's x faces as in the test world.
//
// Navigation rule (the trap that hid the 2026-09-16 wrong-normal bug): placing
// the navigator AT the hit point is not enough. GetGlobalExitNormal then
// returns the nearest face of the entered solid, and a photon leaving the crack
// through its side wall would be reflected about the slab face, so a test built
// that way would pin the bug. Shoot therefore drives the tracking navigator
// through a real step as G4Transportation does: LocateGlobalPointAndSetup(start),
// ComputeStep(start, k, kInfinity, safety), SetGeometricallyLimitedStep(), then
// a relative LocateGlobalPointAndSetup(hit, &k, true) before the post-step
// touchable is built.
//
// Needs G4ENSDFSTATEDATA (G4OpticalPhoton builds the particle table): register
// the test program with NEEDS_ENSDF.
#ifndef WrapperWorld_hh
#define WrapperWorld_hh

#include "BBRTestSupport.hh"

#include "BBRMaterials.hh"
#include "BBSimOpBoundaryProcess.hh"

#include "G4Box.hh"
#include "G4DynamicParticle.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4NistManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleChange.hh"
#include "G4RotationMatrix.hh"
#include "G4SafetyHelper.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Transform3D.hh"
#include "G4TransportationManager.hh"
#include "Randomize.hh"
#include "geomdefs.hh"

namespace wrapperworld {

struct Result {
  BBSimOpBoundaryProcess::BBRBoundaryStatus status;   // GetLastBBRStatus()
  G4ThreeVector dir, pol, pos;                         // proposed by the particle change
  G4double freqGHz;                                    // GetLastHFSSFrequencyGHz()
  G4int uniformsUsed;                                  // draws from the scripted engine
  G4double edep;                                       // proposed local energy deposit
  G4TrackStatus trackStatus;                           // proposed track status
  G4ThreeVector hit;                                   // where the navigated step ended
};

class World {
 public:
  // crackHalf: half-lengths (x along the crack axis, y long, z gap). rot is the
  // ACTIVE rotation of the crack (G4Transform3D): Rx(+90 deg) puts local y
  // (long) along world +z and local z (gap) along world -y.
  explicit World(G4ThreeVector crackHalf = {2 * mm, 5 * mm, 0.026 * mm},
                 const G4RotationMatrix& rot = {}, G4Material* slabMat = nullptr) {
    auto* gal = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    auto* wlv = new G4LogicalVolume(new G4Box("World", 50 * mm, 50 * mm, 50 * mm), gal, "World");
    fWorld = new G4PVPlacement(nullptr, {}, wlv, "World", nullptr, false, 0);
    G4LogicalVolume* mother = wlv;
    if (slabMat) {
      auto* slv = new G4LogicalVolume(new G4Box("Slab", crackHalf.x(), 20 * mm, 20 * mm), slabMat, "Slab");
      new G4PVPlacement(nullptr, {}, slv, "Slab", wlv, false, 0);
      mother = slv;
    }
    auto* clv = new G4LogicalVolume(new G4Box("crack", crackHalf.x(), crackHalf.y(), crackHalf.z()),
                                    BBRMaterials::GetVacuumWG(), "crack");
    new G4PVPlacement(G4Transform3D(rot, G4ThreeVector()), clv, "crack", mother, false, 0);

    auto* tm = G4TransportationManager::GetTransportationManager();
    fNav = tm->GetNavigatorForTracking();
    fNav->SetWorldVolume(fWorld);
    tm->GetSafetyHelper()->InitialiseHelper();

    // The wrapped stock process is used only by pass-through steps. Both are
    // left to the end of the process (no destruction-order dependence at exit).
    fProc = new BBSimOpBoundaryProcess();
    fProc->RegisterProcess(new G4OpBoundaryProcess());
  }

  BBSimOpBoundaryProcess& Process() { return *fProc; }

  // Photon at start with direction k (unit), polarization pol and energy E_eV.
  // Drives the navigator through a real step from start to the next boundary
  // along k, builds the G4Step (post-step status postStatus, fGeomBoundary by
  // default) and calls the wrapper's PostStepDoIt with eng as the random engine.
  Result Shoot(G4ThreeVector start, G4ThreeVector k, G4ThreeVector pol, G4double E_eV,
               bbrtest::ScriptedEngine& eng, G4StepStatus postStatus = fGeomBoundary) {
    G4Random::setTheEngine(&eng);
    k = k.unit();
    fNav->LocateGlobalPointAndSetup(start, &k, false, false);
    G4TouchableHandle preTh(fNav->CreateTouchableHistory());
    G4double safety = 0.;
    const G4double s = fNav->ComputeStep(start, k, kInfinity, safety);
    const G4ThreeVector hit = start + s * k;
    fNav->SetGeometricallyLimitedStep();
    fNav->LocateGlobalPointAndSetup(hit, &k, true, false);
    G4TouchableHandle postTh(fNav->CreateTouchableHistory());

    auto* dp = new G4DynamicParticle(G4OpticalPhoton::Definition(), k, E_eV * eV);
    dp->SetPolarization(pol);
    G4Track track(dp, 0., hit);
    track.SetTouchableHandle(postTh);
    G4Step step;
    step.SetStepLength(s);
    G4StepPoint* pre = step.GetPreStepPoint();
    pre->SetPosition(start);
    pre->SetTouchableHandle(preTh);
    pre->SetMaterial(preTh->GetVolume()->GetLogicalVolume()->GetMaterial());
    track.SetStep(&step);
    track.SetStepLength(s);
    G4StepPoint* post = step.GetPostStepPoint();
    post->SetStepStatus(postStatus);
    post->SetPosition(hit);
    post->SetTouchableHandle(postTh);
    post->SetMaterial(postTh->GetVolume()->GetLogicalVolume()->GetMaterial());

    const std::size_t before = eng.Draws();
    auto* pc = static_cast<G4ParticleChange*>(fProc->PostStepDoIt(track, step));
    Result r;
    r.status = fProc->GetLastBBRStatus();
    r.freqGHz = fProc->GetLastHFSSFrequencyGHz();
    r.uniformsUsed = G4int(eng.Draws() - before);
    r.dir = *pc->GetMomentumDirection();
    r.pol = *pc->GetPolarization();
    r.pos = *pc->GetPosition();
    r.edep = pc->GetLocalEnergyDeposit();
    r.trackStatus = pc->GetTrackStatus();
    r.hit = hit;
    return r;
  }

 private:
  G4VPhysicalVolume* fWorld = nullptr;
  G4Navigator* fNav = nullptr;
  BBSimOpBoundaryProcess* fProc = nullptr;
};

// Photon energy [eV] of frequency nu_GHz.
inline G4double EnergyEV(G4double nu_GHz) { return CLHEP::h_Planck * nu_GHz * 1e9 * CLHEP::hertz / eV; }

}  // namespace wrapperworld

#endif
