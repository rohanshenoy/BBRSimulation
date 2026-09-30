// testPassthrough — bit-exact pass-through regression for BBSimOpBoundaryProcess.
//
//   testPassthrough <stock|wrapped|mutant> <nEvents> <out.txt>
//
// A geometry with stock boundary optics only (no vacuum_wg, no material
// REFLECTIVITY): two touching dielectric boxes A and B with RINDEX, a ground
// dielectric_dielectric border surface A->B (with surface REFLECTIVITY and
// TRANSMITTANCE) and a ground dielectric_metal border surface B->World, a
// polished metal skin with EFFICIENCY on box C, and box D without RINDEX.
// Sequential run manager, fixed seed, optical photons only. One line per
// optical-photon step with every double as a hex float (%a), so two runs can be
// compared byte for byte (tests/CMakeLists.txt: cmake -E compare_files).
//
// Modes: stock = G4OpticalPhysics only; wrapped = + BBSimPhysics; mutant = a
// test-only wrapper that draws one extra random number per PostStepDoIt (the
// comparison must see it). PT_WLS=1 in the environment makes OpWLS active in
// box B, which is what exposes the wrapper's position in the process list.
//
// The program fails (exit 1) if any of the boundary statuses 1-10, 13, 14 never
// appears, or if a wrapper intercepted a step; exit 2 is a usage error.
#include "BBRTestSupport.hh"

#include "BBSimOpBoundaryProcess.hh"
#include "BBSimPhysics.hh"

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace {

// n <= 0: no RINDEX (the stock process then kills with NoRINDEX).
G4Material* MakeMat(const char* name, G4double n, G4double absLen, G4double rayLen = -1.)
{
  auto* m = new G4Material(name, 1., 1.008 * g / mole, 1e-3 * g / cm3);
  const std::vector<G4double> e = {1e-5 * eV, 1e-1 * eV};
  auto* mpt = new G4MaterialPropertiesTable();
  if (n > 0.) mpt->AddProperty("RINDEX", e, {n, n});
  if (absLen > 0.) mpt->AddProperty("ABSLENGTH", e, {absLen, absLen});
  if (rayLen > 0.) mpt->AddProperty("RAYLEIGH", e, {rayLen, rayLen});
  m->SetMaterialPropertiesTable(mpt);
  return m;
}

class Geometry : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override
  {
    auto* vac = MakeMat("PTVacuum", 1.0, -1.);
    auto* glA = MakeMat("GlassA", 1.5, 20. * cm);
    auto* glB = MakeMat("GlassB", 2.0, 50. * cm, 30. * cm);
    auto* glC = MakeMat("GlassC", 1.3, -1.);
    if (std::getenv("PT_WLS")) {
      // OpWLS active in B: the one process whose place relative to the boundary
      // process decides which random number each of them draws.
      auto* mptB = glB->GetMaterialPropertiesTable();
      mptB->AddProperty("WLSABSLENGTH", std::vector<G4double>{1e-5 * eV, 1e-1 * eV}, {40. * cm, 40. * cm});
      mptB->AddProperty("WLSCOMPONENT", std::vector<G4double>{1e-5 * eV, 1e-3 * eV}, {1., 1.});
      mptB->AddConstProperty("WLSTIMECONSTANT", 0.5 * ns);
    }
    auto* noRI = MakeMat("NoRindex", -1., -1.);

    auto* worldLV = new G4LogicalVolume(new G4Box("World", 50 * cm, 50 * cm, 50 * cm), vac, "World");
    auto* world = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0);
    auto* aLV = new G4LogicalVolume(new G4Box("A", 5 * cm, 5 * cm, 5 * cm), glA, "A");
    auto* bLV = new G4LogicalVolume(new G4Box("B", 5 * cm, 5 * cm, 5 * cm), glB, "B");
    auto* a = new G4PVPlacement(nullptr, {-5 * cm, 0, 0}, aLV, "A", worldLV, false, 0);
    auto* b = new G4PVPlacement(nullptr, {+5 * cm, 0, 0}, bLV, "B", worldLV, false, 0);

    const std::vector<G4double> e = {1e-5 * eV, 1e-1 * eV};
    // A -> B: ground dielectric_dielectric with a surface REFLECTIVITY and
    // TRANSMITTANCE, so the stock Transmission and surface-Absorption branches run.
    auto* sAB = new G4OpticalSurface("AB", unified, ground, dielectric_dielectric, 0.1);
    auto* mAB = new G4MaterialPropertiesTable();
    mAB->AddProperty("SPECULARLOBECONSTANT", e, {0.3, 0.3});
    mAB->AddProperty("SPECULARSPIKECONSTANT", e, {0.2, 0.2});
    mAB->AddProperty("BACKSCATTERCONSTANT", e, {0.5, 0.5});
    mAB->AddProperty("REFLECTIVITY", e, {0.9, 0.9});
    mAB->AddProperty("TRANSMITTANCE", e, {0.05, 0.05});
    sAB->SetMaterialPropertiesTable(mAB);
    new G4LogicalBorderSurface("AB", a, b, sAB);

    // B -> World: ground dielectric_metal; REFLECTIVITY on the SURFACE, not a material.
    auto* sBW = new G4OpticalSurface("BW", unified, ground, dielectric_metal, 0.2);
    auto* mBW = new G4MaterialPropertiesTable();
    mBW->AddProperty("REFLECTIVITY", e, {0.7, 0.7});
    mBW->AddProperty("SPECULARLOBECONSTANT", e, {0.5, 0.5});
    mBW->AddProperty("SPECULARSPIKECONSTANT", e, {0.3, 0.3});
    sBW->SetMaterialPropertiesTable(mBW);
    new G4LogicalBorderSurface("BW", b, world, sBW);

    // C: polished metal skin with a detection efficiency (Detection status).
    auto* cLV = new G4LogicalVolume(new G4Box("C", 5 * cm, 5 * cm, 5 * cm), glC, "C");
    new G4PVPlacement(nullptr, {-5 * cm, +11 * cm, 0}, cLV, "C", worldLV, false, 0);
    auto* sC = new G4OpticalSurface("Cskin", unified, polished, dielectric_metal);
    auto* mC = new G4MaterialPropertiesTable();
    mC->AddProperty("REFLECTIVITY", e, {0.5, 0.5});
    mC->AddProperty("EFFICIENCY", e, {0.6, 0.6});
    sC->SetMaterialPropertiesTable(mC);
    new G4LogicalSkinSurface("Cskin", cLV, sC);

    // D: material without RINDEX, no surface (NoRINDEX status).
    auto* dLV = new G4LogicalVolume(new G4Box("D", 5 * cm, 5 * cm, 5 * cm), noRI, "D");
    new G4PVPlacement(nullptr, {-5 * cm, -11 * cm, 0}, dLV, "D", worldLV, false, 0);
    return world;
  }
};

// Isotropic 500 GHz photons with random linear polarization from the centre of box A.
class Gun : public G4VUserPrimaryGeneratorAction {
 public:
  Gun() : fGun(1)
  {
    fGun.SetParticleDefinition(G4OpticalPhoton::Definition());
    fGun.SetParticleEnergy(2.07e-3 * eV);
    fGun.SetParticlePosition({-5 * cm, 0, 0});
  }
  void GeneratePrimaries(G4Event* ev) override
  {
    const G4double ct = 2. * G4UniformRand() - 1., st = std::sqrt(1. - ct * ct);
    const G4double ph = CLHEP::twopi * G4UniformRand();
    const G4ThreeVector k(st * std::cos(ph), st * std::sin(ph), ct);
    const G4double ang = CLHEP::twopi * G4UniformRand();
    const G4ThreeVector e1 = k.orthogonal().unit(), e2 = k.cross(e1);
    fGun.SetParticleMomentumDirection(k);
    fGun.SetParticlePolarization(std::cos(ang) * e1 + std::sin(ang) * e2);
    fGun.GeneratePrimaryVertex(ev);
  }
 private:
  G4ParticleGun fGun;
};

// The wrapper is named "BBSimOpBoundary" + the wrapped process's name
// (G4WrapperProcess::RegisterProcess appends it): print both as "OpBoundary".
const char* Canon(const G4VProcess* p)
{
  if (!p) return "none";
  if (dynamic_cast<const G4OpBoundaryProcess*>(p) || dynamic_cast<const BBSimOpBoundaryProcess*>(p))
    return "OpBoundary";
  return p->GetProcessName().c_str();
}

class Stepping : public G4UserSteppingAction {
 public:
  explicit Stepping(FILE* out) : fOut(out) {}
  std::map<int, long> statuses;   // boundary status on fGeomBoundary steps
  long intercepted = 0, steps = 0;
  G4OpBoundaryProcess* fBoundary = nullptr;
  BBSimOpBoundaryProcess* fWrapper = nullptr;

  void UserSteppingAction(const G4Step* s) override
  {
    if (!fBoundary) Resolve();
    const auto* post = s->GetPostStepPoint();
    const auto* trk = s->GetTrack();
    const int st = fBoundary->GetStatus();
    if (post->GetStepStatus() == fGeomBoundary) ++statuses[st];
    if (fWrapper && fWrapper->GetLastBBRStatus() != BBSimOpBoundaryProcess::kBBRNone) ++intercepted;
    ++steps;
    const G4ThreeVector& x = post->GetPosition();
    const G4ThreeVector& k = post->GetMomentumDirection();
    const G4ThreeVector& p = post->GetPolarization();
    const auto* pv = post->GetPhysicalVolume();
    std::fprintf(fOut, "%d %d %d %s %d %s %d %d %a %a %a %a %a %a %a %a %a %a\n",
                 G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(),
                 trk->GetTrackID(), trk->GetCurrentStepNumber(),
                 pv ? pv->GetName().c_str() : "OutOfWorld", int(post->GetStepStatus()),
                 Canon(post->GetProcessDefinedStep()), st, int(trk->GetTrackStatus()),
                 x.x(), x.y(), x.z(), k.x(), k.y(), k.z(), p.x(), p.y(), p.z(), trk->GetGlobalTime());
  }

 private:
  void Resolve()
  {
    auto* pl = G4OpticalPhoton::Definition()->GetProcessManager()->GetProcessList();
    for (std::size_t i = 0; i < pl->size(); ++i) {
      if (auto* w = dynamic_cast<BBSimOpBoundaryProcess*>((*pl)[i])) {
        fWrapper = w;
        fBoundary = w->GetWrappedProcess();
      } else if (auto* b = dynamic_cast<G4OpBoundaryProcess*>((*pl)[i])) {
        fBoundary = b;
      }
    }
  }
  FILE* fOut;
};

// Sensitivity check for the harness: one extra random number per PostStepDoIt.
class MutantWrapper : public BBSimOpBoundaryProcess {
 public:
  G4VParticleChange* PostStepDoIt(const G4Track& t, const G4Step& s) override
  {
    (void)G4UniformRand();
    return BBSimOpBoundaryProcess::PostStepDoIt(t, s);
  }
};
class MutantPhysics : public G4VPhysicsConstructor {
 public:
  MutantPhysics() : G4VPhysicsConstructor("MutantPhysics") {}
  void ConstructParticle() override {}
  void ConstructProcess() override
  {
    auto* pm = G4OpticalPhoton::Definition()->GetProcessManager();
    auto* pl = pm->GetProcessList();
    G4VProcess* b = nullptr;
    for (std::size_t i = 0; i < pl->size() && !b; ++i)
      if (dynamic_cast<G4OpBoundaryProcess*>((*pl)[i])) b = (*pl)[i];
    pm->RemoveProcess(b);
    auto* w = new MutantWrapper;
    w->RegisterProcess(b);
    pm->AddDiscreteProcess(w);
  }
};

}  // namespace

int main(int argc, char** argv)
{
  const std::string mode = argc > 1 ? argv[1] : "";
  if (argc != 4 || (mode != "stock" && mode != "wrapped" && mode != "mutant")) {
    std::fprintf(stderr, "usage: %s <stock|wrapped|mutant> <nEvents> <out.txt>  (PT_WLS=1: active WLS in box B)\n",
                 argc > 0 ? argv[0] : "testPassthrough");
    return 2;
  }
  bbrtest::Handler();   // a fatal G4Exception is reported below, not an abort
  const int nEv = std::atoi(argv[2]);
  FILE* out = std::fopen(argv[3], "w");
  if (!out) { std::fprintf(stderr, "cannot write %s\n", argv[3]); return 1; }

  Stepping* stepping = nullptr;
  try {
    auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
    rm->SetUserInitialization(new Geometry);
    // Optical photons only: Cerenkov and scintillation need charged particles
    // (G4EmSaturation creates e- after PreInit and G4OpticalPhysics then aborts).
    G4OpticalParameters::Instance()->SetProcessActivation("Cerenkov", false);
    G4OpticalParameters::Instance()->SetProcessActivation("Scintillation", false);
    auto* phys = new G4VModularPhysicsList;
    phys->RegisterPhysics(new G4OpticalPhysics);
    if (mode == "wrapped") phys->RegisterPhysics(new BBSimPhysics);
    if (mode == "mutant") phys->RegisterPhysics(new MutantPhysics);
    rm->SetUserInitialization(phys);
    rm->SetUserAction(new Gun);
    stepping = new Stepping(out);
    rm->SetUserAction(stepping);
    rm->Initialize();
    G4Random::setTheSeed(20260928);
    rm->BeamOn(nEv);
  } catch (const bbrtest::G4ExceptionCaught& e) {
    std::fprintf(stderr, "FAIL G4Exception %s from %s: %s\n", e.code.c_str(), e.origin.c_str(), e.what());
    std::fclose(out);
    return 1;
  }
  std::fclose(out);

  auto* pm = G4OpticalPhoton::Definition()->GetProcessManager();
  auto* pv = pm->GetPostStepProcessVector(typeDoIt);
  std::printf("mode=%s%s post-step DoIt vector:", mode.c_str(), std::getenv("PT_WLS") ? " (PT_WLS)" : "");
  for (std::size_t i = 0; i < pv->size(); ++i) std::printf(" %s", (*pv)[i]->GetProcessName().c_str());
  std::printf("\nsteps=%ld wrapper-intercepted=%ld boundary statuses:", stepping->steps, stepping->intercepted);
  for (const auto& [k, v] : stepping->statuses) std::printf(" %d:%ld", k, v);
  std::printf("\n");

  bool ok = true;
  for (int st : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14})
    if (!stepping->statuses.count(st)) { std::fprintf(stderr, "FAIL boundary status %d never appeared\n", st); ok = false; }
  if (mode != "stock" && !stepping->fWrapper) { std::fprintf(stderr, "FAIL no wrapper in the process list\n"); ok = false; }
  if (stepping->intercepted != 0) { std::fprintf(stderr, "FAIL the wrapper intercepted %ld steps\n", stepping->intercepted); ok = false; }
  return ok ? 0 : 1;
}
