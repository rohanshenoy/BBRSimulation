#ifndef BBSimOpBoundaryProcess_hh
#define BBSimOpBoundaryProcess_hh

#include "G4OpBoundaryProcess.hh"
#include "G4ParticleChange.hh"
#include "G4WrapperProcess.hh"
#include "globals.hh"
#include <atomic>
#include <utility>

// Wrapper around G4OpBoundaryProcess. PostStepDoIt intercepts photons that
// cross into a volume filled with material "vacuum_wg" and routes them through
// the HFSS diffraction model. The dataset is selected by (volume name, photon
// frequency): the volume name is the HFSS dataset ID and BBRCrackLibrary picks
// the grid frequency nearest in log space. Everything else is a pure
// pass-through.
//
// Geometry rule: a `vacuum_wg` volume is entered through one of its two ±x
// faces (the crack axis). Entry through a side face is not detected and is
// treated as an axial entry.
class BBSimOpBoundaryProcess : public G4WrapperProcess
{
 public:
  explicit BBSimOpBoundaryProcess(const G4String& name = "BBSimOpBoundary");
  ~BBSimOpBoundaryProcess() override = default;

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;

  // Forwards to the wrapped process, then runs
  // BBRCrackLibrary::ValidatePlacedCracks() (once per process), so a missing or
  // wrong HFSS dataset, or a crack too small for its exit grid, stops the run at
  // the first run initialization (/run/initialize with an MT or task run
  // manager, whose Initialize calls BeamOn(0); /run/beamOn with a sequential
  // one), before any event. A geometry without vacuum_wg volumes sees nothing
  // new.
  void BuildPhysicsTable(const G4ParticleDefinition& particle) override;

  // What the wrapper itself did on the last PostStepDoIt invocation.
  // kBBRNone means the step was passed through to the stock process, whose
  // own GetStatus() is then the authoritative status.
  enum BBRBoundaryStatus {
    kBBRNone,
    kBBRDiffractionTransmit,
    kBBRDiffractionReflect,
    kBBRReflect,
    kBBRAbsorb
  };
  BBRBoundaryStatus GetLastBBRStatus() const { return fLastBBRStatus; }
  G4String GetLastBBRStatusString() const;

  // HFSS grid frequency [GHz] used by the last intercepted crack crossing, or
  // -1 if the last PostStepDoIt did not run the diffraction handler. Read by
  // an example's stepping action for the crossings column hfss_freq_GHz.
  G4double GetLastHFSSFrequencyGHz() const { return fLastHFSSFreqGHz; }

  // Transmissions whose start on the ray through the sampled exit point was not valid and
  // whose start at the normal inset was not valid either, so the photon started at the
  // normal inset anyway (see ExitStart). Also printed as exit_start_fallback= on the
  // "[BBR] diffraction events=" line.
  G4int GetExitStartFallbackCount() const { return fNExitFallback; }

  // Process-wide totals over every thread since the last reset: HFSS transmissions and, of those,
  // exit-start fallbacks. The run action of the master resets them at the start of a run and
  // prints them at its end (BBRAnalysis), so the count is one number per run, not one per worker.
  static void ResetRunTotals();
  static std::pair<G4long, G4long> RunTotals();   // {transmissions, fallbacks}

  G4OpBoundaryProcessStatus GetStatus() const;
  // Forward SetInvokeSD to the wrapped G4OpBoundaryProcess (no getter in G4 API).
  void SetInvokeSD(G4bool flag);
  // Public accessor so external observers (e.g. an example's stepping action) can
  // query the inner process status directly.
  G4OpBoundaryProcess* GetWrappedProcess() const;

 private:
  G4VParticleChange* HandleDiffractionBoundary(const G4Track& aTrack,
                                               const G4Step& aStep);
  G4VParticleChange* HandleReflectanceBoundary(const G4Track& aTrack,
                                               const G4Step& aStep);

  // Where a transmitted photon starts inside the crack, at normal distance inset from the exit
  // plane. xf is the touchable's global-to-local transform, pos_out the sampled exit point on the
  // exit plane, dir_out the sampled direction and normal_hat the exit face's outward axis. First
  // choice: on the ray through pos_out, pos_out - (inset/n) dir_out with n = dir_out.normal_hat, so
  // the photon crosses the exit plane at its sampled point. Fallback: the normal inset
  // pos_out - inset normal_hat. A start is valid if it is inside the solid and its ray reaches the
  // exit plane before leaving the solid (DistanceToOut, to the surface tolerance). If neither is
  // valid the normal inset is returned with fellBack = true.
  G4ThreeVector ExitStart(const G4VTouchable& touch, const G4AffineTransform& xf,
                          const G4ThreeVector& pos_out, const G4ThreeVector& dir_out,
                          const G4ThreeVector& normal_hat, G4double inset,
                          G4bool& fellBack) const;

  G4ParticleChange fParticleChange;

  BBRBoundaryStatus fLastBBRStatus = kBBRNone;
  G4double fLastHFSSFreqGHz = -1.;

  // Per-instance counters (safe in MT — each thread gets its own process clone).
  G4int fNDiffraction         = 0;
  G4int fNDiffractionTransmit = 0;
  G4int fNReflectance         = 0;
  G4int fNReflectanceAbsorb   = 0;
  G4int fNExitFallback        = 0;

  static std::atomic<G4long> sRunTransmissions;
  static std::atomic<G4long> sRunExitFallbacks;
};

#endif
