#ifndef BBSimOpBoundaryProcess_hh
#define BBSimOpBoundaryProcess_hh

#include "G4OpBoundaryProcess.hh"
#include "G4ParticleChange.hh"
#include "G4WrapperProcess.hh"

// Wrapper around G4OpBoundaryProcess. PostStepDoIt intercepts photons that
// cross into a volume filled with material "vacuum_wg" and routes them through
// the HFSS diffraction model. The dataset is selected by (volume name, photon
// frequency): the volume name is the HFSS dataset ID and BBRCrackLibrary picks
// the grid frequency nearest in log space. Everything else is a pure
// pass-through.
class BBSimOpBoundaryProcess : public G4WrapperProcess
{
 public:
  explicit BBSimOpBoundaryProcess(const G4String& name = "BBSimOpBoundary");
  ~BBSimOpBoundaryProcess() override = default;

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;

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

  G4ParticleChange fParticleChange;

  BBRBoundaryStatus fLastBBRStatus = kBBRNone;
  G4double fLastHFSSFreqGHz = -1.;

  // Per-instance counters (safe in MT — each thread gets its own process clone).
  G4int fNDiffraction         = 0;
  G4int fNDiffractionTransmit = 0;
  G4int fNReflectance         = 0;
  G4int fNReflectanceAbsorb   = 0;
};

#endif
