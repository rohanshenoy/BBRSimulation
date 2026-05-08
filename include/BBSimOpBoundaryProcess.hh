#ifndef BBSimOpBoundaryProcess_hh
#define BBSimOpBoundaryProcess_hh

#include "G4OpBoundaryProcess.hh"
#include "G4ParticleChange.hh"
#include "G4WrapperProcess.hh"

// Wrapper around G4OpBoundaryProcess. PostStepDoIt intercepts photons that
// cross into a volume filled with material "vacuum_wg" and routes them through
// the HFSS diffraction model. The volume name is the HFSS dataset ID used to
// look up the correct CSV data via BBRCrackLibrary. Everything else is a
// pure pass-through.
class BBSimOpBoundaryProcess : public G4WrapperProcess
{
 public:
  explicit BBSimOpBoundaryProcess(const G4String& name = "BBSimOpBoundary");
  ~BBSimOpBoundaryProcess() override = default;

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;

  G4OpBoundaryProcessStatus GetStatus() const;
  // Forward SetInvokeSD to the wrapped G4OpBoundaryProcess (no getter in G4 API).
  void SetInvokeSD(G4bool flag);

 private:
  G4VParticleChange* HandleDiffractionBoundary(const G4Track& aTrack,
                                               const G4Step& aStep);
  G4VParticleChange* HandleReflectanceBoundary(const G4Track& aTrack,
                                               const G4Step& aStep);

  G4OpBoundaryProcess* GetWrappedProcess() const;

  G4ParticleChange fParticleChange;
};

#endif
