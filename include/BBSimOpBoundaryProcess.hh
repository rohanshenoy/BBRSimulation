#ifndef BBSimOpBoundaryProcess_hh
#define BBSimOpBoundaryProcess_hh

#include "G4OpBoundaryProcess.hh"
#include "G4WrapperProcess.hh"

// Wrapper around G4OpBoundaryProcess. PostStepDoIt is where custom boundary
// physics (diffraction, vacuum->copper) will be injected. Currently a pure
// pass-through — results must be identical to the unwrapped simulation.
class BBSimOpBoundaryProcess : public G4WrapperProcess
{
 public:
  explicit BBSimOpBoundaryProcess(
    const G4String& name = "BBSimOpBoundary");
  ~BBSimOpBoundaryProcess() override = default;

  // Core override — pure pass-through for now; conditions added later.
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;

  // Call-throughs for G4OpBoundaryProcess-specific API used by SteppingAction.
  G4OpBoundaryProcessStatus GetStatus() const;
  G4bool                    GetInvokeSD() const;

 private:
  // Static cast is safe: we always wrap exactly one G4OpBoundaryProcess.
  G4OpBoundaryProcess* GetWrappedProcess() const;
};

#endif
