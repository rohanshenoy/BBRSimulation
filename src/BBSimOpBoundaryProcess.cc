#include "BBSimOpBoundaryProcess.hh"

BBSimOpBoundaryProcess::BBSimOpBoundaryProcess(const G4String& name)
  : G4WrapperProcess(name)
{}

// Pure pass-through. Future structure:
//   if (condition1: vacuum->copper)    return HandleVacuumToCopperBoundary(...)
//   else if (condition2: diffraction)  return HandleDiffractionBoundary(...)
//   else                               return pRegProcess->PostStepDoIt(...)
G4VParticleChange* BBSimOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack,
                                                        const G4Step& aStep)
{
  return pRegProcess->PostStepDoIt(aTrack, aStep);
}

G4OpBoundaryProcessStatus BBSimOpBoundaryProcess::GetStatus() const
{
  return GetWrappedProcess()->GetStatus();
}

G4OpBoundaryProcess* BBSimOpBoundaryProcess::GetWrappedProcess() const
{
  return static_cast<G4OpBoundaryProcess*>(pRegProcess);
}
