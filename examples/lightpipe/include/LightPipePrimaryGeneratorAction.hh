#ifndef LightPipePrimaryGeneratorAction_hh
#define LightPipePrimaryGeneratorAction_hh

#include "BBRPrimarySource.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class LightPipePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  LightPipePrimaryGeneratorAction() = default;
  ~LightPipePrimaryGeneratorAction() override = default;
  void GeneratePrimaries(G4Event* event) override;

 private:
  BBRPrimarySource fSource;
};

#endif
