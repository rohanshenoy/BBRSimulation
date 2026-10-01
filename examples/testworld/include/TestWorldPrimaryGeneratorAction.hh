#ifndef TestWorldPrimaryGeneratorAction_hh
#define TestWorldPrimaryGeneratorAction_hh

#include "BBRPrimarySource.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class TestWorldPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  TestWorldPrimaryGeneratorAction() = default;
  ~TestWorldPrimaryGeneratorAction() override = default;
  void GeneratePrimaries(G4Event* event) override;

 private:
  BBRPrimarySource fSource;
};

#endif
