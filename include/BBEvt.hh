#ifndef BBEvt_h
#define BBEvt_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

class BBEvt {
 public:
  BBEvt();
  virtual ~BBEvt();

  G4ThreeVector position;
  G4ThreeVector direction;
  G4double      energy;    // raw eV — NOT Geant4 internal units
};

#endif
