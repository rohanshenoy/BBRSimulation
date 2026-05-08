#ifndef BBRReflectanceDetectorConstruction_hh
#define BBRReflectanceDetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"

// World: 100 mm vacuum_wg cube.
// Daughter: OFHC_Cu slab (halfX=1 mm, halfY=halfZ=5 mm), centre at x=+1 mm.
// Front Cu face is at x=0.  Gun (BBRDiffractionPGA default) fires from
// (−20 mm, 0, 0) along +x → hits Cu at x=0 at normal incidence.
class BBRReflectanceDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  G4VPhysicalVolume* Construct() override;
};

#endif
