# OpNovice2 reference (not built)

These are the **unmodified** Geant4 `OpNovice2` example sources, macros, and
expected output, kept here as a known-good reference for stock optical-photon
boundary behaviour. BBRsim was originally scaffolded from this example.

**Nothing in this directory is compiled into `BBRSim`.** The build only globs
`src/*.cc` and `include/*.hh` at the repository root, so these files are inert.
They are retained purely for reading and comparison.

The canonical, authoritative copy ships with every Geant4 install at:

```
$G4INSTALL/examples/extended/optical/OpNovice2
```

Prefer that copy when you need a guaranteed-current reference; this snapshot
reflects the version BBRsim was forked from.

## Contents

- `OpNovice2.cc` — original entry point (separate `main`).
- `src/`, `include/` — stock `DetectorConstruction`, `PrimaryGeneratorAction`,
  `SteppingAction`, `Run`, `RunAction`, `HistoManager`, `TrackingAction`,
  `TrackInformation`, `ActionInitialization`, and the three messenger classes.
- `macros/` — stock test macros (`boundary.mac`, `fresnel.mac`, `wls.mac`, …)
  and their reference `.out` files.

## Note on wrapper pass-through verification

`BBSimOpBoundaryProcess` falls through to the stock `G4OpBoundaryProcess` for
any geometry without `vacuum_wg` volumes or `REFLECTIVITY` materials. These
OpNovice2 sources are *not* wired into the `BBRSim` executable, so they cannot
be used as a live pass-through regression as-is. A future automated regression
should use a purpose-built minimal geometry exercising stock boundary optics
rather than resurrecting this example.
