# bbrsimTestWorld — the BBRsim test world

## Introduction

`bbrsimTestWorld` is the reference application for the BBRsim library: two
copper plates with a machining-gap crack each, lit by thermal (Planck) radiation
or by a fixed photon gun. It exercises the BBR physics the library provides:
tabulated full-Drude Cu reflectance at the plates, HFSS-tabulated transmission through the
cracks, and the Planck emitter. Most of the validation fixtures run it.

## Geometry

`TestWorldDetectorConstruction`: a 50 cm vacuum world holding two touching Cu
plates, each as thick as its crack, so that every crack volume equals its HFSS
model exactly. `CuPlateCrack1` is 1.0 mm thick (z from −25 to 1.5 mm) and
`CuPlateCrack2` is 1.5 mm thick (z from 1.5 to 25 mm); both are 50 mm wide in y
and have their front faces at x = 0. The plate material comes from
`BBRMaterials::GetCopper(RRR, T)`, set with `/bbr/det/setCuMaterial`,
`/bbr/det/setCuRRR` and `/bbr/det/setCuStageT` before `/run/initialize`
(default OFHC_Cu, RRR 100, 4 K). Each plate holds one crack daughter, a
`vacuum_wg` volume named by its HFSS dataset ID and filling a slot through the
plate:
`InfParallelPlate_crack1Rohan` (1.0 × 10 × 0.050 mm along the propagation,
long and gap directions, centred at z = 0) and
`InfParallelPlate_crack2` (1.5 × 10 × 0.100 mm, centred at z = 3 mm). For each,
the library reads `<dataDir>/waveguides/<id>_<freq>GHz_Ephi={0,1}`. Since
OR-020 the solids equal the HFSS models in length and section (before, the
slots were 4 mm long, 10.2 mm wide and 2 µm taller than HFSS); at load, exit
points that HFSS placed on a wall are clamped 1 nm inside the declared section
(see *Dataset sidecars* in `validation/README.md`). `/bbr/testworld/cracks false`
(PreInit only, default `true`) leaves out both plates and their cracks, for a
crack-free world.

Opt-in round gap: `/bbr/testworld/roundGap true` (PreInit only, default
`false`) also places `RoundGapPlate`, a 0.4 mm Cu plate centred at z = −80 mm
with its front face at x = 0, holding `RoundGap_r50um`, a straight round
`vacuum_wg` hole along x (a `G4Tubs` turned onto local x inside a
`G4DisplacedSolid`; radius 50 µm, equal to the HFSS radius). With the
default `false` the world, and so the pinned regression macros, are unchanged.
A run with the round gap on needs a data root holding a `RoundGap_r50um`
dataset, for example the mock written by
`validation/Scripts/make_mock_round_gap.py`, selected with `/bbr/dataDir`.

Geometry rule: a `vacuum_wg` volume is entered through one of its two ±x faces
(the crack axis). Entry through a side face is not detected and is treated as
an axial entry.

## Primary event

`TestWorldPrimaryGeneratorAction` has two modes:

- **Planck emitter** (default): a `ThermalSurface` box, 1×20×20 mm centred at
  x = −50 mm by default (`/bbr/thermal/emitterCenter`, `/bbr/thermal/emitterSize`),
  emitting outward from all six faces. Energies follow the Planck photon-number
  spectrum at `/bbr/thermal/setT` (default 4 K) over 10 GHz–20 THz; directions
  are uniform in θ over the outward hemisphere.
- **Fixed gun** (`/bbr/gun/mode true`): one photon per event from
  `/bbr/gun/posX|posY|posZ` along `/bbr/gun/dirX|dirY|dirZ` at
  `/bbr/gun/energy_eV` (default 500 GHz), with random polarization unless
  `/bbr/gun/pol` sets one. z = 0 aims at crack1, z = 3 mm at crack2, z ≳ 5 mm at
  solid copper.

## Commands

All are issued before `/run/initialize` unless noted.

| Command | Owner | Effect |
|---|---|---|
| `/bbr/det/setCuMaterial`, `/bbr/det/setCuRRR`, `/bbr/det/setCuStageT` | library (`BBRConfigMessenger`) | Cu RRR and temperature stage |
| `/bbr/dataDir` | library | data root holding `waveguides/` |
| `/bbr/thermal/*`, `/bbr/gun/*` | library | Planck emitter and fixed gun; also valid after `/run/initialize` |
| `/bbr/testworld/roundGap <bool>` | this example (`TestWorldMessenger`) | also place the round gap (default `false`); refused after `/run/initialize` |
| `/bbr/testworld/cracks <bool>` | this example (`TestWorldMessenger`) | place the two crack plates (default `true`); `false` gives a crack-free world, the round gap unaffected; refused after `/run/initialize` |

## Execution & output

Build against an installed BBRsim (see the [top-level README](../../README.md)):

```bash
. <prefix>/share/BBRsim/bbrsim_env.sh
cmake --preset clang-release          # build/; finds the sourced prefix, else ../../install
cmake --build --preset clang-release
cd build && ./bbrsimTestWorld planck.mac
```

Outside the preset, configure with `-DCMAKE_PREFIX_PATH=<prefix>`. The preset's
install prefix stays `../../install`; use `cmake --install build --prefix
<prefix>` to install the example beside another library prefix. With no
macro argument the binary opens an interactive session and runs `vis.mac`. The
macros live in `G4Macros/` and are copied beside the binary:

To build this example with the current source library, configure from the
repository root with `-DBUILD_BBRSIM_EXAMPLES=ON`; the binary is then
`build/examples/testworld/bbrsimTestWorld`. Set `BBRSIMDATA` to the repository's
`data/` directory when running it without installing. A batch macro error,
including one in a nested macro, makes the executable exit nonzero.

| Macro | Run |
|---|---|
| `planck.mac` | Planck emitter at 4 K, 10 000 events |
| `planck_10K.mac` | Planck emitter at 10 K, 1 000 000 events |
| `planck_5M.mac` | Planck emitter at 4 K, 5 000 000 events (production) |
| `planck_50M.mac` | Planck emitter at 4 K, 50 000 000 events |
| `reflectance.mac` | fixed gun into OFHC_Cu (RRR 100) at z = 10 mm, 500 GHz, 10 000 events |
| `reflectance_OF_Cu.mac` | the same into OF_Cu (RRR 3) at z = 5 mm, 2 000 events |
| `reflectance_HP_Cu.mac` | the same into HP_Cu (RRR 6) at z = 5 mm, 2 000 events |
| `config_mt.mac` | two runs in one session (4 K, then 10 K), one file each (`output/bbr_mt_r0.root`, `output/bbr_mt_r1.root`): the worker config clones follow the broadcast |
| `vis.mac` | visualization setup for the interactive session |

Each run writes `output/bbr.root` (ntuples `crossings` and `abspoints`) and
`output/bbr.metadata.json` under the directory it runs in; a later `/run/beamOn`
overwrites them unless the macro names a new file with `/analysis/setFileName`.
Read them with the `bbrsim` Python package (`from bbrsim import io`).

If you `cmake --install` the example into the library's prefix, a shell that has sourced `bbrsim_env.sh` loads the example library from `<prefix>/lib`, because DYLD_LIBRARY_PATH is searched before RPATH. A build-tree binary then runs the installed copy. Re-install after every rebuild, or run the build-tree binary as `env -u DYLD_LIBRARY_PATH ./bbrsimTestWorld <macro>`.

## Testing

The regression runner runs this executable on three of the macros here,
`reflectance.mac`, `planck.mac` and `config_mt.mac`, and on six fixtures that
exist only in `validation/G4Macros/` (`Validation_CrackWall`,
`Validation_WorldExit`, `Validation_CrackTransmit`, `Validation_CrackOblique`,
`Validation_CrackFrequency`, `Validation_RoundGap`). The three macros are therefore regression
inputs: the drift guard `regression macros pinned` fails when one of their
command lines changes (only full-line comment and blank-line edits
pass), so change them only together with their validators and then update the
pin in `validation/Scripts/drift_guards.sh`. Run everything with
`validation/Scripts/run_regression.sh`; see
[validation/README.md](../../validation/README.md).

## History

The example classes were renamed from the pre-reorg test-world classes in the
reorg; follow their history with `git log --follow -M30% -- <file>`.
