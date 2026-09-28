# bbrsimLightPipe — a 4 K → mixing-chamber light pipe

## Introduction

`bbrsimLightPipe` models a light pipe that carries thermal radiation from a
4 K stage toward the mixing chamber: a copper tube with a Planck emitter just
upstream of its warm aperture. It links the same BBRsim library as the test
world and carries its own copies of the test world's primary-generator,
stepping, run and action-initialization classes; only the detector
construction and its messenger differ. It is a scaffold: the differential
4 K → MXC measurement it is meant to support is not implemented yet.

## Geometry

`LightPipeDetectorConstruction` places the wall in a 50 cm vacuum world, in one
of two modes selected by `/bbr/lightpipe/mode`:

- **`parametric`** (default) — a tube along +x with its warm aperture at
  x = −50 mm: bore radius 5 mm, length 100 mm, wall thickness 2 mm.
- **`cad`** — an ASCII `.STL` mesh from `/bbr/lightpipe/stlPath`, loaded with
  the bundled header-only CADMesh (`library/include/CADMesh.hh`). The built-in
  reader is ASCII-only; binary STL needs assimp. A sample mesh ships at
  `data/cad/box_sample.stl` (installed as
  `<prefix>/share/BBRsim/data/cad/box_sample.stl`).

The wall is the library's Drude copper at the `/bbr/det/` RRR and stage
temperature, or a perfect reflector.

| Command | Argument | Default |
|---|---|---|
| `/bbr/lightpipe/mode` | `parametric \| cad` | `parametric` |
| `/bbr/lightpipe/bore` | length + unit | 5 mm (inner bore radius) |
| `/bbr/lightpipe/length` | length + unit | 100 mm |
| `/bbr/lightpipe/wallThickness` | length + unit | 2 mm |
| `/bbr/lightpipe/wallMaterial` | `Cu \| reflector` | `Cu` |
| `/bbr/lightpipe/stlPath` | path | none (required in `cad` mode) |

All `/bbr/lightpipe/*` commands are **PreInit only** and are not broadcast to
worker threads: issue them before `/run/initialize`. Errors are fatal `LP001`
(a non-positive dimension), `LP010` (`cad` mode without a path), `LP011` (STL
not found) and `LP012` (STL not parsable).

## Primary event

`LightPipePrimaryGeneratorAction` is the test world's generator: the Planck
emitter by default, or the fixed gun with `/bbr/gun/mode true`. The test-world
default emitter box (1×20×20 mm centred at x = −50 mm) overlaps the tube wall,
so photons would be created inside the copper. `lightpipe.mac` sets
`/bbr/thermal/emitterCenter -51 0 0 mm` and `/bbr/thermal/emitterSize 1 7 7 mm`:
the emitting face sits 0.5 mm upstream of the aperture with its corners inside
the bore. In parametric mode the construction warns (`LP002`) if the configured
emitter reaches into the wall.

## Execution & output

Build against an installed BBRsim (see the [top-level README](../../README.md)):

```bash
. <prefix>/share/BBRsim/bbrsim_env.sh
cmake --preset clang-release          # build/; finds the sourced prefix, else ../../install
cmake --build --preset clang-release
cd build && ./bbrsimLightPipe lightpipe.mac
```

Outside the preset, configure with `-DCMAKE_PREFIX_PATH=<prefix>`. The preset's
install prefix stays `../../install`; use `cmake --install build --prefix
<prefix>` to install the example beside another library prefix.

| Macro | Run |
|---|---|
| `lightpipe.mac` | parametric Cu tube, 4 K Planck emitter sized to the bore, 50 000 events |
| `vis.mac` | visualization setup for the interactive session (no macro argument) |

Each run writes `output/bbr.root` and `output/bbr_legend.json` under the
directory it runs in, with the same ntuples as the test world.

If you `cmake --install` the example into the library's prefix, a shell that has sourced `bbrsim_env.sh` loads the example library from `<prefix>/lib`, because DYLD_LIBRARY_PATH is searched before RPATH. A build-tree binary then runs the installed copy. Re-install after every rebuild, or run the build-tree binary as `env -u DYLD_LIBRARY_PATH ./bbrsimLightPipe <macro>`.

## Testing

`validation/G4Macros/Validation_LightPipe.mac` is a frozen copy of
`lightpipe.mac`; the regression runner checks its output with
`check_no_photons_in_metal.py` and `check_term_status.py`. See
[validation/README.md](../../validation/README.md).

## History

The example classes were renamed from the pre-reorg light-pipe classes in the
reorg, and the four action classes were copied from the test world; follow
their history with `git log --follow -M30% -- <file>`.
