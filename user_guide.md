# BBRsim User Guide

How to run the test world, what every `/bbr/` command does, what the output
contains, and which script checks what. Building, installing and the Python
setup are in the [README](README.md); this guide assumes you have done its
[Quick start](README.md#quick-start).

## Running the test world

`bbrsimTestWorld` runs from its build directory, `examples/testworld/build/`,
in a shell that has sourced the env script. The macros sit beside the binary
there, and output goes to `output/` under the directory you run from.

```bash
./bbrsimTestWorld                 # interactive: opens the UI and runs vis.mac
./bbrsimTestWorld <macro.mac>     # batch
```

The geometry: a 50 cm vacuum world, a 4 mm Cu slab with its front face at
x = 0, and two `vacuum_wg` crack volumes in the slab (crack1: 52 µm gap at
z = 0; crack2: 102 µm gap at z = 3 mm).

In the commands below, `./bbrsimTestWorld` runs in `examples/testworld/build/`
and the `conda run` commands run from the repository root. The light-pipe
example has its own [README](examples/lightpipe/README.md).

## Test cases

### 1. Planck thermal emitter

The default mode. A 1×20×20 mm box at x = −50 mm emits photons with energies
drawn from the Planck photon-number spectrum (10 GHz–20 THz) toward the slab
and the cracks.

```bash
./bbrsimTestWorld planck.mac         # 10 000 events at 4 K
./bbrsimTestWorld planck_5M.mac      # 5 000 000 events at 4 K
./bbrsimTestWorld planck_10K.mac     # 1 000 000 events at 10 K
```

Check the emitted spectrum:

```bash
conda run -n bbrsim python validation/check_planck_spectrum.py examples/testworld/build/output/bbr.root --temp 4
```

The crack2/crack1 entry ratio should match the ratio of their apertures. Both
cracks need entries, so use at least ~1M events from the default emitter:

```bash
./bbrsimTestWorld planck_10K.mac
conda run -n bbrsim python validation/check_crack_ratio.py examples/testworld/build/output/bbr.root
```

The regression runner checks the same ratio on the second run of
`validation/G4Macros/Validation_CrackTransmit.mac`, whose emitter sits just off
the copper face over both cracks, so 200 000 events are enough.

### 2. HFSS crack diffraction

A photon entering a crack volume is handled by the HFSS lookup, which decides
whether it transmits and samples its exit direction and position.
Geometry rule: a `vacuum_wg` volume is entered through one of its two ±x faces
(the crack axis). Entry through a side face is not detected and is treated as
an axial entry. At 500 GHz and normal incidence:

| Crack | Gap | Observed T (40 000 photons) |
|---|---|---|
| crack1 | 52 µm | 50.0 ± 0.25 % |
| crack2 | 102 µm | 50.3 ± 0.25 % |

The expected value is exactly 50 %. A parallel-plate gap narrower than λ/2 is
a perfect polarization filter: only the TEM mode, which has no cutoff,
transmits, so unpolarized light transmits half. The raw HFSS power ratios come
out slightly above 1 (1.0545 for crack1, a port-normalization artefact), so
they are capped at 1 when the tables load; a `[BBR] HFSS ... capped to 1` line
reports it.

To aim the fixed gun at a crack:

```mac
/run/initialize
/bbr/gun/mode true
/bbr/gun/posZ 0.0        # crack1 (3.0 for crack2)
/run/beamOn 10000
```

Plot the exit angles against the HFSS far field:

```bash
conda run -n bbrsim python tools/plot_crack_angular.py examples/testworld/build/output/bbr.root --iwt 180 --iwp 0
```

**One frequency of data.** The lookup gives each photon the HFSS dataset
nearest its frequency (in log frequency) and records the choice in the
`hfss_freq_GHz` column. Only a 500 GHz dataset exists, so every photon of a
broadband Planck run gets the 500 GHz tables, and broadband crack results are
indicative, not quantitative. Adding real data means adding
`<id>_<freq>GHz_Ephi=N` directories under `waveguides/`. There is no
interpolation between frequencies; a photon beyond the grid uses the nearest
edge, with one `BBR008` warning per crack and side.

### 3. Copper reflectance

A fixed gun fires 500 GHz photons at the solid Cu face, away from the cracks,
and counts the fraction absorbed.

```bash
./bbrsimTestWorld reflectance.mac         # OFHC_Cu (RRR 100), 10 000 events
./bbrsimTestWorld reflectance_OF_Cu.mac   # OF_Cu   (RRR 3),    2 000 events
./bbrsimTestWorld reflectance_HP_Cu.mac   # HP_Cu   (RRR 6),    2 000 events
```

Expected absorptance at 500 GHz and 4 K: 4.9×10⁻⁵ for RRR 100 (the Drude
relaxation plateau, so only 0–2 absorptions in 10 000 events), 1.0×10⁻³ for
RRR 3, 6.3×10⁻⁴ for RRR 6. Check the absorbed count against the model
(Poisson test):

```bash
conda run -n bbrsim python validation/check_reflectance.py --root examples/testworld/build/output/bbr.root
```

For a Planck run, compare the absorptance over the whole spectrum with the
Planck-weighted Drude value. The script reads RRR and temperature from the
material name; `--rrr` and `--temp` only override them:

```bash
./bbrsimTestWorld planck_5M.mac
conda run -n bbrsim python validation/check_cu_absorptance.py examples/testworld/build/output/bbr.root
```

Plot the reflectance curves for the three grades:

```bash
conda run -n bbrsim python tools/plot_cu_reflectance.py --out output/cu_reflectance_plots.png
```

Each worker thread also prints a running tally,
`[BBR] reflectance mat=... N=... A_obs=... R_theory=...`, every 1000 Cu hits
on that thread. In a multithreaded run a thread may never reach 1000 hits, so
the line can be missing; add `/run/numberOfThreads 1` before
`/run/initialize` to see it.

## Command reference

All `/bbr/` commands exist from the first macro line. They differ in when they
are allowed and whether they reach the worker threads:

| Commands | Allowed | Broadcast to workers | Notes |
|---|---|---|---|
| `/bbr/dataDir` | before `/run/initialize` | no | runtime data root |
| `/bbr/det/` | before `/run/initialize` | no | geometry is built once, on the master |
| `/bbr/thermal/` | any time | yes | emitter rebuilt at the next event |
| `/bbr/gun/` | any time | yes | read every event |
| `/bbr/config/print` | any time | — | prints every current setting |

There is no runtime geometry change: after `/run/initialize`, changing a
`/bbr/det/` value requires a new session. `config_mt.mac` runs twice in one
session (4 K, then 10 K) to show that settings reach every worker between runs.

### `/bbr/dataDir`

| Command | Argument | Default |
|---|---|---|
| `/bbr/dataDir` | path | `$BBRSIMDATA`, else `<install prefix>/share/BBRsim/data` |

The directory must contain `waveguides/` with the HFSS datasets. Quote a path
that contains spaces.

### `/bbr/det/`: the Cu wall

Three ways to set the copper, all before `/run/initialize`:

```mac
/bbr/det/setCuMaterial OFHC_Cu    # named grade: OFHC_Cu (RRR 100), OF_Cu (RRR 3), HP_Cu (RRR 6); resets T to 4 K
/bbr/det/setCuRRR 250             # any integer RRR >= 1, when you have a measured value
/bbr/det/setCuStageT 40 K         # temperature stage, for warm shield layers
```

For a warm layer, set the temperature, then the RRR. The run confirms the
material at startup:

```
[BBR] Cu wall material: Cu_RRR50_T40K  (RRR=50, T=40 K)
```

RRR is the only property you supply. The room-temperature conductivity,
5.96×10⁷ S/m, is the same for every copper grade; at 4 K the conductivity is
RRR times that. Above about 50 K a phonon term is added (Matthiessen's rule).
The Drude model then gives the reflectance, tabulated from 10 GHz to 20 THz.
The [copper notebook](notebooks/copper_reflectance.ipynb) derives it and
compares it with Serov (2016).

### `/bbr/thermal/`: the Planck emitter

| Command | Argument | Default |
|---|---|---|
| `/bbr/thermal/setT` | value + unit, default K; mK accepted | `4 K` |
| `/bbr/thermal/emitterCenter` | x y z + unit | `-50 0 0 mm` |
| `/bbr/thermal/emitterSize` | full extents Wx Wy Wz + unit | `1 20 20 mm` |

The emitter is a box radiating outward from all six faces. The default suits
the test world; other geometries need their own. The light pipe uses
`-51 0 0 mm` / `1 7 7 mm` to keep the emitter inside its 5 mm bore.

Energies follow the Planck photon-number spectrum, ∝ ν²/(e^{hν/kT} − 1), the
right weighting when each event is one photon. It peaks at hν ≈ 1.59 kT
(133 GHz at 4 K), not at the energy-spectrum peak of 2.82 kT (235 GHz).
The band is fixed at 10 GHz–20 THz. Between 2.2 K and 117 K it holds at least
99 % of the spectrum; outside that range the run prints one `BBR021` warning
(the band misses 13.8 % of the spectrum at 0.5 K and 3.9 % at 150 K).
Directions are uniform in θ over the outward hemisphere, following Chang's
convention; the approximation washes out after a few reflections.

### `/bbr/gun/`: the photon gun

```mac
/bbr/gun/mode true         # true = fixed gun; false = Planck emitter (default)
/bbr/gun/posX -20.0        # position in mm; defaults -20 0 0
/bbr/gun/posY   0.0
/bbr/gun/posZ   0.0        # z = 0: crack1; z = 3: crack2; z > ~5: solid Cu
/bbr/gun/dirX 1.0          # direction, normalized internally; default 1 0 0
/bbr/gun/dirY 0.0
/bbr/gun/dirZ 0.0
/bbr/gun/energy_eV 2.07e-3 # photon energy (2.07e-3 eV = 500 GHz)
/bbr/gun/pol 0 0 0         # 0 0 0 = random polarization (default); any other
                           # vector is projected perpendicular to the direction
```

## Output files

`output/` is relative to the directory the run starts in.

| File | Contents |
|---|---|
| `output/bbr.root` | Two ntuples. **`crossings`**: one row per optical-photon boundary crossing (run and event IDs, position, energy, momentum before and after, incidence angles, volume, material, status and event-type codes, crossing count, and `hfss_freq_GHz`: the HFSS frequency chosen at a crack entry, −1 on every other row). **`abspoints`**: one row per photon death (position, energy, final momentum, reflection count, final volume and status). |
| `output/bbr_legend.json` | Maps the integer code columns (status, event type, volume, material) back to names. |

`abspoints.term_status` is `WorldExit` for a photon that left the world,
`BulkAbsorption` for absorption inside a material, and otherwise the boundary
status that killed it (for example `BBRAbsorb`).

Read the output with the `bbrsim` package, which decodes the codes:

```python
from bbrsim import io
crossings, abspoints = io.load("examples/testworld/build/output/bbr.root")
```

**Selecting crack entries.** Every crossing is logged from the world side:
`mat_pre` and `vol_pre` are always `G4_Galactic` and `World`, and the entered
material and volume are `mat_post` and `vol_post`. So a crack entry is
`mat_post == "vacuum_wg"`, split by `vol_post`, and a first copper hit is a
`mat_post` starting with `Cu_RRR`.

## Scripts

Run every script as `conda run -n bbrsim python <script>`. Each reads BBRsim
output (default `output/bbr.root` in the current directory) through the
`bbrsim` package.

**Validators** (`validation/`), run by the regression runner:

| Script | Checks |
|---|---|
| `check_physics.py` | the `bbrsim` formulas against reference values (no input) |
| `check_reflectance.py` | absorbed count vs the Drude model (`--root`, `--RRR`, `--T_K`, `--freq`) |
| `check_cu_serov.py` | Drude loss for `OF_Cu` / `HP_Cu` vs Serov (2016) within ±10 %; `HP_Cu` is the known expected failure |
| `check_planck_spectrum.py` | the emitted spectrum's peak (`--temp`) |
| `check_nreflect.py` | the per-photon reflection-count distribution |
| `check_angle_distribution.py` | incidence angles at the copper vs uniform-in-θ emission (KS test, sized for the 10 000-event Planck run) |
| `check_crack_transmittance.py` | T = 0.50 at normal incidence, no exits along the crack face |
| `check_crack_ratio.py` | crack2/crack1 entry ratio vs the aperture ratio, within 3 σ |
| `check_crack_wall_reflection.py` | reflection off a crack's side wall flips only p_z |
| `check_crack_oblique.py` | 45° incidence: 138 checks against the HFSS tables and the Python model |
| `check_crack_frequency.py` | dataset choice per photon on the mock five-frequency tree: 89 checks (`--data-dir`) |
| `check_invariants.py` | no photon ever travels inside a metal, and every photon death is labelled correctly (`--allow-no-crossings` for a run that crosses no boundary) |

**Validator run by hand**, because it needs a large run:
`check_cu_absorptance.py` (`planck_5M.mac`), shown above.

**Plot scripts** (`tools/`): `plot_cu_reflectance.py` (writes to the current
directory, or `--out`), `plot_crack_angular.py` (`--iwt`, `--iwp`; writes
next to its input), `plot_crack_frequency.py` (writes into its output
directory, or `--out`), `plot_test_output.py` (overview of one output file;
writes next to its input).

The fixtures these run on, and the runner, are described in
[validation/README.md](validation/README.md).

## Troubleshooting

**`[BBR] det/setCuMaterial: unknown alias`**: only `OFHC_Cu`, `OF_Cu` and
`HP_Cu` exist. Use `/bbr/det/setCuRRR <N>` for any other grade.

**`[BBR] det/setCuRRR: RRR must be >= 1`**: RRR is a positive integer.

**A crack macro aborts with `BBR011`**: the HFSS data were not found. Run
`cmake --install build` in the repository, or point `/bbr/dataDir` or
`$BBRSIMDATA` at a directory containing `waveguides/`. The message names the
path it tried.

**No `[BBR] reflectance` lines**: expected with many threads; see
[Copper reflectance](#3-copper-reflectance).

**`check_cu_absorptance.py` fails** (A_obs / A_theory outside [0.3, 3]): too
few events (at RRR 100 a 10 000-event run gives 0–2 absorptions), or an
`--rrr` / `--temp` override that does not match the run. It takes the path to
a ROOT file.

**`ModuleNotFoundError` in a Python script**: run it with
`conda run -n bbrsim python`, not `python3`, and make `bbrsim` importable:
source the env script, or run `conda run -n bbrsim pip install -e tools/python`
once.
