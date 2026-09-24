# BBRsim User Guide

This guide covers everything needed to build BBRsim, run its existing test
cases, interpret the output, and configure material properties.

---

## Prerequisites

| Requirement | Version | Notes |
|---|---|---|
| Geant4 | 11.x | Must be built with optical physics (`-DGEANT4_USE_OPENGL_X11=ON` optional) |
| CMake | ≥ 3.16 | |
| C++ compiler | C++17 | Clang or GCC |
| Python | 3.x via conda | `bbrsim` environment — see below |

**Python environment setup** (one-time):

```bash
conda create -n bbrsim python=3.11 numpy scipy matplotlib pandas
conda activate bbrsim
```

All Python scripts must be run as `conda run -n bbrsim python <script>` — not
`python3` directly, even if packages appear installed in the base environment.

---

## Build

```bash
cd BBRSimulation
cmake --preset clang-release          # Apple Clang, Release, binaryDir build/
cmake --build --preset clang-release
```

`CMakePresets.json` pins `/usr/bin/clang`/`clang++` and `Release`; IDE CMake
integrations (VSCode CMake Tools) pick the preset up automatically instead of
offering their own compiler kits. The manual equivalent is
`cmake -S . -B build -DCMAKE_C_COMPILER=/usr/bin/clang
-DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DCMAKE_BUILD_TYPE=Release`; a
`clang-debug` preset builds into `build-debug/`.

On macOS, always configure with Apple Clang explicitly. A bare `cmake ..` can
select Homebrew GCC (libstdc++); it compiles, then fails at link against an
Apple-Clang (libc++) Geant4 install with undefined symbols on every API taking
`std::` types. That is a compiler/ABI mismatch, not a code bug — check
`CMAKE_CXX_COMPILER` in `build/CMakeCache.txt`.

Batch-only build (no UI or visualization — faster, no display required):

```bash
cmake -DWITH_GEANT4_UIVIS=OFF ..
make
```

CMake copies all `.mac` files to `build/` alongside the `BBRSim` executable.

---

## Running BBRsim

The executable is `BBRSim`, always run from the `build/` directory.

**Interactive** (requires UI + visualization in your Geant4 build):

```bash
./BBRSim
```

**Batch** (pass a macro file):

```bash
./BBRSim <macro.mac>
```

---

## Test cases

All test cases run in the consolidated test world
(`BBRTestDetectorConstruction`): a 50 cm vacuum world, a 4 mm Cu slab with
its front face at x = 0, and two `vacuum_wg` crack daughters (crack1: 52 µm
gap at z = 0; crack2: 102 µm gap at z = 3 mm). Every optical-photon boundary
crossing is logged to the `crossings` ntuple in `output/bbr.root` (see
*Output files*).

### 1. Planck thermal emitter

The default mode (`/bbr/gun/mode false`): a 1×20×20 mm box surface at
x = −50 mm emits photons with energies drawn from the Planck photon-number
spectrum (10 GHz–20 THz) toward the Cu slab and cracks.

```bash
./BBRSim planck.mac        # 10 000 events at 4 K
./BBRSim test.mac          # 5 000 000 events at 4 K
./BBRSim test_10K.mac      # 1 000 000 events at 10 K
```

Validate the emitted spectrum:

```bash
conda run -n bbrsim python scripts/check_planck_spectrum.py build/output/bbr.root --temp 4
```

### 2. HFSS crack diffraction

Photons entering a `vacuum_wg` crack volume are routed through the HFSS
S-parameter lookup. The `[BBR] diffraction` stdout line reports the running
observed transmittance. At 500 GHz normal incidence:

| Crack | Gap | Observed T |
|---|---|---|
| crack1 | 52 µm | 50.0 ± 0.25% |
| crack2 | 102 µm | 50.3 ± 0.25% |

The benchmark is **50%**: a parallel-plate gap thinner than λ/2 is a perfect
polarization filter — only the cutoff-free TEM mode transmits — so an
unpolarized beam transmits exactly half. Both observations (40k photons each;
`crack_transmit.mac` + `check_crack_transmittance.py`) are consistent with that
ideal. Raw HFSS power ratios above 1 — a port-normalization artefact, 1.0545
for crack1 — are capped at 1 when the tables load (a `[BBR] HFSS ... capped to
1` line is printed). Before that cap the observed values were 51.9% / 50.5%,
and an older "52.7%" figure in archived plans carried the same artefact; do not
quote either.

> **Scope limit — one HFSS frequency of data.** The lookup itself is keyed by
> frequency: `BBRCrackLibrary` reads the grid from the dataset directory names
> (`<id>_<freq>GHz_Ephi=N` under `<dataDir>/waveguides`) and gives each photon
> the grid point nearest in log frequency, recording it in `hfss_freq_GHz`.
> Only 500 GHz exists, so Planck-mode runs spanning 10 GHz–20 THz still apply
> the 500 GHz transmittance and angular PDFs to every photon, and broadband
> crack results remain indicative rather than quantitative. Adding real exports
> is a matter of adding directories; there is no interpolation between grid
> points, and a photon beyond the grid uses the nearest edge (one `BBR008`
> warning per crack and side; a one-point grid never warns). The
> quarter-symmetry azimuth fold/unfold is validated at normal incidence and at
> 45° off normal (`crack_oblique.mac`), the only oblique row in the HFSS data.

To aim the fixed gun at a crack:

```mac
/run/initialize
/bbr/gun/mode true
/bbr/gun/posZ 0.0        # crack1 (use 3.0 for crack2)
/run/beamOn 10000
```

Crack analyses read the `crossings` ntuple from `output/bbr.root`:

```bash
conda run -n bbrsim python scripts/check_crack_ratio.py build/output/bbr.root
conda run -n bbrsim python scripts/plot_crack_angular.py build/output/bbr.root --iwt 180 --iwp 0
```

Crack entries are selected as `mat_post == "vacuum_wg"` and split by `vol_post`.
Every crossing is logged World-side, so `mat_pre`/`vol_pre` are always
`G4_Galactic`/`World`; a `mat_pre == "vacuum_wg"` filter matches zero rows.

### 3. Copper reflectance

Fires 500 GHz photons at the solid Cu wall (away from the cracks) and
measures the fraction absorbed. The `[BBR] reflectance` output line is
printed every 1000 events.

```bash
./BBRSim reflectance.mac   # OFHC_Cu (RRR=100), 10 000 events at z=10 mm
./BBRSim test_of_cu.mac    # OF_Cu   (RRR=3),   2 000 events at z=5 mm
./BBRSim test_hp_cu.mac    # HP_Cu   (RRR=6),   2 000 events at z=5 mm
```

The output line format is:

```
[BBR] reflectance mat=Cu_RRR100_T4K N=10000 A_obs=... R_theory=0.999951
```

- `A_obs` = fraction of photons absorbed (1 − R_obs)
- `R_theory` = full complex Drude reflectance at the gun frequency

With the full complex Drude model, OFHC Cu at 4 K sits on the relaxation
plateau: D ≈ 4.9×10⁻⁵ at 500 GHz, i.e. only ~0–2 absorptions per 10 000
events. Lower grades absorb more (OF_Cu RRR=3: D ≈ 1.0×10⁻³; HP_Cu RRR=6:
D ≈ 6.3×10⁻⁴). Statistical check (Poisson test on the absorbed count):

```bash
./BBRSim reflectance.mac
conda run -n bbrsim python scripts/check_reflectance.py
```

For Planck-mode runs, compare the aggregate absorptance against the
Planck-weighted Drude integral. The script reads the ROOT file and computes
absorptance from decoded crossings; it defaults to `build/output/bbr.root` and
infers RRR and temperature from the run, so `--rrr` / `--temp` are overrides
only:

```bash
./BBRSim test.mac
conda run -n bbrsim python scripts/check_cu_absorptance.py build/output/bbr.root
```

Plot reflectance curves for all three grades:

```bash
conda run -n bbrsim python scripts/plot_cu_reflectance.py
# Saves: build/cu_reflectance_plots.png
```

---

## Configuring the Cu wall material

The Cu material is set by messenger commands, all issued **before**
`/run/initialize`. There are three ways to specify it.

### Option A — Named alias

```mac
/bbr/det/setCuMaterial OFHC_Cu    # RRR=100, 4 K
/bbr/det/setCuMaterial OF_Cu      # RRR=3,   4 K
/bbr/det/setCuMaterial HP_Cu      # RRR=6,   4 K
/run/initialize
```

### Option B — Direct RRR (recommended when you have a measured spec)

```mac
/bbr/det/setCuRRR 250             # any integer >= 1
/run/initialize
```

### Option C — RRR + temperature stage (warm shield layers)

```mac
/bbr/det/setCuStageT 40 K         # set temperature first
/bbr/det/setCuRRR 50
/run/initialize
```

The simulation confirms the resolved material at startup:

```
[BBR] Cu wall material: Cu_RRR50_T40K  (RRR=50, T=40 K)
```

### Physics of the RRR → σ mapping

σ_RT = 5.96×10⁷ S/m is universal for all copper grades (you never set this).
At 4 K the impurity term dominates and σ_DC = RRR × σ_RT is an excellent
approximation. At higher temperatures, Matthiessen's rule adds a phonon
contribution. The full complex Drude model (σ(ω) = σ_DC/(1−iωτ) inserted into
ε̃ = 1 + iσ/(ε₀ω)) then gives R(ω), tabulated over 10 GHz–20 THz to match the
Planck emitter range. See
[docs/physics/copper_reflectance_model.md](docs/physics/copper_reflectance_model.md)
for the derivation.

---

## Command availability (PreInit vs Idle)

All `/bbr/...` commands are owned by `BBRConfigMessenger`, which
`BBRConfigManager` registers at startup, so every command exists from the first
macro line. What differs is *when* each is allowed and whether it reaches worker
threads:

| Directory | Valid states | Broadcast to workers? | Notes |
|---|---|---|---|
| `/bbr/dataDir` | `PreInit` **only** | No | Root of the runtime data tree (must contain `waveguides/`); default `$BBRSIMDATA`, else `../data` |
| `/bbr/det/` | `PreInit` **only** | No | Geometry is built on the master in `Construct()`; issue before `/run/initialize` |
| `/bbr/gun/` | `PreInit` and `Idle` | Yes | Read fresh each event |
| `/bbr/thermal/` | `PreInit` and `Idle` | Yes | Planck CDF / emitter box rebuilt on the next event |
| `/bbr/config/print` | `PreInit` and `Idle` | — | Dumps all current settings |

`BBRConfigManager::Instance()` is a thread-local clone: the master builds from
compiled defaults and each worker copy-constructs from the master, so broadcast
settings propagate without shared mutable state. `verify_config_mt.mac`
exercises this with two runs (4 K then 10 K) in one session.

There is **no runtime geometry reinitialization**. `/run/reinitializeGeometry`
is not supported by any BBRsim detector construction; changing a `/bbr/det/`
value after `/run/initialize` requires a new session.

---

## Planck emitter configuration

```mac
/bbr/thermal/setT 10.0                    # emitter temperature [K], default 4.0
/bbr/thermal/emitterCenter -50 0 0 mm     # emitter box centre, world frame (default)
/bbr/thermal/emitterSize   1 20 20 mm     # emitter box FULL extents Wx Wy Wz (default)
/run/initialize
/run/beamOn 10000
```

The emitter is a box radiating outward from all six faces; its centre and
extents are read from `BBRConfigManager` and the box is rebuilt on the next
event when either changes. The defaults suit the `BBRSim` test world (a patch
50 mm in front of the Cu slab). Other geometries must size it themselves — the
light pipe uses `-51 0 0 mm` / `1 7 7 mm` so the emitting face sits 0.5 mm
upstream of the aperture with its corners inside the 5 mm bore (see
`lightpipe.mac`); the default box would put the emitting face inside the tube
wall.

Photon energies are drawn from the Planck **photon-number** spectrum
∝ ν²/(e^{hν/kT}−1) over a fixed 10 GHz–20 THz range — the correct weighting for
an unweighted photon Monte Carlo, where each event is one photon. Note this
peaks at u = hν/kT ≈ 1.5936 (≈133 GHz at 4 K), *not* at the familiar
energy-spectrum peak hν = 2.82 kT (≈235 GHz at 4 K).

Direction is sampled uniformly over the outward hemisphere (θ uniform in
[0°, 90°], φ uniform in [0°, 360°]). This is Chang's original convention and is
intentional — a known approximation whose effect washes out after multiple
reflections inside a cavity.

---

## Gun configuration

The particle gun is a 500 GHz optical photon. Gun commands are read fresh each
event and are valid both before and after `/run/initialize`.

```mac
/run/initialize

/bbr/gun/mode true         # true  = fixed particle gun
                           # false = Planck thermal emitter (default)

/bbr/gun/posX -20.0        # gun X position [mm]
/bbr/gun/posY   0.0
/bbr/gun/posZ   0.0        # z=0 → crack1, z=3 → crack2, z>~5 → solid Cu

/bbr/gun/dirX 1.0          # momentum direction (normalized internally)
/bbr/gun/dirY 0.0
/bbr/gun/dirZ 0.0

/bbr/gun/energy_eV 2.07e-3 # photon energy in eV (500 GHz = 2.07e-3 eV)
/bbr/gun/pol 0 0 0         # polarization: 0 0 0 = random (default); any other
                           # vector is projected perpendicular to the direction

/run/beamOn 1000
```

---

## Output files

| File | Written by | Contents |
|---|---|---|
| `output/bbr.root` | `BBRRunAction` / `BBRTestSteppingAction` (via `G4AnalysisManager`) | Two ntuples. **`crossings`** — one row per optical-photon boundary crossing: run_id, event_id, position, energy, pre/post momentum, incidence angles, volume/material/status/event-type codes, per-track crossing count, and `hfss_freq_GHz` (the HFSS grid frequency selected for a crack entry, `-1` on every other row). **`abspoints`** — one row per photon termination (`fStopAndKill`): position, energy, final momentum, n_reflect, terminating volume + status codes. `term_status` is `WorldExit` for a photon that left the world, `BulkAbsorption` for a `G4OpAbsorption` kill, otherwise the boundary status of the killing step (e.g. `BBRAbsorb`). Categorical fields are integer codes; runs are multithreaded and the per-thread ntuples are merged into this one file. |
| `output/bbr_legend.json` | `BBRRunAction` (master thread) | `{category: {code: name}}` dictionary decoding the integer code columns (status / event_type / volume / material). Consumed by `analysis/bbrsim/io.py`. |
| `bbrsim stdout` | redirect from stdout | `[BBR] reflectance` and `[BBR] diffraction` running tallies |
| `build/cu_reflectance_plots.png` | `plot_cu_reflectance.py` | 3-panel reflectance / absorptance / temperature-dependence plot |

Read the ROOT output in Python via the shared loader:

```python
import sys; sys.path.insert(0, "analysis")
from bbrsim.io import load
crossings, abspoints = load("build/output/bbr.root")   # decoded DataFrames
```

---

## Analysis scripts

All scripts live in `scripts/` and must be run from the repo root with
`conda run -n bbrsim python scripts/<name>.py`. Every one of them reads BBRsim
output as `build/output/bbr.root` through the `analysis/bbrsim` loader — none
read a BBRsim-produced CSV. (`plot_crack_angular.py` and `plot_cu_reflectance.py`
additionally load external reference CSVs: the HFSS far-field tables and the
Palik / Serov / Geant4-IR copper comparison sets.)

Physics formulas are not duplicated per script: `analysis/bbrsim/physics.py` is
the single source of truth for the complex-Drude reflectance, the Planck
photon-number spectrum, and Hagen-Rubens, mirroring the C++ implementation and
self-tested by `check_physics.py`.

| Script | Purpose | Key flags |
|---|---|---|
| `check_physics.py` | Self-test of `bbrsim.physics` against C++ anchors | — |
| `check_reflectance.py` | Poisson test of absorbed count vs full-Drude theory (reflectance.mac) | `--root`, `--RRR`, `--T_K`, `--freq` |
| `check_cu_absorptance.py` | Compare A_obs from decoded crossings to Planck-weighted Drude theory; PASS/FAIL | positional path, `--rrr`, `--temp` |
| `check_cu_serov.py` | Full-Drude D for the OF_Cu / HP_Cu RRR aliases vs Serov (2016) measured points (±10%). HP_Cu = RRR 6 is currently 13% low — an open decision | — |
| `check_crack_wall_reflection.py` | crack_wall.mac: crack→Cu wall reflections flip pz, keep px/py | positional path |
| `check_no_photons_in_metal.py` | Invariant: no crossing starts inside a Cu / perfect-metal material | positional path |
| `check_term_status.py` | Invariant: no `unknown`; world exits are `WorldExit`; absorptions have a volume; BBRAbsorb counts match | positional path |
| `check_crack_transmittance.py` | crack_transmit.mac: T_obs = 0.50 ± 3σ, no tangential exits | positional path |
| `check_crack_oblique.py` | crack_oblique.mac: 16 runs at 45° off normal — k_y-sign conservation, T vs table, exit marginals vs the `bbrsim.hfss` mirror, ±y/±z mirror KS, specular reflection (138 checks) | positional output dir (default `build/output`) |
| `check_crack_frequency.py` | crack_frequency.mac against the mock five-frequency tree: per-photon selection vs the `bbrsim.hfss` mirror, T_obs vs the selected dataset, far-field signature, clamping, Planck bin coverage (89 checks) | positional output dir, `--data-dir MOCK_ROOT` |
| `notebooks/copper_reflectance.ipynb` | Narrative walkthrough of the copper Drude model: bounce counts, the corner frequency, Hagen-Rubens breakdown, the relaxation plateau, the 24-point Geant4 table, and validation against Serov 2016. Needs `pip install ipykernel` in the `bbrsim` env | open in VS Code / Jupyter |
| `plot_crack_frequency.py` | Four-panel overview of the frequency-keyed lookup: selection staircase, observed vs predicted transmittance, the Planck spectrum coloured by dataset, and the far-field exit signature | positional output dir, `--data-dir MOCK_ROOT`, `--out PATH` |
| `make_mock_hfss_frequencies.py` | Build the mock multi-frequency HFSS tree from the real 500 GHz data (scaled transmittance + truncated far field per frequency); refuses to write inside `data/` | `--src --dst --ids [--freqs --scales]` |
| `run_regression.sh` | Build + eight regression macros (parallel, under `<build>/regression/`) + every validator, one PASS/FAIL line each; exit code = unexpected failures, known reds reported as XFAIL | positional build dir (default `build`) |
| `check_planck_spectrum.py` | Validate emitted spectrum against Planck photon-number peak | positional path, `--temp <K>` |
| `check_nreflect.py` | Per-track reflection-count distribution sanity checks | — |
| `check_angle_distribution.py` | KS test of Cu incidence angles | — |
| `check_crack_ratio.py` | crack2/crack1 rate ratio vs aperture ratio | positional path |
| `plot_cu_reflectance.py` | Reflectance/absorptance curves vs frequency, temperature panel | `--out <path>` |
| `plot_crack_angular.py` | Outgoing crack angular distributions vs HFSS far-field theory | `--iwt`, `--iwp` |
| `plot_test_output.py` | Overview plots of the boundary-crossing output | `--temp <K>` |

`check_angle_distribution.py` is a KS test and is therefore N-sensitive: it
passes on its designated 10k-event `planck.mac` workload but over-rejects on a
5M-event run. That is a property of the test, not a regression.

---

## Troubleshooting

**`[BBR] setCuMaterial: unknown alias`**
Only `OFHC_Cu`, `OF_Cu`, and `HP_Cu` are valid aliases. Use `/bbr/det/setCuRRR <N>`
for any other grade.

**`[BBR] setCuRRR: RRR must be >= 1`**
RRR must be a positive integer. RRR < 1 has no physical meaning.

**No `[BBR] reflectance` lines in output**
The reflectance tally prints every 1000 Cu hits. Check that the gun position
and direction point at the Cu face (x = 0 plane) away from the cracks
(`/bbr/gun/mode true`, `/bbr/gun/posZ 10.0`), or that enough Planck-mode
events have run.

**`RESULT: FAIL` from `check_cu_absorptance.py`**
A_obs/A_theory outside [0.3, 3.0]. Common causes: too few events (< 1000 give
high statistical noise — at RRR=100 the absorptance is ~4.9×10⁻⁵, so a 10k run
yields ~0–2 absorptions), or an `--rrr`/`--temp` override that contradicts the
Cu material actually used in the run. Without overrides the script parses RRR
and T from the material name, so a mismatch is usually a stale override.

**`check_cu_absorptance.py` reports "no Cu boundary crossings found"**
You passed the wrong file, or the run never reached the Cu slab. The script
takes a **ROOT path** (positional, default `build/output/bbr.root`). The old
`check_cu_absorptance.py out.txt --rrr 100` stdout-parsing form no longer
exists.

**Python scripts crash with `ModuleNotFoundError`**
Run as `conda run -n bbrsim python scripts/<name>.py`, not `python3 scripts/<name>.py`.
