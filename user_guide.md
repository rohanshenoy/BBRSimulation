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
mkdir build && cd build
cmake ..
make
```

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
observed transmittance. Expected at 500 GHz normal incidence:

| Crack | Gap | Expected T |
|---|---|---|
| crack1 | 52 µm | ≈ 52.7% |
| crack2 | 102 µm | ≈ 50.4% |

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

> **Note (Phase A):** `check_crack_ratio.py` and `plot_crack_angular.py` still
> read the legacy CSV via pandas and are being ported to the ROOT loader
> (`analysis/bbrsim/io.py`). Only `check_planck_spectrum.py` and
> `check_reflectance.py` read `output/bbr.root` today.

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
Planck-weighted Drude integral:

```bash
./BBRSim test.mac >out.txt 2>&1
conda run -n bbrsim python scripts/check_cu_absorptance.py out.txt --rrr 100
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

## Gun configuration

The particle gun is a 500 GHz optical photon. All gun commands must appear
after `/run/initialize`.

```mac
/run/initialize

/bbr/gun/mode true         # true = direct-hit mode (solid Cu, z=5 mm)
                           # false = crack-test mode (default)

/bbr/gun/posX -20.0        # gun X position [mm]
/bbr/gun/posY   0.0
/bbr/gun/posZ   0.0        # z=0 → crack1, z=3 → crack2, z>~5 → solid Cu

/bbr/gun/dirX 1.0          # momentum direction (normalized internally)
/bbr/gun/dirY 0.0
/bbr/gun/dirZ 0.0

/bbr/gun/energy_eV 2.07e-3 # photon energy in eV (500 GHz = 2.07e-3 eV)

/run/beamOn 1000
```

---

## Output files

| File | Written by | Contents |
|---|---|---|
| `output/bbr.root` | `BBRRunAction` / `BBRTestSteppingAction` (via `G4AnalysisManager`) | Two ntuples. **`crossings`** — one row per optical-photon boundary crossing: run_id, event_id, position, energy, pre/post momentum, incidence angles, volume/material/status/event-type codes, per-track crossing count. **`abspoints`** — one row per photon termination (`fStopAndKill`): position, energy, final momentum, n_reflect, terminating volume + status codes. Categorical fields are integer codes; runs are multithreaded and the per-thread ntuples are merged into this one file. |
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
`conda run -n bbrsim python scripts/<name>.py`. `check_planck_spectrum.py` and
`check_reflectance.py` read `build/output/bbr.root` via the `analysis/bbrsim`
loader; the remaining `check_*`/`plot_*` scripts still default to the legacy
CSV and are being ported to the ROOT loader. All accept a path argument.

| Script | Purpose | Key flags |
|---|---|---|
| `check_reflectance.py` | Poisson test of absorbed count vs full-Drude theory (reflectance.mac) | `--RRR`, `--T_K`, `--freq` |
| `check_cu_absorptance.py` | Compare stdout A_obs to Planck-weighted Drude theory; PASS/FAIL | `--rrr <N>` |
| `check_cu_serov.py` | Verify σ_eff values reproduce Serov (2016) reference points | — |
| `check_planck_spectrum.py` | Validate emitted spectrum against Planck photon-number peak | `--temp <K>` |
| `check_nreflect.py` | Per-track reflection-count distribution sanity checks | — |
| `check_angle_distribution.py` | KS test of Cu incidence angles | — |
| `check_crack_ratio.py` | crack2/crack1 rate ratio vs aperture ratio | — |
| `plot_cu_reflectance.py` | Reflectance/absorptance curves vs frequency, temperature panel | `--out <path>` |
| `plot_crack_angular.py` | Outgoing crack angular distributions vs HFSS far-field theory | `--iwt`, `--iwp` |
| `plot_test_output.py` | Overview plots of the boundary-crossing output (legacy CSV; ROOT port pending) | `--temp <K>` |

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
high statistical noise), wrong `--rrr` flag (must match the mac file's RRR),
or a photon energy far from 500 GHz (Planck-weighted theory assumes a 4 K
spectrum).

**Python scripts crash with `ModuleNotFoundError`**
Run as `conda run -n bbrsim python scripts/<name>.py`, not `python3 scripts/<name>.py`.
