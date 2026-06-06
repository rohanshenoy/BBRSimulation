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

### 1. Wrapper regression

Verifies that the `BBSimOpBoundaryProcess` wrapper is transparent — identical
output to the stock Geant4 `G4OpBoundaryProcess` on an ordinary geometry.

```bash
./BBRSim verify_wrapper.mac
```

Run this before and after any change to `BBSimOpBoundaryProcess` or `BBSimPhysics`
and confirm the output is bit-identical. This is the core regression guard.

---

### 2. HFSS diffraction

Fires a 500 GHz optical photon at a narrow parallel-plate crack and measures
the fraction that transmits. The HFSS S-parameter tables in `data/` encode the
diffraction physics.

```bash
./BBRSim diffraction.mac          # crack1: b = 52 µm at z = 0
./BBRSim diffraction_crack2.mac   # crack2: b = 102 µm at z = 3 mm
./BBRSim validation.mac           # both cracks in one session
```

Expected transmittances at 500 GHz normal incidence:

| Crack | Gap | Expected T |
|---|---|---|
| crack1 | 52 µm | ≈ 52.7% |
| crack2 | 102 µm | ≈ 50.4% |

Output is written to `diffraction_output.csv`:

```
crack_id,x,y
InfParallelPlate_crack1Rohan_500GHz,3.14,-1.02
...
```

Plot the exit-position distribution and per-crack transmittance:

```bash
conda run -n bbrsim python ../scripts/plot_diffraction.py
```

---

### 3. Copper reflectance

Fires photons at a solid Cu wall and measures the fraction absorbed.
The `[BBR] reflectance` output line is printed every 1000 events.

```bash
./BBRSim test.mac          # OFHC_Cu (RRR=100), 500 GHz, 10 000 events
./BBRSim test_of_cu.mac    # OF_Cu   (RRR=3),   500 GHz, 2 000 events
./BBRSim test_hp_cu.mac    # HP_Cu   (RRR=6),   500 GHz, 2 000 events
```

The output line format is:

```
[BBR] reflectance mat=Cu_RRR100_T4K N=10000 A_obs=0.001600 R_theory=0.998377
```

- `A_obs` = fraction of photons absorbed (1 − R_obs)
- `R_theory` = Drude model reflectance at the gun frequency

Compare `A_obs` against the Planck-weighted Drude theory:

```bash
./BBRSim test.mac >out.txt 2>&1
conda run -n bbrsim python ../scripts/check_cu_absorptance.py out.txt --rrr 100
```

Expected output:

```
[BBR] log  mat=Cu_RRR100_T4K  N=10000  A_obs=1.600000e-03
RRR                        = 100
sigma_DC                   = 5.9600e+09 S/m
tau                        = 2.497 ps
A_theory (Planck+Drude)    = 1.XXXXXXe-03
ratio A_obs/A_theory       = X.XXX  (expected [0.3, 3.0])
RESULT: PASS
```

Plot reflectance curves for all three grades across the full 50 GHz–20 THz range:

```bash
conda run -n bbrsim python ../scripts/plot_cu_reflectance.py
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
contribution. The full Drude model (Griffiths §9.4) then gives R(ω) across
the 50 GHz–20 THz range. See [docs/physics/copper_reflectance_model.md](physics/copper_reflectance_model.md)
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
/bbr/gun/posZ   0.0        # or use the shortcut below for crack-test
/bbr/gun/setZ   3.0 mm     # shortcut: set gun Z (crack2 is at z=3 mm)

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
| `diffraction_output.csv` | `BBRDiffractionSteppingAction` | Crack ID + exit position (x, y) per transmitted photon |
| `bbrsim_stdout.txt` | redirect from stdout | `[BBR] reflectance` lines for `check_cu_absorptance.py` |
| `build/cu_reflectance_plots.png` | `plot_cu_reflectance.py` | 3-panel reflectance / absorptance / temperature-dependence plot |

---

## Analysis scripts

All scripts live in `scripts/` and must be run from the repo root with
`conda run -n bbrsim python scripts/<name>.py`.

| Script | Purpose | Key flags |
|---|---|---|
| `plot_diffraction.py` | 3-panel diffraction plots (transmittance vs event count, exit position) | `--out <path>` |
| `plot_cu_reflectance.py` | Reflectance/absorptance curves vs frequency, temperature panel | `--out <path>` |
| `check_cu_absorptance.py` | Compare simulated A_obs to Planck+Drude theory; PASS/FAIL | `--rrr <N>` |
| `check_cu_serov.py` | Verify σ_eff values reproduce Serov (2016) reference points | — |
| `check_reflectance.py` | Standalone Drude reflectance at a single frequency | — |
| `plot_test_output.py` | Quick plot of test.mac CSV output | — |

---

## Troubleshooting

**`[BBR] setCuMaterial: unknown alias`**
Only `OFHC_Cu`, `OF_Cu`, and `HP_Cu` are valid aliases. Use `/bbr/det/setCuRRR <N>`
for any other grade.

**`[BBR] setCuRRR: RRR must be >= 1`**
RRR must be a positive integer. RRR < 1 has no physical meaning.

**No `[BBR] reflectance` lines in output**
The reflectance logger only fires when a photon hits the Cu slab face-on in
`test.mac` gun mode. Check that `/bbr/gun/mode true` is set and that the gun
position and direction point at the Cu face (x = 0 plane).

**`RESULT: FAIL` from `check_cu_absorptance.py`**
A_obs/A_theory outside [0.3, 3.0]. Common causes: too few events (< 1000 give
high statistical noise), wrong `--rrr` flag (must match the mac file's RRR),
or a photon energy far from 500 GHz (Planck-weighted theory assumes a 4 K
spectrum).

**Python scripts crash with `ModuleNotFoundError`**
Run as `conda run -n bbrsim python scripts/<name>.py`, not `python3 scripts/<name>.py`.
