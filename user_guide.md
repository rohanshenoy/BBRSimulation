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

The geometry: a 50 cm vacuum world and two touching Cu plates with their front
faces at x = 0, each as thick as its crack, holding two `vacuum_wg` crack
volumes that equal the HFSS models (crack1: 1.0 mm long, 10 mm wide, 50 µm gap
at z = 0, in a 1.0 mm plate; crack2: 1.5 mm long, 10 mm wide, 100 µm gap at
z = 3 mm, in a 1.5 mm plate).

In the commands below, `./bbrsimTestWorld` runs in `examples/testworld/build/`
and the `conda run` commands run from the repository root. The light pipe
application moved out of this repository on 2026-10-26; BBRSimulation keeps the
test world.

## Test cases

### 1. Planck thermal emitter

The default mode. A 1×20×20 mm box at x = −50 mm emits photons with energies
drawn from the Planck photon-number spectrum (10 GHz–20 THz) toward the plates
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
| crack1 | 50 µm | 50.0 ± 0.25 % |
| crack2 | 100 µm | 50.3 ± 0.25 % |

The expected value is exactly 50 %. A parallel-plate gap narrower than λ/2 is
a perfect polarization filter: only the TEM mode, which has no cutoff,
transmits, so unpolarized light transmits half. The raw HFSS power ratios come
out slightly above 1 (1.0545 for crack1, a port-normalization artefact). When
the tables load, an incidence key whose largest transmittance over linear
polarizations exceeds 1 has both polarizations' transmittances divided by it,
so the maximum is 1 and their ratio is kept; a `[BBR] HFSS ... normalized to 1`
line reports each such key.

The crack volumes equal the HFSS models, so HFSS exit points that lie on a wall
sit on the Geant4 wall too. When a frequency's tables load, every exit point
closer than 1 nm to the declared cross-section's boundary, or past it by at most
the 1e-6 relative tolerance, is moved to 1 nm inside it, with its weight and
its Ephi pairing unchanged, and one line reports the count:
`[BBR] HFSS <id>_<freq>GHz: N exit point(s) clamped 1 nm inside the declared
cross-section` (4500 for crack1 and 6000 for crack2 at 500 GHz, 15 keys each).
After an HFSS transmission the photon is placed 2e-8 mm (20 times Geant4's
surface tolerance) inside the crack's exit face, so Geant4 11.1's boundary
process does not warn `OpBoun06` on its first step into the world.

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
`<id>_<freq>GHz_Ephi=N` directories under `waveguides/`, each pair with its
`<id>_<freq>GHz.dataset.json` sidecar (see `validation/README.md`, Dataset
sidecars). There is no interpolation between frequencies; a photon beyond the
grid uses the nearest edge, with one `BBR008` warning per crack and side.
When a photon's frequency and the dataset serving it lie on opposite sides of
the guide's lowest-mode cutoff (`modes.cutoff_ghz` in the sidecar: TE10 at
14.99 GHz for the slot cracks, TE11 at 1757 GHz for the 50 µm round gap),
BBRsim warns once with `BBR026`. The 2000 GHz round-gap data are therefore for
fixed-energy gun runs only until a dataset below cutoff exists. A schema-2.x
sidecar lists every regime cutoff (`modes.regime_cutoffs_ghz`) and the warning
names the mode it crosses (`mode=<label>`); a 2.x table at or above its
cross-term limit also raises `BBR028` once.

After a transmission the photon starts just inside the crack on the ray through
its sampled exit point. `exit_start_normal` in the `[BBR] run <id> ...` line counts
the transmissions that start beside a wall and graze away from it (the on-ray start
would cross the wall behind); they never hit a wall, and their exit-plane crossing
moves inward by at most inset x t/n. Photons grazing toward a wall use the on-ray
start and are counted by neither counter, and the tree check's `exit direction` row
proves none of them hits a wall; `exit_start_fallback` counts starts valid in
neither way (see `validation/README.md`, Dataset sidecars).

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
| `/bbr/testworld/roundGap` | before `/run/initialize` | no | test world only, default `false`: `true` also places the opt-in straight round gap `RoundGap_r50um`, which needs a data root holding a `RoundGap_r50um` dataset (`/bbr/dataDir`) |
| `/bbr/testworld/cracks` | before `/run/initialize` | no | test world only, default `true`: `false` leaves out the two crack plates and their cracks, for a crack-free world; the round gap is unaffected |
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
the test world; other geometries need their own.

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

**Selecting crack entries.** Every crack entry is logged from the world side:
`mat_pre` and `vol_pre` are `G4_Galactic` and `World`, and the entered
material and volume are `mat_post` and `vol_post`. So a crack entry is
`mat_post == "vacuum_wg"`, split by `vol_post`, and a first copper hit is a
`mat_post` starting with `Cu_RRR`. The rows with `mat_pre == "vacuum_wg"` are
photons leaving a crack: into the world (after an HFSS transmission or, for a
photon that starts inside the crack as in the crack-wall fixture, after
bouncing between the walls), or into the copper when the photon starts inside
the crack.

## Scripts

Run every script as `conda run -n bbrsim python <script>`. Each reads BBRsim
output (default `output/bbr.root` in the current directory) through the
`bbrsim` package.

**Validators** (`validation/`), run by the regression runner:

| Script | Checks |
|---|---|
| `check_reflectance.py` | absorbed count vs the Drude model (`--root`, `--RRR`, `--T_K`, `--freq`) |
| `check_cu_serov.py` | Drude loss for `OF_Cu` / `HP_Cu` vs Serov (2016) within ±10 %; `HP_Cu` is the known expected failure |
| `check_planck_spectrum.py` | the emitted spectrum at `--temp`: its peak, and a KS test against the Planck photon-number CDF truncated to the emitter band |
| `check_nreflect.py` | the per-photon reflection-count distribution |
| `check_angle_distribution.py` | incidence angles at the copper vs uniform-in-θ emission (KS test, sized for the 10 000-event Planck run) |
| `check_crack_transmittance.py` | T = 0.50 at normal incidence, no exits along the crack face |
| `check_crack_ratio.py` | crack2/crack1 entry ratio vs the aperture ratio, within 3 σ |
| `check_crack_wall_reflection.py` | reflection off a crack's side wall flips only p_z |
| `check_crack_oblique.py` | 45° incidence: 138 checks against the HFSS tables and the Python model |
| `check_crack_frequency.py` | dataset choice per photon on the mock five-frequency tree: 89 checks (`--data-dir`) |
| `check_round_gap.py` | the opt-in straight round gap at 2000 GHz on its mock dataset: transmittance per polarization, exit positions inside the HFSS radius, the radial exit profile per polarization and the mean exit direction, 49 checks over nine runs, one of them the startup crack line in the run log (`--data-dir`, `--log`) |
| `check_invariants.py` | no photon ever travels inside a metal, and every photon death is labelled correctly (`--allow-no-crossings` for a run that crosses no boundary) |
| `check_dataset_sidecars.py` | reads HFSS data trees, not output (the runner passes `data/waveguides` and the three mock trees): every dataset is a complete trio with a valid schema-1 or 2.x sidecar (the mode lists or regime cutoffs, CSV checksums, C1-C5), and the frequencies of one ID agree on the frequency-independent physics; see [Dataset sidecars](validation/README.md#dataset-sidecars) |
| `check_tree_run.py` | reads a run that `validation/Scripts/make_tree_check.py` generated for an HFSS data tree: every photon used the right frequency table, and transmission, reflection and exit directions follow it, against `bbrsim.hfss` reading the same tree (`--log`); see [Real-tree check](validation/README.md#real-tree-check) |

**Validator run by hand**, because it needs a large run:
`check_cu_absorptance.py` (`planck_5M.mac`), shown above.

**Checking a new data tree**: `validation/Scripts/make_tree_check.py <root> --out <dir>`
writes a test-world macro for every dataset under `<root>/waveguides` that the
test world can place, and `check_tree_run.py` checks the run; the recipe is in
[Real-tree check](validation/README.md#real-tree-check).

**Python tests** (`tools/python/tests/`, pytest): the `bbrsim` formulas against
analytic limits and reference values, the HFSS mirror, the ROOT loader, and the
validators on synthetic bad outputs. Run
`conda run -n bbrsim python -m pytest -q -p no:cacheprovider tools/python/tests`
from the repository root (pytest comes with `pip install -e "tools/python[test]"`).

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

**`/run/initialize` aborts with `BBR024`** (`HFSS dataset sidecar ... cannot
be opened`, or a sidecar field that is missing, mistyped or names another
dataset or frequency): every `<id>_<freq>GHz_Ephi={0,1}` pair needs its
`<id>_<freq>GHz.dataset.json` beside it, and an old data tree or a hand-copied
dataset has none. Fetch the sidecar the HFSS runner wrote with the dataset, or
write one (`bbrsim.sidecar.build_from_csvs`; the mock generators show how),
then run `conda run -n bbrsim python validation/check_dataset_sidecars.py
<dir>/waveguides` before the Geant4 run. The message names the absolute
sidecar path and the field. See [Dataset
sidecars](validation/README.md#dataset-sidecars). A `BBR024` that names
`geometry.extent_mm.l` or `.g` and the full width of the declared section means
the sidecar contradicts itself (for a rectangle the two must agree to 1e-6
relative): regenerate the sidecar from the producer's job, do not edit one number.
Under schema 1, `provenance.material` must be a string or null; any other type is `BBR024`.
Under schema 2.x it is required: null, an absent key or no `provenance` is
`BBR024`, and the type message then reads "must be a string under schema 2.x".
Other 2.x causes: a `schema_version` that is neither 1.x nor 2.x ("BBRsim reads
1.x and 2.x"); two frequencies of one ID that differ in schema major or in the
incidence grid (the message names the block); a malformed or contradictory
per-key record; a `dark` basis whose `t_flux` is not 0; a rectangle axis with
fewer than 2 nodes; an `exit_field.grid` step that is not (max - min)/(count - 1);
a negative per-key `t_flux`; or, in the Python mirror only, a 16-column 2.x
`waveguide.csv` read without its sidecar. In the C++ library a 16-column table
loaded without a 2.x sidecar is `BBR013` (wrong field count; through the crack
library a dataset without any sidecar stops at discovery as `BBR024`).

A 2.x `BBR012` at load that says "has no per-key record" means a CSV key is off
its declared incidence angle by more than 1e-9 degree (inside the 0.01 degree
rounding of the key checks): fix the angle in the export, never the sidecar. One
that says "OutgoingPower ... is negative" means a row holds a negative exit flux;
the producer stops on that, so the file was altered or came from another tool.

**`BBR025`**: the dataset declares a convention or a geometry BBRsim does not
implement. The message names the field (a frame, incidence, far-field or
polarization definition, the boundary model, or the modes block) and the one
convention the sampler implements; or, at `/run/initialize`, it names the crack volume whose
solid does not equal the declared exit cross-section (the Geant4 section must
equal the HFSS one to within 1 nm, or 1 nm plus 1e-6 relative where larger, as
the test-world cracks do since OR-020; check also that local x is the
propagation axis: a bare `G4Tubs` has it on z); or, at the first CSV load of a
frequency, it counts the exit points outside the declared cross-section. Two further messages name the
length and the material. `the crack volume's length ... differs from the HFSS model's ... mm` means the
solid's extent along local x is not `geometry.extent_mm.p` (within 1 nm, or 1 nm plus 1e-6 relative where
longer): the Geant4 crack must be exactly as long as the HFSS model, as the test-world plates make it.
`provenance.material is "..."; the gap must be vacuum` means the dataset was solved with a gap material
other than vacuum, which the sampler does not model (a null material in a legacy sidecar is accepted). A dataset produced outside the
canonical pose is rejected, not remapped. See [Dataset
sidecars](validation/README.md#dataset-sidecars). Under schema 2.x it also names a fixed string, a recomputed value
or a mode the reader's own derivation does not give: `transmittance.cell_90`,
a `c_k`, `modes.regime_cutoffs_ghz`, `cross_term.valid_below_ghz`,
`cell_90_assumes` (an emitter other than `uniform_theta`), `exit_field.columns`
or a `reference.passed` that is not true; a declared theta list without a 90 degree
entry is refused the same way.

**`BBR012`**: a CSV does not match its sidecar (a far-field grid, exit lattice or row count
different from the declaration). Under schema 2.x it also reports `InterceptFactor`
that is not constant within a key or differs from the key's c_k or from the closed
form of the declared grid; a dark basis with a nonzero `OutgoingPower` on any row; and
a key whose T, `OutgoingPower` / (`IngoingPower` x `InterceptFactor`), differs from the
sidecar's `t_flux` / c_k (1e-6 relative). The message names the key, the polarization
and the value. Regenerate the dataset; BBRsim never resets T to make a row fit.

**`BBR027`**: a dataset gives a transmitting polarization nowhere to go. At
the first load of a frequency, a key's basis polarization has T > 0 but its
far field has no outward intensity (every row zero or in the exit-face plane)
or its exit field is zero at every exit point; mid-run, the two basis far
fields or exit fields cancel for one photon's mixed polarization. BBRsim no
longer draws a uniform direction or exit point instead. The export is
inconsistent: check it with `check_dataset_sidecars.py` and ask the HFSS side
for a re-export. See [Dataset sidecars](validation/README.md#dataset-sidecars).

**A `BBR026` banner in a Planck run**: a warning, not an error. A photon's
frequency and the grid frequency serving it lie on opposite sides of a regime
cutoff of the guide (the lowest mode, `modes.cutoff_ghz`, under schema 1; every
`modes.regime_cutoffs_ghz` entry under 2.x), so the table was computed in the
other propagation regime; the photon is served unchanged and the warning
appears once per dataset, mode and direction (`mode=<label>` in the text). Expected in a broadband run on a
single-frequency grid (the 2000 GHz round gap is above its 1757 GHz cutoff,
which is why it is for gun runs only); unexpected in a fixed-energy gun run
whose energy was meant to match the dataset, where it signals the wrong
dataset or frequency. The runner tolerates it only where it is expected
([Dataset sidecars](validation/README.md#dataset-sidecars)).

**A `BBR028` banner**: a warning, not an error, for schema-2.x data only. A table at
or above the dataset's polarization-filter limit (`cross_term.valid_below_ghz`)
served a photon, so the polarization cross term, declared valid only below that
limit, is not trusted for it; the table is used unchanged and the warning appears
once per dataset. Expected for the 2000 GHz round gap (limit 1757 GHz); in a run
meant to stay below the limit it means the wrong frequency was served. The
real-tree check expects it exactly where the manifest implies it
([Real-tree check](validation/README.md#real-tree-check)).

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
