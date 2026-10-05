# BBRsim validation

Fixed-seed fixtures and PASS/FAIL validators that gate every change to BBRsim.

```
G4Macros/Validation_*.mac   the six validation-only fixtures
check_*.py                  the 12 validators
Scripts/run_regression.sh   build + fixtures + validators, one command
Scripts/drift_guards.sh     source-tree consistency checks (run by the runner, or alone)
Scripts/make_mock_hfss_frequencies.py   the mock HFSS tree for Validation_CrackFrequency
Scripts/consumer_smoke/     external find_package(BBRsim) + link smoke test
Scripts/numbers.baseline    the three fixed-seed numbers that BBR_PIN=1 compares
Scripts/tests/              make_bad_output.py (synthetic outputs for the validator negative tests),
                            test_env.sh (the env scripts in bash, zsh, tcsh, dash), check_links.sh
```

## Running

```bash
validation/Scripts/run_regression.sh [BUILD_DIR]     # default: build
```

The runner builds and installs the library with the C++ tests
(`BUILD_BBRSIM_TESTS=ON`; prefix `BBR_PREFIX`, default `install/`), sources
the installed env script, builds both examples, runs the drift guards, the
consumer smoke test and the `bbrsim data default` and `version stamp` checks,
then step 1c: the C++
tests (`ctest`), the Python tests (`pytest`, from `/` on the installed
`bbrsim`; pytest must be installed, `pip install -e "tools/python[test]"`) and
`Scripts/tests/test_env.sh` (`env scripts`). It runs the ten cases in parallel,
each in `BUILD_DIR/regression/<case>/`, with `G4FORCENUMBEROFTHREADS` pinned to
8 (`BBR_THREADS` overrides; Geant4 warns when a run has more threads than the
square root of its events, and the log scan would fail on a many-core machine),
and runs the validators on their output. Six cases run the fixtures here; four
run example macros directly (see [Regression inputs](#regression-inputs)).
Last come two rows: `installed examples` (both examples installed into
`BUILD_DIR/expfx`, never the prefix, run without `BBRSIMDATA` or
`DYLD_LIBRARY_PATH`: `Validation_CrackTransmit.mac` must pass
`check_crack_transmittance.py`, `lightpipe.mac` must write its output) and
`tools smoke` (the four `tools/plot_*.py` and the notebook's code cells on the
fixture output, in `BUILD_DIR/tools_smoke`). With `BBR_PIN=1` the runner also
compares the three fixed-seed numbers (reflectance pull, crack T_obs, Planck
peak ratio) exactly with `Scripts/numbers.baseline` (row `fixed-seed
numbers`), for refactors that must not change behaviour. The exit code is the
number of unexpected failures; a green run ends with the line
`fail=0  xfail=1  xpass=0`; the pass count includes the optional fixed-seed row
when `BBR_PIN=1`. New rows cover batch errors and the optional in-tree example
build. The application provenance row edits and rebuilds a standalone example
copy, checking that its recorded fingerprint changes while the installed library
fingerprint stays fixed; pytest also exercises real redirected output and per-result metadata.

| Case | Macro | Executable | Validators |
|---|---|---|---|
| `refl` | `examples/testworld/G4Macros/reflectance.mac` | `bbrsimTestWorld` | `check_reflectance.py --root`, `check_invariants.py` |
| `planck` | `examples/testworld/G4Macros/planck.mac` | `bbrsimTestWorld` | `check_planck_spectrum.py --temp 4`, `check_nreflect.py`, `check_angle_distribution.py`, `check_invariants.py` |
| `config_mt` | `examples/testworld/G4Macros/config_mt.mac` | `bbrsimTestWorld` | `check_planck_spectrum.py --temp 4` on run 0 (`bbr_mt_r0.root`) and `--temp 10` on run 1 (`bbr_mt_r1.root`); `check_invariants.py` on each |
| `wall` | `Validation_CrackWall.mac` | `bbrsimTestWorld` | `check_crack_wall_reflection.py`, `check_invariants.py` |
| `exit` | `Validation_WorldExit.mac` | `bbrsimTestWorld` | `check_invariants.py --allow-no-crossings` (the fixture crosses no boundary) |
| `transmit` | `Validation_CrackTransmit.mac` | `bbrsimTestWorld` | `check_crack_transmittance.py` on run 1 (`bbr.root`); `check_crack_ratio.py` on run 2 (`bbr_ratio.root`); `check_invariants.py` on both (one row) |
| `oblique` | `Validation_CrackOblique.mac` | `bbrsimTestWorld` | `check_crack_oblique.py` on `output/`; `check_invariants.py` on each of the 16 `bbr_oblique_rNN.root` |
| `frequency` | `Validation_CrackFrequency.mac` | `bbrsimTestWorld` | `check_crack_frequency.py --data-dir mock_hfss` on `output/`; `check_invariants.py` on each of the 16 `bbr_freq_rNN.root` |
| `lp` | `examples/lightpipe/G4Macros/lightpipe.mac` | `bbrsimLightPipe` | `check_invariants.py` |
| `lp_cad` | `Validation_LightPipeCAD.mac` | `bbrsimLightPipe` | `check_invariants.py` (cad mode, the bundled `box_sample.stl` through `BBRSIMDATA`) |
| — | (no ROOT input) | — | `check_cu_serov.py` (XFAIL) |

A case passes when the binary exits 0, writes its ROOT output, and its log
holds no `GeomNav`, `G4Exception`, `BBR0xx` or `LP002` line (the frequency case
tolerates its `BBR008` clamp warnings, and exactly one per side is required).

## Drift guards

`Scripts/drift_guards.sh` prints one PASS/FAIL row per guard and exits with the
FAIL count; the runner adds its rows to the summary.

- `library file lists`, `example file lists`, `tools file lists` — every library source and header, example source and macro, `bbrsim` module and plot script is named in its `CMakeLists.txt`.
- `library names no example` — no library file names an example class.
- `every fixture is run` — the runner calls every `G4Macros/Validation_*.mac`.
- `regression macros pinned` — see [Regression inputs](#regression-inputs).
- `notebooks without outputs` — tracked notebooks carry no outputs or execution counts.
- `markdown links` — every relative link and in-page `#anchor` in the tracked `*.md` resolves (`Scripts/tests/check_links.sh`).
- `banned names` — no tracked file (outside `ChangeHistory`, `data/` and the guard itself) edits the Python module search path, reads the HFSS data by a relative path instead of through `bbrsim.paths`, or names one of the executables, script paths, property keys or output files the reorganization retired. The list is in the guard; quoting it anywhere else trips it.
- `action copies identical` — the four `TestWorld*` / `LightPipe*` action classes match once the prefix is normalised, so a fix to one copy reaches the other.

In an installed copy of `validation/` (not the top of a git checkout, no
`examples/`) the guards that need them report themselves skipped.

`Scripts/tests/test_env.sh REPO PREFIX` sources the env scripts of the source
tree and of the install in bash, zsh, tcsh and dash (22 cases; a shell that is
not installed is skipped) and exits with its FAIL count.

## PASS criteria

- `check_reflectance.py` — the absorbed count is within 5 σ (Poisson) of N·D from the full Drude model.
- `check_planck_spectrum.py` — two rows, both required: the peak of E/kT of the emitted photons is within [0.65, 1.35] of the photon-number peak 1.5936, and a Kolmogorov-Smirnov test of the first-crossing energies against the photon-number Planck CDF truncated to the emitter band (4.14e-5 to 8.27e-2 eV) gives p > 0.01. The peak row alone passes 10 K data analysed at 9–12 K; the KS row fails a 10 % temperature error (and so a worker that missed the `config_mt` broadcast).
- `check_nreflect.py` — schema-2 per-track boundary ordinals must be contiguous and unique; cumulative reflection counts must match decoded statuses, and final termination totals must match the crossing rows. Every track must have exactly one termination. Zero-crossing world exits pass with zero totals. Legacy schemas are explicitly unsupported by this check. The plot shows final actual reflections per track, including zero. Runs on the Planck, crack-wall, world-exit and light-pipe fixtures.
- `check_angle_distribution.py` — first-hit Cu incidence angles below 15° match uniform-in-θ emission (KS p > 0.01). The test is N-sensitive: it is meant for the 10k-event `planck` case (`planck.mac`).
- `check_crack_wall_reflection.py` — every crack→Cu reflection flips p_z and keeps p_x, p_y to 1e-9.
- `check_crack_transmittance.py` — T_obs within 3 σ (binomial) of 0.50 and no tangential exits.
- `check_crack_ratio.py` — the crack2/crack1 entry ratio is within 3 σ (Poisson) of the aperture ratio A2/A1 = 1.962.
- `check_crack_oblique.py`, `check_crack_frequency.py` — all 138 and 89 checks respectively.
- `check_invariants.py` — both invariants, each printed in its own section: no photons in metal (no crossing starts inside a `Cu_RRR*` or `BBR_Perfect*` material, and the file holds at least one crossing; `--allow-no-crossings` waives only the latter, for the world-exit fixture), and termination labels (the file holds at least one `abspoints` row, no `unknown` label, every world exit is `WorldExit`, every absorption has a volume, and the `BBRAbsorb` counts agree between the two ntuples). A code with no legend entry fails the section that reads that column (`legend lacks code(s) …`), so a legend gap cannot make a check pass vacuously.
- `check_cu_serov.py` — full-Drude loss for the `OF_Cu` (RRR 3) and `HP_Cu` (RRR 6) aliases within ±10 % of Serov et al. (2016). **XFAIL:** `HP_Cu` comes out 13 % low at 230 GHz, because RRR 6 was derived with Hagen-Rubens. Whether to move `HP_Cu` to RRR 5 or accept a wider tolerance is an open decision; the runner reports the check as XFAIL, and as XPASS (a failure) if it starts passing.

One validator needs an output the fixtures do not produce, so it is run by
hand: `check_cu_absorptance.py` (0.3 < A_obs/A_theory < 3 against the
Planck-weighted Drude absorptance, on a run with enough Cu absorptions, such as
`planck_5M.mac`).

## Regression inputs

Four cases run example macros directly: `refl` runs
`examples/testworld/G4Macros/reflectance.mac`, `planck` runs
`examples/testworld/G4Macros/planck.mac`, `config_mt` runs
`examples/testworld/G4Macros/config_mt.mac` and `lp` runs
`examples/lightpipe/G4Macros/lightpipe.mac`. These macros double as regression
inputs, and their validators rely on the seeds, event counts and settings in
them. The drift guard `regression macros pinned` records, for each, the sha256
of its lines with full-line comments and blank lines stripped. It fails, naming
the macro, when those lines change or when the runner no longer runs the macro.
Only edits to full-line comments and blank lines pass: an inline `# …` after a
command, or trailing whitespace, changes the hash although Geant4 ignores it.
Change such a macro only together with its validators, then update its pin in
`Scripts/drift_guards.sh`.

The other six cases run the `Validation_*.mac` fixtures, which exist only
here. Their validators rely on the seeds and event counts too, and
`check_crack_oblique.py` and `check_crack_frequency.py` index the runs by
position, so change a fixture only together with its validators. Every
`Validation_*.mac` must also be called by the runner (drift guard
`every fixture is run`).

The library install copies this directory to `<prefix>/validation`
(`INSTALL_VALIDATION`, on by default). The four example macros are installed
only with `INSTALL_EXAMPLES`, so without it the installed `drift_guards.sh`
reports the pin guard as skipped. `cmake --install` never deletes: a file
removed from `validation/` stays in a `<prefix>/validation` installed before,
so remove that directory before reinstalling.

To run one fixture by hand, use an empty run directory, since every fixture
writes `output/` in its CWD:

```bash
mkdir -p run/transmit && cd run/transmit
<repo>/examples/testworld/build/bbrsimTestWorld <repo>/validation/G4Macros/Validation_CrackTransmit.mac >run.log 2>&1
conda run -n bbrsim python <repo>/validation/check_crack_transmittance.py output/bbr.root
```

## Dataset sidecars

**Placement and timing.** Every `<id>_<freq>GHz_Ephi={0,1}` directory pair
under `waveguides/` has a sidecar `<id>_<freq>GHz.dataset.json` beside it
(JSON, schema 1.x, agreed with Blackbody-Simulations on 2026-10-05). BBRsim
reads it when it discovers the dataset, and checks every placed crack at the
first `/run/beamOn`, before any event. A missing or inconsistent sidecar is the
fatal `BBR024`; a convention BBRsim does not implement, or a declared
cross-section that does not fit the crack solid, is the fatal `BBR025`. The
frequency-independent blocks are repeated in every frequency's file, so a
copied frequency stays self-describing, and all frequencies of one ID must
agree on them (`modes.propagating_count`, the only per-frequency field of
`modes`, excepted).

**Frames.** Canonical (p, l, g) = (propagation from entrance to exit, long,
gap), p × l = g, which is the Geant4 crack-local (x, y, z). The sampler
implements one frame: HFSS global (X, Y, Z) = (−g, +l, +p), incidence angles
of the arrival direction (k = −r), and the exit CS (Z, Y, −X) = (p, l, g). A
sidecar that declares any other frame is rejected (F4).

**Fields (schema 1.0).**

| Block | Fields |
|---|---|
| top level | `schema_version`, `dataset_id`, `frequency_ghz`, `frequency_label` (the CSV `Freq` token; the directory token plus `GHz`) |
| `provenance` | `producer`, `hand_written`, `git_commit`, `created_utc`, `job_id`, `manifest_sha256`, `aedt_version`, `project`, `design`, `object`, `pyaedt_version`, `project_sha256`, `run_design`, `material`, `model_units` |
| `frames` | `canonical`; `hfss_global_axes_in_canonical` {x, y, z}; `entrance_face` and `exit_face` {axis, side}; `exit_cs` {name, origin, x, y, z, origin_mm_global}, all in HFSS global; `exit_cs_axes_in_canonical` {x, y, z}; `pose_rule`; `entrance_outward_normal_global`; `exit_outward_normal_global` |
| `excitation` | `coordinate_system`, `incidence_convention`, `normal_entry_theta_deg`, `incident_phi_deg`, `incident_theta_deg`, `ei_v_per_m`, `incoming_power_w`, `polarization_convention` (verbatim `"Ephi=0: E_theta=1; Ephi=1: E_phi=1"`), `plane_wave_origin`, `origin_mm_global` |
| `far_field` | `coordinate_system`, `definition`, `component_basis`, `theta_deg` and `phi_deg` {min, max, count, step}, `columns`, `radiation_surface`, `points_per_key` |
| `exit_field` | `coordinate_system`, `points_in_si`, `field_in_ref_cs`, `field_components_frame`, `plane`, `resolution_mm`, `grid` {y_e, z_e}, `cross_section`, `outside_points`, `bounds_method`, `rim_points`, `points_per_key_retained`, `columns` |
| `exit_field.cross_section` | `{"shape": "rectangle", "y_e_half_m", "z_e_half_m"}`, `{"shape": "disc", "radius_m"}`, or `{"shape": "polygon", "vertices_m": [[y, z], ...]}` |
| `transmittance` | `definition`, `outgoing_power`, `incoming_includes_cos_theta` |
| `symmetry` | `mirror_l`, `mirror_g`, `end_to_end`, `rotational` |
| `boundaries` | `entrance`, `exit`, `walls` |
| `modes` (required for rectangle and disc; Rohan 2026-10-04: computed, no nulls; **final**, agreed with the HFSS session 2026-10-05) | `cutoff_ghz` and `mode`: the lowest mode, which `BBR026` compares against. `basis`: the closed-PEC model string, numbers as Python repr in metres; recorded only. `list_limit_ghz` = 20000.0, inclusive. `polarization_filter_limit_ghz`: TE01 = c/(2b) for a rectangle, the TE11 cutoff for a disc. For a rectangle, `gap_family_onsets` [{n, mode, cutoff_ghz}] (n = 0 is TE10, then TE0k at k·c/(2b)) plus `mode_count_below_limit`. For a disc, `cutoffs` [{mode, cutoff_ghz, degeneracy}], ascending. `propagating_count`: entries with cutoff ≤ this file's `frequency_ghz`, the only per-frequency field. Index pairs are counted once each (TE and TM separately; disc degeneracy not doubled). Axes: a = 2·y_e_half_m with m along l, b = 2·z_e_half_m with n along g; the lowest mode is TE10 if a ≥ b, else TE01. Ties within 1e-12 relative go TE before TM, then by index. Labels are `TE{n}{p}` when both indices are below 10, otherwise `TE{n},{p}`. All derived from the declared solid; omitted for a polygon |
| `geometry` | `shape`, `extent_mm` {p, l, g}, `bounding_box_mm` (HFSS global, [xmin, ymin, zmin, xmax, ymax, zmax]) |
| `files` | keyed by the path relative to the sidecar, e.g. `"<dataset_id>_<label>_Ephi=0/waveguide.csv"`, each {sha256, bytes, rows} |

**Checks.** A field that is missing or has the wrong JSON type is `BBR024`. A
later 1.x sidecar may carry fields this version does not know; they are
ignored. At discovery and initialization, sidecar only (no CSV is read):

- F1 `schema_version`: the major version is 1 (`BBR024`).
- F2 `dataset_id`: equals the directory ID and the stripped physical-volume name (`BBR024`).
- F3 `frequency_label`, `frequency_ghz`: the label equals the directory token plus `GHz`; `frequency_ghz` and the label's parsed value are within 0.1 % of the directory frequency (`BBR024`).
- F4 `frames`: `hfss_global_axes_in_canonical` is a right-handed orthonormal basis, `exit_cs.z` = x × y, `exit_cs` mapped through it equals `exit_cs_axes_in_canonical`, the exit and entrance outward normals map to +p and −p, and the frame is the one the sampler implements (`BBR025`).
- F5 `excitation`: `"global"`, `"arrival_direction"`, `normal_entry_theta_deg` 180 and the verbatim polarization string (`BBR025`); the incidence lists are non-empty lists of numbers (`BBR024`).
- F6 `far_field`: `"exit_cs"`, `"Theta-Phi"`, `"spherical_in_exit_cs"`, and `columns` in the positional order `BBRHFSSData` reads (`BBR025`).
- F7 `exit_field`: `"exit_cs"`, `points_in_si` true, `"x_e=0"`, `columns` in the order `BBRHFSSData` reads, and `field_in_ref_cs` agreeing with `field_components_frame` (false with `"hfss_global"`, true with `"exit_cs"`; the frame is otherwise recorded only, because BBRsim uses the components only through \|E\|² and Σ E₀·E₁*) (`BBR025`); `outside_points` is `none`, `omitted` or `zero` (`BBR024`).
- F8 `exit_field.cross_section`: `rectangle` or `disc` with positive dimensions (`BBR024` for an unknown shape or a dimension ≤ 0); `polygon` is reserved and not supported yet (`BBR025`).
- F9 `transmittance`: the supported `definition` string and `incoming_includes_cos_theta` false (`BBR025`).
- F10 `symmetry`: `mirror_l`, `mirror_g` and `end_to_end` all true, because the sampler folds by both mirrors and serves both ends from one table; `rotational` is recorded only (`BBR025`).
- F11 `geometry.extent_mm`: p, l, g > 0 (`BBR024`); compared with the placed solid's bounding limits along local x, y, z, printed as one `[BBR]` info line and stored in `metadata.json` (info only, not a warning).
- F12 `exit_field.cross_section` against the placed solid, at the first `/run/beamOn`, for every placement of the ID: the solid's local origin is inside it, and the declared section's boundary, mapped onto both exit faces (±x) under both transverse mirrors with the wrapper's axial inset, is strictly inside the solid. With `rim_points` `"included"`, the recorded Geant4 margin is what makes this pass (`BBR025`).
- F13 `modes.mode`, `.cutoff_ghz`, `.polarization_filter_limit_ghz`: re-derived from `cross_section` to 1e-6 relative with closed forms (rectangle f_mn = (c/2)·√((m/a)² + (n/b)²), c = 299792458 m/s; disc TE11 with x′₁₁ = 1.8411837813). The full lists are re-derived only by the Python validator, with scipy (`BBR025`).

At the CSV load, per frequency, in `BBRHFSSData`:

- C1 the CSV headers equal `far_field.columns` and `exit_field.columns` (`BBR013`).
- C2 every incidence key of the CSVs lies on the declared `incident_phi_deg` × `incident_theta_deg` grid. This is a subset check, not an equality: HFSS grids are complete, test grids may be sparse (`BBR007`).
- C3 per key, the far-field row count equals `points_per_key` and every Phi and Theta lies in the declared [min, max]; over the whole file the distinct Phi and Theta counts equal the declared counts. The `step` is checked by the Python validator only (`BBR012`).
- C4 X = 0 on every exit row; the exit rows per key and polarization equal `points_per_key_retained`; every Y and Z lies in the declared `grid` range, and the distinct counts equal the declared counts (`BBR012`).
- C5 every exit point, at any key and either polarization, lies inside the declared cross-section; the message gives the count and an offending point (`BBR025`). C5 checks the data against the sidecar, F12 the sidecar against the geometry; together they put every exit point inside the placed solid.

The containment tolerance is 1e-6 relative (the HFSS runner's own; never use a
stricter one). Grid and angle comparisons use 1e-9 (degrees for angles,
relative to the largest declared |value| for exit coordinates).

At run time, `BBR026` (a `JustWarning`, once per dataset and direction) fires
when a photon's frequency and the grid frequency serving it lie on opposite
sides of `modes.cutoff_ghz`; the runner tolerates it only where it is expected.

Recorded but not checked: `provenance`, `boundaries`, `geometry.shape`,
`geometry.bounding_box_mm`, `symmetry.rotational`, `pose_rule` and the origins
are copied into `<stem>.metadata.json`. The `files` checksums are verified by
the runner's `check_dataset_sidecars.py`; BBRsim itself does not verify them.

**Tools.** `validation/check_dataset_sidecars.py` re-derives the full mode
lists (scipy) and the CSV checksums. `bbrsim.sidecar` (`tools/python`) is the
Python twin of the C++ checks and holds the builders the mock generators use.
`validation/Scripts/write_legacy_sidecars.py` wrote the sidecars of the two
500 GHz datasets.

## Mock HFSS tree

`Validation_CrackFrequency.mac` tests the frequency-keyed HFSS lookup against a
mock five-frequency tree (50 / 150 / 500 / 1500 / 5000 GHz) that
`Scripts/make_mock_hfss_frequencies.py` builds from the real 500 GHz data; the
script refuses to write inside `data/`. The runner builds the tree in
`BUILD_DIR/mock_hfss` when it is missing or older than the data or the script
(about 1.3 GB), and links it into the regression directory, where the fixture's
`/bbr/dataDir ../mock_hfss` finds it. The fixture header gives the by-hand
recipe.

## Leak guard

The real `data/waveguides` must hold only the 500 GHz datasets; the runner
checks this with `bbrsim.hfss.discover_frequencies`, so mock data written into
`data/` by mistake fails the run. The runner also stops if the installed data
copy under `<prefix>/share/BBRsim/data` differs from `data/`, because
`cmake --install` never deletes stale files.

## Data note

`data/materials/cu_reflectance_measured.csv` is literature data that no script
reads; it is kept for reference.
