# BBRsim validation

Fixed-seed fixtures and PASS/FAIL validators that gate every change to BBRsim.

```
G4Macros/Validation_*.mac   the seven validation-only fixtures
check_*.py                  the 14 validators
Scripts/run_regression.sh   build + fixtures + validators, one command
Scripts/drift_guards.sh     source-tree consistency checks (run by the runner, or alone)
Scripts/make_mock_hfss_frequencies.py   the mock HFSS tree for Validation_CrackFrequency
Scripts/make_mock_round_gap.py          the mock RoundGap_r50um dataset at 2000 GHz for Validation_RoundGap, with links to the real crack datasets
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
`Scripts/tests/test_env.sh` (`env scripts`). It runs the eleven cases in parallel,
each in `BUILD_DIR/regression/<case>/`, with `G4FORCENUMBEROFTHREADS` pinned to
8 (`BBR_THREADS` overrides; Geant4 warns when a run has more threads than the
square root of its events, and the log scan would fail on a many-core machine),
and runs the validators on their output. Seven cases run the fixtures here; four
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
| `round` | `Validation_RoundGap.mac` | `bbrsimTestWorld` | `check_round_gap.py --data-dir mock_round_gap --log run.log` on `output/`; `check_invariants.py` on each of the 9 `bbr_round_rNN.root` |
| `lp` | `examples/lightpipe/G4Macros/lightpipe.mac` | `bbrsimLightPipe` | `check_invariants.py` |
| `lp_cad` | `Validation_LightPipeCAD.mac` | `bbrsimLightPipe` | `check_invariants.py` (cad mode, the bundled `box_sample.stl` through `BBRSIMDATA`) |
| — | (no ROOT input) | — | `check_cu_serov.py` (XFAIL); `check_dataset_sidecars.py` on `data/waveguides`, `BUILD_DIR/mock_hfss/waveguides` and `BUILD_DIR/mock_round_gap/waveguides` |

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
- `check_round_gap.py` — all 49 checks over the nine runs of `Validation_RoundGap.mac` (2000 GHz, against the mock `RoundGap_r50um` table): one gap entry per event along the gun direction; T_obs within 3 σ (binomial) of the `bbrsim.hfss` prediction; every exit on the exit face and within the 50 µm HFSS radius, although the Geant4 hole is 51 µm; for the three fixed polarizations, the share of exits inside R/2 within 4 σ of the table's prediction, which tells the Ephi=0 and Ephi=1 radial profiles apart; and the mean exit direction within 4 standard errors of the prediction in each component. The last run, at the diagonal key (45°, 135°) with a fixed polarization of equal θ̂ and φ̂ components, makes the polarization cross term 2 E_θ E_φ √(T₀T₁) Re ρ large (Re ρ = 0.288 in the mock): T = 0.511 against 0.289 with the term's sign flipped, 28 σ apart, so it pins that sign end to end; the random-polarization diagonal run averages the term out. With `--log` (as the runner calls it) one more check: the run log holds exactly one `[BBR] crack RoundGap_r50um:` line reporting its sidecar fit, the startup check of the placed solid against the sidecar; without it, 48 checks.
- `check_invariants.py` — both invariants, each printed in its own section: no photons in metal (no crossing starts inside a `Cu_RRR*` or `BBR_Perfect*` material, and the file holds at least one crossing; `--allow-no-crossings` waives only the latter, for the world-exit fixture), and termination labels (the file holds at least one `abspoints` row, no `unknown` label, every world exit is `WorldExit`, every absorption has a volume, and the `BBRAbsorb` counts agree between the two ntuples). A code with no legend entry fails the section that reads that column (`legend lacks code(s) …`), so a legend gap cannot make a check pass vacuously.
- `check_cu_serov.py` — full-Drude loss for the `OF_Cu` (RRR 3) and `HP_Cu` (RRR 6) aliases within ±10 % of Serov et al. (2016). **XFAIL:** `HP_Cu` comes out 13 % low at 230 GHz, because RRR 6 was derived with Hagen-Rubens. Whether to move `HP_Cu` to RRR 5 or accept a wider tolerance is an open decision; the runner reports the check as XFAIL, and as XPASS (a failure) if it starts passing.
- `check_dataset_sidecars.py` — every dataset in the real tree (`data/waveguides`) and the two mock trees (`BUILD_DIR/mock_hfss/waveguides`, `BUILD_DIR/mock_round_gap/waveguides`) has a complete trio and a sidecar that passes `bbrsim.sidecar.check_full` (checks, full mode lists, CSV checksums, C1-C5), and the frequencies of each ID agree on the frequency-independent blocks; each directory must hold at least one dataset. See [Dataset sidecars](#dataset-sidecars).

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

The other seven cases run the `Validation_*.mac` fixtures, which exist only
here. Their validators rely on the seeds and event counts too, and
`check_crack_oblique.py`, `check_crack_frequency.py` and `check_round_gap.py`
index the runs by position, so change a fixture only together with its
validators. Every `Validation_*.mac` must also be called by the runner (drift guard
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
first run initialization (`/run/initialize` with a multithreaded or task run
manager, as in the examples; `/run/beamOn` with a sequential one), before any
event. A missing or inconsistent sidecar is the fatal `BBR024`; a convention
BBRsim does not implement, or a declared cross-section that does not fit the
crack solid, is the fatal `BBR025`; a missing `far_field.csv` or
`waveguide.csv` in either Ephi directory is the fatal `BBR001` or `BBR002`
(the CSVs themselves load at the first photon that selects the frequency).
Every placed `vacuum_wg` volume is validated, whether or not a photon reaches
it, so a data tree named with `/bbr/dataDir` or `BBRSIMDATA` must hold a
complete dataset with sidecars for each placed crack; the test world always
places both `InfParallelPlate_crack1Rohan` and `InfParallelPlate_crack2`. The
frequency-independent blocks are repeated in every frequency's file, so a
copied frequency stays self-describing. Every frequency of one ID must agree on
the frame mapping (`frames.hfss_global_axes_in_canonical` and
`frames.exit_cs_axes_in_canonical`), `symmetry`, `boundaries`, `geometry.shape`
and `geometry.extent_mm`, the `modes` block without `basis` and
`propagating_count`, and `exit_field.cross_section`: numbers to 1e-9 relative,
everything else exactly, over the fields both sidecars carry. Descriptive and
pose-dependent fields (`frames.canonical`, `pose_rule`, the face selectors, the
origins, `geometry.bounding_box_mm`) may differ, so a hand-written legacy
sidecar and one written by the HFSS runner can serve one ID.

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
| `modes` (required for rectangle and disc; computed, no nulls) | `cutoff_ghz` and `mode`: the lowest mode, which `BBR026` compares against. `basis`: the closed-PEC model string, numbers as Python repr in metres; recorded only. `list_limit_ghz` = 20000.0, inclusive. `polarization_filter_limit_ghz`: TE01 = c/(2b) for a rectangle, the TE11 cutoff for a disc. For a rectangle, `gap_family_onsets` [{n, mode, cutoff_ghz}] (n = 0 is TE10, then TE0k at k·c/(2b)) plus `mode_count_below_limit`. For a disc, `cutoffs` [{mode, cutoff_ghz, degeneracy}], ascending. `propagating_count`: entries with cutoff ≤ this file's `frequency_ghz`, the only per-frequency field. Index pairs are counted once each (TE and TM separately; disc degeneracy not doubled). Axes: a = 2·y_e_half_m with m along l, b = 2·z_e_half_m with n along g; the lowest mode is TE10 if a ≥ b, else TE01. Ties within 1e-12 relative go TE before TM, then by index. Labels are `TE{n}{p}` when both indices are below 10, otherwise `TE{n},{p}`. All derived from the declared solid; omitted for a polygon |
| `geometry` | `shape`, `extent_mm` {p, l, g}, `bounding_box_mm` (HFSS global, [xmin, ymin, zmin, xmax, ymax, zmax]) |
| `files` | keyed by the path relative to the sidecar, e.g. `"<dataset_id>_<label>_Ephi=0/waveguide.csv"`, each {sha256, bytes, rows}; `rows` is the number of non-blank lines after the header (a blank line holds only whitespace), which is the pandas data-row count `len(read_csv(...))` with or without a trailing newline |

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
- F11 `geometry.extent_mm`: p, l, g > 0 (`BBR024`); compared with the placed solid's bounding limits along local x, y, z, printed as one `[BBR]` info line and stored as `hfss_extent_mm` and `geant4_extent_mm` in `data.hfss_datasets` of each `<stem>.metadata.json` (info only, not a warning).
- F12 `exit_field.cross_section` against the placed solid, at the first run initialization, for every placement of the ID: the solid's local origin is inside it, and the declared section's boundary, mapped onto both exit faces (±x) under both transverse mirrors with the wrapper's axial inset, is strictly inside the solid. With `rim_points` `"included"`, the recorded Geant4 margin is what makes this pass (`BBR025`).
- F13 `modes.mode`, `.cutoff_ghz`, `.polarization_filter_limit_ghz`: re-derived from `cross_section` to 1e-6 relative with closed forms (rectangle f_mn = (c/2)·√((m/a)² + (n/b)²), c = 299792458 m/s; disc TE11 with x′₁₁ = 1.8411837813). The full lists are re-derived only by the Python validator, with scipy (`BBR025`).

At the CSV load, per frequency, in `BBRHFSSData`:

- C1 the CSV headers equal `far_field.columns` and `exit_field.columns` (`BBR013`).
- C2 every incidence key of the CSVs lies on the declared `incident_phi_deg` × `incident_theta_deg` grid. This is a subset check, not an equality: HFSS grids are complete, test grids may be sparse. Within one polarization, `far_field.csv` and `waveguide.csv` carry the same key set (`BBR007` for either).
- C3 per key, the far-field row count equals `points_per_key` and every Phi and Theta lies in the declared [min, max]; over the whole file the distinct Phi and Theta counts equal the declared counts (`BBR012`). Only the Python validator checks the `step`: when `count` > 1 it must equal (max − min)/(count − 1) of the data to 1e-9 relative (`BBR012`).
- C4 X = 0 on every exit row; the exit rows per key and polarization equal `points_per_key_retained`; every Y and Z lies in the declared `grid` range, and the distinct counts equal the declared counts (`BBR012`).
- C5 every exit point, at any key and either polarization, lies inside the declared cross-section; the message gives the count and an offending point (`BBR025`). With `outside_points` `"zero"`, a point outside the section is exempt when its six field components are exactly zero in both polarizations; with `"none"` or `"omitted"` every point must lie inside. C5 checks the data against the sidecar, F12 the sidecar against the geometry; together they put every exit point that carries field inside the placed solid.

C1-C5 run on both polarizations, and the two must agree: the Ephi=1 CSVs have
the same incidence keys, far-field and exit grids and per-key row counts as the
Ephi=0 CSVs (`BBR012`; `BBRHFSSData` also pairs the rows by position). A CSV
row with the wrong number of fields, or an empty field, is `BBR013`.

The containment tolerance is 1e-6 relative (the HFSS runner's own; never use a
stricter one). Incidence keys are rounded to 0.01° before the C2 subset test,
as `BBRHFSSData` rounds them. The far-field ranges use 1e-9 degrees and the
exit-grid ranges 1e-9 relative to the largest declared |value|.

At run time, `BBR026` (a `JustWarning`, once per dataset and direction) fires
when a photon's frequency and the grid frequency serving it lie on opposite
sides of `modes.cutoff_ghz`; the runner tolerates it only where it is expected.

Recorded but not checked against the conventions BBRsim implements (of these,
only `boundaries`, `geometry.shape` and `symmetry.rotational` enter the
agreement across frequencies): `provenance`, `boundaries`, `geometry.shape`,
`geometry.bounding_box_mm`, `symmetry.rotational`, `pose_rule` and the origins
(`frames.exit_cs.origin_mm_global`, `excitation.origin_mm_global`) are copied
verbatim into `<stem>.metadata.json`, a missing one as `null`. They go in
`data.hfss_datasets`, one entry per placed crack volume in the order they were
validated: `volume`, `dataset_id`, `geant4_extent_mm` (the placed solid's
extent along local x, y, z), `hfss_extent_mm` (`geometry.extent_mm`) and
`frequencies`, each with its `label`, `frequency_ghz`, `sidecar` file name and
`recorded` fields. `recorded` is one object, `{"provenance", "boundaries",
"geometry": {"shape", "bounding_box_mm"}, "symmetry": {"rotational"},
"frames": {"pose_rule", "exit_cs_origin_mm_global"}, "excitation":
{"origin_mm_global"}}`, each value taken from the sidecar field of the same
path except `frames.exit_cs_origin_mm_global`, a flattened key holding the
sidecar's `frames.exit_cs.origin_mm_global`. A geometry without cracks records
`[]`. The `files` checksums are verified by the runner's
`check_dataset_sidecars.py`; BBRsim itself does not verify them.

**Tools.** `validation/check_dataset_sidecars.py <waveguides-dir> [...]`
(default: `waveguides/` under the data root of `bbrsim.paths.data_dir()`)
checks every dataset in the given directories. A dataset `<id>_<freq>GHz` is
the trio `<stem>_Ephi=0/`, `<stem>_Ephi=1/` and `<stem>.dataset.json`; other
files (a tree manifest, `SHA256SUMS`) are ignored, and a stem missing any part
of the trio fails. A complete dataset is checked with
`bbrsim.sidecar.check_full`: F1-F11 and F13, the full mode lists re-derived
with scipy, the sha256, size and row count of the four CSVs, and C1-C5 on both
polarizations. The frequencies of one ID in one directory must also agree on
the frequency-independent blocks (`bbrsim.sidecar.invariant_diff` names the
blocks that differ). It prints one `PASS` or `FAIL` line per dataset, after the
agreement check, and a `FAIL` line for a directory argument that cannot be
listed or holds no dataset; then the count and `RESULT: PASS` or
`RESULT: FAIL`, and exits 0 or 1. `bbrsim.sidecar` (`tools/python`) is the Python twin of the C++
checks and holds the builders the mock generators use.
`validation/Scripts/write_legacy_sidecars.py` wrote the sidecars of the two
500 GHz datasets in `data/waveguides/` from their CSVs and the Blackbody-Simulations
reference configs; only `boundaries.walls` (PEC) is inferred, as
`provenance.inferred` records. It refuses to write when a CSV differs from the
checksums it records; it verifies every checksum and checks every built
sidecar before it writes any, so a refusal leaves `data/` untouched.

## Mock HFSS tree

`Validation_CrackFrequency.mac` tests the frequency-keyed HFSS lookup against a
mock five-frequency tree (50 / 150 / 500 / 1500 / 5000 GHz) that
`Scripts/make_mock_hfss_frequencies.py` builds from the real 500 GHz data; the
script refuses to write inside `data/`. Each mock frequency gets its own
sidecar `<id>_<f>GHz.dataset.json`, a copy of the 500 GHz one with the
frequency, `modes.propagating_count`, the CSV checksums and the provenance
(`mock_of`, `transmittance_scale`, `far_field_theta_above_90_zeroed`)
rewritten, so the tree passes BBRsim's discovery checks and
`check_dataset_sidecars.py`. The runner builds the tree in
`BUILD_DIR/mock_hfss` when it is missing or older than the data, the script,
`tools/python/bbrsim/sidecar.py` or `tools/python/bbrsim/hfss.py` (about
1.3 GB), and links it into the
regression directory, where the fixture's `/bbr/dataDir ../mock_hfss` finds it.
The fixture header gives the by-hand recipe.

## Mock round-gap tree

`Validation_RoundGap.mac` turns on the opt-in straight round gap
(`/bbr/testworld/roundGap true`: `RoundGap_r50um`, a 51 µm-radius `vacuum_wg`
hole along x through a 0.4 mm Cu plate at z = −80 mm) and feeds it the mock
dataset `RoundGap_r50um_2000GHz` that `Scripts/make_mock_round_gap.py` writes:
a 50 µm disc in sidecar schema 1.0 with synthetic tables whose transmittance
and radial exit profile differ between the two polarizations and whose
far-field lobe moves with the incidence key, so each mapping is visible to
`check_round_gap.py`. The test world places both slab cracks as well and
BBRsim checks every placed crack before the first event, so the tree also links the real crack datasets of
`data/waveguides`, with copies of their sidecars. The script refuses to write
inside the data tree it links. The runner builds the tree (a few MB) in
`BUILD_DIR/mock_round_gap` when it is missing or older than the data, the
script, `tools/python/bbrsim/sidecar.py` or `tools/python/bbrsim/hfss.py`, and links it into the regression
directory, where the fixture's `/bbr/dataDir ../mock_round_gap` finds it. By
hand:

```bash
conda run -n bbrsim python validation/Scripts/make_mock_round_gap.py \
    --real data/waveguides --dst D/../mock_round_gap
cd D && bbrsimTestWorld <repo>/validation/G4Macros/Validation_RoundGap.mac > run.log 2>&1
conda run -n bbrsim python <repo>/validation/check_round_gap.py output --data-dir ../mock_round_gap --log run.log
```

## Leak guard

The real `data/waveguides` must hold only the 500 GHz datasets; the runner
checks this with `bbrsim.hfss.discover_frequencies`, so mock data written into
`data/` by mistake fails the run. The runner also stops if the installed data
copy under `<prefix>/share/BBRsim/data` differs from `data/`, because
`cmake --install` never deletes stale files.

## Data note

`data/materials/cu_reflectance_measured.csv` is literature data that no script
reads; it is kept for reference.
