# BBRsim validation

Fixed-seed fixtures and PASS/FAIL validators that gate every change to BBRsim.

```
G4Macros/Validation_*.mac   the six validation-only fixtures
check_*.py                  the 12 validators
Scripts/run_regression.sh   build + fixtures + validators, one command
Scripts/drift_guards.sh     source-tree consistency checks (run by the runner, or alone)
Scripts/make_mock_hfss_frequencies.py   the mock HFSS tree for Validation_CrackFrequency
Scripts/consumer_smoke/     external find_package(BBRsim) + link smoke test
Scripts/tests/              make_bad_output.py (synthetic outputs for the validator negative tests),
                            test_env.sh (the env scripts in bash, zsh, tcsh, dash), check_links.sh
```

## Running

```bash
validation/Scripts/run_regression.sh [BUILD_DIR]     # default: build
```

The runner builds and installs the library (prefix `BBR_PREFIX`, default
`install/`), sources the installed env script, builds both examples, runs the
drift guards and the consumer smoke test, runs the ten cases in parallel, each
in `BUILD_DIR/regression/<case>/`, and runs the validators on their output.
Six cases run the fixtures here; four run example macros directly (see
[Regression inputs](#regression-inputs)). The exit code is the number of
unexpected failures; a green run ends with the line
`pass=47  fail=0  xfail=1  xpass=0` (two spaces between fields, then the output
directory).

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
- `check_nreflect.py` — n_reflect = 1 is the modal bin and the counts do not rise for n = 1…10.
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
