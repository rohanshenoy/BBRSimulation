# BBRsim validation

Fixed-seed fixtures and PASS/FAIL validators that gate every change to BBRsim.

```
G4Macros/Validation_*.mac   the eight fixtures
check_*.py                  the 14 validators
Scripts/run_regression.sh   build + fixtures + validators, one command
Scripts/drift_guards.sh     source-tree consistency checks (run by the runner, or alone)
Scripts/make_mock_hfss_frequencies.py   the mock HFSS tree for Validation_CrackFrequency
Scripts/consumer_smoke/     external find_package(BBRsim) + link smoke test
```

## Running

```bash
validation/Scripts/run_regression.sh [BUILD_DIR]     # default: build
```

The runner builds and installs the library (prefix `BBR_PREFIX`, default
`install/`), sources the installed env script, builds both examples, runs the
drift guards and the consumer smoke test, runs the fixtures in parallel, each in
`BUILD_DIR/regression/<case>/`, and runs the validators on their output. The
exit code is the number of unexpected failures; a green run ends with the line
`pass=41  fail=0  xfail=1  xpass=0` (two spaces between fields, then the output
directory).

| Case | Fixture | Executable | Validators |
|---|---|---|---|
| `refl` | `Validation_Reflectance.mac` | `bbrsimTestWorld` | `check_reflectance.py --root`, `check_term_status.py`, `check_no_photons_in_metal.py` |
| `planck` | `Validation_Planck.mac` | `bbrsimTestWorld` | `check_planck_spectrum.py --temp 4`, `check_nreflect.py`, `check_angle_distribution.py`, `check_term_status.py`, `check_no_photons_in_metal.py` |
| `wall` | `Validation_CrackWall.mac` | `bbrsimTestWorld` | `check_crack_wall_reflection.py`, `check_no_photons_in_metal.py`, `check_term_status.py` |
| `exit` | `Validation_WorldExit.mac` | `bbrsimTestWorld` | `check_term_status.py` |
| `transmit` | `Validation_CrackTransmit.mac` | `bbrsimTestWorld` | `check_crack_transmittance.py`, `check_no_photons_in_metal.py`, `check_term_status.py` |
| `oblique` | `Validation_CrackOblique.mac` | `bbrsimTestWorld` | `check_crack_oblique.py` on `output/`; `check_no_photons_in_metal.py` and `check_term_status.py` on each of the 16 `bbr_oblique_rNN.root` |
| `frequency` | `Validation_CrackFrequency.mac` | `bbrsimTestWorld` | `check_crack_frequency.py --data-dir mock_hfss` on `output/`; the two invariants on each of the 16 `bbr_freq_rNN.root` |
| `lp` | `Validation_LightPipe.mac` | `bbrsimLightPipe` | `check_no_photons_in_metal.py`, `check_term_status.py` |
| — | (no ROOT input) | — | `check_physics.py`, `check_cu_serov.py` (XFAIL) |

A fixture run passes when the binary exits 0, writes its ROOT output, and its
log holds no `GeomNav`, `G4Exception`, `BBR0xx` or `LP002` line (the frequency
case tolerates its `BBR008` clamp warnings, and exactly one per side is
required).

## PASS criteria

- `check_reflectance.py` — the absorbed count is within 5 σ (Poisson) of N·D from the full Drude model.
- `check_planck_spectrum.py` — the peak of E/kT of the emitted photons is within [0.65, 1.35] of the photon-number peak 1.5936.
- `check_nreflect.py` — n_reflect = 1 is the modal bin and the counts do not rise for n = 1…10.
- `check_angle_distribution.py` — first-hit Cu incidence angles below 15° match uniform-in-θ emission (KS p > 0.01). The test is N-sensitive: it is meant for the 10k-event Planck fixture.
- `check_crack_wall_reflection.py` — every crack→Cu reflection flips p_z and keeps p_x, p_y to 1e-9.
- `check_crack_transmittance.py` — T_obs within 3 σ (binomial) of 0.50 and no tangential exits.
- `check_crack_oblique.py`, `check_crack_frequency.py` — all 138 and 89 checks respectively.
- `check_no_photons_in_metal.py` — no crossing starts inside a `Cu_RRR*` or `BBR_Perfect*` material, and the file holds at least one crossing.
- `check_term_status.py` — no `unknown` label, every world exit is `WorldExit`, every absorption has a volume, and the `BBRAbsorb` counts agree between the two ntuples.
- `check_physics.py` — the `bbrsim` formulas (Drude absorptance, Hagen-Rubens, the Planck peak, the HFSS frequency-selection rule, `paths.data_dir`) match documented reference values.
- `check_cu_serov.py` — full-Drude loss for the `OF_Cu` (RRR 3) and `HP_Cu` (RRR 6) aliases within ±10 % of Serov et al. (2016). **XFAIL:** `HP_Cu` comes out 13 % low at 230 GHz, because RRR 6 was derived with Hagen-Rubens. Whether to move `HP_Cu` to RRR 5 or accept a wider tolerance is an open decision; the runner reports the check as XFAIL, and as XPASS (a failure) if it starts passing.

Two validators need an output the fixtures do not produce, so they are run by
hand: `check_crack_ratio.py` (crack2/crack1 entry ratio within 3 σ of the
aperture ratio, on a Planck run of at least 1M events, such as
`planck_10K.mac`) and `check_cu_absorptance.py`
(0.3 < A_obs/A_theory < 3 against the Planck-weighted Drude absorptance, on a
run with enough Cu absorptions, such as `planck_5M.mac`).

## Fixtures are frozen copies

`Validation_Reflectance.mac`, `Validation_Planck.mac` and
`Validation_LightPipe.mac` are frozen copies of the example macros
`reflectance.mac`, `planck.mac` and `lightpipe.mac`. A change to an example
macro does not change a fixture. The validators rely on the fixtures' seeds and
event counts, and `check_crack_oblique.py` and `check_crack_frequency.py` index
the runs by position, so change a fixture only together with its validators.
Every `Validation_*.mac` must also be called by the runner (drift guard
`every fixture is run`).

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
