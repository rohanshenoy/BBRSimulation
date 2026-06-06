# Cu Material Library — Current State & Remaining Work

> **Status (updated 2026-06-05):** Drude model fully implemented in `BBRMaterials.hh`.
> Remaining work: update verification script and plot to use Drude theory instead of
> Hagen-Rubens, and update test-mac comments.

---

## Physics basis (revised)

**RRR is the sole user-facing parameter** for copper optical properties.  
Users supply an integer RRR per cryostat component; the simulation derives everything else.

| Quantity | Formula | Note |
|----------|---------|------|
| σ(4 K) | RRR × σ_RT | σ_RT = 5.96×10⁷ S/m universal for Cu |
| τ | σ_DC × mₑ / (nₑ e²) | Drude scattering time |
| σ(ω) | σ_DC / (1 − iωτ) | Full Drude AC conductivity |
| R(ω) | \|(ñ−1)/(ñ+1)\|² | Fresnel, normal incidence, Griffiths §9.4 |

**Matthiessen's rule & σ_phonon:**  
σ_phonon ∝ 1/T is only the simple power-law limit, valid for T > ~50 K.  
Below ~50 K the phonon contribution is governed by umklapp (Bloch-Grüneisen regime)
and becomes negligible compared to the impurity term at 4 K for any RRR > 1.  
→ For all BBRsim targets (detector at mK, copper shield at 4 K), use σ_DC = RRR × σ_RT.  
→ For warm shield layers (40 K, 77 K), umklapp makes the 1/T approximation unreliable;
  create a separate `G4MaterialPropertiesTable` per temperature stage using `GetCopper(RRR, T_K)`.

**Do NOT ask users to supply σ_DC directly** — it is administratively tedious and they will
not have measured values. Accept RRR as the sole required input.

**Named aliases** (OFHC_Cu, OF_Cu, HP_Cu) are convenience mappings to typical RRR values:

| Alias | RRR | Typical grade |
|-------|-----|---------------|
| OFHC_Cu | 100 | Standard OFHC copper |
| OF_Cu | 3 | Low-purity or cold-worked |
| HP_Cu | 6 | Hydrogen-annealed, 99.999% |

These were back-calculated to roughly reproduce Serov et al. (2016) measured
absorptance at 4 K, but the RRR is the physical parameter — not the σ.

---

## ✅ DONE — implemented in BBRMaterials.hh

- `BuildDrudeMaterial(name, Z, A, ρ, RRR, T_K)` — internal helper; full Drude reflectance table
- `GetCopper(RRR, T_K=4.0)` — public API; name = `Cu_RRR{N}_T{T}K`; cached by G4 material store
- `GetCopperByName(name)` — maps OFHC_Cu / OF_Cu / HP_Cu → `GetCopper(RRR, 4.0)`
- `GetVacuumWG()`, `GetPerfectAbsorber()`, `GetPerfectReflector()` — ancillary materials
- Below 50 K: σ_DC = RRR × σ_RT (impurity dominated)
- Above 50 K: Matthiessen's rule with σ_phonon ≈ σ_RT × 273/T (simple 1/T approximation)

## ✅ DONE — BBRTestDetectorConstruction

`SetCuMaterial(name)` messenger wired; `/bbr/det/setCuMaterial <name>` selects material
before `/run/initialize`.

## ✅ DONE — test macs

`test.mac`, `test_of_cu.mac`, `test_hp_cu.mac` all exist and exercise the messenger.
Comments still reference old Serov/Hagen-Rubens framing — see Task 1 below.

---

## Remaining tasks

---

### Task 1 — Update check_cu_absorptance.py: --rrr flag, Drude theory

The current script uses `--sigma` and computes theory via Hagen-Rubens:
`A_HR(ω) = 2√(2ε₀ω/σ)`.
The simulation uses the full Drude model, so the theory comparison must also use Drude.
Mismatch causes the ratio A_obs/A_theory to be wrong for frequencies above
the H-R validity cutoff (~64 GHz for RRR=100 at 4K).

**Files:** `scripts/check_cu_absorptance.py`

- [ ] Replace `--sigma` argument with `--rrr` (integer, default 100).
      Derive σ_DC = RRR × σ_RT inside the script.
- [ ] Replace `A_HR(omega)` integrand with `A_Drude(omega)` = 1 − R_Drude(omega),
      where R_Drude is the full Griffiths §9.4 formula (same as in BBRMaterials.hh).
- [ ] Update usage docs and expected outputs.
- [ ] Verify: run against `of_cu_out.txt` with `--rrr 3` and `test` output with `--rrr 100`.

---

### Task 2 — Update plot_cu_reflectance.py: replace H-R aliases with Drude curves

The plot currently shows Hagen-Rubens curves labelled "OF_Cu (Serov 2016, Hagen-Rubens)"
and "HP_Cu (Serov 2016, Hagen-Rubens)". These should be replaced with Drude curves for
the corresponding RRR values so the plot accurately shows what the simulation does.

**Files:** `scripts/plot_cu_reflectance.py`

- [ ] Replace `materials_HR` dict entries for OF_Cu/HP_Cu with Drude entries
      (RRR=3, RRR=6) using the same `drude_R(freq_Hz, sigma_drude(4.0, RRR))` calls
      already used for OFHC_Cu.
- [ ] Update legend labels: "Cu_RRR3_T4K (Drude)" etc.
- [ ] Add an annotation note on panel 3 (temperature dependence) marking the ~50 K
      boundary where the simple 1/T phonon approximation breaks down (umklapp regime below).
- [ ] Keep Serov reference data points as experimental validation markers.

---

### Task 3 — Update test mac comments

`test_of_cu.mac` and `test_hp_cu.mac` still reference Serov σ_eff values.
Update the header comments to reflect the Drude/RRR model.

**Files:** `test_of_cu.mac`, `test_hp_cu.mac`

- [ ] Replace "sigma_eff = 1/rho_0" comment with "RRR=3, Drude model at 4K".
- [ ] Replace "HP annealed, sigma_eff back-calc" with "RRR=6, Drude model at 4K".

---

## Self-review checklist

- [ ] `check_cu_absorptance.py` theory uses Drude (matches simulation)
- [ ] `--rrr` flag derives σ from RRR × σ_RT (not hardcoded)
- [ ] `plot_cu_reflectance.py` shows Drude curves for each alias
- [ ] Serov data points remain as experimental reference markers
- [ ] Test mac comments no longer reference Hagen-Rubens or Serov σ_eff
- [ ] Build passes with no changes to C++ (plan is Python-only)
