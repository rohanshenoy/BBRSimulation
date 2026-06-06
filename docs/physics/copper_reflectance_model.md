# Copper Reflectance Model in BBRsim

BBRsim simulates free-space blackbody radiation from 50 GHz to 20 THz. In this
band, copper reflectance is controlled by the Drude model of free-electron
conduction, not by room-temperature optical constants or simplified empirical
fits. This document describes the physics, the parameterization, and the
BBRsim implementation.

---

## Why Not Hagen-Rubens?

The Hagen-Rubens formula is widely used as a first estimate of metal reflectance:

```
D = 1 − R ≈ 2√(2ε₀ω/σ)
```

It is the **low-frequency limit** of the Drude model, valid only when ωτ << 1,
i.e., when the photon frequency is well below the electron scattering rate 1/τ.

The scattering time τ is set by the DC conductivity:

```
τ = σ_DC · mₑ / (nₑ · e²)
```

For OFHC copper (RRR = 100) at 4 K:

| Quantity | Value |
|---|---|
| σ_DC = 100 × σ_RT | 5.96×10⁹ S/m |
| τ | ≈ 2.5 ps |
| H-R breakdown frequency 1/(2πτ) | **≈ 64 GHz** |

BBRsim's simulation range starts at 50 GHz. Hagen-Rubens is wrong for
essentially the entire range. The full Drude model must be used.

---

## The Drude–Griffiths Model

From Griffiths, *Introduction to Electrodynamics* §9.4, for a plane wave in a
conductor with conductivity σ(ω):

**Complex wave vector:**

```
k = ω√(ε₀μ₀/2) · √[ √(1 + (σ_r/ε₀ω)²) + 1 ]
κ = ω√(ε₀μ₀/2) · √[ √(1 + (σ_r/ε₀ω)²) − 1 ]
```

where σ_r = Re[σ(ω)] is the real part of the AC conductivity.

**Complex refractive index:**

```
ñ = n + iκ = c(k + iκ)/ω
```

**Normal-incidence reflectance:**

```
R = |(ñ − 1)/(ñ + 1)|²  =  [(n−1)² + κ²] / [(n+1)² + κ²]
```

**Full Drude AC conductivity** (replaces static σ in Griffiths's equations):

```
σ(ω) = σ_DC / (1 − iωτ)
```

At low frequency (ωτ << 1): σ(ω) → σ_DC and the Hagen-Rubens limit is
recovered. At high frequency (ωτ >> 1): σ(ω) → σ_DC/(−iωτ) → purely imaginary,
reducing reflectance below the H-R prediction.

---

## Temperature and RRR Parameterization

### RRR — the primary user parameter

The Residual Resistance Ratio is the single number characterizing copper quality:

```
RRR = ρ(273 K) / ρ(4 K)  =  σ(4 K) / σ(273 K)
```

At room temperature, resistivity is dominated by phonon scattering — the same
for all copper grades. At cryogenic temperatures, phonons are largely frozen out
and resistivity is set by impurity and defect scattering, which varies with
purity. RRR captures this ratio directly.

**RRR is the only parameter users need to supply.** The simulation derives
σ_DC, τ, and the full REFLECTIVITY table internally. Users never supply σ
directly.

### Conductivity vs temperature: Matthiessen's rule

```
1/σ(T) = 1/σ_impurity + 1/σ_phonon(T)
```

where:
- `σ_impurity = RRR × σ_RT` — fixed by impurities/defects, independent of T
- `σ_phonon(T) ≈ σ_RT × 273/T` — simple power-law approximation, valid for T > ~50 K

**At T = 4 K:** phonon term is negligible, so:

```
σ(4 K) ≈ RRR × σ_RT
```

This is an excellent approximation for all copper grades at the detector stage.

**At 10–50 K** (cryostat shield layers): the simple 1/T phonon model breaks
down because phonon-phonon umklapp scattering enters the Bloch-Grüneisen
regime. The code switches to the impurity-dominated formula below 50 K rather
than extrapolating the inaccurate 1/T approximation.

**At T ≥ 50 K** (outer shields, room temperature): Matthiessen's rule with the
1/T phonon term is used. For warm shield layers users should create a separate
material with the appropriate stage temperature (`/bbr/det/setCuStageT <T> K`).

---

## Copper Material Constants

| Constant | Symbol | Value | Units |
|---|---|---|---|
| Room-temp conductivity | σ_RT | 5.96×10⁷ | S/m |
| Free electron density | nₑ | 8.49×10²⁸ | m⁻³ |
| Electron mass | mₑ | 9.109×10⁻³¹ | kg |
| Electron charge | e | 1.602×10⁻¹⁹ | C |
| Permittivity of free space | ε₀ | 8.854×10⁻¹² | F/m |

σ_RT is universal for all copper grades. It does not depend on purity and users
should never override it.

---

## Typical RRR Values

| Grade | RRR | σ(4 K) [S/m] | τ(4 K) [ps] | H-R valid below |
|---|---|---|---|---|
| Cold-worked / disordered | ~1–5 | ~0.6–3×10⁸ | ~0.03–0.13 ps | > 1 THz |
| Commercial Cu | ~10 | 5.96×10⁸ | 0.25 ps | ~640 GHz |
| OFHC Cu | ~100 | 5.96×10⁹ | 2.5 ps | ~64 GHz |
| Ultra-pure crystal | ~500 | 3.0×10¹⁰ | 12.5 ps | ~13 GHz |

For all grades the Hagen-Rubens validity cutoff is inside or below the BBRsim
range (50 GHz–20 THz). **The full Drude model is always required.**

---

## BBRsim Implementation

The implementation lives in `include/BBRMaterials.hh` (header-only, compiled
into any translation unit that includes it).

### Internal builder

```cpp
// Internal — not part of the public API.
// Builds a G4Material with a 20-point log-spaced REFLECTIVITY table.
G4Material* BuildDrudeMaterial(const G4String& name,
                               G4double Z, G4double A_g_mol, G4double density_g_cm3,
                               G4int RRR, G4double T_K);
```

Physics path inside `BuildDrudeMaterial`:
1. Compute σ_DC: uses RRR × σ_RT if T_K < 50 K; Matthiessen's rule otherwise.
2. Compute τ = σ_DC × mₑ / (nₑ e²).
3. For each of 20 log-spaced energies (50 GHz → 20 THz):
   - Evaluate Re[σ(ω)] = σ_DC / (1 + ω²τ²)
   - Compute k, κ from Griffiths §9.4
   - Compute R = |(ñ−1)/(ñ+1)|²
4. Store as `REFLECTIVITY` on the material's `G4MaterialPropertiesTable`.

### Public API

```cpp
// Parameterized by RRR and temperature stage.
// T_K defaults to 4.0 K (cryogenic baseline).
// Material name is deterministic ("Cu_RRR{N}_T{T}K") — cached by G4 material store.
G4Material* BBRMaterials::GetCopper(G4int RRR, G4double T_K = 4.0);

// Named aliases for convenience — resolve to GetCopper(RRR, 4.0).
// Valid names: "OFHC_Cu" (RRR=100), "OF_Cu" (RRR=3), "HP_Cu" (RRR=6).
G4Material* BBRMaterials::GetCopperByName(const G4String& name);
```

Usage examples:

```cpp
// Standard OFHC copper at the 4K stage
G4Material* mat = BBRMaterials::GetCopper(100);

// User's specific sample with measured RRR=250
G4Material* mat = BBRMaterials::GetCopper(250);

// 4K plate shield at a warmer stage (40K radiation shield)
G4Material* mat = BBRMaterials::GetCopper(50, 40.0);

// Named alias
G4Material* mat = BBRMaterials::GetCopperByName("OFHC_Cu");
```

### Messenger interface

In a Geant4 macro file, before `/run/initialize`:

```mac
# Named alias — most common
/bbr/det/setCuMaterial OFHC_Cu

# Direct RRR — when you know your sample's spec sheet value
/bbr/det/setCuRRR 250

# Warm shield layer
/bbr/det/setCuStageT 40 K
/bbr/det/setCuRRR 50
```

The simulation prints the resolved material name at geometry construction time:

```
[BBR] Cu wall material: Cu_RRR250_T4K  (RRR=250, T=4 K)
```

---

## Frequency Regimes

| Regime | Condition | Model | Notes |
|---|---|---|---|
| Hagen-Rubens | f << 1/(2πτ) | D ≈ 2√(2ε₀ω/σ) | Low-frequency limit only; plotted as reference |
| Drude transition | f ~ 1/(2πτ) | Full Griffiths §9.4 | BBRsim operating range (50 GHz–20 THz) |
| Near-infrared | f > ~10 THz | Interband transitions | d-band electrons dominate; Drude overestimates R |

BBRsim focuses on the Drude transition regime. Above ~10 THz, interband effects
begin to contribute; the Drude model overestimates R slightly, consistent with
comparisons against Palik room-temperature data.

---

## Validation

Three independent data sources are compared in `scripts/plot_cu_reflectance.py`:

- **Serov et al. (2016)** — direct reflectance measurements at 4 K for OF copper
  (150 GHz) and hydrogen-annealed copper (230 GHz). Primary cryogenic validation points.
- **Palik Handbook Vol. 1, Table 1** — room-temperature n, k optical constants;
  R computed from Fresnel. Valid above ~1 THz where k > 1 (anomalous 1–2 THz points excluded).
- **Geant4 IR reflectivity table** — 56-point digitized table from the Geant4
  optical physics package.

The `check_cu_absorptance.py` script compares simulated absorptance against a
Planck-weighted Drude integral at 4 K and flags deviations beyond a factor of 3:

```bash
./BBRSim test.mac >out.txt 2>&1
conda run -n bbrsim python scripts/check_cu_absorptance.py out.txt --rrr 100
```

---

## References

- Griffiths, D.J. (2017). *Introduction to Electrodynamics*, 4th ed. §9.4.
- Serov, Y.L. et al. (2016). *IEEE Trans. Microwave Theory Tech.* **64**(11), 3828.
- Palik, E.D. (1985). *Handbook of Optical Constants of Solids*, Vol. 1, Table 1.
- Chang, Y.-Y. (2023). *SuperCDMS HVeV Run 2 … Blackbody Radiation in Cryogenic
  Experiments*. PhD thesis, Caltech, Chapter 5.
- Golwala & Mirabolfathi, *BBRsim* NSF QIS proposal (2025).
