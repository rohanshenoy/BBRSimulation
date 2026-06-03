# Copper Reflectance Model in BBRsim

BBRsim simulates free-space blackbody radiation from 50 GHz to 20 THz. In this band, copper reflectance is controlled by the Drude model of free-electron conduction, not by room-temperature optical constants or simplified empirical fits. This document describes the physics, the parameterization, and the BBRsim implementation.

---

## Why Not Hagen-Rubens?

The Hagen-Rubens formula is widely used as a first estimate of metal reflectance:

```
D = 1 − R ≈ 2√(2ε₀ω/σ)
```

It is the **low-frequency limit** of the Drude model, valid only when ωτ << 1, i.e. when the photon frequency is well below the electron scattering rate 1/τ.

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

BBRsim's simulation range starts at 50 GHz. Hagen-Rubens is wrong for essentially the entire range. The full Drude model must be used.

---

## The Drude–Griffiths Model

From Griffiths, *Introduction to Electrodynamics* §9.4, for a plane wave in a conductor with conductivity σ(ω):

**Complex wave vector:**

```
k = ω√(ε₀μ₀/2) · √[ √(1 + (σ_r/ε₀ω)²) + 1 ]
κ = ω√(ε₀μ₀/2) · √[ √(1 + (σ_r/ε₀ω)²) − 1 ]
```

where σ_r = Re[σ(ω)] is the real part of the AC conductivity.

**Complex refractive index:**

```
ñ = n + ik = c(k + iκ)/ω
```

**Normal-incidence reflectance:**

```
R = |(ñ − 1)/(ñ + 1)|²  =  [(n−1)² + k²] / [(n+1)² + k²]
```

**Full Drude AC conductivity** (replaces static σ in Griffiths's equations):

```
σ(ω) = σ_DC / (1 − iωτ)
```

At low frequency (ωτ << 1): σ(ω) → σ_DC and the Hagen-Rubens limit is recovered. At high frequency (ωτ >> 1): σ(ω) → σ_DC/(−iωτ) → pure imaginary, reducing reflectance below H-R.

---

## Temperature and RRR Parameterization

### RRR definition

The Residual Resistance Ratio is the single number characterizing copper quality:

```
RRR = ρ(273 K) / ρ(4 K)  =  σ(4 K) / σ(273 K)
```

At room temperature, resistivity is dominated by phonon scattering — the same for all copper grades. At cryogenic temperatures, phonons freeze out and resistivity is set by impurity/defect scattering, which varies with purity. RRR captures this ratio.

### Conductivity vs temperature: Matthiessen's rule

```
1/σ(T) = 1/σ_impurity + 1/σ_phonon(T)
```

where:
- `1/σ_impurity = 1/(RRR × σ_RT)` — fixed by defects, independent of T
- `1/σ_phonon(T) ≈ T / (273 × σ_RT)` — linear in T above ~50 K, freezes out below

At T = 4 K: σ_phonon >> σ_impurity, so σ(4 K) ≈ RRR × σ_RT. At T = 273 K: both terms contribute and the result is σ_RT regardless of RRR, consistent with the RRR definition.

---

## Copper Material Constants

| Constant | Symbol | Value | Units |
|---|---|---|---|
| Room-temp conductivity | σ_RT | 5.96×10⁷ | S/m |
| Free electron density | nₑ | 8.49×10²⁸ | m⁻³ |
| Electron mass | mₑ | 9.109×10⁻³¹ | kg |
| Electron charge | e | 1.602×10⁻¹⁹ | C |
| Debye temperature | Θ_D | 343 | K |

σ_RT is universal for all copper grades. It does not depend on purity.

---

## Typical RRR Values

| Grade | RRR | σ(4K) [S/m] | τ(4K) [ps] | H-R valid below |
|---|---|---|---|---|
| Commercial Cu | ~10 | 5.96×10⁸ | 0.25 ps | ~640 GHz |
| OFHC Cu | ~100 | 5.96×10⁹ | 2.5 ps | ~64 GHz |
| HP (H₂-annealed) Cu | ~57 | 3.4×10⁹ | 1.4 ps | ~115 GHz |
| Ultra-pure crystal | ~500 | 3.0×10¹⁰ | 12.5 ps | ~13 GHz |

For all grades, the Hagen-Rubens validity cutoff is inside or below the BBRsim range (50 GHz–20 THz). **The full Drude model is always required.**

---

## BBRsim API

The current implementation (`BBRMaterials.hh`) uses Hagen-Rubens with hardcoded conductivities. The planned upgrade replaces this with the full Drude model:

```cpp
// Planned API — replaces GetOFHCCopper(), GetOFCopperSerov(), GetHPCopperSerov():
G4Material* mat = BBRMaterials::GetCopper(/*RRR=*/100, /*T_K=*/4.0);

// Messenger-compatible name lookup (maps grade name to (RRR, T)):
G4Material* mat = BBRMaterials::GetCopperByName("OFHC_Cu");   // RRR=100, 4K
G4Material* mat = BBRMaterials::GetCopperByName("OF_Cu");     // RRR=57,  4K
G4Material* mat = BBRMaterials::GetCopperByName("HP_Cu");     // RRR=57,  4K (H2-annealed)
```

The function generates a 20-point log-spaced REFLECTIVITY table from 50 GHz to 20 THz, stored on the material's `G4MaterialPropertiesTable` as the `REFLECTIVITY` property read by `BBSimOpBoundaryProcess::HandleReflectanceBoundary`.

---

## Frequency Regimes

| Regime | Condition | Model | Notes |
|---|---|---|---|
| Hagen-Rubens | f << 1/(2πτ) | D ≈ 2√(2ε₀ω/σ) | Valid only below ~64 GHz for OFHC Cu at 4K |
| Drude transition | f ~ 1/(2πτ) | Full Griffiths §9.4 | BBRsim range (50 GHz–20 THz) |
| Near-infrared | f > ~10 THz | Interband transitions dominate | Surface finish and passivation matter; Drude becomes less accurate |

BBRsim focuses on the Drude transition regime. Above ~10 THz, interband effects (d-band electrons in Cu) begin to contribute; the Drude model overestimates R slightly in this range, consistent with what is seen when comparing against Palik room-temperature data.

---

## Validation Data

Three independent data sources are plotted in `scripts/plot_cu_reflectance.py`:

- **Serov et al. (2016)** — direct reflectance measurements at 4 K for OF_Cu (150 GHz) and HP_Cu (230 GHz). These are the primary cryogenic validation points.
- **Palik Handbook Vol. 1, Table 1** — room-temperature n, k optical constants; R computed from Fresnel. Valid above ~1 THz where k > 1 (anomalous data at 1–2 THz excluded).
- **Geant4 IR reflectivity table** (`Geant4_copper_IR_reflectivity.ods`) — digitized 56-point table used by the Geant4 optical physics package.

---

## References

- Griffiths, D.J. (2017). *Introduction to Electrodynamics*, 4th ed. §9.4 — Electromagnetic waves in conductors; complex wave vector; reflectance.
- Serov, Y.L. et al. (2016). Reflectivity of technical copper grades at cryogenic temperatures. *Cryogenics*.
- Palik, E.D. (1985). *Handbook of Optical Constants of Solids*, Vol. 1, Table 1 — n, k for copper at room temperature.
- Chang, Y.-Y. (2023). *SuperCDMS HVeV Run 2 … and the Blackbody Radiation in Cryogenic Experiments*. PhD thesis, Caltech, Chapter 5.
- Golwala & Mirabolfathi, *BBRsim* NSF QIS proposal (2025).
