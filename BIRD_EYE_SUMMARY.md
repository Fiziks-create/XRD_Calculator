# CharacterizationCalc XRD — Bird's Eye Summary

> **What the package does, why each step matters, and how the outputs connect**

---

## 30-Second Picture

```
Raw diffractogram (.xlsx / .csv)
         │
         ▼
┌─────────────────────────────────────┐
│   1. Background removal             │
│   2. Peak detection                 │  ◄── scipy.signal.find_peaks
│   3. Pseudo-Voigt fitting → FWHM   │  ◄── Thompson et al. (1987)
└─────────────────────────────────────┘
         │
         ├──► Bragg's Law           → d-spacing (Å) per peak
         │
         ├──► Scherrer Equation     → mean crystallite size D̄ (nm)
         │
         ├──► Williamson–Hall (UDM) → D_WH (nm)  +  ε_WH (dimensionless)
         │                                           size & strain together
         ├──► Size–Strain Plot      → D_SSP (nm) +  ε_rms
         │
         ├──► Stokes–Wilson         → ε per peak
         │
         ├──► Dislocation density   → δ = 1/D² [lines m⁻²]
         │
         ├──► Crystallinity %       → area-ratio method (SNIP BG)
         │
         ├──► Phase comparison      → d-spacing match vs reference
         │
         └──► Texture coefficient   → Harris TC per peak
```

---

## Analysis Chain — Step by Step

### Step 1 — Data Ingestion
The tool reads a two-column file (2θ in degrees, intensity in arbitrary units).  
Columns are auto-detected by name so any reasonable header works.  
An optional second file (reference standard) enables phase comparison.

### Step 2 — Peak Detection
`scipy.signal.find_peaks` scans the pattern for local maxima.  
Three tunable thresholds (height %, prominence %, minimum distance in 2°) filter noise peaks from real Bragg reflections.  
Typically 3–20 peaks are identified for a polycrystalline sample.

### Step 3 — Pseudo-Voigt Peak Fitting
Each detected peak is individually fitted with a **Pseudo-Voigt** profile  
`η × Lorentzian + (1-η) × Gaussian`  
to extract the **FWHM (β)** in degrees, which is then converted to radians.  
This is the most critical step — all size and strain calculations depend on β.

> **Why Pseudo-Voigt?**  
> Real XRD peaks are neither purely Gaussian (thermal broadening limit) nor purely Lorentzian (strain broadening limit). The mixing parameter η captures both contributions simultaneously.

### Step 4 — Bragg's Law → d-spacings
```
d (Å) = λ / (2 sin θ)
```
Every peak position gives a d-spacing — the fingerprint for phase identification.  
Standard d-spacings from JCPDS/PDF-4 cards are compared when a reference is uploaded.

### Step 5 — Scherrer Equation → Crystallite Size
```
D (nm) = K λ / (β cos θ)
```
Applied **per peak** → distribution of sizes across reflections.  
Optional instrument-broadening correction via Warren's formula removes the diffractometer's intrinsic contribution.  
**Limitation:** Scherrer treats all broadening as size-only; it overestimates D when microstrain is present.

### Step 6 — Williamson–Hall Analysis → Size + Strain (simultaneously)
```
β cos θ = (Kλ/D) + 4ε sin θ
```
Plotting β·cosθ vs 4·sinθ gives a straight line:
- **y-intercept** → volume-weighted crystallite size D_WH (reciprocal)
- **slope** → upper-bound microstrain ε_WH

> **Key insight:** If the WH plot is truly linear with good R², size and strain are the dominant broadening mechanisms. Non-linearity suggests anisotropic strain or stacking faults.

### Step 7 — Size–Strain Plot (SSP) — Independent Cross-Check
```
(d · β · cos θ)² = (Kλ/D_SSP) · (d²·β·cosθ) + (ε_rms/0.9)²
```
An alternative linearisation (Mote et al. 2012) that gives:
- **Area-weighted** size D_SSP (less sensitive to large crystallite tails)
- **rms microstrain** ε_rms

> **Use together:** If D_WH ≈ D_SSP and ε_WH ≈ ε_rms, the analysis is internally consistent.  
> Large divergence suggests non-uniform strain or mixed phases.

### Step 8 — Stokes–Wilson Microstrain — Per-Peak View
```
ε = β / (4 tan θ)
```
A quick per-peak estimate. Scatter across peaks is normal and reflects anisotropic strain.  
The WH/SSP values are statistically more robust but this plot identifies outlier peaks.

### Step 9 — Dislocation Density
```
δ = 1 / D²   [lines m⁻²]
```
Connects crystallite size to defect density.  
Typical ranges:
- Annealed metals: 10¹⁰ – 10¹² m⁻²
- Cold-worked / ball-milled: 10¹⁴ – 10¹⁶ m⁻²
- Nanoparticles: 10¹⁵ – 10¹⁷ m⁻²

### Step 10 — Crystallinity Index
```
Xc = A_crystalline / (A_crystalline + A_amorphous)  × 100 %
```
Background is estimated by SNIP (Statistics-sensitive Non-linear Iterative Peak-clipping).  
The ratio of crystalline peak area to total pattern area gives the crystallinity fraction.  
**Use as relative comparison** between batches — absolute values depend on background model.

### Step 11 — Phase Comparison (if reference supplied)
Measured d-spacings are matched to reference d-values within a tolerance (default ±0.05 Å).  
Match % = (matched peaks) / (total measured peaks) × 100.  
A match > 70 % strongly suggests the correct phase; < 50 % may indicate impurities or incorrect reference.

### Step 12 — Texture Coefficient (Harris, 1952)
```
TC(hkl) = [I(hkl) / I₀(hkl)] / [ (1/N) Σ I(hkl)/I₀(hkl) ]
```
TC = 1 for a randomly oriented (powder) sample.  
TC > 1 means that plane is *preferentially oriented* (texture/preferred growth direction).  
TC < 1 means that plane is *underrepresented*.  
Requires a reference for I₀ values.

---

## Output Quick-Reference

| Output file | What to look at first |
|-------------|----------------------|
| `00_summary_dashboard.png` | 4-panel overview — start here |
| `01_raw_pattern.png` | Check that all peaks were detected (red dots) |
| `02_peak_fits.png` | Verify FWHM fits are sensible (not clipped / too wide) |
| `03_scherrer.png` | Size distribution and histogram |
| `04_williamson_hall.png` | R² value — is the fit linear? |
| `05_size_strain_plot.png` | Cross-check WH; compare D and ε |
| `06_microstrain.png` | Any peaks with anomalously high strain? |
| `07_dislocation.png` | Which reflections carry the most defects? |
| `08_reference.png` | Which peaks matched / were missed? |
| `XRD_Analysis.xlsx` | All numbers, 5 sheets, ready for your report |
| `Analysis_Report.txt` | Narrative summary with references |

---

## Common Interpretation Guide

| Observation | Likely meaning |
|-------------|----------------|
| **D_Scherrer ≪ D_WH** | Significant microstrain; Scherrer alone underestimates |
| **D_WH ≈ D_SSP** | Analysis self-consistent; broadening is size + strain dominated |
| **WH slope ≈ 0** | Negligible microstrain; size broadening only |
| **WH R² < 0.7** | Non-uniform strain, mixed phases, or too few peaks |
| **Crystallinity < 40 %** | Predominantly amorphous — consider amorphous halo fitting |
| **TC >> 1 for one plane** | Strong preferred orientation / epitaxial growth |
| **Phase match < 50 %** | Wrong reference, impurity phase, or non-standard stoichiometry |
| **ε > 1 × 10⁻²** | Heavy mechanical deformation or very small nanoparticles |

---

## Computational Stack

| Layer | Library | Role |
|-------|---------|------|
| UI | `gradio 4.x` | Browser-based interface, no frontend code |
| Numerics | `numpy`, `scipy` | Peak finding, curve fitting, regression |
| Visualisation | `matplotlib` | 9 publication-quality plots |
| Data I/O | `pandas`, `openpyxl` | Excel/CSV read-write |
| Packaging | `zipfile`, `json` | ZIP output assembly |

---

## Module Roadmap

| Module | Status |
|--------|--------|
| **XRD** | ✅ Complete |
| FTIR | 🔜 Next |
| Raman | 🔜 Planned |
| UV-vis | 🔜 Planned |

---

*CharacterizationCalc · XRD Module v1.0*
