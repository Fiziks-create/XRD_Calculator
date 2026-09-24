# CharacterizationCalc — XRD Analysis Module

> **Comprehensive X-ray Diffraction data analysis web application**  
> Built with Python · Gradio · Deployable on Render.com or Hugging Face Spaces

---

## Table of Contents
1. [Overview](#overview)
2. [Scientific Analyses](#scientific-analyses)
3. [Input File Format](#input-file-format)
4. [Quick Start (Local)](#quick-start-local)
5. [Deploy on Render.com](#deploy-on-rendercom)
6. [Deploy on Hugging Face Spaces](#deploy-on-hugging-face-spaces)
7. [Output Files](#output-files)
8. [Scientific Methods & References](#scientific-methods--references)
9. [Troubleshooting](#troubleshooting)

---

## Overview

CharacterizationCalc XRD is a **browser-based analysis tool** that processes raw X-ray diffraction data and returns a full suite of structural characterisation metrics — no coding required.

| Feature | Detail |
|---------|--------|
| Input formats | `.xlsx`, `.xls`, `.csv` |
| X-ray sources | Cu Kα · Co Kα · Mo Kα · Cr Kα · Fe Kα · Ag Kα · Custom |
| Analyses | 12 independent methods |
| Output | ZIP with 9 plots · Excel · CSV · JSON · text report |
| Deployment | Render.com (free) · Hugging Face Spaces · Local |

---

## Scientific Analyses

| # | Method | What you get |
|---|--------|-------------|
| 1 | **Peak Detection** | Position (2θ), height, prominence for all peaks |
| 2 | **Pseudo-Voigt Fitting** | FWHM (β) per peak; η mixing parameter |
| 3 | **Integral Breadth** | β = Area / Height alternative broadening measure |
| 4 | **Bragg's Law** | d-spacing (Å) → interplanar spacing |
| 5 | **Scherrer Equation** | Mean crystallite size D (nm) ± instrument correction |
| 6 | **Williamson–Hall (UDM)** | Volume-weighted D + microstrain ε (simultaneous) |
| 7 | **Size–Strain Plot (SSP)** | Area-weighted D + rms strain; Mote et al. (2012) |
| 8 | **Microstrain (Stokes–Wilson)** | Per-peak ε = β / (4 tan θ) |
| 9 | **Dislocation Density** | δ = 1/D² [lines m⁻²]; Williamson & Smallman (1956) |
| 10 | **Crystallinity Index** | % crystallinity via SNIP background subtraction |
| 11 | **Phase Comparison** | d-spacing match against optional reference standard |
| 12 | **Texture Coefficient** | Harris TC per peak; preferred orientation detection |

---

## Input File Format

### Sample Data File (required)

Your file must contain **exactly two numeric columns** — the tool auto-detects them by name.

**Accepted column name variants:**

| Column | Accepted headers (case-insensitive) |
|--------|-------------------------------------|
| 2θ angle | `2theta`, `2θ`, `2-theta`, `twotheta`, `angle`, `2th`, `degrees` |
| Intensity | `intensity`, `counts`, `signal`, `i`, `int`, `y`, `cps` |

**Example CSV:**
```
2theta,intensity
10.00,120
10.05,118
10.10,135
...
45.23,8500
45.28,22300
45.33,9100
```

**Example Excel:**  
Column A = `2θ (°)`, Column B = `Intensity (a.u.)` — both numeric, no merged cells.

### Reference File (optional)

A second file with the standard d-spacings for phase identification:
```
2theta,relative_intensity
31.77,100
45.48,42
56.60,20
...
```
**Source:** Download standard JCPDS/PDF cards, or generate from COD/ICSD data.

---

## Quick Start (Local)

### Prerequisites
- Python 3.9 – 3.12
- `pip`

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/YOUR_USERNAME/characterization-calc-xrd.git
cd characterization-calc-xrd

# 2. Create virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run the app
python app.py
```

Open your browser at **http://localhost:7860**

---

## Deploy on Render.com

This repository ships with `render.yaml` for one-click Blueprint deployment.

### Steps

1. **Fork** this repository to your GitHub account.

2. Go to [render.com](https://render.com) → **New** → **Blueprint**.

3. Connect your GitHub account and select the forked repo.

4. Render reads `render.yaml` automatically — click **Apply**.

5. Wait ~2 minutes for the build to complete.

6. Your app is live at `https://characterization-calc-xrd.onrender.com`  
   *(free tier spins down after inactivity; first request takes ~30 s to wake)*

### Manual Render Setup (alternative)

If you prefer manual configuration:

| Setting | Value |
|---------|-------|
| **Runtime** | Python 3 |
| **Build Command** | `pip install -r requirements.txt` |
| **Start Command** | `python app.py` |
| **Environment Variable** | `PORT = 7860` |
| **Plan** | Free |

---

## Deploy on Hugging Face Spaces

1. Create a new Space at [huggingface.co/spaces](https://huggingface.co/spaces).
2. Choose **Gradio** as the SDK.
3. Upload `app.py`, `backend.py`, `requirements.txt`.
4. The Space builds automatically — no `render.yaml` needed.

> **Note:** Rename `app.py` → keep as `app.py` (Spaces expects this).

---

## Output Files

After clicking **Run Analysis**, a ZIP file is generated with:

```
XRD_<MaterialName>_<Timestamp>.zip
├── plots/
│   ├── 00_summary_dashboard.png   ← 4-panel overview (publication-ready)
│   ├── 01_raw_pattern.png         ← Diffractogram + detected peaks
│   ├── 02_peak_fits.png           ← Individual Pseudo-Voigt fits
│   ├── 03_scherrer.png            ← Crystallite size vs 2θ + histogram
│   ├── 04_williamson_hall.png     ← WH (UDM) linear regression plot
│   ├── 05_size_strain_plot.png    ← SSP linear regression plot
│   ├── 06_microstrain.png         ← Per-peak microstrain bar chart
│   ├── 07_dislocation.png         ← Dislocation density bar chart
│   └── 08_reference.png           ← Phase comparison (if ref uploaded)
├── results/
│   ├── XRD_Analysis.xlsx          ← Full data (5 Excel sheets)
│   ├── peak_data.csv              ← Per-peak data (CSV)
│   ├── Analysis_Report.txt        ← Complete text report + references
│   └── analysis_data.json         ← Raw numerical output (JSON)
└── README.txt                     ← Quick-reference inside ZIP
```

### Excel Sheets

| Sheet | Content |
|-------|---------|
| `Peak Analysis` | 2θ, d, FWHM, Scherrer D, microstrain, dislocation density per peak |
| `Summary` | Averaged/global quantities (WH size, SSP size, crystallinity, …) |
| `Williamson-Hall` | x = 4sinθ, y = β·cosθ raw data for custom plotting |
| `Size-Strain Plot` | x = d²·β·cosθ, y = (d·β·cosθ)² raw data |
| `Reference Match` | Measured d vs reference d, deviation, match status |

---

## Scientific Methods & References

### Bragg's Law [R1]
```
d = λ / (2 sin θ)
```
Converts 2θ peak positions to interplanar d-spacings.

### Scherrer Equation [R2, R3]
```
D = K λ / (β cos θ)
```
- **K** = shape factor (default 0.9 for near-spherical crystallites)
- **β** = FWHM in radians (instrument-broadening corrected if σ_instr > 0)
- Instrument correction: β²_sample = β²_measured − β²_instrument  [R7]

### Williamson–Hall (UDM) [R5]
```
β cos θ = (K λ / D) + 4 ε sin θ
```
Linear fit of β·cosθ vs 4·sinθ:  
- **Intercept** → volume-weighted crystallite size D_WH  
- **Slope** → upper-bound microstrain ε_WH

### Size–Strain Plot (SSP) [R10]
```
(d · β · cos θ)² = (K λ / D²) · (d² · β · cos θ) + (ε/0.9)²
```
- **Slope** → area-weighted crystallite size D_SSP  
- **Intercept** → rms microstrain ε_rms  
- Reference: Mote, Purushotham & Dole (2012) *J. Theor. Appl. Phys.* 6, 6

### Microstrain — Stokes & Wilson [R4]
```
ε = β / (4 tan θ)
```
Per-peak estimate; scatter between peaks is normal.

### Dislocation Density [R6]
```
δ = 1 / D²   [lines m⁻²]
```
Higher crystallinity / larger D → lower dislocation density.

### Crystallinity Index
```
Xc = A_crystalline / (A_crystalline + A_amorphous)
```
Background estimated by SNIP algorithm; use as relative comparison.

### Texture Coefficient — Harris [R8]
```
TC(hkl) = [I(hkl) / I₀(hkl)] / [(1/N) Σ I(hkl) / I₀(hkl)]
```
TC > 1 indicates preferred orientation in that plane.

### Full Reference List

| Code | Citation |
|------|----------|
| [R1] | Bragg, W. L. (1913). *Proc. R. Soc. A* 88, 428–438 |
| [R2] | Scherrer, P. (1918). *Nachr. Ges. Wiss. Göttingen*, 98–100 |
| [R3] | Patterson, A. L. (1939). *Phys. Rev.* 56, 978–982 |
| [R4] | Stokes, A. R. & Wilson, A. J. C. (1944). *Proc. Phys. Soc.* 56, 174 |
| [R5] | Williamson, G. K. & Hall, W. H. (1953). *Acta Metall.* 1, 22–31 |
| [R6] | Williamson, G. K. & Smallman, R. E. (1956). *Phil. Mag.* 1, 34–46 |
| [R7] | Warren, B. E. (1969). *X-ray Diffraction*. Addison-Wesley |
| [R8] | Harris, G. B. (1952). *Phil. Mag.* 43, 113–123 |
| [R9] | Thompson, P. et al. (1987). *J. Appl. Crystallogr.* 20, 79–83 |
| [R10] | Mote, V. D. et al. (2012). *J. Theor. Appl. Phys.* 6, 6 |
| [R11] | Cullity, B. D. & Stock, S. R. (2001). *Elements of X-ray Diffraction*, 3rd Ed. |
| [R12] | Venkateswarlu, K. et al. (2010). *Trans. Nonferrous Met. Soc. China* 20, 941 |

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| **"Could not auto-detect 2θ column"** | Rename your columns to `2theta` and `intensity` |
| **"No peaks detected"** | Lower *Min Height %* or *Min Prominence %* in the UI |
| **Williamson–Hall shows negative D** | Fewer than 3 peaks or extremely narrow peaks — increase data range |
| **All sizes identical** | Data may be a single sharp peak — need ≥ 3 peaks for WH/SSP |
| **App not waking on Render (free)** | Free tier sleeps after 15 min; first request wakes it in ~30 s |
| **Excel file opens blank** | Ensure `openpyxl` is installed: `pip install openpyxl` |
| **Crystallinity is 0%** | Background estimation failed; try a wider 2θ range |

---

## Project Structure

```
characterization-calc-xrd/
├── app.py              ← Gradio web interface
├── backend.py          ← Scientific analysis engine (12 methods)
├── requirements.txt    ← Python dependencies
├── render.yaml         ← Render.com Blueprint deployment config
├── .gitignore          ← Git exclusions
└── README.md           ← This file
```

---

## Licence

MIT — free to use, modify, and distribute with attribution.

---

*CharacterizationCalc XRD · Part of the CharacterizationCalc suite*  
*Other modules: FTIR · Raman · UV-vis (coming soon)*
