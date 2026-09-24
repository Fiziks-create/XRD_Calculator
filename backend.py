"""
================================================================================
  XRD Analysis Backend — CharacterizationCalc
================================================================================
Comprehensive X-ray Diffraction data analysis engine.

Scientific Methods:
  1.  Peak Detection           — scipy.signal.find_peaks (height/prominence/distance)
  2.  Pseudo-Voigt Fitting     — Thompson et al. (1987)
  3.  Integral Breadth         — alternative broadening measure
  4.  Bragg's Law              — d-spacing from 2θ
  5.  Scherrer Equation        — crystallite size (with instrument-broadening correction)
  6.  Williamson-Hall (UDM)    — size + microstrain separation
  7.  Size-Strain Plot (SSP)   — alternative WH analysis (Mote et al., 2012)
  8.  Microstrain (Stokes-Wilson) — per-peak strain estimate
  9.  Dislocation Density      — Williamson & Smallman (1956)
  10. Crystallinity Estimate   — area-ratio method
  11. Phase Comparison         — d-spacing match with optional reference
  12. Texture Coefficient      — Harris (1952)

Key References:
  [R1]  Bragg (1913) Proc. R. Soc. A 88, 428–438
  [R2]  Scherrer (1918) Nachr. Ges. Wiss. Göttingen, 98–100
  [R3]  Patterson (1939) Phys. Rev. 56, 978–982
  [R4]  Stokes & Wilson (1944) Proc. Phys. Soc. 56, 174
  [R5]  Williamson & Hall (1953) Acta Metall. 1, 22–31
  [R6]  Williamson & Smallman (1956) Phil. Mag. 1, 34–46
  [R7]  Warren (1969) X-ray Diffraction, Addison-Wesley
  [R8]  Harris (1952) Phil. Mag. 43, 113–123
  [R9]  Thompson et al. (1987) J. Appl. Crystallogr. 20, 79–83
  [R10] Mote et al. (2012) J. Theor. Appl. Phys. 6, 6
  [R11] Cullity & Stock (2001) Elements of X-ray Diffraction, 3rd Ed.
  [R12] Venkateswarlu et al. (2010) Trans. Nonferrous Met. Soc. China 20, 941
================================================================================
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
import matplotlib.ticker as ticker

from scipy.signal import find_peaks, peak_widths, savgol_filter
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.ndimage import uniform_filter1d

import warnings
import io
import os
import json
import zipfile
from datetime import datetime

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
WAVELENGTHS = {
    "Cu Kα (λ = 1.5406 Å)": 1.54060,
    "Co Kα (λ = 1.7902 Å)": 1.79020,
    "Mo Kα (λ = 0.7107 Å)": 0.71073,
    "Cr Kα (λ = 2.2909 Å)": 2.29090,
    "Fe Kα (λ = 1.9373 Å)": 1.93730,
    "Ag Kα (λ = 0.5594 Å)": 0.55940,
}

SCHERRER_K_DEFAULT = 0.9          # Scherrer shape factor for spherical crystallites
INSTRUMENT_FWHM_DEFAULT = 0.0     # Default instrument broadening (°), 0 = no correction


# ─────────────────────────────────────────────────────────────────────────────
# Utility: Column Auto-Detection
# ─────────────────────────────────────────────────────────────────────────────
_TWO_THETA_ALIASES = [
    '2theta', '2θ', '2-theta', 'twotheta', 'two_theta', 'angle', '2th',
    '2t', 'degrees', 'deg', 'position', 'pos', '2θ(°)', '2θ (°)'
]
_INTENSITY_ALIASES = [
    'intensity', 'counts', 'count', 'int', 'i(', 'cts', 'signal',
    'y', 'i_net', 'i_', 'peak', 'yobs', 'cps', 'intensity (a.u.)',
    'intensity(a.u.)', 'int (a.u.)', 'a.u.'
]

def _detect_columns(df: pd.DataFrame):
    """
    Auto-detect the 2θ and Intensity columns from a DataFrame.
    Returns (theta_col, intensity_col).
    Falls back to columns [0, 1] if detection fails.
    """
    theta_col = intensity_col = None
    for col in df.columns:
        low = col.lower().strip()
        if any(alias in low for alias in _TWO_THETA_ALIASES):
            theta_col = col
        if any(alias in low for alias in _INTENSITY_ALIASES):
            intensity_col = col
    # fallback
    if theta_col is None:
        theta_col = df.columns[0]
    if intensity_col is None:
        intensity_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
    return theta_col, intensity_col


def load_xrd_data(file_path: str):
    """
    Load XRD data from Excel (.xlsx / .xls) or CSV (.csv).
    Returns (two_theta, intensity, col_names_dict, preview_html).
    Raises ValueError with user-friendly messages on failure.
    """
    ext = os.path.splitext(file_path)[-1].lower()
    try:
        if ext in ('.xlsx', '.xls'):
            df = pd.read_excel(file_path)
        elif ext == '.csv':
            # try common separators
            for sep in (',', '\t', ';', ' '):
                try:
                    df_try = pd.read_csv(file_path, sep=sep)
                    if df_try.shape[1] >= 2:
                        df = df_try
                        break
                except Exception:
                    continue
            else:
                df = pd.read_csv(file_path)
        else:
            raise ValueError(f"Unsupported file type: '{ext}'. Upload .xlsx, .xls, or .csv.")
    except Exception as exc:
        raise ValueError(f"Could not read file: {exc}")

    if df.shape[1] < 2:
        raise ValueError("File must have at least 2 columns (2θ and Intensity).")

    # Strip whitespace from column names
    df.columns = [str(c).strip() for c in df.columns]
    theta_col, int_col = _detect_columns(df)

    # Drop NaN rows, coerce to float
    sub = df[[theta_col, int_col]].dropna()
    try:
        two_theta = sub[theta_col].astype(float).values
        intensity = sub[int_col].astype(float).values
    except Exception:
        raise ValueError("Non-numeric values found in data columns. Check your file.")

    if len(two_theta) < 10:
        raise ValueError("Fewer than 10 valid data points found — check your file.")

    # Ensure sorted ascending by 2θ
    sort_idx = np.argsort(two_theta)
    two_theta = two_theta[sort_idx]
    intensity = intensity[sort_idx]

    col_info = {"theta_col": theta_col, "intensity_col": int_col,
                "n_rows": len(two_theta), "columns": list(df.columns)}
    return two_theta, intensity, col_info


def load_reference(file_path: str):
    """
    Load reference/standard XRD pattern (CSV or Excel).
    Returns DataFrame with first two columns as [2θ, RelIntensity].
    """
    ext = os.path.splitext(file_path)[-1].lower()
    try:
        if ext in ('.xlsx', '.xls'):
            ref = pd.read_excel(file_path)
        else:
            ref = pd.read_csv(file_path)
        ref.columns = [str(c).strip() for c in ref.columns]
        ref = ref.dropna().iloc[:, :2]
        ref.columns = ['two_theta_ref', 'intensity_ref']
        ref = ref.astype(float)
        ref = ref.sort_values('two_theta_ref').reset_index(drop=True)
        return ref
    except Exception as exc:
        raise ValueError(f"Could not read reference file: {exc}")


# ─────────────────────────────────────────────────────────────────────────────
# Preprocessing
# ─────────────────────────────────────────────────────────────────────────────
def smooth_signal(intensity: np.ndarray, window: int = 11, poly: int = 3) -> np.ndarray:
    """Savitzky-Golay smoothing. Falls back to uniform filter if array too short."""
    if len(intensity) < window:
        return uniform_filter1d(intensity, size=max(3, len(intensity) // 5))
    try:
        return savgol_filter(intensity, window_length=window, polyorder=poly)
    except Exception:
        return uniform_filter1d(intensity, size=window)


def estimate_background(two_theta: np.ndarray, intensity: np.ndarray,
                         n_iter: int = 40) -> np.ndarray:
    """
    SNIP (Statistics-sensitive Non-linear Iterative Peak-clipping) background.
    Widely used for XRD background subtraction.

    Ryan & Cousins (1988) Nuclear Instruments and Methods.
    """
    bg = intensity.copy().astype(float)
    n = len(bg)
    for p in range(1, n_iter + 1):
        bg_new = bg.copy()
        for i in range(p, n - p):
            v = (bg[i - p] + bg[i + p]) / 2.0
            bg_new[i] = min(bg[i], v)
        bg = bg_new
    return bg


# ─────────────────────────────────────────────────────────────────────────────
# Peak Detection
# ─────────────────────────────────────────────────────────────────────────────
def detect_peaks(two_theta: np.ndarray, intensity: np.ndarray,
                 min_height_frac: float = 0.05,
                 min_prominence_frac: float = 0.03,
                 min_distance_deg: float = 0.5) -> tuple:
    """
    Detect XRD peaks using scipy.signal.find_peaks with scaled thresholds.

    Parameters
    ----------
    two_theta         : 2θ array (°)
    intensity         : Intensity array
    min_height_frac   : Minimum height as fraction of max intensity
    min_prominence_frac : Minimum prominence as fraction of max intensity
    min_distance_deg  : Minimum separation between peaks (°)

    Returns
    -------
    peaks : indices of detected peaks
    props : dict with scipy peak properties
    """
    step = float(np.median(np.diff(two_theta)))
    min_dist_pts = max(1, int(min_distance_deg / step))
    max_int = np.max(intensity)

    peaks, props = find_peaks(
        intensity,
        height=min_height_frac * max_int,
        prominence=min_prominence_frac * max_int,
        distance=min_dist_pts,
        width=1
    )
    return peaks, props


# ─────────────────────────────────────────────────────────────────────────────
# Peak Fitting (Pseudo-Voigt)
# ─────────────────────────────────────────────────────────────────────────────
def _pseudo_voigt(x, amplitude, center, fwhm, eta, baseline):
    """
    Pseudo-Voigt profile (Thompson, Cox & Hastings, 1987, J. Appl. Crystallogr.).
    Both Gaussian and Lorentzian components share the same FWHM parameter.

    PV(x) = A * [η·L(x) + (1−η)·G(x)] + baseline
    G(x)  = exp(−(x−c)²/(2σ_G²)),  σ_G = fwhm/(2√(2ln2))
    L(x)  = σ_L²/((x−c)²+σ_L²),   σ_L = fwhm/2
    η ∈ [0,1]: 0 = pure Gaussian, 1 = pure Lorentzian
    """
    sigma_g = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    sigma_l = fwhm / 2.0
    gauss = np.exp(-((x - center) ** 2) / (2.0 * sigma_g ** 2))
    lorentz = (sigma_l ** 2) / ((x - center) ** 2 + sigma_l ** 2)
    return amplitude * (eta * lorentz + (1.0 - eta) * gauss) + baseline


def _gaussian(x, amplitude, center, fwhm, baseline):
    """Pure Gaussian profile (fallback)."""
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return amplitude * np.exp(-((x - center) ** 2) / (2.0 * sigma ** 2)) + baseline


def _r2(y_obs, y_fit):
    ss_res = np.sum((y_obs - y_fit) ** 2)
    ss_tot = np.sum((y_obs - np.mean(y_obs)) ** 2)
    return 0.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot


def fit_single_peak(two_theta: np.ndarray, intensity: np.ndarray,
                    peak_idx: int, window_deg: float = 3.5) -> dict:
    """
    Fit a single XRD peak with a Pseudo-Voigt profile.
    Falls back to Gaussian if Pseudo-Voigt fails.

    Returns dict with keys:
      type, center, fwhm_deg, fwhm_rad, amplitude, eta, baseline,
      r_squared, x_fit, y_fit, x_region, y_region, integral_breadth_deg
    """
    step = float(np.median(np.diff(two_theta)))
    half_win = max(40, int(window_deg / step))
    start = max(0, peak_idx - half_win)
    end   = min(len(two_theta), peak_idx + half_win + 1)

    xr = two_theta[start:end]
    yr = intensity[start:end].copy()

    peak_max   = float(intensity[peak_idx])
    baseline_0 = float(np.percentile(yr, 10))
    amplitude_0 = peak_max - baseline_0
    center_0   = float(two_theta[peak_idx])

    # Initial FWHM from scipy peak widths
    try:
        w, _, _, _ = peak_widths(intensity, [peak_idx], rel_height=0.5)
        fwhm_0 = float(w[0]) * step
    except Exception:
        fwhm_0 = 0.3

    fwhm_0 = max(fwhm_0, 2 * step)

    result = {"type": "failed", "fwhm_deg": fwhm_0, "fwhm_rad": np.radians(fwhm_0),
              "r_squared": 0.0, "center": center_0, "amplitude": amplitude_0,
              "eta": 0.5, "baseline": baseline_0,
              "x_region": xr, "y_region": yr,
              "x_fit": xr, "y_fit": yr,
              "integral_breadth_deg": fwhm_0}

    # ── Try Pseudo-Voigt ──────────────────────────────────────────────────────
    try:
        p0 = [amplitude_0, center_0, fwhm_0, 0.5, baseline_0]
        lo = [0, center_0 - 1.5, step, 0.0, 0.0]
        hi = [amplitude_0 * 3, center_0 + 1.5, 10.0, 1.0, amplitude_0 * 0.5]
        popt, pcov = curve_fit(_pseudo_voigt, xr, yr, p0=p0,
                               bounds=(lo, hi), maxfev=8000)
        amp, cen, fwhm, eta, base = popt
        xf = np.linspace(xr.min(), xr.max(), 600)
        yf = _pseudo_voigt(xf, *popt)
        r2 = _r2(yr, _pseudo_voigt(xr, *popt))

        # Integral breadth β_I = Area / Height
        area = float(np.trapezoid(np.maximum(yr - base, 0), xr))
        ib   = area / max(amp, 1e-12)

        if fwhm > 0 and r2 > 0.5:
            result.update({"type": "pseudo_voigt",
                           "center": cen, "fwhm_deg": fwhm,
                           "fwhm_rad": np.radians(fwhm),
                           "amplitude": amp, "eta": eta, "baseline": base,
                           "r_squared": r2,
                           "x_fit": xf, "y_fit": yf,
                           "integral_breadth_deg": ib})
            return result
    except Exception:
        pass

    # ── Fallback: Gaussian ────────────────────────────────────────────────────
    try:
        p0g = [amplitude_0, center_0, fwhm_0, baseline_0]
        lo  = [0, center_0 - 1.5, step, 0.0]
        hi  = [amplitude_0 * 3, center_0 + 1.5, 10.0, amplitude_0 * 0.5]
        popt_g, _ = curve_fit(_gaussian, xr, yr, p0=p0g,
                              bounds=(lo, hi), maxfev=8000)
        amp_g, cen_g, fwhm_g, base_g = popt_g
        xf = np.linspace(xr.min(), xr.max(), 600)
        yf = _gaussian(xf, *popt_g)
        r2_g = _r2(yr, _gaussian(xr, *popt_g))
        area = float(np.trapezoid(np.maximum(yr - base_g, 0), xr))
        ib   = area / max(amp_g, 1e-12)
        if fwhm_g > 0 and r2_g > 0.3:
            result.update({"type": "gaussian",
                           "center": cen_g, "fwhm_deg": fwhm_g,
                           "fwhm_rad": np.radians(fwhm_g),
                           "amplitude": amp_g, "eta": 0.0, "baseline": base_g,
                           "r_squared": r2_g,
                           "x_fit": xf, "y_fit": yf,
                           "integral_breadth_deg": ib})
    except Exception:
        pass

    return result


def correct_instrument_broadening(fwhm_meas_rad: np.ndarray,
                                  beta_instr_deg: float = 0.0) -> np.ndarray:
    """
    Warren correction for instrument broadening (Gaussian convolution model):
        β_sample² = β_measured² − β_instrument²

    References: Warren (1969), Cullity & Stock (2001) p.96
    """
    if beta_instr_deg <= 0.0:
        return fwhm_meas_rad
    b_i = np.radians(beta_instr_deg)
    corrected = np.sqrt(np.maximum(fwhm_meas_rad ** 2 - b_i ** 2, (0.01 * fwhm_meas_rad) ** 2))
    return corrected


# ─────────────────────────────────────────────────────────────────────────────
# Crystallographic Analysis
# ─────────────────────────────────────────────────────────────────────────────
def bragg_law(two_theta_deg: np.ndarray, wavelength: float) -> np.ndarray:
    """
    Apply Bragg's Law to obtain d-spacing.
        nλ = 2d·sin θ  →  d = λ / (2·sin θ)

    Reference: [R1] Bragg (1913)

    Parameters
    ----------
    two_theta_deg : 2θ in degrees
    wavelength    : X-ray wavelength in Å

    Returns
    -------
    d_spacing : d-spacings in Å
    """
    theta_rad = np.radians(np.asarray(two_theta_deg, dtype=float) / 2.0)
    return wavelength / (2.0 * np.sin(theta_rad))


def scherrer_size(fwhm_rad: np.ndarray, two_theta_deg: np.ndarray,
                  wavelength: float, K: float = 0.9) -> np.ndarray:
    """
    Scherrer equation for volume-averaged crystallite size.
        D = Kλ / (β·cos θ)

    References: [R2] Scherrer (1918), [R3] Patterson (1939)

    Parameters
    ----------
    fwhm_rad      : FWHM values in radians (already instrument-corrected)
    two_theta_deg : 2θ positions in degrees
    wavelength    : X-ray wavelength in Å
    K             : shape factor (0.9 default; 0.94 for cubic, 1.0 for sphere)

    Returns
    -------
    D_nm : crystallite sizes in nm
    """
    theta_rad = np.radians(np.asarray(two_theta_deg, dtype=float) / 2.0)
    beta = np.asarray(fwhm_rad, dtype=float)
    # Guard against zero or negative
    beta = np.maximum(beta, 1e-10)
    D_angstrom = (K * wavelength) / (beta * np.cos(theta_rad))
    return D_angstrom / 10.0  # Å → nm


# ─────────────────────────────────────────────────────────────────────────────
# Williamson-Hall Analysis (UDM)
# ─────────────────────────────────────────────────────────────────────────────
def williamson_hall(two_theta_deg: np.ndarray, fwhm_rad: np.ndarray,
                    wavelength: float, K: float = 0.9) -> dict:
    """
    Uniform Deformation Model (UDM) Williamson-Hall analysis.
        β·cos θ = (Kλ/D) + 4ε·sin θ

    Plot: Y = β·cos θ  vs  X = 4·sin θ
    → Slope     = ε   (microstrain, signed: + tensile, − compressive)
    → Intercept = Kλ/D → D = Kλ / intercept

    Minimum 3 peaks required for meaningful linear regression.

    Reference: [R5] Williamson & Hall (1953)
    """
    two_theta_deg = np.asarray(two_theta_deg, dtype=float)
    fwhm_rad      = np.asarray(fwhm_rad,      dtype=float)

    if len(two_theta_deg) < 3:
        return None

    theta_rad = np.radians(two_theta_deg / 2.0)
    X = 4.0 * np.sin(theta_rad)          # x-axis: 4 sin θ
    Y = fwhm_rad * np.cos(theta_rad)     # y-axis: β cos θ  [rad]

    slope, intercept, r, p_val, se = linregress(X, Y)

    D_wh_nm = None
    if intercept > 1e-10:
        D_wh_nm = float((K * wavelength) / intercept) / 10.0  # Å → nm

    # Residuals for error estimate
    Y_pred   = slope * X + intercept
    residuals = Y - Y_pred

    return {
        "x_wh"             : X,
        "y_wh"             : Y,
        "slope"            : float(slope),
        "intercept"        : float(intercept),
        "r_squared"        : float(r ** 2),
        "p_value"          : float(p_val),
        "std_error"        : float(se),
        "residuals"        : residuals,
        "D_wh_nm"          : D_wh_nm,
        "microstrain"      : float(slope),          # ε from slope
        "microstrain_pct"  : float(abs(slope) * 100),
        "strain_sign"      : "tensile" if slope >= 0 else "compressive",
        "valid"            : (intercept > 1e-10)    # physical check
    }


# ─────────────────────────────────────────────────────────────────────────────
# Size-Strain Plot (SSP) — alternative to Williamson-Hall
# ─────────────────────────────────────────────────────────────────────────────
def size_strain_plot(two_theta_deg: np.ndarray, fwhm_rad: np.ndarray,
                     wavelength: float, K: float = 0.9) -> dict:
    """
    Size-Strain Plot (SSP) analysis.
        (d_hkl · β · cos θ)² = (K/D_s) · (d_hkl² · β · cos θ) + (1.5 ε)²

    Y = (d · β · cos θ)²
    X = d² · β · cos θ
    → Slope     = K/D_s    (area-weighted size)
    → Intercept = (1.5 ε)² → ε_rms = √intercept / 1.5

    Reference: [R10] Mote et al. (2012)
    """
    two_theta_deg = np.asarray(two_theta_deg, dtype=float)
    fwhm_rad      = np.asarray(fwhm_rad,      dtype=float)

    if len(two_theta_deg) < 3:
        return None

    theta_rad = np.radians(two_theta_deg / 2.0)
    d         = bragg_law(two_theta_deg, wavelength)   # Å
    bc        = fwhm_rad * np.cos(theta_rad)           # β cos θ  [rad]

    X = d ** 2 * bc             # Å² · rad
    Y = (d * bc) ** 2           # Å² · rad²

    slope, intercept, r, p_val, se = linregress(X, Y)

    D_ssp_nm  = None
    eps_rms   = None
    if slope > 1e-20:
        D_ssp_nm = float((K * wavelength) / slope) / 10.0  # Å → nm
    if intercept >= 0:
        eps_rms  = float(np.sqrt(intercept) / 1.5)

    return {
        "x_ssp"   : X,
        "y_ssp"   : Y,
        "slope"   : float(slope),
        "intercept": float(intercept),
        "r_squared": float(r ** 2),
        "D_ssp_nm": D_ssp_nm,
        "eps_rms" : eps_rms,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Per-Peak Microstrain (Stokes & Wilson)
# ─────────────────────────────────────────────────────────────────────────────
def microstrain_stokes_wilson(fwhm_rad: np.ndarray,
                               two_theta_deg: np.ndarray) -> np.ndarray:
    """
    Per-peak microstrain estimate:
        ε = β / (4·tan θ)

    This is equivalent to  ε = β·cos θ / (4·sin θ).

    Reference: [R4] Stokes & Wilson (1944)
    """
    theta_rad = np.radians(np.asarray(two_theta_deg, dtype=float) / 2.0)
    return np.asarray(fwhm_rad, dtype=float) / (4.0 * np.tan(theta_rad))


# ─────────────────────────────────────────────────────────────────────────────
# Dislocation Density
# ─────────────────────────────────────────────────────────────────────────────
def dislocation_density(D_nm: np.ndarray) -> np.ndarray:
    """
    Dislocation density from crystallite size (Williamson-Smallman formula):
        δ = 1 / D²  [m⁻²]

    Reference: [R6] Williamson & Smallman (1956)
    """
    D_m = np.asarray(D_nm, dtype=float) * 1e-9   # nm → m
    return 1.0 / (D_m ** 2)


# ─────────────────────────────────────────────────────────────────────────────
# Crystallinity Estimation
# ─────────────────────────────────────────────────────────────────────────────
def estimate_crystallinity(two_theta: np.ndarray, intensity: np.ndarray,
                            peaks: np.ndarray) -> tuple:
    """
    Estimate degree of crystallinity via area-ratio method.
        CI (%) = (Σ Crystalline Peak Areas) / (Total Corrected Area) × 100

    Background is estimated using the SNIP method.

    Reference: Alexander (1969) X-ray Diffraction Methods in Polymer Science.
    """
    bg = estimate_background(two_theta, intensity, n_iter=30)
    intensity_corr = np.maximum(intensity - bg, 0.0)
    total_area = float(np.trapezoid(intensity_corr, two_theta))

    if total_area < 1e-12:
        return 0.0, 100.0

    step = float(np.median(np.diff(two_theta)))
    cryst_area = 0.0

    for pk in peaks:
        hw = max(25, int(2.0 / step))
        s  = max(0, pk - hw)
        e  = min(len(intensity_corr), pk + hw + 1)
        xr = two_theta[s:e]
        yr = intensity_corr[s:e]
        # Local linear baseline inside this window
        local_bg = np.linspace(yr[0], yr[-1], len(yr))
        cryst_area += float(np.trapezoid(np.maximum(yr - local_bg, 0.0), xr))

    ci = min(100.0, (cryst_area / total_area) * 100.0)
    return round(ci, 2), round(100.0 - ci, 2)


# ─────────────────────────────────────────────────────────────────────────────
# Phase Comparison with Reference
# ─────────────────────────────────────────────────────────────────────────────
def compare_with_reference(d_spacings_meas: np.ndarray,
                            intensities_meas: np.ndarray,
                            ref_df: pd.DataFrame,
                            wavelength: float,
                            tolerance_A: float = 0.05) -> dict:
    """
    Compare measured d-spacings with a reference pattern.
    Returns match list and a match-quality score (0–1).
    """
    ref_2th  = ref_df['two_theta_ref'].values
    ref_int  = ref_df['intensity_ref'].values
    ref_int_norm = ref_int / (ref_int.max() + 1e-12) * 100.0
    ref_d    = bragg_law(ref_2th, wavelength)

    matches = []
    matched_ref = set()
    for i, d_m in enumerate(d_spacings_meas):
        diffs = np.abs(ref_d - d_m)
        j     = int(np.argmin(diffs))
        if diffs[j] < tolerance_A and j not in matched_ref:
            matched_ref.add(j)
            matches.append({
                "peak_no"      : i + 1,
                "d_measured_A" : round(float(d_m), 4),
                "d_ref_A"      : round(float(ref_d[j]), 4),
                "Δd_A"         : round(float(d_m - ref_d[j]), 4),
                "ref_2theta"   : round(float(ref_2th[j]), 3),
                "ref_int_pct"  : round(float(ref_int_norm[j]), 1),
                "meas_int"     : round(float(intensities_meas[i]), 1),
            })

    quality = len(matches) / max(len(d_spacings_meas), 1)
    return {"matches": matches,
            "n_matches": len(matches),
            "n_measured": len(d_spacings_meas),
            "match_quality_pct": round(quality * 100, 1)}


# ─────────────────────────────────────────────────────────────────────────────
# Texture Coefficient (Harris Analysis)
# ─────────────────────────────────────────────────────────────────────────────
def texture_coefficient(intensities_meas: np.ndarray,
                         intensities_ref: np.ndarray) -> np.ndarray:
    """
    Harris Texture Coefficient:
        TC(hkl) = [I(hkl)/I₀(hkl)] / [(1/N)·Σ(I(hkl)/I₀(hkl))]

    TC = 1 → random orientation
    TC > 1 → preferred orientation in that direction

    Reference: [R8] Harris (1952)
    """
    I    = np.asarray(intensities_meas,  dtype=float)
    I0   = np.asarray(intensities_ref,   dtype=float)
    I0   = np.maximum(I0, 1e-12)
    ratios = I / I0
    return ratios / (np.mean(ratios) + 1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# Plot Helpers
# ─────────────────────────────────────────────────────────────────────────────
STYLE = {
    "bg"      : "#FAFAFA",
    "grid_c"  : "#E0E0E0",
    "data_c"  : "#2563EB",   # blue
    "peak_c"  : "#DC2626",   # red
    "fit_c"   : "#D97706",   # amber
    "line1"   : "#7C3AED",   # violet
    "line2"   : "#059669",   # emerald
    "line3"   : "#0891B2",   # cyan
    "title_fs": 13,
    "label_fs": 11,
    "tick_fs" : 9,
}

def _apply_style(ax, title="", xlabel="", ylabel=""):
    ax.set_facecolor(STYLE["bg"])
    ax.grid(True, color=STYLE["grid_c"], linewidth=0.6, linestyle="--", alpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=STYLE["tick_fs"])
    if title:  ax.set_title(title,  fontsize=STYLE["title_fs"], fontweight="bold", pad=8)
    if xlabel: ax.set_xlabel(xlabel, fontsize=STYLE["label_fs"])
    if ylabel: ax.set_ylabel(ylabel, fontsize=STYLE["label_fs"])


def plot_raw_pattern(two_theta, intensity, peaks, material_name="Sample"):
    """Full XRD diffractogram with annotated detected peaks."""
    fig, ax = plt.subplots(figsize=(13, 5), dpi=140)
    fig.patch.set_facecolor("white")

    ax.plot(two_theta, intensity, color=STYLE["data_c"], lw=0.9,
            label="XRD Pattern", zorder=2)

    if len(peaks):
        ax.scatter(two_theta[peaks], intensity[peaks],
                   color=STYLE["peak_c"], s=55, zorder=5,
                   label=f"Detected Peaks  (n = {len(peaks)})", marker="v")
        for pk in peaks:
            ax.annotate(f"{two_theta[pk]:.2f}°",
                        xy=(two_theta[pk], intensity[pk]),
                        xytext=(0, 12), textcoords="offset points",
                        ha="center", fontsize=7, color=STYLE["peak_c"],
                        rotation=60, clip_on=True)

    _apply_style(ax,
                 title=f"XRD Diffractogram — {material_name}",
                 xlabel="2θ (°)",
                 ylabel="Intensity (a.u.)")
    ax.set_xlim(two_theta.min(), two_theta.max())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.legend(fontsize=9)
    plt.tight_layout()
    return fig


def plot_peak_fits(two_theta, intensity, peaks, fit_results, wavelength, material_name="Sample"):
    """Grid of individual peak fits."""
    n = min(len(peaks), 12)
    if n == 0:
        return None
    cols = min(4, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.5 * cols, 3.5 * rows), dpi=130)
    fig.patch.set_facecolor("white")
    axes_flat = np.array(axes).flatten() if n > 1 else [axes]

    for i in range(n):
        ax  = axes_flat[i]
        fr  = fit_results[i]
        xr  = fr["x_region"]
        yr  = fr["y_region"]
        xf  = fr["x_fit"]
        yf  = fr["y_fit"]
        d   = bragg_law(fr["center"], wavelength)

        ax.plot(xr, yr, ".", ms=3, color=STYLE["data_c"], alpha=0.8, label="Data")
        ax.plot(xf, yf, "-", lw=1.8, color=STYLE["fit_c"],
                label=f"{fr['type'].replace('_',' ').title()}\nR²={fr['r_squared']:.4f}")

        # FWHM bar
        half = fr["amplitude"] / 2.0 + fr.get("baseline", 0)
        ax.axhline(y=half, color=STYLE["line2"], ls=":", lw=0.9, alpha=0.8)

        ax.set_title(f"2θ={fr['center']:.2f}°  d={d:.3f} Å\nFWHM={fr['fwhm_deg']:.4f}°",
                     fontsize=7.5, fontweight="bold")
        ax.set_xlabel("2θ (°)", fontsize=7)
        ax.set_ylabel("Intensity", fontsize=7)
        ax.legend(fontsize=6, loc="upper right")
        _apply_style(ax)

    # Hide unused subplots
    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle(f"Peak Fitting (Pseudo-Voigt) — {material_name}",
                 fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    return fig


def plot_scherrer(two_theta_peaks, D_nm, material_name="Sample"):
    """Crystallite size vs 2θ with distribution histogram."""
    if len(D_nm) == 0:
        return None
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), dpi=140)
    fig.patch.set_facecolor("white")

    D_mean = float(np.mean(D_nm))
    D_std  = float(np.std(D_nm))

    ax1.plot(two_theta_peaks, D_nm, "o-", color=STYLE["data_c"],
             ms=8, mfc="white", mew=2, lw=1.5)
    ax1.axhline(D_mean, color=STYLE["peak_c"], ls="--", lw=1.5,
                label=f"Mean D = {D_mean:.2f} nm")
    ax1.fill_between(two_theta_peaks,
                     D_mean - D_std, D_mean + D_std,
                     alpha=0.15, color=STYLE["peak_c"],
                     label=f"±1σ = {D_std:.2f} nm")
    _apply_style(ax1, "Scherrer Crystallite Size vs 2θ", "2θ (°)", "D (nm)")
    ax1.legend(fontsize=9)

    bins = min(max(3, len(D_nm)), 15)
    ax2.hist(D_nm, bins=bins, color=STYLE["data_c"], edgecolor="white",
             alpha=0.8, linewidth=0.8)
    ax2.axvline(D_mean, color=STYLE["peak_c"], ls="--", lw=2,
                label=f"Mean = {D_mean:.2f} nm")
    _apply_style(ax2, "Crystallite Size Distribution", "D (nm)", "Count")
    ax2.legend(fontsize=9)

    fig.suptitle(f"Scherrer Analysis — {material_name}",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    return fig


def plot_williamson_hall(wh: dict, material_name="Sample"):
    """Williamson-Hall (UDM) plot with linear regression."""
    if wh is None:
        return None
    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)
    fig.patch.set_facecolor("white")

    X, Y = wh["x_wh"], wh["y_wh"]
    ax.scatter(X, Y, s=90, color=STYLE["data_c"], edgecolors="white",
               lw=1.5, zorder=5, label=f"Data (n={len(X)})")

    xf = np.linspace(X.min() * 0.9, X.max() * 1.1, 300)
    yf = wh["slope"] * xf + wh["intercept"]
    ax.plot(xf, yf, color=STYLE["peak_c"], lw=2,
            label=f"Fit  slope={wh['slope']:.6f}  "
                  f"R²={wh['r_squared']:.4f}")

    # Equation box
    D_str = f"{wh['D_wh_nm']:.2f} nm" if wh.get("D_wh_nm") else "N/A (negative intercept)"
    box_txt = (
        f"UDM Williamson–Hall\n"
        f"β·cosθ = (Kλ/D) + 4ε·sinθ\n\n"
        f"D$_{{WH}}$ = {D_str}\n"
        f"ε  = {wh['microstrain']:.5f}  ({wh['strain_sign']})\n"
        f"R² = {wh['r_squared']:.4f}"
    )
    ax.text(0.04, 0.96, box_txt, transform=ax.transAxes,
            fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="#ccc", alpha=0.9))

    _apply_style(ax,
                 title=f"Williamson–Hall Plot (UDM) — {material_name}",
                 xlabel="4 · sin θ",
                 ylabel="β · cos θ  (rad)")
    ax.legend(fontsize=9, loc="lower right")
    plt.tight_layout()
    return fig


def plot_ssp(ssp: dict, material_name="Sample"):
    """Size-Strain Plot (SSP)."""
    if ssp is None:
        return None
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=140)
    fig.patch.set_facecolor("white")

    X, Y = ssp["x_ssp"], ssp["y_ssp"]
    ax.scatter(X, Y, s=80, color=STYLE["line1"], edgecolors="white",
               lw=1.5, zorder=5, label=f"Data (n={len(X)})")
    xf = np.linspace(X.min() * 0.9, X.max() * 1.1, 200)
    ax.plot(xf, ssp["slope"] * xf + ssp["intercept"],
            color=STYLE["peak_c"], lw=2, label=f"Fit  R²={ssp['r_squared']:.4f}")

    D_str   = f"{ssp['D_ssp_nm']:.2f} nm" if ssp.get("D_ssp_nm") else "N/A"
    eps_str = f"{ssp['eps_rms']:.5f}" if ssp.get("eps_rms") is not None else "N/A"
    box = (f"Size–Strain Plot\n"
           f"D$_{{SSP}}$ = {D_str}\n"
           f"ε$_{{rms}}$ = {eps_str}\n"
           f"R² = {ssp['r_squared']:.4f}")
    ax.text(0.04, 0.96, box, transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc="#eef", ec="#aaa", alpha=0.9))

    _apply_style(ax, f"Size–Strain Plot — {material_name}",
                 "d² · β · cosθ  (Å² · rad)",
                 "(d · β · cosθ)²  (Å² · rad²)")
    ax.legend(fontsize=9)
    plt.tight_layout()
    return fig


def plot_microstrain(two_theta_peaks, epsilon, material_name="Sample"):
    """Per-peak microstrain bar chart."""
    if len(epsilon) == 0:
        return None
    fig, ax = plt.subplots(figsize=(11, 4.5), dpi=140)
    fig.patch.set_facecolor("white")

    eps_k = epsilon * 1e3
    colors = [STYLE["line2"] if e >= 0 else STYLE["peak_c"] for e in eps_k]
    bars = ax.bar(range(len(two_theta_peaks)), eps_k,
                  color=colors, edgecolor="white", lw=0.8, alpha=0.85)
    ax.axhline(np.mean(eps_k), color="#666", ls="--", lw=1.5,
               label=f"Mean ε = {np.mean(eps_k):.4f} ×10⁻³")
    ax.axhline(0, color="black", lw=0.8)

    ax.set_xticks(range(len(two_theta_peaks)))
    ax.set_xticklabels([f"{t:.2f}°" for t in two_theta_peaks],
                       rotation=45, ha="right", fontsize=8)
    for bar, val in zip(bars, eps_k):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2,
                h + 0.01 * abs(eps_k.max() - eps_k.min()),
                f"{val:.3f}", ha="center", va="bottom", fontsize=7)

    _apply_style(ax, f"Microstrain (Stokes–Wilson) — {material_name}",
                 "2θ Position", "Microstrain ε (×10⁻³)")
    ax.legend(fontsize=9)
    plt.tight_layout()
    return fig


def plot_dislocation_density(two_theta_peaks, delta, material_name="Sample"):
    """Dislocation density per peak (semi-log)."""
    if len(delta) == 0:
        return None
    fig, ax = plt.subplots(figsize=(11, 4.5), dpi=140)
    fig.patch.set_facecolor("white")

    ax.semilogy(two_theta_peaks, delta, "s-",
                color=STYLE["line1"], ms=9, mfc="white", mew=2, lw=1.5)
    ax.axhline(np.mean(delta), color=STYLE["peak_c"], ls="--", lw=1.5,
               label=f"Mean δ = {np.mean(delta):.3e} m⁻²")

    _apply_style(ax, f"Dislocation Density — {material_name}",
                 "2θ (°)", "δ (m⁻²)")
    ax.legend(fontsize=9)
    plt.tight_layout()
    return fig


def plot_reference_comparison(two_theta, intensity, peaks, ref_df,
                               wavelength, material_name="Sample"):
    """Overlay measured pattern with reference stick pattern."""
    if ref_df is None:
        return None
    fig, ax = plt.subplots(figsize=(13, 5), dpi=140)
    fig.patch.set_facecolor("white")

    int_norm = intensity / intensity.max() * 100
    ax.plot(two_theta, int_norm, color=STYLE["data_c"], lw=0.9,
            label="Measured (normalised)", zorder=2)

    ref_2th = ref_df["two_theta_ref"].values
    ref_int = ref_df["intensity_ref"].values
    ref_int_n = ref_int / ref_int.max() * 100
    ax.vlines(ref_2th, 0, -ref_int_n * 0.4,
              colors=STYLE["peak_c"], lw=1.8, alpha=0.75,
              label="Reference (stick, inverted)", zorder=3)

    if len(peaks):
        ax.scatter(two_theta[peaks], int_norm[peaks],
                   marker="^", s=60, color=STYLE["line2"], zorder=6,
                   label=f"Detected Peaks (n={len(peaks)})")

    ax.axhline(0, color="#888", lw=0.8)
    _apply_style(ax, f"Measured vs Reference — {material_name}",
                 "2θ (°)", "Normalised Intensity (%)")
    ax.legend(fontsize=9)
    ax.set_xlim(two_theta.min(), two_theta.max())
    plt.tight_layout()
    return fig


def plot_summary_dashboard(results: dict, material_name="Sample"):
    """4-panel summary dashboard."""
    fig = plt.figure(figsize=(14, 9), dpi=130)
    fig.patch.set_facecolor("white")
    gs  = gridspec.GridSpec(2, 2, hspace=0.45, wspace=0.35)

    peaks      = results["peaks"]
    two_theta  = results["two_theta"]
    intensity  = results["intensity"]
    D_nm       = results["D_nm"]
    epsilon    = results["epsilon"]
    wh         = results.get("wh")

    # ── Panel 1: Raw pattern ──────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(two_theta, intensity, color=STYLE["data_c"], lw=0.8)
    if len(peaks):
        ax1.scatter(two_theta[peaks], intensity[peaks],
                    color=STYLE["peak_c"], s=45, marker="v", zorder=5)
    _apply_style(ax1, f"XRD Pattern — {material_name}", "2θ (°)", "Intensity (a.u.)")
    ax1.set_xlim(two_theta.min(), two_theta.max())

    # ── Panel 2: Crystallite size ─────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[1, 0])
    if len(D_nm):
        ax2.bar(range(len(D_nm)), D_nm, color=STYLE["line2"],
                edgecolor="white", alpha=0.85)
        ax2.axhline(np.mean(D_nm), color=STYLE["peak_c"], ls="--", lw=1.5,
                    label=f"Mean = {np.mean(D_nm):.2f} nm")
        ax2.set_xticks(range(len(results["two_theta_peaks"])))
        ax2.set_xticklabels([f"{t:.1f}°" for t in results["two_theta_peaks"]],
                             rotation=45, ha="right", fontsize=7)
        ax2.legend(fontsize=8)
    _apply_style(ax2, "Crystallite Size D (Scherrer)", "Peak (2θ)", "D (nm)")

    # ── Panel 3: Williamson-Hall ──────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 1])
    if wh and len(wh.get("x_wh", [])) >= 3:
        X, Y = wh["x_wh"], wh["y_wh"]
        ax3.scatter(X, Y, s=70, color=STYLE["data_c"], edgecolors="white", lw=1.2, zorder=5)
        xf = np.linspace(X.min(), X.max(), 200)
        ax3.plot(xf, wh["slope"] * xf + wh["intercept"],
                 color=STYLE["peak_c"], lw=1.8)
        ax3.set_title(f"Williamson–Hall  R²={wh['r_squared']:.3f}",
                      fontsize=STYLE["title_fs"] - 1, fontweight="bold")
        ax3.set_xlabel("4 · sin θ", fontsize=STYLE["label_fs"])
        ax3.set_ylabel("β · cos θ (rad)", fontsize=STYLE["label_fs"])
        _apply_style(ax3)
    else:
        ax3.text(0.5, 0.5, "≥ 3 peaks required for\nWilliamson–Hall",
                 transform=ax3.transAxes, ha="center", va="center", fontsize=11,
                 color="#888")
        _apply_style(ax3, "Williamson–Hall")

    fig.suptitle(f"XRD Analysis Summary — {material_name}",
                 fontsize=14, fontweight="bold")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Master Analysis Function
# ─────────────────────────────────────────────────────────────────────────────
def run_xrd_analysis(
    xrd_file_path   : str,
    ref_file_path   : str   = None,
    wavelength      : float = 1.5406,
    material_name   : str   = "Sample",
    min_ht_frac     : float = 0.05,
    min_prom_frac   : float = 0.03,
    min_dist_deg    : float = 0.5,
    K               : float = 0.9,
    beta_instr_deg  : float = 0.0,
    two_theta_min   : float = None,
    two_theta_max   : float = None,
) -> tuple:
    """
    Run the complete XRD analysis pipeline.

    Returns
    -------
    results : dict — all numerical results
    figures : dict — all matplotlib figures
    summary : str  — human-readable text summary
    """
    # ── 1. Load ──────────────────────────────────────────────────────────────
    two_theta, intensity, col_info = load_xrd_data(xrd_file_path)

    # ── 2. Trim range ─────────────────────────────────────────────────────────
    mask = np.ones(len(two_theta), dtype=bool)
    if two_theta_min is not None:
        mask &= two_theta >= two_theta_min
    if two_theta_max is not None:
        mask &= two_theta <= two_theta_max
    two_theta = two_theta[mask]
    intensity = intensity[mask]

    if len(two_theta) < 10:
        raise ValueError("After 2θ range filter, fewer than 10 data points remain.")

    # ── 3. Detect peaks ───────────────────────────────────────────────────────
    peaks, _props = detect_peaks(two_theta, intensity,
                                  min_ht_frac, min_prom_frac, min_dist_deg)

    # ── 4. Fit each peak ──────────────────────────────────────────────────────
    fit_results = [fit_single_peak(two_theta, intensity, pk) for pk in peaks]

    # Extract centre positions and FWHM
    centers     = np.array([fr["center"]    for fr in fit_results])
    fwhm_deg    = np.array([fr["fwhm_deg"]  for fr in fit_results])
    fwhm_rad_m  = np.array([fr["fwhm_rad"]  for fr in fit_results])
    ib_deg      = np.array([fr.get("integral_breadth_deg", fr["fwhm_deg"]) for fr in fit_results])
    fit_types   = [fr["type"] for fr in fit_results]
    r_sq        = np.array([fr["r_squared"] for fr in fit_results])
    int_at_peak = intensity[peaks]

    # ── 5. Instrument-broadening correction ───────────────────────────────────
    fwhm_rad_corr = correct_instrument_broadening(fwhm_rad_m, beta_instr_deg)

    # ── 6. Bragg's law → d-spacings ──────────────────────────────────────────
    d_spacings = bragg_law(centers, wavelength)

    # ── 7. Scherrer crystallite size ──────────────────────────────────────────
    D_nm = scherrer_size(fwhm_rad_corr, centers, wavelength, K)

    # ── 8. Williamson-Hall ───────────────────────────────────────────────────
    wh = williamson_hall(centers, fwhm_rad_corr, wavelength, K)

    # ── 9. Size-Strain Plot ───────────────────────────────────────────────────
    ssp = size_strain_plot(centers, fwhm_rad_corr, wavelength, K)

    # ── 10. Microstrain per peak (Stokes-Wilson) ──────────────────────────────
    epsilon = microstrain_stokes_wilson(fwhm_rad_corr, centers)

    # ── 11. Dislocation density ───────────────────────────────────────────────
    delta = dislocation_density(D_nm)

    # ── 12. Crystallinity ─────────────────────────────────────────────────────
    ci, ai = estimate_crystallinity(two_theta, intensity, peaks)

    # ── 13. Reference comparison ─────────────────────────────────────────────
    ref_df       = None
    ref_compare  = None
    tc_values    = None
    if ref_file_path:
        try:
            ref_df = load_reference(ref_file_path)
            ref_compare = compare_with_reference(
                d_spacings, int_at_peak, ref_df, wavelength
            )
            # Texture coefficient (if matching peaks exist)
            if ref_compare["matches"] and len(ref_compare["matches"]) > 0:
                match_meas = [m["meas_int"]     for m in ref_compare["matches"]]
                match_ref  = [m["ref_int_pct"]  for m in ref_compare["matches"]]
                if all(v > 0 for v in match_ref):
                    tc_values = texture_coefficient(
                        np.array(match_meas),
                        np.array(match_ref) * int_at_peak.max() / 100.0
                    )
        except Exception as exc:
            ref_compare = {"error": str(exc)}

    # ── Assemble results dict ─────────────────────────────────────────────────
    results = {
        "two_theta"       : two_theta,
        "intensity"       : intensity,
        "col_info"        : col_info,
        "peaks"           : peaks,
        "two_theta_peaks" : centers,
        "int_at_peaks"    : int_at_peak,
        "fit_results"     : fit_results,
        "fit_types"       : fit_types,
        "r_squared_fits"  : r_sq,
        "fwhm_deg"        : fwhm_deg,
        "fwhm_rad"        : fwhm_rad_corr,
        "integral_breadth_deg" : ib_deg,
        "d_spacings"      : d_spacings,
        "D_nm"            : D_nm,
        "D_mean_nm"       : float(np.mean(D_nm)) if len(D_nm) else None,
        "D_std_nm"        : float(np.std(D_nm))  if len(D_nm) else None,
        "wh"              : wh,
        "ssp"             : ssp,
        "epsilon"         : epsilon,
        "delta"           : delta,
        "crystallinity_pct"  : ci,
        "amorphous_pct"      : ai,
        "ref_df"          : ref_df,
        "ref_compare"     : ref_compare,
        "tc_values"       : tc_values,
        "wavelength"      : wavelength,
        "K"               : K,
        "beta_instr_deg"  : beta_instr_deg,
        "material_name"   : material_name,
        "analysis_time"   : datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    # ── Build figures ─────────────────────────────────────────────────────────
    figures = {}
    figures["01_raw_pattern"]     = plot_raw_pattern(two_theta, intensity, peaks, material_name)
    figures["02_peak_fits"]       = plot_peak_fits(two_theta, intensity, peaks, fit_results, wavelength, material_name)
    figures["03_scherrer"]        = plot_scherrer(centers, D_nm, material_name)
    figures["04_williamson_hall"] = plot_williamson_hall(wh, material_name)
    figures["05_size_strain_plot"] = plot_ssp(ssp, material_name)
    figures["06_microstrain"]     = plot_microstrain(centers, epsilon, material_name)
    figures["07_dislocation"]     = plot_dislocation_density(centers, delta, material_name)
    if ref_df is not None:
        figures["08_reference"]   = plot_reference_comparison(two_theta, intensity, peaks, ref_df, wavelength, material_name)
    figures["00_summary_dashboard"] = plot_summary_dashboard(results, material_name)

    summary = _build_summary(results)
    return results, figures, summary


# ─────────────────────────────────────────────────────────────────────────────
# Results Table
# ─────────────────────────────────────────────────────────────────────────────
def build_results_dataframe(results: dict) -> pd.DataFrame:
    """Build the main per-peak results DataFrame."""
    peaks      = results["peaks"]
    two_theta  = results["two_theta"]
    intensity  = results["intensity"]
    centers    = results["two_theta_peaks"]
    d_sp       = results["d_spacings"]
    fwhm_d     = results["fwhm_deg"]
    fwhm_r     = results["fwhm_rad"]
    ib         = results["integral_breadth_deg"]
    D          = results["D_nm"]
    eps        = results["epsilon"]
    delta      = results["delta"]
    r_sq       = results["r_squared_fits"]
    ftypes     = results["fit_types"]
    int_pk     = results["int_at_peaks"]

    rows = []
    for i in range(len(peaks)):
        rows.append({
            "Peak #"                   : i + 1,
            "2θ (°)"                  : round(float(centers[i]), 4),
            "θ (°)"                   : round(float(centers[i]) / 2, 4),
            "Intensity (a.u.)"        : round(float(int_pk[i]), 2),
            "d-spacing (Å)"           : round(float(d_sp[i]), 4),
            "FWHM (°)"                : round(float(fwhm_d[i]), 5),
            "FWHM (rad)"              : round(float(fwhm_r[i]), 7),
            "Integral Breadth (°)"    : round(float(ib[i]), 5),
            "Fit Type"                : ftypes[i],
            "Fit R²"                  : round(float(r_sq[i]), 4),
            "Crystallite Size D (nm)" : round(float(D[i]), 3),
            "Microstrain ε"           : f"{float(eps[i]):.6e}",
            "Microstrain ε ×10⁻³"    : round(float(eps[i]) * 1e3, 4),
            "Dislocation δ (m⁻²)"    : f"{float(delta[i]):.4e}",
        })
    return pd.DataFrame(rows)


def build_summary_dataframe(results: dict) -> pd.DataFrame:
    """One-row statistical summary."""
    wh  = results.get("wh") or {}
    ssp = results.get("ssp") or {}
    rows = [
        ["Number of Detected Peaks",        len(results["peaks"])],
        ["Wavelength (Å)",                  results["wavelength"]],
        ["Scherrer K Factor",               results["K"]],
        ["Instrument FWHM Correction (°)",  results["beta_instr_deg"]],
        ["─── Scherrer ─────────────────",  ""],
        ["Mean Crystallite Size D (nm)",    f"{results['D_mean_nm']:.3f}" if results['D_mean_nm'] else "N/A"],
        ["Std Dev D (nm)",                  f"{results['D_std_nm']:.3f}"  if results['D_std_nm']  else "N/A"],
        ["─── Williamson–Hall ───────────",  ""],
        ["WH Crystallite Size D_WH (nm)",   f"{wh.get('D_wh_nm', 'N/A'):.3f}" if wh.get("D_wh_nm") else "N/A"],
        ["WH Microstrain ε",                f"{wh.get('microstrain', 'N/A'):.6e}" if wh.get("microstrain") is not None else "N/A"],
        ["WH Strain (%)",                   f"{wh.get('microstrain_pct', 0):.4f}" if wh.get("microstrain_pct") is not None else "N/A"],
        ["WH Strain Character",             wh.get("strain_sign", "N/A")],
        ["WH R²",                           f"{wh.get('r_squared', 0):.4f}" if wh.get("r_squared") is not None else "N/A"],
        ["─── Size–Strain Plot ──────────",  ""],
        ["SSP Crystallite Size D_SSP (nm)", f"{ssp.get('D_ssp_nm', 'N/A'):.3f}" if ssp.get("D_ssp_nm") else "N/A"],
        ["SSP RMS Strain ε_rms",            f"{ssp.get('eps_rms', 'N/A'):.6e}"  if ssp.get("eps_rms") else "N/A"],
        ["SSP R²",                          f"{ssp.get('r_squared', 0):.4f}"     if ssp.get("r_squared") is not None else "N/A"],
        ["─── Microstructure ────────────",  ""],
        ["Mean Microstrain ε (Stokes–Wilson)",f"{float(np.mean(results['epsilon'])):.6e}" if len(results['epsilon']) else "N/A"],
        ["Mean Dislocation Density δ (m⁻²)", f"{float(np.mean(results['delta'])):.4e}"   if len(results['delta']) else "N/A"],
        ["─── Crystallinity ─────────────",  ""],
        ["Crystallinity (%)",               results["crystallinity_pct"]],
        ["Amorphous Fraction (%)",          results["amorphous_pct"]],
    ]
    if results.get("ref_compare") and "match_quality_pct" in (results["ref_compare"] or {}):
        rc = results["ref_compare"]
        rows += [
            ["─── Reference Match ───────────", ""],
            ["Matching Peaks",                f"{rc['n_matches']} / {rc['n_measured']}"],
            ["Phase Match Quality (%)",       rc["match_quality_pct"]],
        ]
    return pd.DataFrame(rows, columns=["Parameter", "Value"])


# ─────────────────────────────────────────────────────────────────────────────
# Text Summary
# ─────────────────────────────────────────────────────────────────────────────
def _build_summary(results: dict) -> str:
    wh  = results.get("wh") or {}
    ssp = results.get("ssp") or {}
    rc  = results.get("ref_compare") or {}
    n_pk = len(results["peaks"])
    D_nm = results["D_nm"]
    eps  = results["epsilon"]
    dlt  = results["delta"]

    L  = []
    hr = "═" * 68
    L += [hr,
          f"  XRD ANALYSIS REPORT  ·  {results['material_name'].upper()}",
          hr,
          f"  Date       : {results['analysis_time']}",
          f"  Wavelength : {results['wavelength']} Å",
          f"  K factor   : {results['K']}",
          f"  Instr. β   : {results['beta_instr_deg']}°",
          ""]

    L += ["── DATA OVERVIEW " + "─"*50,
          f"  2θ range   : {results['two_theta'].min():.2f}° – {results['two_theta'].max():.2f}°",
          f"  Step size  : {float(np.median(np.diff(results['two_theta']))):.4f}°",
          f"  Data pts   : {len(results['two_theta'])}",
          ""]

    L += [f"── PEAKS DETECTED  ({n_pk} peaks) " + "─"*34]
    for i in range(n_pk):
        L.append(f"  #{i+1:>2}  2θ={results['two_theta_peaks'][i]:.3f}°"
                 f"   d={results['d_spacings'][i]:.4f} Å"
                 f"   FWHM={results['fwhm_deg'][i]:.4f}°"
                 f"   R²={results['r_squared_fits'][i]:.4f}")
    L.append("")

    L += ["── SCHERRER CRYSTALLITE SIZE " + "─"*38]
    if len(D_nm):
        L += [f"  Mean D       : {np.mean(D_nm):.3f} nm",
              f"  Std Dev      : {np.std(D_nm):.3f} nm",
              f"  Range        : {np.min(D_nm):.3f} – {np.max(D_nm):.3f} nm"]
    L.append("")

    L += ["── WILLIAMSON–HALL ANALYSIS (UDM) " + "─"*33]
    if wh:
        D_wh_str = f"{wh['D_wh_nm']:.3f} nm" if wh.get("D_wh_nm") else "N/A (non-positive intercept)"
        L += [f"  D_WH         : {D_wh_str}",
              f"  Microstrain ε: {wh.get('microstrain', 0):.6e}  ({wh.get('strain_sign', 'N/A')})",
              f"  Strain (%)   : {wh.get('microstrain_pct', 0):.4f}",
              f"  R²           : {wh.get('r_squared', 0):.4f}",
              f"  Slope        : {wh.get('slope', 0):.6e}",
              f"  Intercept    : {wh.get('intercept', 0):.6e}"]
        if not wh.get("valid"):
            L.append("  ⚠ Negative intercept — WH size estimate unreliable;"
                     " use Scherrer value.")
    else:
        L.append("  Requires ≥ 3 peaks.")
    L.append("")

    L += ["── SIZE–STRAIN PLOT (SSP) " + "─"*41]
    if ssp and ssp.get("D_ssp_nm"):
        L += [f"  D_SSP        : {ssp['D_ssp_nm']:.3f} nm",
              f"  ε_rms        : {ssp.get('eps_rms', 'N/A'):.5e}" if ssp.get("eps_rms") else "  ε_rms      : N/A",
              f"  R²           : {ssp.get('r_squared', 0):.4f}"]
    else:
        L.append("  Requires ≥ 3 peaks or intercept < 0.")
    L.append("")

    L += ["── MICROSTRAIN (Stokes–Wilson) " + "─"*36]
    if len(eps):
        L += [f"  Mean ε       : {np.mean(eps):.6e}",
              f"  Std Dev ε    : {np.std(eps):.6e}"]
    L.append("")

    L += ["── DISLOCATION DENSITY " + "─"*44]
    if len(dlt):
        L += [f"  Mean δ       : {np.mean(dlt):.4e} m⁻²",
              f"  Std Dev δ    : {np.std(dlt):.4e} m⁻²"]
    L.append("")

    L += ["── CRYSTALLINITY " + "─"*50,
          f"  Crystalline  : {results['crystallinity_pct']:.2f} %",
          f"  Amorphous    : {results['amorphous_pct']:.2f} %",
          ""]

    if rc and "n_matches" in rc:
        L += ["── REFERENCE COMPARISON " + "─"*43,
              f"  Matches      : {rc['n_matches']} / {rc['n_measured']}",
              f"  Match Quality: {rc['match_quality_pct']} %",
              ""]

    L += [hr,
          "  SCIENTIFIC REFERENCES",
          hr,
          "  [R1]  Bragg (1913) Proc. R. Soc. A 88, 428–438",
          "  [R2]  Scherrer (1918) Nachr. Ges. Wiss. Göttingen, 98–100",
          "  [R3]  Patterson (1939) Phys. Rev. 56, 978–982",
          "  [R4]  Stokes & Wilson (1944) Proc. Phys. Soc. 56, 174",
          "  [R5]  Williamson & Hall (1953) Acta Metall. 1, 22–31",
          "  [R6]  Williamson & Smallman (1956) Phil. Mag. 1, 34–46",
          "  [R7]  Warren (1969) X-ray Diffraction, Addison-Wesley",
          "  [R8]  Harris (1952) Phil. Mag. 43, 113–123",
          "  [R9]  Thompson et al. (1987) J. Appl. Crystallogr. 20, 79–83",
          "  [R10] Mote et al. (2012) J. Theor. Appl. Phys. 6, 6",
          "  [R11] Cullity & Stock (2001) Elements of X-ray Diffraction",
          "  [R12] Venkateswarlu et al. (2010) Trans. Nonferrous Met. 20, 941",
          hr]

    return "\n".join(L)


# ─────────────────────────────────────────────────────────────────────────────
# ZIP Export
# ─────────────────────────────────────────────────────────────────────────────
def create_results_zip(results: dict, figures: dict, summary: str) -> str:
    """
    Pack all analysis outputs into a downloadable ZIP file.
    Returns the path to the created ZIP.
    """
    mat   = results["material_name"].replace(" ", "_")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    zip_path = f"/tmp/XRD_{mat}_{stamp}.zip"

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:

        # ── Plots ─────────────────────────────────────────────────────────────
        for name, fig in figures.items():
            if fig is None:
                continue
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=180,
                        bbox_inches="tight", facecolor="white")
            buf.seek(0)
            zf.writestr(f"plots/{name}.png", buf.read())
            plt.close(fig)

        # ── Excel workbook ────────────────────────────────────────────────────
        df_peaks = build_results_dataframe(results)
        df_summ  = build_summary_dataframe(results)

        buf_xl = io.BytesIO()
        with pd.ExcelWriter(buf_xl, engine="openpyxl") as writer:
            df_peaks.to_excel(writer, sheet_name="Peak Analysis", index=False)
            df_summ.to_excel(writer,  sheet_name="Summary",       index=False)
            # WH data
            wh = results.get("wh")
            if wh and len(wh.get("x_wh", [])):
                pd.DataFrame({"4·sinθ": wh["x_wh"],
                              "β·cosθ (rad)": wh["y_wh"]}).to_excel(
                    writer, sheet_name="Williamson-Hall", index=False)
            # SSP data
            ssp = results.get("ssp")
            if ssp and len(ssp.get("x_ssp", [])):
                pd.DataFrame({"d²·β·cosθ": ssp["x_ssp"],
                              "(d·β·cosθ)²": ssp["y_ssp"]}).to_excel(
                    writer, sheet_name="Size-Strain Plot", index=False)
            # Reference match
            rc = results.get("ref_compare")
            if rc and rc.get("matches"):
                pd.DataFrame(rc["matches"]).to_excel(
                    writer, sheet_name="Reference Match", index=False)
        buf_xl.seek(0)
        zf.writestr("results/XRD_Analysis.xlsx", buf_xl.read())

        # ── CSV ────────────────────────────────────────────────────────────────
        zf.writestr("results/peak_data.csv",
                    df_peaks.to_csv(index=False))

        # ── Summary text ───────────────────────────────────────────────────────
        zf.writestr("results/Analysis_Report.txt", summary)

        # ── JSON numerical results ─────────────────────────────────────────────
        json_out = _serialize_results(results)
        zf.writestr("results/analysis_data.json",
                    json.dumps(json_out, indent=2))

        # ── README inside ZIP ─────────────────────────────────────────────────
        zf.writestr("README.txt", _zip_readme(mat, stamp))

    return zip_path


def _serialize_results(results: dict) -> dict:
    """Convert numpy/pandas types to JSON-serialisable."""
    skip = {"two_theta", "intensity", "fit_results", "ref_df", "col_info"}
    out  = {}
    for k, v in results.items():
        if k in skip:
            continue
        if isinstance(v, np.ndarray):
            out[k] = v.tolist()
        elif isinstance(v, (np.floating, np.integer)):
            out[k] = float(v)
        elif isinstance(v, dict):
            inner = {}
            for kk, vv in v.items():
                if isinstance(vv, np.ndarray):
                    inner[kk] = vv.tolist()
                elif isinstance(vv, (np.floating, np.integer)):
                    inner[kk] = float(vv)
                elif isinstance(vv, pd.DataFrame):
                    continue
                else:
                    try:
                        json.dumps(vv)
                        inner[kk] = vv
                    except (TypeError, ValueError):
                        pass
            out[k] = inner
        elif isinstance(v, pd.DataFrame):
            continue
        else:
            try:
                json.dumps(v)
                out[k] = v
            except (TypeError, ValueError):
                pass
    return out


def _zip_readme(mat, stamp):
    return f"""XRD Analysis Results — {mat}
Generated : {stamp}
By        : CharacterizationCalc XRD Module

STRUCTURE
=========
plots/
  00_summary_dashboard.png   ─ 4-panel overview
  01_raw_pattern.png         ─ Diffractogram + detected peaks
  02_peak_fits.png           ─ Individual Pseudo-Voigt fits
  03_scherrer.png            ─ Scherrer size vs 2θ + histogram
  04_williamson_hall.png     ─ WH (UDM) plot
  05_size_strain_plot.png    ─ Size-Strain Plot
  06_microstrain.png         ─ Per-peak microstrain bar chart
  07_dislocation.png         ─ Dislocation density
  08_reference.png           ─ vs reference (if uploaded)

results/
  XRD_Analysis.xlsx          ─ Full data in Excel (5 sheets)
  peak_data.csv              ─ Per-peak CSV
  Analysis_Report.txt        ─ Text summary with references
  analysis_data.json         ─ Raw numerical output (JSON)

ANALYSES
========
 1. Peak detection          (scipy.signal.find_peaks)
 2. Pseudo-Voigt fitting    (Thompson et al. 1987)
 3. Integral breadth        (area / height)
 4. Bragg's Law             (d-spacing)
 5. Scherrer equation       (crystallite size D)
 6. Williamson–Hall UDM     (D + microstrain ε)
 7. Size–Strain Plot        (alternative size + strain)
 8. Microstrain per peak    (Stokes & Wilson 1944)
 9. Dislocation density     (Williamson & Smallman 1956)
10. Crystallinity estimate  (area-ratio method)
11. Phase comparison        (if reference provided)
12. Texture coefficient     (Harris 1952, if reference matches)

KEY REFERENCES
==============
Bragg (1913) · Scherrer (1918) · Patterson (1939)
Stokes & Wilson (1944) · Williamson & Hall (1953)
Williamson & Smallman (1956) · Harris (1952)
Thompson et al. (1987) · Mote et al. (2012)
Cullity & Stock (2001)
"""
