"""
================================================================================
  CharacterizationCalc — XRD Analysis Web Application
  Built with Gradio · Deployable on Render / HuggingFace Spaces
================================================================================
"""

import os
import traceback

import gradio as gr

# ── Gradio 4.x / pydantic >= 2.11 compatibility patch ────────────────────────
# Newer pydantic emits JSON-schema nodes like {"additionalProperties": true}.
# gradio_client 1.3.x assumes every schema node is a dict and crashes with
# "TypeError: argument of type 'bool' is not iterable" when building the API
# info, which makes the startup self-check fail (-> 502 on Render).
# requirements.txt also pins pydantic; this patch is a second safety net.
try:
    import gradio_client.utils as _gcu

    _orig_get_type = _gcu.get_type
    _orig_j2p = _gcu._json_schema_to_python_type

    def _safe_get_type(schema):
        if isinstance(schema, bool):
            return "any"
        return _orig_get_type(schema)

    def _safe_j2p(schema, defs=None):
        if isinstance(schema, bool):
            return "Any"
        return _orig_j2p(schema, defs)

    _gcu.get_type = _safe_get_type
    _gcu._json_schema_to_python_type = _safe_j2p
except Exception:  # pragma: no cover
    pass
import numpy as np
import pandas as pd

# ── Startup diagnostics: print library versions into the Render log ─────────
try:
    from importlib.metadata import version as _v
    print("[startup] " + ", ".join(
        f"{p}={_v(p)}" for p in ("gradio", "gradio_client", "fastapi",
                                 "starlette", "pydantic", "numpy",
                                 "scipy", "pandas")), flush=True)
except Exception:
    pass

from backend import (
    WAVELENGTHS,
    load_xrd_data,
    run_xrd_analysis,
    build_results_dataframe,
    create_results_zip,
)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
def _wavelength_from_selection(source: str, custom: float) -> float:
    if source == "Custom":
        try:
            w = float(custom)
            assert 0.01 < w < 10
            return w
        except Exception:
            return 1.5406
    return WAVELENGTHS.get(source, 1.5406)


def inspect_file(file):
    """Return a quick preview of the uploaded file so the user can confirm columns."""
    if file is None:
        return "No file uploaded."
    try:
        ext = os.path.splitext(file)[-1].lower()
        if ext in ('.xlsx', '.xls'):
            df = pd.read_excel(file, nrows=8)
        else:
            df = pd.read_csv(file, nrows=8, sep=None, engine='python')

        cols = list(df.columns)
        info = (
            f"✅  File loaded\n"
            f"Columns  : {cols}\n"
            f"Rows preview (first 8):\n"
            f"{df.to_string(index=False)}"
        )
        return info
    except Exception as exc:
        return f"⚠ Could not preview file: {exc}"


# ─────────────────────────────────────────────────────────────────────────────
# Main Analysis Callback
# ─────────────────────────────────────────────────────────────────────────────
def run_analysis(
    xrd_file,
    ref_file,
    material_name,
    source_choice,
    custom_wl,
    min_height_pct,
    min_prom_pct,
    min_dist_deg,
    scherrer_k,
    instr_fwhm,
    two_theta_min,
    two_theta_max,
    progress=gr.Progress(track_tqdm=True),
):
    """
    Central callback: runs the full XRD pipeline and returns all Gradio outputs.
    """
    # ── Validation ────────────────────────────────────────────────────────────
    if xrd_file is None:
        msg = "❌  Please upload an XRD data file before running analysis."
        return (None,)*9 + (msg, None, None)

    try:
        progress(0.05, desc="Reading XRD file …")
        wavelength = _wavelength_from_selection(source_choice, custom_wl)
        mat = (material_name or "Sample").strip() or "Sample"

        t_min = float(two_theta_min) if (two_theta_min not in (None, "")) else None
        t_max = float(two_theta_max) if (two_theta_max not in (None, "")) else None

        progress(0.12, desc="Detecting peaks …")
        results, figures, summary = run_xrd_analysis(
            xrd_file_path  = xrd_file,
            ref_file_path  = ref_file if ref_file else None,
            wavelength     = wavelength,
            material_name  = mat,
            min_ht_frac    = float(min_height_pct) / 100.0,
            min_prom_frac  = float(min_prom_pct)   / 100.0,
            min_dist_deg   = float(min_dist_deg),
            K              = float(scherrer_k),
            beta_instr_deg = float(instr_fwhm),
            two_theta_min  = t_min,
            two_theta_max  = t_max,
        )

        progress(0.70, desc="Generating plots …")
        n_peaks = len(results["peaks"])
        if n_peaks == 0:
            summary = ("⚠ No peaks detected with current settings.\n"
                       "Try lowering 'Minimum Peak Height' or 'Minimum Prominence'.\n\n"
                       + summary)

        progress(0.88, desc="Building results table …")
        df_peaks = build_results_dataframe(results)

        progress(0.95, desc="Packaging ZIP …")
        zip_path = create_results_zip(results, figures, summary)

        progress(1.00, desc="Done ✓")

        # Map figures (keys are like "01_raw_pattern")
        figs = figures
        return (
            figs.get("00_summary_dashboard"),
            figs.get("01_raw_pattern"),
            figs.get("02_peak_fits"),
            figs.get("03_scherrer"),
            figs.get("04_williamson_hall"),
            figs.get("05_size_strain_plot"),
            figs.get("06_microstrain"),
            figs.get("07_dislocation"),
            figs.get("08_reference"),
            summary,
            df_peaks,
            zip_path,
        )

    except Exception as exc:
        err = (
            f"❌ Analysis failed:\n{exc}\n\n"
            f"Traceback:\n{traceback.format_exc()}"
        )
        return (None,)*9 + (err, None, None)


# ─────────────────────────────────────────────────────────────────────────────
# Gradio UI
# ─────────────────────────────────────────────────────────────────────────────
_CSS = """
/* ── Global ────────────────────────────────────────────────────────────── */
body { font-family: 'Inter', sans-serif; }
.app-header { text-align:center; padding:24px 0 8px; }
.app-header h1 { font-size:2rem; font-weight:800;
                  background:linear-gradient(135deg,#1d4ed8,#7c3aed);
                  -webkit-background-clip:text; -webkit-text-fill-color:transparent; }
.app-header p  { color:#6b7280; font-size:0.95rem; margin-top:4px; }
.section-label { font-size:0.8rem; font-weight:700;
                  text-transform:uppercase; letter-spacing:.07em;
                  color:#6b7280; margin:16px 0 4px; }
.badge { display:inline-block; padding:2px 8px; border-radius:999px;
         font-size:0.7rem; font-weight:700; background:#dbeafe; color:#1e40af;
         margin:2px; }
.run-btn { background:linear-gradient(135deg,#1d4ed8,#7c3aed) !important;
           color:white !important; font-weight:700 !important;
           font-size:1rem !important; border-radius:8px !important; }
"""

SOURCE_CHOICES = ["Custom"] + list(WAVELENGTHS.keys())
SOURCE_DEFAULT = "Cu Kα (λ = 1.5406 Å)"

_METHODS_MD = """
### Analyses Performed

| # | Method | Formula | Reference |
|---|--------|---------|-----------|
| 1 | **Peak Detection** | scipy height/prominence/distance | scipy.signal |
| 2 | **Pseudo-Voigt Fit** | PV = A[η·L + (1−η)·G] + bg | Thompson et al. 1987 |
| 3 | **Integral Breadth** | β_I = Area / Height | — |
| 4 | **Bragg's Law** | d = λ / (2 sinθ) | Bragg 1913 |
| 5 | **Scherrer Eqn.** | D = Kλ / (β cosθ) | Scherrer 1918 / Patterson 1939 |
| 6 | **Williamson–Hall (UDM)** | β cosθ = Kλ/D + 4ε sinθ | Williamson & Hall 1953 |
| 7 | **Size–Strain Plot** | (d·β·cosθ)² = f(d²·β·cosθ) | Mote et al. 2012 |
| 8 | **Microstrain (S–W)** | ε = β / (4 tan θ) | Stokes & Wilson 1944 |
| 9 | **Dislocation Density** | δ = 1/D² | Williamson & Smallman 1956 |
|10 | **Crystallinity** | CI = area_cryst / area_total | Alexander 1969 |
|11 | **Phase Comparison** | d-spacing match | — |
|12 | **Texture Coefficient** | TC = (I/I₀) / ⟨I/I₀⟩ | Harris 1952 |

> **Instrument Broadening Correction** (optional):  
> β²_sample = β²_measured − β²_instrument (Warren formula)
"""

_FORMAT_MD = """
### Required Input Formats

**XRD Data File** (`.xlsx` / `.xls` / `.csv`)  
Two columns — headers are auto-detected:

| Column 1 | Column 2 |
|----------|----------|
| 2θ (°) | Intensity (a.u.) |

Column names may vary; the tool accepts common aliases such as `2theta`, `angle`, `counts`, `I`, etc.

---

**Reference Standard File** (`.csv` / `.xlsx`) — *optional*  
Two columns:

| Column 1 | Column 2 |
|----------|----------|
| 2θ (°) | Relative Intensity (%) |

Source: JCPDS/ICDD card, published pattern, or simulation.  
*Analysis runs completely without this file — it only enables phase comparison.*
"""


def build_app() -> gr.Blocks:
    with gr.Blocks(
        title="CharacterizationCalc — XRD",
        theme=gr.themes.Soft(
            primary_hue=gr.themes.colors.blue,
            secondary_hue=gr.themes.colors.violet,
            font=[gr.themes.GoogleFont("Inter"), "sans-serif"],
        ),
        css=_CSS,
    ) as demo:

        # ── Header ────────────────────────────────────────────────────────────
        gr.HTML("""
        <div class="app-header">
          <h1>🔬 CharacterizationCalc — XRD Analysis</h1>
          <p>Complete X-ray Diffraction analysis platform · Peak finding · Scherrer · Williamson-Hall · Strain · More</p>
          <div style="margin-top:8px;">
            <span class="badge">Bragg's Law</span>
            <span class="badge">Scherrer</span>
            <span class="badge">Williamson-Hall</span>
            <span class="badge">Size-Strain Plot</span>
            <span class="badge">Microstrain</span>
            <span class="badge">Dislocation Density</span>
            <span class="badge">Crystallinity</span>
            <span class="badge">Phase Comparison</span>
          </div>
        </div>
        """)

        with gr.Row(equal_height=False):
            # ── LEFT PANEL: Inputs ─────────────────────────────────────────
            with gr.Column(scale=1, min_width=360):

                # ── File Uploads ─────────────────────────────────────────────
                gr.Markdown("## 📁 Data Upload")

                with gr.Group():
                    gr.Markdown("### XRD Data File *(required)*")
                    gr.Markdown(
                        "Upload `.xlsx` / `.xls` / `.csv` with two columns:\n"
                        "- **Column 1** → `2θ (°)` (angle)\n"
                        "- **Column 2** → `Intensity (a.u.)`\n\n"
                        "*Column headers are auto-detected.*"
                    )
                    xrd_file = gr.File(
                        label="Upload XRD Data",
                        file_types=[".xlsx", ".xls", ".csv"],
                        type="filepath",
                    )
                    file_preview = gr.Textbox(
                        label="File Preview",
                        lines=6, interactive=False,
                        placeholder="Upload a file to see a preview …",
                    )
                    xrd_file.change(fn=inspect_file, inputs=xrd_file, outputs=file_preview)

                with gr.Group():
                    gr.Markdown("### Reference Standard File *(optional)*")
                    gr.Markdown(
                        "Upload a JCPDS/ICDD reference CSV for phase comparison.\n"
                        "- **Column 1** → `2θ (°)`\n"
                        "- **Column 2** → `Relative Intensity (%)`\n\n"
                        "*The full analysis runs even without this file.*"
                    )
                    ref_file = gr.File(
                        label="Upload Reference Standard (optional)",
                        file_types=[".csv", ".xlsx"],
                        type="filepath",
                    )

                # ── Sample Information ────────────────────────────────────────
                gr.Markdown("## ⚗️ Sample & Instrument")

                with gr.Group():
                    material_name = gr.Textbox(
                        label="Material / Sample Name",
                        placeholder="e.g. ZnO nanoparticles, TiO₂, CuO …",
                        value="Sample",
                    )

                    source_choice = gr.Dropdown(
                        choices=SOURCE_CHOICES,
                        value=SOURCE_DEFAULT,
                        label="X-ray Source",
                        info="Select radiation source; most labs use Cu Kα",
                    )
                    custom_wl = gr.Number(
                        value=1.5406, label="Custom Wavelength (Å)",
                        precision=5, visible=False,
                        info="Enter wavelength if 'Custom' selected",
                    )
                    source_choice.change(
                        fn=lambda s: gr.update(visible=(s == "Custom")),
                        inputs=source_choice, outputs=custom_wl,
                    )

                    instr_fwhm = gr.Number(
                        value=0.0, label="Instrument FWHM β_instr (°)",
                        precision=4,
                        info="FWHM from a standard (LaB₆/Si). 0 = no correction.",
                    )

                # ── Peak Detection ────────────────────────────────────────────
                gr.Markdown("## 🎯 Peak Detection")

                with gr.Group():
                    min_height_pct = gr.Slider(
                        1, 50, value=5, step=1,
                        label="Min Peak Height  (% of max intensity)",
                        info="Lower → detect smaller/weaker peaks",
                    )
                    min_prom_pct = gr.Slider(
                        1, 30, value=3, step=1,
                        label="Min Peak Prominence  (% of max intensity)",
                        info="Filters shoulder peaks — raise if too many peaks detected",
                    )
                    min_dist_deg = gr.Slider(
                        0.1, 10.0, value=0.5, step=0.1,
                        label="Min Peak Separation  (°)",
                        info="Minimum 2θ gap between adjacent peaks",
                    )

                # ── Scherrer ──────────────────────────────────────────────────
                gr.Markdown("## 📐 Scherrer Parameters")

                with gr.Group():
                    scherrer_k = gr.Number(
                        value=0.9, label="Shape Factor K",
                        precision=3,
                        info="0.90 general / 0.94 cubic / 1.00 perfect sphere",
                    )

                # ── 2θ Range ──────────────────────────────────────────────────
                gr.Markdown("## 📏 2θ Analysis Range *(optional)*")

                with gr.Group():
                    with gr.Row():
                        two_theta_min = gr.Number(
                            value=None, label="Min 2θ (°)  [e.g. 10]",
                            info="Leave blank to use full range",
                        )
                        two_theta_max = gr.Number(
                            value=None, label="Max 2θ (°)  [e.g. 80]",
                        )

                # ── Run Button ────────────────────────────────────────────────
                run_btn = gr.Button(
                    "🚀  Run Full XRD Analysis",
                    variant="primary",
                    elem_classes=["run-btn"],
                    size="lg",
                )

            # ── RIGHT PANEL: Results ───────────────────────────────────────
            with gr.Column(scale=2):
                gr.Markdown("## 📊 Analysis Results")

                with gr.Tabs():
                    with gr.TabItem("🗺 Summary Dashboard"):
                        plot_dashboard = gr.Plot(label="4-Panel Summary")

                    with gr.TabItem("📈 XRD Pattern"):
                        plot_raw = gr.Plot(label="Diffractogram + Detected Peaks")

                    with gr.TabItem("🔍 Peak Fitting"):
                        plot_fits = gr.Plot(label="Pseudo-Voigt Peak Fits")

                    with gr.TabItem("🔬 Scherrer Size"):
                        plot_scherrer = gr.Plot(label="Scherrer Crystallite Size")

                    with gr.TabItem("📐 Williamson-Hall"):
                        plot_wh = gr.Plot(label="Williamson–Hall Plot (UDM)")

                    with gr.TabItem("🔷 Size-Strain Plot"):
                        plot_ssp = gr.Plot(label="Size–Strain Plot (Mote 2012)")

                    with gr.TabItem("⚡ Microstrain"):
                        plot_ms = gr.Plot(label="Microstrain (Stokes–Wilson)")

                    with gr.TabItem("🔩 Dislocation Density"):
                        plot_dd = gr.Plot(label="Dislocation Density")

                    with gr.TabItem("📌 Reference"):
                        plot_ref = gr.Plot(label="Measured vs Reference Pattern")

                gr.Markdown("### 📋 Analysis Summary Report")
                summary_box = gr.Textbox(
                    label="Report", lines=28, interactive=False,
                    placeholder="Run analysis to see the full report …",
                )

                gr.Markdown("### 📊 Per-Peak Results Table")
                results_df = gr.DataFrame(
                    label="Peak-by-Peak Analysis",
                    wrap=True,
                )

                gr.Markdown("### ⬇️ Download All Results")
                download_file = gr.File(
                    label="Download ZIP (plots + Excel + CSV + JSON + report)",
                    interactive=False,
                )

        # ── Accordion: Methods & Formats ──────────────────────────────────────
        with gr.Accordion("📚 Methods, Formulas & References", open=False):
            gr.Markdown(_METHODS_MD)

        with gr.Accordion("📋 File Format Guide", open=False):
            gr.Markdown(_FORMAT_MD)

        with gr.Accordion("💡 Interpretation Guide", open=False):
            gr.Markdown("""
### How to interpret the results

**Scherrer Crystallite Size (D)**
- Gives the *volume-averaged* coherently scattering domain size
- Smaller D → more nanocrystalline / poorly crystallised
- *Caution*: Scherrer conflates size AND strain — use WH if strain is present

**Williamson–Hall (UDM)**
- Separates size (intercept) from microstrain (slope)
- Positive slope → **tensile** microstrain; negative → **compressive**
- If R² < 0.7, results should be treated cautiously (few peaks or noise)
- Non-positive intercept means the WH size is unphysical; use Scherrer instead

**Size–Strain Plot (SSP)**
- Independent alternative to WH; gives area-weighted size + rms strain
- Cross-check with WH — agreement validates the analysis

**Microstrain (Stokes–Wilson, ε = β/4tanθ)**
- Per-peak estimate; scatter is normal for experimental data
- Average ε from WH is more statistically robust

**Dislocation Density (δ = 1/D²)**
- Higher crystallinity / larger D → lower δ
- Typical annealed metals: 10¹⁰–10¹² m⁻²; heavily deformed: 10¹⁴–10¹⁶ m⁻²

**Crystallinity Index**
- Sensitive to baseline method; use as *relative* comparison between samples
- Below ~40 % → predominantly amorphous material

**Phase Match Quality**
- 100 % → every measured peak matches the reference
- >70 % generally indicates the correct phase
""")

        # ── Wire up the Run button ─────────────────────────────────────────────
        run_btn.click(
            fn=run_analysis,
            inputs=[
                xrd_file, ref_file, material_name,
                source_choice, custom_wl,
                min_height_pct, min_prom_pct, min_dist_deg,
                scherrer_k, instr_fwhm,
                two_theta_min, two_theta_max,
            ],
            outputs=[
                plot_dashboard,
                plot_raw,
                plot_fits,
                plot_scherrer,
                plot_wh,
                plot_ssp,
                plot_ms,
                plot_dd,
                plot_ref,
                summary_box,
                results_df,
                download_file,
            ],
        )

    return demo


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app = build_app()
    app.launch(
        server_name="0.0.0.0",
        server_port=int(os.environ.get("PORT", 7860)),
        share=False,
        show_error=True,
        favicon_path=None,
    )
