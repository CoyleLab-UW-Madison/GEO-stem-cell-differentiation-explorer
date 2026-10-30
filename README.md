# GEO Stem Cell Differentiation Explorer

**Interactive Phase-Vocoder Frequency Modulation (FM), Decoded Physiological State ($\Delta F_{\text{stem}}$), Colony-Edge Zonation ($d_{\text{edge}}$), ROC Mesoderm/EMT ($\Delta s_\perp^{\text{stem}}$), and Stoichiometry-Matched Oscillation Probability ($\Delta P_{\text{osc}}^{\text{stem}}$) Atlas across 420 Fields of View ($N = 448,054$ Oscillating Single Cells)**

🔗 **Live Interactive Explorer (GitHub Pages):**  
**[https://coylelab-uw-madison.github.io/GEO-stem-cell-differentiation-explorer/](https://coylelab-uw-madison.github.io/GEO-stem-cell-differentiation-explorer/)**

---

## Overview

This repository hosts the interactive 420-FOV spatial and single-cell explorer for genetically encoded MinDE oscillation microscopy during H9 human embryonic stem cell (hESC) trilineage differentiation (`Control`, `Ectoderm`, `Endoderm`, `Mesoderm` sampled every `12 hours` from `Day 0.0` to `Day 3.0`, with `3 biological replicate wells × 5 FOVs/well = 420 FOVs`).

Each FOV is rendered as a pixel-exact **12-Panel (`2 × 6`) Tile Sheet** (`fov_atlas/*.png`, web-optimized `2880 × 1050 px` on GitHub Pages and full-resolution `4320 × 1575 px` on the Research Drive fileserver) where all 11 spatial channels share an identical relative bounding box (`0.00 px` registration offset). This enables:
* **Full 12-Panel (`2 × 6`) Inspection** or single-channel zoom (`Tiles 1–12`),
* **Interactive Side-by-Side & Pixel-Exact Live Blend / Overlay** of any two channels (`Base A` + `Overlay B`) with a live **Opacity Slider (`0–100%`)**, **CSS Blend Formula Selector** (`screen`, `normal`, `hard-light`, `overlay`, `lighten`), **Curtain-Wipe Comparison Slider (`0–100%`)**, and optional **Emerald Green Tint** on Channel 1,
* **Browser-Persisted FOV Bookmarking & Research Notes (`localStorage` + 1-Click CSV Export)** for tagging and exporting notes across all 420 FOVs.

---

## Dataset Funnel & Physical Architecture

| Dataset Tier | Count Across 420 FOVs | Median Area (`px`) | Role in Biophysical Pipeline |
| :--- | ---: | ---: | :--- |
| **Full Expressing Tissue (`in_tissue`)** | `420` FOV maps (`E ≥ 560 \| D ≥ 188 \| mask > 0`) | Continuous tissue | **Primary domain for Axis 2 ($P_{\text{ctrl}}, \Delta P_{\text{osc}}$, Panels `7`, `8A–8C`, `9`)**, capturing both segmented cells and unsegmented (`mask == 0`) silent/quenched colony regions. |
| **Raw Amplitude Masks (`mask > 0`)** | `1,817,010` objects | `9 px` | Raw Cellpose output on the oscillation amplitude projection (includes `959,046` `0`-osc-pixel fragments of median area `6 px`). |
| **Real-Sized Segmented Masks (`area ≥ 100 px`)** | `285,421` masks | `198 px` | Cell-sized amplitude masks reported as `Real Masks ≥100px` in the interactive explorer. |
| **Masks with $\ge 1$ Strict Osc Pixel (`df1`)** | `857,964` objects | `24 px` | `2026_10_20_stemcell_all_single_cells_vocoder_physiology.parquet` (includes `409,910` small fragments with `1–9` oscillating pixels). |
| **Oscillating Single Cells (`n_osc_pixels ≥ 10`, `df2`)** | **`448,054` cells** | **`175 px`** | `2026_10_24_single_cell_posterior_confidence_audit.parquet`. **Primary single-cell cohort for Axis 1 ($\Delta F_{\text{stem}}, \Delta s_\parallel, \Delta s_\perp$, and radial colony-edge zonation).** |
| **High-Confidence Sail Core (`g_strict`)** | **`131,806` cells** | **`265 px`** | `n_osc_pixels ≥ 35`, `in_pivot_frac ≤ 0.15`, `cnt0_ge50_frac ≥ 0.80`, `chi2_min ≤ 0.85`, `px_nopr_voc_std ≤ 0.20`. **Used to calibrate the $(\Delta s_\parallel^{\text{stem}}, \Delta s_\perp^{\text{stem}})$ Sail manifold.** |

---

## The 12 Panels Bundled per Field of View (`2 × 6` Grid)

| Tile | Panel Name | Biophysical Quantity |
| :---: | :--- | :--- |
| **1** | **1. Grayscale Still (`*_g.tif`)** | Time-averaged `MinE-GFP` fluorescence (`1st–99.5th` percentile clean grayscale LUT) showing the full colony morphology (`100%` of expressing tissue) for effortless visual inspection and live blending. |
| **2** | **2. Dual Still (`MinD` + `MinE`)** | Composite `Red = mCherry-MinD`, `Green = MinE-GFP` time-averaged fluorescence still with subtle colony boundaries. |
| **3** | **3. Stoichiometry `E/D`** | Dimensionless carrier ratio `MinE / MinD` (`[2.0, 5.6]`) across the full colony footprint (where cell thickness / optical path length $h_{\text{cell}}$ cancels out identically: $\frac{[E]\cdot h_{\text{cell}}}{[D]\cdot h_{\text{cell}}} = \frac{[E]}{[D]}$). |
| **4** | **4. Raw FM Barcode `f_mHz`** | Raw phase-vocoder oscillation frequency (`18.5–49.5 mHz`) showing clonal salt-and-pepper heterogeneity driven by carrier `E/D`. |
| **5** | **5. Decoded State $\Delta F_{\text{stem}}$** | Non-parametric stoichiometry-decoded frequency shift vs. homeostatic `Control` (`-10.5` to `+10.5 mHz`, grey steady-state band `[-1.5, +1.5] mHz`)—revealing radial colony-edge zonation and multi-cell physiological domains that transcend clonal `E/D` barcodes. |
| **6** | **6. ROC Meso/EMT $\Delta s_\perp^{\text{stem}}$** | ROC-calibrated orthogonal departure from the differentiation highway ($\tau^* = +0.85\text{ mHz}$ grey-to-pink cutoff, `[-1.0, +2.8] mHz`). |
| **7** | **7. Expected Control $P_{\text{ctrl}}(D, E)$** | Expected oscillation probability (`0–82%`, `viridis`) predicted from per-pixel $(D, E)$ expression using the homeostatic `Control (d1.5–d3.0)` surface. |
| **8** | **8A. Surface Shift $\Delta P_{\text{osc}}^{\text{surf}}(D, E)$** | Population-level surface shift $P_{\text{cond}}(D, E) - P_{\text{ctrl}}(D, E)$ (`[-65%, +65%]`) mapped onto each pixel's $(D, E)$ coordinates. |
| **9** | **8B. Local Gated $\Delta P_{\text{osc}}^{\text{stem}}(y, x)$** | Smoothed local oscillation fraction minus expected Control competence $P_{\text{local}}(\sigma=6\text{ px}) - P_{\text{ctrl}}(D, E)$ (`Blue / Grey / Gold ±12%`), gated to competent/active tissue. |
| **10** | **8C. Dual-ROC Surprise (`Blue / Grey / Pink`)** | Calibrated surprise map (`[-100, +100]`) keeping `~80–91%` of `Control` tissue **Flat Slate Grey** (`|S_roc| < 28.5`) while separating **Competent-Domain Quenching (`0-on-1` Blue, $P_{\text{dom}} - P_{\text{ctrl}} \le -22\%$)** in Ectoderm/Endoderm interiors and **Low-Threshold Active Cells (`1-on-0` Pink, $P_{\text{ctrl}} \le 32.5\%$)** at colony edges and in Mesoderm EMT. |
| **11** | **9. 5-Class Dead-Zone Map** | Pixel-level classification across the full expressing colony: **Pearl-Grey** (`Strict-Z Oscillating`), **Gold Class 1** (`Relaxed-Z Low-Amplitude Oscillators, 4.0 ≤ D_z < 5.0`), **Cyan Class 2** (`High E/D > 4.6 Over-Driven`), **Violet Class 3** (`Balanced-E/D Damped Interior`), **Coral Class 4** (`Low E/D < 2.5 Arrested`), and **Crimson Class 5** (`Dim Reporter Tissue`). |
| **12** | **10. Cells on 2D Stem-Cell Sail** | Single cells plotted on the $(\Delta s_\parallel^{\text{stem}}, \Delta s_\perp^{\text{stem}})$ Sail manifold with the homeostatic `Control` bullseye `(0, 0)` and `ROC Meso+` threshold ($\Delta s_\perp^{\text{stem}} = +0.85\text{ mHz}$). |

---

## Key Summary Figures

1. **[Figure 50 — Out-of-Sample Reference Robustness (`LOWO`, Day-Matched, Edge-Matched $P_{\text{ctrl}}(D, E \mid d_{\text{edge}})$), Knot-Free Domain LLR, & Competent-Mask Audit](figures/Fig50_Reference_Robustness_EdgeMatched_DayMatched_LOMO_Audit.png)**
2. **[Figure 49 — Colony-Edge Physiological Zonation (`24/24` Wells) & Edge-Conditioned Non-Clonal Neighbor Wave Coherence (`Δφ` Phase Independence & `±1 px` Border Step)](figures/Fig49_Colony_Edge_Zonation_and_Edge_Conditioned_Neighbor_Coherence.png)**
3. **[Figure 48 — Master Representative Contact Sheet: All Spatial Views & Candidate $\Delta P_{\text{osc}}$ Formulations Across All 4 Lineages (`Day 3.0`)](figures/Fig48_Representative_Contact_Sheet_All_Spatial_Views_4Lineages.png)**
4. **[Figure 47 — Dual-ROC Calibrated Oscillation Surprise Map: Control Baseline (`1-on-1` Grey) vs. Quenched (`0-on-1` Blue) & Low-Threshold Active (`1-on-0` Pink)](figures/Fig47_Dual_ROC_Surprise_Flat_Control_Blue_Quench_Pink_Meso.png)**
5. **[Figure 46 — Marrying Bulk Trajectories with the Fit-Free Dual Sail: $(D, E)$ Phase-Portrait Bins in $(\Delta F_{\text{stem}}, \Delta P_{\text{osc}}^{\text{stem}})$](figures/Fig46_StoichiometryBin_FitFree_Dual_Sail_dF_vs_dPosc.png)**
6. **[Figure 45 — Fit-Free Dual Encoding & Excitability Space: Stoichiometry-Matched $\Delta F_{\text{stem}}$ vs. $\Delta P_{\text{osc}}^{\text{stem}}$ in 2D, 3D, and $4\times 7$ Landscapes](figures/Fig45_FitFree_Dual_Space_dFstem_vs_dPosc_2D_3D_and_Sail.png)**
7. **[Figure 44 — Control-Competent Gated $\Delta P_{\text{osc}}^{\text{stem}}$ & Coupling to ROC Meso/EMT ($\Delta s_\perp^{\text{stem}}$)](figures/Fig44_ControlCompetent_Gated_dPosc_and_Meso_Coupling.png)**
8. **[Figure 43 — Stoichiometry-Matched Oscillation Probability Surface $P_{\text{osc}}^{\text{Control}}(D, E)$ & Lineage $\Delta P_{\text{osc}}^{\text{stem}}$](figures/Fig43_Stoichiometry_Matched_Oscillation_Probability.png)**
9. **[Figure 42 — Transcending the Clonal FM Barcode & Diagnosing Colony Dead Zones in $(D, E)$ Expression Space](figures/Fig42_Barcode_Transcendence_and_DeadZone_Diagnostics.png)**
10. **[Figure 41 — 4-Channel FOV Decoding: Green Still, Raw FM Barcode, Decoded State Shift, and ROC Meso/EMT Axis](figures/Fig41_4Channel_GreenStill_FMBarcode_StateShift_ROCMeso.png)**
11. **[Figure 40 — Dataset-Wide (419 FOVs) Spatial Wave Coherence & Patchy Domain Remodeling](figures/Fig40_Spatial_Wave_Coherence_419FOVs.png)**
12. **[Figure 34 — Global 2D Stem-Cell-Centered Sail Manifold Trajectories](figures/Fig34_Global_2D_Sail_Trajectories.png)**

---

## Repository Structure

* `index.html` / `H9_FOV_Physiology_Atlas_Viewer.html`: Self-contained interactive 420-FOV 12-panel & live-blend explorer with browser-persisted FOV bookmarking (`localStorage` + CSV export).
* `fov_atlas/*.png`: All 420 12-panel (`2 × 6`) FOV tile sheets (`2880 × 1050 px` web-optimized; `4320 × 1575 px` full-resolution on the Research Drive fileserver).
* `figures/*.png`: Summary figures (`Figs 34, 40–50`).
* `2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv` & `.json`: Per-FOV quantitative metrics across all 420 FOVs (including `n_cells`, `n_seg_100`, `df_edge_med`, `df_int_med`, and `8A–8C` metrics).
* `FRESH_EYES.md`, `FRESH_EYES_REVIEWED.md`, `FRESH_EYES_CONSENSUS.md`, `Fresh_Eyes_AgreedUpon.md`: Independent methodological review, empirical mask-funnel audit, and consensus execution plan.
* `build_h9_fov_atlas_and_html_viewer.py`: Python build pipeline that computes the reference surfaces, renders all 420 12-panel FOV sheets, and compiles the interactive HTML viewer (`--html-only` for fast HTML/metadata rebuilds).
