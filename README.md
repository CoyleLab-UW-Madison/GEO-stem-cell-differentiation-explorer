# GEO Stem Cell Differentiation Explorer

**Interactive Phase-Vocoder Frequency Modulation (FM), Decoded Physiological State ($\Delta F_{\text{stem}}$), ROC Mesoderm/EMT ($\Delta s_\perp^{\text{stem}}$), and Stoichiometry-Matched Oscillation Probability ($\Delta P_{\text{osc}}^{\text{stem}}$) Atlas across 420 Fields of View ($N = 857,964$ Oscillating Cells)**

🔗 **Live Interactive Explorer (GitHub Pages):**  
**[https://coylelab-uw-madison.github.io/GEO-stem-cell-differentiation-explorer/](https://coylelab-uw-madison.github.io/GEO-stem-cell-differentiation-explorer/)**

---

## Overview

This repository hosts the interactive 420-FOV spatial and single-cell explorer for genetically encoded MinDE oscillation microscopy during H9 human embryonic stem cell (hESC) trilineage differentiation (`Control`, `Ectoderm`, `Endoderm`, `Mesoderm` sampled every `12 hours` from `Day 0.0` to `Day 3.0`, with `3 replicate wells × 5 FOVs/well = 420 FOVs`).

Each FOV is rendered as a pixel-exact **10-Panel (`2 × 5`) Tile Sheet** (`fov_atlas/*.png`) that can be inspected as a full 10-panel sheet, cropped dynamically to any single channel (`1–10`), or compared side-by-side in any pair of channels across single FOVs, 7-day timecourses, 4-lineage comparisons, or 15-FOV condition montages.

### The 10 Channels Bundled per Field of View

| Panel | Channel Name | Biophysical Quantity |
| :---: | :--- | :--- |
| **1** | **Green Still (`*_g.tif`)** | Time-averaged `MinE-GFP` fluorescence (`1st–99.5th` percentile emerald LUT) showing the entire colony footprint (both oscillating and silent tissue). |
| **2** | **Dual Expression Still (`MinD` + `MinE`)** | Composite `Red = MinD`, `Green = MinE` time-averaged still with white boundaries overlaid around oscillating cell segments. |
| **3** | **Stoichiometry `E/D` (Full Colony)** | Carrier expression ratio `MinE / MinD` evaluated across the full expressing colony footprint. |
| **4** | **Raw Frequency `f_mHz` (`The FM Barcode`)** | Raw phase-vocoder oscillation frequency (`18–50 mHz`, `turbo` LUT)—showing salt-and-pepper heterogeneity driven by clonal carrier `E/D` differences. |
| **5** | **Decoded State Shift $\Delta F_{\text{stem}}^{\text{surface}}$** | Non-parametric stoichiometry-matched frequency shift vs. steady-state `Control` (`-10.5` to `+10.5 mHz`, grey steady-state band `[-1.5, +1.5] mHz`)—revealing coherent multi-cell physiological domains that transcend the raw carrier barcode. |
| **6** | **ROC Meso/EMT Axis $\Delta s_\perp^{\text{stem}}$** | ROC-calibrated orthogonal departure (`grey` $< +0.85\text{ mHz} \to$ `pink/crimson` $\ge +0.85\text{ mHz}$) isolating Mesoderm/EMT trajectory cells ($4.7\times$ odds enrichment). |
| **7** | **Expected Control $P_{\text{osc}}^{\text{Ctrl}}(D, E)$** | Expected oscillation probability (`0–82%`, `viridis`) predicted purely from each pixel's $(D, E)$ expression using the steady-state `Control (d1.5–3.0)` competence surface. |
| **8** | **Control-Competent Gated $\Delta P_{\text{osc}}^{\text{stem}}$** | Stoichiometry-matched oscillation probability shift $\Delta P_{\text{osc}}^{\text{stem}} = P_{\text{local}} - P_{\text{ctrl}}(D, E)$ gated to **Control-Competent territory ($P_{\text{ctrl}} \ge 20\%$)**. Non-competent `"0-on-0"` regions ($P_{\text{ctrl}} < 20\%$) are masked dark charcoal; **Deep Royal Blue** (`<-15%`) marks intrinsic biophysical quenching (`Ectoderm`/`Endoderm`), and **Orange/Gold** (`>+15%`) marks `"More Alive Than Control"` hyper-excitability (`Mesoderm/EMT`). |
| **9** | **5-Class Dead-Zone Diagnostic Map** | Classifies every expressing colony pixel into **White** (`Oscillating`), **Bright Cyan** (`High E/D > 5.0` saturation dead zone), **Gold** (`Optimal E/D & D`, likely `Z ≥ 5` filter drop), **Crimson Red** (`Dim D < 210` sub-threshold dead zone), and **Purple** (`Low E/D < 2.2`). |
| **10** | **Cells on 2D Stem-Cell Sail** | Single cells plotted on the $(\Delta s_\parallel^{\text{stem}}, \Delta s_\perp^{\text{stem}})$ Sail and colored by their cell-level $\Delta P_{\text{osc}}^{\text{cell}}$ (`Blue` $\to$ `Grey` $\to$ `Orange/Gold`), showing how Orange/Gold `+\Delta P_{\text{osc}}` cells sit directly above the dashed pink `ROC Meso+` threshold ($\Delta s_\perp^{\text{stem}} = +0.85\text{ mHz}$). |

---

## Key Summary Figures

1. **[Figure 44 — Control-Competent Gated $\Delta P_{\text{osc}}^{\text{stem}}$ & Coupling of `+ΔP_osc` to ROC Meso/EMT](figures/Fig44_ControlCompetent_Gated_dPosc_and_Meso_Coupling.png)**
2. **[Figure 43 — Stoichiometry-Matched Oscillation Probability Surface $P_{\text{osc}}^{\text{Control}}(D, E)$ & Lineage $\Delta P_{\text{osc}}^{\text{stem}}$](figures/Fig43_Stoichiometry_Matched_Oscillation_Probability.png)**
3. **[Figure 42 — Transcending the Clonal FM Barcode & Diagnosing Colony Dead Zones in $(D, E)$ Expression Space](figures/Fig42_Barcode_Transcendence_and_DeadZone_Diagnostics.png)**
4. **[Figure 41 — 4-Channel FOV Decoding: Green Still, Raw FM Barcode, Decoded State Shift, and ROC Meso/EMT Axis](figures/Fig41_4Channel_GreenStill_FMBarcode_StateShift_ROCMeso.png)**
5. **[Figure 40 — Dataset-Wide (419 FOVs) Spatial Wave Coherence & Patchy Domain Remodeling](figures/Fig40_Spatial_Wave_Coherence_419FOVs.png)**
6. **[Figure 34 — Global 2D Stem-Cell-Centered Sail Manifold Trajectories](figures/Fig34_Global_2D_Sail_Trajectories.png)**

---

## Repository Structure

* `index.html` / `H9_FOV_Physiology_Atlas_Viewer.html`: Self-contained interactive 420-FOV explorer.
* `fov_atlas/*.png`: All 420 10-panel (`2 × 5`) FOV tile sheets (`2625 × 1102 px`).
* `figures/*.png`: High-resolution summary figures (`Figs 34, 40–44`).
* `2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv` & `.json`: Per-FOV quantitative metrics across all 420 FOVs.
* `build_h9_fov_atlas_and_html_viewer.py`: Python build pipeline that computes the 2D `Control` reference surfaces, renders all 420 10-panel FOV sheets, and compiles the interactive HTML viewer.
