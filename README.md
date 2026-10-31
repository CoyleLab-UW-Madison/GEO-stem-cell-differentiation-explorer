# GEO Stem Cell Differentiation Explorer

Interactive field-of-view (FOV) atlas for genetically encoded oscillator (GEO) microscopy during human embryonic stem cell (H9 hESC) trilineage differentiation:

**[https://coylelab-uw-madison.github.io/GEO-stem-cell-differentiation-explorer/](https://coylelab-uw-madison.github.io/GEO-stem-cell-differentiation-explorer/)**

---

## Dataset

Time-lapse dual-color fluorescence microscopy (`mCherry-MinD` and `MinE-GFP`) across **420 fields of view** ($N = 448,054$ segmented oscillating cells):
* **Lineages (4):** Pluripotent hESC Control (`mTeSR`), Ectoderm, Endoderm, and Mesoderm.
* **Timepoints (7):** Sampled every 12 hours following singularization (`Day 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0`).
* **Replicates:** `3` biological replicate wells per condition $\times$ `5` fields of view per well (`15` FOVs per condition $\times$ `28` conditions = `420` FOVs).

---

## Spatial Views & Decoded Physiology Metrics

Each field of view is rendered as a 5-panel horizontal strip (`fov_atlas/*.png`) with pixel-exact registration across all channels:

| Panel | Metric | Scale & Colormap | Description |
| :---: | :--- | :--- | :--- |
| **1** | **Fluorescence Still** | Grayscale (`1st–99.2th` pct) | Time-averaged `MinE-GFP` fluorescence showing colony morphology across the full expressing tissue footprint. |
| **2** | **FM Barcode ($f$)** | `18 – 52 mHz` (`Turbo`) | Raw single-cell oscillation frequency determined by phase-vocoder demodulation, reflecting clonal `(MinD, MinE)` expression stoichiometry. |
| **3** | **Frequency Shift ($\Delta F_{\text{stem}}$)** | `-12 to +12 mHz` (`Blue ↔ White ↔ Red`) | Stoichiometry-demodulated carrier frequency shift $\Delta F_{\text{stem}} = f - f_{\text{ref}}(D, E)$ relative to the empirical pluripotent hESC reference surface $f_{\text{ref}}(D, E)$. |
| **4** | **Competence Shift ($\Delta P_{\text{osc}}$)** | `-50% to +50%` (`Teal ↔ Grey ↔ Gold`) | Shift in local oscillation competence $\Delta P_{\text{osc}} = P_{\text{osc}} - P_{\text{ref}}(D, E)$ relative to the pluripotent hESC competence surface $P_{\text{ref}}(D, E)$, with neutral silver-grey spanning the control interquartile range (`±13%`). |
| **5** | **Stoichiometric Tilt ($\Delta s_\perp$)** | `-1.0 to +2.8 mHz` (`Grey → Purple`) | Rotation of the local $(D, E)$ frequency gradient relative to the pluripotent hESC reference surface, with neutral silver-grey below the empirical hESC 95th-percentile ROC cutoff (`≤ +0.85 mHz`) and purple highlighting the supra-threshold tail (`> +0.85 mHz`). |

---

## Repository Files

* `index.html`: Self-contained interactive explorer for browsing all 420 FOVs by lineage, day, well, or replicate montage, with single-panel zoom and two-channel overlay blending.
* `fov_atlas/*.png`: 5-panel spatial strips for all 420 FOVs.
* `figures/`: Summary figures (`Figure5_GEO_Physiology_Encoding_Endpoints.png` and `Figure6_GEO_Physiology_Timecourse_Trajectories.png`).
* `2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv` / `.json`: Per-FOV summary statistics across all 420 FOVs.
* `build_h9_fov_atlas_and_html_viewer.py`: Script used to render the 5-panel FOV strips and generate `index.html`.
