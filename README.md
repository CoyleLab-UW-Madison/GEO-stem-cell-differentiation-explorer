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

Each field of view is rendered as a 7-panel horizontal strip (`fov_atlas/*.png`) with pixel-exact registration across all channels:

| Panel | Metric | Scale & Colormap | Description |
| :---: | :--- | :--- | :--- |
| **1** | **Fluorescence Still** | Grayscale (`2nd–99.2nd` pct) | Time-averaged `MinE-GFP` fluorescence showing colony morphology across the full expressing tissue footprint. |
| **2** | **FM Barcode ($f$)** | `18 – 52 mHz` (`Turbo`) | Raw single-cell oscillation frequency determined by phase-vocoder demodulation, reflecting clonal `(MinD, MinE)` expression stoichiometry. |
| **3** | **Frequency Shift ($\Delta F_{\text{stem}}$)** | `-12 to +12 mHz` (`Blue ↔ White ↔ Red`) | Stoichiometry-demodulated carrier frequency shift $\Delta F_{\text{stem}} = f - f_{\text{ref}}(D, E)$ relative to the empirical pluripotent hESC reference surface $f_{\text{ref}}(D, E)$. |
| **4** | **Competence Shift ($\Delta P_{\text{osc}}$)** | `-50% to +50%` (`Teal ↔ Grey ↔ Gold`) | Shift in local oscillation competence $\Delta P_{\text{osc}} = P_{\text{osc}} - P_{\text{ref}}(D, E)$ relative to the pluripotent hESC competence surface $P_{\text{ref}}(D, E)$, with neutral silver-grey spanning the control interquartile range (`±13%`). |
| **5** | **Stoichiometric Tilt ($\Delta s_\perp$)** | `-1.0 to +2.8 mHz` (`Grey → Purple`) | Rotation of the local $(D, E)$ frequency gradient relative to the pluripotent hESC reference surface, with neutral silver-grey below the empirical hESC 95th-percentile ROC cutoff (`≤ +0.85 mHz`) and purple highlighting the supra-threshold tail (`> +0.85 mHz`). |
| **6** | **Aclonal Coherence ($\mathcal{C}_{\text{nc}}$)** | `0.00 to 0.58` (`TV-Static → Hot-Pink`) | Mesoscale physiological coherence ($r \le 45\;\mu\text{m}$) restricted strictly to **non-clonal** neighbor pairs ($|\Delta \log_2(E/D)| > 0.38$), rendering uncoupled/incoherent cells as TV-static white noise and multiclonal synchronized territories in solid hot-pink/magenta. |
| **7** | **Aclonal Scale ($r_{1/2}$)** | `8 to 58 µm` (`TV-Static → Spring-Emerald`) | Single-cell radial fall-off correlation lengthscale integrated across expanding concentric non-clonal shells ($22, 38, 56, 78\;\mu\text{m}$), highlighting interior multiclonal domain cores in vivid spring-emerald and carving crisp TV-static boundaries between adjacent opposing domains. |

---

## Repository Files

* `index.html`: Self-contained interactive explorer for browsing all 420 FOVs by lineage, day, well, or replicate montage, with single-panel zoom and two-channel overlay blending.
* `fov_atlas/*.png`: 7-panel spatial strips for all 420 FOVs.
* `figures/`: Summary figures (`Figure5_GEO_Physiology_Encoding_Endpoints.png`, `Figure6_GEO_Physiology_Timecourse_Trajectories.png`, `Diagnostic_NonClonal_Radial_Decay_and_SingleCell_Coherence_Maps.png`, and `Diagnostic_3State_Coherence_Composite_Dissection.png`).
* `2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv` / `.json`: Per-FOV summary statistics across all 420 FOVs.
* `build_h9_fov_atlas_and_html_viewer.py`: Script used to render the 7-panel FOV strips and generate `index.html`.
