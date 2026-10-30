# Fresh_Eyes_AgreedUpon.md — Final Consensus & Execution Plan (Opus 5.5 × Argon)

*Prepared 2026-10-29 after synthesizing [`FRESH_EYES.md`](file:///home/smcoyle/.gemini/antigravity-cli/brain/aba364ca-b51f-459e-82e0-cad8dfe277e5/FRESH_EYES.md), [`FRESH_EYES_REVIEWED.md`](file:///home/smcoyle/.gemini/antigravity-cli/brain/aba364ca-b51f-459e-82e0-cad8dfe277e5/FRESH_EYES_REVIEWED.md), and [`FRESH_EYES_CONSENSUS.md`](file:///home/smcoyle/.gemini/antigravity-cli/brain/aba364ca-b51f-459e-82e0-cad8dfe277e5/FRESH_EYES_CONSENSUS.md).*

---

## 1. Locked-In Scientific & Methodological Consensus

### 1.1 The Exact Dataset Funnel (420 FOVs: 4 Lineages × 7 Days × 3 Wells × 5 FOVs/Well)
Opus's empirical check ([`consensus_mask_funnel_check.py`](file:///home/smcoyle/h9_vocoder_work/workspace_archive/consensus_mask_funnel_check.py)) and our full-parquet verification settle the exact meaning of every tier in the dataset:

| Tier | Count Across 420 FOVs | Median Area (`px`) | Physical Meaning & Role in Pipeline |
|---|---:|---:|---|
| **Full Expressing Tissue (`in_tissue`)** | `420` FOV maps (`E ≥ 560 \| D ≥ 188 \| mask > 0`) | Continuous tissue | **Primary domain for Axis 2 ($P_{\text{osc}}, \Delta P_{\text{osc}}$, Panels `7`, `8A–8C`, `9`)**, because silent/quenched cells lack oscillation amplitude and mostly sit in `mask == 0`. |
| **All Amplitude Masks (`mask > 0`)** | `1,817,010` objects | `9 px` | Raw Cellpose output on the amplitude projection. Includes `959,046` (`52.8%`) `0`-osc-pixel fragments (`median 6 px`, `98.8% < 100 px`). **Do not use raw `n_seg_total` as a cell denominator.** |
| **Masks with $\ge 1$ Strict Osc Pixel (`df1`)** | `857,964` objects | `24 px` | `2026_10_20_stemcell_all_single_cells_vocoder_physiology.parquet`. Includes `409,910` small fragments with `1–9` osc pixels (`median 11 px`, `96.4% < 100 px`). |
| **Oscillating Single Cells (`n_osc_pixels ≥ 10`, `df2`)** | **`448,054` cells** | **`175 px`** | `2026_10_24_single_cell_posterior_confidence_audit.parquet`. **Primary population for Axis 1 ($\Delta F_{\text{stem}}$ and colony-edge vs. interior spatial physiology).** |
| **High-Confidence Sail Core (`g_strict`)** | **`131,806` cells** | **`265 px`** | `n_osc_pixels ≥ 35`, `in_pivot_frac ≤ 0.15`, `cnt0_ge50_frac ≥ 0.80`, `chi2_min ≤ 0.85`, `px_nopr_voc_std ≤ 0.20`. **Used to calibrate the $(\Delta s_\parallel, \Delta s_\perp)$ Sail manifold.** |

> [!IMPORTANT]
> **Why Axis 1 Is Cell-Level and Axis 2 Is Pixel/Domain-Level:**
> Because Cellpose segmentation was generated from the **oscillation amplitude image**, quiescent cells in quenched Ectoderm/Endoderm domains rarely receive a segmentation mask (`mask == 0`). Therefore, **Axis 1 ($\Delta F_{\text{stem}}, \Delta s_\parallel, \Delta s_\perp$)** is evaluated on the `448,054` oscillating cells (`n_osc_pixels ≥ 10`), whereas **Axis 2 ($P_{\text{ctrl}}, \Delta P_{\text{osc}}$, `8A–8C`)** is evaluated on the full expressing tissue footprint (`in_tissue`).

---

### 1.2 Agreed Physical & Biological Principles
1. **3D ER Waves & Optical Path-Length ($h_{\text{cell}}$) Cancellation:**
   * Mammalian MinDE waves run on the ER membrane as a 3D volume-filling oscillation whose frequency is insensitive to cell shape.
   * In widefield imaging, per-pixel intensity scales with $[\text{concentration}] \times h_{\text{cell}}$. In the ratio $E/D = \frac{[E] \cdot h_{\text{cell}}}{[D] \cdot h_{\text{cell}}}$, cell height $h_{\text{cell}}$ **cancels identically**, making $\Delta F_{\text{stem}}(E/D)$, $\Delta s_\parallel$, and $\Delta s_\perp$ invariant to cell flattening.
   * By contrast, $P_{\text{ctrl}}(D, E)$ uses per-pixel ADU $(D, E)$, which scales with $h_{\text{cell}}$. To control for thinner/flatter cells at colony edges without new wet-lab imaging, we benchmark edge cells against an **edge-distance-stratified Control reference** $P_{\text{ctrl}}(D, E \mid d_{\text{edge}})$ (as well as day-matched and leave-one-well-out Control surfaces).
2. **The Colony-Edge vs. Interior Physiological Gradient (`24 / 24` Wells):**
   * Across all 24 wells on Days 2.0 and 3.0, colony-edge cells (`d_edge < 40 px`) are up-shifted in $\Delta F_{\text{stem}}$ by `+2.0` to `+2.9 mHz` in pluripotent `Control` (`+3.80 mHz` between `0–20 px` and `> 80 px`) and by `+0.5` to `+6.6 mHz` across all lineages at matched `E/D`.
   * Conversely, the deep frequency slowdowns in **Ectoderm (`-10.22 mHz`)** and **Endoderm (`-5.34 mHz`)** are concentrated in the **colony interior (`> 80 px`)**, while the outer `0–20 px` rim lags near `-1.91 mHz` and `-1.35 mHz`. This explains why the Ectoderm rosette-rim cells (`Point0045`) remain active and less shifted along the highway.
3. **Refined Framing for Mesoderm Axis 2 ("Broadened / Decoupled Competence" Rather Than Blanket "Hyper-Excitability"):**
   * Mesoderm oscillates at low per-pixel $(D, E)$ where Control is silent (`Pink` in `8C`, including `28–32%` in the deep interior `> 40 px` vs. `5–7%` in Control), while at high $(D, E)$ some real-sized Mesoderm masks are quieter than Control (`71–82%` vs. `97%`). We therefore replace promotional language ("hyper-excitable", "more alive") with precise biophysical terminology: **flattened / low-threshold $(D, E)$ competence envelope**.
4. **Well-Level Replication (`n = 3` Wells per Condition-Day):**
   * All headline summary claims ($\Delta F_{\text{stem}}$, $\Delta s_\parallel$, $\Delta s_\perp$, $\Delta P_{\text{osc}}$, `% Quenched`, `% Low-Threshold Active`, and neighbor correlation $r$) will be reported with **well-level means and CIs (`n = 3` biological replicate wells)**.

---

## 2. Step-by-Step Execution Plan

### Phase 1: Figure 49 — Colony-Edge Zonation & Edge-Conditioned Non-Clonal Neighbor Coherence
**Script to create:** `/home/smcoyle/h9_vocoder_work/workspace_archive/plot_fig49_colony_edge_zonation_and_conditioned_neighbor_coherence.py`
**Output artifact:** `Fig49_H9_Colony_Edge_Zonation_and_Edge_Conditioned_Neighbor_Coherence.png` (saved to artifact dir and synced to `Z:` drive `/mnt/researchdrive/General/Scott/stemcellcodes/vocoder_2026_10/figures/`).

* **Panel A (Radial $\Delta F_{\text{stem}}$ & $\Delta s_\parallel$ Profiles Across Time & Lineage):**
  * Continuous radial profiles of $\Delta F_{\text{stem}}(d_{\text{edge}})$ and $\Delta s_\parallel(d_{\text{edge}})$ from `0 px` (colony boundary) to `> 100 px` (deep interior) across `d0.0 → d3.0` for all 4 lineages, with **well-level ribbons (`n = 3` wells)**.
  * Highlights the steady pluripotent `Control` edge gradient (`+2.90 mHz` at `0–20 px` vs. `-0.90 mHz` at `> 80 px`) and the **inside-out differentiation downshift** in `Ectoderm` and `Endoderm`.
* **Panel B (Well-by-Well Reproducibility & Within-`E/D`-Bin Invariance):**
  * Shows all `24 / 24` wells (`< 40 px` edge vs. `≥ 40 px` interior on Days 2.0 and 3.0) and verifies that the radial gradient holds within narrow `0.5`-wide `E/D` stoichiometry bins and after accounting for per-pixel `log D`.
* **Panel C (Spatial Map of $d_{\text{edge}}$ vs. $\Delta F_{\text{stem}}$ on Representative Day 3.0 Colonies):**
  * Visual side-by-side of colony boundary distance contours (`20, 40, 80 px`) overlaid on $\Delta F_{\text{stem}}$ for `Control`, `Ectoderm`, `Endoderm`, and `Mesoderm`.
* **Panel D (Non-Clonal Neighbor Coherence vs. Joint `(E/D, d_edge)`-Matched Permutation Null):**
  * Tests touching cell pairs with discordant clonal barcodes ($|\Delta(E/D)| > 0.5$) against three hierarchy-of-rigor nulls:
    1. `E/D`-matched random pairs,
    2. **Joint `(E/D, d_edge)`-matched permutation null within the same FOV** (removing the shared colony-edge gradient!), and
    3. Spatial decay of residual $\Delta F_{\text{stem}}$ correlation as a function of centroid distance (`15 px → 250 px`) within the same radial zone.
* **Panel E (Cross-Boundary Phase Independence $\Delta\phi$ vs. Frequency Coherence):**
  * Computes the instantaneous oscillator phase difference $|\Delta\phi|$ across touching cell boundaries alongside $\Delta F_{\text{stem}}$ correlation, testing that neighboring cells share **physiological frequency state without phase-locking across membranes** (ruling out optical/FFT pixel bleed).

---

### Phase 2: Figure 50 — Reference Robustness Audit & Likelihood-Ratio Surprise Calibration
**Script to create:** `/home/smcoyle/h9_vocoder_work/workspace_archive/plot_fig50_reference_robustness_and_llr_surprise_audit.py`
**Output artifact:** `Fig50_H9_Reference_Robustness_EdgeMatched_DayMatched_LOMO_Audit.png` (saved to artifact dir and synced to `Z:` drive).

* **Panel A (Unifying `P_ctrl` + Leave-One-Well-Out & Day-Matched Stability):**
  * Compares the canonical pooled homeostatic `P_ctrl(D, E)` (`d1.5–d3.0`) against:
    1. **Leave-one-well-out (LOWO)** $P_{\text{ctrl}}^{(-w)}(D, E)$, and
    2. **Day-matched** $P_{\text{ctrl}}^{(t)}(D, E)$ (`Control` from the exact same imaging day $t$).
  * Reports well-level (`n = 3`) means and CIs for `% Grey`, `% Blue`, and `% Pink` across all 4 lineages.
* **Panel B (Edge-Distance-Stratified Reference $P_{\text{ctrl}}(D, E \mid d_{\text{edge}})$):**
  * Evaluates $\Delta P_{\text{osc}}$ and `8C` when edge pixels/cells (`0–20`, `20–40 px`) are benchmarked against **Control edge pixels/cells at the same radial distance**, quantifying exactly how much of the Mesoderm and Ectoderm-rim low-$(D, E)$ survival is intrinsic vs. optical path-length flattening.
* **Panel C (Principled Domain Log-Likelihood Ratio / Binomial-$Z$ vs. ROC Knots):**
  * Replaces the hard-coded piecewise interpolation knots in `8C` with a continuous **domain-integrated signed binomial log-likelihood ratio ($\text{LLR}$)** whose flat-Control deadband (`[q_05, q_95]`) is derived automatically from leave-one-well-out Control wells.
* **Panel D (Resolution of the Mesoderm Competent-Mask Question from `Consensus §4`):**
  * Evaluates real-sized masks (`area ≥ 100 px`) across **both Strict-Z and Relaxed-Z (`Class 1 Gold`)** oscillation tiers, stratified by $d_{\text{edge}}$ and $(D, E)$ zone, with `n = 3` well CIs.
* **Panel E (Variance Decomposition Summary Table):**
  * Fits a hierarchical binomial/linear mixed-effects decomposition across pixel tiles / cells (`Lineage`, `Day`, `d_edge`, `E/D` or `(D, E)`, plus `(1 | well)` and `(1 | FOV)`) to report the exact `% variance explained` by each biological factor.

---

### Phase 3: Codebase & Documentation Cleanup (Zero-Re-Render Updates + [`Fig 47`](file:///home/smcoyle/.gemini/antigravity-cli/brain/aba364ca-b51f-459e-82e0-cad8dfe277e5/Fig47_H9_Dual_ROC_Surprise_Flat_Control_Blue_Quench_Pink_Meso.png) Fix)
1. **Fix [`Fig 47`](file:///home/smcoyle/.gemini/antigravity-cli/brain/aba364ca-b51f-459e-82e0-cad8dfe277e5/Fig47_H9_Dual_ROC_Surprise_Flat_Control_Blue_Quench_Pink_Meso.png) Threshold Annotations Automatically in Code:**
   * Update [`plot_fig47_dual_roc_surprise_flat_control_blue_quench_pink_meso.py`](file:///home/smcoyle/h9_vocoder_work/workspace_archive/plot_fig47_dual_roc_surprise_flat_control_blue_quench_pink_meso.py) so $\tau_-^*$ and $\tau_+^*$ are derived programmatically from the ROC optimization and formatted directly into the figure labels without hard-coded string mismatches (`-22%` raw domain shift vs. `-28.5` rescaled score).
2. **Upgrade the Interactive HTML Explorer (`index.html` & `H9_FOV_Physiology_Atlas_Viewer.html`):**
   * **Fix the segmented-mask denominator:** Report real-sized segmented masks (`area ≥ 100 px`) alongside `n_cells` (`n_osc_pixels ≥ 10`) in the metric card so the denominator is not inflated by `6 px` fragments.
   * **Add the Interactive FOV Bookmark / Tagging Tool (`Consensus §2H`):** Add a lightweight `★ Star / Tag / Notes` bar in the viewer (persisted in `localStorage` with one-click `Export Tagged FOVs CSV`) so you can record notes while browsing the 420 FOVs.
   * **Neutralize promotional strings** in the viewer headers/captions ("Transcends Barcode!", "Vivid Pink Meso!") in favor of clean biophysical descriptions.
3. **Update [`README.md`](file:///home/smcoyle/h9_vocoder_work/README.md) (Local, GitHub, and `Z:` Drive):**
   * Update [`README.md`](file:///home/smcoyle/h9_vocoder_work/README.md) to document the **12-panel (`2 × 6`)** layout, the live **Overlay / Blend & Curtain-Wipe** controls, the exact `1.82M → 858k → 448k → 132k` data funnel, and the new **Figures 47–50**.

---

### Phase 4 (Optional Full-Atlas Re-Render When Ready)
Because re-rendering all 420 full-resolution (`4320 × 1575 px`) and web-optimized (`2880 × 1050 px`) 12-panel PNGs takes ~25 minutes on 12 workers, we will **hold off on re-rendering the 420 PNGs** until you review Figures 49 and 50. When you give the green light for a final atlas render pass, [`build_h9_fov_atlas_and_html_viewer.py`](file:///home/smcoyle/h9_vocoder_work/build_h9_fov_atlas_and_html_viewer.py) will also include:
* **Single unified `P_ctrl` reference** across Panels `7`, `8A`, `8B`, and `8C`, and
* **Panel 10 (`10. Cells on Sail`)** colored by $\Delta s_\perp^{\text{stem}}$ (or `8C` class) instead of noisy single-cell `osc_pct - P_ctrl`.
