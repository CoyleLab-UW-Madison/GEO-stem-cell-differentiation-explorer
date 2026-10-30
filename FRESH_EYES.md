# FRESH_EYES.md — Independent Review of the H9 MinDE Phase-Vocoder Project

*Written 2026-10-29 after re-reading the conversation, the build script, the `workspace_archive/` figure scripts, the README, and running three new sanity checks (Sections 3 and 3b; Section 3 revised after your ER-wave feedback).*

---

## 1. What we did (my own summary)

**Experiment.** H9 hESCs expressing a synthetic MinD/MinE oscillator were singularized and taken through trilineage differentiation (`control`, `ecto`, `endo`, `meso`). They were imaged at 7 timepoints (`d0.0 → d3.0`, every 12 h), with 3 wells × 5 FOVs per condition-day, giving **420 FOVs**.

**Signal extraction.** We ran a per-pixel FFT / phase-vocoder on each time series to get the expression levels (`D`, `E`), amplitude metrics (Z-scores), and frequency (peak bin + vocoder estimate). Amplitude-based Cellpose masks gave **1,817,010 segmented objects**. After a posterior confidence audit, **448,054** of these were kept as high-confidence oscillating cells.

**Axis 1: "how" cells oscillate (frequency encoding).**
- Raw frequency mostly reflects each clone's `E/D` stoichiometry (the "FM barcode").
- We subtracted the Control frequency surface at matched stoichiometry. That gives **ΔF_stem**.
- Phase-portrait fits split each cell's shift into **Δs∥** (along the main differentiation path) and **Δs⊥** (orthogonal to it). Δs⊥ is enriched in Mesoderm and in freshly singularized d0 cells. Its ROC threshold is τ* = +0.85 mHz.
- In the colony maps, ΔF_stem forms multi-cell patches that span neighbours with different `E/D`.

**Axis 2: "whether" cells oscillate (the dual).**
- **P_ctrl(D, E)** is the probability that a steady-state Control pixel at that expression level oscillates.
- **ΔP_osc** compares each condition against P_ctrl. Three spatial versions exist:
  - **8A**, the population surface shift
  - **8B**, a smoothed local observed − expected
  - **8C**, a dual-ROC "surprise" map
- Main picture: Ecto/Endo interiors are quenched (Blue). Mesoderm has many active cells sitting at stoichiometries where Control would usually be silent (Pink).

**Deliverables.**
- A 12-panel atlas for all 420 FOVs.
- A GitHub Pages viewer with blend/overlay.
- Figures 34 and 40–48.
- Full-resolution copies on the Z: drive.

Your own distillation from the conversation was accurate. The rest of this doc covers what I would challenge before any of it goes into a manuscript.

---

## 2. Corrections I owe you (things said earlier that don't hold up)

1. **"Every segmented mask is oscillating, so a mask can't be Blue in 8C." This is wrong.**
   - The masks are amplitude-segmented. Only **448k of 1.82M (25%)** made it into the audited cell set.
   - In 8C, Grey/Pink is assigned only to audited cells. Every other tissue pixel can turn Blue, including the ~1.37M segmented objects that failed the audit.
   - So Blue means "competent tissue that is not a confident oscillator". That is not the same as "tissue with no mask".

2. **Ectoderm "edge spindle" mechanism. My explanation was contradicted by my own table.**
   - I said the Pink rim cells survive by a *deeper* highway downshift.
   - The numbers I showed say the opposite. Pink cells are *less* shifted than Grey cells: Δs∥ = +0.76 vs +3.70, and ΔF = −6.26 vs −7.79 mHz.
   - The table also mixed two populations: n = 727 from the full parquet and n = 228 from the audited set.
   - The "1-D confinement / longer period" story was speculation, not something we found. See Section 3b for what the data actually suggests.

3. **The Fig 47 thresholds are reported inconsistently.**
   - τ₋* appears as −22% in the script's log and as −28% in the figure text and HTML.
   - τ₊* has appeared as P_ctrl ≤ 32%, ≤ 35.5% and ≤ 38.5% at different times.
   - The cause is hard-coded interpolation knots in `compute_unified_roc_surprise_map` that are not tied to the ROC output. Pick one set, derive it in code, and print it into the figure automatically.

4. **We have two different "Control references".**
   - Panels 7 / 8B / 8C use `P_ctrl_2d` from the atlas bundle: control d1.5–3.0, `max_fovs=4` per day, geomspace bins.
   - 8A uses `P_ref` from the Fig 46 cache, which has different edges and a different sample.
   - They should be one object.

5. **Panel 10 still colours cells by single-cell `osc_fraction − P_ctrl`.** We agreed that quantity is mostly mask-boundary noise (±26% smear). Colour by Δs⊥, or by the 8C class, instead.

---

## 3. Sanity check: what drives 8C "Pink"? (d3.0, 2 FOVs/well, audited cells)

**Question:** is 8C "Pink" just a colony-edge / low-intensity effect, or an artifact of the σ = 2 smoothing blurring cell intensity into the background?

**Definitions:**
- *Edge* = mean distance to the tissue boundary < 40 px.
- *Smoothed* = the current pipeline (σ = 2 Gaussian, full mask).
- *Raw core* = unsmoothed intensity inside a 2-px-eroded mask.

| Lineage | Position | n cells | Pink % (smoothed) | Pink % (raw core) | median D (smoothed) | median D (raw core) |
|---|---|---:|---:|---:|---:|---:|
| Control | interior | 15,829 | 5.2 | 6.8 | 333 | 341 |
| Control | **edge** | 5,733 | **15.8** | 16.3 | 269 | 276 |
| Ecto | interior | 1,564 | 18.4 | 17.1 | 280 | 290 |
| Ecto | **edge** | 1,010 | **61.8** | 57.5 | 203 | 212 |
| Endo | interior | 6,545 | 8.3 | 9.4 | 316 | 322 |
| Endo | **edge** | 817 | **43.6** | 41.4 | 225 | 232 |
| Meso | interior | 3,473 | 31.8 | 27.9 | 241 | 251 |
| Meso | **edge** | 4,352 | **60.3** | 53.9 | 199 | 207 |

**What this says:**
- **Smoothing / partial-volume is *not* the explanation.** Eroded raw-core intensities give nearly the same Pink rates, and D rises by only ~8 ADU.
- **Colony-edge position strongly shapes Pink in every lineage, including Control.** Control edge cells are Pink 3× as often as interior cells, and edge cells have genuinely lower per-pixel D. This is not a nuisance to remove (see 3b). It is a spatial variable to model explicitly.
- **Over half of Mesoderm's audited cells (4,352 of 7,825) sit within 40 px of a tissue boundary**, so lineage and position are entangled there.
- **A lineage effect survives in colony interiors.** Interior Meso is 28–32% Pink and interior Ecto 17–18%, versus 5–7% for interior Control.

**Caveat on what "low D" means (revised after your ER point).**
- I originally framed this as the oscillator responding to cell shape. You're right that this is unlikely here. In mammalian cells the waves run on the ER, as a 3D volume-filling process, and are robust to shape.
- The narrower, still-valid point concerns the **readout**, not the wave physics. P_ctrl(D, E) uses *intensity per pixel*, which is concentration × optical path length.
- A flatter cell with the *same* concentration reads as lower D and lower E (same E/D). It is therefore looked up in a lower-P_ctrl bin and can be called Pink without anything unusual about its oscillation.
- **ΔF_stem is immune to this** because it is referenced on E/D. Only the absolute-(D, E) metrics (panels 7, 8A–C, and the dead-zone classes) need a height or density covariate.

---

## 3b. Edge cells are frequency-shifted, and it's robust (follow-up check)

ΔF_stem for all audited cells on Days 2.0 and 3.0, binned by distance to the tissue boundary. Medians in mHz; d3.0 shown:

| Lineage | 0–20 px | 20–40 | 40–80 | >80 (interior) |
|---|---:|---:|---:|---:|
| Control | **+2.90** | +1.12 | +0.12 | −0.90 |
| Ecto | **−1.91** | −5.26 | −7.91 | **−10.22** |
| Endo | **−1.35** | −1.85 | −3.80 | **−5.34** |
| Meso | **−0.43** | +0.19 | −1.08 | −2.99 |

- **Every well shows it: 24 of 24** (4 lineages × 3 wells × 2 days). Edge (<40 px) minus interior ΔF_stem is positive in all of them, ranging from +0.5 to +6.6 mHz. Control is +2.0 to +2.9; Ecto at d3.0 is +5.7 to +6.6.
- **It is not E/D composition.** In Control, the edge up-shift persists inside every 0.5-wide E/D bin, from +1.3 to +3.6 mHz.
- **It is mostly not explained by absolute D.** Within E/D bins in Control, the slope is ΔF/∂log D = −3.0 mHz per e-fold. The edge-vs-interior difference in D (≈335 → 270) therefore predicts only **+0.65 mHz** of the observed +2–3 mHz. Cell area has no effect (slope −0.07).
- **Δs∥ follows the same gradient**: edge negative, interior positive, in every lineage. Interior cells are further along the highway.
- **The lineage downshifts are largely interior phenomena.** Ecto interior is −10.2 mHz, but Ecto edge cells sit near Control (−1.9). Endo shows the same pattern, more mildly.

**My read:** this is a genuine, very consistent spatial phenomenon, and an interesting one.
- hESC colony edges are a known distinct niche: different mechanics, YAP/TAZ state, BMP/TGF-β receptor accessibility, and edge-first fate in micropatterned colonies.
- Our oscillator reports that the edge niche has a distinct physiological state *in all four conditions, including pluripotent Control*.
- Differentiation-driven downshifts then propagate from the interior outward.

So yes, disentangling is the goal, but not to remove the edge. The aim is to estimate **lineage**, **position (edge ↔ interior)** and **their interaction** as separate effects, with the absolute-intensity readout caveat handled separately for the P_osc metrics.

---

## 4. Methodological issues to fix before the manuscript

**A. Circularity / in-sample evaluation**
- P_ctrl is built from Control FOVs, and "Control stays 91% grey" is then measured on Control FOVs from the same pool. Use **leave-one-well-out** surfaces: build P_ctrl without well *k*, then evaluate on well *k*.
- The ROC thresholds (Fig 41, Fig 47) were chosen and reported on the same d2–3 data. Report them with cross-validation, or calibrate on one replicate set and test on the others.
- The showcase FOVs were swapped during the session (for good reasons). For any claim drawn from them, also show the **distribution across all 15 FOVs**. The viewer's table already has the numbers.

**B. Unit of replication**
- Pixels and cells are not independent. The replicate unit is the **well** (n = 3 per condition-day).
- Headline effects (ΔF, ΔP, % Pink, % Blue) need well-level summaries with CIs, or mixed models with random effects for well and FOV.
- With N ≈ 448k cells, any cell-level p-value is meaningless on its own.

**C. Absolute-intensity coordinates across imaging days**
- P_ctrl, the tissue gates (`E ≥ 560`, `D ≥ 188` ADU) and the dead-zone classes are all in **absolute ADU**.
- If illumination, detector or medium drifted between imaging days, every ΔP shifts.
- Control is imaged on every day. Use **day-matched Control** (or per-day intensity normalization against Control) and show the result is stable.

**D. What 8C's "Pink" actually measures**
- For audited cells, Pink depends *only* on `P_ctrl(D_cell, E_cell)`. It is a stoichiometry/intensity classifier applied to cells already selected as oscillating.
- AUROC = 0.83 therefore partly means "Meso oscillators have different (D, E)", which mixes expression, colony position, and the intensity readout, not necessarily dynamics.
- This is fine as a visualization, but the paper's language ("hyper-excitable", "more alive than control") should not get ahead of it.

**E. Magic numbers**
- Z thresholds (5/4), tissue gates, smoothing σ (2, 6, 10, 12), competence gate (20 vs 22%), ±12/15% bands, and the interpolation knots in the surprise map all need either a principled derivation or a robustness sweep (a supplementary figure that shows the conclusions survive ±50% changes).

**F. Spatial-coherence claim (the most exciting one) needs a proper null**
- "Coherent ΔF_stem across non-clonal neighbours" is the most exciting result, so it deserves the hardest test.
- Sources of artificial neighbour correlation in our own pipeline:
  - σ = 2 smoothing of D/E
  - FFT/vocoder pixel bleed
  - segmentation boundary sharing
  - shared local environment (density, edge distance)
- The edge gradient in 3b is itself a large source of neighbour correlation: neighbours share distance-to-edge. The coherence claim should be tested **after conditioning on distance-to-edge**.
- Tests:
  - a permutation null that shuffles ΔF among cells **matched on E/D and edge distance**
  - correlation vs distance *beyond* the smoothing scale
- Because MinDE is intracellular, we expect **no phase-locking** between neighbours but **frequency correlation** if they share cell state. Showing exactly that pattern would make a strong figure.

---

## 5. New directions

1. **Replace binned ΔP and ROC fixes with a cell-level model.** This uses *all* 1.82M segmented objects, including the 0s:
   ```
   logit P(osc_i) = s(log D_i, log E_i) + s(edge_dist_i) + s(area_i) + s(local_density_i)
                    + lineage × day + (1 | well) + (1 | FOV)
   ```
   - This gives calibrated ΔP with CIs and **decomposes the effect** into expression, colony position and lineage.
   - It directly separates lineage effects from edge-niche effects and from the absolute-intensity readout.
   - It resolves the "discrete 0/1 state" worry you raised: the model is the per-cell probability, and each observed 0/1 is a draw from it.

2. **Information-theoretic surprise in place of the ROC-tuned map.**
   - Per-cell surprise = −log₂ P_ctrl(observed state | covariates). It treats 0s and 1s symmetrically.
   - Summed over a domain it becomes a log-likelihood ratio, so "this domain is quenched" gets a real test statistic.
   - It needs no hand-tuned thresholds and should naturally give a flat Control.

3. **Formalize "Δs⊥ ~ ΔP_osc" as a hurdle model.**
   - Model P(osc) and E[f | osc] jointly with shared covariates.
   - Then ask: does Δs⊥ predict oscillation competence *after* adjusting for (D, E, edge distance)?
   - That is the rigorous version of your intuition.

4. **The colony-edge niche as a first-class spatial axis** (see 3b).
   - Model ΔF_stem, Δs∥, Δs⊥ and P_osc as smooth functions of distance-to-edge × lineage × day.
   - Ask whether lineage effects are mostly interior effects (Ecto/Endo look that way), and whether the edge "front" moves inward over days.
   - Link the edge state to known niche biology: YAP/TAZ, BMP4 response, pSMAD1/5 staining at edges.
   - Separately, add a density or height covariate (nuclei per area, or a z-stack on a subset) **only** for the absolute-intensity P_osc metrics. This is a readout correction, not a claim about wave physics, since the ER-supported waves are shape-robust.

5. **Ground truth and perturbations.**
   - Post-hoc immunofluorescence on the same FOVs (OCT4, PAX6, SOX17, T/BRA) to validate Δs⊥ as an EMT readout cell by cell.
   - Targeted perturbations to map which cellular axis drives ΔF and which drives P_osc:
     - ATP: 2-DG / oligomycin
     - volume: osmotic shifts
     - spreading/contractility: Y-27632, blebbistatin

6. **Day 0 singularization as a built-in positive control.**
   - Singularization should produce the Δs⊥ / Pink signature, which then relaxes.
   - That time course is a clean, quantifiable internal benchmark for every metric above.

---

## 6. Engineering & reproducibility

| Issue | Suggestion |
|---|---|
| 34 figure-numbered scripts importing each other (`plot_fig47` → `fig43` → `fig42` …), absolute paths, `sys.path` hacks, symlinks | Small package (`h9vocoder/`: `io.py`, `reference.py`, `metrics.py`, `render.py`) plus one `config.yaml` holding every threshold |
| `.pkl` caches with undocumented provenance | Store a hash of inputs and parameters with each cache; rebuild automatically when they change |
| 1,650-line build script with the HTML in a Python f-string | Move the HTML/JS into a template file; generate only the JSON |
| README is stale (still describes the 10-panel / 2×5 layout and quotes N = 857,964 instead of the 448,054 audited cells) | Regenerate the README from the same metadata as the viewer |
| 200 MB of PNGs re-committed on every render, plus amend/force-push history | Deploy from an orphan `gh-pages` branch that is overwritten each time (or use Git LFS / release assets). Keep `main` for code only |
| Thresholds and colormaps baked into 420 PNGs (any change = full re-render) | Export per-FOV label maps + per-cell tables (npz/JSON) and colour client-side (canvas/WebGL), so thresholds can be adjusted live in the viewer |
| `credential.helper store` saves your token in plaintext (`~/.git-credentials`) | Switch to `gh auth login`, or `git-credential-libsecret` / `cache --timeout` |
| No tests | Add a few golden-value tests (P_ctrl surface checksum, cell counts per condition, one FOV's metrics) to catch silent drift |
| Figure captions/HTML use promotional language ("Transcends Barcode!", "Vivid Pink Meso!") | Neutral wording before sharing with collaborators or reviewers |

**For your manual browse of the 420 FOVs:** a small viewer feature would help. Let you tag or annotate FOVs (stored in `localStorage`, exportable as CSV), so your qualitative impressions become data you can compare against the metrics. A **distance-to-edge overlay** panel would also make the 3b gradient visible while you browse.

---

## 7. Suggested priority order

1. Leave-one-well-out P_ctrl + day-matched normalization. Cheap, and it protects every ΔP claim.
2. The cell-level GLMM (Section 5.1) with distance-to-edge × lineage. This separates lineage, niche and readout effects for both ΔF and P_osc. The edge gradient (3b) is a strong candidate headline result in its own right.
3. Spatial-coherence null (conditioned on edge distance) + phase-locking test. This protects the most exciting claim.
4. Well-level statistics for all headline numbers.
5. Refactor + `gh-pages` deployment.
6. Then the wet-lab extensions (IF ground truth, edge-niche markers, perturbations).
