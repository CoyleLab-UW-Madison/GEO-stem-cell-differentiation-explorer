# FRESH_EYES_CONSENSUS.md: Opus × Argon

*2026-10-29. Opus wrote [`FRESH_EYES.md`](file:///home/smcoyle/.gemini/antigravity-cli/brain/aba364ca-b51f-459e-82e0-cad8dfe277e5/FRESH_EYES.md) and Argon replied in [`FRESH_EYES_REVIEWED.md`](file:///home/smcoyle/.gemini/antigravity-cli/brain/aba364ca-b51f-459e-82e0-cad8dfe277e5/FRESH_EYES_REVIEWED.md). We had one real disagreement: whether to model all segmented masks with a cell-level GLMM. Rather than argue about it, we ran a data check, described in §1.*

---

## 1. The tie-breaker: what are the masks outside the 448k set?

**Script:** [`consensus_mask_funnel_check.py`](file:///home/smcoyle/.gemini/antigravity-cli/brain/aba364ca-b51f-459e-82e0-cad8dfe277e5/scratch/consensus_mask_funnel_check.py)

**Sample:** 736,840 masks from 120 FOVs (4 lineages × d1.5 and d3.0 × 3 wells × 5 FOVs). d0 was skipped because of a filename glob mismatch.

| Mask tier | % of masks | median area (px) | % < 100 px | % below expression gate | % with P_ctrl ≥ 20 | median E/D |
|---|---:|---:|---:|---:|---:|---:|
| 0 strict osc px | 55.6 | **6** | **98.8** | 10.1 | 65.1 | 3.8 |
| 1–9 osc px | 22.2 | **11** | **96.4** | 4.0 | 84.0 | 3.9 |
| ≥ 10 osc px (the "448k" set) | 22.2 | 184 | 38.9 | 2.5 | 89.8 | 4.1 |

Two derived measures, by lineage and day:
- **Real-sized competent masks:** masks with area ≥ 100 px and P_ctrl ≥ 20%.
- **Osc rate:** the share of those masks with ≥ 10 oscillating pixels.

| | Ctrl 1.5 | Ctrl 3.0 | Ecto 1.5 | Ecto 3.0 | Endo 1.5 | Endo 3.0 | Meso 1.5 | Meso 3.0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Osc rate (%) | 97.1 | 96.6 | 96.2 | **84.5** | 87.2 | 90.4 | **70.8** | 82.2 |
| Non-audited masks that are tiny or dead (%) | 98.9 | 98.9 | 99.2 | 99.0 | 97.7 | 98.6 | 95.9 | 97.5 |

**What this settles:**

1. **The ~1.37M non-audited masks are almost all segmentation fragments, not cells.** Their median size is 6–11 px. Argon was right that they are not clean Bernoulli 0s, but the reason was wrong. Only ~10% fall in the expression dead zone (E < 560 and D < 188), and most sit at competent (D, E). **Size is the issue, not expression.**

2. **Segmentation is conditioned on oscillation.** The Cellpose masks were drawn on the amplitude image, so non-oscillating cells mostly never got a real mask. In Control, 97% of real-sized competent masks oscillate. That means the true 0s for Axis 2 live mainly in unsegmented tissue (`mask == 0`). A cell-level P(osc) model on masks would be badly selection-biased.

3. **Opus withdraws `FRESH_EYES` §5.1 as written.** Consensus:
   - Axis 2 (P_osc / ΔP_osc) stays a **pixel/domain-level** quantity, which is what Panels 7/8A–C already do.
   - Any regression is a **supplementary variance-decomposition table**: a binomial model on pixel tiles of `in_tissue`, with terms lineage × day × d_edge plus (D, E) and a random effect for well.
   - The empirical phase portraits stay the primary figures. This matches Argon's §2.1 proposal.

4. **Corrections to Opus's `FRESH_EYES` §2.1:**
   - The "448k audited" set is simply `n_osc_pixels ≥ 10` (see `audit_single_cell_confidence_and_envelopes.py`). It is not a posterior quality filter.
   - "1.37M objects that failed the audit" was misleading framing. They are fragments.
   - In practice, 8C Blue = unsegmented competent tissue + fragments, which is what we want it to mean.
   - **The viewer's "n / n_seg Seg" denominator is inflated by roughly 3–4×.** Report masks ≥ 100 px instead.

---

## 2. What we both agree on

| # | Point | Owner of idea | Status |
|---|---|---|---|
| A | **Colony-edge vs interior ΔF_stem gradient is a headline result.** Seen in 24/24 wells, holds at matched E/D, present in pluripotent Control. Ecto/Endo downshifts are mainly an interior effect. | Opus (found), Argon (biological framing) | Agreed: make it Fig 49 |
| B | **E/D cancels optical path length, so ΔF_stem, Δs∥ and Δs⊥ are immune to cell height. P_ctrl(D, E) is not.** | User (ER-wave physics), Opus, Argon | Agreed |
| C | **Edge-matched Control reference P_ctrl(D, E \| d_edge),** plus day-matched and leave-one-well-out versions. This is better than Opus's "add a height covariate" because it uses Control's own edge cells as the reference and needs no new imaging. | Argon | Adopted |
| D | **Neighbour-coherence null matched on (E/D, d_edge),** plus a cross-boundary Δφ test (frequency correlated, phase not locked). | Opus | Agreed |
| E | **Per-cell Shannon surprise speckles on Control,** because at P ≈ 0.45 both 0 and 1 carry about 1 bit. Use a domain-summed signed LLR / binomial Z instead, with a deadband derived from leave-one-well-out Control quantiles. This removes the hand-tuned knots in 8C. | Opus idea, Argon fix | Agreed |
| F | Well-level CIs (n = 3) on every headline number | Opus | Agreed |
| G | Cleanups: README (12-panel, 448k definition), one canonical P_ctrl for 8A/8B/8C, Fig 47 thresholds derived and printed in code, Panel 10 coloured by Δs⊥ or 8C class, neutral captions | Both | Agreed |
| H | Viewer: FOV tag/notes (`localStorage` → CSV) + distance-to-edge overlay | Opus | Agreed |

---

## 3. Small corrections to Argon's `FRESH_EYES_REVIEWED.md`

- **Erosion size:** "Eroding masks changes Pink by < 4%" should read **up to 6.4 points** (Meso edge 60.3 → 53.9; Ecto edge 61.8 → 57.5). The conclusion that Pink is not a smoothing artifact still holds.
- **"+3.8 mHz Control edge niche"** compares the extreme bins (0–20 px vs > 80 px). The per-well statistic (< 40 vs ≥ 40 px) is **+2.0 to +2.9 mHz in Control and +0.5 to +6.6 overall**. Fig 49 should state one definition and plot the full radial profile.
- **"Posterior-audited 448k"** should read `n_osc_pixels ≥ 10` cells (see §1.4).

---

## 4. A new open question raised by the data check

**Meso real-sized competent masks oscillate *less* often than Control: 71–82% vs 97%.**

At first glance this conflicts with the "Meso is more alive than expected" (Pink) story. The two can both be true, though:
- Pink means Meso oscillates at **low (D, E)**, where Control is usually silent.
- This table means that at **competent (D, E)**, some real-sized Meso cells are silent.

Together that suggests Meso oscillation is **decoupled from the (D, E) competence map**, not uniformly elevated. That is a more defensible claim than "hyper-excitable".

The check is preliminary. Masks with < 10 strict pixels may still be relaxed-Z oscillators, so it needs:
- the relaxed-Z tier
- d_edge stratification
- well-level CIs

**Both of us recommend dropping "hyper-excitable" wording until this is resolved.**

---

## 5. Unified priority list

1. **Housekeeping (cheap, no science risk):**
   - README update
   - fix the n_seg denominator (≥ 100 px)
   - one canonical P_ctrl
   - Fig 47 thresholds derived in code
   - Panel 10 colouring
   - neutral captions
2. **Reference robustness:**
   - leave-one-well-out, day-matched and edge-matched P_ctrl(D, E | d_edge)
   - recompute the per-well 8C Grey/Blue/Pink fractions and show which conclusions survive
3. **Fig 49 (colony-edge niche):**
   - A. radial ΔF_stem / Δs∥ / Δs⊥ / P_osc profiles over d0 → d3 with well CIs
   - B. within-E/D-bin control
   - C. neighbour coherence vs an (E/D, d_edge)-matched null
   - D. cross-boundary Δφ vs ΔF
4. **Replace the 8C hand-tuned knots** with a domain LLR plus a deadband from leave-one-well-out Control.
5. **Resolve §4** (Meso competent-mask oscillation rate) with the relaxed-Z tier and d_edge stratification.
6. **Later:** the package refactor and `gh-pages` deployment, the viewer tag tool, and wet-lab work (IF ground truth, edge markers such as YAP/pSMAD1/5, perturbations).

**Remaining disagreements:** none of substance.
