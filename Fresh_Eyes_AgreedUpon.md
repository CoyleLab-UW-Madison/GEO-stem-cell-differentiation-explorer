# Fresh Eyes — Final Locked Consensus Plan (`Fresh_Eyes_AgreedUpon.md`)
**Authors:** Argon (Gemini) & Opus (Claude) — Unified Methodological & Biological Consensus  
**Dataset:** 420 FOVs (`4 lineages × 7 days × 3 biological replicate wells × 5 FOVs/well`), `448,054` audited oscillating cells (`n_osc_pixels ≥ 10`).

---

## 1. Locked Methodological Specification for $\Delta F_{\text{stem}}^{\text{2D}}(D, E)$

> [!IMPORTANT]
> **Unified 2D $(D, E)$ Architecture Across the Entire Manuscript:**
> Every decoded metric—$\Delta F_{\text{stem}}^{\text{2D}}(D, E)$ (`Panel 5`), $\Delta P_{\text{osc}}^{\text{stem}}(D, E)$ (`Panels 7, 8A–C`), and the 28-landscape Sail coordinates $(\Delta s_\parallel^{\text{stem}}, \Delta s_\perp^{\text{stem}})$ (`Panels 6 & 10`)—now shares the **exact same 2D $(D, E)$ reference architecture** on the `Control Days 1.5–3.0` homeostatic pool (`32 × 32` geomspace grid, $D \in [150, 950]\text{ ADU}$, $E \in [440, 3200]\text{ ADU}$), requiring **zero** ad-hoc $(B_D, B_E)$ background subtraction parameters in $\Delta F_{\text{stem}}$.

| Specification Item | Locked Implementation | Out-of-Sample Validation (`Control d1.5–3.0` LOWO) |
| :--- | :--- | :--- |
| **1. Primary $\Delta F_{\text{stem}}^{\text{2D}}$ (`Cell-2D`)** | Cell-level 2D lookup: each cell's mean $(D_{\text{cell}}, E_{\text{cell}})$ is matched to the median frequency $f_{\text{Ctrl, 2D}}(i_D, i_E)$ of `Control Days 1.5–3.0` cells in that exact $32 \times 32$ 2D bin. | **`r(D, resid) = -0.002`**, **`r(E, resid) = -0.053`**, **`RMSE = 8.27 mHz`** (`98.2%` cell support; outperforms `1D_BG` which had `r(D, resid) = +0.227`, `RMSE = 8.56 mHz`). |
| **2. Support Rule & Out-of-Support Flagging** | Each 2D $(i_D, i_E)$ bin requires **$\ge 20$ distinct `Control` cells** (no silent nearest-neighbor extrapolation). Cells in sparse bins (`1.9%` Control, `1.4%` Endo, `7.1%` Ecto, `13.4%` Meso on Day 3.0) are flagged (`#384459` cool slate in `Panel 5`), excluded from primary cell-level medians, and their % is explicitly reported. | Eliminates low-$D$ extrapolation bias where `Mesoderm` enters the `1-on-0` regime ($D < 185\text{ ADU}$). |
| **3. Sensitivity Bound (`Pixel-2D`)** | Pixel-level 2D lookup ($f_i - \mu_{\text{Ctrl, 2D}}(D_i, E_i)$ per strict oscillating pixel $i$, then averaged per cell, requiring $\ge 20$ distinct `Control` cells per bin; `99.9–100.0%` cell support). | **`r(D, resid) = -0.054`**, **`r(E, resid) = +0.154`**, **`RMSE = 8.39 mHz`**. Reported alongside `Cell-2D` so key effects are quoted with their full method range. |
| **4. Illumination-Dome Removal & `d_macro`** | Unmasked pixels (`mask == 0`) require `local_std_G >= 12.0 ADU` in addition to intensity above background, eliminating the circular empty-media illumination dome in 75 sparse FOVs (`Endoderm d2.0–3.0`). Macroscopic colony boundary `d_macro` uses edge-padded (`mode='edge'`) + pinhole-filled (`< 3,500` sub-px) distance transforms. | 0/420 FOVs retain the empty-media illumination dome; FOV frame edges and tiny inter-nuclear gaps are never misclassified as colony boundaries. |

---

## 2. Locked Quantitative Table (Day 3.0, $n = 3$ Wells/Condition)

| Lineage (Day 3.0) | **Primary `Cell-2D` $\Delta F_{\text{stem}}^{\text{2D}}$** (`Overall` \| `Edge <40px` \| `Int ≥40px` \| `Edge−Int`) | `Cell-2D` Support | **Sensitivity `Pixel-2D` $\Delta F_{\text{stem}}^{\text{2D}}$** (`Overall` \| `Edge <40px` \| `Int ≥40px` \| `Edge−Int`) | `Pixel-2D` Support | **Raw Physical $f_{\text{mHz}}$** (`Edge` \| `Int` \| `Edge−Int`) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Control** | **`+0.25 ± 1.01`** (`+1.36` \| `+0.24` \| **`+1.12 ± 0.53`**) | `98.1%` | `+0.26 ± 0.98` (`+0.53` \| `+0.26` \| `+0.27 ± 0.38`) | `100.0%` | `32.83` \| `35.14` \| `-2.30 ± 1.27` |
| **Ectoderm** | **`-6.27 ± 0.71`** (`-1.16` \| `-7.65` \| **`+6.49 ± 0.42`**) | `92.9%` | `-6.21 ± 0.78` (`-3.52` \| `-7.07` \| `+3.55 ± 0.38`) | `100.0%` | `29.65` \| `29.36` \| `+0.29 ± 1.61` |
| **Endoderm** | **`-4.68 ± 0.26`** (`-1.91` \| `-4.85` \| **`+2.94 ± 0.46`**) | `98.6%` | `-4.35 ± 0.25` (`-2.96` \| `-4.46` \| `+1.50 ± 0.55`) | `100.0%` | `32.07` \| `31.49` \| `+0.58 ± 1.94` |
| **Mesoderm** | **`-1.62 ± 0.49`** (`-0.85` \| `-1.78` \| **`+0.94 ± 0.37`**) | **`86.6%`** | `-3.45 ± 0.29` (`-3.87` \| `-3.27` \| `-0.60 ± 0.20`) | `99.9%` | `29.62` \| `30.34` \| `-0.72 ± 0.22` |

---

## 3. Three Framing Guardrails for Figures & Text
1. **Quote primary (`Cell-2D`) per-zone numbers first, with `Pixel-2D` as the sensitivity range:**
   - *Ectoderm (`d3.0`):* overall `-6.27 ± 0.71 mHz` (`-6.21 mHz` pixel-2D), with interior `-7.65 mHz` vs. rim `-1.16 mHz` (`Edge − Int = +6.49 ± 0.42 mHz` cell-2D; `-7.07` vs. `-3.52 mHz`, `+3.55 ± 0.38 mHz` pixel-2D; `+3.3` to `+5.4 mHz` net of Control's own edge difference).
   - *Endoderm (`d3.0`):* overall `-4.68 ± 0.26 mHz` (`-4.35 mHz` pixel-2D), with interior `-4.85 mHz` vs. rim `-1.91 mHz` (`Edge − Int = +2.94 ± 0.46 mHz` cell-2D; `+1.50 ± 0.55 mHz` pixel-2D; `+1.2` to `+1.8 mHz` net of Control).
   - *Mesoderm (`d3.0`):* modest clock slowdown of `-1.62 ± 0.49 mHz` (`Cell-2D`, `86.6%` support) to `-3.45 ± 0.29 mHz` (`Pixel-2D`, `99.9%` support), with its primary differentiation signature residing in **`Tilt` / $\Delta s_\perp^{\text{stem}}$** and **`1-on-0` Vivid Pink** persistence at low `MinD`.
2. **Flatness in `Control` and `Mesoderm` is consistent with 3D ER volume-wave robustness, not a standalone mechanistic proof:**
   - State clearly that `Control` (`Edge − Int = +0.27` to `+1.12 mHz`) and post-EMT `Mesoderm` (`-0.60` to `+0.94 mHz`, where cells form a dispersed mesenchymal lawn rather than cohesive epithelial colonies) show no systematic decoded-state edge zonation, consistent with mammalian 3D ER-supported volume waves being insensitive to boundary flattening per se.
3. **Present the biological interpretation of `Ectoderm` / `Endoderm` edge zonation as a testable spatial hypothesis:**
   - The decoded-state edge zonation in `Ectoderm` (`+3.55` to `+6.49 mHz`) and `Endoderm` (`+1.50` to `+2.94 mHz`) arises because dense interior rosette/sheet cores maintain higher $(D, E)$ expression (`Expected f_Ctrl ≈ 36.4–38.1 mHz`) while oscillating at `29.4–31.5 mHz`, whereas outer colony margins have lower $(D, E)$ (`Expected f_Ctrl ≈ 33.8–35.1 mHz`). Frame the reading that interior rosette/sheet cores are further along in differentiation while colony rims lag behind as a **spatial hypothesis testable by rim-vs-core fate-marker co-staining (`PAX6` / `SOX2` in Ectoderm, `SOX17` / `FOXA2` in Endoderm)**.
