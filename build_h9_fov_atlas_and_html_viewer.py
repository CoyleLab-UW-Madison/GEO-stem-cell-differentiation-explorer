#!/usr/bin/env python3
"""
Builds the 5-Panel (1x5 Horizontal Strip) 420-FOV GEO Stem Cell Differentiation Atlas
and Interactive HTML Explorer (`index.html` & `H9_FOV_Physiology_Atlas_Viewer.html`).

Each of the 420 FOVs is rendered on a strict, pixel-exact 1-Row x 5-Column grid
(`22.5 x 4.70 in` @ `100 DPI` = `2250 x 470 px`) where all 5 spatial panels place their
`512 x 512` colony image at the exact same relative sub-box (`[x0 + 0.007, 0.115, 0.162, 0.7755]`)
inside their `1/5`-width column (guaranteeing 0.00-pixel registration offset in overlay mode):

  1. Fluorescence Still (`MinE-GFP`, grayscale + 50 µm scale bar)
  2. FM Barcode `f` (`[18, 52] mHz`, `Turbo`)
  3. Frequency Shift `ΔF_stem` (`[-12, +12] mHz`, `Blue <-> Crisp White (0) <-> Red`)
  4. Competence Shift `ΔP_osc` (`[-50%, +50%]`, `Teal <-> Silver-Grey (±13% IQR) <-> Gold`)
  5. Stoichiometric Tilt `Δs_⊥` (`[-1.0, +2.8] mHz`, `Silver-Grey (<= +0.85 ROC) -> Vivid Purple`)
"""

import io
import json
import pickle
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter, binary_dilation
import tifffile

LOCAL_WORK_DIR = Path("/home/smcoyle/h9_vocoder_work")
ARCHIVE_DIR = LOCAL_WORK_DIR / "workspace_archive"
sys.path.insert(0, str(LOCAL_WORK_DIR))
sys.path.insert(0, str(ARCHIVE_DIR))

from plot_mockup_fig5_and_fig6 import (
    make_df_blue_red_cmap,
    make_dposc_blue_grey_gold_cmap,
    make_sperp_slate_magenta_cmap,
)

WELL_CKPT_DIR = LOCAL_WORK_DIR / "well_checkpoints"
LOCAL_MASK_DIR = LOCAL_WORK_DIR / "consolidated_masks"
STILLS_DIR = LOCAL_WORK_DIR / "consolidated_stills"
ATLAS_DIR = LOCAL_WORK_DIR / "fov_atlas"

REMOTE_OUT_DIR = Path("/mnt/researchdrive/General/Scott/stemcellcodes/vocoder_2026_10")
REMOTE_ATLAS_DIR = REMOTE_OUT_DIR / "fov_atlas"
GEO_REPO_DIR = Path("/home/smcoyle/GEO-stem-cell-differentiation-explorer")

TREATMENTS = ["control", "ecto", "endo", "meso"]
TREAT_LABELS = {
    "control": "hESC Control (mTeSR)",
    "ecto": "Ectoderm",
    "endo": "Endoderm",
    "meso": "Mesoderm",
}

_G_CELL_DICT = None
_G_P_CTRL_2D = None
_G_D_EDGES = None
_G_E_EDGES = None


def make_grayscale_fluorescence_cmap():
    return LinearSegmentedColormap.from_list(
        "grayscale_fluo",
        [
            (0.00, "#070a0f"),
            (0.15, "#1e2430"),
            (0.45, "#646e7d"),
            (0.75, "#b3bcc8"),
            (1.00, "#f4f6f9"),
        ],
        N=256,
    )


def _init_worker(bundle_path_str: str):
    global _G_CELL_DICT, _G_P_CTRL_2D, _G_D_EDGES, _G_E_EDGES
    with open(bundle_path_str, "rb") as f:
        b = pickle.load(f)
    _G_CELL_DICT = b["cell_dict"]
    _G_P_CTRL_2D = b["P_ctrl_2d"]
    _G_D_EDGES = b["d_edges"]
    _G_E_EDGES = b["e_edges"]


def process_single_well(ckpt_path_str: str) -> list:
    ckpt_path = Path(ckpt_path_str)
    stem = ckpt_path.stem
    day_folder, cond_folder = stem.split("__", 1)
    parts = cond_folder.split("_")
    treat = parts[0]
    day_val = float(parts[2])
    well_val = int(parts[4])

    df_w = pd.read_parquet(
        ckpt_path,
        columns=["y", "x", "filename", "D_z", "E_z", "bin_diff", "F_bin", "freq_mHz"]
    )

    cmap_gray = make_grayscale_fluorescence_cmap()
    cmap_barcode = plt.get_cmap("turbo")
    cmap_df = make_df_blue_red_cmap()
    cmap_dposc = make_dposc_blue_grey_gold_cmap()
    cmap_sperp = make_sperp_slate_magenta_cmap()

    norm_barcode = Normalize(vmin=18.0, vmax=52.0)
    norm_df = Normalize(vmin=-12.0, vmax=12.0)
    norm_dposc = Normalize(vmin=-50.0, vmax=50.0)
    norm_sperp = Normalize(vmin=-1.0, vmax=2.8)

    filenames = sorted(df_w["filename"].unique())
    fov_records = []

    for fov_idx, fn in enumerate(filenames):
        fov_short = fn.replace(".nd2", "")
        mp = LOCAL_MASK_DIR / day_folder / cond_folder / f"{fov_short}_masks.tif"
        gp = STILLS_DIR / day_folder / cond_folder / f"{fov_short}_g.tif"
        rp = STILLS_DIR / day_folder / cond_folder / f"{fov_short}_r.tif"
        if not (mp.exists() and gp.exists() and rp.exists()):
            continue

        mask = tifffile.imread(mp)
        g_im = tifffile.imread(gp).astype(np.float32)
        r_im = tifffile.imread(rp).astype(np.float32)

        df_w_fov = df_w[df_w["filename"] == fn]
        m_st = (
            (df_w_fov["D_z"] >= 5.0)
            & (df_w_fov["E_z"] >= 4.0)
            & (df_w_fov["bin_diff"] <= 1)
            & (df_w_fov["F_bin"] >= 3)
            & (df_w_fov["freq_mHz"] >= 14.0)
            & (df_w_fov["freq_mHz"] <= 80.0)
        )
        strict_yx = df_w_fov.loc[m_st, ["y", "x"]].to_numpy()

        step = 2
        g_raw_sub = g_im[::step, ::step]
        g_sm = gaussian_filter(g_im, sigma=2.0)[::step, ::step]
        r_sm = gaussian_filter(r_im, sigma=2.0)[::step, ::step]
        m_sub = mask[::step, ::step]

        local_std_g = np.sqrt(
            np.maximum(0.0, gaussian_filter((g_raw_sub - gaussian_filter(g_raw_sub, sigma=5.0)) ** 2, sigma=5.0))
        )
        in_tissue = (m_sub > 0) | (((g_sm >= 560.0) | (r_sm >= 188.0)) & (local_std_g >= 12.0))

        # 1. Grayscale Still (*_g.tif : MinE-GFP, pure unannotated microscopy still)
        p_lo, p_hi = np.percentile(g_raw_sub, [2.0, 99.2])
        g_norm_still = np.clip((g_raw_sub - p_lo) / max(p_hi - p_lo, 1.0), 0.0, 1.0) ** 0.80
        lum = 0.04 + 0.92 * g_norm_still
        rgb_gray = np.stack([lum, lum, lum], axis=-1).astype(np.float32)

        # Shared dark microscopy background for panels 2..5 (exact match to Figures 5 and 6)
        p02, p995 = np.percentile(g_raw_sub, [2.0, 99.5])
        g_norm = np.clip((g_raw_sub - p02) / max(p995 - p02, 1.0), 0.0, 1.0) ** 0.65
        bg_fluo = np.stack([0.50 * g_norm, 0.51 * g_norm, 0.54 * g_norm], axis=-1).astype(np.float32)
        bg_fluo[~in_tissue] *= 0.14

        m_r = np.roll(m_sub, -1, axis=1)
        m_d = np.roll(m_sub, -1, axis=0)
        cell_border = ((m_sub > 0) & (m_sub != m_r)) | ((m_sub > 0) & (m_sub != m_d))
        border_mult = 0.52

        rgb_bar = bg_fluo.copy()
        rgb_df = bg_fluo.copy()
        rgb_dposc = bg_fluo.copy()
        rgb_sperp = bg_fluo.copy()

        di = np.clip(np.searchsorted(_G_D_EDGES, r_sm) - 1, 0, len(_G_D_EDGES) - 2)
        ei = np.clip(np.searchsorted(_G_E_EDGES, g_sm) - 1, 0, len(_G_E_EDGES) - 2)
        p_exp_map = 100.0 * _G_P_CTRL_2D[di, ei]
        sigma_nb = 32.0 / step
        tis_nb = np.maximum(gaussian_filter(in_tissue.astype(np.float32), sigma=sigma_nb), 0.15)

        c_info = _G_CELL_DICT.get((treat, day_val, well_val, fn), None)
        max_lbl = int(mask.max())
        in_aud_cell = np.zeros(m_sub.shape, dtype=bool)

        if c_info is not None and len(c_info["label_id"]) > 0:
            c_ids = c_info["label_id"]
            valid_c = (c_ids > 0) & (c_ids <= max_lbl)
            c_v = c_ids[valid_c]
            f_v = c_info["freq_mHz"][valid_c]
            df_v = c_info["df_stem"][valid_c]
            sp_v = c_info["ds_perp_roc"][valid_c]
            dp_v = c_info["dP_osc_cell"][valid_c]
            ed_v = c_info["ED_ratio"][valid_c]

            flat_m = m_sub.ravel()
            flat_p = p_exp_map.ravel()
            cnts = np.bincount(flat_m, minlength=max_lbl + 1)
            sums = np.bincount(flat_m, weights=flat_p, minlength=max_lbl + 1)
            mean_p_per_lbl = np.where(cnts > 0, sums / np.maximum(cnts, 1), 55.0)
            p_cell_v = mean_p_per_lbl[c_v]
            dp_eff_v = np.where(
                np.isfinite(dp_v),
                np.where((dp_v >= 0.0) & (p_cell_v < 24.0), np.maximum(dp_v, 36.0 - p_cell_v), dp_v),
                0.0,
            )

            has_cell_lut = np.zeros(max_lbl + 1, dtype=bool)
            has_cell_lut[c_v] = True
            in_aud_cell = has_cell_lut[m_sub]

            lut_bar_c = np.zeros((max_lbl + 1, 3), dtype=np.float32)
            lut_df_c = np.zeros((max_lbl + 1, 3), dtype=np.float32)
            lut_dp_c = np.zeros((max_lbl + 1, 3), dtype=np.float32)
            lut_sp_c = np.zeros((max_lbl + 1, 3), dtype=np.float32)

            lut_bar_c[c_v] = cmap_barcode(norm_barcode(f_v))[:, :3].astype(np.float32)
            oos_m = ~np.isfinite(df_v)
            sup_m = np.isfinite(df_v)
            lut_df_c[c_v[oos_m]] = np.array([0.90, 0.92, 0.95], dtype=np.float32)
            if np.any(sup_m):
                lut_df_c[c_v[sup_m]] = cmap_df(norm_df(df_v[sup_m]))[:, :3].astype(np.float32)
            lut_dp_c[c_v] = cmap_dposc(norm_dposc(dp_eff_v))[:, :3].astype(np.float32)
            lut_sp_c[c_v] = cmap_sperp(norm_sperp(sp_v))[:, :3].astype(np.float32)

            tex_w = (0.76 + 0.24 * g_norm)[..., None]
            rgb_bar[in_aud_cell] = np.clip(lut_bar_c[m_sub[in_aud_cell]] * tex_w[in_aud_cell] * 1.18, 0.0, 1.0)
            rgb_df[in_aud_cell] = np.clip(lut_df_c[m_sub[in_aud_cell]] * tex_w[in_aud_cell] * 1.18, 0.0, 1.0)
            rgb_dposc[in_aud_cell] = np.clip(lut_dp_c[m_sub[in_aud_cell]] * tex_w[in_aud_cell] * 1.18, 0.0, 1.0)
            rgb_sperp[in_aud_cell] = np.clip(lut_sp_c[m_sub[in_aud_cell]] * tex_w[in_aud_cell] * 1.18, 0.0, 1.0)

            rgb_bar[cell_border & in_aud_cell] *= border_mult
            rgb_df[cell_border & in_aud_cell] *= border_mult
            rgb_dposc[cell_border & in_aud_cell] *= border_mult
            rgb_sperp[cell_border & in_aud_cell] *= border_mult
        else:
            c_v = np.array([], dtype=np.int32)
            f_v = np.array([0.0], dtype=np.float32)
            df_v = np.array([0.0], dtype=np.float32)
            sp_v = np.array([0.0], dtype=np.float32)
            dp_v = np.array([0.0], dtype=np.float32)
            ed_v = np.array([0.0], dtype=np.float32)

        # Multi-cell neighborhood competence on non-oscillating colony tissue in control-competent (D, E) space
        p_dom = 100.0 * 1.45 * (gaussian_filter(in_aud_cell.astype(np.float32), sigma=sigma_nb) / tis_nb)
        quenched_dom = in_tissue & (~in_aud_cell) & (p_exp_map >= 24.0)
        if np.any(quenched_dom):
            dp_dom = np.minimum(0.0, p_dom[quenched_dom] - p_exp_map[quenched_dom])
            dom_cols = cmap_dposc(norm_dposc(dp_dom))[:, :3].astype(np.float32)
            tex_q = (0.76 + 0.24 * g_norm[quenched_dom])[..., None]
            rgb_dposc[quenched_dom] = np.clip(dom_cols * tex_q * 1.16, 0.0, 1.0)

        # Summary statistics for this FOV
        n_cells = int(len(c_v))
        med_f = float(np.nanmedian(f_v)) if n_cells > 0 else 0.0
        iqr_f = float(np.nanpercentile(f_v, 75) - np.nanpercentile(f_v, 25)) if n_cells > 1 else 0.0
        med_df_stem = float(np.nanmedian(df_v)) if np.any(np.isfinite(df_v)) else 0.0
        med_dposc = float(np.nanmedian(dp_v)) if np.any(np.isfinite(dp_v)) else 0.0
        med_s_perp = float(np.nanmedian(sp_v)) if np.any(np.isfinite(sp_v)) else 0.0
        p90_s_perp = float(np.nanpercentile(sp_v, 90)) if np.any(np.isfinite(sp_v)) else 0.0
        n_meso_roc = int(np.sum(sp_v > 0.85)) if n_cells > 0 else 0
        pct_meso_roc = 100.0 * n_meso_roc / max(n_cells, 1)
        med_ed = float(np.nanmedian(ed_v)) if n_cells > 0 else 0.0

        rel_png = f"fov_atlas/{treat}_d{day_val:.1f}_w{well_val}_fov{fov_idx}_{fov_short}.png"
        out_png = LOCAL_WORK_DIR / rel_png

        # =====================================================================
        # RENDER CLEAN 1x5 HORIZONTAL STRIP (22.5 x 4.70 in @ 100 DPI = 2250 x 470 px)
        # Every Tile (c in {0..4}) occupies x in [c*0.20, (c+1)*0.20].
        # Every 512x512 colony image occupies [x0 + 0.007, 0.115, 0.162, 0.7755]
        # (exact 1:1 square = 364.5 x 364.5 px, 0.00 px overlay offset!).
        # =====================================================================
        fig = plt.figure(figsize=(22.5, 4.70), dpi=100)
        fig.patch.set_facecolor("#0f172a")

        fov_pt = fov_short.split("_")[0]
        tiles = [
            (
                0, rgb_gray, cmap_gray, Normalize(0, 100),
                "1. Fluorescence Still (MinE-GFP)",
                f"{treat.upper()} d{day_val:.1f} · well {well_val} ({fov_pt}) | {n_cells} cells",
                "#f1f5f9", None
            ),
            (
                1, rgb_bar, cmap_barcode, norm_barcode,
                "2. FM Barcode f (18–52 mHz)",
                f"median f = {med_f:.1f} mHz (IQR {iqr_f:.1f} mHz)",
                "#38bdf8", [20, 30, 40, 50]
            ),
            (
                2, rgb_df, cmap_df, norm_df,
                "3. Frequency Shift ΔF_stem (mHz)",
                f"median ΔF_stem = {med_df_stem:+.2f} mHz",
                "#f87171", [-10, -5, 0, 5, 10]
            ),
            (
                3, rgb_dposc, cmap_dposc, norm_dposc,
                "4. Competence Shift ΔP_osc (%)",
                f"median ΔP_osc = {med_dposc:+.1f}% (±13% IQR grey)",
                "#fbbf24", [-40, -20, 0, 20, 40]
            ),
            (
                4, rgb_sperp, cmap_sperp, norm_sperp,
                "5. Stoichiometric Tilt Δs_⊥ (mHz)",
                f"tail > +0.85 mHz: {pct_meso_roc:.1f}% (n={n_meso_roc}) | p90={p90_s_perp:+.2f}",
                "#c084fc", [-0.5, 0.5, 1.5, 2.5]
            ),
        ]

        w_col = 0.20
        for c_t, rgb_arr, cm, norm_obj, title_s, sub_s, t_col, cb_ticks in tiles:
            x0 = c_t * w_col
            ax = fig.add_axes([x0 + 0.007, 0.115, 0.162, 0.7755])
            ax.imshow(rgb_arr, interpolation="nearest", aspect="auto")
            # 50 µm scale bar (100 px at step=2) in bottom-right corner
            h_im, w_im = rgb_arr.shape[:2]
            ax.plot([w_im - 118, w_im - 18], [h_im - 18, h_im - 18], color="#ffffff", lw=2.2, solid_capstyle="butt")
            ax.set_xlim(0, w_im)
            ax.set_ylim(h_im, 0)

            ax.set_title(title_s, fontsize=10.0, fontweight="bold", color=t_col, pad=6)
            ax.set_xlabel(sub_s, fontsize=8.2, color="#cbd5e1", labelpad=5)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_edgecolor("#334155")
                sp.set_linewidth(1.0)

            if cm is not None and norm_obj is not None:
                cax = fig.add_axes([x0 + 0.1715, 0.115, 0.0055, 0.7755])
                cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm_obj, cmap=cm), cax=cax)
                if cb_ticks is not None:
                    cb.set_ticks(cb_ticks)
                cb.ax.tick_params(colors="#cbd5e1", labelsize=7.2, length=2.5)
                for sp in cax.spines.values():
                    sp.set_edgecolor("#475569")
                    sp.set_linewidth(0.6)

        # Save full-res to research drive if available, and web-optimized PNG to fov_atlas/
        remote_png = REMOTE_ATLAS_DIR / out_png.name
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=140, facecolor=fig.get_facecolor())
        plt.close(fig)
        if REMOTE_ATLAS_DIR.exists():
            remote_png.write_bytes(buf.getvalue())

        buf.seek(0)
        im_pil = Image.open(buf).convert("RGB")
        im_web = im_pil.resize((2250, 470), Image.Resampling.LANCZOS)
        im_q = im_web.quantize(colors=256, method=Image.Quantize.MAXCOVERAGE, dither=Image.Dither.NONE)
        im_q.save(out_png, optimize=True)

        fov_records.append(
            {
                "treatment": treat,
                "treatment_label": TREAT_LABELS[treat],
                "day": day_val,
                "well": well_val,
                "fov_index": fov_idx,
                "filename": fn,
                "fov_short": fov_short,
                "img_rel_path": rel_png,
                "n_pixels": int(len(strict_yx)),
                "n_cells": n_cells,
                "med_ED_ratio": round(med_ed, 2),
                "med_freq_mHz": round(med_f, 1),
                "iqr_freq_mHz": round(iqr_f, 1),
                "med_df_stem_mHz": round(med_df_stem, 2),
                "med_dposc_pct": round(med_dposc, 1),
                "med_s_perp_roc_mHz": round(med_s_perp, 2),
                "p90_s_perp_roc_mHz": round(p90_s_perp, 2),
                "pct_meso_roc": round(pct_meso_roc, 1),
                "n_meso_roc": n_meso_roc,
            }
        )

    return fov_records


def build_interactive_html_viewer(fov_records: list):
    for r in fov_records:
        for k, v in list(r.items()):
            if isinstance(v, float) and not np.isfinite(v):
                r[k] = 0.0
    records_sorted = sorted(
        fov_records,
        key=lambda r: (TREATMENTS.index(r["treatment"]), r["day"], r["well"], r["fov_index"]),
    )
    json_payload = json.dumps(records_sorted, allow_nan=False)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>GEO Stem Cell Differentiation Explorer — 420-FOV Atlas</title>
<style>
  :root {{
    --bg: #0b0f19;
    --panel: #111827;
    --panel-alt: #1e293b;
    --border: #2e3a4e;
    --text: #f1f5f9;
    --muted: #94a3b8;
    --accent: #38bdf8;
    --control: #475569;
    --ecto: #dc2626;
    --endo: #0d9488;
    --meso: #7e22ce;
    --blend-opacity: 0.65;
    --blend-mode: screen;
    --wipe-right: 15.5%;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    padding: 0;
    font-family: -apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background: var(--bg);
    color: var(--text);
    -webkit-font-smoothing: antialiased;
  }}
  header {{
    background: #0f172a;
    border-bottom: 1px solid var(--border);
    padding: 12px 24px;
    display: flex;
    flex-wrap: wrap;
    justify-content: space-between;
    align-items: center;
    gap: 12px;
    position: sticky;
    top: 0;
    z-index: 100;
  }}
  .title-block h1 {{ margin: 0; font-size: 16px; font-weight: 600; letter-spacing: -0.01em; color: #f8fafc; }}
  .title-block p {{ margin: 3px 0 0 0; font-size: 12px; color: var(--muted); }}
  .mode-tabs {{ display: flex; gap: 6px; flex-wrap: wrap; }}
  .mode-btn {{
    background: var(--panel-alt);
    color: #cbd5e1;
    border: 1px solid var(--border);
    padding: 6px 12px;
    border-radius: 5px;
    cursor: pointer;
    font-size: 12px;
    font-weight: 500;
  }}
  .mode-btn:hover {{ border-color: #64748b; color: #fff; }}
  .mode-btn.active {{ background: #e2e8f0; border-color: #e2e8f0; color: #0f172a; font-weight: 600; }}
  .controls-bar {{
    background: var(--panel);
    border-bottom: 1px solid var(--border);
    padding: 9px 24px;
    display: flex;
    flex-wrap: wrap;
    gap: 16px;
    align-items: center;
  }}
  .focus-bar {{
    background: #0d1322;
    border-bottom: 1px solid var(--border);
    padding: 8px 24px;
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
    align-items: center;
  }}
  .ctrl-group {{ display: flex; align-items: center; gap: 5px; flex-wrap: wrap; }}
  .ctrl-label {{
    font-size: 11px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.03em;
    color: var(--muted);
    margin-right: 4px;
  }}
  .pill-btn {{
    background: var(--panel-alt);
    color: #cbd5e1;
    border: 1px solid var(--border);
    padding: 4px 10px;
    border-radius: 4px;
    cursor: pointer;
    font-size: 11.5px;
    font-weight: 500;
  }}
  .pill-btn:hover {{ border-color: #64748b; color: #fff; }}
  .pill-btn.active {{ background: #38bdf8; color: #0f172a; border-color: #38bdf8; font-weight: 600; }}
  .pill-btn.treat-control.active {{ background: #e2e8f0; color: #0f172a; border-color: #e2e8f0; }}
  .pill-btn.treat-ecto.active {{ background: var(--ecto); color: #fff; border-color: var(--ecto); }}
  .pill-btn.treat-endo.active {{ background: var(--endo); color: #fff; border-color: var(--endo); }}
  .pill-btn.treat-meso.active {{ background: var(--meso); color: #fff; border-color: var(--meso); }}
  .kb-hint {{
    margin-left: auto;
    font-size: 11px;
    color: var(--muted);
  }}
  main {{ padding: 16px 24px 48px 24px; max-width: 1920px; margin: 0 auto; }}
  .metrics-strip {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(175px, 1fr));
    gap: 10px;
    margin-bottom: 14px;
  }}
  .metric-card {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 9px 13px;
  }}
  .metric-card .m-label {{ font-size: 10.5px; color: var(--muted); text-transform: uppercase; letter-spacing: 0.02em; font-weight: 600; }}
  .metric-card .m-val {{ font-size: 16px; font-weight: 600; margin-top: 3px; color: #f8fafc; }}
  .metric-card .m-sub {{ font-size: 11px; color: var(--muted); margin-top: 2px; }}
  .fov-card {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 12px 14px;
    margin-bottom: 14px;
  }}
  .fov-card-header {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 8px;
    margin-bottom: 10px;
    padding-bottom: 7px;
    border-bottom: 1px solid var(--border);
    font-size: 12.5px;
  }}
  .badge {{ padding: 2px 8px; border-radius: 4px; font-weight: 600; font-size: 11px; }}
  .badge-control {{ background: rgba(148,163,184,0.18); color: #e2e8f0; border: 1px solid #64748b; }}
  .badge-ecto {{ background: rgba(220,38,38,0.18); color: #fca5a5; border: 1px solid #dc2626; }}
  .badge-endo {{ background: rgba(13,148,136,0.18); color: #5eead4; border: 1px solid #0d9488; }}
  .badge-meso {{ background: rgba(126,34,206,0.22); color: #d8b4fe; border: 1px solid #9333ea; }}
  .fov-img-wrap img {{ width: 100%; height: auto; border-radius: 5px; display: block; background: #0f172a; }}

  .single-tile-grid {{
    display: flex;
    flex-wrap: wrap;
    gap: 16px;
    justify-content: center;
    align-items: flex-start;
  }}
  .tile-col-wrap {{
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 5px;
  }}
  .tile-col-tag {{
    font-size: 11.5px;
    font-weight: 600;
    color: #cbd5e1;
    background: #1e293b;
    border: 1px solid var(--border);
    padding: 3px 10px;
    border-radius: 4px;
  }}
  /* Each tile is 1/5 of 2250x470 = 450x470 px */
  .tile-zoom-box {{
    overflow: hidden;
    border-radius: 6px;
    border: 1px solid var(--border);
    background: #0f172a;
    position: relative;
    aspect-ratio: 450 / 470;
  }}
  .tile-zoom-box img.base-layer {{
    width: 500%;
    max-width: none;
    height: auto;
    display: block;
  }}
  /* Overlay layer clipped strictly to the 512x512 colony image box:
     left = 0.007/0.20 = 3.5%, right = 100% - (0.007+0.162)/0.20 = 15.5%,
     top = 10.95%, bottom = 11.5% */
  .tile-zoom-box .overlay-clip {{
    position: absolute;
    inset: 0;
    pointer-events: none;
    clip-path: inset(10.95% var(--wipe-right, 15.5%) 11.5% 3.5%);
    opacity: var(--blend-opacity, 0.65);
    mix-blend-mode: var(--blend-mode, screen);
  }}
  .tile-zoom-box .overlay-clip img {{
    width: 500%;
    max-width: none;
    height: auto;
    display: block;
  }}
  select.pair-select {{
    background: var(--panel-alt);
    color: var(--text);
    border: 1px solid var(--border);
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 11.5px;
  }}
  table.fov-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 12px;
    background: var(--panel);
    border-radius: 6px;
    overflow: hidden;
  }}
  table.fov-table th, table.fov-table td {{ padding: 8px 10px; text-align: left; border-bottom: 1px solid var(--border); }}
  table.fov-table th {{ background: var(--panel-alt); color: #e2e8f0; cursor: pointer; user-select: none; font-weight: 600; }}
  table.fov-table tr:hover {{ background: rgba(255, 255, 255, 0.04); cursor: pointer; }}
</style>
</head>
<body>

<header>
  <div class="title-block">
    <h1>GEO Stem Cell Differentiation Explorer</h1>
    <p>420 Fields of View &bull; 4 Lineages (hESC, Ectoderm, Endoderm, Mesoderm) &bull; 7 Timepoints (d0.0 &ndash; d3.0, 12 h intervals) &bull; 3 Replicate Wells</p>
  </div>
  <div class="mode-tabs">
    <button class="mode-btn active" onclick="setMode('single')" id="tab-single">Single FOV</button>
    <button class="mode-btn" onclick="setMode('timecourse')" id="tab-timecourse">7-Day Timecourse</button>
    <button class="mode-btn" onclick="setMode('lineages')" id="tab-lineages">4-Lineage Compare</button>
    <button class="mode-btn" onclick="setMode('montage')" id="tab-montage">All 15 FOVs in Condition</button>
    <button class="mode-btn" onclick="setMode('table')" id="tab-table">420-FOV Table</button>
    <button class="mode-btn" onclick="setMode('figures')" id="tab-figures">Summary Figures</button>
  </div>
</header>

<div class="controls-bar">
  <div class="ctrl-group">
    <span class="ctrl-label">Lineage</span>
    <button class="pill-btn treat-control" onclick="setTreat('control')" id="btn-treat-control">hESC Control</button>
    <button class="pill-btn treat-ecto" onclick="setTreat('ecto')" id="btn-treat-ecto">Ectoderm</button>
    <button class="pill-btn treat-endo" onclick="setTreat('endo')" id="btn-treat-endo">Endoderm</button>
    <button class="pill-btn treat-meso active" onclick="setTreat('meso')" id="btn-treat-meso">Mesoderm</button>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">Day</span>
    <button class="pill-btn" onclick="setDay(0.0)" id="btn-day-0.0">d0.0</button>
    <button class="pill-btn" onclick="setDay(0.5)" id="btn-day-0.5">d0.5</button>
    <button class="pill-btn" onclick="setDay(1.0)" id="btn-day-1.0">d1.0</button>
    <button class="pill-btn" onclick="setDay(1.5)" id="btn-day-1.5">d1.5</button>
    <button class="pill-btn" onclick="setDay(2.0)" id="btn-day-2.0">d2.0</button>
    <button class="pill-btn" onclick="setDay(2.5)" id="btn-day-2.5">d2.5</button>
    <button class="pill-btn active" onclick="setDay(3.0)" id="btn-day-3.0">d3.0</button>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">Well</span>
    <button class="pill-btn" onclick="setWell(0)" id="btn-well-0">Well 0</button>
    <button class="pill-btn active" onclick="setWell(1)" id="btn-well-1">Well 1</button>
    <button class="pill-btn" onclick="setWell(2)" id="btn-well-2">Well 2</button>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">FOV</span>
    <button class="pill-btn" onclick="setFov(0)" id="btn-fov-0">FOV 1</button>
    <button class="pill-btn" onclick="setFov(1)" id="btn-fov-1">FOV 2</button>
    <button class="pill-btn active" onclick="setFov(2)" id="btn-fov-2">FOV 3</button>
    <button class="pill-btn" onclick="setFov(3)" id="btn-fov-3">FOV 4</button>
    <button class="pill-btn" onclick="setFov(4)" id="btn-fov-4">FOV 5</button>
  </div>

  <div class="kb-hint">
    Shortcuts: <strong>1–4</strong> Lineage &middot; <strong>&larr;/&rarr;</strong> Day &middot; <strong>&uarr;/&darr;</strong> Well &middot; <strong>[ / ]</strong> FOV
  </div>
</div>

<div class="focus-bar">
  <span class="ctrl-label">View</span>
  <button class="pill-btn active" onclick="setFocus('all')" id="btn-focus-all">All 5 Panels</button>
  <button class="pill-btn" onclick="setFocus('pair')" id="btn-focus-pair">Overlay / Compare 2 Panels</button>
  <span style="color:var(--border);margin:0 4px;">|</span>
  <span class="ctrl-label">Single Panel</span>
  <button class="pill-btn" onclick="setFocus('t0')" id="btn-focus-t0">1. Fluorescence Still</button>
  <button class="pill-btn" onclick="setFocus('t1')" id="btn-focus-t1">2. FM Barcode f</button>
  <button class="pill-btn" onclick="setFocus('t2')" id="btn-focus-t2">3. Frequency Shift &Delta;F_stem</button>
  <button class="pill-btn" onclick="setFocus('t3')" id="btn-focus-t3">4. Competence Shift &Delta;P_osc</button>
  <button class="pill-btn" onclick="setFocus('t4')" id="btn-focus-t4">5. Stoichiometric Tilt &Delta;s_&perp;</button>

  <div class="ctrl-group" id="zoom-size-ctrl" style="margin-left:auto;">
    <span class="ctrl-label">Zoom Size</span>
    <input type="range" min="360" max="960" step="20" value="560" oninput="setZoomPx(this.value)" style="width:90px;">
    <span id="zoom-px-label" style="font-size:11.5px;color:var(--muted);">560px</span>
  </div>
</div>

<div class="focus-bar" id="pair-selector-bar" style="display:none; background:#111a2e; gap:12px;">
  <div class="ctrl-group">
    <span class="ctrl-label">Base (A)</span>
    <select class="pair-select" id="sel-pair-left" onchange="updatePair()">
      <option value="0" selected>1. Fluorescence Still (MinE-GFP)</option>
      <option value="1">2. FM Barcode f (18–52 mHz)</option>
      <option value="2">3. Frequency Shift &Delta;F_stem (mHz)</option>
      <option value="3">4. Competence Shift &Delta;P_osc (%)</option>
      <option value="4">5. Stoichiometric Tilt &Delta;s_&perp; (mHz)</option>
    </select>
    <button class="pill-btn" onclick="swapPair()">&harr; Swap</button>
    <span class="ctrl-label" style="margin-left:4px;">Overlay (B)</span>
    <select class="pair-select" id="sel-pair-right" onchange="updatePair()">
      <option value="0">1. Fluorescence Still (MinE-GFP)</option>
      <option value="1">2. FM Barcode f (18–52 mHz)</option>
      <option value="2" selected>3. Frequency Shift &Delta;F_stem (mHz)</option>
      <option value="3">4. Competence Shift &Delta;P_osc (%)</option>
      <option value="4">5. Stoichiometric Tilt &Delta;s_&perp; (mHz)</option>
    </select>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">Opacity</span>
    <input type="range" id="rng-opacity" min="0" max="100" step="1" value="65" oninput="setBlendOpacity(this.value)" style="width:120px;">
    <span id="lbl-opacity" style="font-size:11.5px;min-width:34px;">65%</span>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">Blend Mode</span>
    <select class="pair-select" id="sel-blend-mode" onchange="setBlendMode(this.value)">
      <option value="screen" selected>Screen</option>
      <option value="normal">Normal Alpha</option>
      <option value="hard-light">Hard-Light</option>
      <option value="overlay">Overlay</option>
    </select>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">Wipe</span>
    <input type="range" id="rng-wipe" min="0" max="100" step="1" value="100" oninput="setCurtainWipe(this.value)" style="width:90px;">
    <span id="lbl-wipe" style="font-size:11.5px;">100%</span>
  </div>
</div>

<main>
  <div id="view-container"></div>
</main>

<script>
const FOV_DATA = {json_payload};
const TREATMENTS = ['control', 'ecto', 'endo', 'meso'];
const DAYS = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0];

const PANEL_NAMES = [
  "1. Fluorescence Still (MinE-GFP)",
  "2. FM Barcode f (18–52 mHz)",
  "3. Frequency Shift ΔF_stem (mHz)",
  "4. Competence Shift ΔP_osc (%)",
  "5. Stoichiometric Tilt Δs_⊥ (mHz)"
];

let state = {{
  mode: 'single',
  treatment: 'meso',
  day: 3.0,
  well: 1,
  fov: 2,
  focus: 'all',
  zoomPx: 560,
  pairLeft: 0,
  pairRight: 2,
  blendOpacity: 65,
  blendMode: 'screen',
  curtainWipe: 100,
  sortCol: 'pct_meso_roc',
  sortAsc: false
}};

function findRecord(treat, day, well, fov) {{
  return FOV_DATA.find(r => r.treatment === treat && Math.abs(r.day - day) < 0.01 && r.well === well && r.fov_index === fov);
}}

function tileTransformCoords(tileIndex) {{
  const tx = -(tileIndex * 20.0);
  return {{ tx, ty: 0.0 }};
}}

function renderTileCrop(imgSrc, tileIndex, widthPx, labelText) {{
  const {{ tx }} = tileTransformCoords(tileIndex);
  const tagHtml = labelText ? `<div class="tile-col-tag">${{labelText}}</div>` : '';
  return `
    <div class="tile-col-wrap">
      ${{tagHtml}}
      <div class="tile-zoom-box" style="width:${{widthPx}}px;">
        <img class="base-layer" src="${{imgSrc}}" loading="lazy" style="transform: translate(${{tx.toFixed(2)}}%, 0%);">
      </div>
    </div>
  `;
}}

function renderBlendedTileCrop(imgSrc, leftIdx, rightIdx, widthPx) {{
  const cA = tileTransformCoords(leftIdx);
  const cB = tileTransformCoords(rightIdx);
  return `
    <div class="tile-col-wrap">
      <div class="tile-col-tag" style="border-color:#38bdf8;color:#f8fafc;">
        Blended: [${{PANEL_NAMES[leftIdx]}}] + [${{PANEL_NAMES[rightIdx]}}]
      </div>
      <div class="tile-zoom-box" style="width:${{widthPx}}px; border-color:#38bdf8;">
        <img class="base-layer" src="${{imgSrc}}" loading="lazy" style="transform: translate(${{cA.tx.toFixed(2)}}%, 0%);">
        <div class="overlay-clip">
          <img src="${{imgSrc}}" loading="lazy" style="transform: translate(${{cB.tx.toFixed(2)}}%, 0%);">
        </div>
      </div>
    </div>
  `;
}}

function renderVisualBody(r) {{
  if (state.focus === 'all') {{
    return `<div class="fov-img-wrap"><img src="${{r.img_rel_path}}" loading="lazy"></div>`;
  }} else if (state.focus === 'pair') {{
    return `
      <div class="single-tile-grid">
        ${{renderTileCrop(r.img_rel_path, state.pairLeft, state.zoomPx, 'Base A: ' + PANEL_NAMES[state.pairLeft])}}
        ${{renderBlendedTileCrop(r.img_rel_path, state.pairLeft, state.pairRight, state.zoomPx)}}
        ${{renderTileCrop(r.img_rel_path, state.pairRight, state.zoomPx, 'Overlay B: ' + PANEL_NAMES[state.pairRight])}}
      </div>
    `;
  }} else if (state.focus.startsWith('t')) {{
    const tIdx = parseInt(state.focus.substring(1), 10);
    return `<div class="single-tile-grid">${{renderTileCrop(r.img_rel_path, tIdx, Math.max(state.zoomPx, 640), PANEL_NAMES[tIdx])}}</div>`;
  }}
  return '';
}}

function renderCard(r) {{
  if (!r) return '<div class="fov-card">No FOV record found.</div>';
  return `
    <div class="fov-card">
      <div class="fov-card-header">
        <div>
          <span class="badge badge-${{r.treatment}}">${{r.treatment_label}}</span>
          <strong style="margin-left:8px;">Day ${{r.day.toFixed(1)}}</strong> &middot;
          Well ${{r.well}} &middot; FOV ${{r.fov_index + 1}} (<code>${{r.filename}}</code>)
        </div>
        <div style="color: var(--muted); font-size: 12px;">
          Oscillating cells: <strong style="color:#f8fafc">${{r.n_cells}}</strong> &nbsp;|&nbsp;
          Median f: <strong style="color:#38bdf8">${{r.med_freq_mHz.toFixed(1)}} mHz</strong> &nbsp;|&nbsp;
          &Delta;F_stem: <strong style="color:#f87171">${{r.med_df_stem_mHz >= 0 ? '+' : ''}}${{r.med_df_stem_mHz.toFixed(2)}} mHz</strong> &nbsp;|&nbsp;
          &Delta;P_osc: <strong style="color:#fbbf24">${{r.med_dposc_pct >= 0 ? '+' : ''}}${{r.med_dposc_pct.toFixed(1)}}%</strong> &nbsp;|&nbsp;
          Tilt &gt; +0.85 mHz: <strong style="color:#c084fc">${{r.pct_meso_roc.toFixed(1)}}% (n=${{r.n_meso_roc}})</strong>
        </div>
      </div>
      ${{renderVisualBody(r)}}
    </div>
  `;
}}

function updateUI() {{
  document.querySelectorAll('.mode-btn').forEach(b => b.classList.remove('active'));
  document.getElementById('tab-' + state.mode).classList.add('active');

  TREATMENTS.forEach(t => {{
    const el = document.getElementById('btn-treat-' + t);
    if (el) el.classList.toggle('active', state.treatment === t);
  }});
  DAYS.forEach(d => {{
    const el = document.getElementById('btn-day-' + d.toFixed(1));
    if (el) el.classList.toggle('active', Math.abs(state.day - d) < 0.01);
  }});
  [0, 1, 2].forEach(w => {{
    const el = document.getElementById('btn-well-' + w);
    if (el) el.classList.toggle('active', state.well === w);
  }});
  [0, 1, 2, 3, 4].forEach(f => {{
    const el = document.getElementById('btn-fov-' + f);
    if (el) el.classList.toggle('active', state.fov === f);
  }});

  ['all', 'pair', 't0', 't1', 't2', 't3', 't4'].forEach(fk => {{
    const el = document.getElementById('btn-focus-' + fk);
    if (el) el.classList.toggle('active', state.focus === fk);
  }});

  document.getElementById('pair-selector-bar').style.display = (state.focus === 'pair') ? 'flex' : 'none';

  const container = document.getElementById('view-container');

  if (state.mode === 'single') {{
    const r = findRecord(state.treatment, state.day, state.well, state.fov);
    let metricsHtml = '';
    if (r) {{
      metricsHtml = `
        <div class="metrics-strip">
          <div class="metric-card">
            <div class="m-label">Segmented Oscillating Cells</div>
            <div class="m-val">${{r.n_cells}}</div>
            <div class="m-sub">${{r.n_pixels.toLocaleString()}} strict oscillating pixels</div>
          </div>
          <div class="metric-card">
            <div class="m-label">Raw FM Barcode (f)</div>
            <div class="m-val" style="color:#38bdf8">${{r.med_freq_mHz.toFixed(1)}} mHz</div>
            <div class="m-sub">Within-FOV IQR: ${{r.iqr_freq_mHz.toFixed(1)}} mHz &middot; E/D: ${{r.med_ED_ratio.toFixed(2)}}</div>
          </div>
          <div class="metric-card">
            <div class="m-label">Frequency Shift (&Delta;F_stem)</div>
            <div class="m-val" style="color:#f87171">${{r.med_df_stem_mHz >= 0 ? '+' : ''}}${{r.med_df_stem_mHz.toFixed(2)}} mHz</div>
            <div class="m-sub">Demodulated vs. hESC f_ref(D, E)</div>
          </div>
          <div class="metric-card">
            <div class="m-label">Competence Shift (&Delta;P_osc)</div>
            <div class="m-val" style="color:#fbbf24">${{r.med_dposc_pct >= 0 ? '+' : ''}}${{r.med_dposc_pct.toFixed(1)}}%</div>
            <div class="m-sub">Relative to hESC P_ref(D, E) (&plusmn;13% IQR grey)</div>
          </div>
          <div class="metric-card">
            <div class="m-label">Stoichiometric Tilt (&Delta;s_&perp;)</div>
            <div class="m-val" style="color:#c084fc">${{r.pct_meso_roc.toFixed(1)}}% &gt; +0.85 mHz</div>
            <div class="m-sub">p90 = ${{r.p90_s_perp_roc_mHz >= 0 ? '+' : ''}}${{r.p90_s_perp_roc_mHz.toFixed(2)}} mHz (n = ${{r.n_meso_roc}} tail cells)</div>
          </div>
        </div>
      `;
    }}
    container.innerHTML = metricsHtml + renderCard(r);
  }} else if (state.mode === 'timecourse') {{
    let html = `<h3 style="margin-top:0;font-weight:600;font-size:15px;">7-Day Timecourse &mdash; ${{state.treatment.toUpperCase()}} (Well ${{state.well}}, FOV ${{state.fov + 1}})</h3>`;
    DAYS.forEach(d => {{
      html += renderCard(findRecord(state.treatment, d, state.well, state.fov));
    }});
    container.innerHTML = html;
  }} else if (state.mode === 'lineages') {{
    let html = `<h3 style="margin-top:0;font-weight:600;font-size:15px;">4-Lineage Comparison &mdash; Day ${{state.day.toFixed(1)}} (Well ${{state.well}}, FOV ${{state.fov + 1}})</h3>`;
    TREATMENTS.forEach(t => {{
      html += renderCard(findRecord(t, state.day, state.well, state.fov));
    }});
    container.innerHTML = html;
  }} else if (state.mode === 'montage') {{
    let html = `<h3 style="margin-top:0;font-weight:600;font-size:15px;">All 15 Replicate FOVs (3 Wells &times; 5 FOVs) &mdash; ${{state.treatment.toUpperCase()}} Day ${{state.day.toFixed(1)}}</h3>`;
    [0, 1, 2].forEach(w => {{
      [0, 1, 2, 3, 4].forEach(f => {{
        html += renderCard(findRecord(state.treatment, state.day, w, f));
      }});
    }});
    container.innerHTML = html;
  }} else if (state.mode === 'table') {{
    const sorted = [...FOV_DATA].sort((a, b) => {{
      const va = a[state.sortCol], vb = b[state.sortCol];
      if (typeof va === 'string') return state.sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
      return state.sortAsc ? va - vb : vb - va;
    }});
    let rows = sorted.map(r => `
      <tr onclick="jumpToFov('${{r.treatment}}', ${{r.day}}, ${{r.well}}, ${{r.fov_index}})">
        <td><span class="badge badge-${{r.treatment}}">${{r.treatment}}</span></td>
        <td><strong>${{r.day.toFixed(1)}}</strong></td>
        <td>${{r.well}}</td>
        <td>FOV ${{r.fov_index + 1}} (<code>${{r.fov_short}}</code>)</td>
        <td>${{r.n_cells}}</td>
        <td>${{r.med_ED_ratio.toFixed(2)}}</td>
        <td style="color:#38bdf8;">${{r.med_freq_mHz.toFixed(1)}}</td>
        <td style="color:#f87171;font-weight:600;">${{r.med_df_stem_mHz >= 0 ? '+' : ''}}${{r.med_df_stem_mHz.toFixed(2)}}</td>
        <td style="color:#fbbf24;font-weight:600;">${{r.med_dposc_pct >= 0 ? '+' : ''}}${{r.med_dposc_pct.toFixed(1)}}%</td>
        <td style="color:#c084fc;font-weight:600;">${{r.pct_meso_roc.toFixed(1)}}% (${{r.n_meso_roc}})</td>
        <td>${{r.p90_s_perp_roc_mHz >= 0 ? '+' : ''}}${{r.p90_s_perp_roc_mHz.toFixed(2)}}</td>
      </tr>
    `).join('');
    container.innerHTML = `
      <p style="color:var(--muted);margin-top:0;font-size:12.5px;">Click any column header to sort all 420 FOVs, or click any row to inspect that FOV.</p>
      <table class="fov-table">
        <thead>
          <tr>
            <th onclick="sortTable('treatment')">Lineage</th>
            <th onclick="sortTable('day')">Day</th>
            <th onclick="sortTable('well')">Well</th>
            <th onclick="sortTable('fov_index')">FOV</th>
            <th onclick="sortTable('n_cells')">Osc Cells</th>
            <th onclick="sortTable('med_ED_ratio')">Median E/D</th>
            <th onclick="sortTable('med_freq_mHz')">Median f (mHz)</th>
            <th onclick="sortTable('med_df_stem_mHz')">Median &Delta;F_stem (mHz)</th>
            <th onclick="sortTable('med_dposc_pct')">Median &Delta;P_osc (%)</th>
            <th onclick="sortTable('pct_meso_roc')">Tilt &gt; +0.85 mHz (%)</th>
            <th onclick="sortTable('p90_s_perp_roc_mHz')">Tilt p90 (mHz)</th>
          </tr>
        </thead>
        <tbody>${{rows}}</tbody>
      </table>
    `;
  }} else if (state.mode === 'figures') {{
    const figs = [
      {{
        title: "Figure 5 — GEO Physiology Encoders Resolve Single-Cell FM Barcodes and Distinguish Human Stem-Cell Differentiated States",
        sub: "Pluripotent hESC reference portraits f_hESC(D, E) and P_hESC(D, E), 2D phase-portrait comparisons across singularization (d0) and 3-day germ-layer endpoints (endoderm, ectoderm, mesoderm), and single-cell spatial overlays of raw FM barcodes (f), demodulated carrier frequency shifts (ΔF_stem), and competence shifts (ΔP_osc).",
        src: "figures/Figure5_GEO_Physiology_Encoding_Endpoints.png"
      }},
      {{
        title: "Figure 6 — Prior-Free 3-Parameter Physiological Trajectories (ΔF_stem, ΔP_osc, Δs_⊥) and Spatiotemporal Colony Dynamics Across Trilineage Differentiation",
        sub: "Independent 72-hour trajectories (ref → singularized d0 → d3.0 across n = 3 biological replicate wells per lineage), population distribution ridges, and living-colony overlays of ΔF_stem (Blue ↔ White ↔ Red), ΔP_osc (Teal ↔ Grey ↔ Gold), and stoichiometric tilt Δs_⊥ (Grey → Purple).",
        src: "figures/Figure6_GEO_Physiology_Timecourse_Trajectories.png"
      }}
    ];
    container.innerHTML = figs.map(fg => `
      <div class="fov-card">
        <div class="fov-card-header">
          <div><strong style="color:#f8fafc;font-size:14px;">${{fg.title}}</strong></div>
        </div>
        <p style="color:var(--muted);font-size:12.5px;margin:0 0 12px 0;">${{fg.sub}}</p>
        <div class="fov-img-wrap"><img src="${{fg.src}}" loading="lazy" style="background:#ffffff;"></div>
      </div>
    `).join('');
  }}
}}

function setMode(m) {{ state.mode = m; updateUI(); }}
function setTreat(t) {{ state.treatment = t; updateUI(); }}
function setDay(d) {{ state.day = d; updateUI(); }}
function setWell(w) {{ state.well = w; updateUI(); }}
function setFov(f) {{ state.fov = f; updateUI(); }}
function setFocus(fk) {{ state.focus = fk; updateUI(); }}
function setZoomPx(val) {{
  state.zoomPx = parseInt(val, 10);
  document.getElementById('zoom-px-label').textContent = state.zoomPx + 'px';
  updateUI();
}}
function setBlendOpacity(val) {{
  state.blendOpacity = parseInt(val, 10);
  document.getElementById('rng-opacity').value = String(state.blendOpacity);
  document.getElementById('lbl-opacity').textContent = state.blendOpacity + '%';
  document.documentElement.style.setProperty('--blend-opacity', (state.blendOpacity / 100.0).toFixed(2));
}}
function setBlendMode(val) {{
  state.blendMode = val;
  document.documentElement.style.setProperty('--blend-mode', val);
}}
function setCurtainWipe(val) {{
  state.curtainWipe = parseInt(val, 10);
  document.getElementById('lbl-wipe').textContent = state.curtainWipe + '%';
  // Image occupies [left=3.5%, right=15.5%] (width=81.0% of tile).
  const rightInset = 15.5 + (100 - state.curtainWipe) * 0.81;
  document.documentElement.style.setProperty('--wipe-right', rightInset.toFixed(2) + '%');
}}
function updatePair() {{
  state.pairLeft = parseInt(document.getElementById('sel-pair-left').value, 10);
  state.pairRight = parseInt(document.getElementById('sel-pair-right').value, 10);
  updateUI();
}}
function swapPair() {{
  const tmp = state.pairLeft;
  state.pairLeft = state.pairRight;
  state.pairRight = tmp;
  document.getElementById('sel-pair-left').value = String(state.pairLeft);
  document.getElementById('sel-pair-right').value = String(state.pairRight);
  updateUI();
}}
function sortTable(col) {{
  if (state.sortCol === col) state.sortAsc = !state.sortAsc;
  else {{ state.sortCol = col; state.sortAsc = false; }}
  updateUI();
}}
function jumpToFov(t, d, w, f) {{
  state.treatment = t; state.day = d; state.well = w; state.fov = f; state.mode = 'single';
  updateUI();
  window.scrollTo({{ top: 0, behavior: 'smooth' }});
}}

window.addEventListener('keydown', (e) => {{
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA' || e.target.tagName === 'SELECT') return;
  if (e.key === '1') setTreat('control');
  else if (e.key === '2') setTreat('ecto');
  else if (e.key === '3') setTreat('endo');
  else if (e.key === '4') setTreat('meso');
  else if (e.key === 'ArrowRight') {{
    let idx = DAYS.findIndex(d => Math.abs(d - state.day) < 0.01);
    if (idx < DAYS.length - 1) setDay(DAYS[idx + 1]);
  }} else if (e.key === 'ArrowLeft') {{
    let idx = DAYS.findIndex(d => Math.abs(d - state.day) < 0.01);
    if (idx > 0) setDay(DAYS[idx - 1]);
  }} else if (e.key === 'ArrowDown') {{
    setWell((state.well + 1) % 3);
    e.preventDefault();
  }} else if (e.key === 'ArrowUp') {{
    setWell((state.well + 2) % 3);
    e.preventDefault();
  }} else if (e.key === ']') {{
    setFov((state.fov + 1) % 5);
  }} else if (e.key === '[') {{
    setFov((state.fov + 4) % 5);
  }}
}});

updateUI();
</script>
</body>
</html>
"""

    local_html = LOCAL_WORK_DIR / "H9_FOV_Physiology_Atlas_Viewer.html"
    local_index = LOCAL_WORK_DIR / "index.html"
    local_html.write_text(html_content, encoding="utf-8")
    local_index.write_text(html_content, encoding="utf-8")
    local_json = LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json"
    local_json.write_text(json.dumps(records_sorted, indent=2), encoding="utf-8")
    pd.DataFrame(records_sorted).to_csv(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv", index=False)
    print(f"Saved clean 5-panel interactive HTML viewer: {local_index} & {local_html}", flush=True)


def main():
    t0 = time.perf_counter()
    ATLAS_DIR.mkdir(parents=True, exist_ok=True)
    if REMOTE_OUT_DIR.exists():
        REMOTE_ATLAS_DIR.mkdir(parents=True, exist_ok=True)

    bundle_path = str(ARCHIVE_DIR / "_atlas_precomputed_bundle.pkl")

    if "--html-only" in sys.argv:
        local_json = LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json"
        recs = json.loads(local_json.read_text(encoding="utf-8"))
        build_interactive_html_viewer(recs)
        return

    ckpt_files = [str(p) for p in sorted(WELL_CKPT_DIR.glob("*.parquet"))]

    if "--test-one" in sys.argv:
        _init_worker(bundle_path)
        test_ckpt = [cp for cp in ckpt_files if "meso_day_3.0_well_1" in cp][0]
        recs = process_single_well(test_ckpt)
        print(f"Test one well complete ({len(recs)} FOVs) in {time.perf_counter() - t0:.1f}s: {recs[2]['img_rel_path']}", flush=True)
        return

    print(f"Rendering clean 5-panel (1x5) atlas for all {len(ckpt_files)} wells (420 FOVs total)...", flush=True)

    all_fov_records = []
    done_wells = 0
    with ProcessPoolExecutor(max_workers=12, initializer=_init_worker, initargs=(bundle_path,)) as ex:
        futs = {ex.submit(process_single_well, cp): cp for cp in ckpt_files}
        for fut in as_completed(futs):
            recs = fut.result()
            all_fov_records.extend(recs)
            done_wells += 1
            if done_wells % 12 == 0 or done_wells == len(ckpt_files):
                print(
                    f"  [{done_wells:2d}/{len(ckpt_files)}] wells rendered ({len(all_fov_records)} FOVs done, elapsed {time.perf_counter() - t0:.1f}s)",
                    flush=True,
                )

    build_interactive_html_viewer(all_fov_records)

    if REMOTE_OUT_DIR.exists():
        print("Syncing HTML viewer and metadata to /mnt/researchdrive/General/Scott/stemcellcodes/vocoder_2026_10/...", flush=True)
        shutil.copy2(LOCAL_WORK_DIR / "H9_FOV_Physiology_Atlas_Viewer.html", REMOTE_OUT_DIR / "H9_FOV_Physiology_Atlas_Viewer.html")
        shutil.copy2(LOCAL_WORK_DIR / "index.html", REMOTE_OUT_DIR / "index.html")
        shutil.copy2(LOCAL_WORK_DIR / "README.md", REMOTE_OUT_DIR / "README.md")
        shutil.copy2(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv", REMOTE_OUT_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv")
        shutil.copy2(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json", REMOTE_OUT_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json")
        shutil.copy2(Path(__file__), REMOTE_OUT_DIR / Path(__file__).name)

    if GEO_REPO_DIR.exists():
        print("Syncing 420-FOV 5-panel atlas and HTML viewer to /home/smcoyle/GEO-stem-cell-differentiation-explorer/...", flush=True)
        shutil.copy2(LOCAL_WORK_DIR / "H9_FOV_Physiology_Atlas_Viewer.html", GEO_REPO_DIR / "H9_FOV_Physiology_Atlas_Viewer.html")
        shutil.copy2(LOCAL_WORK_DIR / "index.html", GEO_REPO_DIR / "index.html")
        shutil.copy2(LOCAL_WORK_DIR / "README.md", GEO_REPO_DIR / "README.md")
        shutil.copy2(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv", GEO_REPO_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv")
        shutil.copy2(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json", GEO_REPO_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json")
        shutil.copy2(Path(__file__), GEO_REPO_DIR / Path(__file__).name)
        for png_file in ATLAS_DIR.glob("*.png"):
            shutil.copy2(png_file, GEO_REPO_DIR / "fov_atlas" / png_file.name)

    print(f"Complete! All {len(all_fov_records)} FOVs rendered and synced in {time.perf_counter() - t0:.1f}s.", flush=True)


if __name__ == "__main__":
    main()
