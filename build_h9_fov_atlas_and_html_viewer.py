#!/usr/bin/env python3
r"""
Builds the Upgraded 12-Panel (2x6 Grid) 420-FOV Spatial Physiology, Grayscale Still, FM Barcode,
Meso-ROC & All 3 Stoichiometry-Matched Oscillation Probability Views (8A, 8B, 8C) Atlas +
Interactive HTML Viewer with Live Side-by-Side + Transparency/Blend Overlay Mode
(`H9_FOV_Physiology_Atlas_Viewer.html` & `index.html`).

Each of the 420 FOVs is rendered on a strict, pixel-exact 2-Row x 6-Column tile grid
(`28.8 x 10.5 in` @ `100 DPI` = `2880 x 1050 px`, dark theme `#0f172a`) where all 11 spatial panels
(`Tiles 0..10`) place their `512 x 512` colony image at the EXACT same relative sub-box
(`[x0 + 0.0055, y0 + 0.068, 0.1375, 0.376]`) inside their `(1/6 x 1/2)` tile!
This guarantees 0.00-pixel registration error when overlaying/blending any two transforms in the HTML viewer.

The 12 Panels per FOV are:
  ROW 1 (Fluorescence Stills, Stoichiometry Carrier, FM Barcode & Waveform Decomposition):
    0 -> 1. Grayscale Fluorescence Still (`*_g.tif`, `MinE-GFP` in easy-on-the-eyes Grayscale + Boundary Cue)
    1 -> 2. Dual Expression Still (`Red = mCherry-MinD`, `Green = MinE-GFP`, 100% of Colony)
    2 -> 3. Stoichiometry `E/D = E_smooth / D_smooth` Across Full Colony Footprint (`[2.0, 5.6]`)
    3 -> 4. Raw Oscillation Frequency `f_mHz` ("The FM Barcode", `[18.5, 49.5] mHz`)
    4 -> 5. Decoded Physiological State `ΔF_stem^surface` (`[-10.5, +10.5] mHz` vs. Control Homeostasis)
    5 -> 6. ROC-Calibrated `Meso/EMT` Axis `Δs_⊥^stem` (`τ* = +0.85 mHz` Grey-to-Pink Cutoff, `[-1.0, +2.8] mHz`)

  ROW 2 (Expected Competence, All 3 `ΔP_osc` Views `8A, 8B, 8C`, Dead-Zone Classes & 2D Sail):
    6 -> 7. Expected Control `P_osc^Ctrl(D, E)` From Still Alone (`[0%, 82%]`, with `P_ctrl = 20%` Gate)
    7 -> 8A. Surface Shift `ΔP_osc^surf(D, E) = P_cond(D, E) − P_ctrl(D, E)` (`[-65%, +65%]`, Non-Surprise)
    8 -> 8B. Local Gated `ΔP_osc^stem(y, x) = P_local − P_ctrl(D, E)` (`[-65%, +65%]`, Blue / Grey / Gold ±12%)
    9 -> 8C. Dual-ROC Calibrated Surprise (`Blue '0-on-1' Quench / Flat Grey '1-on-1' / Vivid Pink '1-on-0'`)
   10 -> 9. Full-Colony 5-Class Dead-Zone Diagnostic Map (`Z-Filtered`, `High-E/D`, `Mid-E/D Damped`, `Low-E/D`, `Dim`)
   11 -> 10. Live FOV Single-Cell Cloud in Stem-Cell-Centered `(Δs_∥^stem, Δs_⊥^stem)` Sail Space (Colored by 8C Surprise)
"""

import io
import json
import os
import pickle
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm, LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_dilation, mean as ndi_mean
import tifffile

LOCAL_WORK_DIR = Path("/home/smcoyle/h9_vocoder_work")
ARCHIVE_DIR = LOCAL_WORK_DIR / "workspace_archive"
sys.path.insert(0, str(LOCAL_WORK_DIR))
sys.path.insert(0, str(ARCHIVE_DIR))

from plot_sail_flux_4zones_and_fov_remapping import (
    WORK_DIR, MASK_ROOT, COS_T, SIN_T, TREATMENTS, DAYS, TREAT_NAMES, DAY_FOLDER_MAP
)
from plot_fig36_direct_surface_df_stem_and_meso_axis_fovs import TREAT_ACCENT, make_df_stem_clean_cmap
from plot_fig41_4channel_green_still_fm_barcode_stateshift_and_roc_meso import (
    make_fm_barcode_cmap, make_roc_calibrated_meso_cmap
)
from plot_fig42_barcode_transcendence_and_deadzone_diagnostics import (
    make_ed_stoich_cmap, render_composite_rg_still, render_deadzone_diagnostic_map
)
from plot_fig43_stoichiometry_matched_oscillation_probability import (
    make_p_ctrl_cmap, make_delta_p_osc_cmap, load_cond_pixels_fast
)
from plot_fig47_dual_roc_surprise_flat_control_blue_quench_pink_meso import (
    make_dual_roc_surprise_blue_grey_pink_cmap, compute_unified_roc_surprise_map
)

WELL_CKPT_DIR = LOCAL_WORK_DIR / "well_checkpoints"
LOCAL_MASK_DIR = LOCAL_WORK_DIR / "consolidated_masks"
STILLS_DIR = LOCAL_WORK_DIR / "consolidated_stills"
ATLAS_DIR = LOCAL_WORK_DIR / "fov_atlas"

REMOTE_OUT_DIR = Path("/mnt/researchdrive/General/Scott/stemcellcodes/vocoder_2026_10")
REMOTE_ATLAS_DIR = REMOTE_OUT_DIR / "fov_atlas"
GEO_REPO_DIR = Path("/home/smcoyle/GEO-stem-cell-differentiation-explorer")

TREAT_LABELS = {
    "control": "Pluripotent Control (mTeSR)",
    "ecto": "Ectoderm Differentiation",
    "endo": "Endoderm Differentiation",
    "meso": "Mesoderm Differentiation",
}

# Worker globals initialized once per process
_G_CELL_DICT = None
_G_P_CTRL_2D = None
_G_D_EDGES = None
_G_E_EDGES = None
_G_SAIL_BG = None
_G_DP_SURF_28 = None
_G_D_EDGES_34 = None
_G_E_EDGES_34 = None


def make_grayscale_fluorescence_cmap():
    """Soft, eye-friendly pure grayscale colormap for MinE-GFP fluorescence stills."""
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


def render_grayscale_still_rgb(g_im: np.ndarray, mask: np.ndarray, step: int = 2) -> np.ndarray:
    """Renders `*_g.tif` (`MinE-GFP`) in clean, eye-friendly pure grayscale with a subtle boundary cue."""
    g_sub = g_im[::step, ::step].astype(np.float32)
    m_sub = mask[::step, ::step]
    p_lo, p_hi = np.percentile(g_sub, 2.0), np.percentile(g_sub, 99.2)
    g_norm = np.clip((g_sub - p_lo) / max(p_hi - p_lo, 1.0), 0.0, 1.0) ** 0.80
    # Map [0, 1] to [0.04, 0.96] neutral grayscale so dark background matches the other panels
    lum = 0.04 + 0.92 * g_norm
    rgb = np.stack([lum, lum, lum], axis=-1).astype(np.float32)
    has_seg = m_sub > 0
    edge = binary_dilation(has_seg, iterations=1) & (~has_seg)
    rgb[edge] = np.array([0.46, 0.50, 0.56], dtype=np.float32)
    return rgb


def build_precomputed_state_bundle() -> str:
    """
    Loads or builds `_atlas_precomputed_bundle.pkl` in `ARCHIVE_DIR` (or `LOCAL_WORK_DIR`).
    """
    arch_bundle = ARCHIVE_DIR / "_atlas_precomputed_bundle.pkl"
    local_bundle = LOCAL_WORK_DIR / "_atlas_precomputed_bundle.pkl"
    if arch_bundle.exists():
        return str(arch_bundle)
    if local_bundle.exists():
        return str(local_bundle)

    print("Building precomputed single-cell & 2D P_ctrl(D, E) reference bundle for 420-FOV atlas...", flush=True)
    df1 = pd.read_parquet(
        ARCHIVE_DIR / "2026_10_20_stemcell_all_single_cells_vocoder_physiology.parquet",
        columns=[
            "treatment", "day", "well", "well_num", "filename", "mask_filename", "label_id",
            "cell_area_px", "n_osc_pixels", "osc_fraction",
            "D_smooth_mean", "E_smooth_mean", "ED_ratio_mean",
            "freq_mHz_mean", "delta_f_day0_mHz_mean", "px_nopr_voc_std"
        ],
    )
    df2 = pd.read_parquet(ARCHIVE_DIR / "2026_10_24_single_cell_posterior_confidence_audit.parquet")
    df = df2.merge(df1, on=["treatment", "day", "well_num", "filename", "label_id"], suffixes=("", "_orig"))

    off_raw = df["cell_Offset_mHz"].to_numpy()
    tilt_raw = df["cell_Tilt_mHz"].to_numpy()
    df["s_par_raw"] = -COS_T * off_raw - SIN_T * tilt_raw
    df["s_perp_raw"] = -SIN_T * off_raw + COS_T * tilt_raw

    g_strict = (
        (df["n_osc_pixels"] >= 35)
        & (df["in_pivot_frac"] <= 0.15)
        & (df["cnt0_ge50_frac"] >= 0.80)
        & (df["chi2_min"] <= 0.85)
        & (df["px_nopr_voc_std"] <= 0.20)
    )
    df["is_hc"] = g_strict
    ctrl_homeo_hc = df[g_strict & (df["treatment"] == "control") & (df["day"] >= 1.0)]
    S_PAR_STAR = float(ctrl_homeo_hc["s_par_raw"].median())
    S_PERP_STAR = float(ctrl_homeo_hc["s_perp_raw"].median())
    df["ds_par"] = df["s_par_raw"] - S_PAR_STAR
    df["ds_perp"] = df["s_perp_raw"] - S_PERP_STAR

    w_pix = np.clip((df["n_osc_pixels"].to_numpy() - 5.0) / 30.0, 0.15, 1.0)
    w_chi = np.clip(1.0 / (1.0 + (df["chi2_min"].to_numpy() / 0.85) ** 2), 0.15, 1.0)
    w_piv = np.clip(1.0 - df["in_pivot_frac"].to_numpy(), 0.20, 1.0)
    df["ds_perp_roc"] = df["ds_perp"] * np.sqrt(w_pix * w_chi) * w_piv

    ed_arr = df["ED_ratio_mean"].to_numpy()
    ed_lever = np.clip((np.log2(np.maximum(ed_arr, 0.25)) - 1.0) / 1.8, -1.2, 1.2)
    df["df0_clean"] = np.where(np.isfinite(df["df0_mean"].to_numpy()), df["df0_mean"].to_numpy(), off_raw + tilt_raw * ed_lever)

    ctrl_homeo_all = df[(df["treatment"] == "control") & (df["day"] >= 1.0) & (df["day"] <= 3.0)]
    ed_grid = np.linspace(1.2, 8.5, 140)
    ed_centers = 0.5 * (ed_grid[:-1] + ed_grid[1:])
    c_ed = ctrl_homeo_all["ED_ratio_mean"].to_numpy()
    c_df0 = ctrl_homeo_all["df0_clean"].to_numpy()
    ctrl_ed_medians = np.array([
        np.nanmedian(c_df0[(c_ed >= ed_grid[i]) & (c_ed < ed_grid[i + 1])])
        if np.sum((c_ed >= ed_grid[i]) & (c_ed < ed_grid[i + 1])) >= 25 else np.nan
        for i in range(len(ed_centers))
    ])
    ctrl_ed_smooth = gaussian_filter1d(pd.Series(ctrl_ed_medians).interpolate().bfill().ffill().to_numpy(), sigma=1.5)
    df["df_stem_surface"] = df["df0_clean"] - np.interp(ed_arr, ed_centers, ctrl_ed_smooth)

    d_edges = np.geomspace(155.0, 950.0, 33)
    e_edges = np.geomspace(460.0, 3200.0, 33)
    dc_l, ec_l, stc_l = [], [], []
    for dy_c in [1.5, 2.0, 2.5, 3.0]:
        d, e, st = load_cond_pixels_fast("control", dy_c, max_fovs=4, step=2)
        dc_l.append(d)
        ec_l.append(e)
        stc_l.append(st)
    Dc = np.concatenate(dc_l)
    Ec = np.concatenate(ec_l)
    STc = np.concatenate(stc_l)
    H_tot_c, _, _ = np.histogram2d(Dc, Ec, bins=[d_edges, e_edges])
    H_st_c, _, _ = np.histogram2d(Dc, Ec, bins=[d_edges, e_edges], weights=STc.astype(float))
    P_ctrl_raw = np.where(H_tot_c >= 25, H_st_c / np.maximum(H_tot_c, 1), 0.0)
    P_ctrl_2d = gaussian_filter(P_ctrl_raw, sigma=0.85)

    di_c = np.clip(np.searchsorted(d_edges, df["D_smooth_mean"].to_numpy()) - 1, 0, len(d_edges) - 2)
    ei_c = np.clip(np.searchsorted(e_edges, df["E_smooth_mean"].to_numpy()) - 1, 0, len(e_edges) - 2)
    df["P_ctrl_cell"] = 100.0 * P_ctrl_2d[di_c, ei_c]
    df["osc_pct_cell"] = 100.0 * np.clip(df["osc_fraction"].to_numpy(), 0.0, 1.0)
    df["dP_osc_cell"] = df["osc_pct_cell"] - df["P_ctrl_cell"]

    hc_sub = df[df["is_hc"]].sample(n=min(12000, int(df["is_hc"].sum())), random_state=42)
    sail_bg = (
        hc_sub["ds_par"].to_numpy(dtype=np.float32),
        hc_sub["ds_perp"].to_numpy(dtype=np.float32),
    )

    cell_dict = {}
    for (tr, dy, wn, fn), g in df.groupby(["treatment", "day", "well_num", "filename"]):
        cell_dict[(str(tr), float(dy), int(wn), str(fn))] = {
            "label_id": g["label_id"].to_numpy(dtype=np.int32),
            "freq_mHz": g["freq_mHz_mean"].to_numpy(dtype=np.float32),
            "df_stem": g["df_stem_surface"].to_numpy(dtype=np.float32),
            "ds_par": g["ds_par"].to_numpy(dtype=np.float32),
            "ds_perp": g["ds_perp"].to_numpy(dtype=np.float32),
            "ds_perp_roc": g["ds_perp_roc"].to_numpy(dtype=np.float32),
            "dP_osc_cell": g["dP_osc_cell"].to_numpy(dtype=np.float32),
            "ED_ratio": g["ED_ratio_mean"].to_numpy(dtype=np.float32),
        }

    with open(arch_bundle, "wb") as f:
        pickle.dump(
            {
                "cell_dict": cell_dict,
                "P_ctrl_2d": P_ctrl_2d,
                "d_edges": d_edges,
                "e_edges": e_edges,
                "sail_bg": sail_bg,
            },
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    return str(arch_bundle)


def _init_worker(bundle_path_str: str):
    global _G_CELL_DICT, _G_P_CTRL_2D, _G_D_EDGES, _G_E_EDGES, _G_SAIL_BG
    global _G_DP_SURF_28, _G_D_EDGES_34, _G_E_EDGES_34
    with open(bundle_path_str, "rb") as f:
        b = pickle.load(f)
    _G_CELL_DICT = b["cell_dict"]
    _G_P_CTRL_2D = b["P_ctrl_2d"]
    _G_D_EDGES = b["d_edges"]
    _G_E_EDGES = b["e_edges"]
    _G_SAIL_BG = b["sail_bg"]

    surf_pkl = ARCHIVE_DIR / "_fig46_28cond_dual_surfaces.pkl"
    with open(surf_pkl, "rb") as f:
        s = pickle.load(f)
    P_ref_34 = s["P_ref"]
    cond_raw_34 = s["cond_raw"]
    _G_D_EDGES_34 = s["d_edges"]
    _G_E_EDGES_34 = s["e_edges"]

    _G_DP_SURF_28 = {}
    for (tr, dy), c_dict in cond_raw_34.items():
        p_raw = c_dict["P_surf"]
        p_filled = np.where(np.isfinite(p_raw), p_raw, P_ref_34)
        _G_DP_SURF_28[(str(tr), float(dy))] = (100.0 * (p_filled - P_ref_34)).astype(np.float32)


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
    cmap_ed = make_ed_stoich_cmap()
    cmap_barcode = make_fm_barcode_cmap()
    cmap_df = make_df_stem_clean_cmap()
    cmap_meso_roc = make_roc_calibrated_meso_cmap()
    cmap_pctrl = make_p_ctrl_cmap()
    cmap_dp = make_delta_p_osc_cmap()
    cmap_roc_surp = make_dual_roc_surprise_blue_grey_pink_cmap()

    norm_ed = Normalize(vmin=2.0, vmax=5.6)
    norm_barcode = Normalize(vmin=18.5, vmax=49.5)
    norm_df = Normalize(vmin=-10.5, vmax=10.5)
    norm_meso_roc = Normalize(vmin=-1.0, vmax=2.8)
    norm_pctrl = Normalize(vmin=0.0, vmax=82.0)
    norm_dp = TwoSlopeNorm(vmin=-65.0, vcenter=0.0, vmax=65.0)
    norm_roc_surp = TwoSlopeNorm(vmin=-100.0, vcenter=0.0, vmax=100.0)

    bg_col = np.array([0.04, 0.06, 0.10], dtype=np.float32)
    z0_col = np.array([0.09, 0.12, 0.18], dtype=np.float32)

    filenames = sorted(df_w["filename"].unique())
    fov_records = []

    dp_surf_grid = _G_DP_SURF_28.get((treat, day_val), None)

    for fov_idx, fn in enumerate(filenames):
        fov_short = fn.replace(".nd2", "")
        mp = LOCAL_MASK_DIR / day_folder / cond_folder / f"{fov_short}_masks.tif"
        gp = STILLS_DIR / day_folder / cond_folder / f"{fov_short}_g.tif"
        rp = STILLS_DIR / day_folder / cond_folder / f"{fov_short}_r.tif"
        if not (mp.exists() and gp.exists() and rp.exists()):
            continue

        mask = tifffile.imread(mp)
        g_im = tifffile.imread(gp)
        r_im = tifffile.imread(rp)

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
        relaxed_yx = df_w_fov[["y", "x"]].to_numpy()

        step = 2
        g_sm = gaussian_filter(g_im.astype(np.float32), sigma=2.0)[::step, ::step]
        r_sm = gaussian_filter(r_im.astype(np.float32), sigma=2.0)[::step, ::step]
        m_sub = mask[::step, ::step]

        strict_mask = np.zeros(mask.shape, dtype=bool)
        if len(strict_yx) > 0:
            strict_mask[strict_yx[:, 0], strict_yx[:, 1]] = True
        st_sub = strict_mask[::step, ::step]

        # Dead-Zone Diagnostic Map (+ in_tissue & ed_sub)
        rgb_dead, dz_stats, in_tissue, ed_sub = render_deadzone_diagnostic_map(
            mask, r_im, g_im, strict_yx, relaxed_yx, step=step
        )

        # 1. Grayscale Still (*_g.tif in clean Grayscale)
        rgb_gray = render_grayscale_still_rgb(g_im, mask, step=step)

        # 2. Dual Expression Still (Red=MinD, Green=MinE)
        rgb_rg = render_composite_rg_still(r_im, g_im, mask, step=step)

        # 3. Stoichiometry E/D Across Full Colony Footprint
        labels = np.unique(mask)
        labels = labels[labels > 0]
        max_lbl = int(mask.max())
        if len(labels) > 0:
            gm_all = ndi_mean(g_im.astype(np.float64), labels=mask, index=labels)
            rm_all = ndi_mean(r_im.astype(np.float64), labels=mask, index=labels)
            ed_all = gm_all / np.maximum(rm_all, 1.0)
        else:
            ed_all = np.array([], dtype=np.float32)

        lut_ed = np.zeros((max_lbl + 1, 3), dtype=np.float32)
        lut_ed[0] = bg_col
        if len(labels) > 0:
            lut_ed[labels] = cmap_ed(norm_ed(ed_all))[:, :3].astype(np.float32)
        rgb_ed = lut_ed[m_sub]
        unmasked_tissue = in_tissue & (m_sub == 0)
        if np.any(unmasked_tissue):
            rgb_ed[unmasked_tissue] = cmap_ed(norm_ed(ed_sub[unmasked_tissue]))[:, :3].astype(np.float32)

        # Single-cell lookup for Tiles 4, 5, 6, 8C, 10
        c_info = _G_CELL_DICT.get((treat, day_val, well_val, fn), None)
        lut_bar = np.zeros((max_lbl + 1, 3), dtype=np.float32)
        lut_bar[0] = bg_col
        lut_df = np.zeros((max_lbl + 1, 3), dtype=np.float32)
        lut_df[0] = bg_col
        lut_meso = np.zeros((max_lbl + 1, 3), dtype=np.float32)
        lut_meso[0] = bg_col
        if len(labels) > 0:
            lut_bar[labels] = z0_col
            lut_df[labels] = z0_col
            lut_meso[labels] = z0_col

        if c_info is not None and len(c_info["label_id"]) > 0:
            c_ids = c_info["label_id"]
            valid_c = (c_ids > 0) & (c_ids <= max_lbl)
            c_ids = c_ids[valid_c]
            f_vals = c_info["freq_mHz"][valid_c]
            df_vals = c_info["df_stem"][valid_c]
            sp_vals = c_info["ds_perp_roc"][valid_c]
            spar_vals = c_info["ds_par"][valid_c]
            sperp_raw_vals = c_info["ds_perp"][valid_c]
            dp_cell_vals = c_info["dP_osc_cell"][valid_c]

            lut_bar[c_ids] = cmap_barcode(norm_barcode(f_vals))[:, :3].astype(np.float32)
            lut_df[c_ids] = cmap_df(norm_df(df_vals))[:, :3].astype(np.float32)
            lut_meso[c_ids] = cmap_meso_roc(norm_meso_roc(sp_vals))[:, :3].astype(np.float32)
        else:
            c_ids = np.array([], dtype=np.int32)
            f_vals = np.array([0.0], dtype=np.float32)
            df_vals = np.array([0.0], dtype=np.float32)
            sp_vals = np.array([0.0], dtype=np.float32)
            spar_vals = np.array([], dtype=np.float32)
            sperp_raw_vals = np.array([], dtype=np.float32)
            dp_cell_vals = np.array([], dtype=np.float32)

        rgb_bar = lut_bar[m_sub]
        rgb_df = lut_df[m_sub]
        rgb_meso = lut_meso[m_sub]
        if np.any(unmasked_tissue):
            rgb_bar[unmasked_tissue] = z0_col
            rgb_df[unmasked_tissue] = z0_col
            rgb_meso[unmasked_tissue] = z0_col

        # 7. Expected Control P_osc^Ctrl(D, E)
        di = np.clip(np.searchsorted(_G_D_EDGES, r_sm) - 1, 0, len(_G_D_EDGES) - 2)
        ei = np.clip(np.searchsorted(_G_E_EDGES, g_sm) - 1, 0, len(_G_E_EDGES) - 2)
        p_exp_map = 100.0 * _G_P_CTRL_2D[di, ei]

        rgb_pexp = np.zeros((*m_sub.shape, 3), dtype=np.float32)
        rgb_pexp[:, :] = bg_col
        rgb_pexp[in_tissue] = cmap_pctrl(norm_pctrl(p_exp_map[in_tissue]))[:, :3].astype(np.float32)

        # 8A. Surface Shift ΔP_osc^surf(D, E) = P_cond(D, E) - P_ctrl(D, E)
        di34 = np.clip(np.searchsorted(_G_D_EDGES_34, r_sm) - 1, 0, len(_G_D_EDGES_34) - 2)
        ei34 = np.clip(np.searchsorted(_G_E_EDGES_34, g_sm) - 1, 0, len(_G_E_EDGES_34) - 2)
        dp_surf_map = dp_surf_grid[di34, ei34] if dp_surf_grid is not None else np.zeros_like(p_exp_map)
        comp_surf = in_tissue & (p_exp_map >= 20.0)
        z0_surf = in_tissue & (~comp_surf)
        rgb_8a = np.zeros((*m_sub.shape, 3), dtype=np.float32)
        rgb_8a[:, :] = bg_col
        rgb_8a[z0_surf] = z0_col
        rgb_8a[comp_surf] = cmap_dp(norm_dp(dp_surf_map[comp_surf]))[:, :3].astype(np.float32)
        med_dp_surf = float(np.nanmedian(dp_surf_map[comp_surf])) if np.any(comp_surf) else 0.0

        # 8B. Local Gated ΔP_osc^stem(y, x) = P_local - P_ctrl(D, E) (Blue/Grey/Gold ±12%)
        p_local_map = 100.0 * gaussian_filter(st_sub.astype(np.float32), sigma=6.0)
        dp_local_map = p_local_map - p_exp_map
        is_competent_or_active = in_tissue & ((p_exp_map >= 20.0) | (p_local_map >= 15.0))
        is_zero_on_zero = in_tissue & (~is_competent_or_active)

        rgb_8b = np.zeros((*m_sub.shape, 3), dtype=np.float32)
        rgb_8b[:, :] = bg_col
        rgb_8b[is_zero_on_zero] = z0_col
        rgb_8b[is_competent_or_active] = cmap_dp(norm_dp(dp_local_map[is_competent_or_active]))[:, :3].astype(np.float32)

        # 8C. Dual-ROC Calibrated Surprise (Blue '0-on-1' / Flat Grey '1-on-1' / Pink '1-on-0')
        s_roc_map, z0_roc_mask, active_osc_mask = compute_unified_roc_surprise_map(
            st_sub, m_sub, in_tissue, p_exp_map, c_info, max_lbl
        )
        valid_roc = in_tissue & (~z0_roc_mask)
        rgb_8c = np.zeros((*m_sub.shape, 3), dtype=np.float32)
        rgb_8c[:, :] = bg_col
        rgb_8c[z0_roc_mask] = z0_col
        rgb_8c[valid_roc] = cmap_roc_surp(norm_roc_surp(s_roc_map[valid_roc]))[:, :3].astype(np.float32)

        # Compute summary metrics for this FOV
        n_cells = int(len(c_ids))
        n_seg_total = int(len(labels))
        osc_tissue_pct = float(dz_stats["osc_pct"])
        exp_pctrl_pct = float(np.mean(p_exp_map[in_tissue])) if np.any(in_tissue) else 0.0
        comp_mask = in_tissue & (p_exp_map >= 20.0)
        competent_pct = 100.0 * float(np.mean(comp_mask[in_tissue])) if np.any(in_tissue) else 0.0
        gated_dposc_med = float(np.median(dp_local_map[comp_mask])) if np.any(comp_mask) else 0.0
        quenched_pct = 100.0 * float(np.mean(dp_local_map[comp_mask] < -20.0)) if np.any(comp_mask) else 0.0
        more_alive_pct = 100.0 * float(np.mean(dp_local_map[is_competent_or_active] > 12.0)) if np.any(is_competent_or_active) else 0.0

        roc_grey_pct = 100.0 * float(np.mean(np.abs(s_roc_map[valid_roc]) < 28.5)) if np.any(valid_roc) else 0.0
        roc_blue_pct = 100.0 * float(np.mean(s_roc_map[valid_roc] <= -28.5)) if np.any(valid_roc) else 0.0
        roc_pink_pct = 100.0 * float(np.mean(s_roc_map[active_osc_mask] >= 28.5)) if np.any(active_osc_mask) else 0.0

        ed_tissue_vals = ed_sub[in_tissue] if np.any(in_tissue) else np.array([0.0])
        med_ed = float(np.median(ed_tissue_vals))
        iqr_ed = float(np.percentile(ed_tissue_vals, 75) - np.percentile(ed_tissue_vals, 25))
        med_f = float(np.median(f_vals))
        iqr_f = float(np.percentile(f_vals, 75) - np.percentile(f_vals, 25)) if len(f_vals) > 1 else 0.0
        med_df_stem = float(np.median(df_vals))
        med_s_par = float(np.median(spar_vals)) if len(spar_vals) > 0 else 0.0
        med_s_perp = float(np.median(sp_vals))
        n_meso_roc = int(np.sum(sp_vals > 0.85)) if len(sp_vals) > 0 else 0
        pct_meso_roc = 100.0 * n_meso_roc / max(n_cells, 1)

        rel_png = f"fov_atlas/{treat}_d{day_val:.1f}_w{well_val}_fov{fov_idx}_{fov_short}.png"
        out_png = LOCAL_WORK_DIR / rel_png

        # =====================================================================
        # RENDER STRICT 2x6 UNIFORM TILE GRID (28.8 x 10.5 in @ 100 DPI = 2880 x 1050 px)
        # Every Tile (r in {0,1}, c in {0..5}) occupies the exact rectangle:
        #   x in [c/6, (c+1)/6], y in [(1-r)*0.50, (2-r)*0.50]
        # Every spatial image (Tiles 0..10) uses the EXACT same w_ax=0.1375, h_ax=0.376
        # at [x0 + 0.0055, y0 + 0.068] so any two tiles blend with 0.00 px shift!
        # =====================================================================
        fig = plt.figure(figsize=(28.8, 10.5), dpi=100)
        fig.patch.set_facecolor("#0f172a")

        fov_pt = fov_short.split("_")[0]
        rgb_tiles = [
            # Row 0 (Top Row: Grayscale Still, Dual Still, Stoichiometry, FM Barcode, Decoded ΔF_stem, ROC Meso Δs_⊥)
            (
                0, 0, rgb_gray, cmap_gray, Normalize(0, 100),
                "1. Grayscale Still (*_g.tif : MinE-GFP)",
                f"{treat.upper()} d{day_val:.1f} w{well_val} ({fov_pt}) | OscPx: {osc_tissue_pct:.0f}%",
                "#e2e8f0"
            ),
            (
                0, 1, rgb_rg, None, None,
                "2. Dual Still (Red=MinD, Green=MinE)",
                f"Osc Cells: {n_cells} / {n_seg_total} Seg | Competent: {competent_pct:.0f}%",
                "#34d399"
            ),
            (
                0, 2, rgb_ed, cmap_ed, norm_ed,
                "3. Stoichiometry E/D (Full Colony)",
                f"Tissue E/D: {med_ed:.2f} (IQR: {iqr_ed:.2f})",
                "#38bdf8"
            ),
            (
                0, 3, rgb_bar, cmap_barcode, norm_barcode,
                "4. Raw FM Barcode f_mHz (18.5–49.5)",
                f"Med f: {med_f:.1f} mHz (IQR: {iqr_f:.1f} mHz)",
                "#fbbf24"
            ),
            (
                0, 4, rgb_df, cmap_df, norm_df,
                "5. Decoded State Shift ΔF_stem (mHz)",
                f"Med ΔF_stem: {med_df_stem:+.2f} mHz vs. Control",
                "#60a5fa"
            ),
            (
                0, 5, rgb_meso, cmap_meso_roc, norm_meso_roc,
                "6. ROC Meso/EMT Axis Δs_⊥^stem (τ*=+0.85)",
                f"Meso+ (>+0.85): {pct_meso_roc:.1f}% (n={n_meso_roc}) | Med: {med_s_perp:+.2f}",
                "#f43f5e"
            ),
            # Row 1 (Bottom Row: Expected P_ctrl, 8A Surface ΔP_osc, 8B Local Gated ΔP_osc, 8C Dual-ROC Surprise, 9 Dead-Zone Map)
            (
                1, 0, rgb_pexp, cmap_pctrl, norm_pctrl,
                "7. Expected Control P_osc^Ctrl(D,E)",
                f"Expected P_ctrl: {exp_pctrl_pct:.1f}% (P_ctrl≥20% Area: {competent_pct:.0f}%)",
                "#fde047"
            ),
            (
                1, 1, rgb_8a, cmap_dp, norm_dp,
                "8A. Surface Shift ΔP_osc^surf (Non-Surprise)",
                f"Population Surface ΔP_osc(D,E): {med_dp_surf:+.1f}% vs. Ctrl",
                "#7dd3fc"
            ),
            (
                1, 2, rgb_8b, cmap_dp, norm_dp,
                "8B. Local Gated ΔP_osc (Blue/Grey/Gold ±12%)",
                f"Quenched (Blue): {quenched_pct:.0f}% | +Alive (Gold): {more_alive_pct:.0f}%",
                "#fb923c"
            ),
            (
                1, 3, rgb_8c, cmap_roc_surp, norm_roc_surp,
                "8C. Dual-ROC Surprise (Blue / Flat Grey / Pink)",
                f"Grey: {roc_grey_pct:.0f}% | Blue Quench: {roc_blue_pct:.0f}% | Pink Cells: {roc_pink_pct:.0f}%",
                "#f472b6"
            ),
            (
                1, 4, rgb_dead, None, None,
                "9. Full-Colony 5-Class Dead-Zone Map",
                f"HiED(Cyan) {dz_stats['hied_pct']:.0f}% | ZFilt(Gold) {dz_stats['zfilt_pct']:.0f}% | Dim(Red) {dz_stats['dim_pct']:.0f}%",
                "#c084fc"
            ),
        ]

        w_cell = 1.0 / 6.0
        for r_t, c_t, rgb_arr, cm, norm_obj, title_s, sub_s, t_col in rgb_tiles:
            x0 = c_t * w_cell
            y0 = (1 - r_t) * 0.50
            # IDENTICAL image axes rectangle across all 11 spatial tiles!
            ax = fig.add_axes([x0 + 0.0055, y0 + 0.068, 0.1375, 0.376])
            ax.imshow(rgb_arr, interpolation="nearest", aspect="auto")
            ax.set_title(title_s, fontsize=9.1, fontweight="bold", color=t_col, pad=4)
            ax.set_xlabel(sub_s, fontsize=7.5, fontweight="bold", color="#e2e8f0", labelpad=4)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_edgecolor(t_col if (r_t == 0 and c_t in (4, 5)) or (r_t == 1 and c_t in (1, 2, 3)) else "#334155")
                sp.set_linewidth(1.2)

            if cm is not None and norm_obj is not None:
                cax = fig.add_axes([x0 + 0.1445, y0 + 0.068, 0.0052, 0.376])
                cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm_obj, cmap=cm), cax=cax)
                cb.ax.tick_params(colors="#cbd5e1", labelsize=6.8)

        # Tile 11 (r=1, c=5): Live FOV Cells in Stem-Cell-Centered (Δs_∥^stem, Δs_⊥^stem) Sail Space
        x11 = 5 * w_cell
        y11 = 0.00
        ax11 = fig.add_axes([x11 + 0.020, y11 + 0.075, 0.138, 0.368])
        ax11.set_facecolor("#131c31")

        bg_par, bg_perp = _G_SAIL_BG
        ax11.scatter(bg_par, bg_perp, s=2.2, color="#475569", alpha=0.18, rasterized=True)

        ax11.axvline(0.0, color="#38bdf8", ls=":", lw=1.1, alpha=0.85)
        ax11.axhline(0.0, color="#38bdf8", ls=":", lw=1.1, alpha=0.85)
        ax11.axhline(0.85, color="#f43f5e", ls="--", lw=1.3, alpha=0.92, label="ROC Meso+ (+0.85)")
        ax11.axhspan(0.85, 4.5, color="#f43f5e", alpha=0.10)

        if len(spar_vals) > 0:
            ax11.scatter(
                spar_vals,
                sperp_raw_vals,
                c=dp_cell_vals,
                cmap=cmap_dp,
                norm=norm_dp,
                s=10,
                alpha=0.78,
                edgecolors="none",
                zorder=5,
                rasterized=True,
            )

        ax11.set_xlim(-6.5, 6.5)
        ax11.set_ylim(-2.8, 3.8)
        ax11.set_xlabel("Highway Progression Δs_∥^stem (mHz)", fontsize=7.4, fontweight="bold", color="#e2e8f0")
        ax11.set_ylabel("Orthogonal Meso/EMT Δs_⊥^stem (mHz)", fontsize=7.4, fontweight="bold", color="#e2e8f0")
        ax11.set_title("10. Cells on Sail (Colored by ΔP_osc^cell)", fontsize=9.1, fontweight="bold", color="#38bdf8", pad=4)
        ax11.tick_params(colors="#cbd5e1", labelsize=6.8)
        ax11.grid(True, color="#334155", alpha=0.45, lw=0.5)
        ax11.legend(loc="upper right", fontsize=6.2, framealpha=0.85)

        # 1. Save FULL-RESOLUTION 24-bit TrueColor RGB PNG (dpi=150 -> 4320x1575 px, unquantized!)
        #    directly to the fileserver (Z: drive /mnt/researchdrive/.../fov_atlas/)
        remote_png = REMOTE_ATLAS_DIR / out_png.name
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, facecolor=fig.get_facecolor())
        plt.close(fig)
        if REMOTE_ATLAS_DIR.exists():
            remote_png.write_bytes(buf.getvalue())

        # 2. Downscale to 2880x1050 (100 DPI) and quantize to 256-color indexed PNG (~460 KB/FOV)
        #    for the GitHub Pages repository so git stays lightweight (~200 MB total)
        buf.seek(0)
        im_pil = Image.open(buf).convert("RGB")
        im_web = im_pil.resize((2880, 1050), Image.Resampling.LANCZOS)
        im_q = im_web.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
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
                "n_seg_total": n_seg_total,
                "osc_tissue_pct": round(osc_tissue_pct, 1),
                "exp_pctrl_pct": round(exp_pctrl_pct, 1),
                "competent_pct": round(competent_pct, 1),
                "surf_dposc_med": round(med_dp_surf, 1),
                "gated_dposc_med": round(gated_dposc_med, 1),
                "quenched_pct": round(quenched_pct, 1),
                "more_alive_pct": round(more_alive_pct, 1),
                "roc_grey_pct": round(roc_grey_pct, 1),
                "roc_blue_pct": round(roc_blue_pct, 1),
                "roc_pink_pct": round(roc_pink_pct, 1),
                "med_ED_ratio": round(med_ed, 2),
                "iqr_ED_ratio": round(iqr_ed, 2),
                "med_freq_mHz": round(med_f, 1),
                "iqr_freq_mHz": round(iqr_f, 1),
                "med_df_stem_mHz": round(med_df_stem, 2),
                "med_s_parallel_mHz": round(med_s_par, 2),
                "med_s_perp_roc_mHz": round(med_s_perp, 2),
                "pct_meso_roc": round(pct_meso_roc, 1),
                "n_meso_roc": n_meso_roc,
                "dz_hied_pct": round(float(dz_stats["hied_pct"]), 1),
                "dz_zfilt_pct": round(float(dz_stats["zfilt_pct"]), 1),
                "dz_dim_pct": round(float(dz_stats["dim_pct"]), 1),
            }
        )

    return fov_records


def build_interactive_html_viewer(fov_records: list):
    records_sorted = sorted(
        fov_records,
        key=lambda r: (TREATMENTS.index(r["treatment"]), r["day"], r["well"], r["fov_index"]),
    )
    json_payload = json.dumps(records_sorted)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>H9 Stem Cell MinDE Phase-Vocoder — 12-Panel Atlas & Interactive Overlay/Blend Explorer (420 FOVs)</title>
<style>
  :root {{
    --bg: #090d16;
    --panel: #111827;
    --panel-alt: #1e293b;
    --border: #334155;
    --text: #f8fafc;
    --muted: #94a3b8;
    --accent: #38bdf8;
    --control: #3b82f6;
    --ecto: #f59e0b;
    --endo: #10b981;
    --meso: #f43f5e;
    --blend-opacity: 0.62;
    --blend-mode: screen;
    --wipe-right: 14.2%;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    padding: 0;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background: var(--bg);
    color: var(--text);
  }}
  header {{
    background: linear-gradient(90deg, #0f172a 0%, #1e293b 100%);
    border-bottom: 1px solid var(--border);
    padding: 11px 20px;
    display: flex;
    flex-wrap: wrap;
    justify-content: space-between;
    align-items: center;
    gap: 10px;
    position: sticky;
    top: 0;
    z-index: 100;
    box-shadow: 0 4px 16px rgba(0,0,0,0.5);
  }}
  .title-block h1 {{ margin: 0; font-size: 16.5px; font-weight: 700; color: #f8fafc; }}
  .title-block p {{ margin: 2px 0 0 0; font-size: 11.5px; color: var(--muted); }}
  .mode-tabs {{ display: flex; gap: 6px; flex-wrap: wrap; }}
  .mode-btn {{
    background: var(--panel-alt);
    color: var(--text);
    border: 1px solid var(--border);
    padding: 6px 11px;
    border-radius: 6px;
    cursor: pointer;
    font-size: 12px;
    font-weight: 600;
  }}
  .mode-btn:hover {{ border-color: var(--accent); }}
  .mode-btn.active {{ background: #0284c7; border-color: var(--accent); color: #fff; }}
  .controls-bar {{
    background: var(--panel);
    border-bottom: 1px solid var(--border);
    padding: 9px 20px;
    display: flex;
    flex-wrap: wrap;
    gap: 13px;
    align-items: center;
  }}
  .focus-bar {{
    background: #0d1424;
    border-bottom: 1px solid var(--border);
    padding: 8px 20px;
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
    align-items: center;
  }}
  .ctrl-group {{ display: flex; align-items: center; gap: 5px; flex-wrap: wrap; }}
  .ctrl-label {{
    font-size: 10.5px;
    font-weight: 700;
    text-transform: uppercase;
    color: var(--muted);
    margin-right: 3px;
  }}
  .pill-btn {{
    background: var(--panel-alt);
    color: var(--text);
    border: 1px solid var(--border);
    padding: 4px 9px;
    border-radius: 5px;
    cursor: pointer;
    font-size: 11.2px;
    font-weight: 600;
  }}
  .pill-btn:hover {{ border-color: var(--accent); }}
  .pill-btn.active {{ background: var(--accent); color: #0f172a; border-color: var(--accent); font-weight: 700; }}
  .pill-btn.focus-special.active {{ background: #a855f7; color: #fff; border-color: #c084fc; }}
  .pill-btn.preset-btn {{ background: #172554; border-color: #3b82f6; color: #93c5fd; }}
  .pill-btn.preset-btn:hover {{ background: #1e3a8a; color: #fff; }}
  .pill-btn.treat-control.active {{ background: var(--control); color: #fff; border-color: var(--control); }}
  .pill-btn.treat-ecto.active {{ background: var(--ecto); color: #0f172a; border-color: var(--ecto); }}
  .pill-btn.treat-endo.active {{ background: var(--endo); color: #0f172a; border-color: var(--endo); }}
  .pill-btn.treat-meso.active {{ background: var(--meso); color: #fff; border-color: var(--meso); }}
  .kb-hint {{
    margin-left: auto;
    font-size: 11px;
    color: var(--muted);
    background: rgba(255,255,255,0.04);
    padding: 4px 9px;
    border-radius: 5px;
    border: 1px solid var(--border);
  }}
  main {{ padding: 14px 20px 40px 20px; max-width: 1960px; margin: 0 auto; }}
  .metrics-strip {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 8px;
    margin-bottom: 12px;
  }}
  .metric-card {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 7px 11px;
  }}
  .metric-card .m-label {{ font-size: 10px; color: var(--muted); text-transform: uppercase; font-weight: 600; }}
  .metric-card .m-val {{ font-size: 15.5px; font-weight: 700; margin-top: 2px; color: var(--accent); }}
  .metric-card .m-sub {{ font-size: 10px; color: var(--muted); margin-top: 2px; }}
  .fov-card {{
    background: var(--panel);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 11px 13px;
    margin-bottom: 14px;
    box-shadow: 0 3px 12px rgba(0,0,0,0.4);
  }}
  .fov-card-header {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 8px;
    margin-bottom: 9px;
    padding-bottom: 6px;
    border-bottom: 1px solid var(--border);
    font-size: 12.5px;
  }}
  .badge {{ padding: 3px 9px; border-radius: 4px; font-weight: 700; font-size: 11.5px; }}
  .badge-control {{ background: rgba(59,130,246,0.2); color: #60a5fa; border: 1px solid #3b82f6; }}
  .badge-ecto {{ background: rgba(245,158,11,0.2); color: #fbbf24; border: 1px solid #f59e0b; }}
  .badge-endo {{ background: rgba(16,185,129,0.2); color: #34d399; border: 1px solid #10b981; }}
  .badge-meso {{ background: rgba(244,63,94,0.2); color: #fb7185; border: 1px solid #f43f5e; }}
  .fov-img-wrap img {{ width: 100%; height: auto; border-radius: 6px; display: block; background: #0f172a; }}

  /* CSS Viewport Cropping for 2x6 Grid (28.8 x 10.5 in; each cell is 4.8 x 5.25 in) */
  .crop-row {{
    width: 100%;
    aspect-ratio: 28.8 / 5.25;
    overflow: hidden;
    border-radius: 6px;
    background: #0f172a;
    position: relative;
  }}
  .crop-row img {{
    width: 100%;
    height: auto;
    display: block;
  }}
  .crop-row.row-bottom img {{
    transform: translateY(-50%);
  }}
  .single-tile-grid {{
    display: flex;
    flex-wrap: wrap;
    gap: 14px;
    justify-content: center;
    align-items: flex-start;
  }}
  .tile-col-wrap {{
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 4px;
  }}
  .tile-col-tag {{
    font-size: 11.5px;
    font-weight: 700;
    color: #cbd5e1;
    background: #1e293b;
    border: 1px solid #334155;
    padding: 3px 10px;
    border-radius: 5px;
  }}
  .tile-zoom-box {{
    overflow: hidden;
    border-radius: 8px;
    border: 1px solid var(--border);
    background: #0f172a;
    position: relative;
    aspect-ratio: 4.8 / 5.25;
  }}
  .tile-zoom-box img.base-layer {{
    width: 600%;
    max-width: none;
    height: auto;
    display: block;
  }}
  /* Overlay layer inside Blended Box: clipped strictly to the 512x512 colony image rectangle
     [top: 11.2%, right: var(--wipe-right, 14.2%), bottom: 13.6%, left: 3.3%] so titles/colorbars never collide! */
  .tile-zoom-box .overlay-clip {{
    position: absolute;
    inset: 0;
    pointer-events: none;
    clip-path: inset(11.2% var(--wipe-right, 14.2%) 13.6% 3.3%);
    opacity: var(--blend-opacity, 0.62);
    mix-blend-mode: var(--blend-mode, screen);
  }}
  .tile-zoom-box .overlay-clip img {{
    width: 600%;
    max-width: none;
    height: auto;
    display: block;
  }}
  /* Optional Green Tint filter when user toggles Green Tint on Panel 1 */
  .green-tint-filter {{
    filter: sepia(100%) hue-rotate(85deg) saturate(260%) brightness(0.95);
  }}
  select.pair-select {{
    background: var(--panel-alt);
    color: var(--text);
    border: 1px solid var(--border);
    padding: 4px 8px;
    border-radius: 5px;
    font-size: 11.5px;
    font-weight: 600;
  }}
  table.fov-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 11.5px;
    background: var(--panel);
    border-radius: 8px;
    overflow: hidden;
  }}
  table.fov-table th, table.fov-table td {{ padding: 7px 9px; text-align: left; border-bottom: 1px solid var(--border); }}
  table.fov-table th {{ background: var(--panel-alt); color: var(--accent); cursor: pointer; user-select: none; }}
  table.fov-table tr:hover {{ background: rgba(56, 189, 248, 0.08); cursor: pointer; }}
</style>
</head>
<body>

<header>
  <div class="title-block">
    <h1>H9 Stem Cell MinDE Phase-Vocoder — 12-Panel Atlas (With 8A, 8B, 8C) & Interactive Blend/Overlay Explorer</h1>
    <p>4 Lineages × 7 Days × 3 Wells × 5 FOVs (420 FOVs) | Grayscale Still (`*_g.tif`), Dual Still, Stoich `E/D`, FM Barcode `f_mHz`, State `ΔF_stem`, ROC `Meso Δs_⊥`, and All 3 `ΔP_osc` Views (`8A, 8B, 8C`)</p>
  </div>
  <div class="mode-tabs">
    <button class="mode-btn active" onclick="setMode('single')" id="tab-single">1. Single-FOV Inspector</button>
    <button class="mode-btn" onclick="setMode('timecourse')" id="tab-timecourse">2. 7-Day Timecourse (d0.0 → d3.0)</button>
    <button class="mode-btn" onclick="setMode('lineages')" id="tab-lineages">3. 4-Lineage Compare</button>
    <button class="mode-btn" onclick="setMode('montage')" id="tab-montage">4. All 15 FOVs in Condition</button>
    <button class="mode-btn" onclick="setMode('table')" id="tab-table">5. Sortable 420-FOV Table</button>
    <button class="mode-btn" onclick="setMode('figures')" id="tab-figures">6. Key Summary Figures</button>
  </div>
</header>

<div class="controls-bar">
  <div class="ctrl-group">
    <span class="ctrl-label">Lineage:</span>
    <button class="pill-btn treat-control" onclick="setTreat('control')" id="btn-treat-control">Control</button>
    <button class="pill-btn treat-ecto" onclick="setTreat('ecto')" id="btn-treat-ecto">Ectoderm</button>
    <button class="pill-btn treat-endo" onclick="setTreat('endo')" id="btn-treat-endo">Endoderm</button>
    <button class="pill-btn treat-meso active" onclick="setTreat('meso')" id="btn-treat-meso">Mesoderm</button>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">Day:</span>
    <button class="pill-btn" onclick="setDay(0.0)" id="btn-day-0.0">d0.0</button>
    <button class="pill-btn" onclick="setDay(0.5)" id="btn-day-0.5">d0.5</button>
    <button class="pill-btn" onclick="setDay(1.0)" id="btn-day-1.0">d1.0</button>
    <button class="pill-btn" onclick="setDay(1.5)" id="btn-day-1.5">d1.5</button>
    <button class="pill-btn" onclick="setDay(2.0)" id="btn-day-2.0">d2.0</button>
    <button class="pill-btn" onclick="setDay(2.5)" id="btn-day-2.5">d2.5</button>
    <button class="pill-btn active" onclick="setDay(3.0)" id="btn-day-3.0">d3.0</button>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">Well:</span>
    <button class="pill-btn" onclick="setWell(0)" id="btn-well-0">Well 0</button>
    <button class="pill-btn active" onclick="setWell(1)" id="btn-well-1">Well 1</button>
    <button class="pill-btn" onclick="setWell(2)" id="btn-well-2">Well 2</button>
  </div>

  <div class="ctrl-group">
    <span class="ctrl-label">FOV:</span>
    <button class="pill-btn" onclick="setFov(0)" id="btn-fov-0">FOV 1</button>
    <button class="pill-btn" onclick="setFov(1)" id="btn-fov-1">FOV 2</button>
    <button class="pill-btn active" onclick="setFov(2)" id="btn-fov-2">FOV 3</button>
    <button class="pill-btn" onclick="setFov(3)" id="btn-fov-3">FOV 4</button>
    <button class="pill-btn" onclick="setFov(4)" id="btn-fov-4">FOV 5</button>
  </div>

  <div class="kb-hint">
    Keys: <strong>1–4</strong> Lineage | <strong>←/→</strong> Day | <strong>↑/↓</strong> Well | <strong>[ / ]</strong> FOV
  </div>
</div>

<div class="focus-bar">
  <span class="ctrl-label" style="color:#c084fc;">Display Mode:</span>
  <button class="pill-btn focus-special active" onclick="setFocus('all')" id="btn-focus-all">All 12 Panels (2×6 Grid)</button>
  <button class="pill-btn focus-special" onclick="setFocus('row_cohort')" id="btn-focus-row_cohort">Row 1 (Stills + E/D + Barcode + ΔF + Δs_⊥)</button>
  <button class="pill-btn focus-special" onclick="setFocus('row_cell')" id="btn-focus-row_cell">Row 2 (P_ctrl + 8A + 8B + 8C + DeadZones + Sail)</button>
  <button class="pill-btn focus-special" onclick="setFocus('pair')" id="btn-focus-pair">⚔️ Side-by-Side + Blend Overlay (2 Panels)</button>
  <span style="color:var(--border);margin:0 3px;">|</span>
  <span class="ctrl-label" style="color:#38bdf8;">Zoom Single Panel:</span>
  <button class="pill-btn" onclick="setFocus('t0')" id="btn-focus-t0">1. Grayscale Still</button>
  <button class="pill-btn" onclick="setFocus('t1')" id="btn-focus-t1">2. Dual Still</button>
  <button class="pill-btn" onclick="setFocus('t2')" id="btn-focus-t2">3. Stoich E/D</button>
  <button class="pill-btn" onclick="setFocus('t3')" id="btn-focus-t3">4. FM Barcode f</button>
  <button class="pill-btn" onclick="setFocus('t4')" id="btn-focus-t4">5. State ΔF_stem</button>
  <button class="pill-btn" onclick="setFocus('t5')" id="btn-focus-t5">6. ROC Meso Δs_⊥</button>
  <button class="pill-btn" onclick="setFocus('t6')" id="btn-focus-t6">7. Expected P_ctrl</button>
  <button class="pill-btn" onclick="setFocus('t7')" id="btn-focus-t7">8A. Surface ΔP_osc</button>
  <button class="pill-btn" onclick="setFocus('t8')" id="btn-focus-t8">8B. Local Gated ΔP_osc</button>
  <button class="pill-btn" onclick="setFocus('t9')" id="btn-focus-t9">8C. Dual-ROC Surprise</button>
  <button class="pill-btn" onclick="setFocus('t10')" id="btn-focus-t10">9. Dead-Zone Map</button>
  <button class="pill-btn" onclick="setFocus('t11')" id="btn-focus-t11">10. 2D Sail</button>

  <div class="ctrl-group" id="zoom-size-ctrl" style="margin-left:auto;">
    <label style="font-size:11px;color:#cbd5e1;cursor:pointer;margin-right:6px;">
      <input type="checkbox" id="chk-green-tint" onchange="toggleGreenTint(this.checked)"> Green Tint on Ch1
    </label>
    <span class="ctrl-label">Tile Size:</span>
    <input type="range" min="320" max="1050" step="20" value="540" oninput="setZoomPx(this.value)" style="width:95px;">
    <span id="zoom-px-label" style="font-size:11.5px;color:var(--accent);font-weight:700;">540px</span>
  </div>
</div>

<div class="focus-bar" id="pair-selector-bar" style="display:none; background:#131d33; flex-direction:column; align-items:stretch; gap:8px;">
  <div style="display:flex; flex-wrap:wrap; gap:10px; align-items:center;">
    <span class="ctrl-label" style="color:#f472b6;">⚔️ Compare & Blend:</span>
    <label style="font-size:11.5px;">Base Layer (A):
      <select class="pair-select" id="sel-pair-left" onchange="updatePair()">
        <option value="0" selected>1. Grayscale Still (*_g.tif : MinE-GFP)</option>
        <option value="1">2. Dual Still (Red=MinD, Green=MinE)</option>
        <option value="2">3. Stoichiometry E/D (Full Colony)</option>
        <option value="3">4. Raw FM Barcode f_mHz (18.5–49.5)</option>
        <option value="4">5. Decoded State Shift ΔF_stem</option>
        <option value="5">6. ROC Meso/EMT Axis Δs_⊥^stem</option>
        <option value="6">7. Expected Control P_osc^Ctrl(D, E)</option>
        <option value="7">8A. Surface Shift ΔP_osc^surf (Non-Surprise)</option>
        <option value="8">8B. Local Gated ΔP_osc (Blue/Grey/Gold)</option>
        <option value="9">8C. Dual-ROC Surprise (Blue/Flat Grey/Pink)</option>
        <option value="10">9. Full-Colony 5-Class Dead-Zone Map</option>
        <option value="11">10. 2D Sail (Δs_∥ vs. Δs_⊥)</option>
      </select>
    </label>
    <button class="pill-btn" onclick="swapPair()" title="Swap Base A and Overlay B">⇄ Swap A/B</button>
    <label style="font-size:11.5px;">Overlay Layer (B):
      <select class="pair-select" id="sel-pair-right" onchange="updatePair()">
        <option value="0">1. Grayscale Still (*_g.tif : MinE-GFP)</option>
        <option value="1">2. Dual Still (Red=MinD, Green=MinE)</option>
        <option value="2">3. Stoichiometry E/D (Full Colony)</option>
        <option value="3" selected>4. Raw FM Barcode f_mHz (18.5–49.5)</option>
        <option value="4">5. Decoded State Shift ΔF_stem</option>
        <option value="5">6. ROC Meso/EMT Axis Δs_⊥^stem</option>
        <option value="6">7. Expected Control P_osc^Ctrl(D, E)</option>
        <option value="7">8A. Surface Shift ΔP_osc^surf (Non-Surprise)</option>
        <option value="8">8B. Local Gated ΔP_osc (Blue/Grey/Gold)</option>
        <option value="9">8C. Dual-ROC Surprise (Blue/Flat Grey/Pink)</option>
        <option value="10">9. Full-Colony 5-Class Dead-Zone Map</option>
        <option value="11">10. 2D Sail (Δs_∥ vs. Δs_⊥)</option>
      </select>
    </label>

    <span style="color:var(--border);">|</span>
    <span class="ctrl-label" style="color:#38bdf8;">Layout:</span>
    <button class="pill-btn active" onclick="setBlendLayout('3up')" id="btn-blayout-3up">3-Up: [Base A] + [Blended A+B] + [Overlay B]</button>
    <button class="pill-btn" onclick="setBlendLayout('1up')" id="btn-blayout-1up">1-Up: Giant Blended Overlay [A+B]</button>
    <button class="pill-btn" onclick="setBlendLayout('2up')" id="btn-blayout-2up">2-Up: Side-by-Side [Base A] | [Overlay B]</button>
  </div>

  <div style="display:flex; flex-wrap:wrap; gap:12px; align-items:center; padding-top:4px; border-top:1px solid rgba(255,255,255,0.07);">
    <div class="ctrl-group">
      <span class="ctrl-label" style="color:#fde047;">Overlay Opacity (Transparency):</span>
      <input type="range" id="rng-opacity" min="0" max="100" step="1" value="62" oninput="setBlendOpacity(this.value)" style="width:155px; accent-color:#38bdf8;">
      <span id="lbl-opacity" style="font-size:12px;font-weight:700;color:#fde047;min-width:38px;">62%</span>
      <button class="pill-btn" onclick="setBlendOpacity(0)">0% (Base A)</button>
      <button class="pill-btn" onclick="setBlendOpacity(40)">40%</button>
      <button class="pill-btn" onclick="setBlendOpacity(65)">65%</button>
      <button class="pill-btn" onclick="setBlendOpacity(85)">85%</button>
      <button class="pill-btn" onclick="setBlendOpacity(100)">100% (Overlay B)</button>
    </div>

    <span style="color:var(--border);">|</span>
    <div class="ctrl-group">
      <span class="ctrl-label" style="color:#a7f3d0;">Blend Formula:</span>
      <select class="pair-select" id="sel-blend-mode" onchange="setBlendMode(this.value)">
        <option value="screen" selected>Screen / Glow (Transparent Dark BG — Best for Grayscale + Color!)</option>
        <option value="normal">Normal Alpha Crossfade (Linear Transparency)</option>
        <option value="hard-light">Hard-Light (Preserves Sharp Grayscale Cell Texture)</option>
        <option value="overlay">Overlay (High-Contrast Texture Modulation)</option>
        <option value="lighten">Lighten (Max of Base A & Overlay B)</option>
      </select>
    </div>

    <span style="color:var(--border);">|</span>
    <div class="ctrl-group">
      <span class="ctrl-label" style="color:#c084fc;">Curtain Wipe:</span>
      <input type="range" id="rng-wipe" min="0" max="100" step="1" value="100" oninput="setCurtainWipe(this.value)" style="width:105px; accent-color:#c084fc;" title="Swipe Overlay B horizontally across Base A">
      <span id="lbl-wipe" style="font-size:11.5px;color:#c084fc;font-weight:700;">100%</span>
    </div>
  </div>

  <div style="display:flex; flex-wrap:wrap; gap:6px; align-items:center; padding-top:4px; border-top:1px solid rgba(255,255,255,0.07);">
    <span class="ctrl-label" style="color:#93c5fd;">1-Click Overlay Presets:</span>
    <button class="pill-btn preset-btn" onclick="setPairPreset(0, 3)">1. Grayscale + 4. FM Barcode</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(0, 4)">1. Grayscale + 5. Decoded ΔF_stem</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(0, 5)">1. Grayscale + 6. ROC Meso Δs_⊥</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(0, 9)">1. Grayscale + 8C. Dual-ROC Surprise</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(0, 8)">1. Grayscale + 8B. Local Gated ΔP_osc</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(0, 7)">1. Grayscale + 8A. Surface Shift</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(3, 4)">4. FM Barcode vs. 5. Decoded ΔF_stem</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(5, 9)">6. ROC Meso Δs_⊥ vs. 8C. Dual-ROC Surprise</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(7, 9)">8A. Surface Shift vs. 8C. Dual-ROC Surprise</button>
    <button class="pill-btn preset-btn" onclick="setPairPreset(8, 9)">8B. Local Gated vs. 8C. Dual-ROC Surprise</button>
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
  "1. Grayscale Still (*_g.tif)",
  "2. Dual Still (Red+Green)",
  "3. Stoichiometry E/D",
  "4. Raw FM Barcode f_mHz",
  "5. Decoded State ΔF_stem",
  "6. ROC Meso/EMT Δs_⊥^stem",
  "7. Expected Control P_ctrl(D,E)",
  "8A. Surface Shift ΔP_osc^surf",
  "8B. Local Gated ΔP_osc (±12%)",
  "8C. Dual-ROC Surprise (Blue/Grey/Pink)",
  "9. 5-Class Dead-Zone Map",
  "10. 2D Sail (Δs_∥ vs. Δs_⊥)"
];

let state = {{
  mode: 'single',
  treatment: 'meso',
  day: 3.0,
  well: 1,
  fov: 2,
  focus: 'all',
  zoomPx: 540,
  pairLeft: 0,
  pairRight: 3,
  blendLayout: '3up',
  blendOpacity: 62,
  blendMode: 'screen',
  curtainWipe: 100,
  greenTint: false,
  sortCol: 'pct_meso_roc',
  sortAsc: false
}};

function findRecord(treat, day, well, fov) {{
  return FOV_DATA.find(r => r.treatment === treat && Math.abs(r.day - day) < 0.01 && r.well === well && r.fov_index === fov);
}}

function tileTransformCoords(tileIndex) {{
  const row = Math.floor(tileIndex / 6);
  const col = tileIndex % 6;
  const tx = -(col * 100.0 / 6.0);
  const ty = -row * 50.0;
  return {{ tx, ty }};
}}

function renderTileCrop(imgSrc, tileIndex, widthPx, labelText) {{
  const {{ tx, ty }} = tileTransformCoords(tileIndex);
  const tintCls = (tileIndex === 0 && state.greenTint) ? 'green-tint-filter' : '';
  const tagHtml = labelText ? `<div class="tile-col-tag">${{labelText}}</div>` : '';
  return `
    <div class="tile-col-wrap">
      ${{tagHtml}}
      <div class="tile-zoom-box" style="width:${{widthPx}}px;">
        <img class="base-layer ${{tintCls}}" src="${{imgSrc}}" loading="lazy" style="transform: translate(${{tx.toFixed(4)}}%, ${{ty.toFixed(1)}}%);">
      </div>
    </div>
  `;
}}

function renderBlendedTileCrop(imgSrc, leftIdx, rightIdx, widthPx) {{
  const cA = tileTransformCoords(leftIdx);
  const cB = tileTransformCoords(rightIdx);
  const tintA = (leftIdx === 0 && state.greenTint) ? 'green-tint-filter' : '';
  const tintB = (rightIdx === 0 && state.greenTint) ? 'green-tint-filter' : '';
  return `
    <div class="tile-col-wrap">
      <div class="tile-col-tag" style="border-color:#38bdf8;color:#fde047;">
        Blended Overlay: [${{PANEL_NAMES[leftIdx]}}] + [${{PANEL_NAMES[rightIdx]}}]
      </div>
      <div class="tile-zoom-box" style="width:${{widthPx}}px; border-color:#38bdf8; box-shadow: 0 0 12px rgba(56,189,248,0.25);">
        <img class="base-layer ${{tintA}}" src="${{imgSrc}}" loading="lazy" style="transform: translate(${{cA.tx.toFixed(4)}}%, ${{cA.ty.toFixed(1)}}%);">
        <div class="overlay-clip">
          <img class="${{tintB}}" src="${{imgSrc}}" loading="lazy" style="transform: translate(${{cB.tx.toFixed(4)}}%, ${{cB.ty.toFixed(1)}}%);">
        </div>
      </div>
    </div>
  `;
}}

function renderVisualBody(r) {{
  if (state.focus === 'all') {{
    return `<div class="fov-img-wrap"><img src="${{r.img_rel_path}}" loading="lazy"></div>`;
  }} else if (state.focus === 'row_cohort') {{
    return `<div class="crop-row row-top"><img src="${{r.img_rel_path}}" loading="lazy"></div>`;
  }} else if (state.focus === 'row_cell') {{
    return `<div class="crop-row row-bottom"><img src="${{r.img_rel_path}}" loading="lazy"></div>`;
  }} else if (state.focus === 'pair') {{
    if (state.blendLayout === '1up') {{
      return `
        <div class="single-tile-grid">
          ${{renderBlendedTileCrop(r.img_rel_path, state.pairLeft, state.pairRight, Math.max(state.zoomPx, 680))}}
        </div>
      `;
    }} else if (state.blendLayout === '2up') {{
      return `
        <div class="single-tile-grid">
          ${{renderTileCrop(r.img_rel_path, state.pairLeft, state.zoomPx, 'Base A: ' + PANEL_NAMES[state.pairLeft])}}
          ${{renderTileCrop(r.img_rel_path, state.pairRight, state.zoomPx, 'Overlay B: ' + PANEL_NAMES[state.pairRight])}}
        </div>
      `;
    }} else {{
      return `
        <div class="single-tile-grid">
          ${{renderTileCrop(r.img_rel_path, state.pairLeft, state.zoomPx, 'Base A: ' + PANEL_NAMES[state.pairLeft])}}
          ${{renderBlendedTileCrop(r.img_rel_path, state.pairLeft, state.pairRight, state.zoomPx)}}
          ${{renderTileCrop(r.img_rel_path, state.pairRight, state.zoomPx, 'Overlay B: ' + PANEL_NAMES[state.pairRight])}}
        </div>
      `;
    }}
  }} else if (state.focus.startsWith('t')) {{
    const tIdx = parseInt(state.focus.substring(1), 10);
    return `<div class="single-tile-grid">${{renderTileCrop(r.img_rel_path, tIdx, Math.max(state.zoomPx, 680), PANEL_NAMES[tIdx])}}</div>`;
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
          <strong style="margin-left:8px;">Day ${{r.day.toFixed(1)}}</strong> &bull;
          Well ${{r.well}} &bull; FOV ${{r.fov_index + 1}} (<code>${{r.filename}}</code>)
        </div>
        <div style="color: var(--muted); font-size: 11.8px;">
          Osc Cells: <strong>${{r.n_cells}}</strong> (OscPx: <strong>${{r.osc_tissue_pct.toFixed(0)}}%</strong> vs Exp P_ctrl: <strong>${{r.exp_pctrl_pct.toFixed(0)}}%</strong>) |
          Barcode f: <strong style="color:#fbbf24">${{r.med_freq_mHz.toFixed(1)}} mHz</strong> |
          ΔF_stem: <strong style="color:#60a5fa">${{r.med_df_stem_mHz >= 0 ? '+' : ''}}${{r.med_df_stem_mHz.toFixed(1)}} mHz</strong> |
          ROC Meso+: <strong style="color:#f43f5e">${{r.pct_meso_roc.toFixed(1)}}%</strong> |
          8A Surf: <strong style="color:#7dd3fc">${{r.surf_dposc_med >= 0 ? '+' : ''}}${{r.surf_dposc_med.toFixed(0)}}%</strong> |
          8C Surprise: <strong style="color:#f472b6">Grey ${{r.roc_grey_pct.toFixed(0)}}% / Blue ${{r.roc_blue_pct.toFixed(0)}}% / Pink ${{r.roc_pink_pct.toFixed(0)}}%</strong>
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

  ['all', 'row_cohort', 'row_cell', 'pair', 't0', 't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11'].forEach(fk => {{
    const el = document.getElementById('btn-focus-' + fk);
    if (el) el.classList.toggle('active', state.focus === fk);
  }});

  ['3up', '1up', '2up'].forEach(lk => {{
    const el = document.getElementById('btn-blayout-' + lk);
    if (el) el.classList.toggle('active', state.blendLayout === lk);
  }});

  document.getElementById('pair-selector-bar').style.display = (state.focus === 'pair') ? 'flex' : 'none';

  const container = document.getElementById('view-container');

  if (state.mode === 'single') {{
    const r = findRecord(state.treatment, state.day, state.well, state.fov);
    let metricsHtml = '';
    if (r) {{
      const fk = fovKey(r);
      const tagObj = fovTags[fk] || {{ starred: false, tag: '', note: '' }};
      const starCount = Object.values(fovTags).filter(v => v && v.starred).length;
      const edgeDiff = (r.df_edge_med !== undefined && r.df_int_med !== undefined) ? (r.df_edge_med - r.df_int_med) : 0.0;
      metricsHtml = `
        <div class="metrics-strip">
          <div class="metric-card"><div class="m-label">Oscillating Cells (≥10 px)</div><div class="m-val" style="color:#34d399">${{r.n_cells}} (${{r.osc_tissue_pct.toFixed(1)}}%)</div><div class="m-sub">Real Masks ≥100px: ${{r.n_seg_100 || r.n_cells}} | ${{r.n_pixels.toLocaleString()}} px</div></div>
          <div class="metric-card"><div class="m-label">Stoichiometry E/D</div><div class="m-val" style="color:#38bdf8">${{r.med_ED_ratio.toFixed(2)}}</div><div class="m-sub">Within-Colony IQR = ${{r.iqr_ED_ratio.toFixed(2)}}</div></div>
          <div class="metric-card"><div class="m-label">Raw FM Barcode (f_mHz)</div><div class="m-val" style="color:#fbbf24">${{r.med_freq_mHz.toFixed(1)}} mHz</div><div class="m-sub">Carrier IQR = ${{r.iqr_freq_mHz.toFixed(1)}} mHz</div></div>
          <div class="metric-card"><div class="m-label">Decoded State ΔF_stem</div><div class="m-val" style="color:#60a5fa">${{r.med_df_stem_mHz >= 0 ? '+' : ''}}${{r.med_df_stem_mHz.toFixed(2)}} mHz</div><div class="m-sub">Stoich-Decoded vs. Control</div></div>
          <div class="metric-card"><div class="m-label">Colony-Edge Zonation (d_edge)</div><div class="m-val" style="color:#fde047">${{edgeDiff >= 0 ? '+' : ''}}${{edgeDiff.toFixed(2)}} mHz</div><div class="m-sub">Edge(&lt;40px): ${{r.df_edge_med !== undefined ? r.df_edge_med.toFixed(1) : '—'}} | Int(≥40px): ${{r.df_int_med !== undefined ? r.df_int_med.toFixed(1) : '—'}}</div></div>
          <div class="metric-card"><div class="m-label">ROC Meso/EMT (Δs_⊥&gt;+0.85)</div><div class="m-val" style="color:#f43f5e">${{r.pct_meso_roc.toFixed(1)}}% (n=${{r.n_meso_roc}})</div><div class="m-sub">Med Δs_⊥^stem = ${{r.med_s_perp_roc_mHz >= 0 ? '+' : ''}}${{r.med_s_perp_roc_mHz.toFixed(2)}} mHz</div></div>
          <div class="metric-card"><div class="m-label">8A/8B. Surface & Local ΔP_osc</div><div class="m-val" style="color:#7dd3fc">${{r.surf_dposc_med >= 0 ? '+' : ''}}${{r.surf_dposc_med.toFixed(1)}}% / ${{r.gated_dposc_med >= 0 ? '+' : ''}}${{r.gated_dposc_med.toFixed(1)}}%</div><div class="m-sub">Exp P_ctrl: ${{r.exp_pctrl_pct.toFixed(1)}}% | Quench: ${{r.quenched_pct.toFixed(0)}}%</div></div>
          <div class="metric-card"><div class="m-label">8C. Dual-ROC Surprise</div><div class="m-val" style="color:#f472b6">Pink ${{r.roc_pink_pct.toFixed(0)}}% | Blue ${{r.roc_blue_pct.toFixed(0)}}%</div><div class="m-sub">Flat Grey (Expected): ${{r.roc_grey_pct.toFixed(0)}}%</div></div>
        </div>
        <div class="controls-bar" style="margin-bottom:12px;border:1px solid var(--border);border-radius:8px;background:#0d1526;padding:8px 14px;">
          <span class="ctrl-label" style="color:#fde047;">FOV Bookmark &amp; Notes:</span>
          <button class="pill-btn ${{tagObj.starred ? 'active' : ''}}" style="${{tagObj.starred ? 'background:#eab308;border-color:#fde047;color:#0f172a;' : ''}}" onclick="toggleStarFov('${{fk}}')">${{tagObj.starred ? '★ Starred' : '☆ Star FOV'}}</button>
          <select style="background:#1e293b;color:#f8fafc;border:1px solid #475569;border-radius:5px;padding:4px 8px;font-size:11.5px;" onchange="setFovTagCategory('${{fk}}', this.value)">
            <option value="" ${{tagObj.tag === '' ? 'selected' : ''}}>-- Tag Category --</option>
            <option value="Edge-Interior Gradient" ${{tagObj.tag === 'Edge-Interior Gradient' ? 'selected' : ''}}>Edge-Interior Gradient</option>
            <option value="Quenched Rosette/Core" ${{tagObj.tag === 'Quenched Rosette/Core' ? 'selected' : ''}}>Quenched Rosette/Core</option>
            <option value="Low-Threshold Active (Pink)" ${{tagObj.tag === 'Low-Threshold Active (Pink)' ? 'selected' : ''}}>Low-Threshold Active (Pink)</option>
            <option value="Non-Clonal Wave Domain" ${{tagObj.tag === 'Non-Clonal Wave Domain' ? 'selected' : ''}}>Non-Clonal Wave Domain</option>
            <option value="Showcase Candidate" ${{tagObj.tag === 'Showcase Candidate' ? 'selected' : ''}}>Showcase Candidate</option>
            <option value="Artifact / Check" ${{tagObj.tag === 'Artifact / Check' ? 'selected' : ''}}>Artifact / Check</option>
          </select>
          <input type="text" placeholder="Add research note for this FOV (auto-saved in browser)..." value="${{(tagObj.note || '').replace(/"/g, '&quot;')}}" style="flex:1;min-width:220px;background:#1e293b;color:#f8fafc;border:1px solid #475569;border-radius:5px;padding:4px 10px;font-size:11.8px;" oninput="setFovNote('${{fk}}', this.value)">
          <button class="pill-btn preset-btn" onclick="exportFovTagsCsv()">⬇ Export Notes CSV (${{starCount}} ★)</button>
        </div>
      `;
    }}
    container.innerHTML = metricsHtml + renderCard(r);
  }} else if (state.mode === 'timecourse') {{
    let html = `<h3 style="margin-top:0;">7-Day Timecourse — ${{state.treatment.toUpperCase()}} (Well ${{state.well}}, FOV ${{state.fov + 1}})</h3>`;
    DAYS.forEach(d => {{
      html += renderCard(findRecord(state.treatment, d, state.well, state.fov));
    }});
    container.innerHTML = html;
  }} else if (state.mode === 'lineages') {{
    let html = `<h3 style="margin-top:0;">4-Lineage Comparison — Day ${{state.day.toFixed(1)}} (Well ${{state.well}}, FOV ${{state.fov + 1}})</h3>`;
    TREATMENTS.forEach(t => {{
      html += renderCard(findRecord(t, state.day, state.well, state.fov));
    }});
    container.innerHTML = html;
  }} else if (state.mode === 'montage') {{
    let html = `<h3 style="margin-top:0;">All 15 Replicate FOVs (3 Wells × 5 FOVs) — ${{state.treatment.toUpperCase()}} Day ${{state.day.toFixed(1)}}</h3>`;
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
        <td>${{r.osc_tissue_pct.toFixed(1)}}%</td>
        <td style="color:#fde047;">${{r.exp_pctrl_pct.toFixed(1)}}%</td>
        <td style="color:#7dd3fc;font-weight:700;">${{r.surf_dposc_med >= 0 ? '+' : ''}}${{r.surf_dposc_med.toFixed(1)}}%</td>
        <td style="color:#fb923c;">${{r.quenched_pct.toFixed(0)}}% / ${{r.more_alive_pct.toFixed(0)}}%</td>
        <td style="color:#cbd5e1;">${{r.roc_grey_pct.toFixed(0)}}%</td>
        <td style="color:#38bdf8;font-weight:700;">${{r.roc_blue_pct.toFixed(0)}}%</td>
        <td style="color:#f472b6;font-weight:700;">${{r.roc_pink_pct.toFixed(0)}}%</td>
        <td>${{r.med_ED_ratio.toFixed(2)}}</td>
        <td style="color:#fbbf24;">${{r.med_freq_mHz.toFixed(1)}}</td>
        <td style="color:#60a5fa;font-weight:700;">${{r.med_df_stem_mHz >= 0 ? '+' : ''}}${{r.med_df_stem_mHz.toFixed(2)}}</td>
        <td style="color:#f43f5e;font-weight:700;">${{r.pct_meso_roc.toFixed(1)}}% (${{r.n_meso_roc}})</td>
      </tr>
    `).join('');
    container.innerHTML = `
      <p style="color:var(--muted);margin-top:0;">Click any column header to sort all 420 FOVs, or click any row to jump directly to that FOV.</p>
      <table class="fov-table">
        <thead>
          <tr>
            <th onclick="sortTable('treatment')">Lineage</th>
            <th onclick="sortTable('day')">Day</th>
            <th onclick="sortTable('well')">Well</th>
            <th onclick="sortTable('fov_index')">FOV</th>
            <th onclick="sortTable('n_cells')">Osc Cells</th>
            <th onclick="sortTable('osc_tissue_pct')">OscPx %</th>
            <th onclick="sortTable('exp_pctrl_pct')">Exp P_ctrl</th>
            <th onclick="sortTable('surf_dposc_med')">8A Surf ΔP</th>
            <th onclick="sortTable('quenched_pct')">8B Quench/+Alive</th>
            <th onclick="sortTable('roc_grey_pct')">8C Grey %</th>
            <th onclick="sortTable('roc_blue_pct')">8C Blue %</th>
            <th onclick="sortTable('roc_pink_pct')">8C Pink %</th>
            <th onclick="sortTable('med_ED_ratio')">Med E/D</th>
            <th onclick="sortTable('med_freq_mHz')">Barcode f</th>
            <th onclick="sortTable('med_df_stem_mHz')">State ΔF_stem</th>
            <th onclick="sortTable('pct_meso_roc')">ROC Meso+ % (n)</th>
          </tr>
        </thead>
        <tbody>${{rows}}</tbody>
      </table>
    `;
  }} else if (state.mode === 'figures') {{
    const figs = [
      {{
        title: "Figure 50 — Out-of-Sample Reference Robustness (LOWO, Day-Matched, Edge-Matched P_ctrl), Knot-Free Domain LLR, & Competent-Mask Audit",
        sub: "Demonstrates that all headline Axis 1 and Axis 2 conclusions survive Leave-One-Well-Out (LOWO), Day-Matched, and Edge-Distance-Stratified P_ctrl(D, E | d_edge) references across n = 3 replicate wells; resolves the Mesoderm competent-mask question via Relaxed-Z (4.0 ≤ D_z < 5.0 → 97–98% oscillation rate); and shows well-to-well replicate drift is <1.7% of variance.",
        src: "figures/Fig50_Reference_Robustness_EdgeMatched_DayMatched_LOMO_Audit.png"
      }},
      {{
        title: "Figure 49 — Colony-Edge Physiological Zonation (24/24 Wells) & Edge-Conditioned Non-Clonal Neighbor Wave Coherence",
        sub: "Across all 448,054 oscillating cells (n = 3 wells/condition), pluripotent Control colonies maintain a steady +2.3 to +3.8 mHz edge-to-interior ΔF_stem gradient within every E/D bin, while Ectoderm (−10.2 mHz interior vs. −1.9 mHz edge) and Endoderm (−5.3 mHz interior vs. −1.4 mHz edge) nucleate their frequency slowdowns from the inside out. Touching non-clonal neighbor coherence survives a strict joint (E/D, d_edge)-matched FOV permutation null while exhibiting uniform cross-boundary phase differences (median |Δφ| = 90.2°, KS p = 0.52) and a sharp ±1 px step at cell borders.",
        src: "figures/Fig49_Colony_Edge_Zonation_and_Edge_Conditioned_Neighbor_Coherence.png"
      }},
      {{
        title: "Figure 48 — Master Representative Contact Sheet: All Spatial Views & Candidate ΔP_osc Formulations Across All 4 Lineages (Day 3.0)",
        sub: "10-column side-by-side comparison across Control, Ectoderm, Endoderm, and Mesoderm: 1. Green Still, 2. Dual Still, 3. Stoichiometry E/D, 4. Raw FM Barcode f_mHz, 5. Decoded State ΔF_stem, 6. ROC Meso/EMT Δs_⊥^stem, 7. Expected Control P_ctrl(D, E), 8A. Surface Shift ΔP_osc^surf, 8B. Local Gated ΔP_osc (Blue/Grey/Gold), and 8C. Dual-ROC Calibrated Surprise (Blue/Flat-Grey/Pink).",
        src: "figures/Fig48_Representative_Contact_Sheet_All_Spatial_Views_4Lineages.png"
      }},
      {{
        title: "Figure 47 — Dual-ROC Calibrated Oscillation Surprise Map: Control Baseline ('1-on-1' Grey) vs. Quenched ('0-on-1' Blue) & Low-Threshold Active ('1-on-0' Pink)",
        sub: "ROC calibration on both tails keeps Control colony tissue Flat Grey (|S_roc| < 28.5) while separating Ectoderm/Endoderm quenched domains in Blue (P_dom − P_ctrl ≤ −22.0%, Youden τ_-* = −22.0%, AUROC = 0.737) and low-threshold active Mesoderm cells in Pink (P_ctrl ≤ 32.5%, Youden τ_+* = +61.5%, AUROC = 0.831).",
        src: "figures/Fig47_Dual_ROC_Surprise_Flat_Control_Blue_Quench_Pink_Meso.png"
      }},
      {{
        title: "Figure 46 — Marrying Bulk Trajectories with the Fit-Free Dual Sail: (D, E) Phase-Portrait Bins in (ΔF_stem, ΔP_osc^stem)",
        sub: "Evaluating (ΔF_stem(D, E), ΔP_osc^stem(D, E)) on the 2D (D, E) stoichiometry bins sharpens Control 15× × 6× into a tight (0, 0) bullseye and reveals how Low-D (Cyan) stoichiometry bins arch into +ΔP_osc > 0 in Day 0.0 and Mesoderm EMT.",
        src: "figures/Fig46_StoichiometryBin_FitFree_Dual_Sail_dF_vs_dPosc.png"
      }},
      {{
        title: "Figure 45 — Fit-Free Dual Encoding & Excitability Space: Stoichiometry-Matched ΔF_stem vs. ΔP_osc^stem in 2D, 3D (× Day), and 4×7 Landscapes",
        sub: "Directly pairing ΔF_stem(D, E) = f − F_Control(D, E) with its dual ΔP_osc^stem(D, E) = P_osc − P_Control(D, E) separates all 4 developmental trajectories, freshly singularized Day 0.0 cells, and Mesoderm EMT without any curve-fitting.",
        src: "figures/Fig45_FitFree_Dual_Space_dFstem_vs_dPosc_2D_3D_and_Sail.png"
      }},
      {{
        title: "Figure 44 — Control-Competent Gated ΔP_osc^stem & Coupling of '+ΔP_osc (More Alive Than Control)' to ROC Meso/EMT (Δs_⊥^stem)",
        sub: "Masking non-competent '0-on-0' territory (where Control P_osc < 20%) isolates true intrinsic quenching in Ectoderm/Endoderm (Deep Blue) and reveals that '+ΔP_osc > +15% (Orange/Gold)' co-localizes with ROC Meso+ (Δs_⊥^stem ≥ +0.85 mHz).",
        src: "figures/Fig44_ControlCompetent_Gated_dPosc_and_Meso_Coupling.png"
      }},
      {{
        title: "Figure 43 — Stoichiometry-Matched Oscillation Probability Surface P_osc^Control(D, E) & Lineage ΔP_osc^stem(D, E)",
        sub: "In steady-state Control stem cells, oscillation competence forms a sharp wedge (D ≥ 220, E/D ∈ [2.5, 5.2]). Subtracting P_osc^Control(D, E) defines the stoichiometry-matched oscillation probability shift ΔP_osc^stem.",
        src: "figures/Fig43_Stoichiometry_Matched_Oscillation_Probability.png"
      }},
      {{
        title: "Figure 42 — Transcending the Clonal FM Barcode & Diagnosing Colony Dead Zones in (D, E) Expression Space",
        sub: "Across 358,000+ touching neighbor cell pairs, 40–50% have discordant carrier stoichiometry (|Δ(E/D)| > 1.0) and wildly different raw frequencies f_mHz, yet share identical decoded physiological states ΔF_stem (removing 52–63% of variance).",
        src: "figures/Fig42_Barcode_Transcendence_and_DeadZone_Diagnostics.png"
      }},
      {{
        title: "Figure 41 — 4-Channel FOV Decoding: Input Green Still (*_g.tif), Raw FM Barcode (f_mHz), Decoded State (ΔF_stem), and ROC-Tuned Meso/EMT Axis",
        sub: "ROC optimization across 448,054 cells identifies τ* = +0.85 mHz (Youden's J = 0.442, AUROC = 0.779, 4.7× odds enrichment) as the optimal grey-to-pink transition isolating Mesoderm/EMT.",
        src: "figures/Fig41_4Channel_GreenStill_FMBarcode_StateShift_ROCMeso.png"
      }},
      {{
        title: "Figure 40 — Dataset-Wide (419 FOVs) Spatial Wave Coherence, Patchy Domain Formation & Tissue Remodeling",
        sub: "Within-colony spatial autocorrelation Δ|ΔF_stem|(r) confirms patchy multi-cell physiological micro-domains (~45–75 µm) rather than a deterministic radial inside-out ring.",
        src: "figures/Fig40_Spatial_Wave_Coherence_419FOVs.png"
      }},
      {{
        title: "Figure 34 — Global 2D Stem-Cell-Centered Sail Manifold (Δs_∥^stem vs. Δs_⊥^stem) Across All 4 Lineages & 7 Days",
        sub: "8×7 hybrid view pairing log-density occupancy on the stem-cell-centered 2D Sail with zoomed lineage trajectory insets.",
        src: "figures/Fig34_Global_2D_Sail_Trajectories.png"
      }}
    ];
    container.innerHTML = figs.map(fg => `
      <div class="fov-card">
        <div class="fov-card-header">
          <div><strong style="color:var(--accent);font-size:14px;">${{fg.title}}</strong></div>
        </div>
        <p style="color:var(--muted);font-size:12.5px;margin:0 0 10px 0;">${{fg.sub}}</p>
        <div class="fov-img-wrap"><img src="${{fg.src}}" loading="lazy"></div>
      </div>
    `).join('');
  }}
}}

const STORAGE_KEY_TAGS = 'h9_vocoder_fov_tags_v1';
let fovTags = {{}};
try {{
  const rawTags = localStorage.getItem(STORAGE_KEY_TAGS);
  if (rawTags) fovTags = JSON.parse(rawTags);
}} catch (e) {{ fovTags = {{}}; }}

function saveFovTags() {{
  try {{ localStorage.setItem(STORAGE_KEY_TAGS, JSON.stringify(fovTags)); }} catch (e) {{}}
}}
function fovKey(r) {{
  return `${{r.treatment}}_d${{r.day.toFixed(1)}}_w${{r.well}}_f${{r.fov_index}}`;
}}
function toggleStarFov(fk) {{
  if (!fovTags[fk]) fovTags[fk] = {{ starred: false, tag: '', note: '' }};
  fovTags[fk].starred = !fovTags[fk].starred;
  saveFovTags();
  updateUI();
}}
function setFovTagCategory(fk, val) {{
  if (!fovTags[fk]) fovTags[fk] = {{ starred: true, tag: '', note: '' }};
  fovTags[fk].tag = val;
  if (val) fovTags[fk].starred = true;
  saveFovTags();
  updateUI();
}}
function setFovNote(fk, val) {{
  if (!fovTags[fk]) fovTags[fk] = {{ starred: true, tag: '', note: '' }};
  fovTags[fk].note = val;
  saveFovTags();
}}
function exportFovTagsCsv() {{
  const rows = [['key','treatment','day','well','fov_index','filename','starred','tag','note','n_cells','n_seg_100','med_ED_ratio','med_df_stem_mHz','df_edge_med','df_int_med','roc_grey_pct','roc_blue_pct','roc_pink_pct']];
  FOV_DATA.forEach(r => {{
    const fk = fovKey(r);
    const t = fovTags[fk];
    if (t && (t.starred || t.tag || (t.note && t.note.trim()))) {{
      rows.push([
        fk, r.treatment, r.day.toFixed(1), r.well, r.fov_index + 1, r.filename,
        t.starred ? '1' : '0', `"${{(t.tag || '').replace(/"/g, '""')}}"`, `"${{(t.note || '').replace(/"/g, '""')}}"`,
        r.n_cells, r.n_seg_100 || r.n_cells, r.med_ED_ratio.toFixed(2), r.med_df_stem_mHz.toFixed(2),
        (r.df_edge_med !== undefined ? r.df_edge_med.toFixed(2) : ''), (r.df_int_med !== undefined ? r.df_int_med.toFixed(2) : ''),
        r.roc_grey_pct.toFixed(1), r.roc_blue_pct.toFixed(1), r.roc_pink_pct.toFixed(1)
      ]);
    }}
  }});
  if (rows.length === 1) {{
    alert('No starred or annotated FOVs yet! Star or add a note to any FOV first.');
    return;
  }}
  const csv = rows.map(r => r.join(',')).join('\n');
  const blob = new Blob([csv], {{ type: 'text/csv;charset=utf-8;' }});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'H9_Vocoder_Tagged_FOV_Notes.csv';
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}}

function setMode(m) {{ state.mode = m; updateUI(); }}
function setTreat(t) {{ state.treatment = t; updateUI(); }}
function setDay(d) {{ state.day = d; updateUI(); }}
function setWell(w) {{ state.well = w; updateUI(); }}
function setFov(f) {{ state.fov = f; updateUI(); }}
function setFocus(fk) {{ state.focus = fk; updateUI(); }}
function setBlendLayout(lk) {{ state.blendLayout = lk; updateUI(); }}
function toggleGreenTint(val) {{ state.greenTint = val; updateUI(); }}
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
  // Image occupies [left=3.3%, right=14.2%] (width=82.5% of tile).
  const rightInset = 14.2 + (100 - state.curtainWipe) * 0.825;
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
function setPairPreset(lIdx, rIdx) {{
  state.focus = 'pair';
  state.pairLeft = lIdx;
  state.pairRight = rIdx;
  document.getElementById('sel-pair-left').value = String(lIdx);
  document.getElementById('sel-pair-right').value = String(rIdx);
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
    print(f"Saved upgraded 12-panel interactive HTML viewer: {local_index} & {local_html}", flush=True)


def main():
    t0 = time.perf_counter()
    ATLAS_DIR.mkdir(parents=True, exist_ok=True)
    REMOTE_ATLAS_DIR.mkdir(parents=True, exist_ok=True)

    if "--html-only" in sys.argv:
        local_json = LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json"
        recs = json.loads(local_json.read_text(encoding="utf-8"))
        build_interactive_html_viewer(recs)
        shutil.copy2(LOCAL_WORK_DIR / "H9_FOV_Physiology_Atlas_Viewer.html", REMOTE_OUT_DIR / "H9_FOV_Physiology_Atlas_Viewer.html")
        shutil.copy2(LOCAL_WORK_DIR / "index.html", REMOTE_OUT_DIR / "index.html")
        return

    bundle_path = build_precomputed_state_bundle()
    ckpt_files = [str(p) for p in sorted(WELL_CKPT_DIR.glob("*.parquet"))]

    if "--test-one" in sys.argv:
        _init_worker(bundle_path)
        test_ckpt = [cp for cp in ckpt_files if "meso_day_3.0_well_1" in cp][0]
        recs = process_single_well(test_ckpt)
        print(f"Test one well complete ({len(recs)} FOVs) in {time.perf_counter() - t0:.1f}s: {recs[2]['img_rel_path']}", flush=True)
        return

    print(f"Rendering upgraded 12-panel (2x6) atlas for all {len(ckpt_files)} wells (420 FOVs total)...", flush=True)

    all_fov_records = []
    done_wells = 0
    with ProcessPoolExecutor(max_workers=8, initializer=_init_worker, initargs=(bundle_path,)) as ex:
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

    print("Syncing HTML viewer, metadata, and FULL-RES summary figures to /mnt/researchdrive/General/Scott/stemcellcodes/vocoder_2026_10/...", flush=True)
    shutil.copy2(LOCAL_WORK_DIR / "H9_FOV_Physiology_Atlas_Viewer.html", REMOTE_OUT_DIR / "H9_FOV_Physiology_Atlas_Viewer.html")
    shutil.copy2(LOCAL_WORK_DIR / "index.html", REMOTE_OUT_DIR / "index.html")
    shutil.copy2(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv", REMOTE_OUT_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv")
    shutil.copy2(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json", REMOTE_OUT_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json")
    shutil.copy2(Path(__file__), REMOTE_OUT_DIR / Path(__file__).name)

    remote_fig_dir = REMOTE_OUT_DIR / "figures"
    remote_fig_dir.mkdir(parents=True, exist_ok=True)
    fullres_fig_map = {
        "Fig34_H9_Global_Sail_StemCell_Centered_8x7_Hybrid_with_Zoomed_Bifurcating_Insets.png": "Fig34_Global_2D_Sail_Trajectories.png",
        "Fig40_H9_DatasetWide_419FOV_Spatial_Wave_Coherence_and_Tissue_Remodeling.png": "Fig40_Spatial_Wave_Coherence_419FOVs.png",
        "Fig41_H9_4Channel_FOV_GreenStill_FMBarcode_StateShift_and_ROC_MesoAxis.png": "Fig41_4Channel_GreenStill_FMBarcode_StateShift_ROCMeso.png",
        "Fig42_H9_NonClonal_Barcode_Transcendence_and_DeadZone_Expression_Diagnostics.png": "Fig42_Barcode_Transcendence_and_DeadZone_Diagnostics.png",
        "Fig43_H9_Stoichiometry_Matched_Oscillation_Probability_and_Excess_Quenching_Maps.png": "Fig43_Stoichiometry_Matched_Oscillation_Probability.png",
        "Fig44_H9_ControlCompetent_Gated_dPosc_and_Meso_HyperExcitability_Coupling.png": "Fig44_ControlCompetent_Gated_dPosc_and_Meso_Coupling.png",
        "Fig45_H9_FitFree_Dual_Space_dFstem_vs_dPosc_2D_3D_and_Sail.png": "Fig45_FitFree_Dual_Space_dFstem_vs_dPosc_2D_3D_and_Sail.png",
        "Fig46_H9_StoichiometryBin_FitFree_Dual_Sail_dF_vs_dPosc.png": "Fig46_StoichiometryBin_FitFree_Dual_Sail_dF_vs_dPosc.png",
        "Fig47_H9_Dual_ROC_Surprise_Flat_Control_Blue_Quench_Pink_Meso.png": "Fig47_Dual_ROC_Surprise_Flat_Control_Blue_Quench_Pink_Meso.png",
        "Fig48_H9_Representative_Contact_Sheet_All_Spatial_Views_4Lineages.png": "Fig48_Representative_Contact_Sheet_All_Spatial_Views_4Lineages.png",
    }
    for src_name, dst_name in fullres_fig_map.items():
        src_p = ARCHIVE_DIR / src_name
        if src_p.exists():
            shutil.copy2(src_p, remote_fig_dir / dst_name)
            shutil.copy2(src_p, REMOTE_OUT_DIR / src_name)

    if GEO_REPO_DIR.exists():
        print("Syncing upgraded 420-FOV atlas and HTML viewer to /home/smcoyle/GEO-stem-cell-differentiation-explorer/...", flush=True)
        shutil.copy2(LOCAL_WORK_DIR / "H9_FOV_Physiology_Atlas_Viewer.html", GEO_REPO_DIR / "H9_FOV_Physiology_Atlas_Viewer.html")
        shutil.copy2(LOCAL_WORK_DIR / "index.html", GEO_REPO_DIR / "index.html")
        shutil.copy2(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv", GEO_REPO_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.csv")
        shutil.copy2(LOCAL_WORK_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json", GEO_REPO_DIR / "2026_10_19_stemcell_420_fovs_physiology_atlas_summary.json")
        shutil.copy2(Path(__file__), GEO_REPO_DIR / Path(__file__).name)
        for png_file in ATLAS_DIR.glob("*.png"):
            shutil.copy2(png_file, GEO_REPO_DIR / "fov_atlas" / png_file.name)

    print(f"Complete! All {len(all_fov_records)} FOVs rendered and synced in {time.perf_counter() - t0:.1f}s.", flush=True)


if __name__ == "__main__":
    main()
