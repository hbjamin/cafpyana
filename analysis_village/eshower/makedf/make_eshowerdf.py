"""
Eshower slice-level dataframe: primary shower + primary track per slice,
with slc (nu_score, barycenterFM), razzled/dazzle PID, and recalc dEdx.
"""
from makedf.makedf import (
    make_slcdf,
    make_mcdf,
    make_potdf_bnb,
    make_hdrdf,
    loadbranches,
)
from makedf.branches import (
    slcbranches,
    barycenterFMbranches,
    trkbranches,
    shwbranches,
    pfpbranch,
)
from pyanalib.pandas_helpers import multicol_merge, multicol_add
from makedf.util import InAV
import pandas as pd
import numpy as np

# Razzled (shower PID) and dazzle (track PID) - per PFP
_razzled_branches = [
    pfpbranch + "razzled.bestScore",
    pfpbranch + "razzled.electronScore",
    pfpbranch + "razzled.muonScore",
    pfpbranch + "razzled.pdg",
    pfpbranch + "razzled.photonScore",
    pfpbranch + "razzled.pionScore",
    pfpbranch + "razzled.protonScore",
]
_dazzle_branches = [
    pfpbranch + "trk.dazzle.bestScore",
    pfpbranch + "trk.dazzle.muonScore",
    pfpbranch + "trk.dazzle.otherScore",
    pfpbranch + "trk.dazzle.pdg",
    pfpbranch + "trk.dazzle.pionScore",
    pfpbranch + "trk.dazzle.protonScore",
]

# Full branch set for eshower PFP (trk + shw + razzled + dazzle)
_eshower_pfp_branches = trkbranches + shwbranches + _razzled_branches + _dazzle_branches

# Fallback: no razzled/dazzle and no bestplane (v10_14_02 and similar CAFs lack rec.slc.reco.pfp.shw.bestplane)
_eshower_trkbranches = [b for b in trkbranches if "bestplane" not in b]
_eshower_shwbranches = [b for b in shwbranches if "bestplane" not in b]
_eshower_pfp_branches_fallback = _eshower_trkbranches + _eshower_shwbranches


def _find_col(df, name):
    """Return first column tuple whose last non-empty level equals name."""
    for c in df.columns:
        t = c if isinstance(c, tuple) else (c,)
        if t[-1] == name or (len(t) > 1 and t[-2] == name):
            return c
    return None


def make_eshowerdf(f):
    """Build slice-level df: one row per slice with slc, primshw, primtrk (MultiIndex columns)."""
    # Slice-level: slc + barycenterFM
    slcdf = make_slcdf(f)
    if slcdf.empty:
        return _empty_eshower_df()

    bcdf = loadbranches(f["recTree"], barycenterFMbranches)
    bcdf = bcdf.rec.slc.barycenterFM
    _flat = [c if isinstance(c, str) else c[-1] or c[0] for c in bcdf.columns]
    bcdf.columns = pd.MultiIndex.from_tuples(
        [("slc", "barycenterFM", c, "", "", "") for c in _flat]
    )
    slcdf = multicol_merge(slcdf, bcdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    if slcdf.empty:
        return _empty_eshower_df()

    # PFP with trk, shw, razzled, dazzle
    try:
        pfpdf = loadbranches(f["recTree"], _eshower_pfp_branches)
    except Exception:
        # Fallback if razzled/dazzle or bestplane missing (e.g. v10_14_02 CAFs lack pfp.shw.bestplane)
        pfpdf = loadbranches(f["recTree"], _eshower_pfp_branches_fallback)
    pfpdf = pfpdf.rec.slc.reco

    # Maxplane and maxplane_energy for showers (same logic as make_pfpdf)
    try:
        nh = pfpdf.loc(axis=1)["pfp", "shw", "plane", :, "nHits"]
        pfpdf[("pfp", "shw", "maxplane", "", "", "")] = nh.idxmax(axis=1).apply(
            lambda x: x[3] if isinstance(x, tuple) and len(x) > 3 else "I2"
        )
        mp = pfpdf[("pfp", "shw", "maxplane", "", "", "")]
        pfpdf[("pfp", "shw", "maxplane_energy", "", "", "")] = np.select(
            [mp == "I2", mp == "I1", mp == "I0"],
            [
                pfpdf[("pfp", "shw", "plane", "I2", "energy", "")],
                pfpdf[("pfp", "shw", "plane", "I1", "energy", "")],
                pfpdf[("pfp", "shw", "plane", "I0", "energy", "")],
            ],
            default=np.nan,
        )
    except Exception:
        pfpdf[("pfp", "shw", "maxplane", "", "", "")] = "I2"
        ene_col = ("pfp", "shw", "plane", "I2", "energy", "")
        pfpdf[("pfp", "shw", "maxplane_energy", "", "", "")] = pfpdf[ene_col] if ene_col in pfpdf.columns else np.nan

    # Best-plane dEdx -> maxplane_dEdx_new (from plane dEdx by maxplane; -999 if invalid)
    try:
        mp = pfpdf[("pfp", "shw", "maxplane", "", "", "")]
        pfpdf[("pfp", "shw", "maxplane_dEdx_new", "", "", "")] = np.select(
            [mp == "I2", mp == "I1", mp == "I0"],
            [
                pfpdf[("pfp", "shw", "plane", "I2", "dEdx", "")],
                pfpdf[("pfp", "shw", "plane", "I1", "dEdx", "")],
                pfpdf[("pfp", "shw", "plane", "I0", "dEdx", "")],
            ],
            default=-999.0,
        )
        invalid = np.isnan(pfpdf[("pfp", "shw", "maxplane_dEdx_new", "", "", "")]) | (
            pfpdf[("pfp", "shw", "maxplane_dEdx_new", "", "", "")] <= 0
        )
        pfpdf.loc[invalid, ("pfp", "shw", "maxplane_dEdx_new", "", "", "")] = -999.0
    except Exception:
        pfpdf[("pfp", "shw", "maxplane_dEdx_new", "", "", "")] = -999.0

    # Primary shower: shower with highest maxplane_energy per slice (trackScore < 0.5)
    tscore_col = _find_col(pfpdf, "trackScore")
    is_shower = (pfpdf[tscore_col] < 0.5) if tscore_col is not None else pd.Series(False, index=pfpdf.index)
    shwdf = pfpdf.loc[is_shower].copy()
    if shwdf.empty:
        slcdf = _add_empty_primshw_primtrk(slcdf)
        return _build_final_eshower(slcdf)

    # Sort by maxplane_energy descending, take first per (entry, slc)
    ene_col = ("pfp", "shw", "maxplane_energy", "", "", "")
    if ene_col not in shwdf.columns:
        shwdf[ene_col] = np.nan
    primshw_df = (
        shwdf.sort_values(ene_col, ascending=False)
        .groupby(level=[0, 1], sort=False)
        .nth(0)
        .copy()
    )
    # Drop track-only columns so we keep slc-like + shw + razzled + trackScore
    keep_shw = [c for c in primshw_df.columns if c[1] != "trk"]
    primshw_df = primshw_df[[c for c in primshw_df.columns if c in keep_shw]]
    primshw_df.columns = pd.MultiIndex.from_tuples(
        [("primshw",) + c[1:] for c in primshw_df.columns],
        names=primshw_df.columns.names,
    )

    # Primary track: longest track per slice (trackScore >= 0.5)
    is_track = (pfpdf[tscore_col] >= 0.5) if tscore_col is not None else pd.Series(False, index=pfpdf.index)
    trkdf = pfpdf.loc[is_track].copy()
    if trkdf.empty:
        slcdf = multicol_merge(slcdf, primshw_df.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
        slcdf = _add_empty_primtrk(slcdf)
        return _build_final_eshower(slcdf)

    len_col = ("pfp", "trk", "len", "", "", "")
    if len_col not in trkdf.columns:
        trkdf[len_col] = 0.0
    primtrk_df = (
        trkdf.sort_values(len_col, ascending=False)
        .groupby(level=[0, 1], sort=False)
        .nth(0)
        .copy()
    )
    keep_trk = [c for c in primtrk_df.columns if c[1] != "shw"]
    primtrk_df = primtrk_df[[c for c in primtrk_df.columns if c in keep_trk]]
    primtrk_df.columns = pd.MultiIndex.from_tuples(
        [("primtrk",) + c[1:] for c in primtrk_df.columns],
        names=primtrk_df.columns.names,
    )

    slcdf = multicol_merge(slcdf, primshw_df.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    slcdf = multicol_merge(slcdf, primtrk_df.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    return _build_final_eshower(slcdf)


def _empty_eshower_df():
    return pd.DataFrame()


def _add_empty_primshw_primtrk(slcdf):
    """Add placeholder primshw/primtrk columns (NaN) when no shower/track."""
    for tag in ("primshw", "primtrk"):
        for stub in [
            ("shw", "maxplane_energy", "", "", ""),
            ("shw", "maxplane_dEdx_new", "", "", ""),
            ("shw", "len", "", "", ""),
            ("trk", "len", "", "", ""),
        ]:
            c = (tag,) + stub
            if c not in slcdf.columns:
                slcdf[c] = np.nan
    return slcdf


def _add_empty_primtrk(slcdf):
    for stub in [
        ("primtrk", "trk", "len", "", "", ""),
        ("primtrk", "trackScore", "", "", "", ""),
    ]:
        if stub not in slcdf.columns:
            slcdf[stub] = np.nan
    return slcdf


def _build_final_eshower(slcdf):
    """Ensure 6-level MultiIndex columns and return."""
    nlevel = 6
    def pad(c):
        if isinstance(c, str):
            c = (c,)
        return tuple(list(c) + [""] * (nlevel - len(c)))
    if slcdf.columns.nlevels < nlevel:
        slcdf.columns = pd.MultiIndex.from_tuples(
            [pad(c) for c in slcdf.columns],
            names=(list(slcdf.columns.names or []) + [""] * nlevel)[:nlevel],
        )
    return slcdf
