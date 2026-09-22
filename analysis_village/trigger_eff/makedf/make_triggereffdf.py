"""
Trigger-efficiency dataframes.

Modeled on analysis_village/hnl_nuee_nupi0 nopreselect + savePfp:
  - slice-led table expanded to all PFPs
  - barycenterFM + correctedOpFlash
  - optional GENIE mcnu truth match via slc.tmatch
  - no analysis preselection in the default config
"""

from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *
import warnings

warnings.filterwarnings("ignore")


def make_mcnudf_triggereff(f, **args):
    mcdf = make_mcnudf(f, **args)
    # drop mcdf columns not relevant for this analysis
    if "mu" in list(zip(*list(mcdf.columns)))[0]:
        mcdf = mcdf.drop("mu", axis=1, level=0)
    if "p" in list(zip(*list(mcdf.columns)))[0]:
        mcdf = mcdf.drop("p", axis=1, level=0)
    if "cpi" in list(zip(*list(mcdf.columns)))[0]:
        mcdf = mcdf.drop("cpi", axis=1, level=0)
    return mcdf


def make_triggereffdf_mcnu_nopreselect_savepfp(f):
    return make_triggereffdf_mcnu(f, applyPreselection=False, savePfp=True)


def make_triggereffdf_mcnu(
    f,
    include_weights=False,
    multisim_nuniv=100,
    slim=True,
    applyPreselection=True,
    savePfp=False,
):
    slcdf = make_triggereffdf(f, applyPreselection=applyPreselection, savePfp=savePfp)
    mcdf = make_mcnudf_triggereff(
        f,
        include_weights=include_weights,
        multisim_nuniv=multisim_nuniv,
        slim=slim,
    )
    mcdf.columns = pd.MultiIndex.from_tuples(
        [tuple(["slc", "truth"] + list(c)) for c in mcdf.columns]
    )
    df = multicol_merge(
        slcdf.reset_index(),
        mcdf.reset_index(),
        left_on=[
            ("entry", "", "", "", "", ""),
            ("slc", "tmatch", "idx", "", "", ""),
        ],
        right_on=[
            ("entry", "", "", "", "", ""),
            ("rec.mc.nu..index", "", ""),
        ],
        how="left",
    )
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df


def make_triggereffdf_data(f, applyPreselection=True, savePfp=False):
    slcdf = make_triggereffdf(f, applyPreselection=applyPreselection, savePfp=savePfp)
    # drop truth cols for data
    slcdf = slcdf.drop("tmatch", axis=1, level=1)  # slc level
    slcdf = slcdf.drop("truth", axis=1, level=2)  # pfp level

    framedf = make_framedf(f)
    timingdf = make_timingdf(f)
    ftdf = multicol_merge(
        framedf,
        timingdf,
        left_index=True,
        right_index=True,
        how="left",
        validate="one_to_one",
    )

    df = multicol_merge(
        slcdf.reset_index(),
        ftdf.reset_index(),
        left_on=[("entry", "", "", "", "", "")],
        right_on=[("entry", "", "", "", "", "")],
        how="left",
    )
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df


def make_triggereffdf_data_nopreselect_savepfp(f):
    return make_triggereffdf_data(f, applyPreselection=False, savePfp=True)


def make_triggereffdf(f, applyPreselection=False, savePfp=True):
    """
    One row per slice, or one row per PFP when savePfp=True.

    Slice branches include barycenterFM and correctedOpFlash. No object
    reduction is performed when savePfp=True; selection is left to analysis.
    """
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    DETECTOR = "SBND"

    pfpdf = make_pfpdf(f)
    slcdf = loadbranches(
        f["recTree"], slcbranches + barycenterFMbranches + correctedflashbranches
    )
    slcdf = slcdf.rec

    pfpdf = pfpdf.drop("pfochar", axis=1, level=1)

    if savePfp:
        # Keep slcdf as the left/base table and expand rows with pfp-level info.
        slcdf = multicol_merge(
            slcdf,
            pfpdf.reset_index(level="rec.slc.reco.pfp..index"),
            left_index=True,
            right_index=True,
            how="left",
            validate="one_to_many",
        )

        pfp_idx_col = next(
            c
            for c in slcdf.columns
            if c == "rec.slc.reco.pfp..index"
            or (
                isinstance(c, tuple)
                and len(c) > 0
                and c[0] == "rec.slc.reco.pfp..index"
            )
        )
        slcdf = slcdf.set_index(pfp_idx_col, append=True)

    if applyPreselection:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic == 0]
        slcdf = slcdf[slcdf.slc.nu_score > 0.5]
        slcdf = slcdf[InFV(df=slcdf.slc.vertex, inzback=0, det=DETECTOR)]

    return slcdf
