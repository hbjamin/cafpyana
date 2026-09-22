from analysis_village.trigger_eff.makedf.make_triggereffdf import *

DFS = [
    make_triggereffdf_mcnu_nopreselect_savepfp,
    make_hdrdf,
    make_triggerdf,
    make_potdf_bnb,
    make_mcnudf,
]
NAMES = ["rec", "hdr", "trigger", "pot", "mcnu"]
