from analysis_village.eshower.makedf.make_eshowerdf import make_eshowerdf
from makedf.makedf import make_hdrdf, make_potdf_bnb, make_mcnudf

DFS = [make_eshowerdf, make_hdrdf, make_potdf_bnb, make_mcnudf]
NAMES = ["eshower", "hdr", "pot", "mcnu"]
