# Copyright (C) 2025 Maxwell J. Campbell
import sys

import MDAnalysis as MDa
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis
from utils import save

u = MDa.Universe(
    "/Users/maxcampbell/mdsasa-bolt/bench/10827_dyn_85.psf",  # 4IAQ: https://www.gpcrmd.org/view/85/
    "/Users/maxcampbell/mdsasa-bolt/bench/10824_trj_85.xtc",
)

selected_atoms = u.select_atoms("protein")


analysis = SASAAnalysis(selected_atoms)
analysis.run()

if len(sys.argv) > 1 and sys.argv[1] == "save":
    save(analysis, "mdakit_results.txt")
