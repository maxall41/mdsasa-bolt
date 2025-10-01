# Copyright (C) 2025 Maxwell J. Campbell
import logging
import sys

import MDAnalysis as MDa
from utils import save

from mdsasa_bolt.analysis import SASAAnalysis

# Configure logging to show INFO level messages
logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s - %(message)s")

u = MDa.Universe(
    "/Users/maxcampbell/mdsasa-bolt/bench/10827_dyn_85.psf",
    "/Users/maxcampbell/mdsasa-bolt/bench/10824_trj_85.xtc",
)

selected_atoms = u.select_atoms("protein")

analysis = SASAAnalysis(selected_atoms)
analysis.run()

if len(sys.argv) > 1 and sys.argv[1] == "save":
    save(analysis, "new_sasa_results.txt")
