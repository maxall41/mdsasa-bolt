# Copyright (C) 2025 Maxwell J. Campbell
import logging
import sys

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PRM, TRJ
from utils import save

from mdsasa_bolt.analysis import SASAAnalysis

# Configure logging to show INFO level messages
logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s - %(message)s")

u = mda.Universe(PRM, TRJ)

analysis = SASAAnalysis(u)
analysis.run()

if len(sys.argv) > 1 and sys.argv[1] == "save":
    save(analysis, "new_sasa_results.txt")
