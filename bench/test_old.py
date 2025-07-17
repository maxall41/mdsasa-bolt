# Copyright (C) 2025 Maxwell J. Campbell
import sys

import MDAnalysis as mda
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis
from MDAnalysis.tests.datafiles import PRM, TRJ
from utils import save

u = mda.Universe(PRM, TRJ)

analysis = SASAAnalysis(u)
analysis.run()

if len(sys.argv) > 1 and sys.argv[1] == "save":
    save(analysis, "old_sasa_results.txt")
