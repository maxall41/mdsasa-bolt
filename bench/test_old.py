# Copyright (C) 2025 Maxwell J. Campbell
import sys

import MDAnalysis as MDa
import numpy as np
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis
from utils import save

u = MDa.Universe(
    "/Users/maxcampbell/mdsasa-bolt/bench/10827_dyn_85.psf",
    "/Users/maxcampbell/mdsasa-bolt/bench/10824_trj_85.xtc",
)

selected_atoms = u.select_atoms("not (resname TIP3 or resname SOD or resname CLA or resname POPC)")
filtered_residues = [r for r in u.residues if r.resname not in {"TIP3", "CLA", "SOD", "POPC"}]
u.residues = u.residues[np.array([r.ix for r in filtered_residues])]


analysis = SASAAnalysis(selected_atoms, select="all")
analysis.run()

if len(sys.argv) > 1 and sys.argv[1] == "save":
    print("SAVING")
    save(analysis, "old_sasa_results.txt")
