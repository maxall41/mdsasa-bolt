# Copyright (C) 2025 Maxwell J. Campbell
import logging
import sys

import MDAnalysis as mda
import numpy as np
from MDAnalysis.tests.datafiles import PDB_xvf, TRR_xvf
from utils import save

from mdsasa_bolt.analysis import SASAAnalysis

# Configure logging to show INFO level messages
logging.basicConfig(level=logging.INFO, format="%(name)s - %(levelname)s - %(message)s")

u = mda.Universe(PDB_xvf, TRR_xvf)

selected_atoms = u.select_atoms("not (resname SOL or resname CL or resname NA)")
filtered_residues = [r for r in u.residues if r.resname not in {"SOL", "CL", "NA"}]
u.residues = u.residues[np.array([r.ix for r in filtered_residues])]

analysis = SASAAnalysis(selected_atoms, select="all")
analysis.run()

if len(sys.argv) > 1 and sys.argv[1] == "save":
    save(analysis, "new_sasa_results.txt")
