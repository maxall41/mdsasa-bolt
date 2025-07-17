from pathlib import Path

import MDAnalysis as MDa
import numpy as np
import pytest
from MDAnalysis.core.topologyattrs import Atomnames, Atomtypes, Resids, Resnames, Resnums, Segids
import math
from mdsasa_bolt.analysis import SASAAnalysis

from ..utils import make_universe

PARENT = Path(__file__).parent.parent


class TestSASAAnalysis:
    # fixtures are helpful functions that set up a test
    # See more at https://docs.pytest.org/en/stable/how-to/fixtures.html
    @pytest.fixture
    def universe(self) -> MDa.Universe:
        """Creates dummy universe object with three frames."""
        u = make_universe(
            n_frames=3,
        )

        for ts in u.trajectory:
            ts.positions[: ts.frame] *= -1

        u.add_TopologyAttr(Atomnames(["H"] * len(u.atoms)))
        u.add_TopologyAttr(Resnames(["GLY"] * len(u.residues)))
        u.add_TopologyAttr(Resids(list(range(len(u.residues)))))
        u.add_TopologyAttr(Resnums(list(range(len(u.residues)))))
        u.add_TopologyAttr(Segids(["A"] * len(u.segments)))
        u.add_TopologyAttr(Atomtypes(["O"] * len(u.atoms)))
        return u

    @pytest.fixture
    def analysis(self, universe: MDa.Universe) -> SASAAnalysis:
        """Create analysis class from dummy universe object."""
        return SASAAnalysis(universe)

    @pytest.mark.parametrize(
        "select, n_atoms",  # argument names
        [  # argument values in a tuple, in order
            ("all", 125),
            ("index 0:9", 10),
            ("segindex 3:4", 50),
        ],
    )
    def test_atom_selection(self, universe: MDa.Universe, select, n_atoms) -> None:
        """Test that we can select atoms for analysis."""
        # `universe` here is the fixture defined above
        analysis = SASAAnalysis(universe, select=select)
        assert analysis.atomgroup.n_atoms == n_atoms

    def test_total_sasa_calculation(self, analysis: SASAAnalysis) -> None:
        """Test that analysis pipeline outputs same number of frames in input."""
        analysis.run(stop=3)
        assert analysis.n_frames == 3

    def test_with_atom_group(self, universe: MDa.Universe) -> None:
        analysis = SASAAnalysis(universe.atoms)
        analysis.run()
        assert analysis.n_frames == 3

    def test_total_sasa_calculation_results(self, analysis: SASAAnalysis) -> None:
        analysis.run(stop=3)
        assert analysis.n_frames == 3
        assert analysis.results["total_area"].dtype == np.dtype("float64")
        assert np.all(analysis.results["total_area"] >= 0)

    def test_residue_sasa_calculation_results(self, analysis: SASAAnalysis) -> None:
        analysis.run(stop=3)
        assert analysis.n_frames == 3
        assert analysis.results["residue_area"].dtype == np.dtype("float64")
        assert np.all(analysis.results["residue_area"] >= 0)
        assert analysis.results["residue_area"].shape == (3, 25)

    def test_calculation(self):
        u = MDa.Universe(PARENT / 'data' / 'cobrotoxin.pdb', PARENT / 'data' / 'cobrotoxin.trr')

        selected_atoms = u.select_atoms("not (resname SOL or resname CL or resname NA)")
        filtered_residues = [r for r in u.residues if r.resname not in {"SOL", "CL", "NA"}]
        u.residues = u.residues[np.array([r.ix for r in filtered_residues])]

        analysis = SASAAnalysis(selected_atoms)
        analysis.run()
        assert np.all(analysis.results["residue_area"] >= 0)
        assert math.isclose(analysis.results["total_area"][0],4580.830870509148,abs_tol=20)
        assert math.isclose(analysis.results["total_area"][1],4811.4414303302765,abs_tol=20)
        assert math.isclose(analysis.results["total_area"][2],4685.583839893341,abs_tol=20)
