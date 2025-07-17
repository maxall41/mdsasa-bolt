# Copyright (C) 2025 Maxwell J. Campbell
import hashlib
import logging
import pickle
from collections.abc import Callable
from typing import TYPE_CHECKING, Union

import freesasa
import numpy as np
from MDAnalysis import NoDataError
from MDAnalysis.analysis.base import AnalysisBase
from rust_sasa_python import Residue, calculate_sasa_internal_at_residue_level

from .inference import get_all_radii_methods, get_atom_element

if TYPE_CHECKING:
    from MDAnalysis.core.universe import AtomGroup, Universe

logger = logging.getLogger(__name__)


class SASAAnalysis(AnalysisBase):
    """SASAAnalysis class.

    This class is used to compute the solvant accessible area of a trajectory.

    Parameters
    ----------
    universe_or_atomgroup: :class:`~MDAnalysis.core.universe.Universe` or :class:`~MDAnalysis.core.groups.AtomGroup`
        Universe or group of atoms to apply this analysis to.
        If a trajectory is associated with the atoms,
        then the computation iterates over the trajectory.
    select: str
        Selection string for atoms to extract from the input Universe or
        AtomGroup

    Attributes
    ----------
    universe: :class:`~MDAnalysis.core.universe.Universe`
        The universe to which this analysis is applied
    atomgroup: :class:`~MDAnalysis.core.groups.AtomGroup`
        The atoms to which this analysis is applied
    results: :class:`~MDAnalysis.analysis.base.Results`
        results of calculation are stored here, after calling
        :meth:`SASAAnalysis.run`
    start: Optional[int]
        The first frame of the trajectory used to compute the analysis
    stop: Optional[int]
        The frame to stop at for the analysis
    step: Optional[int]
        Number of frames to skip between each analyzed frame
    n_frames: int
        Number of frames analysed in the trajectory
    times: numpy.ndarray
        array of Timestep times. Only exists after calling
        :meth:`SASAAnalysis.run`
    frames: numpy.ndarray
        array of Timestep frame indices. Only exists after calling
        :meth:`SASAAnalysis.run`

    """

    def __init__(
        self,
        universe_or_atomgroup: Union["Universe", "AtomGroup"],
        select: str = "all",
        **kwargs,
    ) -> None:
        """Initialize SASAAnalysis."""
        super().__init__(universe_or_atomgroup.universe.trajectory, **kwargs)
        self.universe = universe_or_atomgroup.universe
        self.atomgroup: AtomGroup = universe_or_atomgroup.select_atoms(select)
        self._classifier = freesasa.Classifier().getStandardClassifier("protor")

        # Determine the best radius calculation method for this system
        self._radius_method = self._determine_radius_method()

        # Pre-compute radii for all atoms using the determined method
        self._atom_radii = self._calculate_atom_radii()
        self._atom_resnums = np.array([atom.resnum.item() for atom in self.atomgroup])

    def _determine_radius_method(self) -> Callable:
        """Determine the best radius calculation method for this system."""
        # Try each method with the first 3 atoms to see which one works
        x = 3
        test_atoms = self.atomgroup[:x]
        logger.info("Testing radius calculation methods for first %d atoms", x)

        for method in get_all_radii_methods(self._classifier):
            for test_atom in test_atoms:
                try:
                    radius = method(test_atom)
                    if radius is None or radius <= 0:
                        raise NoDataError("Invalid radius")
                    return method
                except NoDataError:
                    pass

        error_msg = "No radius calculation method worked for this system"
        raise ValueError(error_msg)

    def _get_radius_with_fallback(self, atom) -> float:
        """Get radius for an atom with fallback methods if primary method fails."""
        for method in get_all_radii_methods(self._classifier):
            try:
                radius = method(atom)
                logger.info(f"Using fallback radius method for atom {atom.name} in residue {atom.resname} = {radius}")
                if radius is not None and radius > 0:
                    return radius
            except (NoDataError, Exception):
                pass

        error_msg = "No radius calculation method worked for this system"
        raise ValueError(error_msg)

    def _calculate_atom_radii(self) -> np.ndarray:
        """Calculate radii for all atoms using the determined method."""
        radii = np.zeros(len(self.atomgroup), dtype=float)

        for i, atom in enumerate(self.atomgroup):
            radii[i] = self._get_radius_with_fallback(atom)

        logger.info(f"Pre-computed radii for {len(radii)} atoms")
        return radii

    def _prepare(self) -> None:
        self.results.total_area = np.zeros(
            self.n_frames,
            dtype=float,
        )
        self.results.residue_area = np.zeros(
            (self.n_frames, len(self.universe.residues.resids)),
            dtype=float,
        )

    def _single_frame(self) -> None:
        """Calculate data from a single frame of trajectory."""
        # Get current positions and construct input_atoms efficiently
        positions = self.atomgroup.positions

        input_atoms = [
            (tuple(position), radius, resnum)
            for position, radius, resnum, atom in zip(
                positions,
                self._atom_radii,
                self._atom_resnums,
                self.atomgroup,
                strict=False,
            )
            if get_atom_element(atom) != "H"
        ]

        # Calculate and print md5 hash of input_atoms
        md5_hash = hashlib.md5(pickle.dumps(input_atoms)).hexdigest()
        print(
            f"MD5 hash of input_atoms: {md5_hash}",
        )  # MATCH!!! (0351e1fe72a5d39aed5dde0e8d7c3992, d179a1ce1682e00b157b2aab3318b9ea, e3e9c1a743eaedc31db3d33d7330e8bc) -> (0351e1fe72a5d39aed5dde0e8d7c3992,d179a1ce1682e00b157b2aab3318b9ea,e3e9c1a743eaedc31db3d33d7330e8bc)
        # 11,671 -> 6,629

        residue_sasa_values: list[Residue] = calculate_sasa_internal_at_residue_level(input_atoms, 1.4, 100)

        md5_hash = hashlib.md5(pickle.dumps([r.sasa for r in residue_sasa_values])).hexdigest()
        print(
            f"MD5 hash of residue_sasa_values: {md5_hash}",
        )  # NO MATCH!!! (a3679e0ce19bd353899880974f2a5a43,be3945355682df955a1b234703c256e9,56fc2dff0dce7b4ed9d6794a5361a60b) -> (ca076f949041bb20495f8866e6a105e7,6a7d17d19c414e012d9f3a2d759fff50,723cd476889f52058cafc34d1de0b3e0)

        for sasa in residue_sasa_values:
            logger.info(f"Residue {sasa.residue_number} {sasa.residue_name} has SASA {sasa.sasa}")

        self.results.total_area[self._frame_index] = sum([v.sasa for v in residue_sasa_values])

        # Defend against residue counts mismatch
        if len(self.universe.residues.resids) != len(residue_sasa_values):
            logger.error(
                "Residue count does not match the expectation",
            )
        else:
            self.results.residue_area[self._frame_index] = [r.sasa for r in residue_sasa_values]

    def _conclude(self) -> None:
        self.results.mean_total_area = self.results.total_area.mean()
