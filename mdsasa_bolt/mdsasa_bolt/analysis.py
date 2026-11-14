# Copyright (C) 2025 Maxwell J. Campbell
import logging
from collections.abc import Callable
from typing import TYPE_CHECKING, Union

import freesasa
import numpy as np
from MDAnalysis import NoDataError
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.core.groups import Atom

if TYPE_CHECKING:
    from MDAnalysis.core.universe import AtomGroup, Universe

from . import plumber
from .inference import get_all_radii_methods

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
        radii: Union[dict[str, float], list[float], "np.ndarray", None] = None,
        **kwargs,
    ) -> None:
        """Initialize SASAAnalysis.

        Parameters
        ----------
        universe_or_atomgroup : Universe or AtomGroup
            Universe or group of atoms to apply this analysis to.
        select : str
            Selection string for atoms to extract.
        radii : dict[str, float] or array-like of float, optional
            Custom atomic radii. Can be either:
            - A dict mapping atom types to radii (e.g. ``{"BB": 2.64, "SC1": 2.3}``).
            - A list/array of per-atom radii with length matching the number of
              selected atoms.
            When ``None`` (default), radii are inferred automatically using FreeSASA.

        """
        super().__init__(universe_or_atomgroup.universe.trajectory, **kwargs)

        self.universe = universe_or_atomgroup

        self.atomgroup: AtomGroup = universe_or_atomgroup.select_atoms(select)

        self.probe_radius = kwargs.get("probe_radius", 1.4)
        self.n_points = kwargs.get("n_points", 100)

        if radii is not None:
            self._atom_radii = self._resolve_custom_radii(radii)
        else:
            self._classifier = freesasa.Classifier().getStandardClassifier("protor")
            self._radius_method = self._determine_radius_method()
            self._atom_radii = self._calculate_atom_radii(self._radius_method)

        self._atom_resnums = np.array([atom.resnum.item() for atom in self.atomgroup])

    def _resolve_custom_radii(self, radii: Union[dict[str, float], list[float], "np.ndarray"]) -> np.ndarray:
        """Resolve user-provided radii into a per-atom numpy array."""
        n_atoms = len(self.atomgroup)

        if isinstance(radii, dict):
            result = np.zeros(n_atoms, dtype=float)
            for i, atom in enumerate(self.atomgroup):
                atom_type = atom.type
                if atom_type not in radii:
                    msg = f"Custom radii dict missing entry for atom type '{atom_type}' (atom index {i})"
                    raise ValueError(msg)
                result[i] = radii[atom_type]
            return result

        result = np.asarray(radii, dtype=float)
        if result.shape != (n_atoms,):
            msg = f"Custom radii array length ({result.shape}) does not match number of selected atoms ({n_atoms})"
            raise ValueError(msg)
        return result

    def _determine_radius_method(self) -> Callable:
        """Determine the best radius calculation method for this system."""
        # Try each method with the first 3 atoms to see which one works
        x = 3
        test_atoms: AtomGroup = self.atomgroup[:x]
        logger.info("Testing radius calculation methods for first %d atoms", x)

        for method in get_all_radii_methods(self._classifier):
            worked = 0
            for test_atom in test_atoms:
                try:
                    radius = method(test_atom)
                    if radius is None or radius <= 0:
                        msg = "Invalid radius"
                        raise NoDataError(msg)  # noqa: TRY301
                    worked += 1
                except NoDataError:
                    pass
            if worked == x:
                return method

        error_msg = "No radius calculation method worked for this system"
        raise ValueError(error_msg)

    def _get_radius_with_fallback(self, atom: Atom, method: Callable) -> float:
        """Get radius for an atom with fallback methods if primary method fails."""
        try:
            return method(atom)
        except NoDataError:
            for fallback_method in get_all_radii_methods(self._classifier):
                try:
                    radius = fallback_method(atom)
                    if radius is not None and radius > 0:
                        return radius
                except (NoDataError, Exception):
                    pass

        error_msg = "No radius calculation method worked for this system"
        raise ValueError(error_msg)

    def _calculate_atom_radii(self, method: Callable) -> np.ndarray:
        """Calculate radii for all atoms using the determined method."""
        radii = np.zeros(len(self.atomgroup), dtype=float)

        for i, atom in enumerate(self.atomgroup):
            radii[i] = self._get_radius_with_fallback(atom, method)

        logger.info(f"Pre-computed radii for {len(radii)} atoms")
        return radii

    def run(
        self,
        start: int = 0,
        stop: None | int = None,
        step: int = 1,
        frames: None | list[int] = None,
    ) -> None:
        """Run the analysis."""
        # Update frame parameters if provided
        self.start = start
        self.stop = stop
        self.step = step
        self.frames = frames

        self._setup_frames(
            self._trajectory,
            self.start,
            self.stop,
            self.step,
            self.frames,
        )

        input_atoms_per_frame = []

        # Iterate over trajectory (which now respects the frame parameters)
        for _ in self._sliced_trajectory:
            input_atoms = [
                (tuple(position), radius, resnum)
                for position, radius, resnum in zip(
                    self.atomgroup.positions.copy(),
                    self._atom_radii,
                    self._atom_resnums,
                    strict=False,
                )
            ]
            input_atoms_per_frame.append(input_atoms)

        self.stop = len(input_atoms_per_frame)

        self.results.total_area = np.zeros(
            self.n_frames,
            dtype=float,
        )
        self.results.residue_area = np.zeros(
            (self.n_frames, len(self.universe.residues.resids)),
            dtype=float,
        )

        frame_residues = plumber.frames(input_atoms_per_frame, self.probe_radius, self.n_points)

        for frame_index, residues in enumerate(frame_residues):
            self.results.total_area[frame_index] = sum([v.sasa for v in residues])
            if len(self.universe.residues.resids) != len(residues):
                logger.error(
                    f"Residue count does not match the expectation! Not saving per residue SASA data! universe: {len(self.universe.residues.resids)}, frame: {len(residues)}",  # noqa: E501
                )
            else:
                self.results.residue_area[frame_index] = [r.sasa for r in residues]

    def _conclude(self) -> None:
        self.results.mean_total_area = self.results.total_area.mean()
