import logging
from typing import TYPE_CHECKING, Union

import freesasa
import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.core.selection import NoDataError
from MDAnalysis.guesser.tables import vdwradii
from rust_sasa_python import Residue, calculate_sasa_internal_at_residue_level

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
    ):
        super().__init__(universe_or_atomgroup.universe.trajectory, **kwargs)
        self.universe = universe_or_atomgroup.universe
        self.atomgroup: AtomGroup = universe_or_atomgroup.select_atoms(select)
        self._classifier = freesasa.Classifier()

    def _prepare(self):
        self.results.total_area = np.zeros(
            self.n_frames,
            dtype=float,
        )
        self.results.residue_area = np.zeros(
            (self.n_frames, len(self.universe.residues.resids)),
            dtype=float,
        )

    def _single_frame(self):
        """Calculate data from a single frame of trajectory"""
        input_atoms = []
        for i, atom in enumerate(self.atomgroup):
            x, y, z = atom.position

            try:
                radius = vdwradii.get(atom.element)
            except NoDataError:
                try:
                    radius = self._classifier.radius(atom.resname, atom.name)
                except NoDataError:
                    logger.warning(
                        "Using freesasa classifier with ANY residue type! This may result in incorrect radii and SASA values!",
                    )
                    radius = self._classifier.radius("ANY", atom.type)

            index = atom.resnum.item()
            if radius is None:
                raise ValueError(f"Failed to get radius for atom type: {atom.type}")
            input_atoms.append(((x, y, z), radius, index))

        residue_sasa_values: list[Residue] = calculate_sasa_internal_at_residue_level(input_atoms, 1.4, 100)

        self.results.total_area[self._frame_index] = sum([v.sasa for v in residue_sasa_values])

        # Defend agains residue counts mismatch
        if len(self.universe.residues.resids) != len(residue_sasa_values):
            logger.error(
                f"Residude count do not match the expectation, residue SASA not in results { len(self.universe.residues.resids)} != {len(residue_sasa_values)}",
            )
        else:
            self.results.residue_area[self._frame_index] = [r.sasa for r in residue_sasa_values]

    def _conclude(self):
        self.results.mean_total_area = self.results.total_area.mean()
