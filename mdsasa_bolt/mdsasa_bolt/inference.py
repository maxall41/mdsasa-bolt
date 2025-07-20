# Copyright (C) 2025 Maxwell J. Campbell
from collections.abc import Callable

from freesasa import Classifier
from MDAnalysis.core.groups import Atom
from MDAnalysis.guesser.tables import vdwradii


def get_all_radii_methods(classifier: Classifier) -> list[Callable]:
    """Returns a list of methods to determine the radius of an atom.

    Args:
        classifier: The classifier to use for determining the radius.

    Returns:
        A list of methods to determine the radius of an atom.

    """
    return [
        lambda atom: vdwradii.get(atom.type),
        lambda atom: vdwradii.get(atom.name[0]),
        lambda atom: vdwradii.get(atom.type[0]),
        lambda atom: classifier.radius(atom.resname, atom.name),
        lambda atom: classifier.radius("ANY", atom.type),
    ]


def get_all_element_methods() -> list[Callable]:
    """Returns a list of methods to determine the element of an atom.

    Returns:
        A list of methods to determine the element of an atom.

    """
    return [
        lambda atom: atom.element,
        lambda atom: atom.type[0],
        lambda atom: atom.name[0],
    ]


def get_atom_element(atom: Atom) -> str:
    """Attempts to determine the element of an atom.

    Args:
        atom: The atom to determine the element of.

    Returns:
        The element of the atom.

    Raises:
        ValueError: If the element could not be determined.

    """
    el_methods = get_all_element_methods()
    for method in el_methods:
        try:
            element = method(atom)
            if element:
                return element
        except:
            pass
    raise ValueError(f"Could not determine element for atom {atom}")
