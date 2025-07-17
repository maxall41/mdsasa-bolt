from MDAnalysis.guesser.tables import vdwradii


def get_all_radii_methods(classifier):
    methods = [
        lambda atom: vdwradii.get(atom.type),
        lambda atom: vdwradii.get(atom.name[0]),
        lambda atom: vdwradii.get(atom.type[0]),
        lambda atom: classifier.radius(atom.resname, atom.name),
        lambda atom: classifier.radius("ANY", atom.type),
    ]

    return methods


def get_all_element_methods():
    methods = [
        lambda atom: atom.element,
        lambda atom: atom.type[0],
        lambda atom: atom.name[0],
    ]

    return methods


def get_atom_element(atom):
    el_methods = get_all_element_methods()
    for method in el_methods:
        try:
            element = method(atom)
            if element:
                return element
        except:
            pass
    raise ValueError(f"Could not determine element for atom {atom}")
