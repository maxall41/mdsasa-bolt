import warnings
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure
from Bio.PDB.Atom import PDBConstructionWarning


def save(analysis, filename) -> None:
    """Save analysis results by dumping frames as PDB with SASA in B-factor column."""
    with open("bench/" + filename, "w") as f:
        f.writelines(f"{frame}\n" for frame in analysis.results.total_area)
    # Create output directory
    output_dir = Path("sasa_frames")
    output_dir.mkdir(exist_ok=True)

    # Get universe and atomgroup from analysis
    universe = analysis.universe
    atomgroup = analysis.atomgroup

    # Suppress Biopython warnings
    warnings.filterwarnings("ignore", category=PDBConstructionWarning)

    # Process each frame
    for frame_idx in range(analysis.n_frames):
        # Move to the specific frame
        universe.trajectory[analysis.start + frame_idx * analysis.step]

        # Get residue SASA values for this frame
        residue_sasa = analysis.results.residue_area[frame_idx]

        # Map residue SASA to atoms
        atom_sasa = np.zeros(len(atomgroup))
        for i, atom in enumerate(atomgroup):
            # Use the residue index to get SASA value
            res_idx = atom.resindex
            if res_idx < len(residue_sasa):
                atom_sasa[i] = residue_sasa[res_idx]

        # print("atom_sasa", atom_sasa)

        # Create Biopython structure
        structure = Structure.Structure("SASA_STRUCTURE")
        model = Model.Model(0)
        structure.add(model)

        # Group atoms by chain and residue
        current_chain_id = None
        current_chain = None
        current_resnum = None
        current_residue = None

        for i, atom in enumerate(atomgroup):
            # Handle chain
            chain_id = getattr(atom, "chainid", "A") if hasattr(atom, "chainid") else "A"
            if chain_id != current_chain_id:
                current_chain = Chain.Chain(chain_id)
                model.add(current_chain)
                current_chain_id = chain_id
                current_resnum = None  # Reset residue tracking

            # Handle residue
            resnum = atom.resnum
            if resnum != current_resnum:
                resname = atom.resname
                current_residue = Residue.Residue((" ", resnum, " "), resname, " ")
                current_chain.add(current_residue)
                current_resnum = resnum

            # Create atom
            coord = atom.position
            bfactor = atom_sasa[i]
            occupancy = getattr(atom, "occupancy", 1.0) if hasattr(atom, "occupancy") else 1.0
            element = getattr(atom, "element", atom.name[0]) if hasattr(atom, "element") else atom.name[0]

            bio_atom = Atom.Atom(
                name=atom.name,
                coord=coord,
                bfactor=bfactor,
                occupancy=occupancy,
                altloc=" ",
                fullname=atom.name.ljust(4),
                serial_number=i + 1,
                element=element,
            )
            current_residue.add(bio_atom)

        # for atom in structure.get_atoms():
        #     print("AS", atom.bfactor)

        # Write PDB file
        output_file = output_dir / f"frame_{frame_idx:04d}_sasa.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_file))

