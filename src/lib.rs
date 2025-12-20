use std::collections::HashMap;

use pyo3::prelude::*;
mod simd;
use nalgebra::Point3;
use rayon::prelude::*;
use rust_sasa::Atom;
use rust_sasa::calculate_sasa_internal as calculate_sasa_internal_internal;
use simd::simd_sum;

#[pyclass]
#[derive(Clone)]
pub struct Residue {
    #[pyo3(get)]
    pub chain_id: String,
    #[pyo3(get)]
    pub residue_name: String,
    #[pyo3(get)]
    pub residue_number: u32,
    #[pyo3(get)]
    pub sasa: f32,
}

type ProteinAtoms = Vec<((f32, f32, f32), f32, usize)>;

pub fn calculate_sasa_internal(
    atoms_in: Vec<((f32, f32, f32), f32, usize)>,
    probe_radius: f32,
    n_points: usize,
) -> PyResult<Vec<f32>> {
    let atoms: Vec<Atom> = atoms_in
        .into_iter()
        .enumerate()
        .map(|(index, (pos, radius, parent_id))| Atom {
            position: [pos.0, pos.1, pos.2],
            id: index,
            parent_id: Some(parent_id as isize),
            radius,
        })
        .collect();
    Ok(calculate_sasa_internal_internal(
        atoms.as_slice(),
        probe_radius,
        n_points,
        1,
    ))
}

fn calculate_sasa_internal_at_residue_level(
    atoms_in: &ProteinAtoms,
    probe_radius: f32,
    n_points: usize,
) -> PyResult<Vec<Residue>> {
    let atom_sasa = calculate_sasa_internal(atoms_in.clone(), probe_radius, n_points).unwrap();

    let mut residue_groups: HashMap<usize, Vec<f32>> = HashMap::new();

    for (i, atom) in atoms_in.iter().enumerate() {
        let parent_id = atom.2;
        residue_groups
            .entry(parent_id)
            .or_insert_with(Vec::new)
            .push(atom_sasa[i]);
    }

    let mut residue_sasa = Vec::new();
    for (parent_id, sasa_values) in residue_groups {
        let local_sum = simd_sum(sasa_values.as_slice());

        residue_sasa.push(Residue {
            chain_id: "UNK".to_string(),
            residue_name: "UNK".to_string(),
            residue_number: parent_id as u32,
            sasa: local_sum,
        });
    }

    Ok(residue_sasa)
}

#[pyfunction]
pub fn frames(
    frames: Vec<ProteinAtoms>,
    probe_radius: f32,
    n_points: usize,
) -> PyResult<Vec<Vec<Residue>>> {
    let result: Result<Vec<Vec<Residue>>, _> = frames
        .par_iter()
        .map(|frame| calculate_sasa_internal_at_residue_level(frame, probe_radius, n_points))
        .collect();

    result
}

/// A Python module implemented in Rust.
#[pymodule]
fn plumber(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(frames, m)?)?;
    Ok(())
}
