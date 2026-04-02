//! PyO3 bindings for AMBER force field and energy minimization.

use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;

use crate::forcefield::{energy, minimize, params, topology};
use crate::py_pdb::PyPDB;

/// Compute AMBER force field energy of a structure.
///
/// Returns dict with energy components:
///   bond_stretch, angle_bend, torsion, vdw, electrostatic, total
///   (all in kcal/mol)
#[pyfunction]
pub fn compute_energy(py: Python<'_>, pdb: &PyPDB) -> PyResult<PyObject> {
    let amber = params::amber96();
    let topo = topology::build_topology(&pdb.inner, &amber);
    let coords: Vec<[f64; 3]> = topo.atoms.iter().map(|a| a.pos).collect();

    let result = py.allow_threads(|| energy::compute_energy(&coords, &topo, &amber));

    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("bond_stretch", result.bond_stretch)?;
    dict.set_item("angle_bend", result.angle_bend)?;
    dict.set_item("torsion", result.torsion)?;
    dict.set_item("vdw", result.vdw)?;
    dict.set_item("electrostatic", result.electrostatic)?;
    dict.set_item("total", result.total)?;
    Ok(dict.into_any().unbind())
}

/// Minimize hydrogen positions using AMBER force field.
///
/// Freezes all heavy atoms and optimizes only hydrogen positions
/// using steepest descent with adaptive step size.
///
/// Args:
///     pdb: Structure to minimize.
///     max_steps: Maximum optimization steps (default 500).
///     gradient_tolerance: Convergence criterion in kcal/mol/A (default 0.1).
///
/// Returns dict with:
///     coords: Nx3 optimized coordinates
///     initial_energy: energy before minimization
///     final_energy: energy after minimization
///     energy_components: dict of bond/angle/torsion/vdw/es
///     steps: number of steps taken
///     converged: whether optimization converged
#[pyfunction]
#[pyo3(signature = (pdb, max_steps=500, gradient_tolerance=0.1))]
pub fn minimize_hydrogens<'py>(
    py: Python<'py>,
    pdb: &PyPDB,
    max_steps: usize,
    gradient_tolerance: f64,
) -> PyResult<PyObject> {
    let amber = params::amber96();
    let topo = topology::build_topology(&pdb.inner, &amber);
    let coords: Vec<[f64; 3]> = topo.atoms.iter().map(|a| a.pos).collect();

    let result = py.allow_threads(|| {
        minimize::minimize_hydrogens(&coords, &topo, &amber, max_steps, gradient_tolerance)
    });

    let n = result.coords.len();
    let flat: Vec<f64> = result
        .coords
        .iter()
        .flat_map(|c| c.iter().copied())
        .collect();

    let dict = pyo3::types::PyDict::new(py);
    let coords_arr = PyArray1::from_vec(py, flat)
        .reshape([n, 3])
        .expect("reshape");
    dict.set_item("coords", coords_arr)?;
    dict.set_item("initial_energy", result.initial_energy)?;
    dict.set_item("final_energy", result.energy.total)?;
    dict.set_item("steps", result.steps)?;
    dict.set_item("converged", result.converged)?;

    let components = pyo3::types::PyDict::new(py);
    components.set_item("bond_stretch", result.energy.bond_stretch)?;
    components.set_item("angle_bend", result.energy.angle_bend)?;
    components.set_item("torsion", result.energy.torsion)?;
    components.set_item("vdw", result.energy.vdw)?;
    components.set_item("electrostatic", result.energy.electrostatic)?;
    dict.set_item("energy_components", components)?;

    Ok(dict.into_any().unbind())
}

/// Full structure energy minimization using AMBER force field.
///
/// Args:
///     pdb: Structure to minimize.
///     max_steps: Maximum optimization steps (default 1000).
///     gradient_tolerance: Convergence criterion in kcal/mol/A (default 0.1).
///
/// Returns dict with same format as minimize_hydrogens.
#[pyfunction]
#[pyo3(signature = (pdb, max_steps=1000, gradient_tolerance=0.1))]
pub fn minimize_structure<'py>(
    py: Python<'py>,
    pdb: &PyPDB,
    max_steps: usize,
    gradient_tolerance: f64,
) -> PyResult<PyObject> {
    let amber = params::amber96();
    let topo = topology::build_topology(&pdb.inner, &amber);
    let coords: Vec<[f64; 3]> = topo.atoms.iter().map(|a| a.pos).collect();

    let constrained = vec![false; coords.len()]; // nothing constrained

    let result = py.allow_threads(|| {
        minimize::steepest_descent(
            &coords,
            &topo,
            &amber,
            max_steps,
            gradient_tolerance,
            &constrained,
        )
    });

    let n = result.coords.len();
    let flat: Vec<f64> = result
        .coords
        .iter()
        .flat_map(|c| c.iter().copied())
        .collect();

    let dict = pyo3::types::PyDict::new(py);
    let coords_arr = PyArray1::from_vec(py, flat)
        .reshape([n, 3])
        .expect("reshape");
    dict.set_item("coords", coords_arr)?;
    dict.set_item("initial_energy", result.initial_energy)?;
    dict.set_item("final_energy", result.energy.total)?;
    dict.set_item("steps", result.steps)?;
    dict.set_item("converged", result.converged)?;

    let components = pyo3::types::PyDict::new(py);
    components.set_item("bond_stretch", result.energy.bond_stretch)?;
    components.set_item("angle_bend", result.energy.angle_bend)?;
    components.set_item("torsion", result.energy.torsion)?;
    components.set_item("vdw", result.energy.vdw)?;
    components.set_item("electrostatic", result.energy.electrostatic)?;
    dict.set_item("energy_components", components)?;

    Ok(dict.into_any().unbind())
}

// ---------------------------------------------------------------------------
// Module registration
// ---------------------------------------------------------------------------

#[pymodule]
pub fn py_forcefield(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_energy, m)?)?;
    m.add_function(wrap_pyfunction!(minimize_hydrogens, m)?)?;
    m.add_function(wrap_pyfunction!(minimize_structure, m)?)?;
    Ok(())
}
