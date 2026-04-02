//! Energy minimization algorithms.
//!
//! Steepest descent with adaptive step size. Simple but robust
//! for hydrogen position optimization.

use super::energy::{compute_energy, compute_energy_and_forces, EnergyResult};
use super::params::AmberParams;
use super::topology::Topology;

/// Result of energy minimization.
#[derive(Clone, Debug)]
pub struct MinimizeResult {
    /// Optimized coordinates
    pub coords: Vec<[f64; 3]>,
    /// Final energy breakdown
    pub energy: EnergyResult,
    /// Initial total energy
    pub initial_energy: f64,
    /// Number of steps taken
    pub steps: usize,
    /// Whether minimization converged
    pub converged: bool,
}

/// Minimize energy using steepest descent with line search.
///
/// # Arguments
/// * `coords` — Initial coordinates (modified in place)
/// * `topo` — Topology (bonds, angles, etc.)
/// * `params` — Force field parameters
/// * `max_steps` — Maximum iterations
/// * `gradient_tolerance` — Convergence criterion (kcal/mol/Å)
/// * `constrained` — Indices of atoms that should not move
pub fn steepest_descent(
    coords: &[[f64; 3]],
    topo: &Topology,
    params: &AmberParams,
    max_steps: usize,
    gradient_tolerance: f64,
    constrained: &[bool],
) -> MinimizeResult {
    let n = coords.len();
    let mut pos: Vec<[f64; 3]> = coords.to_vec();

    let initial_e = compute_energy(&pos, topo, params);
    let initial_energy = initial_e.total;

    let mut step_size = 0.01; // initial step size in Å
    let mut prev_energy = initial_energy;
    let mut converged = false;
    let mut steps = 0;

    for step in 0..max_steps {
        steps = step + 1;

        let (energy, forces) = compute_energy_and_forces(&pos, topo, params);

        // Compute max force magnitude
        let mut max_force = 0.0f64;
        for i in 0..n {
            if constrained[i] {
                continue;
            }
            let f2 = forces[i][0] * forces[i][0]
                + forces[i][1] * forces[i][1]
                + forces[i][2] * forces[i][2];
            max_force = max_force.max(f2.sqrt());
        }

        // Check convergence
        if max_force < gradient_tolerance {
            converged = true;
            break;
        }

        // Take step along gradient direction
        let mut new_pos = pos.clone();
        for i in 0..n {
            if constrained[i] {
                continue;
            }
            let f_mag = (forces[i][0] * forces[i][0]
                + forces[i][1] * forces[i][1]
                + forces[i][2] * forces[i][2])
            .sqrt();
            if f_mag < 1e-12 {
                continue;
            }
            // Normalize force direction, step by step_size
            let scale = step_size / f_mag;
            new_pos[i][0] += forces[i][0] * scale;
            new_pos[i][1] += forces[i][1] * scale;
            new_pos[i][2] += forces[i][2] * scale;
        }

        let new_energy = compute_energy(&new_pos, topo, params);

        // Adaptive step size
        if new_energy.total < prev_energy {
            // Accept step, increase step size
            pos = new_pos;
            prev_energy = new_energy.total;
            step_size *= 1.2;
            step_size = step_size.min(0.1); // cap at 0.1 Å
        } else {
            // Reject step, decrease step size
            step_size *= 0.5;
            if step_size < 1e-8 {
                break; // can't make progress
            }
        }
    }

    let final_energy = compute_energy(&pos, topo, params);

    MinimizeResult {
        coords: pos,
        energy: final_energy,
        initial_energy,
        steps,
        converged,
    }
}

/// Minimize only hydrogen positions (freeze all heavy atoms).
pub fn minimize_hydrogens(
    coords: &[[f64; 3]],
    topo: &Topology,
    params: &AmberParams,
    max_steps: usize,
    gradient_tolerance: f64,
) -> MinimizeResult {
    let constrained: Vec<bool> = topo
        .atoms
        .iter()
        .map(|a| !a.is_hydrogen)
        .collect();

    steepest_descent(coords, topo, params, max_steps, gradient_tolerance, &constrained)
}
