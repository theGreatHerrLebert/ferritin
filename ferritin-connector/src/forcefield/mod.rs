//! AMBER force field for energy minimization.
//!
//! Implements the four standard AMBER energy terms:
//! - Bond stretching: E = k(r - r0)²
//! - Angle bending: E = k(θ - θ0)²
//! - Torsion: E = (V/div)(1 + cos(f*φ - φ0))
//! - Nonbonded: Lennard-Jones 12-6 + Coulomb
//!
//! Parameters from AMBER96 (parm96.dat).

pub mod params;
pub mod topology;
pub mod energy;
pub mod minimize;
