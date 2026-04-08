//! AMBER96 force field parameters.
//!
//! Atom types, charges, bond/angle/torsion/LJ parameters parsed from amber96.ini.

use std::collections::{HashMap, HashSet};

/// Atom type assignment: residue_name:atom_name → (amber_type, charge)
#[derive(Clone, Debug)]
pub struct AtomTypeEntry {
    pub amber_type: String,
    pub charge: f64, // in elementary charge units
}

/// Bond stretch parameters: E = k * (r - r0)²
#[derive(Clone, Debug)]
pub struct BondParam {
    pub k: f64,  // kcal/mol/Å²
    pub r0: f64, // Å
}

/// Angle bend parameters: E = k * (θ - θ0)²
#[derive(Clone, Debug)]
pub struct AngleParam {
    pub k: f64,     // kcal/mol/rad²
    pub theta0: f64, // radians
}

/// Single torsion term: E = (V/div) * (1 + cos(f*φ - φ0))
#[derive(Clone, Debug)]
pub struct TorsionTerm {
    pub v: f64,    // barrier height (kcal/mol)
    pub phi0: f64, // phase (radians)
    pub f: f64,    // periodicity
    pub div: f64,  // divisor
}

/// Lennard-Jones parameters per atom type
#[derive(Clone, Debug)]
pub struct LJParam {
    pub r: f64,       // van der Waals radius (Å)
    pub epsilon: f64,  // well depth (kcal/mol)
}

/// Complete AMBER force field parameter set.
#[derive(Clone, Debug)]
pub struct AmberParams {
    /// Residue:Atom → (type, charge)
    pub atom_types: HashMap<String, AtomTypeEntry>,
    /// Wildcard (*:Atom) entries
    pub wildcard_types: HashMap<String, AtomTypeEntry>,
    /// (type_i, type_j) → BondParam (sorted key)
    pub bonds: HashMap<(String, String), BondParam>,
    /// (type_i, type_j, type_k) → AngleParam (sorted outer key)
    pub angles: HashMap<(String, String, String), AngleParam>,
    /// (type_i, type_j, type_k, type_l) → Vec<TorsionTerm>
    pub torsions: HashMap<(String, String, String, String), Vec<TorsionTerm>>,
    /// Improper torsion parameters (separate from proper)
    pub improper_torsions: HashMap<(String, String, String, String), Vec<TorsionTerm>>,
    /// Atoms that should have improper torsions (e.g., "ALA:N", "ALA:C")
    pub residue_impropers: HashSet<String>,
    /// type → LJParam
    pub lj: HashMap<String, LJParam>,
    /// 1-4 electrostatic scaling factor
    pub scee: f64,
    /// 1-4 vdW scaling factor
    pub scnb: f64,
}

fn sorted_pair(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}

fn sorted_triple(a: &str, b: &str, c: &str) -> (String, String, String) {
    if a <= c {
        (a.to_string(), b.to_string(), c.to_string())
    } else {
        (c.to_string(), b.to_string(), a.to_string())
    }
}

impl AmberParams {
    /// Parse AMBER parameters from INI file content.
    pub fn from_ini(content: &str) -> Self {
        let mut params = AmberParams {
            atom_types: HashMap::new(),
            wildcard_types: HashMap::new(),
            bonds: HashMap::new(),
            angles: HashMap::new(),
            torsions: HashMap::new(),
            improper_torsions: HashMap::new(),
            residue_impropers: HashSet::new(),
            lj: HashMap::new(),
            scee: 1.2,
            scnb: 2.0,
        };

        let mut section = String::new();

        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with(';') {
                continue;
            }

            // Section header
            if line.starts_with('[') {
                section = line
                    .trim_start_matches('[')
                    .trim_end_matches(']')
                    .to_string();
                continue;
            }

            // Properties
            if line.starts_with('@') {
                if line.starts_with("@SCEE=") {
                    if let Ok(v) = line[6..].parse::<f64>() {
                        params.scee = v;
                    }
                }
                continue;
            }

            // Header lines
            if line.starts_with("ver:") || line.starts_with("key:") {
                continue;
            }

            // Parse data lines by section
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }

            match section.as_str() {
                "QuadraticBondStretch" => {
                    // ver I J k r0 comment
                    if fields.len() >= 5 {
                        let key = sorted_pair(fields[1], fields[2]);
                        if let (Ok(k), Ok(r0)) =
                            (fields[3].parse::<f64>(), fields[4].parse::<f64>())
                        {
                            params.bonds.insert(key, BondParam { k, r0 });
                        }
                    }
                }
                "QuadraticAngleBend" => {
                    // ver I J K k theta0 comment
                    if fields.len() >= 6 {
                        let key = sorted_triple(fields[1], fields[2], fields[3]);
                        if let (Ok(k), Ok(theta0_deg)) =
                            (fields[4].parse::<f64>(), fields[5].parse::<f64>())
                        {
                            params.angles.insert(
                                key,
                                AngleParam {
                                    k,
                                    theta0: theta0_deg.to_radians(),
                                },
                            );
                        }
                    }
                }
                "Torsions" => {
                    // ver I J K L N div V phi0 f comment
                    if fields.len() >= 10 {
                        let key = (
                            fields[1].to_string(),
                            fields[2].to_string(),
                            fields[3].to_string(),
                            fields[4].to_string(),
                        );
                        if let (Ok(div), Ok(v), Ok(phi0_deg), Ok(f)) = (
                            fields[6].parse::<f64>(),
                            fields[7].parse::<f64>(),
                            fields[8].parse::<f64>(),
                            fields[9].parse::<f64>(),
                        ) {
                            let term = TorsionTerm {
                                v,
                                phi0: phi0_deg.to_radians(),
                                f,
                                div: div.max(1.0),
                            };
                            params
                                .torsions
                                .entry(key)
                                .or_insert_with(Vec::new)
                                .push(term);
                        }
                    }
                }
                "ImproperTorsions" => {
                    // Same format as Torsions but stored separately
                    if fields.len() >= 10 {
                        let key = (
                            fields[1].to_string(),
                            fields[2].to_string(),
                            fields[3].to_string(),
                            fields[4].to_string(),
                        );
                        if let (Ok(div), Ok(v), Ok(phi0_deg), Ok(f)) = (
                            fields[6].parse::<f64>(),
                            fields[7].parse::<f64>(),
                            fields[8].parse::<f64>(),
                            fields[9].parse::<f64>(),
                        ) {
                            let term = TorsionTerm {
                                v,
                                phi0: phi0_deg.to_radians(),
                                f,
                                div: div.max(1.0),
                            };
                            params
                                .improper_torsions
                                .entry(key)
                                .or_insert_with(Vec::new)
                                .push(term);
                        }
                    }
                }
                "ResidueImproperTorsions" => {
                    // Single column: residue:atom names
                    if fields.len() >= 1 {
                        let name = fields[0].trim().to_string();
                        if name.contains(':') {
                            params.residue_impropers.insert(name);
                        }
                    }
                }
                "LennardJones" => {
                    // ver I R epsilon comment
                    if fields.len() >= 4 {
                        if let (Ok(r), Ok(eps)) =
                            (fields[2].parse::<f64>(), fields[3].parse::<f64>())
                        {
                            params.lj.insert(
                                fields[1].to_string(),
                                LJParam { r, epsilon: eps },
                            );
                        }
                    }
                }
                "ChargesAndTypeNames" => {
                    // ver name q type
                    if fields.len() >= 4 {
                        let name = fields[1]; // e.g., "ALA:N"
                        if let Ok(q) = fields[2].parse::<f64>() {
                            let atype = fields[3].to_string();
                            let entry = AtomTypeEntry {
                                amber_type: atype,
                                charge: q,
                            };
                            if name.starts_with("*:") {
                                let atom_name = name[2..].to_string();
                                params.wildcard_types.insert(atom_name, entry);
                            } else {
                                params.atom_types.insert(name.to_string(), entry);
                            }
                        }
                    }
                }
                _ => {}
            }
        }

        params
    }

    /// Look up atom type and charge for a given residue:atom pair.
    pub fn get_atom_type(&self, residue: &str, atom: &str) -> Option<&AtomTypeEntry> {
        let key = format!("{residue}:{atom}");
        self.atom_types
            .get(&key)
            .or_else(|| self.wildcard_types.get(atom))
    }

    /// Look up bond parameters for two atom types.
    pub fn get_bond(&self, type_a: &str, type_b: &str) -> Option<&BondParam> {
        let key = sorted_pair(type_a, type_b);
        self.bonds.get(&key)
    }

    /// Look up angle parameters.
    pub fn get_angle(&self, type_a: &str, type_b: &str, type_c: &str) -> Option<&AngleParam> {
        let key = sorted_triple(type_a, type_b, type_c);
        self.angles.get(&key)
    }

    /// Look up torsion parameters (tries specific then wildcard).
    pub fn get_torsion(
        &self,
        type_a: &str,
        type_b: &str,
        type_c: &str,
        type_d: &str,
    ) -> Option<&Vec<TorsionTerm>> {
        // Try exact match first
        let key = (
            type_a.to_string(),
            type_b.to_string(),
            type_c.to_string(),
            type_d.to_string(),
        );
        if let Some(terms) = self.torsions.get(&key) {
            return Some(terms);
        }
        // Try wildcard on outer atoms
        let key = (
            "*".to_string(),
            type_b.to_string(),
            type_c.to_string(),
            "*".to_string(),
        );
        self.torsions.get(&key)
    }

    /// Look up improper torsion parameters (tries specific then wildcard).
    pub fn get_improper_torsion(
        &self,
        type_a: &str,
        type_b: &str,
        type_c: &str,
        type_d: &str,
    ) -> Option<&Vec<TorsionTerm>> {
        let key = (
            type_a.to_string(),
            type_b.to_string(),
            type_c.to_string(),
            type_d.to_string(),
        );
        if let Some(terms) = self.improper_torsions.get(&key) {
            return Some(terms);
        }
        // Wildcard on outer atoms
        let key = (
            "*".to_string(),
            type_b.to_string(),
            type_c.to_string(),
            type_d.to_string(),
        );
        if let Some(terms) = self.improper_torsions.get(&key) {
            return Some(terms);
        }
        let key = (
            "*".to_string(),
            "*".to_string(),
            type_c.to_string(),
            type_d.to_string(),
        );
        self.improper_torsions.get(&key)
    }

    /// Check if an atom should have improper torsions.
    pub fn is_improper_center(&self, residue: &str, atom: &str) -> bool {
        let key = format!("{residue}:{atom}");
        self.residue_impropers.contains(&key)
    }

    /// Look up LJ parameters for an atom type.
    pub fn get_lj(&self, atype: &str) -> Option<&LJParam> {
        self.lj.get(atype)
    }
}

/// Load the embedded AMBER96 parameter set.
pub fn amber96() -> AmberParams {
    AmberParams::from_ini(include_str!("../../data/amber96.ini"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_amber96() {
        let p = amber96();
        assert!(!p.bonds.is_empty());
        assert!(!p.angles.is_empty());
        assert!(!p.torsions.is_empty());
        assert!(!p.lj.is_empty());
        assert!(!p.atom_types.is_empty());
        assert!((p.scee - 1.2).abs() < 1e-6);
    }

    #[test]
    fn test_atom_type_lookup() {
        let p = amber96();
        let entry = p.get_atom_type("ALA", "CA").unwrap();
        assert_eq!(entry.amber_type, "CT");
        assert!((entry.charge - 0.0337).abs() < 0.001);
    }

    #[test]
    fn test_bond_lookup() {
        let p = amber96();
        let b = p.get_bond("CT", "CT").unwrap();
        assert!((b.k - 310.0).abs() < 1.0);
        assert!((b.r0 - 1.526).abs() < 0.01);
    }

    #[test]
    fn test_lj_lookup() {
        let p = amber96();
        let lj = p.get_lj("CT").unwrap();
        assert!((lj.r - 1.908).abs() < 0.01);
        assert!((lj.epsilon - 0.1094).abs() < 0.01);
    }
}
