//! Bond topology inference and force field atom assignment.
//!
//! Infers bonds from interatomic distances, then builds the lists
//! of bond/angle/torsion interactions needed for energy computation.

use std::collections::HashSet;

use super::params::AmberParams;

/// Per-atom data for force field computation.
#[derive(Clone, Debug)]
pub struct FFAtom {
    pub pos: [f64; 3],
    pub amber_type: String,
    pub charge: f64,
    #[allow(dead_code)]
    pub residue_name: String,
    #[allow(dead_code)]
    pub atom_name: String,
    pub element: String,
    pub residue_idx: usize,
    pub is_hydrogen: bool,
}

/// A bond between two atoms.
#[derive(Clone, Debug)]
pub struct Bond {
    pub i: usize,
    pub j: usize,
}

/// An angle between three atoms (i-j-k, j is central).
#[derive(Clone, Debug)]
pub struct Angle {
    pub i: usize,
    pub j: usize,
    pub k: usize,
}

/// A torsion/dihedral between four atoms (i-j-k-l).
#[derive(Clone, Debug)]
pub struct Torsion {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
}

/// Complete topology for energy computation.
#[derive(Clone, Debug)]
pub struct Topology {
    pub atoms: Vec<FFAtom>,
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub torsions: Vec<Torsion>,
    /// Improper torsions for planarity (aromatic rings, peptide bonds)
    pub improper_torsions: Vec<Torsion>,
    /// Pairs that are 1-2 or 1-3 connected (excluded from nonbonded)
    pub excluded_pairs: HashSet<(usize, usize)>,
    /// 1-4 pairs (scaled nonbonded)
    pub pairs_14: HashSet<(usize, usize)>,
    /// Atoms that could not be assigned a force field type (residue:atom names)
    pub unassigned_atoms: Vec<String>,
}

/// Maximum bond distance for element pairs (Å).
fn max_bond_distance(elem_a: &str, elem_b: &str) -> f64 {
    match (elem_a, elem_b) {
        ("H", _) | (_, "H") => 1.3,
        ("S", _) | (_, "S") => 2.2,
        _ => 1.9, // C-C, C-N, C-O, N-O, etc.
    }
}

/// Build topology from a PDB structure.
pub fn build_topology(pdb: &pdbtbx::PDB, params: &AmberParams) -> Topology {
    // Step 1: Extract atoms with type assignment
    let mut atoms = Vec::new();
    let mut unassigned_atoms = Vec::new();
    let mut res_idx = 0;

    for chain in pdb.chains() {
        for residue in chain.residues() {
            let res_name = residue.name().unwrap_or("UNK").to_string();
            for atom in residue.atoms() {
                let (x, y, z) = atom.pos();
                let atom_name = atom.name().trim().to_string();
                let element = atom
                    .element()
                    .map(|e| e.symbol().to_string())
                    .unwrap_or_else(|| "C".to_string());
                let is_hydrogen = element == "H" || element == "D";

                // Look up AMBER type and charge
                let (amber_type, charge) = match params.get_atom_type(&res_name, &atom_name) {
                    Some(e) => (e.amber_type.clone(), e.charge),
                    None => {
                        unassigned_atoms.push(format!("{}:{}", res_name, atom_name));
                        let t = match element.as_str() {
                            "C" => "CT",
                            "N" => "N",
                            "O" => "O",
                            "S" => "S",
                            "H" => "H",
                            "P" => "P",
                            _ => "CT",
                        };
                        (t.to_string(), 0.0)
                    }
                };

                atoms.push(FFAtom {
                    pos: [x, y, z],
                    amber_type,
                    charge,
                    residue_name: res_name.clone(),
                    atom_name,
                    element,
                    residue_idx: res_idx,
                    is_hydrogen,
                });
            }
            res_idx += 1;
        }
    }

    // Step 2: Infer bonds from distances
    let n = atoms.len();
    let mut bonds = Vec::new();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];

    for i in 0..n {
        for j in (i + 1)..n {
            // Only bond within same or adjacent residues
            let res_diff = (atoms[i].residue_idx as isize - atoms[j].residue_idx as isize).unsigned_abs();
            if res_diff > 1 {
                continue;
            }

            let dx = atoms[i].pos[0] - atoms[j].pos[0];
            let dy = atoms[i].pos[1] - atoms[j].pos[1];
            let dz = atoms[i].pos[2] - atoms[j].pos[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            let max_dist = max_bond_distance(&atoms[i].element, &atoms[j].element);
            if dist < max_dist && dist > 0.4 {
                bonds.push(Bond { i, j });
                neighbors[i].push(j);
                neighbors[j].push(i);
            }
        }
    }

    // Step 3: Build angles from bonds
    let mut angles = Vec::new();
    for j in 0..n {
        let nbrs = &neighbors[j];
        for ni in 0..nbrs.len() {
            for nk in (ni + 1)..nbrs.len() {
                angles.push(Angle {
                    i: nbrs[ni],
                    j,
                    k: nbrs[nk],
                });
            }
        }
    }

    // Step 4: Build torsions from angles
    let mut torsions = Vec::new();
    let mut excluded_pairs = HashSet::new();
    let mut pairs_14 = HashSet::new();

    // Record 1-2 and 1-3 exclusions
    for bond in &bonds {
        excluded_pairs.insert((bond.i.min(bond.j), bond.i.max(bond.j)));
    }
    for angle in &angles {
        excluded_pairs.insert((angle.i.min(angle.k), angle.i.max(angle.k)));
    }

    // Build torsions: for each bond j-k, find i bonded to j and l bonded to k
    for bond in &bonds {
        let j = bond.i;
        let k = bond.j;
        for &i in &neighbors[j] {
            if i == k {
                continue;
            }
            for &l in &neighbors[k] {
                if l == j || l == i {
                    continue;
                }
                torsions.push(Torsion { i, j, k, l });
                let pair = (i.min(l), i.max(l));
                if !excluded_pairs.contains(&pair) {
                    pairs_14.insert(pair);
                }
            }
        }
    }

    // Step 5: Build improper torsions
    // For each atom with exactly 3 bonds that is listed in ResidueImproperTorsions,
    // generate improper torsion entries (central atom is position 3 in the 4-tuple).
    let mut improper_torsions = Vec::new();
    for j in 0..n {
        let nbrs = &neighbors[j];
        if nbrs.len() != 3 {
            continue;
        }
        if !params.is_improper_center(&atoms[j].residue_name, &atoms[j].atom_name) {
            continue;
        }
        // Generate all permutations of the 3 neighbors as (a1, a2, center=j, a4)
        let (n0, n1, n2) = (nbrs[0], nbrs[1], nbrs[2]);
        for &(i, k, l) in &[(n0, n1, n2), (n0, n2, n1), (n1, n0, n2),
                             (n1, n2, n0), (n2, n0, n1), (n2, n1, n0)] {
            let ti = &atoms[i].amber_type;
            let tk = &atoms[j].amber_type; // central atom
            let tj = &atoms[k].amber_type;
            let tl = &atoms[l].amber_type;
            if params.get_improper_torsion(ti, tj, tk, tl).is_some() {
                improper_torsions.push(Torsion { i, j: k, k: j, l });
                break; // only one improper per center atom per neighbor set
            }
        }
    }

    Topology {
        atoms,
        bonds,
        angles,
        torsions,
        improper_torsions,
        excluded_pairs,
        pairs_14,
        unassigned_atoms,
    }
}
