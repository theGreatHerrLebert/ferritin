//! Multi-chain (MM) alignment algorithm.
//!
//! Ported from C++ USAlign `MMalign.h`.
//! Aligns two protein/RNA complexes by determining optimal chain-to-chain
//! correspondence, superposing concatenated chains, and iteratively refining.

pub mod chain_assign;
pub mod complex_score;
pub mod dimer;
pub mod iter;
pub mod trim;

use crate::core::types::{Coord3D, MolType, Transform};

use crate::ext::se::SeResult;

/// Per-chain structure data for multi-chain alignment.
#[derive(Debug, Clone)]
pub struct ChainData {
    /// CA (or C3') coordinates.
    pub coords: Vec<Coord3D>,
    /// One-letter amino acid or nucleotide codes.
    pub sequence: Vec<u8>,
    /// Secondary structure assignment.
    pub sec_structure: Vec<u8>,
    /// Chain identifier.
    pub chain_id: String,
    /// Molecule type (protein vs RNA).
    pub mol_type: MolType,
}

impl ChainData {
    /// Number of residues in this chain.
    pub fn len(&self) -> usize {
        self.coords.len()
    }

    /// Whether this chain has no residues.
    pub fn is_empty(&self) -> bool {
        self.coords.is_empty()
    }
}

/// Result of multi-chain alignment.
#[derive(Debug, Clone)]
pub struct MMAlignResult {
    /// Total TM-score of the complex alignment.
    pub total_score: f64,
    /// Chain assignment pairs: `(chain_i_in_complex1, chain_j_in_complex2)`.
    pub chain_assignments: Vec<(usize, usize)>,
    /// Per-chain SE alignment results (one per assigned chain pair).
    pub per_chain_results: Vec<SeResult>,
    /// Per-chain transforms (one per assigned chain pair).
    pub transforms: Vec<Transform>,
}

/// Count the number of nucleic acid and protein chains.
///
/// Returns `(na_count, aa_count)`.
pub fn count_na_aa_chains(mol_vec: &[MolType]) -> (usize, usize) {
    let mut na = 0usize;
    let mut aa = 0usize;
    for m in mol_vec {
        match m {
            MolType::RNA => na += 1,
            MolType::Protein => aa += 1,
        }
    }
    (na, aa)
}

/// Compute total chain lengths split by molecule type.
///
/// Returns `(len_aa, len_na)`.
pub fn total_chain_lengths(chains: &[ChainData]) -> (usize, usize) {
    let mut len_aa = 0usize;
    let mut len_na = 0usize;
    for c in chains {
        match c.mol_type {
            MolType::Protein => len_aa += c.len(),
            MolType::RNA => len_na += c.len(),
        }
    }
    (len_aa, len_na)
}

/// Concatenate assigned chain pairs into single coordinate/sequence arrays.
///
/// Corresponds to C++ `copy_chain_pair_data`.
/// Returns `(xa, ya, seqx, seqy, secx, secy, mol_type_sum)`.
pub fn concatenate_assigned_chains(
    x_chains: &[ChainData],
    y_chains: &[ChainData],
    assign1: &[i32],
    seqx_a_mat: &[Vec<String>],
    seqy_a_mat: &[Vec<String>],
) -> (Vec<Coord3D>, Vec<Coord3D>, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, i32, Vec<String>) {
    let mut xa = Vec::new();
    let mut ya = Vec::new();
    let mut seqx = Vec::new();
    let mut seqy = Vec::new();
    let mut secx = Vec::new();
    let mut secy = Vec::new();
    let mut mol_type_sum = 0i32;
    let mut sequence = vec![String::new(), String::new()];

    for (i, &j) in assign1.iter().enumerate() {
        if j < 0 {
            continue;
        }
        let ju = j as usize;
        let xc = &x_chains[i];
        let yc = &y_chains[ju];

        xa.extend_from_slice(&xc.coords);
        seqx.extend_from_slice(&xc.sequence);
        secx.extend_from_slice(&xc.sec_structure);

        ya.extend_from_slice(&yc.coords);
        seqy.extend_from_slice(&yc.sequence);
        secy.extend_from_slice(&yc.sec_structure);

        let mol_i = match xc.mol_type {
            MolType::Protein => -1,
            MolType::RNA => 1,
        };
        let mol_j = match yc.mol_type {
            MolType::Protein => -1,
            MolType::RNA => 1,
        };
        mol_type_sum += mol_i + mol_j;

        if i < seqx_a_mat.len() && ju < seqx_a_mat[i].len() {
            sequence[0].push_str(&seqx_a_mat[i][ju]);
        }
        if i < seqy_a_mat.len() && ju < seqy_a_mat[i].len() {
            sequence[1].push_str(&seqy_a_mat[i][ju]);
        }
    }

    (xa, ya, seqx, seqy, secx, secy, mol_type_sum, sequence)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_na_aa_chains() {
        let mol_vec = vec![MolType::Protein, MolType::RNA, MolType::Protein];
        let (na, aa) = count_na_aa_chains(&mol_vec);
        assert_eq!(na, 1);
        assert_eq!(aa, 2);
    }

    #[test]
    fn test_total_chain_lengths() {
        let chains = vec![
            ChainData {
                coords: vec![[0.0; 3]; 10],
                sequence: vec![b'A'; 10],
                sec_structure: vec![b'C'; 10],
                chain_id: "A".to_string(),
                mol_type: MolType::Protein,
            },
            ChainData {
                coords: vec![[0.0; 3]; 5],
                sequence: vec![b'a'; 5],
                sec_structure: vec![b'C'; 5],
                chain_id: "B".to_string(),
                mol_type: MolType::RNA,
            },
        ];
        let (len_aa, len_na) = total_chain_lengths(&chains);
        assert_eq!(len_aa, 10);
        assert_eq!(len_na, 5);
    }

    #[test]
    fn test_chain_data_len() {
        let cd = ChainData {
            coords: vec![[1.0, 2.0, 3.0]; 7],
            sequence: vec![b'G'; 7],
            sec_structure: vec![b'H'; 7],
            chain_id: "X".to_string(),
            mol_type: MolType::Protein,
        };
        assert_eq!(cd.len(), 7);
        assert!(!cd.is_empty());
    }
}
