//! Convert between pdbtbx structures and Arrow RecordBatches.
//!
//! Two directions:
//! - `pdb_to_atom_batch` — PDB → Arrow (export)
//! - `atom_batch_to_pdb` — Arrow → PDB (import)

use std::collections::BTreeMap;

use arrow::array::{
    Array, AsArray, BooleanArray, Float64Array, Int64Array, RecordBatch, StringArray, UInt32Array,
};

use crate::atom::AtomBatchBuilder;
use crate::structure::StructureBatchBuilder;

const BACKBONE_ATOMS: &[&str] = &["N", "CA", "C", "O"];

/// Convert a pdbtbx PDB structure into a per-atom Arrow RecordBatch.
///
/// Iterates over all models, chains, residues, and atoms, emitting
/// one row per atom with full metadata.
pub fn pdb_to_atom_batch(pdb: &pdbtbx::PDB, structure_id: &str) -> anyhow::Result<RecordBatch> {
    let n_atoms = pdb.atom_count();
    let mut builder = AtomBatchBuilder::new(n_atoms);

    for (model_idx, model) in pdb.models().enumerate() {
        for chain in model.chains() {
            let chain_id = chain.id();
            for residue in chain.residues() {
                let res_name = residue
                    .name()
                    .unwrap_or("UNK");
                let res_serial = residue.serial_number() as i64;

                for conformer in residue.conformers() {
                    for atom in conformer.atoms() {
                        let atom_name = atom.name();
                        let element = atom.element().map(|e| e.symbol());
                        let pos = atom.pos();
                        let is_backbone = BACKBONE_ATOMS.contains(&atom_name);

                        builder.append(
                            structure_id,
                            model_idx as u32,
                            chain_id,
                            res_name,
                            res_serial,
                            atom_name,
                            atom.serial_number() as i64,
                            element,
                            pos.0,
                            pos.1,
                            pos.2,
                            atom.b_factor(),
                            atom.occupancy(),
                            atom.hetero(),
                            is_backbone,
                        );
                    }
                }
            }
        }
    }

    builder.finish()
}

/// Convert a pdbtbx PDB structure into a per-structure summary Arrow RecordBatch.
pub fn pdb_to_structure_batch(
    pdb: &pdbtbx::PDB,
    structure_id: &str,
) -> anyhow::Result<RecordBatch> {
    let mut builder = StructureBatchBuilder::new(1);

    let model = pdb.models().next();
    let (chain_count, chains_str) = match model {
        Some(m) => {
            let chains: Vec<&str> = m.chains().map(|c| c.id()).collect();
            (chains.len() as u32, chains.join(","))
        }
        None => (0, String::new()),
    };

    builder.append(
        structure_id,
        pdb.atom_count() as i64,
        pdb.residue_count() as i64,
        chain_count,
        pdb.model_count() as u32,
        &chains_str,
    );

    builder.finish()
}

/// Convert multiple PDB structures into a single atom RecordBatch.
///
/// Each structure is identified by its `structure_id` in the batch.
pub fn pdbs_to_atom_batch(
    pdbs: &[(&pdbtbx::PDB, &str)],
) -> anyhow::Result<RecordBatch> {
    let total_atoms: usize = pdbs.iter().map(|(pdb, _)| pdb.atom_count()).sum();
    let mut builder = AtomBatchBuilder::new(total_atoms);

    for (pdb, structure_id) in pdbs {
        for (model_idx, model) in pdb.models().enumerate() {
            for chain in model.chains() {
                let chain_id = chain.id();
                for residue in chain.residues() {
                    let res_name = residue.name().unwrap_or("UNK");
                    let res_serial = residue.serial_number() as i64;

                    for conformer in residue.conformers() {
                        for atom in conformer.atoms() {
                            let atom_name = atom.name();
                            let element = atom.element().map(|e| e.symbol());
                            let pos = atom.pos();
                            let is_backbone = BACKBONE_ATOMS.contains(&atom_name);

                            builder.append(
                                structure_id,
                                model_idx as u32,
                                chain_id,
                                res_name,
                                res_serial,
                                atom_name,
                                atom.serial_number() as i64,
                                element,
                                pos.0,
                                pos.1,
                                pos.2,
                                atom.b_factor(),
                                atom.occupancy(),
                                atom.hetero(),
                                is_backbone,
                            );
                        }
                    }
                }
            }
        }
    }

    builder.finish()
}

// ============================================================================
// Arrow → PDB (from_arrow)
// ============================================================================

/// Helper to downcast a RecordBatch column by index.
fn col_str(batch: &RecordBatch, idx: usize) -> &StringArray {
    batch.column(idx).as_any().downcast_ref::<StringArray>().unwrap()
}
fn col_f64(batch: &RecordBatch, idx: usize) -> &Float64Array {
    batch.column(idx).as_any().downcast_ref::<Float64Array>().unwrap()
}
fn col_i64(batch: &RecordBatch, idx: usize) -> &Int64Array {
    batch.column(idx).as_any().downcast_ref::<Int64Array>().unwrap()
}
fn col_u32(batch: &RecordBatch, idx: usize) -> &UInt32Array {
    batch.column(idx).as_any().downcast_ref::<UInt32Array>().unwrap()
}
fn col_bool(batch: &RecordBatch, idx: usize) -> &BooleanArray {
    batch.column(idx).as_any().downcast_ref::<BooleanArray>().unwrap()
}

/// Convert an atom-schema Arrow RecordBatch back into pdbtbx PDB structures.
///
/// Groups rows by `structure_id`, then by model/chain/residue to rebuild
/// the hierarchy. Returns a Vec of (structure_id, PDB) pairs.
///
/// This enables round-tripping: PDB → Arrow → filter/transform → PDB.
pub fn atom_batch_to_pdbs(batch: &RecordBatch) -> anyhow::Result<Vec<(String, pdbtbx::PDB)>> {
    let n = batch.num_rows();
    if n == 0 {
        return Ok(vec![]);
    }

    let structure_ids = col_str(batch, 0);
    let models = col_u32(batch, 1);
    let chain_ids = col_str(batch, 2);
    let residue_names = col_str(batch, 3);
    let residue_serials = col_i64(batch, 4);
    let atom_names = col_str(batch, 5);
    let atom_serials = col_i64(batch, 6);
    let elements = batch.column(7).as_string::<i32>();
    let xs = col_f64(batch, 8);
    let ys = col_f64(batch, 9);
    let zs = col_f64(batch, 10);
    let b_factors = col_f64(batch, 11);
    let occupancies = col_f64(batch, 12);
    let is_hetero = col_bool(batch, 13);

    // Group rows by structure_id, preserving order
    let mut structure_order: Vec<String> = Vec::new();
    let mut structure_rows: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let sid = structure_ids.value(i).to_string();
        structure_rows.entry(sid.clone()).or_insert_with(|| {
            structure_order.push(sid.clone());
            Vec::new()
        }).push(i);
    }

    let mut results = Vec::new();

    for sid in &structure_order {
        let rows = &structure_rows[sid];
        let mut pdb = pdbtbx::PDB::new();

        // Group by model
        let mut model_rows: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
        for &row in rows {
            model_rows.entry(models.value(row)).or_default().push(row);
        }

        for (_model_idx, m_rows) in &model_rows {
            let mut model = pdbtbx::Model::new(0);

            // Group by chain
            let mut chain_rows: BTreeMap<String, Vec<usize>> = BTreeMap::new();
            let mut chain_order: Vec<String> = Vec::new();
            for &row in m_rows {
                let cid = chain_ids.value(row).to_string();
                chain_rows.entry(cid.clone()).or_insert_with(|| {
                    chain_order.push(cid.clone());
                    Vec::new()
                }).push(row);
            }

            for cid in &chain_order {
                let c_rows = &chain_rows[cid];
                let mut chain = pdbtbx::Chain::new(cid).unwrap();

                // Group by residue (serial number)
                let mut res_rows: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
                let mut res_order: Vec<i64> = Vec::new();
                for &row in c_rows {
                    let rs = residue_serials.value(row);
                    res_rows.entry(rs).or_insert_with(|| {
                        res_order.push(rs);
                        Vec::new()
                    }).push(row);
                }

                for rs in &res_order {
                    let r_rows = &res_rows[rs];
                    let first = r_rows[0];
                    let res_name = residue_names.value(first);

                    let mut residue = pdbtbx::Residue::new(*rs as isize, None, None)
                        .unwrap();
                    let mut conformer = pdbtbx::Conformer::new(res_name, None, None)
                        .unwrap();

                    for &row in r_rows {
                        let element_str = if elements.is_null(row) {
                            "C"
                        } else {
                            elements.value(row)
                        };

                        let atom = pdbtbx::Atom::new(
                            is_hetero.value(row),
                            atom_serials.value(row) as usize,
                            "",  // conformer/alt-loc ID
                            atom_names.value(row),
                            xs.value(row),
                            ys.value(row),
                            zs.value(row),
                            occupancies.value(row),
                            b_factors.value(row),
                            element_str,
                            0,  // charge
                        ).unwrap();

                        conformer.add_atom(atom);
                    }

                    residue.add_conformer(conformer);
                    chain.add_residue(residue);
                }

                model.add_chain(chain);
            }

            pdb.add_model(model);
        }

        results.push((sid.clone(), pdb));
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load_test_pdb() -> Option<pdbtbx::PDB> {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../test-pdbs/1crn.pdb");
        if !std::path::Path::new(path).exists() {
            return None;
        }
        let (pdb, _) = pdbtbx::ReadOptions::new()
            .set_level(pdbtbx::StrictnessLevel::Loose)
            .read(path)
            .ok()?;
        Some(pdb)
    }

    #[test]
    fn test_pdb_to_atom_batch() {
        let pdb = match load_test_pdb() {
            Some(p) => p,
            None => return,
        };

        let batch = pdb_to_atom_batch(&pdb, "1crn").unwrap();
        assert!(batch.num_rows() > 300, "crambin should have >300 atoms");
        assert_eq!(batch.num_columns(), 15);

        let ids = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(ids.value(0), "1crn");
    }

    #[test]
    fn test_pdb_to_structure_batch() {
        let pdb = match load_test_pdb() {
            Some(p) => p,
            None => return,
        };

        let batch = pdb_to_structure_batch(&pdb, "1crn").unwrap();
        assert_eq!(batch.num_rows(), 1);
    }

    #[test]
    fn test_roundtrip() {
        let pdb = match load_test_pdb() {
            Some(p) => p,
            None => return,
        };

        let original_atoms = pdb.atom_count();

        // PDB → Arrow
        let batch = pdb_to_atom_batch(&pdb, "1crn").unwrap();
        assert_eq!(batch.num_rows(), original_atoms);

        // Arrow → PDB
        let rebuilt = atom_batch_to_pdbs(&batch).unwrap();
        assert_eq!(rebuilt.len(), 1);

        let (sid, rebuilt_pdb) = &rebuilt[0];
        assert_eq!(sid, "1crn");
        assert_eq!(rebuilt_pdb.atom_count(), original_atoms);

        // Spot-check: first atom coordinates should match
        let orig_atom = pdb.models().next().unwrap()
            .chains().next().unwrap()
            .residues().next().unwrap()
            .conformers().next().unwrap()
            .atoms().next().unwrap();
        let rebuilt_atom = rebuilt_pdb.models().next().unwrap()
            .chains().next().unwrap()
            .residues().next().unwrap()
            .conformers().next().unwrap()
            .atoms().next().unwrap();

        let (ox, oy, oz) = orig_atom.pos();
        let (rx, ry, rz) = rebuilt_atom.pos();
        assert!((ox - rx).abs() < 1e-10, "x mismatch: {ox} vs {rx}");
        assert!((oy - ry).abs() < 1e-10, "y mismatch: {oy} vs {ry}");
        assert!((oz - rz).abs() < 1e-10, "z mismatch: {oz} vs {rz}");
    }
}
