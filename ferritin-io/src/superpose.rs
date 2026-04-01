//! Superposition output (RasMol/PyMOL scripts).
//!
//! TODO: Port the full output_superpose from C++ TMalign (lines 3047-3660).
//! This generates PDB files and visualization scripts.

use std::io::Write;

use ferritin_align::core::types::{AlignResult, StructureData};

/// Write superposed structures to files with the given prefix.
///
/// Generates:
/// - `{prefix}` — RasMol script for aligned CA traces
/// - `{prefix}.pml` — PyMOL script for aligned CA traces
pub fn output_superpose<W: Write>(
    w: &mut W,
    _result: &AlignResult,
    _chain1: &StructureData,
    _chain2: &StructureData,
    _prefix: &str,
) -> std::io::Result<()> {
    writeln!(w, "# Superposition output not yet implemented")?;
    writeln!(w, "# Use -m option to output rotation matrix instead")?;
    Ok(())
}
