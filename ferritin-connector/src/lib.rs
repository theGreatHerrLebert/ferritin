use pyo3::prelude::*;
use pyo3::wrap_pymodule;

mod py_align;
mod py_structure;
mod py_transform;

/// ferritin_connector — PyO3 bindings for the ferritin structural bioinformatics toolkit.
#[pymodule]
fn ferritin_connector(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(py_align::py_align))?;
    m.add_wrapped(wrap_pymodule!(py_structure::py_structure))?;
    m.add_wrapped(wrap_pymodule!(py_transform::py_transform))?;
    Ok(())
}
