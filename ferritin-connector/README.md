# ferritin-connector

PyO3 bindings layer for the [ferritin](https://github.com/theGreatHerrLebert/ferritin)
structural bioinformatics toolkit.

This package is the compiled Rust ↔ Python bridge. Most users should install the
`ferritin` Pythonic wrapper instead:

```bash
pip install ferritin
```

`ferritin-connector` is pulled in automatically as a dependency and exposes the raw
`#[pyclass]` types the wrapper composes on top of. Direct use is supported but not
the recommended entry point — see the main repository for documentation, examples,
and the public API.

## Project

- Source: https://github.com/theGreatHerrLebert/ferritin
- License: MIT
