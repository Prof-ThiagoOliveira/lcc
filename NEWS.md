# lcc 4.2.2

- Hardened the bootstrap workflow for `lcc()` and related summaries, ensuring parallel workers reliably load the modelling internals.
- Unified the bootstrap metric outputs so LCC, Pearson correlation, and accuracy share a consistent structure across serial and parallel runs.
- Added documented bootstrap confidence intervals to the user workflow and refreshed package documentation accordingly.

## Compatibility

- No breaking user-facing changes; existing APIs continue to work, with the major version signalling the enhanced bootstrap capabilities and infrastructure improvements.
