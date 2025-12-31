# lcc 4.2.2

- Hardened the bootstrap workflow for `lcc()` and related summaries, ensuring parallel workers reliably load the modelling internals.
- Unified the bootstrap metric outputs so LCC, Pearson correlation, and accuracy share a consistent structure across serial and parallel runs.
- Added documented bootstrap confidence intervals to the user workflow and refreshed package documentation accordingly.

## New Features

- Support for multiple bootstrap schemes, reproducible seeding, and normal, percentile, or BCa confidence intervals through a unified builder.
- Expanded plotting controls, including interval styling, point alpha, and axis scaling configured via plotControl helpers.
- New package options (boot.scheme, ci.method, keep.boot.models, boot.seed, numCore) to fine-tune bootstrap execution and persistence.

## Bug Fixes

- Added stronger input validation, improved numerical stability in degenerate scenarios, and harmonised metric list shapes for downstream consumers.

## Documentation and Tests

- Updated README, reference manuals, and CRAN notes to reflect the enhanced bootstrap and plotting workflows.
- Added focused tests covering bootstrap dimensions, metric integrity, and plot control behaviour.

## Compatibility

- No breaking user-facing changes; existing APIs continue to work, with the major version signalling the enhanced bootstrap capabilities and infrastructure improvements.
