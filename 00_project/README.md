# IBD EcoTyper Project

## Directory policy

- `01_raw_data/`
  Original downloaded data. Treat as immutable.

- `02_metadata/`
  Canonical sample, donor, clinical and ID mapping tables.

- `03_reference/`
  Single-cell reference audits, taxonomy mapping and reference construction.

- `04_CIBERSORTx/`
  Signature matrix, Fractions and HiRes outputs.

- `05_EcoTyper/`
  EcoTyper discovery, cell-state and ecotype analyses.

- `06_validation/`
  Independent and sensitivity validation analyses.

- `07_spatial_ecotyper/`
  Spatial EcoTyper analyses.

- `08_results/`
  Final figures, tables and exports.

- `scripts/`
  Analysis code only.

- `environments/`
  Environment specifications and container definitions.

- `tmp/`
  Temporary/intermediate files that may be deleted.

## Core rule

Never modify or overwrite files under `01_raw_data/`.
Derived data should be written to the corresponding downstream directory.
