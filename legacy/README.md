# Legacy publication implementation

This directory preserves the original MATLAB implementation associated with the 2021 FluTO publication.

It is retained for provenance and comparison, but it is **not** the maintained software interface.

The original code contains characteristics typical of exploratory research scripts, including:

- hard-coded relative paths;
- global solver/path setup in scripts;
- model-specific reaction numbers embedded in functions;
- repeated dynamic array growth;
- limited solver-status checking;
- mixed computation and file I/O;
- legacy naming/typos retained for historical traceability.

Use the maintained package under:

```text
src/+fluto/
```

for new analyses.

Two concrete defects identified during modernization are documented in [../MIGRATION.md](../MIGRATION.md): the duplicated second fully-coupled optimization check and omission of a selected reaction itself when expanding fully-coupled alternatives.
