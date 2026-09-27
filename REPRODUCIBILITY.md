# Reproducibility guide

A FluTO result is defined by more than the MATLAB source code. The feasible flux space depends on the exact metabolic model, bounds, environmental condition, solver stack, and numerical tolerances.

## Minimum experiment record

Archive the following together:

- repository Git commit;
- MATLAB release;
- COBRA Toolbox version/commit;
- LP solver name and version;
- Optimization Toolbox/MATLAB version used by `intlinprog`;
- model filename and SHA-256 checksum;
- model source/database release;
- condition JSON/struct;
- FVA/coupling tolerance;
- trade-off enumeration options;
- result CSV files;
- metadata JSON.

## Solver configuration

The maintained package never calls `changeCobraSolver` on behalf of the user.

Configure COBRA Toolbox explicitly before running an experiment and record that configuration. This avoids hidden global state inside FluTO itself.

## Numerical tolerances

The original scripts often used decimal rounding to compare optimization results. The maintained implementation uses explicit numerical tolerances.

Report the tolerance used for:

- flux-range classification;
- zero-range/blocked detection;
- fully-coupled reaction equality.

Small tolerance changes can alter classifications near zero and therefore change downstream enumeration.

## Condition provenance

Do not encode environment-specific reaction numbers directly inside library functions.

Store them in a condition file or struct and archive it with the run. Verify that reaction numbers refer to the exact model snapshot used in that run.

## Output provenance

The file-driven runner records SHA-256 hashes for the model and condition files.

If the model is modified in memory before calling `fluto.runCondition`, save that exact modified model or otherwise record the transformation steps because a source-file checksum alone will no longer describe the analyzed object.

## Publication-era comparison

The original MATLAB scripts are preserved under `legacy/codes/`.

When a modernization changes results, compare:

1. the exact same model and bounds;
2. the same solver;
3. the same condition;
4. publication-era code;
5. maintained code;
6. the documented bug fixes in `MIGRATION.md`.

Do not attribute a difference to the method before ruling out model, solver, or condition differences.
