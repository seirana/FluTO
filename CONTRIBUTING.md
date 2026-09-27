# Contributing

Changes to the maintained FluTO implementation should follow these rules.

## MATLAB style

- Put reusable implementation in `src/+fluto/`.
- Prefer functions over workspace scripts.
- Validate public-function inputs and fail with namespaced error identifiers.
- Keep numerical tolerances explicit.
- Do not hard-code user paths or mutate the MATLAB path inside library functions.
- Do not silently change the active COBRA solver.
- Check solver status before using optimization results.
- Separate computation from file output.
- Keep random behavior seeded and explicit if randomness is introduced.
- Prefer tables/structs with named fields over position-dependent cell arrays.

## Scientific changes

Any change that can alter identified trade-offs must document:

1. the mathematical assumption being changed;
2. why the change is necessary;
3. a regression or synthetic test;
4. whether publication-era output is expected to change.

Do not label a computational prediction as experimentally validated unless the repository contains evidence supporting that statement.

## Tests

Run:

```matlab
addpath("src")
results = runtests("tests");
assertSuccess(results)
```

Tests in CI intentionally avoid requiring COBRA Toolbox or Optimization Toolbox so core validation logic remains continuously testable on a clean MATLAB runner.
