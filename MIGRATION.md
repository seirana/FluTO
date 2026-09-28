# Modernization notes

The original 2021 code is preserved under `legacy/codes/`. The maintained implementation lives under `src/+fluto/`.

## Architecture change

Historical flow:

```text
script
  -> mutate MATLAB path / solver globally
  -> hard-coded model and condition
  -> FVA
  -> classification
  -> coupling
  -> MILP enumeration
  -> Excel output inside the search loop
```

Maintained flow:

```text
validated model
  -> explicit condition struct
  -> checked FVA
  -> canonical reaction orientation
  -> tolerance-aware classification
  -> symmetric coupling test
  -> checked MILP enumeration
  -> structured result
  -> optional result writer
```

Computation is separated from persistence, and functions return diagnostics instead of relying on console output or workspace state.

## Scientific behavior preserved

The maintained MILP keeps the historical decision-variable layout, integer-variable choice, Big-M default of 1001, coefficient bounds of ±100, and trade-off degree search beginning at 2.

The original source remains available so publication-era runs can be compared directly.

## Concrete bug fixes

### 1. Fully-coupled reaction check

The original `makeFCMatrix.m` performed the same constrained optimization twice and then required the second result to be larger than the first. Because the second block repeated the same objective and fixed reaction, that final comparison could not represent the intended reciprocal coupling test.

The maintained `fluto.buildFullyCoupledMatrix` performs an explicit symmetric check:

1. fix reaction j at an interior feasible value and test whether reaction i is uniquely determined;
2. fix reaction i at an interior feasible value and test whether reaction j is uniquely determined.

### 2. Expansion of fully-coupled alternatives

The original trade-off expansion used only entries found in a coupling-matrix row. The selected reaction itself was not included, so a reaction with no coupled alternatives could produce an empty Cartesian product.

The maintained `fluto.expandCoupledAlternatives` always includes the selected reaction itself and then adds its fully-coupled alternatives.

## Deliberate engineering changes

- hard-coded model condition changes are replaced by `fluto.applyCondition`;
- numerical equality uses tolerances instead of decimal rounding;
- COBRA solver status is checked;
- duplicate trade-off supports are removed deterministically;
- output files are written after computation rather than during every MILP iteration;
- external dependency setup is left to the caller instead of changing global MATLAB state.

## Compatibility

The maintained package expects:

- MATLAB R2021a or newer;
- COBRA Toolbox for model I/O and LP/FVA operations;
- Optimization Toolbox for `intlinprog`.

CI exercises solver-independent code with MATLAB's unit-test framework. Full end-to-end model runs additionally require the COBRA Toolbox and a configured LP solver.
