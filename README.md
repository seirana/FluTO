# FluTO — Flux Trade-Off Identification in Metabolic Networks

FluTO is a MATLAB implementation of a constraint-based method for identifying absolute flux trade-offs in metabolic networks.

The method is associated with:

> Hashemi, S., Razaghi-Moghadam, Z., & Nikoloski, Z. (2021). **Identification of flux trade-offs in metabolic networks.** *Scientific Reports*, 11, 23776. https://doi.org/10.1038/s41598-021-03224-9

## Repository status

The repository now contains two clearly separated implementations:

```text
src/+fluto/       maintained, tested implementation
legacy/codes/     original publication-era MATLAB scripts
```

The publication code is preserved for provenance. New analyses should use the maintained package.

## What FluTO does

FluTO studies reaction sets whose fluxes cannot vary independently under a constrained metabolic model.

A typical workflow is:

```text
COBRA model
   |
   v
explicit biological condition
   |
   v
reaction-wise feasible flux ranges
   |
   v
reaction orientation + flux classification
   |
   v
fully-coupled reaction analysis
   |
   v
MILP trade-off enumeration
   |
   v
structured tables + reproducibility metadata
```

The maintained code separates these stages into testable functions rather than mixing solver configuration, model mutation, optimization, and file output in one script.

## Senior-level architecture

The public API lives in the MATLAB package namespace `fluto`:

```text
src/+fluto/
├── validateModel.m
├── applyCondition.m
├── computeFluxRanges.m
├── canonicalizeReactionDirections.m
├── removeBlockedReactions.m
├── classifyFluxes.m
├── buildFullyCoupledMatrix.m
├── expandCoupledAlternatives.m
├── enumerateTradeoffs.m
├── runCondition.m
└── writeResults.m
```

Design principles:

- explicit inputs and outputs;
- no hard-coded user paths;
- no silent global solver changes;
- namespaced error identifiers;
- input and dimension validation;
- explicit numerical tolerances;
- checked optimization status;
- deterministic support de-duplication;
- sparse MILP matrices where appropriate;
- computation separated from persistence;
- machine-readable diagnostics and metadata;
- publication-era behavior archived instead of silently overwritten.

## Important modernization fixes

Two concrete defects were identified in the original MATLAB implementation.

### Fully-coupled reaction check

The original `makeFCMatrix.m` repeated the same constrained optimization in its second check and then compared the duplicate result against the first result.

The maintained `fluto.buildFullyCoupledMatrix` performs the reciprocal check explicitly:

1. fix reaction *j* and test whether reaction *i* is uniquely determined;
2. fix reaction *i* and test whether reaction *j* is uniquely determined.

### Fully-coupled alternative expansion

The historical trade-off expansion only used reactions returned from the coupling-matrix row. The selected reaction itself was not included, so an uncoupled selected reaction could produce an empty Cartesian product.

The maintained `fluto.expandCoupledAlternatives` always includes the selected reaction itself plus its fully-coupled alternatives.

See [MIGRATION.md](MIGRATION.md) for details.

## Dependencies

For the full analysis:

- MATLAB R2021a or newer;
- COBRA Toolbox;
- a COBRA-compatible LP solver configured by the user;
- Optimization Toolbox for `intlinprog`.

The maintained package **does not** call `changeCobraSolver` automatically. Solver choice is environment-specific and should be configured before running FluTO.

The CI test suite is intentionally solver-independent and therefore does not require COBRA Toolbox or Optimization Toolbox.

## Setup

Clone the repository and add the maintained source folder:

```matlab
addpath("src")
```

Initialize COBRA Toolbox separately using the setup recommended by your installed COBRA version.

## Model input

A FluTO model must contain the COBRA fields:

```text
S
rxns
mets
lb
ub
c
```

`fluto.validateModel` checks dimensions, finite values, bounds, and reaction/metabolite identities. Stable `rxnNumber` and `metNumber` fields are added if absent.

Model snapshots already included in the repository are under `models/`.

For new research, record the model source, release, retrieval date, preprocessing, and exact file checksum. See [models/README.md](models/README.md).

## Conditions are explicit configuration

The original E. coli script embedded reaction numbers and fixed values directly in code.

The maintained API uses a condition struct instead:

```matlab
condition = struct( ...
    "blockReactionNumbers", [101, 102], ...
    "activeReactionNumber", 101, ...
    "activeBounds", [-10, -10], ...
    "fixedReactionNumbers", [200, 716, 721, 723], ...
    "fixedFluxValues", [-0.9476, 3.15, 0.0931, 12.612]);
```

These numbers are examples of the configuration format. They must match the model and biological condition being analyzed.

A JSON template is provided at:

```text
examples/condition.template.json
```

## Programmatic workflow

```matlab
addpath("src")

model = readCbModel("models/Ecoli_iJO1366.mat");

condition = struct( ...
    "blockReactionNumbers", [101, 102], ...
    "activeReactionNumber", 101, ...
    "activeBounds", [-10, -10]);

result = fluto.runCondition(model, condition);
```

The returned struct includes:

- validated/reduced model;
- flux ranges;
- flux classifications;
- FVA/removal diagnostics;
- fully-coupled matrix;
- coupling diagnostics;
- trade-off table;
- support matrix;
- solver-run counts;
- enumeration options.

## File-driven runner

A convenience runner is available:

```matlab
addpath("scripts")

output = runFluTO( ...
    "models/Ecoli_iJO1366.mat", ...
    "examples/condition.template.json", ...
    "artifacts/example");
```

It records SHA-256 hashes of the model and condition files together with MATLAB version information.

## Flux classification

After feasible ranges are computed and negative irreversible reactions are put into a canonical positive orientation, reactions are classified as:

- `fixed`;
- `variable` — fixed-sign variable in the paper;
- `reversible` — sign-variable in the paper.

Classification uses an explicit numerical tolerance rather than decimal rounding.

## Fully-coupled reactions

The maintained coupling routine performs reciprocal constrained LP checks and then computes the transitive closure of the fully-coupled relation.

It returns both the logical coupling matrix and diagnostics describing:

- eligible variable reactions;
- number of pairs tested;
- directly coupled pairs;
- pairs after transitive closure;
- numerical tolerance.

## Trade-off enumeration

`fluto.enumerateTradeoffs` implements the maintained MILP search.

The defaults preserve important publication-era formulation choices:

- minimum degree: 2;
- maximum degree: 9;
- Big-M: 1001;
- coefficient bound: ±100;
- integer decision variables.

Unlike the historical loop, the maintained implementation:

- checks `intlinprog` exit status;
- handles infeasibility separately from unexpected solver failures;
- removes duplicate supports deterministically;
- expands fully-coupled alternatives safely;
- returns a structured result rather than writing Excel files during every search iteration.

## Tests

Run:

```matlab
addpath("src")
results = runtests("tests");
assertSuccess(results)
```

The solver-independent tests cover:

- model validation;
- invalid-bound rejection;
- irreversible-reaction canonicalization;
- tolerance-aware classification;
- reaction-number-based condition application;
- unknown-reaction validation;
- fully-coupled support expansion;
- prevention of duplicate support members.

GitHub Actions runs these tests on a clean MATLAB runner and verifies that the maintained package entry points parse and resolve.

GitHub Actions uses the official `matlab-actions/setup-matlab@v3` and `matlab-actions/run-command@v3` workflow.

## Repository layout

```text
.
├── src/+fluto/          maintained package
├── tests/               MATLAB unit tests
├── scripts/             reproducible runners
├── examples/            explicit configuration templates
├── legacy/codes/        original implementation
├── models/              model snapshots
├── figures/             publication figures
├── supplemantary/       historical supplementary material
├── README.md
├── MIGRATION.md
├── CONTRIBUTING.md
└── CITATION.cff
```

The historical directory name `supplemantary/` is retained to avoid silently moving publication artifacts.

## Scientific interpretation

FluTO identifies trade-offs implied by the supplied stoichiometric model and constraints.

A result therefore depends on:

- model reconstruction;
- reaction bounds;
- nutrient/environment conditions;
- fixed-flux assumptions;
- solver tolerances;
- preprocessing choices.

A computational trade-off should not be described as experimentally validated unless independent evidence supports that claim.

## Reproducibility

For a reported experiment, preserve:

- Git commit SHA;
- MATLAB release;
- COBRA Toolbox version;
- LP/MILP solver and version;
- model source/version/checksum;
- explicit condition configuration;
- numerical tolerances;
- enumeration options;
- generated result tables and metadata.

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff).

## License

No explicit software license file is currently included in this repository. Public visibility alone does not grant reuse rights.
