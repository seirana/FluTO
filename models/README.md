# Model snapshots

This directory contains the model snapshots that were already part of the FluTO repository:

- `ArabidopsisCoreModel.mat`
- `Ecoli_iJO1366.mat`
- `yeastGEM.mat`

The filenames alone are not sufficient to reconstruct complete biological provenance.

For any new published analysis, record:

- original database/source;
- model identifier and release/version;
- retrieval date;
- any modifications applied before FluTO;
- exact Git commit of this repository;
- file checksum.

The maintained FluTO package does not silently modify biomass reactions or nutrient conditions. Those choices belong in an explicit condition/configuration step.
