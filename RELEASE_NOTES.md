# Version 1.0.0

Frozen-input reproduction and fixed-model holdout validation across six kinase and five additional non-kinase targets from ChEMBL 37.

- All 680 prediction cases and all outcomes are included.
- Current retained inputs rebuild byte-for-byte; newly trained RF predictions reproduce the archived RF values.
- RF, Ridge and Tanimoto 1-NN specifications were fixed before test evaluation. Candidate selection did not use model scores.
- Includes frozen API records, eligibility and partition records, per-row predictions, absolute errors, conditional document intervals, fit ledgers and scientific figures.
- The offline entry rebuilds inputs, freshly fits models, verifies predictions and recomputes analysis. Full Windows model runs and a separate two-case clean-directory smoke check are documented.

The complete frozen payload is the tracked `verified_bundle_v1.tar.xz`, with SHA-256 in `BUNDLE_SHA256.txt`. ChEMBL-derived data retain CC BY-SA 3.0; the listed author-created scripts use MIT.
