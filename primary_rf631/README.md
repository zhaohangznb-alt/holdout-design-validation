# Primary RF statistics

Exact archived statistics for the primary RF analysis in Tables S3/S4 of *Holdout Design and Utility Estimates in ChEMBL Bioactivity Prediction*.

- `run_metrics.csv`:288 seed-specific rows, including matching training-mean and training-median baselines, R-squared and per-seed RF errors.
- `conditional_document_bootstrap.csv`:96 primary whole-document conditional intervals. Primary and paired-operational views are retained.
- `FILE_MANIFEST.json`:SHA-256 for these exact copies. `bootstrap_provenance.json` identifies the archived protocols and their seeds; local original paths are historical provenance.

These files extend the public materials while leaving `verified_bundle_v1.tar.xz` unchanged. Cross-model statistics are separately under `bundle/results/analysis/` after extraction. Strict MET year-scaffold S4=[-0.0048733095835785225,0.38998809372401705], S20=[0.0010426538016099756,0.39919870256427065]. Different bootstrap seeds explain different endpoints; neither interval was recomputed for this supplement. The main positive-year-interval statement concerns exact-identity evaluation.
