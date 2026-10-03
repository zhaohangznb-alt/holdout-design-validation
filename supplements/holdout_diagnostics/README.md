# Frozen-prediction holdout diagnostics

Figure 6 / SI Table S21: all 11 broad-cohort targets, four frozen routes, 220 existing prediction files. No new fitting, tuning, bootstrap or target selection. `ANALYSIS_PLAN.json` was frozen before computation. The 51,900 test-context appearances are repeated route evaluations, not 51,900 unique measurements.

## Files

- `context_diagnostics.csv`:44 contexts with training/test counts, label distributions, structural proximity and model/baseline errors.
- `target_route_diagnostics.csv`:33 within-target contrasts; separate log RF and baseline MSE changes, their exact difference and skill changes for all three models.
- `fixed_bin_diagnostics.csv`:660 model/context/bin summaries;39 empty bins keep blank metrics.
- `row_similarity.csv.gz`:all test-context nearest training UID, label and maximum similarity. Standard gzip, uncompressed SHA-256 in `PROVENANCE.json`.
- `Figure_6.pdf` and the Figure/Table data CSVs:manuscript display values and full-precision sources.
- `VERIFICATION.json`:executed independent Windows checks, including132 integer bit-count nearest-neighbour checks.

## Recompute from the public frozen bundle

```sh
python unpack_bundle.py
python -m pip install -r bundle/requirements.txt
python -B bundle/reproduce.py --output bundle/diagnostic_prepared --prepare-only
python -B supplements/holdout_diagnostics/diagnostics.py --bundle bundle --prepared bundle/diagnostic_prepared --output diagnostic_recompute
python -B supplements/holdout_diagnostics/verify_diagnostics.py --bundle bundle --prepared bundle/diagnostic_prepared --output diagnostic_recompute
```

Use the pinned Python3.9.10/numerical dependencies. Preparation is verified on Windows; the Linux representative wrapper separately reports compatibility and actual model execution. The diagnostics entry loads only frozen inputs and predictions. An explicit `--resume` may reuse hash-bound structural checkpoints after interruption; it never fits models. Raw source-label parsing matches the original runner. RF MSE averages row losses per seed before taking a median, including within similarity bins. Equal-context weighting is retained. The error log identity is descriptive and does not estimate causal contributions of similarity or training size.

`EXECUTED_LOCAL_analysis_pipeline.py` records the executed Windows implementation and historical absolute paths. `diagnostics.py` is its portable path adapter, checked separately. To redraw the figure, install Matplotlib3.10.8 in a separate rendering environment and run `make_figure.py --output diagnostic_recompute`. The model environment is unchanged.

Original ChEMBL-derived records retain CC BY-SA3.0 attribution and terms. Supplement scripts use the repository's MIT licence. No publisher source originals are included.
