# Common-test, matched-training-size comparison

This extension fixes the shared source-exact/year-exact test intersection for each of eleven targets and the number of training measurement-context rows. The original pools, eligibility, labels and model settings remain unchanged. The source pool is larger in every target and is sampled five times with predeclared seeds; the smaller year pool and both full-pool references reuse archived predictions. The protocol was fixed before calculating new controlled outcomes in `ANALYSIS_PLAN.json`.

Run with the existing pinned scientific environment (no new dependencies):

```powershell
& E:/paper/kinase_followup_execution/.venv-win/Scripts/python.exe experiment.py --run
& E:/paper/kinase_followup_execution/.venv-win/Scripts/python.exe experiment.py --summarize
& E:/paper/kinase_followup_execution/.venv-win/Scripts/python.exe verify_experiment.py --stage final
```

Preparation has already run and should not be repeated. Completed predictions/fits are hash-bound and reused on continuation. An incomplete fit ledger or unledgered model stops for diagnosis. `RUNNING.lock` identifies the running process; do not remove it while that process is alive. RF and Ridge are fitted once per new ordered sample/model/seed. NN builds a fresh ordered index. The 110 full-pool prediction files are reused, not new fits.

`POOL_COUNTS.csv` includes all eleven targets and both original test populations' coverage. The common test has 4,814 rows across targets, ranging from 23 to 1,416. These are aggregated measurement contexts and not unique molecules. ERBB2, PTGS2 and PDE4D have 23, 35 and 52 test rows. Main inference describes these shared intersections, rather than every row in the original test populations.

Primary outcomes are absolute model MSE on identical rows. Each arm also retains its training-mean MSE and skill. Normalized contrasts use the single common-test label variance as a descriptive scale; its test-mean constant is an oracle reference, not a deployable baseline. RF is the median of three seed-specific mean row losses. The five subsampling medians/minima/maxima are sampling summaries, not confidence intervals. There are no new bootstrap intervals, significance tests, target selection or tuning.

The per-draw full-pool MSE gap equals the source row-count sensitivity plus the matched-size pool contrast, since the year pool remains whole. Independently summarized medians are not required to sum. Residual pool contrasts include structural, label, document and chronological differences and do not identify a pure structural effect. Structural moments are diagnostic outcomes; they are not used for sample selection.

The original manuscript and figures are preserved in the preceding delivery. `02_editable_source` and `01_manuscript_pdf` contain the updated independent revision with Figure 7 and Table S22. `analysis/` contains full-precision numeric results. Model binaries in `models/` remain local and are excluded from delivery archives; their hashes and execution histories are retained in the fit ledgers. The portable numerical snapshot contains all prepared inputs and predictions and can be verified without those binaries. Fresh reproduction writes to a new directory using the supplied bootstrap script, retaining the full-pool archived prediction references.

Author-created code is MIT under the companion repository's code scope. ChEMBL-derived data retain attribution and CC BY-SA 3.0 terms; database release DOI: https://doi.org/10.6019/CHEMBL.database.37. Public companion: https://github.com/zhaohangznb-alt/holdout-design-validation. A Zenodo version DOI is not yet assigned.

Routing: the fixed-Astra Codex adapter selected medium effort with classifier confidence 0.77. It reached its 300-second limit after writing the frozen plan and checking common-test feasibility. The parent continued that plan after the failed attempt; no effort/model escalation or model-fitting repeat occurred.

## Portable numerical snapshot and fresh reproduction

Extract `controlled_comparison_v1.tar.xz` into a new directory. The extracted `controlled_comparison/` contains prepared inputs, all predictions, executed scripts and per-job ledgers. Then run:

```sh
cd controlled_comparison
python verify_snapshot.py
# NEW destination; executes 165 RF fits,55 Ridge fits,55 NN indices
python fresh_reproduce.py --output ../new_fresh_control
```

Use Python3.9.10 and the pinned numerical dependencies. `verify_snapshot.py` only verifies data/hashes and recalculates metrics; it does not fit. `fresh_reproduce.py` copies only the required archived full-pool prediction references and immutable prepared inputs into a new directory, then performs fresh matched-size fits. That portable full fresh command is provided; an additional independent full fresh run has not been executed. The current controlled experiment was executed on Windows. The earlier Linux replay is a separate eight-case check.
