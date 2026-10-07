# Holdout-design validation

Companion to *Holdout Design and Utility Estimates in ChEMBL Bioactivity Prediction*.

## V12 manuscript evidence supplement (package 1.1.0)

The [new evidence archive](supplements/manuscript_evidence_v12/README.md) adds F1 full Windows clean-directory replay records, F2 training-document resampling, F5 label/period controls, F8 weighting, phase-3 extensions, limited-grid and targeted checkpoint sensitivities, and conditional release-forward intervals. It supplies data, protocols, project code and artwork; the manuscript, full upstream inputs and model weights are not deposited. See its component licenses, current file manifest and 81-display index. Earlier archives below retain their original bytes and immutable citations.

The historical API cache is not directly bound to an immutable ChEMBL release, despite a contemporaneous status response naming ChEMBL 37. The separately versioned 36/37 transition has a distinct provenance chain. The new evidence records full Windows clean-directory retraining; full Linux 680-case retraining remains unperformed. No journal submission or manuscript approval is implied. No Zenodo DOI exists in these metadata until an actual record has been verified.

The complete offline package is `verified_bundle_v1.tar.xz` (SHA-256 in `BUNDLE_SHA256.txt`). It contains frozen inputs, all 680 prediction cases, analysis, fit ledgers and verification records. Root-level scripts are provided for browsing; run the pipeline from the extracted `bundle/` directory.

```sh
python unpack_bundle.py
python -m venv .venv
# Activate .venv, then:
python -m pip install -r bundle/requirements.txt
python -B bundle/verify_release.py
python -B bundle/reproduce.py --output bundle/reproduced_v1
```

A shorter integration run is `python -B bundle/reproduce.py --output bundle/smoke_v1 --quick-check`.

Full Windows runs produced 680 prediction cases; a separate clean-directory smoke check rebuilt all 408 input tables byte-for-byte and freshly fitted two Ridge cases, reproducing both prediction files byte-for-byte. This smoke check has a separate scope from the full runs. See `bundle/verification/clean_directory_check.json`.

Publication state is tracked in the repository-root `PUBLISH_STATUS.json`. The copy inside the sealed bundle records the state at packaging time.

![Broad-cohort paired changes](paired_changes_revised.png)

![Operational-cohort paired changes](paired_changes_single_protein.png)

## Subsequent verified supplements

- [Common-test matched-training-size control](supplements/matched_train_control/README.md): all11 targets,4,814 shared test rows,five fixed size-matched source-pool draws. New165 RF fits,55 Ridge fits and55 NN indices;110 full-pool prediction files reused. RF year-pool MSE remains higher in9/11targets (exceptions EGFR/PTGS2). Figure7,TableS22 and complete inputs/predictions/scripts are in the independently verifiable supplement.

- [Primary RF statistics](primary_rf631/README.md):exact288 seed statistics and96 original conditional document intervals used by SI Tables S3/S4; files, protocol versions and hashes are public.
- [Representative replay](supplements/representative_replay/README.md):actual Windows/Linux408-table reconstruction, six fresh RF fits and two1-NN indices per platform; Linux run37128892159 passed.
- [Frozen-prediction diagnostics](supplements/holdout_diagnostics/README.md):44 contexts/11 targets,220 frozen prediction files,51,900 repeated test-context appearances, Figure6 and TableS21 with full-precision data. All11 year-exact targets have smaller training sets and lower median structural proximity. PTGS2 has lower RF absolute MSE but lower skill because its baseline improves more.

These additions have separate manifests. `verified_bundle_v1.tar.xz` retains its original bytes. The diagnostic data commit is3553db3a62e7ff696595cac2a6e6a6f559e4dd43; executed Linux code isab031551cfba56b8982d3563221ca04a2d590923. Detailed run results are retained in the supplement. Zenodo linkage and a version DOI remain pending.

## Main findings

In broad cohorts, median within-target changes in skill (document-year exact minus source-document exact) were−0.431,−0.967,−0.149 for RF,Ridge,Tanimoto1-NN on the six kinases; on the five non-kinase targets they were−0.250,−0.592,−0.314. RF changed in the same direction in6/6 and5/5 targets. Directions vary for the other models. In the non-kinase operational cohort, Ridge's median scaffold-filter change was larger than its year-route change. See `REPORT_CN.md`, `results/analysis/` and `figures/` for all results, absolute errors and exceptions.

The new targets are ACHE,DPP4,PTGS2,HDAC1,PDE4D. CA2 failed the frozen data-count criteria and was not replaced. Selection never used model scores. These are new targets from the same resource; models were trained separately for each target. Data follow database annotations, rather than a new publisher-original audit.

## Contents

- `inputs/current_origin/`: original frozen train/valid/test CSVs, cached assay metadata and initial Policy A.
- `inputs/current_reference/`: current split and training manifests, RF configuration and631-item eligibility policy. Original machine paths are provenance fields; execution resolves package-relative names.
- `inputs/external_raw/`: official API JSON, query URLs, retrieval timestamps and hashes, including all six screened candidates.
- `results/predictions/`:680 prediction cases:480 current-panel cases (including paired-operational sensitivities) and200 new-panel cases.
- `results/analysis/`: row-level losses, target/model summaries, route contrasts, absolute errors, conditional document intervals and paired model intervals.
- `results/fit_ledgers/`: configuration, input signatures, actual fit/index-build times and binary hashes. Model binaries are rebuilt by the training entry and omitted from this bundle.
- `verification/`: input reconstruction and independent numerical checks.
- `PROTOCOL.md`, `EXTERNAL_PROTOCOL.md`, `BOOTSTRAP_PROTOCOL.md`: primary model protocol, score-blind candidate/partition rules, and supplementary conditional interval protocol.
- `FILE_MANIFEST.json`: SHA-256 for the sealed package payload.

## Environment

The full model runs used Windows, Python3.9.10, RDKit2025.09.2, scikit-learn1.6.1, numpy2.0.2, pandas2.3.3, scipy1.13.1 and joblib1.5.3. Dependencies are pinned in `requirements.txt`. Use a dedicated environment. No GPU is required. RF uses two threads; other numerical workers use one thread. The fresh-model runs completed138 RF fits,46 Ridge fits and46 nearest-neighbor indices; repeated contexts reuse identical ordered training inputs/configurations/seeds.

Create an environment, activate it and install the pinned requirements:

```sh
python -m venv .venv
python -m pip install -r requirements.txt
```

The second command must use the activated environment's Python. The validated Python/version details and clean-directory scope are in `verification/clean_directory_check.json`. A bounded Linux replay is verified in the subsequent supplement; full680-case Linux training has not been run.

## Offline end-to-end reproduction

After changing into the extracted `bundle/` directory, use a new output directory for each invocation:

```sh
python -B verify_release.py
python -B reproduce.py --output reproduced_v1
```

The full command rebuilds all288 current input tables from original partitions and policy, prepares the new-target cohorts from cached raw API JSON, checks split integrity, freshly fits RF/Ridge and builds1-NN indices, writes680 prediction cases, independently verifies their metrics and computes2000-replicate conditional whole-document intervals. It uses no network and no archived model binaries.

For faster input validation or a two-case integration smoke check:

```sh
python -B reproduce.py --output preparation_v1 --prepare-only
python -B reproduce.py --output smoke_v1 --quick-check
```

The smoke check freshly fits one Ridge case in each panel and verifies it; it is not a full rerun. Checkpointed direct training may use an explicit `--resume` only after diagnosing an interruption. A live lock or incomplete fit fails rather than silently repeating a fit.

To regenerate the scientific figures, additionally install `matplotlib==3.9.4` and run:

```sh
python -B plot_results.py --analysis results/analysis --output rebuilt_figures
```

## Interpretation and data boundary

Training-mean baselines always use the matching training cohort. RF target results are medians over seeds42,53,67; deterministic models occur once. Panel summaries weight targets equally. Median paired target changes and differences of panel medians are different summaries. Route contrasts compare different test populations; cross-model contrasts within a context use identical rows. Conditional document intervals hold the cohort and fitted models fixed. The Ridge alpha and1-NN rule were fixed before testing, and their results do not represent tuned family optima.

Current-panel reproduction starts from archived partitions and accepted source policies; the new panel starts from frozen official API responses. Full raw ChEMBL database reconstruction and publisher-original certification are outside this package's execution boundary. No biological experiment was performed.

## Licensing and citation

ChEMBL-derived records and metadata retain [CC BY-SA3.0 Unported](https://chembl.gitbook.io/chembl-interface-documentation/about) attribution and terms; see `DATA_LICENSE.md`. Newly authored scripts listed in `CODE_LICENSE_SCOPE.csv` use MIT; see `LICENSE_CODE.txt`. Third-party dependencies retain their own licenses. No publisher PDFs or earlier private review conversations are redistributed.

Data source: ChEMBL37, release DOI [10.6019/CHEMBL.database.37](https://doi.org/10.6019/CHEMBL.database.37). API documentation: [ChEMBL Data Web Services](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services). `CITATION.cff` describes this companion package. A package DOI is added only after the actual Zenodo record exists.

The intended publication route is the author's GitHub repository, with a versioned Zenodo archive. See `PUBLICATION.md` for the release sequence and the actual `PUBLISH_STATUS.json` for publication status.

