# V12 manuscript evidence supplement — package version 1.1.0

Data, author-created code, protocols, source locators, numerical results and artwork supporting *Holdout Design and Utility Estimates in ChEMBL Bioactivity Prediction*. The manuscript and full Supporting Information PDF/TeX are not in this deposit. The data/code publication does not establish final author approval or journal submission of the manuscript.

## What is supplied

Use `DISPLAY_SOURCE_INDEX.csv` (81 manuscript displays), `DIGITAL_TABLE_DICTIONARY_SUBMISSION.md`, `DIGITAL_TABLE_DICTIONARY_V12.md`, `source_provenance/SOURCE_EVIDENCE_STATUS.json` and `PUBLIC_DISTRIBUTION_SCOPE.json`. Numeric table excerpts are supplied as data attachments; three unreferenced manuscript-history prose excerpts are retained locally. Historical dictionaries, indices and manifests describe their named versions, not the current payload. `SHA256_MANIFEST.json` is the current manifest.

Results include full Windows clean-directory replay records (F1), training-document resampling (F2), common-test label/period controls (F5), evaluation-unit weighting (F8), phase-3 D-MPNN/42-target panel/post hoc REML/release-transition analyses, limited-grid and targeted checkpoint sensitivities, and post hoc conditional test-document intervals. These analyses have different inferential scopes; see their protocols. Historical and post hoc records are not converted into preregistered analyses by publication.

The original 680-case offline pipeline, frozen inputs and separately published primary RF statistics, representative Windows/Linux replay, diagnostics and size-matching control remain in the repository at their original immutable commits and original hashes. This new archive supplies additional evidence; it is not a full ChEMBL database, all upstream training inputs, or a model-weight distribution. The full Windows replay is distinct from the bounded representative Linux replay; full 680-case Linux retraining remains unperformed.

## Verify and replay without retraining

Extract `manuscript_evidence_v12.tar.xz`, change into `manuscript_evidence_v12/`, then:

```sh
python VERIFY_PACKAGE.py
python forward_uncertainty/forward_uncertainty.py --root /absolute/path/to/manuscript_evidence_v12
```

The first command uses the Python standard library. The second requires existing NumPy and pandas, uses supplied document-loss sums with fixed fits/cohorts, and writes `data/forward_conditional_intervals.csv` (90 intervals) and `posthoc_descriptive/FORWARD_UNCERTAINTY.json`. Do not add `--prepare`; it requires author-local upstream inputs. Other legacy training/plotting scripts also retain their documented upstream dependency requirements. Their inclusion does not certify portable full retraining. No new fits were performed for this distribution.

The replay was byte-identical for both outputs using NumPy 2.3.5 and pandas 3.0.1. Whole-test-document intervals condition on the fitted models and observed documents, not training or target-population uncertainty.

## Provenance limitations

A contemporaneous API status reports ChEMBL 37, but the historical BRAF retrieval manifest does not directly bind its activity pages to an immutable release. The separately versioned ChEMBL 36/37 extension has a distinct provenance chain. Harris Table 1 numbers match the inspected local original, but biochemical/cellular endpoint ambiguity remains. Publisher PDFs and their page renders are excluded. Initial 27h/27i retention-decision timing and the earliest adopted AI-use date are not authenticated.

Author-local path strings are provenance locators, not public downloadable evidence. A hash does not supply an absent upstream file. Private task correspondence and verbatim acceptance wording are excluded; the previous distribution ledger records those transformations. All scientific CSV fields are preserved. No author approval, structure certification or missing publisher-source verification is inferred.

## Licenses and citation

ChEMBL-derived data retain attribution and CC BY-SA 3.0; see `LICENSE_AND_ATTRIBUTION.md`. Newly authored Python scripts are listed with hashes in `CODE_LICENSE_SCOPE.csv` and use MIT under `LICENSE_CODE.txt`. This does not relicense third-party dependencies, chemical applications, fonts or publisher material. Author-created scientific figures accompany the dataset under CC BY-SA 3.0, subject to third-party font/software rights. No font program or third-party package source is distributed.

Use the immutable repository commit and archive hash to cite this package until an actual Zenodo record exists. The repository-root publication status reports the live archive/DOI state; draft metadata or a pending association is not a DOI.
