# Representative Windows and Linux replay

The Windows replay rebuilt408 partition inputs byte-for-byte, completed six fresh RF fits (seeds42/53/67) and two deterministic Tanimoto1-NN indices for ALK and ACHE under main/revised/source_exact. RF differences were at most2.6645352591003757e-15; both NN CSVs were byte-identical. Earlier two Ridge checks are separate.

The Linux container replay also rebuilt all408 inputs byte-for-byte and completed the same six fresh RF fits and two1-NN indices. RF prediction differences were at most2.6645352591003757e-15; NN numerical predictions matched exactly. Linux output CSV bytes differ from Windows serialization. All numeric output columns met the predeclared1e-12 tolerance. See `linux_replay_checks.json` and `LINUX_RUN_PROVENANCE.json`, and the [executed run](https://github.com/zhaohangznb-alt/holdout-design-validation/actions/runs/37128892159).

The original sealed archive is unchanged. The wrapper copies scripts into its new output directory and makes three explicit serialization adaptations: CRLF ordered-input signatures in both runners, and CRLF in reconstructed external input CSVs. Runtime and original hashes are logged in `COMPATIBILITY_ADAPTER.json`. Model algorithms, configurations, label values, row order and eligibility predicates are unchanged. `external_split_manifest.csv` adds the120 archived input hashes omitted from the sealed external-preparation directory. The source data and predictions remain in the original public bundle.

```sh
python unpack_bundle.py
python -m pip install -r bundle/requirements.txt
python -B supplements/representative_replay/linux_replay.py --bundle bundle --preflight-only
python -B supplements/representative_replay/linux_replay.py --bundle bundle --output bundle/linux_replay_v1
```

Use Python3.9.10 and pinned library versions. The Actions job uses a digest-pinned official Python3.9.10-bullseye container. Output must be new; diagnose a failed or unknown invocation before further model execution. Package installation uses network; input preparation and models use only local frozen data. No bootstrap or model tuning is run. Full680-case clean-directory retraining has not been performed on either platform.
