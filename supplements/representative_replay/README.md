# Representative replay

The prior Windows replay rebuilt408 partition input tables byte-for-byte and executed RF seeds42/53/67 and deterministic Tanimoto1-NN for ALK and ACHE, main/revised/source_exact. ALK has239 test rows and ACHE921. Six fresh RF fits matched archived predictions within2.6645352591003757e-15; both NN prediction CSVs were byte-identical. Previous two Ridge checks are separate. Full680-case clean-directory retraining has not been performed.

`windows_replay_checks.json` records this executed check. The frozen input tables and archived predictions needed for replay are in the existing sealed bundle. The Linux wrapper uses the unchanged bundled runners and explicitly limits execution to these same eight cases.

```sh
python unpack_bundle.py
python -m pip install -r bundle/requirements.txt
python -B supplements/representative_replay/linux_replay.py --bundle bundle --preflight-only
python -B supplements/representative_replay/linux_replay.py --bundle bundle --output bundle/linux_replay_v1
```

Use Python3.9.10 and the pinned library versions. Preflight verifies Python's numeric version, library versions and archived reference cases; the actual platform/build string is recorded independently, allowing genuine Linux verification rather than comparing it to a Windows build string. Output must be new, under the extracted bundle; a failed or unknown invocation must be diagnosed before further execution. Models/preparation use frozen local data without network requests. No bootstrap or model tuning is run. The Actions workflow records actual Linux results and logs; its presence alone is not a passing Linux result.
