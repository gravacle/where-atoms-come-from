# Independent hostile-audit verification transcript

**Date:** 2026-08-27  
**Audited lane:** `LANE_GRA_HUST_PUBLIC_CALIBRATED_SOURCE_IDENTIFIABILITY_V001`

## Fresh source reacquisition

Fresh HTTPS downloads to `/private/tmp` gave:

```text
5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb  official supplement, 2,711,453 bytes
23436d4be7600a7a9dffa02cc4167a20b6eea032a181e77899bb57bb90aa02e9  official Nature Table 1, 193,307 bytes
40756ec0fb8f00c1fde31020b294521a3b220a196bef884a2ea5f3534d77dfaa  public university primary-article mirror, 7,999,835 bytes
```

All three exactly match the frozen source objects.

## Independent audit replay

Command:

```text
python3 audit_calibrated_source.py --check
```

Output:

```text
PASS: AUDIT_RESULT.json matches independent recomputation
```

The executable independently transcribes/checks 80 correction fields, parses
42 official HTML error-budget vectors, reconstructs all ten calibrated partial
forwards and comparator gaps, derives AAF sinc corrections, checks ToS signed
dynamic factors and correction-composition scale, and reconstructs the full
category covariance hierarchy.

## Audit verifier

Command:

```text
python3 verify_independent_audit.py
```

Output:

```text
PASS: 303/303 independent-audit checks passed
```

The verifier additionally validates the repaired frozen builder manifest/seal,
runs the builder verifier (`PASS 192/192`), exact-replays the frozen builder
`RESULT.json`, checks all audit manifests/seal bindings, and executes fourteen
M1/N1/interval repair checks across `THEOREM.md`, `RESULT.json`, and the
analyzer. It confirms the ownership, uncertainty-domain, raw-response, and
claim-ceiling language is intrinsic to the repaired core.

## Syntax and text integrity

```text
PASS: Python compile, JSON parse, and text control-byte scan
```

## Final hashes

```text
7491007370ab60e23c59732c747889e0c87e76b4f115e2c1fb15a69a8e7f6f60  AUDIT_RESULT.json
95693c55b483243034358d29dfe6dbe0b4ac4a4e96a906a0c0e72c4b1db68232  INDEPENDENT_HOSTILE_AUDIT.md
6913aa593990cebfad98a090d9c6fc8a01647b5a39291054698791f58b46ff30  AUDIT_MANIFEST.sha256
0478e1f0ae331ba2bd1f227aaae460a1215b596f90a51100d75e6d7e92c0a634  AUDIT_SEAL.json
c237d48e7e0f4d750baa0a7e5d0a8f66bf2ecef70d47e92f6751d8661be5993b  audit_calibrated_source.py
01ae5a212886ee11387d1bdecac5cd76a4cd8ab1775f6ee832f9e941d6c01db1  verify_independent_audit.py
253c313590ed22b76ddfe25bcf50919dd7819d736274bf1b601db4f2de513764  INDEPENDENT_SOURCE_TRANSCRIPTION.json
a9db20e78353727deb394a3386a77ddd232cffae3d97288785f186407bca5277  SOURCE_REACQUISITION.json
```

Only independent hostile-audit artifacts were edited. The repaired builder
core, canonical files, and git were not edited by this audit.
