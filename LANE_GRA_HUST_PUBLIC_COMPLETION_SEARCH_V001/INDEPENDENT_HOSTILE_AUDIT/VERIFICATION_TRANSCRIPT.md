# HUST public-completion independent hostile-audit transcript

**Date:** 2026-08-28  
**Audited lane:** `LANE_GRA_HUST_PUBLIC_COMPLETION_SEARCH_V001`

## Repaired builder replay

Command:

```text
python3 verify_hust_public_completion_search.py
```

Terminal summary:

```text
SUMMARY 48/48 checks passed
search evidence classes partition all entries
second lead is unverified title only
seven-object inventory dependency is explicitly pinned
world exhaustiveness denied
```

The repaired core separates seven literal endpoint/page entries from four
curator-recorded nonexecutable aggregate entries, pins the actual seven-object
prior custody inventory, and enforces the distinct evidentiary types of the two
acquisition leads.

## Independent endpoint and custody replay

The audit re-queried Nature, Crossref, DataCite, both exact Figshare POST
endpoints, Zenodo and the two literal HUST pages on 2026-08-28.  It also
streamed all seven official Nature release objects without retaining or
promoting the binaries.  MOESM1 through MOESM7 reproduced the seven SHA-256
values in the prior HUST custody inventory.

Observed endpoint summary:

```text
Nature associated objects: 7 (MOESM1--MOESM7)
Crossref relation object: {}
DataCite exact related-identifier records: 3 later References
Figshare exact resource-DOI matches: 0 articles; 0 collections
Zenodo exact DOI-text records: 2 later unrelated publications
New qualifying completion roots: 0
```

Mutable public responses were normalized into
`LIVE_REQUERY_NORMALIZED.json`.  The receipt is query-date evidence, not a
timeless or exhaustive index certificate.

## Independent offline replay

Command:

```text
python3 INDEPENDENT_HOSTILE_AUDIT/audit_public_completion_search.py
```

Terminal summary:

```text
SURFACE_COUNTS literal=7 curator=4 total=11
NATURE_OBJECTS 7
DATACITE_RECORDS 3
ZENODO_RECORDS 2
FIGSHARE_MATCHES 0
VERDICT PASS_AFTER_COMPLETENESS_LEAD_TYPING_AND_CUSTODY_REPAIR
SUMMARY 103/103 independent checks passed
```

The independent executable does not import or invoke the builder verifier.  It
validates the repaired core, dependency hashes, prior seven-object inventory,
independently observed endpoint receipt, dissertation-lead typing, boundedness
ceiling, no-binary/no-accepted-`G` boundary, builder manifest and seal, and byte
hygiene independently.

## Repaired core hashes

```text
ee1e2b795206c56875da176153410506d66f7984afc84469f7430f4aa677c4e1  README.md
b763e8d313aa034368f297ffd49bbad016209f4ad70f84a33082c3b049ec733a  THEOREM.md
e55497a1771af7a6e69501736f3f571e34050ee2a4fb480ca711e182ec8698ef  SEARCH_LEDGER.json
ada9043be0a6fb36614cdcc5589009424aac04c9bb6161a6c7eb6b5c7072c0dd  RESULT.json
8307e9ad91513a927bf61a232caf7275c0694e504de7771d9d65b14f9fb803b4  SOURCE_CUSTODY.json
48b82bb582798582fc45cbad99c8e386fb796425d87f0ed8510894511487f941  verify_hust_public_completion_search.py
013a99d6cec089d56b5047422f750312f70ac9d66c0577ac08a7eb6dad46fe86  VERIFICATION.txt
78a0efe8a7ea5f67986a7275f654915877fd7fe0c545b3d69c2ffa66053ead9f  MANIFEST.sha256
2a0f230deb2d87b251fc36d5d97f28bd3462cb3f55415f2ea23540344b74cbbc  LANE_SEAL.sha256
```

## Independent payload hashes before sealing

```text
ca70eef54f45069053f55a3b09f67a207877edce1847048dec4803559e7f9238  CORE_CUSTODY.sha256
693eacff5a7328d43a3b94ccff27cc54f739a9af817c1c9b59eea940a678ad12  LIVE_REQUERY_NORMALIZED.json
faf1b012878b3b4c7daa335a9a9520e23eaa7cc3688f53d7e2a950a7c971c02c  INDEPENDENT_HOSTILE_AUDIT.md
29bda10929ef75e8fc3f593950db063fb6e2323ddb49fe31c3fdb4802164f34a  audit_public_completion_search.py
950348f35738119d8f6c5b394eefacdeb1dd6afff67e3b733702f9ff54f1e00f  verify_independent_audit.py
```

Only files inside `INDEPENDENT_HOSTILE_AUDIT/` were created or edited by the
auditor.  Builder core, canonical files and git were not edited by the auditor.
