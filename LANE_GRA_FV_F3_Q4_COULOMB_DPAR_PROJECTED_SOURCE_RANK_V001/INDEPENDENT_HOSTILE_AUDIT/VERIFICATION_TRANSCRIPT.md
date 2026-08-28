# FV independent hostile-audit verification transcript

**Date:** 2026-08-27  
**Audited lane:** `LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001`

## Repaired builder replay

Command:

```text
python3 verify_projected_source_rank.py
```

Terminal summary:

```text
SUMMARY 89/89 exact checks passed
PATHS 720; PREFIX_DERIVATIVES 3600; J6 63/8
COULOMB_RING_ROWS A1+T2 rank4; DIRECT_NONIDENTITY E rank2
PROJECTED_OPERATOR_RANK 6; WITNESS_DETERMINANT -4678629417/256
FORMAL_THROUGH_H8 rank6; CTP_WARD_TENSOR_POLE_GRAVITY not claimed
```

The repaired core added and executable-checked the explicit
`S10 / FV-PURE` premise, removed the literal carriage-return byte, rebuilt
the core manifest, and resealed the packet.

## Independent scientific replay

Command:

```text
python3 INDEPENDENT_HOSTILE_AUDIT/audit_projected_source_rank.py
```

Terminal summary:

```text
SUMMARY 83/83 independent hostile checks passed
HISTORIES 720; PROPER_PREFIXES 3600; ENVIRONMENTS 4x64
DIRECT E rank2; H6_RING A1+T2 rank4; OPERATOR rank6
WITNESS_DETERMINANT -4678629417/256; FORMAL h^24 + O(h^26)
DISPOSITION PASS_AFTER_FV_PURE_PREMISE_AND_BYTE_HYGIENE_REPAIR
CEILING off-shell formal rank only; no CTP/Ward/pole/gravity/G
```

The independent executable does not import or invoke the builder verifier.
It reconstructs the exact rational perturbation theory, `G_5`, global ice
completions, operator functionals, determinant, slopes, endpoint topology,
and formal-order argument independently.

## Repaired core hashes

```text
22a22985c96003abc319accc7fe74ee0a89f401b9816fecc8caa0f7b73e62b20  DEPENDENCIES.sha256
651a66b9afd7545b04aa80e5f90952fda9327d011ecd19b973aa80ff51a739f3  MANIFEST.sha256
9f42b8fa93ed03184c858b296f185a41d008d0acbb17166f28bf0e9512611948  README.md
b5d4c3de99aa4e100519c19a9b74de487b47c1a2d3671204e77740bd9094771a  RESULT.md
8301da6bbc026d0e14d985592c5dabe3d91072957c4aaa4b1bebf1f45aadd894  SEAL.sha256
431514695cc31364ed5019ca1b6f99a1cd770f7d9bb4364e59d08e4fa0b84354  SELF_AUDIT.md
6fc221a31151340b91a946d33e442971c1373500e067c354b6c610e3964edb1c  THEOREM.md
4c2195bdaef6ac570bf2c0a8451035bf375838d7651d39081bd2bbbb0d98dbfd  VERIFICATION.txt
0e93d84f9eb7cf7fdd62b5a14d5c6705c74841899dd1676bdd7e7a41eb971a00  verify_projected_source_rank.py
```

## Independent payload hashes before sealing

```text
32c4b0f61b053f12331d016562532989c2007be8769f75b4e4c05c55a10fb055  CORE_CUSTODY.sha256
3801fd9ba6ba3c0fe80c9f4792abfdeb6dd7c37c7145663be05b4d56f8160723  INDEPENDENT_HOSTILE_AUDIT.md
4639a642bd27f20beb4bba5d55af1208e92263f523841847c86dd909ff79dd6c  audit_projected_source_rank.py
3a6f780f63a9baa5eaf269702fc09c32eb58af01624f41e53a34c78ad44ddfec  verify_independent_audit.py
```

Only files inside `INDEPENDENT_HOSTILE_AUDIT/` were created or edited by the
auditor.  Builder core, canonical files, and git were not edited by the
auditor.
