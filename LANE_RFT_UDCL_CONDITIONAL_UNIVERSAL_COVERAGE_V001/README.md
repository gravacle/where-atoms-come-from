# U-DCL conditional universal Coverage-U lane

**Version:** V002 formal clarification of the adopted V001 postulate

**Date:** 2026-08-25

**Status:** independently hostile-audited exact conditional theorem; natural
validity of U-DCL remains unproved and empirically open

## Result first

The program adopted U-DCL as its working universal physical-domain postulate.
This lane now gives the exact typed consequence:

\[
 \operatorname{REC}(r)\land DCL_{\rm phys}(r)
 \Longrightarrow \operatorname{COV}_{\cup}(r),
\]

and therefore

\[
 \operatorname{U\!DCL}\Longrightarrow
 \forall r\in\mathfrak R^{\rm actual,bf}_{\rm FM},
 \operatorname{COV}_{\cup}(r).
\]

The implication is a theorem inside the adopted working axiom system.  U-DCL
itself is still an unproved, falsifiable claim about nature.

## Load-bearing repair

The first author version passed Boolean composition checks but left its local
physical witness under-typed.  The independent hostile audit required one
bounded repair, recorded in `TYPED_CLARIFICATION_V002.md`:

1. `DCL_phys(r)` is an ontic physical predicate; `Cert_DCL(r;P)` is its
   prospective evidential certificate.
2. One common branch `b in {K,W}`, frontier, incidence, census, history state,
   maps, instruments, and future-input object must satisfy D1--D4 together.
3. The external macro-incidence is physically bounded and acyclic after only
   already-well-posed internal objects are contracted to exact K/W maps.
4. D2 supplies one joint arm/context state containing the history register and
   its distribution; no joint object is inferred from history conditionals or
   component marginals.
5. The two-sided maps are well posed on that same state, and all later fresh
   inputs have one arm-common **joint** conditional rule.

These restrictions make the claimed bridge to sealed `C`, `S`, and `J` exact.
They do not add a representation category or weaken a falsifier.

## Imported sealed result

The lane imports, without modifying it,
`LANE_RFT_STANDARD_CAUSAL_URFT_SCOPE_V001/THEOREM.md`:

```text
REC(r) & C(r) & S(r) & J(r)
  -> CTS(r)
  -> OCC_union(r, C_r)
  -> faithful finite-mission encoding Xi with ARCH_union and A1--A4.
```

All external dependencies and their sealed manifests are pinned in
`DEPENDENCIES.sha256`.

## Files

- `THEOREM.md` -- typed predicate, local bridge, per-record theorem, and global
  universal generalization.
- `TYPED_CLARIFICATION_V002.md` -- exact relationship between adopted V001
  prose and the theorem predicate.
- `BOUNDARY.md` -- independent domain and exhaustive typed falsifier burden.
- `RESULT.md` -- earned result and claim ceiling.
- `AUDIT.md` -- repaired author-side scope ledger.
- `INDEPENDENT_AUDIT.md` -- hostile findings and post-repair disposition.
- `DEPENDENCIES.sha256` -- pinned adoption, decision/audit, and sealed theorem
  dependencies.
- `verify_udcl_conditional_coverage.py` -- standard-library custody, text, and
  finite-logic regression.
- `VERIFICATION.txt` -- fresh execution transcript.
- `MANIFEST.sha256` -- final lane content seal.

## Reproduction

From the repository root, run:

```bash
python3 LANE_RFT_UDCL_CONDITIONAL_UNIVERSAL_COVERAGE_V001/verify_udcl_conditional_coverage.py
```

The verifier checks the declared formal composition and artifact custody.  It
does not test whether nature obeys U-DCL.

