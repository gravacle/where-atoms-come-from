# T-51 NO-MEASUREMENT ADJUDICATION FRAMEWORK

**Status:** judge design only; no judgment, measurement, registration, or terminal disposition.

**Controlling text:** `FIELD_INSTRUMENT_V001.md` section 1, read verbatim.  The builder and the
quarantined verifier/judgment drafts are evidence about what has been attempted, not authority to
repair or reinterpret the registered rule after seeing the result.

## 1. Scoreability before any final judge sits

| Item | Verbatim registered clause | Present admissibility status | Reason and later scoring rule |
|---|---|---|---|
| V1 | “`|F|` exceeds the same-table measured control floor beyond contact at the larger separation.” | **PENDING VERIFICATION; builder-positive, not yet judged.** | The sealed builder reports `min |F| = 8.562e-08` at connected `w_enc=5` against `1.0e-11`.  A valid quarantine verifier must independently reproduce the far reading and floor, and must keep the source-absent equality a construction certificate rather than count it as a control.  If verified, this clause is Boolean and scoreable. |
| V2 | “the reading follows earned geometry under the placement swap.” | **UNSCOREABLE on the present lane, not FALSE.** | The registered reading is `F(D)`, but the builder silently makes onset order `k` the scored reading, introduces `TOL_SWAP=0.25`, and treats `F` magnitudes as unscored data.  No registered statistic says whether the swap compares `F`, `|F|`, a winding-adjusted reading, onset order, or a curve; no registered tolerance says what “follows” means.  Moreover, the commissioned far class `(d_gen=2,w_enc=5)` is a singleton, so its advertised swap was replaced by a different `(2,4)` class.  The observed ratios (`0.901` near, `1.238` in the substitute far-side pair) cannot be called a violation because no applicable comparator or threshold was registered.  Missing semantics produce `UNSCOREABLE`, never `FALSE`. |
| V3 | “the onset-order bracket contains **connected** `w_enc` at both placements (scored only after the connectivity gate).” | **PENDING VERIFICATION; builder-positive, not yet judged.** | The gate is the precondition.  The builder reports `2.977` in the bracket about `3` and `4.948` in the bracket about `5`; the draft rebuild is supportive but inadmissible as the final verifier.  A valid verifier must independently close connectivity/minimality and the off-grid fits before a judge may compute this Boolean. |
| V4 | “back-action below the declared tolerance.” | **DESCRIPTIVE/QUALIFIED ONLY on the present lane.** | The measured maxima (`1.191e-04` near, `2.492e-05` far) may be reported.  The `1e-03` threshold was not fixed in the field instrument; the lane audit says it was selected after a scratch survey of the same result scale and only before the sealed run.  That is not a prospective falsifier independent of the observed quantity.  Thus “below the lane's chosen threshold” is a true descriptive comparison, but it is not a clean pre-registered Boolean available to the terminal conjunction.  A fresh holdout with the tolerance and its scale fixed first can make V4 scoreable. |
| V5 | “sign attribution licensed only if the winding-sector sweep separates content sign from winding sign.” | **PENDING VERIFICATION; at most a scoped license.** | The builder reports a stable character `c=(0,0), s0=-1` near and `c=(1,1), s0=-1` far.  The far sign reverses under winding alone, so any eventual positive score licenses the sign only relative to the sealed `zbar1/zbar2` representatives; it does not establish an intrinsic, convention-free far sign.  A valid verifier must test the registered all-lambda stability condition, not one selected lambda. |

`UNSCOREABLE` is an admissibility finding.  It is not the negative value of a Boolean and must not
be converted to `FALSE` to force a registered branch.

## 2. Reachable branch

The field instrument registers two terminal routes: `ALL PASS`, and the `V1 FALSE` route (mediator
family sweep, followed by negative closure only if the null persists).  On the present record:

- `ALL PASS` cannot fire because V2 is unscoreable and V4 is only a qualified descriptive result.
- `V1 FALSE` cannot fire because V1 has not been scored false; the builder result is positive.
- A verifier refutation would make the affected item unscoreable; it would not by itself create the
  `V1 FALSE` branch.

Therefore the only legal present disposition is **PARTIAL POSITIVE / REPAIR REQUIRED / T-51
REMAINS OPEN**.  This is the commission brief's partial case (“name exactly what blocks scoring”),
not a new terminal branch.  It authorizes no C-94 row, no C-77 increment, no `T-51 DONE`, and no
claim that the record surface has passed the registered field-side rule.

## 3. Minimum admissibility packet from a quarantine verifier

A later no-measurement judge must refuse to sit until one sealed verifier packet supplies all of
the following:

1. **Frozen inputs and independence.** Hash every builder input before work; preserve code, raw
   output, seals, execution transcript, and D24 error/correction audit.  Treat all current
   `VERIFY_CODEX` and judgment material as quarantined input.  Import no builder helper or program
   geometry/sector implementation; hard-coded builder numbers may appear only in a comparison
   step performed after the independent result exists.
2. **Default-refuted A--E decisions with dependency mapping.** Each attack stays `REFUTED` until
   its own computed gates pass.  Map any refutation to the affected V-item(s); never score around a
   refutation.
3. **A, sector exactness.** Compute, rather than insert, every mediator/source and
   Hamiltonian/source commutator on 3x2 and the commissioned 3x3 constructions; audit actual
   source-changing transitions or an equivalent full-spectrum decomposition.  The quarantined
   draft's literal `mediator_pairings = [[0 ...]]` is not the commissioned proof-by-computation.
4. **B, connected `w_enc`.** Rebuild endpoint, enclosure, winding, and connectivity machinery from
   written definitions; exhaustively attack the claimed minima at both commissioned placements
   and at a fresh placement; exhibit witnesses and counts below the claimed minimum; test the
   excluded-plaquette and winding-representative conventions.
5. **C, winding attribution.** Recompute all four winding sectors at both placements for every
   commissioned lambda and fresh off-grid lambdas.  Apply the registered requirement that the
   same unique `(c,s0)` survive every usable lambda; carry the convention scope.  The quarantined
   draft computes `0.037` and `0.071` but its verdict predicate checks only `0.037`, only near/far,
   and not its fresh placement; that predicate is incomplete.
6. **D, onset.** Use an independently declared off-grid ladder, estimator, numerical floor, and
   perturbation-of-floor audit.  Report fits and local slopes for both commissioned placements,
   and compare only afterward with the builder and connected minima.
7. **E, D-1/directive audit.** Recursively enumerate every in-scope lane artifact existing at the
   verifier's input cut, report raw hits, apply an explicit human-readable adjudication rule, and
   include a positive scanner control.  No constant can decide the result.  The quarantined draft
   scans only top-level files and returns a literal `violation_count: 0` regardless of its hits;
   its E verdict is therefore inadmissible.
8. **V1 and V4 provenance audit.** Independently reproduce the far `F` values and numerical floor,
   distinguish certificates from controls, reproduce raw back-action, and state that the existing
   V4 threshold was trained on the same-result scale.  The verifier must not promote the qualified
   V4 comparison into a prospective Boolean.

The packet must end with itemized `NOT_REFUTED`/`REFUTED` outcomes and raw values, not an `ALL PASS`
statement.  Scoring remains the later judge's work.

## 4. Minimum admissibility packet from a prospective confirmatory lane

The old data may be used as design/training data, but the confirmatory holdout must not exist when
the following are sealed:

1. **A V2 operational supplement.** Without changing the quoted V2 words, name exactly which
   registered reading is compared (`F`, `|F|`, a predeclared winding-normalized value, onset, or a
   whole-curve statistic), the pairing rule, the statistic, the noise treatment, and the pass
   tolerance.  Choose a venue/placement rule that supplies at least two Gamma-equivalent members
   at the scored larger earned separation; no substitute separation after results.  Both branches
   must be reachable.
2. **A prospective V4 budget.** Freeze the tolerance and its independently owned scale before any
   holdout back-action is computed.  A separate calibration, a fixed fraction of a pre-existing
   gap, or the old T51_A data may train the budget, but none of the new holdout rows may choose it.
3. **A complete choice ledger.** Freeze venue, probe/source placement rule, mediator, lambda grids,
   winding conventions, floors, estimators, exclusions, missing-data rules, and how confirmatory
   results combine with V1/V3/V5.  Prefer re-running all V1--V5 on the same holdout so the terminal
   conjunction is coherent.
4. **Builder and independent verifier packets.** Seal raw tables, code, D24 audit, and manifests;
   then run a fresh default-refuted verifier with its own machinery and fresh off-grid points.  A
   final judge receives both packets and adds no measurement.
5. **Judge contract.** Quote V1--V5 verbatim; label refuted dependencies `UNSCOREABLE`; distinguish
   a scoped V5 license from an intrinsic sign; fire only a branch whose every registered
   precondition has been met.

## 5. Decisions the registrar cannot make for Brian

- The field text calls `F(D)` “the field reading” but leaves V2's “follows earned geometry”
  operator undefined.  Now that the results are visible, a registrar cannot select whichever of
  `F`, magnitude, winding-normalized response, or onset makes the rule pass.  If more than one
  scientifically live operational supplement remains, Brian must choose it prospectively (or
  explicitly delegate the choice) before the holdout is built.
- The same applies to an outcome-sensitive V4 budget: the registrar can preserve the qualified
  old comparison and can prepare prospective alternatives, but Brian must choose among materially
  different scientific budgets if the choice changes the program's terminal outcome.
- Section 0 says whether a reading available only through the priced mediator is the surface's
  field or the mediator's “is a judgment C-77 must eventually absorb” and “no control can decide
  it.”  The registrar may land only the already registered **mediated-response** scope after a
  valid `ALL PASS`; any unqualified ownership ruling or broadening of “field” beyond that scope is
  Brian's semantic/program decision.
- The field instrument contains no terminal meaning for the current partial case.  The registrar
  may keep T-51 open and commission repair, but cannot declare partial evidence equivalent to
  `ALL PASS`, close T-51, or invent a third terminal branch.  Any decision to change that boundary
  must come from Brian.

The principal's existing accumulation directive does not require a new decision: if a valid
`ALL PASS` later fires, the two-source composition computation is next and every computed
composition law remains registrable.
