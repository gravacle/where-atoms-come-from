# D24 audit — independent T53-A VSM verification

- Verdict: `NOT_REFUTED`
- Verdict scope: `raw-integrity-and-measurement-only`
- Default: `REFUTED` until every declared predicate executes and succeeds.
- Predicates executed: 21 of 21
- Predicates succeeded: 21 of 21

## Errors and predicate failures

None observed in this bounded raw-integrity and measurement recomputation.

## Scope audit

The non-refutation verdict, if reached, applies only to exact raw-file integrity, protocol semantics exposed by the raw headers, and the independently recomputed numerical measurements. It is not a record-formation law, a universal result, a gravity result, or a program-status decision.

The dataset is retrospective and restricted to five labeled magnetic sediment specimens. It contains no randomized intended-message assignment, no blind prediction, no no-write control cohort, no long common-hold survival study, no independent laboratory reproduction, no cross-surface coverage, and no gravity observable.

## Reproducibility

Run `PYTHONDONTWRITEBYTECODE=1 python3 -B verify.py` from this directory. The verifier uses only the Python standard library and writes deterministic result artifacts with no clock or machine-specific path fields.
