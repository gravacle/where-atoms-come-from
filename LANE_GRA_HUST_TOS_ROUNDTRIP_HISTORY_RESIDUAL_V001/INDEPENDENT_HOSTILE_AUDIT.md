# Independent hostile audit: HUST ToS round-trip history residual

**Audit date:** 2026-08-27  
**Verdict:** `ACCEPT_AFTER_TWO_SCOPED_CUSTODY_CLARIFICATIONS`

## Independent source replay

The audit did not import the analyzer or first verifier and did not accept
`RESULT.json` as a numerical premise. It reopened the pinned official Figure-2
OOXML workbook with SHA-256
`331ea6ee8a6a2558d4f6a8ffa4cd52e4c5d52d47f2e2a697889d0d3c8a3ad27a`.
A separate read-only import through the configured bundled spreadsheet runtime
showed the same `a!B2:D22` and `b!B2:D22` values and no formulas in either
scored range. The independent OOXML replay recovered ten near and ten far
period summaries per panel.

After merging the two blocks by workbook time, all nine `N-F-N` and all nine
`F-N-F` triples per panel satisfy the claimed chronology. Thus those labels
correctly describe endpoint orientation. They do not identify a causal arrow,
record state, or independently randomized history condition. The worksheet
headers encode time, near, and far; the source-present/source-absent meaning is
inherited from the already-audited article-level custody.

## Independent mathematical and numerical replay

The audit reconstructed the `9 x 10` first-difference matrix without using the
lane implementation. Its four-block return map has exact rational rank 36; its
two-block present-minus-background map has rank 18. Direct multiplication gives
the reported Gram blocks: diagonal/off-diagonal `2/-1` for panel returns and
`4/-2` for panel differentials. These are exact Euclidean weight-overlap
matrices, not empirical covariance matrices.

All 36 panel returns, all 18 present-minus-background values, the nine
orientation-even values, the nine orientation-odd values, every summary field,
and every printed theorem mean were reproduced. In particular,

\[
 \bar d_N=-1.8506327675395133\times10^{-10}\ \mathrm{s}^{-2},\qquad
 \bar d_F=-1.8282181693912770\times10^{-10}\ \mathrm{s}^{-2},
\]

\[
 \bar c=-1.8394254684653953\times10^{-10}\ \mathrm{s}^{-2},\qquad
 \bar h=-1.1207299074118054\times10^{-12}\ \mathrm{s}^{-2}.
\]

Those are reproducible arithmetic outputs, not calibrated significant digits.
The workbook periods lie on a nominal `1e-5 s` grid, while acquisition and
rounding covariance are unavailable. The replay also recovered the largest
same-ordinal present/background duration differences, `0.4445027 day` for
`N-F-N` and `0.0421588 day` for `F-N-F`; the residual is not time-normalized.

Each nine-return sum telescopes to the final minus initial transformed endpoint.
The two differential sums telescope to differences of those panel endpoints.
Accordingly, the mean is not nine independent replications of a history effect.

## Repairs made

Two custody ambiguities were material enough to clarify without changing any
number or enlarging the claim:

1. `d_i`, `c_i`, and `h_i` pair source-present and source-absent loops by the
   same **sequence ordinal** in separate panels. This is not simultaneity,
   randomization, or an authenticated matched-run design. The differential and
   orientation means are unchanged by re-pairing, but individual loopwise
   components and the reported odd-component RMS depend on that ordinal pairing.
2. Long decimal outputs are now explicitly typed as replay digits computed from
   source cells on a nominal `1e-5 s` period grid. They are not measurement
   precision and do not supply an uncertainty model.

The duration language was correspondingly narrowed from an unqualified “paired”
comparison to a same-ordinal bookkeeping difference. No numerical result,
weight, rank, chronology, or scientific conclusion changed.

## Decisive counterfactual failure

Every observed endpoint return contains an intervening opposite-position
excursion. The release supplies no matched equal-duration **no-excursion**
trajectory. Therefore the dataset cannot distinguish an excursion-dependent
hysteresis or memory effect from ordinary time drift, fibre ageing, mechanics,
thermal/controller state, panel offset, source-motion history, or other
unmeasured apparatus evolution. Source-present minus source-absent is also a
cross-panel contrast, not a randomized same-run sham.

That missing counterfactual is decisive. Neither the common negative component
nor the much smaller orientation-mean difference is statistical evidence for
memory, and neither authorizes a zero conclusion. Row covariance, event-level
angles, and an uncertainty model are absent.

## Hostile non-promotion screen

This packet is an exact extraction plus a reproducible history-confound
diagnostic. It does not:

- authenticate record formation, retention, KEEP/BREAK, or lineage;
- estimate `beta_TM` or a hidden gravitational charge;
- identify hysteresis, memory, or causal source-history dependence;
- provide a standard error, confidence interval, significance test, or
  empirical coverage;
- confirm RGRL, Gravity Formation Theory, or gravity emergence;
- execute full `GC16`; or
- determine a new value of \(G\).

The strongest surviving statement is exactly the narrow one: the official
figure-level summaries support a reproducible equal-configuration endpoint
return ledger, with exact endpoint-reuse algebra and an unresolved history
confound.
