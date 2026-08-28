# Stringent self-audit: HUST ToS round-trip history residual

## Source and spreadsheet handling

1. The analyzer reads the parent lane's pinned official Figure-2 workbook and
   verifies its SHA-256 before parsing. No workbook byte is copied or edited.
2. The workbook was independently inspected read-only with the configured
   bundled `@oai/artifact-tool`; sheets `a` and `b`, ranges `B2:D22`, were also
   rendered and visually checked. The observed cells match the analyzer's
   explicit OOXML extraction.
3. Only `a!B3:D22` and `b!B3:D22` are scored. Other HUST tables, derived `G`
   rows, accepted constants, and control workbooks are not inputs.

## Observable typing

4. The linear matrices act on transformed endpoint values
   \(y=(2\pi/T)^2\), not directly on the printed periods. The period transform is
   nonlinear; the endpoint-return operation after transformation is exactly
   linear.
5. The opposite-configuration middle observation certifies the schedule but is
   deliberately not used numerically. This is an equal-endpoint return, not the
   earlier A-B-A middle-versus-interpolated-endpoint gravity contrast.
6. `N-F-N` and `F-N-F` are kept separate before the even/odd decomposition.
   Neither parity component is given a causal name in the theorem.
7. Source-present minus source-absent is a cross-panel difference. The two
   panels are not a randomized same-run sham pair and cannot be called one.
8. Every endpoint return follows an opposite-configuration excursion, but
   there is no matched no-excursion counterfactual. Therefore no excursion,
   hysteresis, or memory effect is identified against ordinary time drift; the
   observable is a history-confound diagnostic only.
9. Present/background values with the same index are paired by sequence ordinal
   only. The panels are separate acquisitions; no simultaneity, randomization,
   matched-run authentication, or causal pairing is inferred. Component means
   survive re-pairing, while individual even/odd values and the odd RMS do not.

## Dependence and statistical ceilings

10. Every weight is frozen explicitly. The exact Gram matrices expose adjacent
   endpoint reuse, but identity endpoint covariance is only the algebra used to
   display overlap. It is not a model of the experiment's actual covariance.
11. Each orientation mean telescopes to the two outer endpoints. Treating nine
   overlapping returns as nine independent replicates would be a material error.
12. No row uncertainty or covariance is released. The packet therefore reports
    no standard error, \(p\)-value, confidence interval, significance, or
    equivalence/null conclusion.
13. Printed periods have finite decimal precision whose acquisition/rounding
    covariance is not supplied. It is not propagated as if the last digit were
    independent uniform noise. Long replay decimals are computational digits,
    not claimed measurement precision.
14. Loop durations are near six days but not identical across panels. The
    largest same-ordinal `N-F-N` difference is `0.4445027 day`; no hidden rate
    normalization or interpolation is applied, and “same ordinal” is not a
    matched-trial assertion.

## Physical and theory ceilings

15. A common negative panel-differential mean in both orientations does not
    establish memory, drift, or any other mechanism. It can contain ordinary
    fibre ageing, mechanics, thermal/controller history, source-motion effects,
    panel offsets, and duration mismatch.
16. The much smaller mean orientation difference is not evidence of zero
    orientation dependence. Its descriptive spread and missing covariance
    forbid that promotion.
17. Source motion is not authenticated record formation. No `M`, `L_T`, `L_D`,
    KEEP/BREAK, formation/sham, or lifecycle label is manufactured from source
    position or time order.
18. The result cannot estimate `beta_TM`, a lineage charge, or a record-only
    gravitational response. It contributes no empirical confirmation of RGRL,
    GFT, or gravity emergence.
19. No accepted value of \(G\) enters, but that fact does not turn a period-return
    diagnostic into a new `G` measurement or a full finite-apparatus fit.
20. This analysis was nominated after the parent HUST source had already been
    inspected. It is not blinded or prospective experimental evidence.

## Program hygiene

21. The lane adds no record definition, force law, interaction, metric field,
    or rescue term. It asks one bounded question of one existing public source.
22. No canonical MODEL, URM, gravity certificate, or experiment register is
    modified by this lane. Promotion, if any, requires a separate independent
    hostile audit.
