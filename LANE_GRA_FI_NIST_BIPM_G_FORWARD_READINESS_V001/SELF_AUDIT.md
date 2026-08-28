# Self-audit - NIST/BIPM public G-forward readiness

**Lane:** `GRA-FI-NIST-BIPM-GFR-V001`

**Date:** 2026-08-27

**Verdict:** `PASS_24_OF_24__PUBLIC_SUMMARY_CEILING_RETAINED__READY_FOR_MANIFEST_FREEZE`

## Load-bearing checks

1. The official committed PDF is the only empirical source; its hash and 31
   pages are checked before extraction.
2. All numerical fields have precise PDF-page, printed-page, table/equation
   custody.  Relevant pages were also visually rendered and inspected.
3. `A_j=16 Gamma_j m_s m_t/R_s` is the paper's finite torque-difference
   coefficient.  It is called a finite contrast, not the infinitesimal GC06
   derivative.
4. Nominal `k_g=0` is limited to the source-mass contribution at ideal torque
   extrema.  No finite-dither or total-apparatus zero is inferred.
5. Table-15 torque is treated as an already calibrated summary.  `C=1` in the
   reduced map does not claim unit gain for raw angle or voltage readout.
6. The eight Table-15 ratios, diagonal Type-A WLS, and four-result Table-18
   GLS are labelled diagnostics, not independent G estimates.
7. Table 16 owns the four displayed central values and whole-ppm uncertainty
   rounding; the tenth-ppm uncertainties used to construct the diagnostic
   covariance are explicitly owned by the Table-17 combined row/Table-18
   diagonal.
8. The paper's CODATA-centered Table-19/final-consensus hierarchy is excluded.
   No accepted, consensus, or synthetic G is an input or validation target.
9. Configuration-level physical remainders exactly alias the source column;
   the public paper does not justify replacing them with only two common mode
   offsets.
10. Table 18 is not miscalled the covariance of the eight torque rows, and
   Table 17 category summaries are not promoted to the missing calibration
   covariance.
11. Missing raw geometry, transfer, covariance, calibration, source-scale,
    and remainder fields are itemized.  No full GC16 fit, lineage test,
    Gravity Formation confirmation, or gravity derivation is claimed.

## Residual boundary

The source coefficient is usable for apparatus planning and for testing a
future raw-data adapter.  It is not a substitute for the authors' unreleased
mass-integration inputs or raw acquisition.  Any later release must be hashed
and mapped prospectively; it must not be back-filled from these summary ratios.

The deterministic analyzer passes `24/24`; its JSON mode reproduces
`RESULT.json` byte-for-byte and its zero-argument output reproduces
`VERIFICATION.txt` byte-for-byte.
