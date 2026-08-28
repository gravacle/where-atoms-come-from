# Self-audit: HUST-2018 dual-method G-forward lane

## Result

`ACCEPT_WITH_PROCESSED_COEFFICIENT_AND_FIGURE_LEVEL_CEILINGS`

## Adversarial checks

1. **Was an accepted value of `G` used?** No.  MOESM4, which contains
   historical/accepted-value comparison data, is not pinned or read by the
   analyzer.  Published HUST row outputs are comparison targets only.
2. **Are the 7 ToS and 29 AAF rows called raw?** No.  They are labelled
   already-derived `G` summaries everywhere.
3. **Are the ToS period points called event-level data?** No.  The Figure-2
   caption identifies each as a three-day extracted period.  The unavailable
   0.5-second angle data remain explicitly missing.
4. **Is A-B-A uncertainty understated?** No uncertainty is assigned to the 18
   overlapping triples.  The comparison uses only the authors' printed row
   uncertainty.
5. **Was polynomial order chosen to force agreement?** No theorem rests on the
   polynomial result.  The method-faithful A-B-A contrast is primary; the
   common-quadratic result is independently labelled a drift diagnostic.
6. **Is the 7,200-second AAF segment promoted to a campaign average?** No.  It
   is used only to demonstrate source/background spectral separation.
7. **Is the cross-method 2.72 ratio called a discovery?** No.  It is explicitly
   conditional on zero cross covariance and retained only as a source-model
   stress diagnostic.
8. **Are processed `Delta C_g/I` and `P_g,l,m` called raw geometry?** No.  The
   missing finite mass-coordinate-density files are the leading blocker.
9. **Are missing corrections or remainders set to zero?** No.  The ToS bracket
   residual and the complete correction/remainder ledger remain unowned.
10. **Is ordinary-gravity agreement claimed to confirm record lineage?** No.
    No RGRL/GFT confirmation or lineage charge is inferred.
11. **Is `GC16` said to be complete?** No.  The exact missing geometry,
    transfer, covariance, stress, nuisance, and holdout fields are enumerated.
12. **Was existing sealed or canonical work modified?** No.  Every change is
    confined to this new lane directory.

## Scientific significance retained

The strongest result is not another average of published `G` values.  It is
the real-data separation of the two finite-apparatus ownership channels: ToS
measures gravity in the dressed stiffness/operator, and AAF measures gravity
in the source/forcing response.  The official release reproduces both at the
processed-coefficient level without double counting.
