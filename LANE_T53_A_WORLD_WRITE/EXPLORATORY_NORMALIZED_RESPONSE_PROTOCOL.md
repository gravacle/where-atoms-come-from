# Exploratory normalized writer-response protocol

Status: frozen before executing `analyze_normalized_response.py`, but after the source
data and the primary writer-off results were already observed.  This is retrospective
physics characterization, not a confirmatory test.

## Physical question

Do five real sediment specimens measured by the same DCD protocol share a useful
dimensionless writer-off response when each curve is expressed in terms of

- normalized reverse writer pulse `x = |H_pulse| / |H_zero|`, where `H_zero` is the
  linearly interpolated field at which that specimen's retained moment changes sign; and
- normalized retained state `r = m_after / m_initial`?

This tests whether initial retained moment and the sign-change field absorb most
within-mechanism specimen variation.  It does not test record creation, a universal law,
or a second mechanism.

## Frozen calculation

1. Use only the five Step 2 remanence curves already pinned in `world_observation.json`.
2. Require writer-off acquisition, monotone time, a nonzero initial moment, and exactly
   one strict sign change in retained moment.
3. Compute `(x,r)` directly from every raw row.  Do not fit or discard rows.
4. Linearly interpolate each specimen on the fixed grid
   `x = 0.00, 0.05, ..., 1.50`; no extrapolation is allowed.
5. At each grid point report all five values, their mean, population standard deviation,
   and range.  Report each specimen's `H_zero`, initial moment, and row count.
6. Report the root-mean-square deviation of all specimen/grid values from the gridwise
   mean and the largest population standard deviation.  No post-hoc PASS tolerance is
   assigned.

## Falsifiers and limits

A missing grid point, non-monotone normalized field, multiple retained-moment zero
crossings, or fewer than five curves invalidates the calculation.  Large dispersion is a
physical result, not an adapter failure.  Even close collapse would be evidence only for
this one mechanism and acquisition protocol; it cannot establish formation, necessity,
sufficiency, universality, gravity, or independent reproduction.
