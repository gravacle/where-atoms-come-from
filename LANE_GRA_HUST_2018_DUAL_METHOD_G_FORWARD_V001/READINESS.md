# HUST-2018 public-data readiness against the finite-apparatus G protocol

## Verdict

`PROCESSED_DUAL_CHANNEL_FORWARD_CLOSED__FULL_GC16_NOT_READY`

The official release goes materially beyond a paper-only readiness audit: one
ToS stiffness response can be recovered from Figure-2 period summaries, the
AAF source-response equation can be reproduced for all three campaigns, and a
representative AAF acceleration stream resolves the driven and lab-background
harmonics.  The release still lacks the geometry, complete transfer,
remainder, and covariance packet required for a new or full-apparatus `G`
cross-check.

## Scored official workbook ranges

| Range | Physical content | Classification | Use in this lane |
|---|---|---|---|
| `a!B3:D22` | 10 near + 10 far periods with source masses present | three-day period summaries; raw-like intermediate | A-B-A and common-quadratic `Delta omega^2` |
| `b!B3:D22` | 10 near + 10 far periods without source masses | three-day background summaries; raw-like intermediate | background subtraction |
| `c!B3:D9` | seven ToS `G` outcomes | already derived | custody and comparison only |
| `c!E10:F10` | combined ToS `G` and uncertainty | already derived | cross-method stress only |
| `d!B4:G129604` | residual/free twist and thermal-noise PSD arrays | frequency-domain intermediate | inventoried, not needed for the bounded forward |
| `e!B3:C10001` | 9,999 one-second angular-acceleration samples | raw-like processed stream | source/background spectral diagnostic |
| `e!B3:C7202` | first 7,200 seconds described by the Figure-2 caption | raw-like processed stream | scored two-tone fit |
| `f!B3:C6` | 4 AAF-I `G` outcomes | already derived | campaign mean custody |
| `f!D9:E18` | 10 AAF-II `G` outcomes | already derived | campaign mean custody |
| `f!F21:G35` | 15 AAF-III `G` outcomes | already derived | campaign mean custody |
| `f!H39:I39` | combined AAF `G` and uncertainty | already derived | cross-method stress only |

## Published-PDF fields used

| PDF custody | Field | GC ownership |
|---|---|---|
| pages 4-6, printed SI 3-5 | exact ToS and AAF equations | ownership rule |
| page 19, Supplementary Table 1 | `I_m`, `K_m`, `I`, `K` | mechanical correction/transfer |
| page 20, Supplementary Table 2 | seven `Delta C_g/I`, `Delta omega^2`, and derived `G` rows | processed stiffness source/response |
| page 21, Supplementary Table 3 | three `sum P_g,l,2`, `alpha_t`, and derived `G` rows | processed forcing source/response |

## GC protocol matrix

| GC requirement | Public status | Consequence |
|---|---|---|
| finite source/detector coordinate-density files | absent | `a` and `k_g` cannot be independently recomputed |
| raw observations with run/configuration labels | limited | one figure-level ToS and AAF diagnostic, not a full fit |
| independent global source scale and covariance | absent | `G` cannot be separated anew from all source calibration ownership |
| complete torsion/support/auxiliary transfer | absent | dressed operator cannot be built row by row |
| readout/controller calibration | absent | pre/post-response ownership cannot be re-audited |
| signed correction/remainder ledger | absent | missing effects cannot be set to zero |
| observation and calibration covariance | absent | no authoritative likelihood or coverage statement |
| complete conserved apparatus stress ledger | absent | full finite-apparatus source premise remains open |
| predeclared null/holdout/nuisance packet | absent | no prospective GC protocol execution |

## Claim ceiling

This lane establishes a real-data **processed-coefficient forward** and a
limited **figure-level response extraction**.  It does not establish a new
value of `G`, an independent finite-geometry calculation, a full `GC16`
execution, RGRL/GFT confirmation, or a record-lineage gravitational charge.
