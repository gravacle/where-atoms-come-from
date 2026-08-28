# Independent hostile audit - NIST/BIPM public G-forward readiness

**Lane:** `GRA-FI-NIST-BIPM-GFR-V001`

**Date:** 2026-08-27

## Verdict

`ACCEPT_AFTER_CUSTODY_LABEL_AND_VERIFICATION_REPAIRS__REDUCED_FORWARD_ONLY`

The reduced source formula, units, all eight torque-mode observations, source
columns, rank/identifiability results, and covariance diagnostics survive
independent hostile recomputation.  The packet remains a public-summary
readiness result, not a real-apparatus execution of GC16 and not an independent
measurement of `G`.

Three custody defects were repaired without changing a scientific result:

1. Table 15 literally has four configuration rows containing eight
   torque-mode observations.  The prose now distinguishes the source-table
   layout from the eight-row analysis representation.
2. The four central `G` values come from Table 16, whose displayed uncertainty
   column is rounded to whole ppm.  The tenth-ppm values used to build the
   diagnostic covariance come from the Table-17 combined row/Table-18
   diagonal.  The result schema now owns those sources separately.
3. Exact source tokens were added for all Table-15 observations, all four
   Table-16 values, the Table-17 combined row, and every displayed Table-18
   entry.  The stale reproduction statement and verification capture were
   updated from earlier counts to `24/24`.

## Primary-source and visual custody

The empirical source remains exactly one committed PDF:

```text
c79552d62f4d4f4e85cfbbb00f135c1d985b596d9cdcde9bee57cfe4618f33dc  nist_bipm_2026.pdf
```

It has 31 PDF pages.  I independently rendered and visually inspected PDF
pages 6, 7, 10, 11, 15, 18, 23, 25, 26, and 27, in addition to replaying text
extraction.  The figures, equation numbers, table columns, units, signs,
superscripts, row order, and printed-page offset agree with the packet's
custody ledger:

- PDF p. 6 contains Figure 3 and equations (1)-(5), including the four-source,
  four-test torque law and the 37.67 degree two-extremum contrast.
- PDF p. 7 states below equation (7) that the source-mass contribution has
  zero second derivative at the torque extrema.
- PDF pp. 10-11 contain Table 1's mode sensitivities and Table 2's masses,
  radii, period, inertia, capacitances, and capacitance gradients.
- PDF pp. 15 and 18 contain the mass-integration formulae, sensitivity budget,
  and full/partial geometry comparators.
- PDF pp. 23 and 25 contain the background/inertia checks and autocollimator
  non-linearity summary.
- PDF pp. 26-27 contain Tables 15-18.  These are the only numerical objects
  used in the eight-row and four-result diagnostics.

The analyzer pins the PDF and all three GC dependencies before extracting any
field.  The added row-level tokens prevent a generic table-heading match from
standing in for custody of the actual numerical observations.

## Source formula and units

The paper's equation (1) is

\[
 N(\Delta\phi)=G{8m_sm_t\over R_s}\Gamma(\Delta\phi).
\]

At the two torque extrema, `Gamma_min` is approximately `-Gamma_max`, giving
the Table-15 relation

\[
 \Delta N=G{16\Gamma_{\max}m_sm_t\over R_s}=GA.       \tag{A1}
\]

This independently reproduces the packet's factor of 16.  Since `Gamma` is
dimensionless,

\[
 [A]={\rm kg^2/m},\qquad
 [G A]={\rm m^3\,kg^{-1}\,s^{-2}}{\rm kg^2/m}
      ={\rm kg\,m^2\,s^{-2}}={\rm N\,m}.             \tag{A2}
\]

The conversion `nN m -> 10^-9 N m` used in the analyzer is therefore correct,
and `Delta N/A` has the SI units of `G`.

The nominal stiffness statement is also correctly scoped.  At a torque
extremum, differentiating the source torque with respect to relative angle
vanishes.  Since GC06's `k_g` is the detector-angle derivative of torque per
`G`, the ideal source contribution maps to `k_g=0`.  This does not set the
finite-dither, miscentering, support, local-gravity, or total-apparatus
stiffness correction to zero.

## Independent eight-observation recomputation

Using only the visually transcribed Table-2 masses and Table-15 fields, a
40-digit decimal replay gives:

| configuration | mode | `A` (`kg^2/m`) | `Delta N/A` (SI) |
|---|---:|---:|---:|
| Sapphire | free | 209.479294634547 | `6.673642864985310e-11` |
| Sapphire | servo | 209.479294634547 | `6.672640379272503e-11` |
| Copper 0 deg | free | 467.458078301505 | `6.673946060223553e-11` |
| Copper 0 deg | servo | 467.458078301505 | `6.673582391248962e-11` |
| Copper 120 deg | free | 467.256225338079 | `6.673897170109821e-11` |
| Copper 120 deg | servo | 467.256225338079 | `6.673597548633620e-11` |
| Copper 240 deg | free | 467.252246069467 | `6.674253631166840e-11` |
| Copper 240 deg | servo | 467.252246069467 | `6.673825596841307e-11` |

All values agree with `RESULT.json` to floating-point precision.  The equal
free/servo `A` within a configuration is required by Table 15: those two
methods observe the same mass geometry but use different transfer/readout
calibrations.  The eight ratios remain deterministic summary reconstructions,
not estimates with a complete public covariance.

## Rank and identifiability audit

Let `A` be the eight-entry source column.

1. `rank(A)=1`, and every entry is nonzero, so fixed calibration and fixed
   remainders algebraically identify the single coefficient multiplying `A`.
2. `[A | I_8]` has rank 8, not 9.  An arbitrary torque remainder in every row
   absorbs the source column exactly.
3. Let `B_config` contain the four configuration indicators, each repeated in
   its free and servo row.  Because `A` is constant within each pair,
   `A=B_config c` exactly.  Hence `rank([A|B_config])=4`, equal to
   `rank(B_config)`, and `G` is not identifiable.
4. With only free and servo common offsets, the two method indicators have
   rank 2 and `rank([A|B_method])=3`; the varying configuration amplitudes make
   `A` independent of those two columns.  Identification here is conditional
   on that much stronger nuisance restriction, which the paper does not earn.
5. For a free source scale `s`, the two Jacobian columns of `G s A` are scalar
   multiples of `A`.  Their rank is one and only `p=Gs` is identified.

The direct numerical singular-value replay agrees with each rank.  The tiny
fifth singular value of `[A|B_config]` is `1.12e-17`; the least-squares alias
residual is `2.45e-13` in ordinary double precision and is algebraically zero
from the duplicated-row construction.

## Covariance diagnostics

The 14 Table-17 category rows reproduce the four RSS totals

\[
 (23.2262782,\ 30.3181464,\ 37.5634131,\ 93.8677261)\ {m ppm},
\]

which round to the displayed `(23.2,30.3,37.5,93.9)` ppm.  These same
tenth-ppm values occupy the Table-18 diagonal.  Table 16 separately displays
their whole-ppm rounding `(23,30,38,94)`.

The visually transcribed Table-18 correlation matrix has eigenvalues

\[
 (0.5000073063,\ 0.7657792989,\ 0.9020189020,\ 1.8321944928),
\]

rank four, and condition number `3.66433544`; it is positive definite.  Using
the four Table-16 central values and the Table-18 standard uncertainties gives
the independently reproduced formal weights

\[
 (0.617817179,\ 0.254323418,\ 0.126684518,\ 0.001174885),
\]

formal value `6.673611063585958e-11`, and formal standard uncertainty
`1.415391478337833e-15` (`21.2088 ppm`).  This remains a four-derived-result
diagnostic.  It is not the absent covariance of the eight torque observations
and excludes the later dark-uncertainty model.

The independent eight-observation diagonal Type-A replay likewise reproduces
`6.673920939386688e-11` with formal uncertainty
`2.477297709643543e-16` (`3.71191 ppm`).  Its deliberately unrealistic
independence, exact-denominator, zero-remainder, and zero-calibration-covariance
assumptions remain explicit.

## Accepted-G and CODATA exclusion

PDF p. 27 contains the authors' later hierarchical consensus model and states
its CODATA-centered prior.  The packet does not import that hierarchy:

- no Table-19 value, consensus posterior, dark-uncertainty scale, or CODATA
  numerical value appears in the result graph;
- the only four `G` values used are the primary Table-16 measured summaries,
  and only in the explicitly labelled Table-18 covariance diagnostic;
- the eight source columns use masses, `R_s`, and `Gamma`, never an accepted
  or fitted `G`; and
- ranks are functions of the design matrices alone.

The verifier now checks both the declared exclusion and the absence of a
Table-19/accepted-value literal from the serialized result.  This does not
pretend that the source PDF lacks CODATA prose; it proves that the excluded
hierarchy is not an input to this forward-readiness result.

## Reproduction and frozen payload

After the repairs, the analyzer reports:

```text
SUMMARY 24/24 checks passed
```

Its zero-argument output is byte-identical to `VERIFICATION.txt`, and its JSON
serialization is byte-identical to `RESULT.json`.  The accepted payload hashes
before adding this audit were:

```text
824fd6ea9dc62e564f18875f90f460a1358b9d9acb84a99a6a04d984c6a6d0ef  READINESS.md
ba77f99731775880b05f0885fe8d32c0d2b268541580c348e978821193d6efae  analyze_nist_bipm_g_readiness.py
a84725aa71ee713f82e49d5a19eef3212b64ea70942064cef2fcd8d637c6d7bb  RESULT.json
51357857edc00b0fd1339673fc7beb99abc16b83dcc7e3983043354c26aee82f  VERIFICATION.txt
95a411abfb6b33c17294cb9d18f1979d0fb5e4db192bd4fea3b1f151607bc805  SELF_AUDIT.md
```

## Final ceiling

The packet proves that the public NIST/BIPM summary supplies a real,
unit-correct finite-contrast source column and a nominal zero source-stiffness
entry at ideal extrema.  It also proves why that summary cannot execute the
full GC16 likelihood: raw transfer, eight-observation covariance, calibrated
remainders, complete mass geometry, independent source-scale covariance, and
the conserved apparatus source ledger remain absent.  It neither estimates
`G` independently nor tests lineage, RGRL, Gravity Formation Theory, or
gravity emergence.
