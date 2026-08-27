# ET/GSGB gamma-length convention crosswalk

**Date:** 2026-08-27

**Status:** post-seal convention clarification; no scientific source bytes
changed

The sealed GSGB and ET packets use equivalent metric normalizations but used
the same informal length label for two differently normalized quantities.
This note makes the conversion explicit and does not alter either result.

In the GSGB/EU QFI convention (denote its conversion length by
\(\ell_F\)),

\[
 s=\ell_F^2 F,
 \qquad F=fP,
 \qquad \ell_F^2 f={4a^2\over3}.
\]

In the ET Bures-line convention (denote its conversion length by
\(\ell_B\)),

\[
 q={\ell_B^2\over4}Q,
 \qquad Q=fP,
 \qquad \ell_B^2 f={16a^2\over3}.
\]

When `F` and `Q` denote the same Fisher/QFI bilinear form, these definitions
give the same physical metric exactly when

\[
 \boxed{\ell_B=2\ell_F}
\]

for positive length conventions.  Thus the factor of four in the squared
scale locks is a convention conversion, not a physical disagreement.

Future joins must retain the typed symbols `ell_F` and `ell_B` (or state the
conversion before using a common symbol).  In particular, no proof may compare
the two raw length labels as though their normalization were identical.  This
crosswalk does not derive either scale, identify Fisher information with a
physical coframe, or promote ET to gravity.
