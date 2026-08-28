# Independent hostile audit — FX complete homogeneous H6 source response

**Audit date:** 2026-08-28  
**Audited lane:** `LANE_GRA_FX_F3_Q4_COULOMB_COMPLETE_H6_SOURCE_RESPONSE_V001`  
**Disposition:** `PASS`  
**Independent replay:** `139/139` checks passed  
**Builder replay:** `109/109` checks passed

## 1. Verdict

The frozen FX theorem survives independent hostile audit in its stated
scope.  Conditional on FV `S1`--`S10`/`FV-PURE`, the selected FO 180-state
winding component, homogeneous momentum, and fixed truncation through H6,
the formerly omitted diagonal words and Brillouin-Wigner/Feshbach folds are
complete.  They reduce exactly to the existing direct pair source plus
Hilbert-space identities:

\[
\begin{aligned}
Q_{\rm diag}^{(2)}&=U_dx^2(-D-40I),\\
Q_{\rm diag}^{(4)}&=U_dx^4\left(-\frac{37}{12}D-20I\right),\\
Q_{\rm diag}^{(6)}&=U_dx^6\left(-\frac{16247}{900}D-
\frac{374}{135}I\right).
\end{aligned}
\]

Thus the complete nonidentity homogeneous source through H6 is the audited
FW source with the exact substitution

\[
\rho\longmapsto\rho_E=\rho f_E(x),\qquad
f_E(x)=1-x^2-\frac{37}{12}x^4-\frac{16247}{900}x^6.
\]

For generic nonzero `f_E`, the exact finite hierarchy is

\[
\boxed{5_{\rm mod\ I}\longrightarrow3_{\operatorname{ad}_H}
\longrightarrow2_{\chi^R_{|0\rangle}}=2_{M_1}.}
\]

No remaining mathematical, custody, or documentary defect was found in the
frozen claim.

## 2. Algorithmically independent reconstruction

The audit executable imports or executes neither the FX builder nor the
FO/FW builders.  It independently:

1. rebuilt the `Z_30` diamond quotient from shifts `(0,1,5,19)`;
2. enumerated the absence of two- and four-link cycles, all 120 elementary
   hexagons, the 180-state winding component, its 420 undirected ring
   transitions, and six free length-thirty translation orbits;
3. used the girth-six alternating-cycle argument to classify every nonempty
   proper prefix through H6 as a `Q` state, while treating an empty-parity
   return as an intermediate-`P` fold;
4. evaluated 1,728,600 nonempty one-, two-, and three-link virtual subsets
   and 13,725 distinct signature/family weights with exact rational
   arithmetic;
5. replaced FX's multiset-permutation and closed-resolvent algorithm by a
   count-state dynamic program that interleaves the derivative of every
   hopping numerator and every endpoint-referenced resolvent; and
6. replaced FX's explicit fold formulas by an iterated dual formal-series
   fixed point in `t=h^2` through `t^3`.

That independent route reproduces

\[
a_2=-60,\qquad a_4=-35,\qquad a_6=-\frac{893}{9},
\]

all eighteen order/orbit derivative rows in FX11, and all six pair/triple
signature counts.  The independently replayed one-edge `(4)` and `(6)`
families vanish after the return-to-`P` exclusion.

## 3. Source scaling and response

The audit inferred each pair coefficient from orbit differences before
checking the residual on every orbit.  The residual is exactly isotropic and
state-independent at every order, establishing the pair-plus-identity
relations without assuming them.  Restoring dimensions gives
`h^m/U_d^(m-1)=U_d x^m`,
`J_6=(63/8)U_d x^6`, and `rho=8/(63x^6)`, so the substitution
`rho -> rho_E` follows exactly.

The response composition was then checked in two independent ways.  First,
integer full-component Grams give exact operator/commutator ranks `5 -> 3`
for a generic source and `4 -> 2` for the ring-only cancellation stratum.
Second, the audit reconstructed the exact zero-momentum block and used
rational polynomial projectors in `Q(sqrt(2))`.  Both responding projectors
have rank-one residues, their sum has rank two, and the exact first
commutator moment has rank two.  The pole gaps and residue vectors are the FW
ones with `rho` replaced by `rho_E`:

\[
\Delta_1=2+2\sqrt2,\qquad \Delta_2=4+2\sqrt2.
\]

No centered source amplitude remains outside those two energy projectors.

## 4. Hostile boundary tests

- **Diagonal-history completeness:** passed.  Even multiplicity partitions
  through length six are exhaustive, and intermediate `P` returns are folds.
- **Endpoint reference:** passed.  Every virtual gap derivative subtracts
  the initial ice endpoint energy before the dynamic program is evaluated.
- **Numerator ownership:** passed.  Each individual hopping factor supplies
  `-D_a/2`; no six-hop numerator derivative is omitted or double counted.
- **Fold ownership:** passed.  The formal fixed point independently restores
  every H4/H6 fold and its derivative, including nested derivative terms.
- **Identity quotient:** passed.  Isotropic Hilbert identities are retained
  in the exact reduction and removed only from nonidentity response ranks.
- **Scale conversion:** passed.  `Q=-2 dH/dj`, the powers of `U_d`, the H6
  normalization `63/8`, `rho`, and `rho_E` all agree exactly.
- **Static/dynamic rank separation:** passed.  Operator, commutator, ground
  spectral, and first-moment ranks are computed separately.
- **Custody:** passed.  All nine FX core files, all fifteen dependencies, the
  seven-file builder manifest, and its single-manifest seal match their
  pinned hashes.
- **Promotion control:** passed.  The finite polynomial zero is not called a
  physical threshold; it is not a physical threshold.  The finite rank-two response is not called a
  graviton, Ward identity, gravity derivation, or calculation of `G`.

## 5. Finite root and exact ceiling

Writing `y=x^2`, the truncated polynomial is strictly decreasing for
`y>=0`, and exact rational signs bracket its unique positive root between
`x=1/2` and `x=27/50`; its decimal checksum is
`x=0.5398271903...`.  At that algebraic stratum the finite ranks are
`4 -> 2 -> 2 -> 2`.  The audit agrees with FX that this is not a physical
threshold: H8 and higher terms are uncontrolled and may shift or remove it.

**PASS.**  The result is complete only for the homogeneous nonidentity
source on the selected FO component, under `FV-PURE`, fixed through H6 and
modulo `H_id`.  It does not establish a local or nonzero-momentum source,
other-component universality, H8 completion, CTP/Ward closure, a
thermodynamic massless tensor, RGRL-B, gravity, or Newton's constant.
