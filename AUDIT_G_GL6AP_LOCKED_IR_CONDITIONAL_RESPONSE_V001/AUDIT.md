# Distinct hostile audit — GL6AP locked-sector conditional response

**Target:** `LANE_CROSS_RFT_GRA_GL6AP_LOCKED_IR_CONDITIONAL_RESPONSE_V001/`  
**Frozen theorem SHA-256:** `aab21a99ecd0c6696084dd0a22ed2489533588f2171bf0c319ffc00a4922033e`  
**Frozen author-manifest SHA-256:** `1c6a3db34bda121f7b7fdc64a85aaed0d146dcb12d8f960449d10d873ff94e1f`  
**Frozen author-seal-file SHA-256:** `f56c24594e041485edc6d14cff91a3317bc400fc92ece595ef590291331a8288`  
**Disposition:** `PASS__LINK_KERNEL_T2_NOT_PAIR_E__HOM_S4_ZERO__ALL24_Q4_AUTOMORPHISMS__NATIVE_LOOP_BREAKS_E_COUNT__MASS_AND_CUBIC_ALLOWED__RECIPROCAL_QUADRATIC_CLASSIFICATION_EXACT__SPECTRAL_AND_PHYSICAL_MOMENTUM_CEILINGS_SOUND`

## 1. Custody and independence

The author froze and sealed eleven GL6AP files before this audit.  Their exact
hashes are pinned in `AUDITED_TARGETS.sha256`; the author manifest and seal
pass.  The seven direct dependency rows are confined to the sealed GL6AN
author snapshot and its distinct sealed hostile audit.  Replaying those packet
gates also verifies GL6AN's transitive sealed GL6AK ancestry.  Mutable GL6AL
and the parallel GL6AO/GL6AQ lanes are not premises.

The independent replay imports no GL6AP or GL6AN verifier.  The frozen author
physics replay passes `434/434` in normal and optimized modes, the author
packet passes `100/100` in both modes, the upstream GL6AN packet passes
`79/79`, and the GL6AN hostile-audit packet passes `58/58`.

## 2. Link constraint and representation mismatch

For

\[
 B(\chi)=\begin{pmatrix}1&1&1&1\\z_1&z_2&z_3&z_4\end{pmatrix},
 \qquad s=\sum_a z_a,
\]

direct multiplication verifies

\[
 \Pi_K=I-B^\dagger(BB^\dagger)^{-1}B,
 \qquad BB^\dagger=\begin{pmatrix}4&\bar s\\s&4\end{pmatrix}.
\]

For `|s|<4`, this is a Hermitian rank-two projector annihilated by `B`.  At the
trivial character the rank drops and

\[
 \ker B(1)=\mathbf1^\perp\cong T_2,
 \qquad \dim T_2=3.
\]

With `sum theta_a=0`, exact Taylor expansion gives

\[
 4-|s|={1\over2}\sum_a\theta_a^2+O(|\theta|^4),
\]

so this is a squared singular value; the singular value is linear.  The
generic rank-two space is a direction-dependent plane in complexified `T2`,
not a fixed local irrep.

The independent `S4` character table gives

\[
 \langle\chi_{T_2},\chi_{T_2}\rangle=1,
 \quad \langle\chi_E,\chi_E\rangle=1,
 \quad \langle\chi_{T_2},\chi_E\rangle=0.
\]

The unsigned `K4` port/pair incidence has rank four and nullity two.  Locked
pair variations lie in that nullspace, the centered opposite-pair plane `E`.
Therefore

\[
 \operatorname{Hom}_{S_4}(T_2,E)=0,
\]

and equality of generic dimensions supplies no equivariant identification.

## 3. Quotient symmetry and loop nonconservation

The replay rebuilds the declared `Q4` with 64 cells, 128 constraint nodes and
256 links.  For every `sigma in S4`, it independently forms the integer matrix

\[
 A_\sigma d_j=d_{\sigma(j)}-d_{\sigma(4)}.
\]

The 24 determinants are twelve `+1` and twelve `-1`.  Reducing modulo four,
each map is a vertex and edge bijection and, with the prescribed child shift,
preserves every parent-child incidence.

An Edmonds--Karp completion, distinct from the author's Dinic completion,
extends the specified alternating hexagon collar to a global degree-two `Q4`
configuration.  Toggling the six cycle links changes the local opposite-pair
type counts by

\[
 (N_1,N_2,N_3):(1,2,3)\longmapsto(3,2,1),
 \qquad \delta N=(2,0,-2).
\]

Independent enumeration of all `6!=720` flip orderings gives

\[
 \sum_{\pi\in S_6}\prod_{j=1}^{5}{-1\over E(S_{\pi,j})}
 =-{63\over8}.
\]

Thus the corresponding sixth-order matrix element is nonzero.  The `S4` orbit
of `(2,0,-2)` contains all six coordinate permutations and spans the full
centered pair-`E` plane.  Since the kernel of
`w -> [w dot N_E,H_eff^(6)]` is `S4` invariant, irreducibility leaves only the
zero kernel.  No nonzero uniform bare pair-`E` count is conserved at sixth
order.  This does not exclude a nonlocal winding label, emergent conservation
in a selected phase, or a statement about an uncomputed all-orders operator.

## 4. Reciprocal analytic inverse-response classification

Exact character products give

\[
 \operatorname{Sym}^2(E)=A_1\oplus E,
 \qquad
 \operatorname{Sym}^2(T_2)=A_1\oplus E\oplus T_2,
 \qquad
 \operatorname{End}(E)=A_1\oplus A_2\oplus E.
\]

Consequently there are exactly two quadratic spatial covariants into symmetric
internal matrices: one scalar `A1` contraction and one traceless `E`
contraction.  There is no linear map from `T2` to `End(E)`.  Reciprocity makes
odd spatial orders internally antisymmetric, hence of type `A2`; neither
`T2` nor `Sym^3(T2)` contains `A2`.  Conditional on existence and reciprocal
character analyticity, the complete form through cubic spatial order is

\[
 \Gamma_E^R(\omega,\theta)=
 [a_0^R(\omega)+c_0^R(\omega)I_2(\theta)]I_2
 +c_2^R(\omega)\mathcal T(Q_E(\theta))+O(|\theta|^4).
\]

The scalar in `End(E)` proves that a mass is symmetry allowed.  Also
`Sym^3(E)=A1+A2+E`, so a cubic order-parameter invariant is allowed.  Symmetry
therefore protects neither masslessness nor a continuous transition.

## 5. Spectral and infrared ceilings

The analytic nondissipative pole alternatives are explicitly conditional on
a selected invariant state, a growing-quotient or thermodynamic completion,
inverse-response regularity, positive inertia/stiffness, and absence of a
lower continuum.  Fixed `Q4` has no `theta -> 0` sequence.

For a chosen ground-state positive-frequency measure, the exact bound is

\[
 \Delta_E(\chi;u)\le {f_u^+(\chi)\over S_u^+(\chi)},
 \qquad S_u^+>0.
\]

The audit explicitly checks that including an elastic atom at zero can make
the full-frequency quotient fall below the least positive support point; the
restriction to `(0,infinity)` is therefore necessary.  Threshold closure is
not an isolated pole.  A pole requires a positive atom, while nonvanishing
infrared visibility additionally requires its residue to stay bounded away
from zero.  Ground-state criteria are not promoted to generic stationary or
finite-temperature Liouvillian spectra.

Finally, the character remains only an authenticated translation label.
Physical momentum additionally requires a calibrated embedding and length,
evidence that the automorphism is physical translation, a controlled growing
or refinement family, and a selected phase with stable response scaling.  No
masslessness, pole, physical momentum, cone, gauge phase, photon, graviton,
stress/Ricci/Einstein law, gravity, or `G` is derived.

**Hostile verdict: PASS.**
