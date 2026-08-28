# Native-support nonzero-momentum complete-H6 response on the FO component

**Lane:** `LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001`  
**Date:** 2026-08-28  
**Status:** exact finite-source theorem plus sampled finite-response result  
**Claim class:** native-support source at FO quotient momentum `m=1`, on the
selected 180-state component, under `FV-PURE`, fixed through H6

## 1. Question and scope

FX completed the homogeneous (`m=0`) source through H6.  Its diagonal part
reduced to the direct pair source with coefficient

\[
 f_E(x)=1-x^2-{37\over12}x^4-{16247\over900}x^6,
 \qquad x={h\over U_d}.
 \tag{FY01}
\]

That identity could not simply be assumed at nonzero momentum: summing a
source over the graph before the Feshbach reduction erases where the source
was inserted.  FY therefore retains each insertion's native support, proves
that their uniform sum recovers FX/FW, and only then Fourier resolves the
source at one nonzero cyclic momentum.

Here **complete** means the direct pair source, every diagonal closed word
and BW/Feshbach fold through H6, and the Hermitian 720+720-history ring source,
on the selected FO 180-state winding component, under `FV-PURE`.  It does not
mean H8 completion, other winding sectors, a thermodynamic limit, continuum
locality, a complete spacetime Ward identity, a massless tensor pole, RGRL-B,
gravity, Newton's constant, or `G`.

## 2. Native-support source

The support set is

\[
 \mathcal S=\{A,B,e_0,e_1,e_2,e_3\}\times\mathbb Z_{30}.       \tag{FY02}
\]

Pair-source and virtual-gap derivatives are assigned to the physical `A/B`
vertex where that Coulomb energy lives.  A hopping-numerator derivative is
assigned to the midpoint of its physical link `e_a`.  A differentiated fold
inherits the support of the differentiated factor.  For a ring matrix
element, the forward and reverse endpoint-referenced ledgers are averaged
before applying `Q=-2 partial_j H`; hence each local source operator is
Hermitian.

For `m=1`, every native phase lies in `Q(zeta_240)`.  With support exponents

\[
 s_\sigma=(0,10,5,9,25,201)\quad
 \hbox{for}\quad(A,B,e_0,e_1,e_2,e_3),                         \tag{FY03}
\]

the normalized source is

\[
 Q^{ij}_{1}={1\over\sqrt{60}}
 \sum_{(\sigma,c)\in\mathcal S}
 \zeta_{240}^{,8c+s_\sigma}\,q^{ij}_{\sigma c}.
 \tag{FY04}
\]

The conjugate mode obeys `Q_29=Q_1^dagger`.  Equation (FY04) is the exact
translation label of the finite FO graph; it is not an assumption of
continuum spatial locality.

## 3. Exact recovery and completeness gates

The verifier independently re-enumerates all diagonal even-multiplicity
families `(2)`, `(2,2)`, `(4,2)`, `(2,4)`, and `(2,2,2)`, retaining only
irreducible proper-prefix histories and restoring every fold through H6.  It
recovers exactly

\[
 a_2=-60,\qquad a_4=-35,\qquad a_6=-{893\over9},               \tag{FY05}
\]

and all eighteen FX derivative rows.  The native ring inventory contains
420 undirected matrix entries in fourteen free translation orbits.  Each
representative includes all 720 forward and 720 reverse histories.

At uniform momentum all support phases coalesce, and

\[
 \sqrt{60}\,Q_0^{\rm FY}=Q_{\rm homogeneous}^{\rm FX/FW}       \tag{FY06}
\]

separately for the direct pair, H2, H4, H6, and all ring entries.  Thus the
nonzero-momentum calculation is a refinement of the frozen complete
homogeneous source, not a different source model.

## 4. Exact `m=1` diagonal-lift theorem

Let `D_1=Q_pair(m=1)/U_d`.  Reducing exact rational source polynomials modulo
the irreducible cyclotomic polynomial `Phi_240` gives

\[
\boxed{\begin{aligned}
 {Q_{\rm diag}^{(2)}(1)\over U_dx^2}&=-D_1,\\
 {Q_{\rm diag}^{(4)}(1)\over U_dx^4}&=-{37\over12}D_1,\\
 {Q_{\rm diag}^{(6)}(1)\over U_dx^6}&=-{16247\over900}D_1.
\end{aligned}}                                                \tag{FY07}
\]

This was a neutral exact test: each order had to return either the displayed
equality or a nonzero exact cross-minor proving nonproportionality.  All three
returned equality.  At H2 the stronger ledger identity is visible directly:
`Qdiag2+Qpair` is the same cell-uniform 120-link density on every FO orbit,
so its nonzero cyclic Fourier modes vanish exactly.

Consequently the complete diagonal nonidentity source at `m=1` carries the
same `f_E(x)` as the homogeneous source:

\[
 {Q_{\rm eff}^{(\le6)}(1)\over J_6}
 =\rho f_E(x)D_1+R_1,
 \qquad \rho={8\over63x^6},                                   \tag{FY08}
\]

where `R_1` is the native Fourier transform of the differentiated Hermitian
ring source.  An exact nonzero cyclotomic remainder on an off-diagonal ring
entry proves `R_1 != 0`; because `D_1` is diagonal, the ring is an independent
operator direction.

## 5. Sampled finite response

The FO Hamiltonian is not changed by the source.  Therefore the source can
expose or suppress existing finite-graph poles but cannot move or create
eigenvalues.  At the two declared samples `x=2/5` and `x=1/2`, independent
180-state commutator and Lehmann replays give

\[
 \boxed{6\ \longrightarrow\ 6\ \longrightarrow\ 6
 \ \longrightarrow\ 6\ \longrightarrow\ 6}                 \tag{FY09}
\]

for coordinate-operator rank, `ad_H` rank, ground spectral rank, the static
positive-frequency spectral-kernel rank, and the first gap-weighted moment
rank.  These are **sampled finite ranks**, not a generic-in-`x` theorem.

The two TT contractions have ground-image rank two at both samples.  Four
finite FO gaps carry response (in `J_6` units):

\[
 3.194109035554332,\quad 3.490165912028476,\quad
 6.166688337463908,\quad 9.139267639373482,                    \tag{FY10}
\]

with residue ranks `(1,3,1,1)`.  There is no zero-energy pole at either
sample.  The lowest pole lies above FO's finite two-one-link threshold proxy
`2.059674505691458`; no continuum stability or particle interpretation is
inferred from that comparison.

## 6. Ward diagnostic and strict ceiling

Using the frozen FO wavevector, FY evaluates the purely spatial contraction

\[
 L^j=k_iQ^{ij}.                                                \tag{FY11}
\]

It is nonzero for the pair source and for both complete sampled sources.  The
norm ratios are approximately `1.2101`, `1.2121`, and `1.2735`.  This rules
out a naive claim that the spatial source is already transverse.  It is not a
complete Ward-theorem failure: FY has not constructed the temporal density,
current operators, contact terms, or a continuum conservation law needed by
a spacetime Ward identity.

### Theorem `NSM1-H6`

Under the frozen FO/FV/FW/FX premises and `FV-PURE`, the native-support source
is exactly composable at quotient momentum `m=1`; it recovers the complete
homogeneous H6 source at `m=0`, preserves the three diagonal `f_E`
coefficients exactly at `m=1`, and retains an independent nonzero ring source.
At `x=2/5` and `x=1/2`, the finite response has full six-channel hierarchy
`6->6->6->6->6`, TT rank two, and four nonzero FO pole gaps.

The theorem establishes a complete source-resolved **finite graph** step.  It
does not establish emergent continuum locality, a Ward identity, a massless
graviton, RGRL-B, gravity emergence, `G`, or Newton's law.
