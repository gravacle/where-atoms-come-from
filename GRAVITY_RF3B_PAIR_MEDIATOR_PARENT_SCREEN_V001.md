# RF3b pair-mediator parent screen

**Record:** `GRAVITY-RF3B-PAIR-MEDIATOR-PARENT-SCREEN-V001`  
**Date:** 2026-08-29  
**Role:** design screen for `RF3b` in
`GRAVITY_RECORD_FIRST_DIRECT_PROOF_PLAN_V001.md`  
**Disposition:** exact algebraic `DPAR` dilation found; proposed oscillator
block alone does not own spatial momentum; smallest physical repair identified;
no separate theorem lane should be promoted before that repair is embedded and
audited

## 1. Question screened

The proposed local completion puts one oscillator on every unordered pair of
the four q4 terminals,

\[
 H_{ab}[F]=\Omega_{ab}(F)
 \left[b_{ab}^\dagger b_{ab}+{1\over2}
 -\eta D_{ab}(b_{ab}+b_{ab}^\dagger)\right],
 \qquad D_{ab}=Z_a-Z_b,                           \tag{S01}
\]

with

\[
 \Omega_{ab}(F)=\Omega_0
 \left({|Fr_{ab}|^2\over r_0^2}\right)^{-1/2}.
                                                               \tag{S02}
\]

Here `eta` denotes the oscillator displacement coupling; it is kept distinct
from FU's `DPAR` slope.  The screen asks four separate questions:

1. Does eliminating the oscillator reproduce the F3 degree-pair energy and
   FU's complete ideal-Coulomb `DPAR` source?
2. Can every F3 flip, including the GA/GD charge and recoil dressing, act
   exactly within the conditional mediator sector?
3. Does the displayed oscillator block own physical spatial momentum and
   stress?
4. What is the smallest autonomous repair that can pass `RF3b` without
   introducing another record mechanism or reverting to a particle-first
   route?

The answers are respectively **yes, yes after one exact common dressing, no,
and a dynamical relational terminal/scaffold completion**.

## 2. Exact oscillator dilation

Put

\[
 P_{ab}=Z_aZ_b,
 \qquad c_{ab}=b_{ab}-\eta D_{ab}.
\]

Since

\[
 D_{ab}^2=2(I-P_{ab}),                            \tag{S03}
\]

completing the square gives the operator identity

\[
 H_{ab}[F]=\Omega_{ab}(F)
 \left[c_{ab}^\dagger c_{ab}+{1\over2}
 -2\eta^2(I-P_{ab})\right].                       \tag{S04}
\]

Within the displayed ansatz, the clean exact choice is

\[
 \boxed{\eta={1\over2},\qquad \Omega_0=U_d.}       \tag{S05}
\]

Then the source-dependent identity in each pair block cancels and

\[
 \boxed{
 H_{ab}[F]=\Omega_{ab}(F)
 \left(c_{ab}^\dagger c_{ab}+{1\over2}P_{ab}\right).}            \tag{S06}
\]

After adding the single source-independent reference `U_d I`, the six-pair
parent is

\[
 H_{\rm PM}[F]=U_dI+
 \sum_{a<b}\Omega_{ab}(F)
 \left(c_{ab}^\dagger c_{ab}+{1\over2}P_{ab}\right).             \tag{S07}
\]

Let `K_0` be the conditional-vacuum sector
`c_ab K_0=0` for every pair.  On all sixteen logical q4 states, including the
ice `P` and off-ice/virtual `Q` states,

\[
 \boxed{
 H_{\rm PM}[F]\big|_{K_0}
 =U_dI+{U_d\over2}\sum_{a<b}
 \left({|Fr_{ab}|^2\over r_0^2}\right)^{-1/2}P_{ab}.}            \tag{S08}
\]

At `F=I`, (S08) is exactly

\[
 U_dI+{U_d\over2}\sum_{a<b}P_{ab}
 =U_d(d-2)^2.                                      \tag{S09}
\]

Thus the oscillator block is an exact dilation of the complete ideal-Coulomb
FU pair family, not merely a fit to its first derivative.  The special choice
(S05) is the unique simple choice inside (S01), up to oscillator parity
`b_ab -> -b_ab`, with only the one constant
reference in (S07), that both matches the nonidentity coefficient and removes
the geometry-dependent zero-point identity.  General `eta` can match the pair
coefficient only by retaining and owning an additional source-dependent
identity or by adding a tuned counterterm; neither may be silently discarded.

This is an **existence construction**, not a deeper derivation of FU's radial
law.  Equation (S02) already inserts the `1/r` dependence that FU derived
conditionally from its grounded pair-resolved electromagnetic solder.
Completing the square realizes that law in an oscillator parent; it does not
derive the law from source-off F3 incidence alone.

## 3. Exact first and second source jets

For a fully frozen source chart

\[
 F(j)=I-{1\over2}j,                               \tag{S10}
\]

write `u_ab=r_ab/r_0`,
`Rhat_ab=u_ab u_ab^T`, and
`s_ab(A)=u_ab^T A u_ab=A:Rhat_ab`.  On `K_0`, exact differentiation of
(S08) gives

\[
 D H_{\rm PM}[A]
 ={U_d\over4}\sum_{a<b}s_{ab}(A)P_{ab},           \tag{S11}
\]

and hence, under FU's convention `Q[A]=-2 D H[A]`,

\[
 \boxed{
 Q[A]=-{U_d\over2}\sum_{a<b}s_{ab}(A)P_{ab}.}     \tag{S12}
\]

This is exactly FU's ideal-Coulomb slope
`lambda_DPAR=-1/2`.  The mixed second derivative is

\[
 \boxed{
 D^2H_{\rm PM}[A,B]
 ={3U_d\over8}\sum_{a<b}s_{ab}(A)s_{ab}(B)P_{ab}
 -{U_d\over16}\sum_{a<b}
 u_{ab}^{\mathsf T}(AB+BA)u_{ab}\,P_{ab}.}        \tag{S13}
\]

Equation (S13) is a lawful contact only after the nonlinear source chart
(S10) has been frozen.  FU fixed
`F(j)=I-j/2+O(j^2)` and therefore fixed the first jet but not the second.
A different `O(j^2)` completion changes the seagull/contact while leaving
(S11)--(S12) intact.  The metric source chart must be selected and carried
consistently through the complete parent; it cannot be chosen after a Ward
test is inspected.

The full oscillator parent also has source derivatives proportional to
`c_ab^\dagger c_ab`.  They vanish on the invariant `K_0` sector, but they must
remain in the full ledger until invariance under every admitted interaction
has been proved.

## 4. Exact conditional-state and flip transport

Define one common conditional displacement

\[
 \mathcal P=
 \exp\!\left[{1\over2}\sum_{a<b}
 D_{ab}(b_{ab}^\dagger-b_{ab})\right].             \tag{S14}
\]

Then

\[
 \mathcal P^\dagger c_{ab}\mathcal P=b_{ab},
 \qquad
 K_0=\mathcal P
 \left(\mathcal H_{q4}\otimes|0\rangle^{\otimes6}\right).       \tag{S15}
\]

For every inherited logical flip, use

\[
 \widehat X_a=\mathcal P(X_a\otimes I)\mathcal P^\dagger.        \tag{S16}
\]

It follows exactly that

\[
 \widehat X_a\mathcal P(|z\rangle\otimes|\boldsymbol n\rangle)
 =\mathcal P(X_a|z\rangle\otimes|\boldsymbol n\rangle),         \tag{S17}
\]

for every logical configuration and every conditional mediator occupation.
The common conjugation preserves the complete Pauli algebra and makes `K_0`
invariant.  A bare `X_a` does not do this: it changes three `D_ab` values and
therefore excites the corresponding conditional oscillator states.

The same construction composes with GA and GD.  If
`mathsf X_a^(GA/GD)` is their already charge- and recoil-conserving flip, the
admitted mediator flip is

\[
 \widehat{\mathsf X}_a
 =\mathcal P\mathsf X_a^{(GA/GD)}\mathcal P^\dagger.              \tag{S18}
\]

Because (S14) is built from the logical `Z` operators and internal oscillator
coordinates, it commutes with GA total charge and, before a physical spatial
solder is added, with GD's total mechanical momentum.  Their exact
commutators and encoded equivalence are therefore preserved by conjugation.
Every flip used in a Feshbach history must receive (S18); dressing only the
external states would not preserve the inherited virtual ledger.

There is one covariance condition that must be explicit.  `D_ab` is oriented.
Under a tetrahedral permutation that reverses the chosen orientation of a
pair, the real oscillator coordinate and momentum must also reverse,
equivalently `b_ab -> -b_(pi(a)pi(b))`.  With this signed-edge action the
linear coupling, oscillator number, polaron, and six-pair family are exactly
`S4` covariant.  Treating the oscillator as an unsigned pair label while
retaining `D_ab` would break the claimed covariance.

## 5. Sharp physical ceiling of the displayed candidate

As written, (S01)--(S02) does **not** pass physical `RF3b`.

- `b_ab+b_ab^dagger` is an internal oscillator coordinate.  No canonical
  spatial position or momentum is assigned to it.
- `Omega_ab(F)` is an externally programmed coefficient on a fixed root.
  It supplies a source derivative but no operator that receives the
  equal-and-opposite force.
- The six oscillators are mutually isolated.  They provide neither spatial
  transport nor a causal propagation law.
- Adding GA charge transfer and GD flip recoil does not create the missing
  momentum owner for the diagonal FU pair energy.  GJ already shows why the
  direct GD half-kick density cannot close the supplied-embedding longitudinal
  pair-source mismatch.

Therefore the candidate proves an exact source-family dilation and exact flip
transport, but it does not yet prove that the pair source is physical stress.
Calling `b_ab` a mediator does not change this operator-type boundary.

## 6. Smallest physical repair

The narrow repair is to give the four terminals an owned relational
configuration.  Introduce positions and conjugate momenta

\[
 [R_a^i,\Pi_b^j]=i\hbar\delta_{ab}\delta^{ij},
\]

and replace (S02) at source off by

\[
 \Omega_{ab}(R)=U_d{r_0\over|R_b-R_a|}.            \tag{S19}
\]

Add only a complete kinetic/support Hamiltonian depending on relative
coordinates, together with the already required GA/GD reservoir, controller,
and boundary factors.  Every active transfer to one of those factors must use
a relative-coordinate recoil operator and every active exterior port must
include its exterior momentum owner.

In polaron variables, define

\[
 A_{ab}=c_{ab}^\dagger c_{ab}+{1\over2}P_{ab},
 \qquad
 H_{ab}(R)=U_d{r_0\over r_{ab}}A_{ab},
 \quad r_{ab}=|R_b-R_a|.                          \tag{S20}
\]

Then the pair sector is genuinely translation invariant,

\[
 \left[\sum_a\Pi_a,H_{\rm PM}(R)\right]=0,        \tag{S21}
\]

and owns exact equal-and-opposite central forces.  With
`rho_ab=R_b-R_a`,

\[
 F_{a\leftarrow ab}
 =-U_dr_0A_{ab}{\rho_{ab}\over r_{ab}^3},
 \qquad
 F_{b\leftarrow ab}=-F_{a\leftarrow ab}.          \tag{S22}
\]

The affine derivative of (S20) is the corresponding bond stress.  On `K_0`
its logical part is proportional to `P_ab Rhat_ab`, exactly the diagonal
operator type that the flip-recoil-only completion lacked.  Thus this repair
can in principle enter the required momentum-density/stress continuity
identity; it does not guarantee that the completed Ward identity will pass.

The support is load bearing.  In sectors with `P_ab=-1`, (S20) contains an
attractive `-1/r` contribution and is unbounded as a pair collapses unless an
owned hard core, binding potential, or relational constraint prevents it.
That stabilizer contributes energy, stress, and source contacts and cannot be
omitted from the complete family.

There is also an exact-custody distinction:

- A smooth tetrahedral binding potential produces position fluctuations.
  Then `Omega_ab` is not exactly `U_d` on the support Hilbert space, the old
  fixed-support F3 Hamiltonian is only an approximation, and the enlarged
  parent requires a fresh complete source and Feshbach calculation.
- The narrow exact witness is a translation-owning **rigid relational q4
  scaffold**, or another exactly invariant regular-shape code, with all six
  rest lengths equal to `r_0`.  Its constraint reactions, orientation/support
  sector, and source contacts must be explicit.  It may not be replaced by a
  fixed c-number grid with no recoil owner.

The rigid witness is enough to test the finite-cell ownership problem.  It
does not by itself prove a common cone or a continuum limit.  Later gluing
must put the bond/scaffold blocks on the same finite q4 family and show that
their local energy and momentum transfers compose without unowned long-range
interactions.

The physical-length and norm in (S19) are also conditional inputs.  They must
be the independently earned relational q4 solder allowed by FU `S2/S9`, not
the gravitational metric whose origin the route is trying to prove.  Using a
background Euclidean norm without stating that custody would be circular.

## 7. Exact finite checks for the embedded parent

The following checks are necessary and remain small enough for an exact
finite verifier.

1. Enumerate all `16` q4 `Z` configurations and verify (S03), (S08), and
   (S09), including every ice and off-ice energy.
2. Verify the six-root source map has rank six and verify (S11)--(S13) on a
   fixed six-element basis of `Sym^2(R^3)`.
3. Verify all `24` tetrahedral permutations with the signed-edge oscillator
   action.
4. Verify the four common-polaron flip intertwiners on all `16` logical
   configurations, and symbolically on arbitrary conditional oscillator
   occupations; verify preservation of the complete flip algebra.
5. Replace each flip by (S18) and verify GA total charge, GD total momentum,
   full `P+Q` encoded equivalence, and invariance of `K_0`.
6. Verify positivity of every oscillator excitation gap and stability of the
   **complete** support-plus-pair Hamiltonian.  The isolated `1/r` pair block
   does not pass this check.
7. In the repaired parent, verify (S21)--(S22), central torque balance, and
   the equality between the affine source derivative and the bond stress.
8. Verify that the regular q4 shape sector is exactly invariant.  If it is
   not, stop inheriting FY and redo the complete source/Feshbach calculation
   with controlled support fluctuations.
9. Include kinetic, constraint/binding, mediator, reservoir, recoil,
   controller, source-work, and boundary terms in one first- and second-source
   ledger.  Verify mixed-partial symmetry and retain every identity/contact.
10. Derive the native node/bond divergence and verify that an active port has
    an explicit exterior charge, energy, and momentum owner.  A controller-off
    response sector must be an invariant sector of one autonomous Hamiltonian.

These checks require no new record criterion and no graviton premise.

## 8. Precise `RF3b` pass conditions

`RF3b` passes only if one embedded physical parent establishes all of the
following simultaneously:

1. **Exact carrier custody:** the full sixteen-state q4 `P+Q` carrier and all
   admitted virtual flips reduce exactly to F3 plus the complete FU `DPAR`
   family, or every difference is explicitly retained and recomputed.
2. **Complete dressing:** every logical flip uses the common polaron together
   with GA charge transfer and GD recoil, with exact algebra, charge, and
   momentum conservation.
3. **Physical pair owner:** the pair energy depends on owned relational
   coordinates and has canonical momentum recipients satisfying an exact
   local equal-and-opposite force law.
4. **Stable support:** the regular q4 realization is dynamically stable and
   its binding or constraint energy and stress are included once.
5. **Complete source:** `h_00`, `h_0i`, and `h_ij` couple to every admitted
   term before projection, with all first and second derivatives, work terms,
   and contacts retained.
6. **Native locality:** the terminal, mediator, reservoir, support,
   controller, and boundary allocation defines the native divergence rather
   than borrowing FY's supplied embedding contraction.
7. **Closed ports:** every active boundary exchange owns charge, energy, and
   recoil; the response hold is an invariant sector of the same autonomous
   parent.
8. **No circular geometry:** the length/coframe input is independently earned
   relational structure and is not imported from the Einstein endpoint.

Passing (1)--(8) closes the finite-cell parent construction.  Transporting
its complete source and contact packet through the elementary-hexagon and
native H6 reductions remains `RF3c--RF3d`; finite-family causal gluing remains
`RF3e`.

## 9. Decision status

No theory decision is required now.  The oscillator construction should not
be adopted as a new postulate or promoted as a gravity theorem.  The next
recommended calculation is the smallest rigid-relational dynamic-node
embedding that tests the eight `RF3b` pass conditions.  A theory decision
would arise only if more than one inequivalent physical support architecture
passes those conditions and the existing record-first premises do not select
between them.
