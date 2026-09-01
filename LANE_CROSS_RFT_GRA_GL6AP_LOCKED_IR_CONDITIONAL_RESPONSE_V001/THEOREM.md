# Locked-sector infrared conditional-response theorem and native E-pole no-go

**Short name:** `GL6AP V001`  
**Date:** 2026-08-31  
**Status:** author frozen and sealed; exact checks pass; independent hostile
audit required before promotion  
**Claim class:** exact representation/constraint decomposition inherited from
sealed GL6AN; exact native-loop nonconservation witness; symmetry-complete
quadratic long-wavelength form under explicit analyticity/state hypotheses;
conditional spectral criteria and no-go

**Not claimed:** a selected ground, Gibbs, KMS, vacuum, or phase; a gauge
theory, photon, graviton, hydrodynamic law, infrared pole, gaplessness,
physical momentum, length, speed, Lorentz/common cone, stress tensor,
Ricci/Einstein response, gravity, or `G`.

## 1. Sealed premise and exact question

Use the sealed and independently hostile-audited `GL6AN V001` lock line.  On
its strict inherited positive subfamily,

\[
 \varepsilon_\star=-6U_d,
 \qquad \Delta=4U_d(d_\star-2)>0,
 \qquad U_d>0,\quad d_\star>2,                             \tag{AP01}
\]

the homogeneous F3 Hamiltonian is, modulo a scalar,

\[
 H_{\rm lock}=-h\sum_eX_e+U_d\sum_vq_v^2,
 \qquad q_v=\sum_{e\ni v}n_e-2.                            \tag{AP02}
\]

On the declared finite period-four girth-six quotient `Q_4`, canonical
Hermitian Kato/Schrieffer--Wolff perturbation theory gives common scalar terms

\[
 H_{\rm eff}^{(2)}=-{Mh^2\over2U_d}P_Q,
 \qquad
 H_{\rm eff}^{(4)}=-{7M\over24}{h^4\over U_d^3}P_Q,         \tag{AP03}
\]

and at least one native alternating hexagon matrix element

\[
 \langle f|H_{\rm eff}^{(6)}|i\rangle
 =-{63\over8}{h^6\over U_d^5}.                             \tag{AP04}
\]

The complete sixth-order diagonal and loop operator is not an input because
GL6AN does not derive it.  The present question is narrower: what infrared
structure follows from the exact constraint, symmetry, and the existence of
the native loop move before a phase or continuum interpretation is assumed?

## 2. Exact link-constraint decomposition

Orient every active link from its parent endpoint to its child endpoint and
put

\[
 e_a(x):=n_{(x,a)}-{1\over2}.                               \tag{AP05}
\]

The strict lock is the oriented divergence constraint: at a parent it is
`sum_a e_a=0`, and at a child it is the same equation with the opposite row
sign.  For a translation character `chi`, use the GL6AN phase gauge

\[
 B(\chi)=
 \begin{pmatrix}1&1&1&1\\z_1&z_2&z_3&z_4\end{pmatrix},
 \qquad |z_a|=1,\qquad s=\sum_a z_a.                       \tag{AP06}
\]

For `|s|<4`, the exact orthogonal projector onto the locked link tangent is

\[
 \boxed{
 \Pi_K(\chi)=I_4-B(\chi)^\dagger
 \begin{pmatrix}4&\bar s\\s&4\end{pmatrix}^{-1}B(\chi).}   \tag{AP07}
\]

It has rank two.  Every projected locked-state link correlation satisfies

\[
C_e(\chi,z)=\Pi_K(\chi)C_e(\chi,z)\Pi_K(\chi).             \tag{AP08}
\]

Equation (AP08) applies on a finite periodic quotient with its ordinary
locked projector, or in any infinite-volume state in which every local
degree constraint is null.  It does not insert a global locked projector
into the infinite quasi-local algebra.

At the trivial character the two rows of `B` coincide and

\[
 \Pi_K(1)=I_4-{1\over4}\mathbf1\mathbf1^T.                 \tag{AP09}
\]

The four-port permutation representation is

\[
 \mathbb R^4\cong A_1\oplus T_2,                           \tag{AP10}
\]

so the rank-three space in (AP09) is the native `T2`, not `E`.  Near the
trivial character, write `z_a=exp(i theta_a)` with `sum_a theta_a=0`.  The two
constraints become

\[
 \sum_a e_a=0,
 \qquad \sum_a\theta_a e_a=O(|\theta|^2\|e\|).             \tag{AP11}
\]

Thus the generic rank-two kernel is the plane transverse to the character
direction inside the three-dimensional `T2`.  Its dimension two does not
turn it into the two-dimensional local pair `E` representation.

The constraint Gram eigenvalue normal to that plane is

\[
 \lambda_-(\chi)=4-|s|
 ={1\over2}\sum_a\theta_a^2+O(\theta^4),                   \tag{AP12}
\]

while the corresponding singular value is its square root and is linear in
`|theta|`.  Equation (AP12) is a static constraint statement only.

## 3. Exact locked pair sector and the representation mismatch

At one degree-four node define the three opposite-pair sums

\[
 p_1=M_{12}+M_{34},\quad
 p_2=M_{13}+M_{24},\quad
 p_3=M_{14}+M_{23},
 \qquad M_{ab}=Z_aZ_b.                                     \tag{AP13}
\]

In the strict lock, `p_1+p_2+p_3=-2`.  Their centered plane is the local
two-dimensional `E` representation.  Equivalently,

\[
 \mathbb R^{\{ab\}}\cong A_1\oplus E\oplus T_2,
 \qquad
 P_Q\,\delta M\,P_Q=P_E\,P_Q\delta M P_Q,                 \tag{AP14}
\]

with `A1` fixed and `T2` zero.  The six `k=2` spin assignments give only
three distinct pair vectors because particle-hole complementation leaves
every `M_ab` invariant.  They are the three vertices of an `E`-plane
triangle, not a small continuous vector supplied by the lock itself.

The exact `S4` character inner product is

\[
 \dim\operatorname{Hom}_{S_4}(T_2,E)=0.                    \tag{AP15}
\]

Therefore there is no background-independent `S4`-equivariant linear map
identifying the trivial-character link kernel with the local locked pair
sector.  At a nontrivial character, the two-dimensional kernel in (AP11) is
defined relative to the moving character direction and is covariant over the
character orbit; it is not a fixed local `E` irrep.  This closes the direct
argument

\[
 \text{quadratic incidence-Gram softness}
 \ \not\Longrightarrow\ 
 \text{soft local pair-}E\text{ response}.                 \tag{AP16}
\]

## 4. Native loop motion supplies no E Ward law

Let `T_l` toggle an alternating locked hexagon `l`.  It preserves every
degree constraint and hence acts within `P_Q`, but it need not preserve the
uniform sum of the local pair-`E` variables.

The verifier constructs a degree-two configuration on the same declared
`Q_4`, with the target hexagon alternating, for which the counts of the three
local opposite-pair types on the six hexagon vertices are

\[
 (N_1,N_2,N_3)_i=(1,2,3),
 \qquad
 (N_1,N_2,N_3)_f=(3,2,1).                                  \tag{AP17}
\]

Consequently

\[
 (\mathbf N_f-\mathbf N_i)=(2,0,-2)\ne0.                   \tag{AP18}
\]

Since the off-diagonal matrix element (AP04) is nonzero,

\[
 \langle f|[\mathbf N,H_{\rm eff}^{(6)}]|i\rangle
 =(\mathbf N_f-\mathbf N_i)
 \langle f|H_{\rm eff}^{(6)}|i\rangle\ne0.                \tag{AP19}
\]

Let `N_E` be the centered part of the three-component count.  The required
symmetry is exact on `Q_4`: with port displacements
`d_1,d_2,d_3,d_4=(e_1,e_2,e_3,0)`, every port permutation `sigma` is
implemented modulo four by the invertible map
`A_sigma d_j=d_{sigma(j)}-d_{sigma(4)}`, together with the corresponding
child-sublattice shift.  The verifier checks all 24 graph automorphisms.
Consequently the space

\[
 \{w\in E^*:[w\mathbin\cdot\mathbf N_E,H_{\rm eff}^{(6)}]=0\}
                                                                    \tag{AP19a}
\]

is `S4` invariant because the Hamiltonian and count transform covariantly.
Since `E` is irreducible, this space is either zero or all of
`E^*`.  The nonzero witness (AP19) excludes the latter.  Hence **no nonzero
uniform linear combination of pair-`E` is conserved** by the native loop
dynamics.  Particle-hole symmetry does not repair this: it acts trivially on
`M_ab`.  GL6AN's linear-degree no-go and (AP19) therefore provide no native
continuous charge whose Ward identity forces a zero-frequency `E` response.

This does not exclude nonlocal winding labels of a truncated local-loop model
or an emergent conservation law in a separately established phase.  It says
that neither is supplied by the sealed premise.

## 5. Exact covariance and the symmetry-complete quadratic form

The fixed `Q_4` has a discrete character set and by itself has no literal
`theta -> 0` infrared limit.  From this section onward, impose the additional
premise of a translation-covariant growing-quotient or thermodynamic
completion that agrees with the local locked rules.  Neither existence of
that completion nor analyticity of its response is inferred from `Q_4`.

Choose, as an additional state premise, a translation- and `S4`-invariant
stationary state of the locked effective dynamics, and let
`G_E^R(omega,chi)` be its connected retarded two-by-two pair-`E` response.
For every `sigma in S4`, exact covariance is

\[
 G_E^R(\omega,\sigma\chi)
 =D_E(\sigma)G_E^R(\omega,\chi)D_E(\sigma)^T.              \tag{AP20}
\]

At `chi=1`, Schur's lemma gives

\[
 G_E^R(\omega,1)=g_E^R(\omega)I_2.                         \tag{AP21}
\]

For the spatial expansion, the centered character tangent `theta` transforms
as `T2`.  Define

\[
 \begin{aligned}
 r_1(\theta)&=\theta_1\theta_2+\theta_3\theta_4,\\
 r_2(\theta)&=\theta_1\theta_3+\theta_2\theta_4,\\
 r_3(\theta)&=\theta_1\theta_4+\theta_2\theta_3.
 \end{aligned}                                             \tag{AP22}
\]

Their sum is `-I_2(theta)/2`, where

\[
 I_2(\theta)=\sum_a\theta_a^2.                             \tag{AP23}
\]

The centered triple

\[
 Q_E(\theta)=\left(r_1,r_2,r_3\right)
 -{r_1+r_2+r_3\over3}(1,1,1)                               \tag{AP24}
\]

is the unique quadratic `E` harmonic up to normalization.  Representation
algebra gives

\[
 \operatorname{Sym}^2(E)=A_1\oplus E,
 \qquad
 \operatorname{Sym}^2(T_2)=A_1\oplus E\oplus T_2.         \tag{AP25}
\]

There are exactly two invariant spatial quadratic contractions and no
invariant term linear in `theta`.  Let
`mathcal T:E -> Sym_0(E)` be the unique equivariant isomorphism up to scale.
If the inverse response exists and is analytic in character near `chi=1`,
its most general form through quadratic spatial order is

\[
 \boxed{
 \Gamma_E^R(\omega,\theta)
 =\big[a_0^R(\omega)+c_0^R(\omega)I_2(\theta)\big]I_2
 +c_2^R(\omega)\mathcal T(Q_E(\theta))
 +O(|\theta|^4).}                                          \tag{AP26}
\]

Here reciprocity means
`Gamma_E^R(omega,-theta)=Gamma_E^R(omega,theta)^T`.  The symmetric internal
matrices transform as `A1+E` and are even, while the antisymmetric internal
matrix transforms as `A2` and is odd.  Exact character algebra gives no `A2`
in either `T2` or `Sym^3(T2)`.  The remainder in (AP26) therefore starts at
fourth order for a reciprocal analytic kernel.  Without that extra premise,
(AP20), rather than an even Taylor series, is the exact statement.

For an isolated time-reversal-invariant phase with an analytic nondissipative
low-frequency expansion, (AP26) reduces to

\[
 \Gamma_E^R(\omega,\theta)
 =\big[r_E-Z_E(\omega+i0)^2+c_0I_2(\theta)\big]I_2
 +c_2\mathcal T(Q_E(\theta))+\cdots .                      \tag{AP27}
\]

For a generic stationary state, `a_0^R` may instead contain damping,
thresholds, or branch cuts.  Those are not excluded by `S4`.

Most importantly, the scalar mass `r_E` is symmetry allowed.  In addition,

\[
 \operatorname{Sym}^3(E)\supset A_1,                       \tag{AP28}
\]

so an `E`-cubic invariant is allowed in a nonlinear effective functional.
Neither a massless point nor a continuous transition is symmetry protected.

## 6. Conditional gapped, gapless, and pole criteria

Within the analytic nondissipative quadratic hypothesis (AP27), let
`rho(hat theta)` be the absolute eigenvalue of
`mathcal T(Q_E(theta))/I_2(theta)` in a fixed normalization.  The two inverse
eigenresponses are

\[
 \gamma_\pm(\omega,\theta)
 =r_E-Z_E(\omega+i0)^2
 +I_2(\theta)\big[c_0\pm c_2\rho(\hat\theta)\big]+cdots . \tag{AP29}
\]

This yields only conditional alternatives:

1. If `r_E>0`, `Z_E>0`, and no lower continuum or zero of the exact inverse
   response intervenes, the quadratic `E` response is gapped.
2. If `r_E=0`, `Z_E>0`, and
   `c_0>|c_2| max rho` with a regular nonzero residue, then

   \[
    \omega_\pm(\theta)^2
    ={I_2(\theta)\over Z_E}
      [c_0\pm c_2\rho(\hat\theta)]+o(|\theta|^2),           \tag{AP30}
   \]

   which is a conditional linearly soft character dispersion.
3. If `r_E<0`, the symmetric quadratic point is unstable.  Because `E` is a
   representation of a finite group and the cubic invariant is allowed, a
   broken discrete phase has no symmetry-mandated Goldstone mode.
4. If the exact retarded kernel is damped or nonanalytic, `r_E=0` need not
   produce an isolated pole; it may produce a continuum, overdamped response,
   or a branch point.

There is also an exact spectral formulation that does not assume (AP27).
For a chosen ground-state GNS representation, a normalized `E` covector `u`,
and character `chi`, let `mu_(u,chi)` be the positive-frequency connected
correlation measure on `[0,infinity)` and define

\[
 \Delta_E(\chi;u)=\inf\big(\operatorname{supp}\mu_{u,\chi}
 \setminus\{0\}\big).                                     \tag{AP31}
\]

A uniform positive lower bound over normalized covectors with strictly
positive inelastic weight, near the trivial character, is an `E` response
gap.  Gaplessness requires a sequence `chi_j -> 1` with nonzero inelastic `E`
weight and `Delta_E(chi_j;u_j)->0`.  An isolated `E` pole at each `chi_j` is
stronger: the measure must contain an atom

\[
 Z_j\,\delta(\nu-\Omega_j),
 \qquad Z_j>0,\qquad \Omega_j\to0.                         \tag{AP32}
\]

or an explicitly controlled asymptotically sharp equivalent.  Calling this
a robust infrared pole with nonvanishing `E` visibility additionally requires
`liminf_j Z_j>0`; a pole sequence can mathematically have residues tending to
zero.  A vanishing threshold alone can be a continuum and is not a pole.

For a ground state, the single-mode quotient supplies a useful sufficient
test.  With

\[
 S_u^+(\chi)=\int_{(0,\infty)}d\mu_{u,\chi}(\nu),
 \qquad
 f_u^+(\chi)=\int_{(0,\infty)}\nu\,d\mu_{u,\chi}(\nu),       \tag{AP33}
\]

one has

\[
 \Delta_E(\chi;u)\le {f_u^+(\chi)\over S_u^+(\chi)},
 \qquad S_u^+(\chi)>0.                                    \tag{AP34}
\]

The restriction to `(0,infinity)` removes any elastic zero-frequency atom;
without it the quotient need not bound the lowest positive support point.
Hence `f_u^+/S_u^+ -> 0` with nonzero inelastic weight is a sufficient
gaplessness criterion.  A conservation law could force such a numerator to
vanish, but
(AP19) proves that native hexagon dynamics does not conserve uniform pair
`E`.  No `O(|theta|^2)` numerator, divergent susceptibility, vanishing mass,
or atomic spectral weight follows from GL6AN alone.

At finite temperature or in a merely stationary invariant state, Liouvillian
transition frequencies and detailed-balance weights differ from the
ground-state formulation.  GL6AK supplies existence of invariant states, not
the state premise needed to decide (AP31)--(AP34).

## 7. Character versus physical momentum

The translation character is an exact label of the authenticated `A3`
automorphism.  It becomes physical momentum only after all of the following
additional gates are supplied:

1. a physical embedding `X:A3 -> R^3` and a calibrated length scale;
2. evidence that the authenticated translations act as physical spatial
   translations in the prepared system;
3. a controlled growing-region/refinement family defining
   `chi(x)=exp(i k dot X(x))` with a stable metric and units;
4. a selected stationary phase whose response has a stable threshold or
   dispersion under that limit.

Even those gates would not establish a common physical cone.  That further
requires calibrated time and agreement of the limiting propagation law with
independently sourced sectors.  No such calibration is present here.

## 8. The theorem/no-go conclusion

The sealed lock supplies an exact divergence constraint, a rank-two generic
character kernel, an exact local pair-`E` filter, and a nonzero native loop
move.  It does **not** supply the missing bridge among them:

\[
\boxed{
 \text{link }T_2\text{ transverse constraint}
 \not\equiv \text{ local pair }E,
 \qquad
 [H_{\rm eff}^{(6)},\mathbf N_E]\ne0,
 \qquad r_E\text{ is allowed}.}                            \tag{AP35}
\]

Thus the bare uniform pair-`E` count fails already at sixth order as a
conserved charge of the formal perturbative series; no all-orders effective
Hamiltonian is being asserted here.  Therefore an `E`-sector pole,
gaplessness, or hydrodynamic law remains a property of a selected phase/state
and the complete effective operator.  Within the GL6AN-only premise, the
shortest decisive gate is to supply an independently verified full connected
sixth-order operator, select or construct a controlled thermodynamic state,
and evaluate the matrix spectral threshold (AP31) together with the
single-mode quotient (AP34) on growing authenticated regions.

## 9. Strict ceiling

Nothing in GL6AP assumes or derives a gauge phase, photon, graviton, physical
momentum, calibrated speed, Lorentz/common cone, complete stress law,
Ricci/Einstein response, gravity, or `G`.
