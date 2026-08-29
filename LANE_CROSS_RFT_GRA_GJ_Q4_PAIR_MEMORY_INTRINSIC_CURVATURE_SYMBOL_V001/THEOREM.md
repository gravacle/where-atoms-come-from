# q4 pair-memory intrinsic-curvature symbol theorem

**Lane ID:** `CROSS-RFT-GRA-GJ-Q4-PAIR-MEMORY-INTRINSIC-CURVATURE-SYMBOL-V001`

**Official short name:** `PMICS`

**Date:** 2026-08-29

**Builder status:** `MUTABLE_BUILDER_REPLAY_122_OF_122_PASS__PENDING_HOSTILE_AUDIT__NOT_SEALED`

**Claim class:** exact finite-dimensional Fourier-symbol theorem; exact
composition of EW's pair-memory metric tangent with three-dimensional
linearized intrinsic curvature; exact pure-gradient kernel and curvature-mode
census; conditional interpretation inside the adopted RGRL working theory

**Not claimed:** that scalar gamma is curvature; that every record produces
curvature; a four-dimensional Ricci scalar or Einstein equation; a solved
Hamiltonian constraint; a propagating scalar or tensor mode; a continuum
pole; a nonlinear or finite-background-curvature kernel theorem; a graviton;
F3 realization of the EW variables; empirical confirmation of RGRL; gravity,
Newton's constant, or the numerical value of `G`

## 1. Question

EW's q4 family contains six observable pair expectations

\[
 C_{ab}:=\langle Y_{ab}\rangle,
 \qquad Y_{ab}=s_as_b.                            \tag{PMICS00}
\]

At `theta=0`, the statewise identity

\[
 XX^{\mathsf T}=I_V+\sum_{a<b}(v_a\odot v_b)Y_{ab}
                                                               \tag{PMICS00a}
\]

gives the exact expectation-memory metric relation

\[
 {\cal F}_\theta(C)
 =I_V+\sum_{a<b}C_{ab}(v_a\odot v_b).             \tag{PMICS00b}
\]

Thus the six observable pair-memory coordinates deform the same rank-three
localization Fisher metric through the exact isomorphism

\[
 D_C:\mathbb R^6\overset{\cong}{\longrightarrow}
 \operatorname{Sym}^2(V),
 \qquad D_C(e_{ab})=v_a\odot v_b.                \tag{PMICS01}
\]

EW's natural parameter `J_ab` is a preparation/control coordinate and is not
itself a retained record.  At the uniform point `D_J C=I_6`, so the older
`J`-tangent matrix equals (PMICS01), but PMICS uses `C` as the physical-memory
coordinate throughout.

Uniform Fisher accumulation can change scale without producing curvature.
The first unresolved kinematic question is therefore whether a nonuniform
pair-memory field has the correct gauge quotient to produce every local
intrinsic spatial-curvature direction.

The answer is yes at every nonzero Fourier momentum.  Three pair-memory
directions are pure spatial coordinate gradients.  The remaining three are
one intrinsic scalar-curvature direction and two trace-free tidal-curvature
directions.  No graviton, particle spectrum, pole, or binding premise enters
the proof.

## 2. Frozen geometry and convention

Use EW's three-dimensional contrast space `V` and a flat reference spatial
metric on one Gaussian slice.  Equivalently, retain only the principal
two-derivative symbol in a locally frozen frame.  Let a small pair-memory field
produce

\[
 \delta h_{ij}(x)
 =\ell_F^2(D_C)_{ij}{}^e\delta C_e\,e^{ik\cdot x},
 \qquad k\ne0.                                   \tag{PMICS02}
\]

Equation (PMICS02) is first an exact candidate-metric tangent.  Inside the
adopted working theory, RGRL-B--C realizes the observable pair-memory
expectation fields `C` as local physical spatial-metric fields and joins them
to qualified retained lineages.  In this repaired typing, RGRL's older word
"controls" is instantiated by `C`, not by the preparation coordinate `J`.
RGRL is an adopted working postulate, not empirical evidence and not a
microscopic F3 derivation.

Freeze the linearized intrinsic Riemann convention

\[
 R^{(1)}_{ijkl}
 ={1\over2}\left(
 \partial_k\partial_jh_{il}
 +\partial_l\partial_ih_{jk}
 -\partial_k\partial_ih_{jl}
 -\partial_l\partial_jh_{ik}\right).             \tag{PMICS03}
\]

Let

\[
 n={k\over|k|},\qquad P=I-nn^{\mathsf T},
 \qquad H=P\delta hP\big|_{k^\perp}.             \tag{PMICS04}
\]

## 3. Theorem PMICS-1 -- exact curvature-symbol quotient

For transverse vectors `a,b in k^perp`, direct substitution of
`partial_i -> i k_i` in (PMICS03) gives

\[
 \boxed{
 R^{(1)}_{a n b n}
 ={|k|^2\over2}H_{ab}.}                           \tag{PMICS05}
\]

Because `D_C` is an isomorphism and restriction

\[
 \operatorname{Sym}^2(V)\longrightarrow
 \operatorname{Sym}^2(k^\perp),\qquad
 h\longmapsto PhP|_{k^\perp}                     \tag{PMICS06}
\]

is onto a three-dimensional target,

\[
 \boxed{\operatorname{rank}(\mathcal C_kD_C)=3,}
 \qquad
 \mathcal C_k(h):=R^{(1)}_{a n b n}.              \tag{PMICS07}
\]

The kernel of (PMICS06) consists exactly of symmetric tensors with zero
transverse--transverse block.  Every such tensor, and only such a tensor, can
be written

\[
 h_{ij}=k_i\xi_j+k_j\xi_i                        \tag{PMICS08}
\]

after an immaterial Fourier-phase convention for `xi`.  Thus

\[
 \boxed{
 \ker(\mathcal C_kD_C)
 =D_C^{-1}\!\left\{k\odot\xi:\xi\in V\right\},
 \qquad\dim\ker=3.}                              \tag{PMICS09}
\]

These are the spatial pure-gradient or linearized coordinate directions.  No
physical curvature is lost in the quotient.

The flat-reference/principal-symbol qualifier is load bearing.  On a curved
background a diffeomorphism perturbation changes coordinate components of the
background curvature by its Lie derivative.  PMICS does not claim that the
full finite-background linearized curvature operator has (PMICS09) as a
literal zero kernel; its invariant content must be stated on the gauge
quotient with those lower-derivative background terms retained.

## 4. Theorem PMICS-2 -- one scalar-curvature plus two tidal directions

The transverse symmetric block decomposes uniquely as

\[
 H={1\over2}P\,\operatorname{tr}H+H^{\rm TF},
 \qquad \dim(\operatorname{tr}H,H^{\rm TF})=(1,2). \tag{PMICS10}
\]

The trace is genuine intrinsic scalar curvature:

\[
 \boxed{
 {}^{(3)}R^{(1)}
 =|k|^2\operatorname{tr}(P\delta hP).}            \tag{PMICS11}
\]

It is not a four-dimensional Ricci scalar, a solved Hamiltonian constraint,
or evidence for a propagating scalar mode.  The two transverse trace-free
directions have zero value of (PMICS11) but nonzero intrinsic Ricci/Riemann
tidal curvature.

Equivalently, the linearized three-dimensional Einstein tensor restricted to
the transverse plane is

\[
 \boxed{
 {}^{(3)}G^{(1)}_{ab}
 ={|k|^2\over2}\left(H_{ab}
 -P_{ab}\operatorname{tr}H\right).}              \tag{PMICS12}
\]

Two-dimensional trace reversal is invertible on
`Sym^2(k^perp)`.  Hence (PMICS12) has the same rank-three quotient as
(PMICS05).

## 5. Exact FZ-direction witness

Take FZ's exact nonzero direction

\[
 r=(7,15,-17),\qquad r^2=563                     \tag{PMICS13}
\]

and the orthogonal transverse frame

\[
 u=(15,-7,0),\qquad w=r\times u=(-119,-255,-274). \tag{PMICS14}
\]

Use EW's realization

\[
 n_1=(1,1,1),\ n_2=(1,-1,-1),\
 n_3=(-1,1,-1),\ n_4=(-1,-1,1),\qquad v_a=n_a/2.
                                                               \tag{PMICS15}
\]

In edge order `(12,13,14,23,24,34)`, the three coordinates
`(u^T h u,u^T h w,w^T h w)` of `D_C delta C` are

\[
 \boxed{
 \begin{pmatrix}
 88&-88&-32&-242&-88&88\\
 -2744&3840&1496&-1496&-270&-826\\
 -132840&-44712&-32400&28290&20500&6900
 \end{pmatrix}\delta C.}                         \tag{PMICS16}
\]

The common curvature factor is `|k|^2 ell_F^2/2`, with the harmless
nonorthonormal-frame normalizations understood.  The `(12,13,23)` minor is

\[
 \boxed{-173782321152\ne0,}                       \tag{PMICS17}
\]

which is an exact rational rank-three witness at the momentum already used by
FZ.

## 6. Gamma and curvature are related but not identical

At EW's symmetric point,

\[
 -\log\gamma_Q(0,\delta C)
 ={1\over4}\|\delta C\|^2+O(\|\delta C\|^3).     \tag{PMICS18}
\]

Equation (PMICS18) follows because `D_JC=I_6` and the pair-score covariance is
`I_6` at the uniform point.  It is a local expectation-coordinate statement,
not permission to call the control `J` a record.

Thus every nonzero infinitesimal pair-memory displacement is distinguishable
to Fisher order.  Equations (PMICS07)--(PMICS11) prove that curvature depends
on the direction and spatial variation of that memory, not merely on the
scalar amount (PMICS18).  Three equally distinguishable directions can be
pure coordinate gradients at a given nonzero momentum; three quotient
directions generate intrinsic curvature.  At `k=0`, every uniform first-order
deformation has zero intrinsic curvature even though it can change the metric
or its overall scale.

Consequently this theorem does not prove that every record creates curvature.
It proves that the q4 pair-memory field has exactly the complete nonuniform
gauge-quotient intrinsic-curvature capacity required of the flat-background
or local principal symbol of a spatial metric.

## 7. Strict separation from the FZ response

FZ separately proves that one finite F3 spatial source

\[
 Q^{ij}=-2{\partial H\over\partial j_{ij}}        \tag{PMICS19}
\]

has an exact two-dimensional TT quotient with active finite response at its
two frozen samples.  That source/polarization is not a metric perturbation.
It is an incomplete spatial part of a future second derivative of the
source-dependent generating functional, before temporal, current, contact,
boundary, and seagull completion.  It may not be substituted for `delta h` in
(PMICS02).

The lawful comparison is narrower: PMICS finds the same two-dimensional TT
type inside the three-dimensional metric-curvature quotient, while FZ proves
that a finite F3 source has two active TT response images.  The no-substitution
same-parent `C/J/j/Q` join remains the next proof gate.

## 8. Physical status and next calculation

Inside adopted RGRL-B--C, PMICS is an axiomatic physical curvature-
susceptibility theorem: qualified retained pair-memory fields act on the
physical spatial metric, and their nonuniform quotient has the rank-three
curvature map above.  Outside that adopted law it is a same-family kinematic
theorem awaiting physical soldering.

The deeper F3 derivation must still:

1. identify EW's qualified pair-memory expectations, the physical F3
   deformation `j`, and the source coordinates without type substitution;
2. derive shared-edge gluing, physical transport, and controlled refinement;
3. build one complete `H[j_00,j_0i,j_ij;ports]` before projection;
4. derive the native stress Ward identity and all current/contact terms; and
5. extract the analytic `k^2` effective-action coefficient and remainder.

No graviton or protected particle pole is required.  A positive complete
Einstein--Hilbert coefficient plus EX/RIEHB gives gravity; tensor quanta are a
possible downstream linearized description.

## 9. Disposition

`EW_PAIR_MEMORY_METRIC_MAP_RANK6__NONZERO_MOMENTUM_INTRINSIC_CURVATURE_SYMBOL_RANK3__KERNEL_EXACTLY_SPATIAL_PURE_GRADIENT_RANK3__CURVATURE_QUOTIENT_TRACE1_PLUS_TIDAL2__UNIFORM_K0_CURVATURE_ZERO__GAMMA_CERTIFIES_DISTINGUISHABILITY_NOT_CURVATURE_MAGNITUDE__RGRL_PHYSICAL_INTERPRETATION_AXIOMATIC__FZ_TT_SOURCE_RESPONSE_TYPE_SEPARATE__SAME_PARENT_F3_DYNAMIC_SOURCE_WARD_STIFFNESS_JOIN_OPEN__NO_GRAVITON_OR_GRAVITY_PROMOTION`
