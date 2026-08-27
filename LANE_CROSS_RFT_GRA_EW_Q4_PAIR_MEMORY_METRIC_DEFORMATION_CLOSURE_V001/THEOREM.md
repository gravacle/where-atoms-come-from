# q4 pair-memory metric-deformation closure theorem

**Lane ID:** CROSS-RFT-GRA-EW-Q4-PAIR-MEMORY-METRIC-DEFORMATION-CLOSURE-V001

**Official short name:** PMMDC

**Date:** 2026-08-27

**Builder status:** MUTABLE_PRESCREEN_READY__BUILDER_REPLAY_PASS__NOT_SEALED

**Claim class:** exact finite-dimensional same-parent exponential-family
theorem; exact symmetric-point rank no-go for the SCGQA one-body softmax
family; exact six-mode pair-memory deformation closure of one rank-three
localization Fisher tensor; exact conditional qualified-record binding

**Not claimed:** that four binary ports are EO's four complete reusable
operations; that pair correlations are records without URFT qualification;
that the family is actual-world matter; that its Fisher tensor is physical
space without localization/soldering; a time metric, lapse, shift, Lorentz
constraint, shared-face gluing, Levi--Civita transport, continuum refinement,
complete stress match, Einstein dynamics, gravity, or outcome selection

## 1. Exact question

SCGQA placed record-query gamma, state fidelity, QFI, and q4 sufficient counts
in one exact family.  At its symmetric point, however, its physical state
coordinate has only three independent directions.  RIEHB's composite-metric
gate requires the retained collective variables to span the complete physical
metric-variation space, not merely to define one positive metric.

This lane asks two focused questions:

1. Can the SCGQA one-body softmax family plus one scalar accumulation mode
   span all six symmetric deformations of its rank-three localization metric?
2. Can qualified pair memory among the same four ports supply the missing
   modes without adding new operation types or a spectator metric factor?

The first answer is no at the symmetric point and to first order.  The second
answer is yes in one exact same-family model.  The result closes a local
spatial deformation-rank obstruction, not the physical gravity theorem.

## 2. Tetrahedral contrast notation

Let

\[
 \mathbf1=(1,1,1,1)^{\mathsf T},\qquad
 P=I_4-\frac14\mathbf1\mathbf1^{\mathsf T},\qquad
 V=\mathbf1^\perp,
 \qquad v_a=Pe_a.                                  \tag{EW01}
\]

Then

\[
 \sum_av_a=0,\qquad
 v_a^{\mathsf T}v_b=\delta_{ab}-\frac14,           \tag{EW02}
\]

and the four vectors span the three-dimensional space \(V\).  For vectors
\(u,v\in V\), define the unhalved symmetric product

\[
 u\odot v:=uv^{\mathsf T}+vu^{\mathsf T}.          \tag{EW03}
\]

Let

\[
 {\cal E}_4=\{\{a,b\}:1\le a<b\le4\}               \tag{EW04}
\]

be the six unordered pairs.  Permuting the four port labels gives the natural
six-dimensional edge representation of \(S_4\).

## 3. Theorem PMMDC-1 -- the SCGQA symmetric-point tangent no-go

For the SCGQA categorical family, write its QFI as a function of the
four-outcome probability vector:

\[
 {\cal F}_{\rm SC}(p)
 =\lambda^2\left(\operatorname{diag}p-pp^{\mathsf T}\right)
 \quad\hbox{on }V.                                 \tag{EW05}
\]

At \(p_0=\mathbf1/4\), an allowed probability variation obeys
\(\delta p\in V\).  The first variation restricted to \(V\) is

\[
 \boxed{
 \delta{\cal F}_{\rm SC}\big|_V
 =\lambda^2 P\operatorname{diag}(\delta p)P\big|_V.}
                                                               \tag{EW06}
\]

The softmax Jacobian at that point is
\(D_\theta p=(\lambda/4)P\), an isomorphism on \(V\).  Thus using the
SCGQA state parameter \(\theta\) rather than \(p\) does not add or remove a
tangent direction.

The linear map

\[
 L_{\rm SC}:V\longrightarrow\operatorname{Sym}^2(V),
 \qquad
 \delta p\longmapsto
 P\operatorname{diag}(\delta p)P\big|_V             \tag{EW07}
\]

is injective, has rank three, and has trace-free image.  It therefore cannot
span the six-dimensional space \(\operatorname{Sym}^2(V)\).

### Proof

The derivative of (EW05) is

\[
 \lambda^2\left[
 \operatorname{diag}\delta p-\delta p\,p_0^{\mathsf T}
 -p_0\delta p^{\mathsf T}\right].
\]

The last two terms vanish when both arguments lie in \(V\), proving (EW06).
If \(L_{\rm SC}(\delta p)=0\), evaluation on every root
\(e_a-e_b\in V\) gives

\[
 \delta p_a+\delta p_b=0\qquad(a\ne b).             \tag{EW08}
\]

Three distinct indices force all components to vanish, so the map is
injective.  Its rank is therefore \(\dim V=3\).  Finally,

\[
 \operatorname{tr}_V L_{\rm SC}(\delta p)
 =\operatorname{Tr}\!\left(P\operatorname{diag}\delta p\right)
 =\frac34\sum_a\delta p_a=0.                        \tag{EW09}
\]

QED.

The map is \(S_4\)-equivariant.  In tetrahedral notation,

\[
 \operatorname{Sym}^2(V)=A_1\oplus E\oplus T_2,
 \qquad
 \operatorname{im}L_{\rm SC}=T_2.                  \tag{EW10}
\]

A scalar change of copy number, information-to-length conversion, or overall
QFI scale adds only the \(A_1\) direction.  Thus the SCGQA one-body coordinate
plus any number of scalar rescalings spans at most
\(A_1\oplus T_2\), of dimension four; the two-dimensional \(E\) shear sector
is missing.  Rotating a coframe does not repair this linear no-go at the
isotropic point because the first variation of an isotropic metric under a
frame rotation is zero.

This theorem is local and first-order at the symmetric point.  It does not
exclude a different family, an anisotropic-background construction, or
higher-order reopening.

## 4. One joint localization/pair-memory family

Let \(s=(s_1,s_2,s_3,s_4)\in\{-1,+1\}^4\).  Define the contrast-valued
one-port statistic and six pair statistics

\[
 X(s):=\sum_{a=1}^4v_as_a=Ps\in V,\qquad
 Y_{ab}(s):=s_as_b.                                \tag{EW11}
\]

For \(\theta\in V\) and
\(J=(J_{ab})_{\{a,b\}\in{\cal E}_4}\in\mathbb R^6\), define

\[
 p_{\theta,J}(s)
 ={1\over Z(\theta,J)}
 \exp\!\left[
 \theta^{\mathsf T}X(s)+
 \sum_{a<b}J_{ab}Y_{ab}(s)\right],                 \tag{EW12}
\]

and the commuting carrier state

\[
 \rho_{\theta,J}
 =\sum_s p_{\theta,J}(s)|s\rangle\langle s|.       \tag{EW13}
\]

Every finite-parameter state has full support.  The three independent
components of \(X\) and the six pair characters \(Y_{ab}\) are linearly
independent modulo constants: they are distinct Walsh characters on the
four-bit cube after the sole linear relation \(\sum_av_a=0\) is removed.
Hence (EW12) is a minimal nine-parameter exponential family.

The computational-basis query \(M_s=|s\rangle\langle s|\) is complete for
this diagonal family.  Therefore, for all two parameter values \(u,u'\),

\[
 \boxed{
 \gamma_Q[p_u,p_{u'}]
 =\gamma_{\rm state}(\rho_u,\rho_{u'})
 =\left(\sum_s\sqrt{p_u(s)p_{u'}(s)}\right)^2.}     \tag{EW14}
\]

The value is strictly below one for distinct parameters.  The SLD QFI equals
the classical Fisher covariance of the complete sufficient-statistic vector
\((X,Y)\).

The family is exactly \(S_4\)-covariant.  A port permutation acts on \(s\),
\(\theta\), \(J\), and the basis labels simultaneously; it preserves
\(\theta\cdot X+\sum J_{ab}Y_{ab}\) and conjugates (EW13) by the corresponding
permutation unitary.

## 5. Theorem PMMDC-2 -- pair memory deforms the same localization metric

Set \(\theta=0\).  Global spin flip \(s\mapsto-s\) leaves every pair statistic
and (EW12) invariant.  Thus, for every finite \(J\),

\[
 \boxed{
 \mathbb E_J[s_a]=0,\qquad
 \Pr_J(s_a=+1)=\Pr_J(s_a=-1)=\frac12.}             \tag{EW15}
\]

The one-port marginals and occupancies remain exactly fixed while pair
memory changes.

The Fisher tensor for the localization coordinate \(\theta\) is

\[
 {\cal F}_\theta(J)
 =\operatorname{Cov}_J(X)
 =\mathbb E_J[XX^{\mathsf T}]
 \quad\hbox{on }V.                                 \tag{EW16}
\]

At \(J=0\), the four spins are independent and uniform, so

\[
 \boxed{{\cal F}_\theta(0)=P\big|_V=I_V.}          \tag{EW17}
\]

The pair-coupling Fisher block is

\[
 {\cal F}_{J,ef}(J)=\operatorname{Cov}_J(Y_e,Y_f).
                                                               \tag{EW18}
\]

At \(J=0\), distinct Walsh pair characters are orthogonal:

\[
 \boxed{{\cal F}_J(0)=I_6.}                        \tag{EW19}
\]

The mixed \(\theta\)-\(J\) Fisher block vanishes at every
\((\theta=0,J)\), because \(X\,Y_e\) is odd under global flip.

Most importantly,

\[
 \boxed{
 {\partial{\cal F}_\theta\over\partial J_{ab}}\bigg|_{J=0}
 =v_a\odot v_b.}                                   \tag{EW20}
\]

The six matrices \(v_a\odot v_b\), \(a<b\), form a basis of
\(\operatorname{Sym}^2(V)\).  Hence

\[
 D_J{\cal F}_\theta(0):
 \mathbb R^6\overset{\cong}{\longrightarrow}
 \operatorname{Sym}^2(V)                           \tag{EW21}
\]

is an exact linear isomorphism.

### Proof

At the uniform point, differentiating an expectation with respect to
\(J_{ab}\) inserts the centered score \(Y_{ab}\), whose mean is zero.
Expanding \(XX^{\mathsf T}\) using (EW11), the uniform expectation of
\(s_cs_ds_as_b\) vanishes unless \((c,d)=(a,b)\) or \((b,a)\).  The two
surviving terms give (EW20).

To prove spanning, suppose a symmetric operator \(q\) on \(V\) is Frobenius-
orthogonal to all six matrices.  Then

\[
 v_a^{\mathsf T}qv_b=0\qquad(a\ne b).              \tag{EW22}
\]

Using \(\sum_bv_b=0\),

\[
 v_a^{\mathsf T}qv_a
 =-\sum_{b\ne a}v_a^{\mathsf T}qv_b=0.             \tag{EW23}
\]

All matrix elements of \(q\) on the spanning tetrahedral frame therefore
vanish, so \(q=0\).  The six matrices are linearly independent in a
six-dimensional target and form a basis.  QED.

Because \(J\mapsto{\cal F}_\theta(J)\) is analytic and its derivative at the
origin is invertible, the inverse-function theorem gives neighborhoods
\({\cal U}\subset\mathbb R^6\) of \(0\) and
\({\cal W}\subset\operatorname{Sym}^2(V)\) of \(I_V\) such that

\[
 J\longmapsto{\cal F}_\theta(J)                    \tag{EW24}
\]

is an analytic local diffeomorphism from \({\cal U}\) onto \({\cal W}\).
After shrinking the neighborhoods if necessary, every tensor in
\({\cal W}\) is positive definite.  Thus one same-parent pair-memory family
realizes every sufficiently small symmetric deformation of its own
rank-three localization Fisher metric.  No direct-product metric spectator
is used.

## 6. Exact \(S_4\) mode census

For an edge vector \(x=(x_{ab})\), define opposite-pair sums and differences

\[
 \begin{aligned}
 p&=(x_{12}+x_{34},\ x_{13}+x_{24},\ x_{14}+x_{23}),\\
 d&=(x_{12}-x_{34},\ x_{13}-x_{24},\ x_{14}-x_{23}).
 \end{aligned}                                     \tag{EW25}
\]

The uniform part of \(p\) is \(A_1\); the subspace
\(p_1+p_2+p_3=0\) is the two-dimensional \(E\) sector; and \(d\) is the
three-dimensional \(T_2\) sector.  Therefore

\[
 \boxed{\mathbb R^{{\cal E}_4}=A_1\oplus E\oplus T_2.}
                                                               \tag{EW26}
\]

The same decomposition holds for \(\operatorname{Sym}^2(V)\).  Equation
(EW20) is \(S_4\)-equivariant and invertible, so it carries each edge sector
onto the corresponding spatial-metric sector.  Pair memory supplies exactly
the \(E\) modes absent from PMMDC-1, while retaining the scale and \(T_2\)
modes in one six-coordinate family.

## 7. Explicit edge-quadratic reconstruction

Choose an isometry \(O:V\to\mathbb R^3\) and set \(n_a=2Ov_a\).  A
tetrahedral coordinate realization is

\[
 \begin{aligned}
 n_1&=(1,1,1),&n_2&=(1,-1,-1),\\
 n_3&=(-1,1,-1),&n_4&=(-1,-1,1).
 \end{aligned}                                     \tag{EW27}
\]

For \(q=q^{\mathsf T}\in\operatorname{Sym}^2(\mathbb R^3)\), define six edge
quadratic forms

\[
 y_{ab}=(n_a-n_b)^{\mathsf T}q(n_a-n_b).           \tag{EW28}
\]

In the coordinate order
\((q_{xx},q_{yy},q_{zz},q_{xy},q_{xz},q_{yz})\) and edge order
\((12,13,14,23,24,34)\), this linear map has determinant

\[
 \det{\cal Q}_{\rm edge}=-2^{19}\ne0.              \tag{EW29}
\]

Its inverse is explicit.  Put

\[
 A=y_{14}+y_{23},\qquad
 B=y_{13}+y_{24},\qquad
 C=y_{12}+y_{34}.                                  \tag{EW30}
\]

Then

\[
 \begin{aligned}
 q_{xx}&={A+B-C\over16},&
 q_{yy}&={A+C-B\over16},&
 q_{zz}&={B+C-A\over16},\\
 q_{xy}&={y_{14}-y_{23}\over16},&
 q_{xz}&={y_{13}-y_{24}\over16},&
 q_{yz}&={y_{12}-y_{34}\over16}.
 \end{aligned}                                     \tag{EW31}
\]

Thus the six pair-memory variations first deform
\({\cal F}_\theta\) through the isomorphism (EW21), and the six edge
quadratic forms of that same tensor determine all of its spatial components
through (EW31).  Equation (EW31) is a reconstruction theorem.  Calling the
reconstructed tensor physical space still requires the same-query
localization and scale premises isolated by GSGB/FERS.

## 8. Product accumulation and conditional qualified-record bind

For \(N\) actual independent copies of (EW13), the complete product query has

\[
 \gamma_N=\gamma_1^N,\qquad
 {\cal F}_{\theta,N}=N{\cal F}_\theta,\qquad
 {\cal F}_{J,N}=N{\cal F}_J.                       \tag{EW32}
\]

The accumulated statistics
\(\sum_rX(s^{(r)})\) and \(\sum_rY_{ab}(s^{(r)})\) are sufficient for the
product exponential family.  This is conditional on actual tensor-product
copies; record redundancy alone does not imply it.

To bind the model to record formation, supply one authenticated same-parent
four-port episode whose retained binary marks, complete query, formation
margin, reference-stable descent, and lineage custody pass the adopted URFT
predicates.  Require its actual joint query law and commuting carrier state to
be (EW12)--(EW13), with every alternative-dependent query port included in
the census.  Under that supplied bind, pair-memory alternatives have the same
query gamma and state gamma (EW14), and all six \(J\) directions are
distinguishable at positive Fisher order.

At \(\theta=0\), a KEEP/BREAK comparison can vary pair memory while preserving
all one-port marginals exactly by (EW15).  This is a strong one-body match,
not a complete stress match.  Pair preparation, controller, work, heat,
reservoir, boundary, failure, and quarantine variables may still differ and
must be exhaustively matched or rendered terminal/reference-stably blind
before any response is attributed specifically to lineage.

## 9. Exact RIEHB interface and remaining physical gates

PMMDC closes one algebraic part of RIEHB's composite-metric rank problem:
the same joint family has one rank-three localization Fisher tensor and six
pair-memory controls whose derivative spans every symmetric deformation of
that tensor near the isotropic point.

The following implications are not licensed by this theorem:

1. **q4 operation ancestry.**  Four binary query ports or sixteen joint
   outcomes are not four complete reusable operations and do not establish
   Q4-MERGE.
2. **Physical localization.**  It remains to prove that clocks, matter,
   electromagnetism, and both probes factor through this complete query and
   read \(q=\ell_F^2{\cal F}_\theta(J)\).
3. **Time and constraints.**  Six spatial symmetric variations do not by
   themselves supply lapse, shift, a Lorentzian clock/null bind, causal
   constraints, or dense span of the full four-metric variation quotient.
4. **Gluing and refinement.**  Shared edges, faces, orientation, physical
   Levi--Civita probe transport, and a shape-regular continuum limit remain
   unproved.
5. **Dynamics.**  A physical law must make the retained collective variables
   \(J=J(\Phi)\) respond and propagate without inserting desired edge values.
   The exact chain still requires the actual \(Dg\) and explicit-force gates.
6. **Complete stress and induced action.**  The common-cone fast spectrum,
   positive matched Ricci coefficient, complete stress census, controlled
   remainders, and RIEHB hypotheses remain separate.

If these same-world gates are later earned, (EW21) supplies the missing local
spatial tangent rank rather than requiring six tensor bits in each individual
record.  Pair memory is then a collective deformation coordinate; classical
curvature and gravitational back-reaction remain macroscopic outputs.

## 10. Controls and exact disposition

1. **Symmetric-point scope.**  Promote PMMDC-1 to a global no-go.  Rejected:
   only the first derivative at the isotropic SCGQA point is excluded.
2. **Scalar repair.**  Use multiple scalar accumulation factors to fill the
   \(E\) sector.  Rejected: all such variations remain \(A_1\).
3. **Outcome/operation confusion.**  Call four ports or sixteen joint
   outcomes Q4-MERGE.  Rejected.
4. **Marginal/stress confusion.**  Infer equal energy or stress from the
   uniform one-port marginals (EW15).  Rejected.
5. **QFI typing.**  Call the six-by-six \(J\)-Fisher block a three-dimensional
   spatial metric.  Rejected: the spatial candidate is
   \({\cal F}_\theta(J)\); \({\cal F}_J\) certifies distinguishable
   deformation coordinates.
6. **Metric/spacetime confusion.**  Promote the Fisher tensor to physical
   space without complete-query localization and scale.  Rejected.
7. **Spatial/full-metric confusion.**  Promote six spatial deformation modes
   directly to the full RIEHB dense-range gate.  Rejected.
8. **Recordhood shortcut.**  Treat a pair correlation as a record without
   formation, retention, complete query, and lineage custody.  Rejected.
9. **Gravity shortcut.**  Relabel deformation closure as curvature, Einstein
   response, \(G\), or gravity.  Rejected.

**Disposition:**

SCGQA_ONE_BODY_SOFTMAX_PLUS_SCALAR_SCALE_SPANS_ONLY_A1_PLUS_T2_AT_THE_SYMMETRIC_POINT__SIX_Q4_PAIR_MEMORY_COUPLINGS_FORM_A1_PLUS_E_PLUS_T2__JOINT_EXPONENTIAL_FAMILY_HAS_UNIFORM_ONE_PORT_MARGINALS_FULL_PAIR_FISHER_RANK_AND_EXACT_QUERY_STATE_GAMMA_JOIN__PAIR_COUPLINGS_DEFORM_THE_SAME_RANK3_LOCALIZATION_FISHER_METRIC_WITH_FULL_SYM2V_TANGENT_RANK__EDGE_QUADRATIC_MAP_TO_ALL_SIX_SPATIAL_COMPONENTS_INVERTIBLE__QUALIFIED_RECORD_BIND_CONDITIONAL__Q4_MERGE_PHYSICAL_SOLDERING_TIME_CONSTRAINTS_GLUING_REFINEMENT_COMPLETE_STRESS_RIEHB_AND_GRAVITY_OPEN__MUTABLE_NOT_SEALED
