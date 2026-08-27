# Independent hostile audit of PMMDC

**Audit date:** 2026-08-27

**Audited lane:**
`CROSS-RFT-GRA-EW-Q4-PAIR-MEMORY-METRIC-DEFORMATION-CLOSURE-V001`

**Frozen builder bytes reviewed:**

- `THEOREM.md`:
  `495e4e99171f4e3e5809f24e5a9a5b68116e996f4f29115bd22f780127d4714e`
- `verify_q4_pair_memory_metric_deformation_closure.py`:
  `5c79c9b24121b33e908c2d21b457b31759dd54dbe395385386621cd23239aeb0`
- `DEPENDENCIES.sha256`:
  `da23908c260f2ae86250cb9991894882bfb4d48a1906d542797329437b0ce044`
- pre-audit `MANIFEST.sha256`:
  `3abed82bc4d4dbaa1bde79883fd67db644e8fc253a547c3f9e5b2b47fe8904aa`

**Verdict:**
`ACCEPT__CLEAN_EXACT_FINITE_SPATIAL_TANGENT_THEOREM__SOURCE_FREEZE_AUTHORIZED`

No builder theorem, verifier, dependency, or pre-audit manifest byte was
edited during this review.

## 1. SCGQA tangent no-go

At `p_0=1/4`, differentiating

\[
 {cal F}_{\rm SC}=\lambda^2(\operatorname{diag}p-pp^{\mathsf T})
\]

and restricting both arguments to `V=1^perp` removes the two rank-one terms
and gives

\[
 \delta{cal F}_{\rm SC}|_V
 =\lambda^2P\operatorname{diag}(\delta p)P|_V.
\]

If this quadratic form vanishes, evaluation on every root `e_a-e_b` gives
`delta p_a+delta p_b=0`.  Three distinct indices force every component to
zero, so the map is injective and has rank three.  Its trace on `V` is
`(3/4) sum_a delta p_a=0`.  Equivariance therefore identifies its image with
the `T_2` copy in
`Sym^2(V)=A_1 direct-sum E direct-sum T_2`.  Any number of independent scalar
rescalings still contributes only the single `A_1` tensor direction.  An
infinitesimal orthogonal frame rotation changes an isotropic tensor by zero.
The stated four-dimensional `A_1+T_2` ceiling and missing two-dimensional
`E` sector are exact at the declared symmetric-point tangent; the theorem
correctly makes no global no-go claim.

## 2. Joint exponential family and uniform marginals

The three independent components of `X=Ps` are degree-one Walsh characters,
while the six `Y_ab=s_a s_b` are distinct degree-two Walsh characters.  They
are linearly independent modulo constants, so the nine-parameter family is
minimal.  The fixed computational-basis query is sufficient for the diagonal
state family, and its squared classical fidelity equals squared Uhlmann
fidelity exactly.

At `theta=0`, every pair Hamiltonian is invariant under the global flip
`s -> -s`.  Full finite-parameter support then gives

\[
 E_J[s_a]=0,\qquad P_J(s_a=+1)=P_J(s_a=-1)=1/2
\]

for every finite `J`, not merely at `J=0`.  The mixed `theta-J` Fisher block
also vanishes by odd parity.  This is exact one-port matching only; it says
nothing about work, heat, controller state, reservoirs, boundaries, or stress.

## 3. Rank-six metric deformation

At the uniform point, `Cov(X)=P|_V=I_V`, the six pair characters have identity
covariance, and their means vanish.  Differentiating `E[XX^T]` with respect to
`J_ab` inserts `s_a s_b`.  In the four-spin average, only the ordered pairs
`(c,d)=(a,b)` and `(b,a)` survive, giving exactly

\[
 \partial_{J_{ab}}{cal F}_\theta|_0
 =v_av_b^{\mathsf T}+v_bv_a^{\mathsf T}.
\]

There is no missing factor of two or covariance-subtraction term.  If a
symmetric operator `q` is orthogonal to all six tensors, then
`v_a^T q v_b=0` for `a!=b`; using `sum_b v_b=0` also gives every diagonal
pairing zero.  Since the tetrahedral vectors span `V`, `q=0`.  The six tensors
therefore form a basis of the six-dimensional `Sym^2(V)`.

The map `J -> F_theta(J)` is analytic, and its derivative at zero is this
isomorphism.  The ordinary finite-dimensional inverse-function theorem thus
gives an analytic local diffeomorphism onto a neighborhood of `I_V`; after
shrinking, that neighborhood is positive definite.  This conclusion is
strictly local in one finite family.  It does not produce compactly supported
fields, a projective continuum, or a physical propagation law.

## 4. Mode and edge reconstruction checks

The six-edge permutation character agrees with `Sym^2(V)` and decomposes as
`A_1+E+T_2`, once each.  The derivative map is equivariant and hence carries
each pair-memory sector to the corresponding tensor sector.

For the explicitly unnormalized tetrahedral vertices in (EW27), the six
edge-quadratic map has determinant `-2^19`.  Direct substitution into (EW30)
recovers all three diagonal and all three off-diagonal entries with the stated
`1/16` factors.  This is an exact finite tensor reconstruction, not a proof
that those six numbers are physical spacetime intervals.

## 5. Physical and record ceilings

The record bind is correctly conditional on an independently qualified
same-parent four-port episode whose complete law is the displayed family.
Four binary ports and sixteen outcomes are not thereby four reusable
operations.  Pair correlations are not records without formation, retention,
complete-query, and lineage custody.

The theorem closes only the finite local rank obstruction for six **spatial**
deformations of one information tensor.  It does not close:

- complete-query physical localization or information-to-length scale;
- lapse, shift, a Lorentzian clock/null bind, or the four constraints;
- shared-edge gluing, common-frame transport, or continuum refinement;
- a compactly supported local-field right inverse;
- complete stress/work matching, a positive induced Ricci coefficient,
  RIEHB stationarity, Einstein dynamics, or gravity.

No direct-product spectator metric is used inside the finite model, but
physical identification of its Fisher tensor with space remains a separate
soldering law.  The theorem and disposition preserve all of these ceilings.

## 6. Replay and custody

- Exact verifier replay: `PASS 3108/3108`.
- Dependency replay: `PASS 7/7`.
- Pre-audit manifest replay: `PASS 9/9`.
- All audited source files have terminal newlines and no CR or NUL bytes.
- The theorem's normalization, representation labels, fidelity convention,
  IFT scope, and physical no-go language were inspected directly.

**Final audit disposition:**

`SCGQA_SYMMETRIC_TANGENT_CEILING_EXACT__PAIR_WALSH_FAMILY_MINIMAL_AND_ONE_PORT_MARGINALS_UNIFORM__SIX_PAIR_DERIVATIVES_BASIS_OF_SYM2V_EXACT__LOCAL_IFT_AND_EDGE_RECONSTRUCTION_EXACT__QUALIFIED_RECORD_BIND_CONDITIONAL__LOCAL_FIELD_SOLDERING_TIME_CONSTRAINTS_GLUING_REFINEMENT_STRESS_RIEHB_AND_GRAVITY_OPEN__ACCEPT`
