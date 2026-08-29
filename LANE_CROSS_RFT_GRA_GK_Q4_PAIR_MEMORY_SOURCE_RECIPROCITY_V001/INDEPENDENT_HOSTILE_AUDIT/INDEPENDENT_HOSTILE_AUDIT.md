# Independent hostile audit -- PMSR

**Audited lane:** `CROSS-RFT-GRA-GK-Q4-PAIR-MEMORY-SOURCE-RECIPROCITY-V001`

**Audit date:** 2026-08-29

**Disposition:** `PASS_209_OF_209__NO_REMAINING_MATERIAL_DEFECT`

## 1. Independence and frozen scope

This audit reconstructed the q4 state algebra, covariance spectrum, source
signs, mixed derivatives, and rank statements from the sixteen spin states.
It did not import, invoke, or copy the builder verifier.  The independent
replay is `verify_hostile_pmsr.py` in this directory.

The final audited target bytes are frozen in `AUDITED_TARGETS.sha256`:

```text
9d79e42d128dd3bae72adc12c26aa27d4e4102bfec917a43d1d6c6688ea49ec7  THEOREM.md
08e0329a2c93bb8e9ace86a819b6b16d8303804004b5e06a303fd4d9a0c94163  RESULT.md
244579cc31ea06179473c7babac87444391632ca7d34290a36d371d2869c7f89  SELF_AUDIT.md
```

All ten declared dependency files exist and match their frozen hashes.

## 2. Defects found and repaired before the freeze

The hostile pass found one real scope defect in the mutable draft.  `GK19`
derives the conjugate of the commuting pair Hamiltonian `H_C`, whereas the
complete FU/F3 strain source also owns the one-edge/dressed-transfer sector
and may own other declared terms.  The draft result called
`Q=-2 partial_j H` the complete physical source.  That exceeded the proved
custody.

The final target repairs this consistently:

- `GK19` is the **pair-sector physical conjugate**
  `Q_A=-2 partial_{j_A}H_C`;
- the reciprocity theorem is explicitly a pair-sector physical-source
  statement; and
- the theorem says that the omitted complete F3 source sectors are not
  supplied.

A packaging warning caused by a trailing blank line in
`DEPENDENCIES.sha256` was also removed.  Neither repair changes the exact
q4 algebra.  No material defect remains in the frozen target.

## 3. Exact algebraic replay

For every one of the sixteen states, the audit independently recovered

\[
 X X^{\mathsf T}=P+\sum_eY_eB_e.
\]

It then proved by exact elimination that the six `B_e` span
`Sym^2(1^perp)`.  Taking expectations at `theta=0` therefore gives the exact
affine identity

\[
 {\cal F}_\theta=I_V+{\cal B}C,
 \qquad C_e=\langle Y_e\rangle,
\]

with `D_C F=B` invertible.  This directly validates the lane's essential
typing: `C` is the observable expectation coordinate; `J` remains the
natural control.  The audit separately proved that the constant and six
pair Walsh characters are independent.  Hence, at finite full support,

\[
 D_JC=G=\operatorname{Cov}(Y,Y)>0.
\]

Away from the uniform point `G` is not the identity, so replacing `C` by
`J` would be a genuine mathematical error, not a harmless relabeling.

The six normalized tetrahedral root dyads were rebuilt directly and their
strain-evaluation matrix `M` has exact nonzero determinant.  Under the
explicit whole-pair-source premise `GK-S4`, exact composition consequently
gives

\[
 D_j{\cal F}_\theta
 ={\beta U_d\lambda_{\rm pair}^{\rm net}\over2}{\cal B}GM.
\]

Because all three maps are invertible at finite full support, the rank is
six.  The stronger `GK-S4` premise is necessary: FU's nonzero `E` component
alone does not prove that the complete six-direction tangent is one scalar
multiple of `M`.  For a general audited tangent `L`, the correct formula is
`(beta U_d/2) B G L`, and its rank equals the rank of `L`.  The independent
replay checked both an invertible and a rank-five `L`.

## 4. Source convention and reciprocity

The audit used a negative test slope so that sign errors could not be hidden
by a positivity convention.  It independently recovered

\[
 \partial_jK=-{U_d\lambda\over2}M,
 \qquad
 \partial_jJ=-\beta\partial_jK
 ={\beta U_d\lambda\over2}M,
\]

and

\[
 Q_A=-2\partial_{j_A}H_C
 =U_d\lambda\sum_eM_{eA}Y_e+q_{{\rm id},A}I.
\]

Exact state sums for every source direction, including a nonzero identity
source control, then reproduced

\[
 \partial_{j_A}{\cal F}_{\theta,mn}
 ={\beta\over2}\partial_{\theta_m}\partial_{\theta_n}
 \langle Q_A\rangle.
\]

For squared complete-query fidelity, the audit separately used

\[
 \gamma(t)={Z(t/2)^2\over Z(0)Z(t)}
\]

along an arbitrary localization direction.  Its second derivative gives the
coefficient `F/4` in `-log gamma`; differentiating that coefficient with
respect to all six source directions exactly reproduced the factor
`beta/8` in `GK24`.

## 5. Thermal rank and boundary controls

Direct sector enumeration gives the partition factor

\[
 D=a^4+4a+3,
 \qquad
 c={a^4-1\over D},
 \qquad
 q={a^4-4a+3\over D}.
\]

For the edge covariance, the independently reduced eigenvalues are

\[
 g_{A_1}=1+4c+q-6c^2
 ={8a(3a^4+4a^3+1)\over D^2},
\]

\[
 g_E=1-2c+q={8\over D},
 \qquad
 g_{T_2}=1-q={8a\over D}.
\]

They are strictly positive for every finite `a>0`, with multiplicities
`1`, `2`, and `3`.  Exact enumeration at five distinct rational values of
`a` replayed every eigenvector and formula.

All three singular controls also pass:

1. at `beta=0`, `G=I_6` but `D_jJ=0`, so the state/Fisher response is zero
   even when the source operator exists;
2. at zero complete pair slope, the nonidentity pair source and response
   vanish; and
3. on the exact six-state ice fiber, `G` and `BGM` have rank two, solely in
   the `E` sector.

These controls rule out the false inference that operator availability by
itself supplies a state response or a six-channel projected memory metric.

## 6. FY and physics ceiling

Dependency inspection confirms only direct pair-source **operator custody**
through FV/FY.  FY does not insert the EW `theta` source, transport the EW
complete query through its noncommuting Feshbach reduction, or equate its
response metric to the commuting squared-fidelity metric.  The theorem
correctly notes that a noncommuting Gibbs Hessian is generally the
Bogoliubov--Kubo--Mori metric and is not automatically the SLD/QFI metric.

Accordingly, PMSR proves an exact microscopic pair-memory/source reciprocity
bridge under its stated same-parent solder.  It does not prove full F3 source
completeness, physical spacetime soldering, a Ward identity, a massless pole,
Einstein dynamics, gravity, or a numerical `G`.

## 7. Final result

The independent replay passes **209/209** exact checks.  The stronger
whole-pair-source premise, the general-`L` fallback, the finite-temperature
rank theorem, all boundary controls, the FY custody ceiling, and the frozen
dependency chain are internally consistent.

**Final hostile disposition:**

`PASS__EXACT_PAIR_EXPECTATION_TO_FISHER_AFFINE_ISOMORPHISM__EXACT_CONDITIONAL_PAIR_SOURCE_RECIPROCITY__FINITE_BETA_RANK6__ICE_RANK2__COMPLETE_F3_AND_GRAVITY_CEILINGS_PRESERVED`
