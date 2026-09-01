# Independent hostile audit — GL6U F3 degree-interaction q4 response V001

**Target:** `LANE_CROSS_RFT_GRA_GL6U_F3_DEGREE_INTERACTION_Q4_RESPONSE_V001/`  
**Frozen theorem SHA-256:** `3f23084cd7e5fa500d331c6e9739b6f0b3e7875bdce3b33f5369eafd54bc7965`  
**Frozen MANIFEST-file SHA-256:** `7297d5b6456c106a7267fcd8fa4d2aa09d72a19f23605e0becde1ad14b2dba04`  
**Disposition:** `PASS_AT_EXACT_FINITE_INTERACTION_OWNED_FACTORIZATION_DEFECT_AND_FULL_PAIR_RESPONSE_SCOPE__NO_COLLECTIVE_STIFFNESS_RICCI_GRAVITY_OR_G_PROMOTION`

## Parent and degree-owner audit

The frozen theorem restores an existing F3 term rather than adding a new
interaction.  BS06 owns

\[
 U_d\sum_v(d_v-d_\star)^2,
\]

and FJ15 already exposes its pairwise link coupling.  GL6T's comparator set
`U_d=0`; GL6U lawfully evaluates the same frozen star with `U_d>0`.

The independent `N=0` census gives one active parent, three guard parents,
four active children, sixteen raw links, four active star links, and twelve
blank nonedges.  For occupation word `n=(n_1,...,n_4)` and
`r=sum_a n_a`, the active parent has degree `r`, child `a` has degree `n_a`,
and each guard has degree zero.  Subtracting the common all-blank scalar gives

\[
 (r-d_\star)^2+
 \sum_a(n_a-d_\star)^2+3d_\star^2-8d_\star^2
 =r(r+1-4d_\star),
\]

so the exact Dicke gap

\[
 \delta_r=r\Delta+U_dr(r+1-4d_\star)
\]

includes both active-parent and child degree owners.  No child, guard, or
nonedge energy has been silently discarded.

## Independent finite-state replay

The audit independently double-counts the edges between adjacent Hamming
shells and recovers

\[
 \langle r+1|-h\sum_aX_a|r\rangle
 =-h\sqrt{(r+1)(4-r)}.
\]

It separately reconstructs every `X_1` and `X_1Z_2Z_3` Dicke matrix element,
including the signs

\[
 (a_0,a_1,a_2,a_3)
 =\left({1\over2},-{1\over2\sqrt6},
        -{1\over2\sqrt6},{1\over2}\right).
\]

More strongly, a fresh symbolic polynomial calculation on all sixteen link
states expands `exp(-iHs)` through fifth order with
`s=tau/hbar`.  Without using author rows it obtains

\[
 x=h\delta_1s^2+O(s^4),\qquad
 z=1-2h^2s^2+O(s^4),
\]

\[
 x-y={4\over3}h^3(2\delta_2-\delta_1)s^4+O(s^6),
\]

and, using `2 delta_1-delta_2=-2U_d`,

\[
 \boxed{y-xz^2=-{16\over3}h^3U_ds^4+O(s^6).}
\]

Reality of the Hamiltonian, blank state, and observables makes the exact
expectations even in `s`; the independent expansion also verifies all odd
coefficients through fifth order vanish.  Thus the displayed leading defect
is nonzero on a sufficiently small punctured interval whenever `hU_d!=0`.

## Pair-response and branch audit

A separate integer-matrix Pauli replay checks all thirty-six double
commutators.  Equal pair labels give `-4h(X_a+X_b)`, labels sharing endpoint
`q` give `-4hX_qZ_rZ_s`, and disjoint labels give zero.  Hence

\[
 D^{\rm KEEP}=-8hxI_6-4hyA_{L(K_4)},
\]

with exact sector factors

\[
 D_{A_1}=-8h(x+2y),\quad
 D_{E_2}=-8h(x-y),\quad
 D_{T_2}=-8hx.
\]

The independent line-graph replay supplies a full six-vector eigenbasis with
eigenvalues `4,-2,0` and multiplicities `1,2,3`.  Under the theorem's open
conditions, the leading sector coefficients are all negative.  The rational
witness `(h,Delta,U_d,d_star)=(1,13,1,2)` gives

\[
 D_{A_1}=-168s^2+O(s^4),\quad
 D_{E_2}=-{800\over3}s^4+O(s^6),\quad
 D_{T_2}=-56s^2+O(s^4).
\]

KEEP and BREAK use the same operator law.  `P_K=1` activates the inherited
transverse term in KEEP; `P_K=0` removes it in BREAK while leaving the same
diagonal degree/detuning law.  The BREAK blank is an exact eigenstate and all
pair queries commute with that diagonal block, so `D_BREAK=0` exactly.

## Physics and scope verdict

At `U_d=0` the four link Hamiltonians factorize and `y=xz^2` exactly.  The
nonzero defect proportional to the restored `U_d` is therefore a genuine
interaction-owned inter-link factorization defect.  Calling it that does not
make it a complete connected cumulant, an autonomous phase, or continuum
stiffness.

The Hamiltonian responds to physical `P_K` and occupation, not to `REC` or
semantic provenance.  An unqualified active `K` word in the same physical
state would reproduce the response.  The six `M_ab` are pair queries, not
records, and the algebraic entrance response is not a normalized physical
source/read CTP Hessian.

Consequently GL6U exactly advances the record-native dynamics from GL6T's
factorized full-rank susceptibility to a finite response carrying an inherited
interaction-owned nonfactorization component.  The `E_2` response itself is
already nonzero at `U_d=0`; what GL6U newly owns is the defect proportional to
`U_d`, not the existence of `E_2` by itself.
It does not supply six normalized physical sources, contacts, the two-time
CTP/Schur quotient, a scalable collective phase, common cone, metric solder,
continuum stiffness, Ricci response, gravity, or `G`.

One nonmaterial notation issue remains in the frozen `RESULT.md`: its summary
writes the remainder as `O(tau^6)` after displaying powers of `tau/hbar`.
The dimensionally explicit form is `O((tau/hbar)^6)=O(s^6)`, as used in the
frozen theorem.  For fixed `hbar` the two have the same asymptotic order, so
this shorthand does not alter any coefficient, sign, domain, or claim.

**Audit verdict: PASS.**
