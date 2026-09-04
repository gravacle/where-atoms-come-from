# Independent hostile audit

**Audit date:** 2026-08-27

**Auditor role:** independent adversarial replay; all builder files were
source-frozen and only this audit file was changed

**Verdict:** `PASS`

**Exact remaining correctness defects:** none found

**Scope of verdict:** GRA-DQ is proof-grade as an exact negative real-time
response theorem for the declared `h=t=0` hold, an exact static Gibbs
conditional-reciprocity calculation, and an exact finite conditional-link
Rabi calculation after restoring the inherited BS06 actuator. It does not
prove reciprocal dynamics in the sealed hold, a matched post-GRA-DN phase
comparison, future FPMH record formation, a nonzero fixed-time local
thermodynamic response, or gravity.

## 1. Frozen replay and custody

The audit read every packet file and the pinned BS, GRA-DJ, and GRA-DN
dependencies. Before replacing the pending audit:

- `verify_reciprocal_backpropagation_screen.py` replayed as `PASS 204/204`;
- all six dependency digests matched; and
- all nine pre-audit packet-manifest entries matched.

The executable reconstructs the four-factor hold, conserved operators,
degree increment, Gibbs conditional, static covariance, Rabi transition,
arbitrary carrier-state reduction, short-time series, and finite-size bound.
Those checks corroborate rather than replace the analytic and scope review
below.

## 2. Exact real-time obstruction

In the declared slice,

\[
 H_{\rm hold}=U_d\sum_v[d_v(n)-4]^2
 +\sum_e n_e[a_r+\lambda_J C_e],
 \qquad C_e=q_u+q_v-2q_uq_v .
\]

All incidence occupations commute with this Hamiltonian. At carrier hopping
`t=0`, every `C_f` is also a polynomial in mutually commuting carrier
occupations, and the retained FPMH `K` factor is absent from the Hamiltonian.
Thus

\[
 [H_{\rm hold},n_e]=[H_{\rm hold},C_f]
 =[H_{\rm hold},P_i^K]=0.
\]

It follows exactly that `n_e(tau)=n_e` and `P_i^K(tau)=P_i^K`. Because every
carrier-only perturbation observable acts on a tensor factor distinct from
`n` and `K`, both retarded commutators vanish for every state and time. This
is not an equilibrium or perturbative approximation. Each finite-volume
kernel is identically zero, so every thermodynamic subsequential limit is
also zero.

The conclusion is correctly limited to the declared hold. In particular,
the current-square term changes an energy conditional on incidence but
contains no incidence-changing operator; it cannot substitute for the BS06
`X_e` actuator. The source also correctly refrains from inferring a change in
a future formation query while the FPMH writer/source ports are off and `K`
is conserved.

## 3. Static reciprocity is not causal response

With all other incidences frozen and endpoint degrees excluding the candidate
edge equal to `d_u,d_v`, direct subtraction gives

\[
 (d+1-4)^2-(d-4)^2=2d-7
\]

and therefore

\[
 \Delta_eE(c)=A_e+\lambda_Jc,
 \qquad A_e=a_r+2U_d(d_u+d_v-7).
\]

The displayed logistic Gibbs conditional follows exactly. For positive
`lambda_J`, aligned endpoints have the larger conditional incidence
probability on the same frozen background. For a prescribed carrier state,
only `p_sigma=Tr(sigma C_e)` enters the mean link-addition energy; an
"ordered" label or magnetization alone does not determine it.

For the explicitly stated plus-source convention
`H -> H+xi n_e+zeta C_f`, differentiation of the finite free energy gives

\[
 \partial_\xi\partial_\zeta F
 =-\beta\operatorname{Cov}(n_e,C_f).
\]

This is a symmetric static Hessian of one Gibbs endpoint. The packet
correctly does not promote it to a retarded kernel, a time arrow, a
transition rate, or dynamic back-propagation.

## 4. Restored-BS06 Rabi law

After restoring only the inherited conditional link flip and freezing every
other contribution, the two-level block is

\[
 H_e(c)=-h_NX_e+F_cn_e,
 \qquad F_c=A_e+\lambda_Jc,
 \qquad h_N={\Omega\over N}.
\]

Removing an irrelevant scalar and exponentiating the remaining Pauli block
gives

\[
 R(F,\tau)={4h_N^2\over F^2+4h_N^2}
 \sin^2\!\left({\tau\sqrt{F^2+4h_N^2}\over2\hbar}\right).
\]

Since `C_e` is a conserved projector, an initially blank pure link is
factorized from an arbitrary prescribed endpoint carrier density operator,
and off-diagonal coherence between the two `C_e` sectors cannot affect a
link-only final probability. Hence

\[
 P_\sigma=(1-p_\sigma)R(A_e,\tau)
 +p_\sigma R(A_e+\lambda_J,\tau)
\]

is exact, including for coherent carrier inputs. Direct matrix
exponentiation independently reproduces this law.

The short-time expansion is also exact through fourth order:

\[
 R(F,\tau)={h_N^2\tau^2\over\hbar^2}
 -{h_N^2(F^2+4h_N^2)\tau^4\over12\hbar^4}+O(\tau^6).
\]

Thus the detuning dependence is absent at quadratic order and first appears
in the fourth-order coefficient. Special detunings or pulse times can make
the two sector probabilities equal, which is why the theorem correctly calls
the finite contrast generic rather than universal.

Finally, `sin^2 x <= x^2` yields the uniform bound

\[
 0\le R(F,\tau)\le {\Omega^2\tau^2\over N^2\hbar^2}.
\]

At fixed time and fixed `Omega`, the local response is therefore
`O(N^-2)`. No collective sum, long-time scaling, relaxation law, or nonzero
thermodynamic response is silently inferred.

## 5. Common-background and GRA-DN composition ceiling

The finite Rabi contrast compares carrier preparations only when they share
the same target incidence word, `K` packet, storage occupations, endpoint
degrees, boundary, clock, and BS06 pulse. Different carrier labels must also
actually give different local `p_sigma`; magnetization alone is insufficient.

GRA-DN does not supply this intervention. Its ordered carrier state is formed
on the KEEP full-incidence sheet, while its uniform product state is formed on
the REBIND blank-incidence sheet. Their degree backgrounds and detunings are
therefore different before any reciprocal test. Substituting those two arms
directly into the Rabi formula would confound carrier input with the already
changed incidence background. A reset, transport, or re-preparation mission
and its port/work ledger remain to be supplied. The packet correctly claims
neither preparation-work equality nor a matched laboratory work
distribution.

The declared hold also supplies no new FPMH mission. Identity action on `K`
preserves its occupation/query statistics; switching a source or writer on
later would be a separately specified mission whose formation probability is
not determined by this lane. This is the correct future-record-formation
ceiling and does not amount to a no-go against all later reciprocal record
experiments.

## 6. Final claim boundary

The earned result is narrow and useful: unchanged BS11 supplies a static
carrier-conditioned incidence energy but zero causal incidence or `K`
response in the sealed hold; restored BS06 supplies the smallest exact finite
conditional dynamic observable. No autonomous full-sheet feedback,
semantic `REC` necessity, unique causation, emergent support, metric, physical
stress response, universal coupling, Newton `G`, or gravity follows.

Within that boundary, no remaining algebraic, logical, provenance, or custody
defect was found.
