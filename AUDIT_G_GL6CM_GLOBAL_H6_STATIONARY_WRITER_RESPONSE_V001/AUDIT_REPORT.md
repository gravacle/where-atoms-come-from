# Independent hostile audit — GL6CM global `H6` stationary writer response

**Target:** `LANE_CROSS_RFT_GRA_GL6CM_GLOBAL_H6_STATIONARY_WRITER_RESPONSE_V001`  
**Verdict:** **PASS** on the declared finite-component, writer-only spectral
surface.

## 1. Independent audit surface

The science replay imports and executes no target program.  It independently
reconstructs:

1. the Perron--Frobenius ground-state and gap scope for a finite connected
   locked flip component;
2. the complete excited-state Gram factorization of the two-writer spectral
   response, including positivity, reciprocity, and the exact kernel;
3. the uniform cycle-rescaling null and the finite writer-map pullback
   `W^T K W`;
4. the isolated-`K2` zero and shared-star strict difference response in exact
   `Q(sqrt(2))` arithmetic; and
5. the physical `h/U_d` coefficient, source normalization, and dimensions.

The executable is `verify_gl6cm_independent.py`; its frozen output is
`INDEPENDENT_RESULT.json`.

## 2. Finite-component stationary state — PASS

Put

\[
 A=\sum_cT_c,\qquad H_0=-JA,\qquad J>0.
\]

On the declared finite connected flip component, `A` is real symmetric,
entrywise nonnegative, and irreducible.  Perron--Frobenius therefore gives a
simple largest eigenvalue `rho(A)` and a unique positive eigenvector.  Because
`H_0=-JA`, this is exactly the unique ground ray of `H_0`.  Finiteness and
simplicity make the separation from the next energy strictly positive.

The audit also checked every one of the 771 labeled connected simple graphs
on two through five vertices.  For each `n`-vertex graph it independently
verified the constructive irreducibility witness

\[
 (I+A)^{n-1}>0
\]

entrywise: a connecting path has length at most `n-1` and can be padded by
the diagonal `I` steps.  Multiple positive flip labels on an edge do not
weaken the general argument.

Writing the excited eigenstates as `|n>` and their gaps as
`Delta_n=E_n-E_0>0`, the reduced resolvent is

\[
 R=\sum_{n>0}{|n\rangle\langle n|\over\Delta_n},
\]

which is strictly positive on the orthogonal complement of the ground ray.
The stated `H_0/PF` scope is therefore mathematically complete for the finite
connected component; it does not imply a thermodynamic phase or bulk gap.

## 3. Spectral response, sign, and exact kernel — PASS

For real cycle amplitudes `y_c`, define

\[
 B_y=\sum_c y_cT_c,
 \qquad S_{nc}=\langle n|T_c|0\rangle .
\]

Direct insertion of the spectral resolution gives

\[
 K_{cc'}
 =2\lambda_T^2\operatorname{Re}
   \sum_{n>0}{S_{nc}^{*}S_{nc'}\over\Delta_n}.
\]

Equivalently, in the real case used by the exact replay,

\[
 \boxed{K=2\lambda_T^2S^T\operatorname{diag}(\Delta_n^{-1})S.}
\]

It follows without a material or mechanical postulate that

\[
 y^TKy
 =2\lambda_T^2\sum_{n>0}
   {\left|\sum_cy_cS_{nc}\right|^2\over\Delta_n}\ge0.
\]

This proves reciprocity and positive semidefiniteness.  Because every gap is
strictly positive, equality holds exactly when every excited-state amplitude
vanishes:

\[
 y^TKy=0
 \iff QB_y|0\rangle=0
 \iff B_y|0\rangle\in\operatorname{span}\{|0\rangle\}.
\]

The replay exercised this equivalence on 625 exact integer vectors using a
rank-three four-cycle transition matrix with three unequal rational gaps.  A
second lower-rank construction verified the important ceiling that additional
dark combinations can exist; the uniform null need not be the whole kernel.

## 4. Uniform rescaling and the writer pullback — PASS

The common mode follows from the parent Hamiltonian itself, not from a fitted
zero:

\[
 \sum_cT_c=-{1\over J}H_0.
\]

For every excited state,

\[
 \sum_cS_{nc}
 =-{1\over J}\langle n|H_0|0\rangle=0.
\]

Therefore `S 1=0` and

\[
 \boxed{K\mathbf 1=0.}
\]

For the microscopic pair-source writer map `y=Wj`, substitution into the
same quadratic form gives, identically,

\[
 (Wj)^TK(Wk)=j^T\boxed{(W^TKW)}k.
\]

The independent replay checked 2,187 exact bilinear comparisons for a
nontrivial four-by-six writer map.  Positivity survives the pullback, while
its kernel is exactly the preimage under `W` of the cycle-response kernel.
Consequently, a local rank-three writer does not by itself establish the
rank of the global response; GL6CM correctly retains that qualification.

## 5. Isolated ring and two-overlap witness — PASS

For an isolated ring, the component is `K2`, and its only toggle is
proportional to `H_0`.  Its matrix element between the ground and excited
states vanishes, hence its spectral curvature is exactly zero.

For the two-arm star, ordered as `(center,arm0,arm1)`, the exact normalized
eigenstates may be chosen as

\[
 |g\rangle={\sqrt2|B\rangle+|A\rangle+|C\rangle\over2},\quad
 |m\rangle={|A\rangle-|C\rangle\over\sqrt2},\quad
 |u\rangle={\sqrt2|B\rangle-|A\rangle-|C\rangle\over2},
\]

with energies `-sqrt(2)J,0,+sqrt(2)J`.  The two arm-toggle transition rows
are

\[
 (\langle m|T_0|g\rangle,\langle m|T_1|g\rangle)
 =\left({1\over2},-{1\over2}\right),\qquad
 (\langle u|T_0|g\rangle,\langle u|T_1|g\rangle)=(0,0).
\]

Thus the exact cycle-amplitude response is

\[
 K_{\rm star}={\lambda_T^2\sqrt2\over4J}
 \begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
\]

The common mode is null and the relative mode is strict.  If `w_0,w_1` are
the derivatives of the two physical arm amplitudes along a scalar source
path, then

\[
 \boxed{-E_0''(0)={\sqrt2\over4J}(w_0-w_1)^2.}
\]

The replay also obtained this expression by twice differentiating the exact
ground branch

\[
 E_0(s)=-\sqrt{(J-sw_0)^2+(J-sw_1)^2},
\]

rather than relying only on the residue calculation.

## 6. Physical coefficient and units — PASS

The locked-ring and source-writer scales are

\[
 J={63\over8}{h^6\over U_d^5},\qquad
 \lambda_T={105\over16}{h^6\over U_d^6}.
\]

At the literal Q4 direction inherited from the independently audited GL6CK
witness, `x=(2,0)`, and therefore

\[
 w_0={105\over8}{h^6\over U_d^6},\qquad w_1=0.
\]

Substitution gives

\[
 \boxed{-E_0''(0)={175\sqrt2\over32}{h^6\over U_d^7}>0.}
\]

The unit-`Theta` direction gives half this value.  The reciprocal of the
literal positive spectral response is

\[
 \boxed{({\cal K}_{\rm spec})^{-1}
 ={16\sqrt2\over175}{U_d^7\over h^6}.}
\]

The dimensions are consistent.  Both `h` and `U_d` have energy units, so
`J` is energy and `lambda_T` is dimensionless.  The physical pair source
`j` and its cycle contraction `x` have energy units, while the cycle kernel
has inverse-energy units.  Along `j=s jhat`, `s` has energy units,
`jhat` and `O=dH/ds` are dimensionless, and `-E''(s)` has inverse-energy
units.  Its inverse therefore has energy units.  No calibration coefficient
has been silently inserted.

## 7. Hostile boundary attacks

### Writer-only versus full source Hamiltonian

**Pass.**  The differentiated family is now explicitly
`H_tilde(s)=H_0+sO_writer`.  The result is the two-writer-vertex spectral
piece, not the second derivative of the full source-dependent Feshbach
Hamiltonian.  A possible diagonal order-six first-source vertex remains
unclassified; its spectral square and cross terms are not included.

### Missing contact

**Pass as a ceiling.**  Every source-second contact, including the known
lower-order contact, remains outside the theorem.  Such terms can alter the
complete physical Hessian.  The packet does not call the positive writer
square the complete Hessian.

### Stationarity and Perron--Frobenius scope

**Pass as typed.**  The response uses one ground state and one resolvent of
the same finite component for both source legs.  It does not select a
thermodynamic state, prove a uniform bulk gap, or establish real-time
causality.

### Record authentication

**Pass as a ceiling.**  The input is a pair-source probe.  The theorem does
not claim that this probe is itself an authenticated record, that it is
autonomously generated by record formation, or that every active ring has a
record lineage.

### Accumulation and collar summation

**Pass as typed.**  The exact two-overlap value is the first strict witness.
It is not added once per dense-parent edge.  Larger components can contain
cross terms through their common resolvent and require an owner-once linked
calculation.

### Legendre-transform sign

**Pass.**  The repaired packet now explicitly defines the expectation
coordinate by `phi=-dE/ds` and the convex transform by
`Gamma(phi)=E(s)+s phi`.  Because the ground energy is concave, this gives

\[
 \Gamma''=(-E'')^{-1}={\cal K}^{-1}>0.
\]

The sign is therefore unambiguous.  The positive response is `K=-E''`, and
its reciprocal on a strict subspace is the displayed positive inverse
response.

### Bulk and gravity promotion

**Pass.**  The packet leaves open the full contact-plus-spectral response,
record authentication, linked clusters, thermodynamic/refinement control,
real-time response, locality and causal continuation, a metric, Ricci or
Einstein form, gravity, and Newton's `G`.

## 8. Verdict

**PASS.**  GL6CM proves a genuine same-state stationary response of the
off-diagonal order-six writer on every finite connected locked component.
That response is a derived positive spectral Gram form with an exact kernel,
an unavoidable common-rescaling null, and an exact pair-source pullback.
The concrete two-overlap component supplies its first strictly positive
relative mode with coefficient `175 sqrt(2)/32 h^6/U_d^7`.  The result is a
valid finite microscopic accumulation mechanism.  It is deliberately not
the complete physical Hessian or a bulk/gravity theorem.
