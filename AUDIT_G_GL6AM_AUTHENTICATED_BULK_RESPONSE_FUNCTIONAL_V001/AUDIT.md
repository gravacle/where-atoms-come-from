# Distinct hostile audit — GL6AM authenticated bulk response V001

**Target:** `LANE_CROSS_RFT_GRA_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/`  
**Frozen theorem SHA-256:** `8407cee5196bfa4240f02159a5f59f941903dcf7a10e2baa18cf52a01ac8f743`  
**Frozen author-manifest SHA-256:** `8b4ac6f6ceda2acb117480201ee96ce22be97fe0a99d8c097d8267100efa8c44`  
**Disposition:** `PASS__FINITE_AUTHENTICATED_PULSE_AND_DEFECT_WINDOWS_HAVE_BOUNDARY_INDEPENDENT_BULK_FUNCTIONALS__RETARDED_THETA_AND_FACTORIAL_TAIL_EXACT_NOT_STRICT_CONE__FINITE_WINDOW_CORRELATION_MEASURE_POSITIVE__A1_E_T2_ONLY_S4_CLOSED__K_DEFECT_COCYCLE_AND_STATE_INDEPENDENT_TAIL_SOUND__NO_STATE_SELECTION_BULK_COEFFICIENT_MOMENTUM_CONE_GRAVITY_OR_G`

## 1. Custody and independence

The author froze all eleven GL6AM files before this replay.  Their exact
hashes are pinned in `AUDITED_TARGETS.sha256`; the author manifest covers its
ten non-manifest files and passes.  GL6AM did not self-seal.

The audit replay imports no author module.  It independently performs 1476
checks.  The author structural replay separately passes `3999/3999` and its
packet/custody check passes `131/131`.  The pinned GL6AK post-freeze audit
passes `87/87`, and its independent reconstruction passes `33398/33398`.

All forty-two direct dependency hashes resolve, including the GL6AK theorem,
lane manifest/seal, original audit manifest/seal, and distinct post-freeze
audit and verification.  The GL6AK audit manifest/seal in turn covers its
post-freeze replay code.  No mutable GL6AL result is in the premise chain.

## 2. Finite-pulse norm limit

For a fixed finite pulse word, write

\[
 W_R=U_{n,R}\cdots U_{1,R},\qquad U_{k,R}=\tau_{s_k}^{(R)}(V_k).
\]

All factors are unitary.  The exact product identity

\[
 \prod_kU_{k,S}-\prod_kU_{k,R}
 =\sum_jU_{n,S}\cdots U_{j+1,S}
 (U_{j,S}-U_{j,R})U_{j-1,R}\cdots U_{1,R}
\]

gives

\[
 \|W_S-W_R\|\le\sum_j\mathcal E_R(V_j,s_j).
\]

One more exact conjugation telescoping gives

\[
 \|W_S^*\tau_t^{(S)}(B)W_S-W_R^*\tau_t^{(R)}(B)W_R\|
 \le \mathcal E_R(B,t)
 +2\|B\|\sum_j\mathcal E_R(V_j,s_j).
\]

The replay independently reconstructs the authenticated `A3 x 4` degree as
three same-parent plus three same-child neighbors.  With `J=2|U_d|` and
`lambda_F3=24J/hbar`, the sealed GL6AK boundary coefficient reduces exactly
as

\[
 {72J/\hbar\over\lambda_{\rm F3}}=3.
\]

The shifted factorial tail beats the cubic shell uniformly on compact time
windows.  Therefore every fixed finite pulse/read word has the stated norm
limit for the allowed locally complete open exhaustions.  The argument does
not use, select, or require a finite-volume state limit.

## 3. Retarded functional and support ceiling

For `V(j)=exp(i j M_alpha/2)`, direct differentiation gives

\[
 {d\over dj}\left[V(j)^*B V(j)\right]_{j=0}
 ={i\over2}[B,M_\alpha].
\]

Thus the GL6AM sign and normalization

\[
 \mathcal G^R_{\beta\alpha}(t)
 ={iE_\star^2\over2\hbar}\Theta(t)
 \omega([\tau_t(M_\beta),M_\alpha])
\]

are correct.  A commutator difference supplies the factor two converting the
AK observable error to `E_star^2 mathcal E_R/hbar`, and the arbitrary-support
AK/AI commutator estimate supplies exactly the double support sum in AM12.

The `Theta` factor makes the operational retarded response exactly zero at
negative response time.  For positive time the factorial tail is generally
nonzero at every finite relational distance.  This is a retarded, quasi-local
influence envelope, not strict microcausality and not a physical cone.

## 4. Positivity and exact `S4` scope

For any chosen joint invariant GL6AK state and every finite list of centered
pair reads,

\[
 \mu_{\alpha\beta}(\mathcal B)
 =\langle\psi_\alpha,P_L(\mathcal B)\psi_\beta\rangle
\]

is a positive matrix measure by the projection-Gram identity.  The distinct
replay checks this with complex coefficients, not only real test vectors.
It also reconstructs the retarded spectral difference and confirms that
positivity belongs to the correlation measure; a merely stationary state
does not make the commutator/dissipative difference positive.

The six-pair permutation representation is reconstructed under all twenty-
four `S4` permutations.  Its mutually orthogonal projectors have ranks
`1,2,3` and resolve the identity, so an `S4`-closed one-cell or common finite
envelope has positive scalar `A1`, `E`, and `T2` measures.  A single
off-origin localized envelope has a nontrivial `S4` orbit and is not closed;
it receives covariance only, not the scalar sector theorem.  A character
label is consequently not promoted to physical momentum.

## 5. Finite `K`-word defect and nonequilibrium ceiling

On the selected all-formed branch, the exact finite perturbation

\[
 V_\kappa=h\sum_{p\in D}(1-\kappa_p)X_p
\]

changes `-hX_p` to `-h kappa_p X_p` and nothing else.  The standard bounded-
perturbation cocycle with

\[
 i\hbar\dot W_R=\tau_t^{(R)}(V_\kappa)W_R,
 \qquad \gamma_t^{R,\kappa}(A)=W_R^*\tau_t^{(R)}(A)W_R
\]

has the sign needed to generate `H+V_kappa`.  Unitary Duhamel comparison and
the same AK norm tail give AM23--AM24.

For two finite words, telescope one changed link at a time.  The route-
uniform AI commutator bound is an operator-norm statement, so it yields the
state-independent factor

\[
 {2h\|B_Y\|\over\hbar}
 \sum_p|\kappa_p-\kappa'_p|\sum_{q\in Y}
 \int_0^{|t|}T_{d_L(p,q)}(\lambda_{\rm F3}u)\,du.
\]

The replay checks the exact shifted-tail integral and the `U_d=0`
distinction: positive-distance influence vanishes while overlapping support
is not forced to vanish.

The original homogeneous invariant state is generally not invariant under a
defect generator.  GL6AM correctly withholds stationary positivity,
homogeneous sectors, and AF/AG numerical coefficients from the defect
contrast.  No hidden state identification occurs.

## 6. Hostile promotion attacks

The audit explicitly attacked and rejected the following possible hidden
promotions:

1. **Infinite mission:** every source, read, and collar remains finite.
2. **State selection:** the stationary response is conditional on any chosen
   member of a nonempty invariant family; no vacuum, ground, Gibbs, or KMS
   state is selected.
3. **Bulk coefficient:** AF/AG prepared-blank finite coefficients authenticate
   ancestry only and are never inserted into the stationary functional.
4. **Translation/sector overreach:** arbitrary envelopes and defect words
   receive covariance, not automatic scalar `A1/E/T2` blocks.
5. **Spectral overreach:** positive correlation measure is not called positive
   absorption or a pole.
6. **Causal overreach:** `Theta` time support plus factorial relational decay
   is not called a strict spatial or Lorentz cone.
7. **Gravity overreach:** no physical momentum calibration, complete stress,
   Ricci/Einstein comparison, gravity identification, or `G` calculation is
   present.

## 7. Verdict

The theorem is sound at the pinned snapshot.  It closes an operational and
thermodynamic-completion gate: every fixed finite authenticated pulse/read or
finite retained-`K` defect mission has one boundary-independent bulk
functional with the stated exact tails.  The next physics gate remains an
independently justified bulk-state preparation/selection and a controlled
growing-window/refinement analysis.  GL6AM itself proves no infrared pole,
physical cone, gravity, or `G`.

**Hostile verdict: PASS.**

