# Independent hostile audit: S4 lineage-to-metric response-kernel reduction theorem

**Theorem ID:** `RGRL-LMRK-V001`

**Date:** 2026-08-27

**Audit mode:** independent no-edit algebraic, representation-theoretic,
response-theoretic, convention, scope, and dependency-custody audit

**Audited theorem SHA-256:**
`49e97e9cd3c9d8c75c65f3717156071bfcc0d88b3be3118aa442f74fb711f50d`

## Verdict

**ACCEPT / CLEAN.**

The theorem proves an exact finite-dimensional kinematic reduction and an
exact point-local/zero-spatial-momentum response count.  It does not promote
that count to a nonzero-momentum dispersive classification, infer dynamics
from Fisher positivity, or derive physical scale, curvature, Einstein--
Hilbert stiffness, or numerical `G` from scalar gamma or record count.

## 1. Projectors and explicit deformation map

The opposite-edge operator satisfies `O^2=I`, preserves the uniform edge
vector, and has three-dimensional `+1` and `-1` eigenspaces.  Consequently

\[
 P_A={\mathbb J_6\over6},\qquad
 P_T={I-O\over2},\qquad
 P_E={I+O\over2}-P_A
\]

are mutually orthogonal idempotents of ranks `1`, `3`, and `2`.  Their sum is
the identity.  This independently reproduces
`E_edges = A1 + E2 + T2` with dimensions `1+2+3`.

In the stated orthonormal tetrahedral coframe,

\[
 v_1={1\over2}(1,1,1),\quad v_2={1\over2}(1,-1,-1),\quad
 v_3={1\over2}(-1,1,-1),\quad v_4={1\over2}(-1,-1,1),
\]

so `v_a.v_b=delta_ab-1/4`.  Directly evaluating
`D(e_ab)=v_a odot v_b` gives equal diagonals and opposite off-diagonal signs
within each opposite-edge pair.  Summing and differencing those matrices gives

\[
 D(P_Ax)=-{s\over6}I_3,
\]

\[
 D(P_Ex)=\operatorname{diag}(p'_1,p'_2,p'_3),
\]

and off-diagonal entries

\[
 (D P_Tx)_{yz}=-{d_1\over2},\quad
 (D P_Tx)_{xz}=-{d_2\over2},\quad
 (D P_Tx)_{xy}=-{d_3\over2}.
\]

Thus the displayed map, its signs, and its normalization are correct.  The
uniform mode maps to trace, while the remaining five modes map to the complete
symmetric-trace-free target.  The source and target representations match and
`D` is an equivariant isomorphism.

## 2. Pullback Gram and conditioning

Using

\[
 \langle u\odot v,x\odot y\rangle_F
 =2[(u\cdot x)(v\cdot y)+(u\cdot y)(v\cdot x)]
\]

reproduces Gram entries `5/4` for the same edge, `-1/4` for adjacent
edges, and `+1/4` for opposite edges.  Hence

\[
 D^*D={5\over4}I-{1\over4}A+{1\over4}O.
\]

On `(A1,E2,T2)`, the line-graph eigenvalues are `(4,-2,0)` and the
opposite-edge eigenvalues are `(+1,+1,-1)`.  Substitution yields exactly

\[
 D^*D={1\over2}P_A+2P_E+P_T,
\]

with singular values `1/sqrt(2)`, `sqrt(2)`, and `1`.  The theorem therefore
correctly distinguishes the raw edge norm from the output Frobenius norm.

## 3. Schur count and rotational compatibility

The final theorem now explicitly premises that the six-component lineage
source tangent carries the same natural `S4` edge action as the pair tangent.
At an `S4`-fixed background, an equivariant derivative-zero or `k=0`
lineage-to-pair convolution kernel acts by one generally complex frequency-
space scalar on each of the three inequivalent, multiplicity-one real
irreducibles.  Therefore

\[
 H^R(\omega,0)=h_A^R P_A+h_E^R P_E+h_T^R P_T
\]

is exact.  The equivalent same/adjacent/opposite matrix census gives
`h_A=a+4b+c`, `h_E=a-2b+c`, and `h_T=a-c`.

The theorem does not infer `O(3)` from `S4`.  It requires the source tangent to
carry the pulled-back action `D^{-1}[R(Dx)R^T]` and the response to intertwine
that action with the output tensor action.  Only under that compatible-source
premise does `STF(V)` become one `l=2` irreducible and force `h_E=h_T`, leaving
trace plus shear.  An output-only rotational action is correctly rejected.

The three-factor statement is restricted throughout to the spatially
point-local/derivative-zero term or the kernel value at `k=0`.  The theorem
explicitly notes that a transforming nonzero wavevector permits longitudinal,
transverse, and other tensor structures and withholds a classification of the
full dispersive kernel.  No prose or disposition silently globalizes the
three-factor count.

## 4. Gamma coefficient and exact dynamical nonidentification

At `J=0`, the six pair Walsh characters have covariance `I_6`.  For squared
complete-query fidelity, the standard local expansion is

\[
 \gamma_Q(0,\delta J)
 =1-{1\over4}\delta J^T F_J(0)\delta J+O(\|\delta J\|^3),
\]

so the coefficients `1/4` in both `-log gamma_Q` and `1-gamma_Q` are correct.
At quadratic order, a fixed scalar decrement leaves the five-sphere of pair
directions; unit vectors in the three sectors have the same decrement but
distinct trace/shear outputs and Gram weights.

The stronger nonidentification is also exact.  Given any nonzero admitted
retarded passive scalar response `f^R`, the two laws

\[
 H_1^R=f^R I,
 \qquad
 H_2^R=f^R(2P_A+3P_E+4P_T)
\]

are causal and `S4` equivariant, remain passive under the declared work
pairing, and are sector-nondegenerate wherever `f^R` is nonzero.  They can be
joined to the identical EW state/query family, complete-query gamma, and
carrier count but yield different response kernels.  The `O(3)`-compatible
two-sector variant preserves the same counterexample.  Because the EW family
does not specify a Hamiltonian, physical source, response observable, or KMS
same-observable axiom, no fluctuation--dissipation premise has been silently
used.  For actual independent equal carriers, Fisher additivity and
`gamma_N=gamma_1^N` fix accumulation only, not response direction or
normalization.

## 5. Causality, passivity, Ward identity, and scale ceilings

With the declared transform convention, retardedness gives upper-half-plane
analyticity subject to the stated stability/temperedness and subtraction
conditions.  Before a physical cone is earned, the support statement is only
mission-time retardedness.

For harmonic input proportional to `exp(-i omega t)`, the work convention
`dot W=<delta L,partial_t delta J>_GD` gives

\[
 \overline{\dot W}={1\over2}\sum_r
 \omega g_r\operatorname{Im}h_r^R\|P_r\delta L\|^2,
 \qquad (g_A,g_E,g_T)=({1\over2},2,1).
\]

The stated positive-frequency sign is therefore correct.  It is explicitly
conditional on independently drivable sectors and this particular conjugate
work pairing; Fisher positivity itself supplies neither passivity nor a
response magnitude.

The Ward discussion accurately imports the EX endpoint ceiling: only a fully
soldered continuum effective theory with the complete on-shell Noether/Ward
system, complete spatial variations, prospective zero initial constraints,
and the associated EX propagation premises can close the remaining normal
equations.  The identity does not identify `D` with a physical metric, choose
the three form factors, supply isotropy, select Einstein--Hilbert dynamics, or
fix its coefficient.

Finally, the convention crosswalk is exact:

\[
 s=\ell_F^2F,
 \qquad q={\ell_B^2\over4}Q,
 \qquad \ell_B=2\ell_F
\]

when `F=Q` denotes the same Fisher/QFI form.  Holding `ell_F` independently
fixed in the variation is necessary and correctly stated.  Its value, the
source normalization, proper-time/cone calibration, dispersive kernel, and
numerical `G` all remain open.

## 6. Dependency custody

All six live dependency hashes were recomputed and match the theorem ledger:

- EW: `495e4e99171f4e3e5809f24e5a9a5b68116e996f4f29115bd22f780127d4714e`;
- EX: `59487e2ce0585e291ecb032215ed3a9d23883e418df2994692448ced4cf5a1f2`;
- induced-EH: `e73c077d1402ebab9a6061b060b9b24e9ac323b25809cccc706f4365da7a5e2f`;
- GSGB: `8721183a72f4b864d06f79f6e68405a5393e37d1a6f4d24d5f4c3c2b79a81075`;
- length crosswalk: `e3037d5fcc0b449b8c46414de94365075e372f8da8ca1b6a87df85c7ef85359b`;
- post-adoption structural theorem:
  `733b18ecaa29c7acd755db6947b790a9ae37240a3c74d199752d5e278280783d`.

No theorem or dependency byte was edited by this audit.

## Final disposition

`CLEAN_EXACT_S4_KINEMATIC_AND_K0_RESPONSE_REDUCTION__PROJECTORS_D_MAP_GRAM_AND_SCHUR_COUNT_REDERIVED__O3_REDUCTION_REQUIRES_COMPATIBLE_SOURCE_ACTION__SCALAR_GAMMA_AND_COUNT_EXACTLY_FAIL_TO_IDENTIFY_DYNAMICS__CAUSAL_PASSIVE_WARD_AND_LENGTH_CEILINGS_PRESERVED__SPATIAL_DISPERSION_PHYSICAL_SOLDERING_ABSOLUTE_SCALE_AND_NUMERICAL_G_REMAIN_OPEN`
