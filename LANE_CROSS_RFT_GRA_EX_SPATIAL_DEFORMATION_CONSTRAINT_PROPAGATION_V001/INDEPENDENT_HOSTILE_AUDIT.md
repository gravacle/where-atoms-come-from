# Independent hostile audit: SDCP

**Lane:** `CROSS-RFT-GRA-EX-SPATIAL-DEFORMATION-CONSTRAINT-PROPAGATION-V001`

**Date:** 2026-08-27

**Audit mode:** independent no-edit review of the theorem and verifier; this
artifact was written by the hostile auditor after the theorem bytes stabilized

**Audited theorem SHA-256:**
`59487e2ce0585e291ecb032215ed3a9d23883e418df2994692448ced4cf5a1f2`

**Audited verifier SHA-256:**
`d926d2490afc14dc8615f1d1621a617a6da92fe2ae2a5da99616b0472f0c06f5`

## Verdict

**ACCEPT / CLEAN at the declared conditional claim class.**

The repaired theorem correctly proves that six arbitrary spatial metric
variations, together with the complete on-shell diffeomorphism Ward identity
and prospectively zero initial constraints, imply the full symmetric metric
residual vanishes on the admitted Gaussian-normal development.

This does not establish that EW's finite family has the required spacetime
functional range, that nature implements the soldering, that the nonmetric
shell and boundary terms vanish, or that the initial constraints, positive
induced coefficient, refinement, and record ancestry are physically realized.
It is a conditional residual-closure theorem, not gravity in nature.

## Hostile findings and repaired points

### 1. Variational index convention

The theorem defines the contravariant residual by

\[
 \delta\Gamma={1\over2}\int\sqrt{-g}\,
 {cal E}^{\mu\nu}\delta g_{\mu\nu}.
\]

It now consistently defines
`M_{Aij}=partial h_{ij}/partial J_A` and uses the invariant dual pairing

\[
 {\delta\Gamma\over\delta J_A}
 ={1\over2}\sqrt{-g}\,{cal E}^{ij}M_{Aij}.
\]

Full rank of the six covariant spatial variations therefore annihilates the
contravariant spatial residual. The earlier upper/lower-index ambiguity is
absent from the audited bytes.

### 2. Pointwise rank versus functional range

The fundamental-lemma step is valid only for six local fields
`J_A(tau,x)` with arbitrary compactly supported variations, or for a nonlocal
deformation map whose adjoint is injective in a declared function space. The
theorem now makes that requirement explicit.

It also correctly limits the EW composition. EW's inverse-function theorem
gives a six-by-six matrix isomorphism at one cell near the symmetric point; it
does not by itself provide six spacetime fields or a dense cellwise
refinement. SDCP separately requires that local-field or functional-range
lift, and requires every admitted variation to recheck all ten signed
intervals, the lapse/null Gram, clock bind, Lorentz signature,
nondegeneracy, and shared-edge gluing.

### 3. Complete off-shell Noether custody

The audited theorem displays the full identity before on-shell reduction:

\[
 \nabla_\mu{cal E}^{\mu}{}_{\nu}
 +\sum_B{cal N}^{B}_{\nu}[{cal E}_B]
 +{cal G}_\nu+{cal B}_\nu=0.
\]

The metric Ward identity follows only after every nonmetric Euler--Lagrange
residual is solved and the gauge, ghost, regulator, anomaly, and boundary
terms vanish, or after a separately proved homogeneous enlarged propagation
theorem. Merely listing an off-shell field is explicitly rejected. A direct
collective explicit force must vanish; a factorized force is removable only
when each residual factor has separately been solved to zero.

### 4. Exact Gaussian-normal propagation

For

\[
 ds^2=-d\tau^2+h_{ij}dx^idx^j,
 \qquad K^i{}_j={1\over2}h^{ik}\partial_\tau h_{kj},
\]

the auditor independently recomputed the divergence of a symmetric
contravariant tensor. Once `E^{ij}=0`, the spatial Ward components are exactly

\[
 \partial_\tau(\sqrt h\,{cal E}^{0j})
 =-2\sqrt h\,K^j{}_i{cal E}^{0i}.
\]

This is a homogeneous linear system along each normal curve. In the declared
uniqueness class, zero initial `E^{0i}` gives `E^{0i}=0`. The normal Ward
component then becomes

\[
 \partial_\tau(\sqrt h\,{cal E}^{00})=0,
\]

so zero initial `E^{00}` gives `E^{00}=0`. Together with the six spatial
equations, all ten symmetric components vanish. No linearization about flat
space is used.

The theorem now explicitly requires the residual and coefficients to lie in a
class for which this homogeneous system has a unique solution; weak or
higher-derivative versions require their own distributional uniqueness
theorem.

### 5. EW algebra and the physical ceiling

The companion EW audit independently accepted the exact pair-family algebra:
the unhalved derivative is

\[
 \left.\partial_{J_{ab}}{cal F}_\theta\right|_0
 =v_av_b^{\mathsf T}+v_bv_a^{\mathsf T},
\]

and the six tetrahedral dyads span `Sym^2(V)`. The SDCP verifier replays this
finite algebra, including determinant `-1/2` in its declared orthonormal
coordinate convention, pair-covariance rank six, spatial Fisher tangent rank
six, and exact uniform one-port marginals.

The theorem does not confuse that cellwise algebra with the missing physical
functional lift. It also limits the Einstein--Hilbert consequence to the
Fierz--Pauli/helicity-two kinetic operator about an admitted solution; an
actual protected retarded massless pole retains low-energy, remainder,
hyperbolicity, causal-boundary, and background conditions.

## Replay and disposition

Fresh replay produced:

```text
PASS SDCP exact checks
tetra_pair_deformation_det=-1/2
pair_covariance_rank=6; spatial_fisher_tangent_rank=6
global_flip_single_port_marginals=1/2
```

No theorem or verifier byte was edited by this audit. No dependency ledger,
verification transcript, or manifest existed in the draft lane at audit time;
those remain builder custody before any source freeze.

**Independent disposition:**

`CLEAN_CONDITIONAL_SIX_SPATIAL_VARIATIONS_PLUS_COMPLETE_ON_SHELL_WARD_IDENTITY_AND_PROSPECTIVELY_ZERO_INITIAL_CONSTRAINTS_CLOSE_THE_FULL_METRIC_RESIDUAL__EW_POINTWISE_PAIR_FISHER_RANK_ACCEPTED_BUT_LOCAL_FIELD_OR_FUNCTIONAL_RANGE_LIFT_REMAINS_AN_EXPLICIT_PREMISE__PHYSICAL_SOLDERING_NONMETRIC_SHELL_INITIAL_CONSTRAINT_INDUCED_COEFFICIENT_REFINEMENT_ANCESTRY_AND_GRAVITY_REMAIN_OPEN`
