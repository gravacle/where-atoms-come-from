# Distinct post-freeze hostile audit — GL6AK A3 quasi-local bulk dynamics V001

**Target:** `LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/`  
**Frozen theorem SHA-256:** `083d5fbb8a48e27e365167075da132ffa23e395587a4c0e40cc572d8b761ad30`  
**Frozen MANIFEST-file SHA-256:** `d38f89c618ea6f77c7b399b005ad0f0abe04d3865e06921f8c765feb44f40620`  
**Frozen SEAL-file SHA-256:** `322bf51a00f8fea3f36a09656dda4ebf89ba56b9a88d60b50e9cc7ab33223987`  
**Disposition:** `PASS__FINITE_AUTHENTICATED_ANCESTRY_AND_A3_DEGREE6_EXACT__BOUNDARY_FACTOR3_AND_CAUCHY_TAIL_EXACT__JOINT_INVARIANT_STATE_EXISTENCE_SOUND__POSITIVE_A1_E_T2_LIOUVILLIAN_MEASURES_SOUND__NO_GLOBAL_RECORD_STATE_SELECTION_POLE_PHYSICAL_MOMENTUM_RICCI_GRAVITY_OR_G`

## Custody

The terminal theorem, manifest, and seal match the pins above.  The author
physics replay passes `6304/6304`; the frozen packet and all dependency and
pre-audit pins pass `104/104`.  The final manifest pins twelve author files,
and the seal pins that manifest.

An initial audit read occurred while the author was still completing seal
custody and correctly rejected that transient snapshot.  No audit was sealed
against it.  After the author declared the final bytes stable, this audit
started again from a fresh read and held the three terminal hashes fixed
through sealing.

The frozen `PRESCREEN_AUDIT.sha256` also resolves all five historical
pre-freeze audit bytes exactly, including the independent `79644/79644`
result.  Those files remain unchanged in this directory; the present
`POSTFREEZE_*` files supply the distinct terminal review.

## Independent physics reconstruction

The distinct spot replay imports neither author verifier nor pre-freeze
replay and passes `33398/33398`.  It reconstructs the two unordered edge
families directly on `A3 x {1,2,3,4}`.  Every active-link site has three
same-parent and three same-child partners, while every cell touches six
internal and twelve shared-child terms.  A non-symmetric finite collar maps
into one strict-interior finite FPSS slab and preserves literal shared-child
equality.  This earns finite authenticated ancestry, not one infinite record.

Assigning each pair term to a minimum-radius endpoint gives no more than
eighteen owned terms per cell.  Combining that census with the inherited
two-endpoint commutator bound gives

\[
 {72J\over\hbar}\int_0^{|t|}T_d(\lambda u)du,
 \qquad \lambda={24J\over\hbar},
\]

and therefore exactly the coefficient `3` and shifted tail

\[
 3\|A\||X|\sum_{r\ge R}(2r+1)^3
 T_{r-r_X+1}(\lambda|t|).
\]

The factorial tail beats the cubic shell uniformly on compact time
intervals.  This supports norm-Cauchy convergence and independence for the
stated locally complete open restrictions of the inherited interaction.  It
does not cover arbitrary added boundary laws.

## State and spectral logic

Weak-* compactness plus continuous time averaging proves existence of a
stationary cluster state.  The displayed `Z3` cubes are Følner; translation
averaging preserves stationarity because time and translation actions
commute.  Finite `S4` averaging then preserves both properties.  This proves
existence only: no state is selected uniquely or operationally prepared, and
no ground, KMS, clustering, or finite-volume state convergence is obtained.

In the invariant state's GNS representation, strong continuity implements
the dynamics by `U(t)=exp(itL)`.  The spectral coordinate is consequently
angular frequency `nu`, with energy `E=hbar nu`.  The matrix measure

\[
 \mu_{AB}(B)=\langle\psi_A,P_L(B)\psi_B\rangle
\]

is positive by the projection norm identity.  The correlation and retarded
sign in AK31--AK32 are consistent with this convention.  Exact reconstruction
of the six-pair representation gives mutually orthogonal projectors of ranks
`1,2,3`, so `S4` invariance yields positive scalar measures in
`A1`, `E`, and `T2` without inserting a pole or gravitational operator.

## Strict ceiling

GL6AK closes a mathematical gate: the selected homogeneous finite parent has
a well-defined stationary infinite-volume response object.  It does not
authenticate one infinite record, autonomously select the all-formed member
or stationary state, prove a gap or mode, identify a displacement character
with physical momentum, supply a physical cone or continuum, or derive
Ricci/Einstein response, gravity, or Newton's constant.

**Post-freeze hostile verdict: PASS.**
