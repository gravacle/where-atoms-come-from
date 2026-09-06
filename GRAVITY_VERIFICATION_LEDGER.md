# Gravity verification ledger

**Status:** active, fail-closed verification ledger  
**Authority:** two-tier UV/IR verification structure, adopted 2026-09-05

## Tier separation

The microscopic requirement is exact discrete ledger conservation, not a
one-cell continuum Ward identity.  For a declared finite region \(R\), the
target law is

\[
\Delta Q_R+\sum_{e\in\partial R}J_e-W_R=0.
\]

`Q` is authenticated retained lineage, `J` is an owned transported graph
current, and `W` is an explicitly owned write/read/source contribution.
On a periodic closed member, paired transported seam currents must telescope
exactly.  Continuum-style Ward decay is a separate infrared question:

\[
\widehat\delta(L)=
 {\|P_{\rm bulk}{\cal R}_L\|\over Z_L|k_L|^3}
 \longrightarrow 0.
\]

The second expression is meaningful only after the same-parent family,
physical source/read map, bulk projector, response scale \(Z_L\), and error
control have been earned.

## L=4 baseline — PASS at finite selected-response scope

`RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003/` and its independent result
audit `AUDIT_G_GL6FJ_M001_RESULT_V001/` pass.  The finite unaliased selected
response has final eigenvalue interval

\[
[3.3147386524,\;177.5395794604],
\]

with independent six-owner reconstruction residual
\(6.73\times10^{-16}\).  This is not a physical 1PI kernel, a Ward null, or
gravity.

## Raw 64-character source classification — PASS, PARTIAL

`AUDIT_G_GL6CY_RAW_SOURCE_COMPOSITION_V001/INDEPENDENT_RESULT.json` records
an independent 983-check reconstruction of the declared raw source
composition

\[
{\cal R}^{\rm raw}_{sr}=H_{,AB}\eta_{A,s}\eta_{B,r}
                         +H_{,A}\eta_{A,sr}.
\]

The nonlinear \(\eta_{sr}\) is nonzero and was explicitly included.  On the
declared path its contraction with \(H_{,A}\) vanishes by the inherited
locked/alternating structure; it was not set to zero.  The source-owned
periodic zero-mode is nevertheless nonzero:

\[
{\cal R}^{h4}_{\rm raw}(0)={32128\over27}S,
\qquad
{\cal R}^{h6,\rm diag+cyc}_{\rm raw}(0)
 ={7212448\over6075}S,
\]

where \(S_{xy,z}=-1\), \(S_{xz,y}=+1\), and all other displayed entries
vanish.

| Component | Disposition | Reason |
| --- | --- | --- |
| Source | `OWNED_PARTIAL` | The complete raw H6 source owner is reconstructed. |
| Boundary/seam | `ALGEBRAIC_TELESCOPING_PROVED__PHYSICAL_FLUX_UNPROVED` | GL6DA proves oriented cochain cancellation and retained periods, not a physical stress/current flux. |
| Lattice | `NOT_PURE_LATTICE` | The raw periodic zero mode above is nonzero. |
| Bulk | `PARTIAL_COHERENT_ZERO_MODE` | A raw source-only periodic bulk component is present; it is not a physical anomaly because the total physical owner composition is absent. |

Thus the raw source residual is neither a demonstrated probability/energy
leak nor an emergent Ward remainder.  It is an exact partial owner balance
that must be combined with state, measure, retained-field, boundary,
matching, constraint, and lawful quotient owners.

## Gate dispositions

### Gate A–UV — OPEN

The two-tier structure is adopted, but Gate A–UV is **not yet certified
closed**.  The following remain required:

1. one complete finite update/CTP ledger that identifies \(Q,J,W\) on the
   same physical history;
2. exact transported pairwise seam cancellation for that physical \(J\),
   not only the existing algebraic overlap cochain;
3. a complete owner-once classification of the residual; and
4. no unowned term in the discrete balance law.

This is a discrete requirement.  It does not require
\(\delta_{\rm dual}(4)=0\), a smooth manifold, or a continuum gauge quotient.

### Gate B–IR — DEFINED, NOT YET LAUNCH-AUTHORIZED

The Gate-B observable and acceptance criterion are now formally fixed by
\(\widehat\delta(L)\) above.  An \(L=8\) production run remains unauthorized
until Gate A–UV has the required physical ledger; otherwise an apparent
defect decay would have no physical normalization or ownership meaning.

## Claim boundary

No statement in this ledger derives a physical Ward identity,
Einstein/Fierz--Pauli response, \(1/r\) exchange, gravity, \(C_R\), or \(G\).
