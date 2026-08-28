# Independent hostile audit — FW FV-WITNESS projected response

**Audit date:** 2026-08-28  
**Audited lane:** `LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001`  
**Disposition:** `PASS_AFTER_PRE_FREEZE_WITNESS_SCOPE_REPAIR`  
**Independent replay:** `144/144` checks passed  
**Builder replay:** `114/114` checks passed

## 1. Bottom line

The frozen FW theorem survives independent hostile audit in its repaired scope.  The result is exactly a response theorem for

\[
Q_{\rm FV-WITNESS}=Q_{\rm pair}^{(0)}+Q_{{\rm ring,irr}}^{(6)},
\]

not for FV's complete fixed-order source.  Conditional on the repaired FV `S10 / FV-PURE` premise and the FO finite-sector construction, the independently recomputed hierarchy is

\[
\boxed{
6_{\rm FV\ family}
\longrightarrow 5_{\rm component/mod\ I}
\longrightarrow 3_{\operatorname{ad}_H}
\longrightarrow 2_{\chi^R_{|0\rangle}}
=2_{M_1}.}
\]

The exact ground response contains two orthogonal rank-one residues at

\[
\Delta_1=2+2\sqrt2,
\qquad
\Delta_2=4+2\sqrt2.
\]

No remaining mathematical or documentary defect was found in the frozen claim.

## 2. Material pre-freeze defect and repair boundary

During the audit window, the builder identified a material scope defect in its candidate packet: FV13 retains uncomputed generated diagonal derivatives and folds `Q_diag^(2,4,6)`, so the two explicitly calculated rank-witness pieces could not be called the complete FV response.  The builder repaired the theorem before freezing it by:

1. defining the load-bearing source as `FV-WITNESS`;
2. stating that all generated `Q_diag^(2,4,6)` terms and their differentiated folds are omitted;
3. withdrawing any claim that rank two upper-bounds the complete fixed-order or CTP response; and
4. making generated-diagonal completion the immediate successor calculation.

This audit pins only the repaired bytes.  It does not retroactively validate the broader pre-freeze wording.

## 3. Independent reconstruction

The audit executable imports or executes neither the FW verifier nor the FO verifier.  It independently:

1. rebuilt the `Z_30` diamond quotient using frozenset occupation states;
2. re-enumerated all 120 elementary hexagons and their four missing-label orientations;
3. regenerated the winding-seed ring component with 180 exact ice states and 420 undirected transitions;
4. reconstructed the six free translation orbits and exact integer `6 x 6` zero-momentum Hamiltonian;
5. recomputed the FV family witness determinant from the two direct `E` rows and four ring `A1+T2` rows;
6. rebuilt the direct pair and irreducible ring sources from normalized tetrahedral/root dyads, including the exact FV11-to-`J6` factor `-8/63`;
7. evaluated the full 180-state operator and commutator Grams with integer arithmetic; and
8. evaluated the zero-momentum response in the exact number field `Q(sqrt(2),sqrt(3))`, using rational polynomial energy projectors rather than degenerate numerical eigenvectors.

The independent FV family witness has

\[
\det W=-\frac{4678629417}{256}\ne0,
\]

so its off-shell operator rank is six.  Restriction to the selected homogeneous component gives the exact full-component identities

\[
Q_A=\frac{60\rho I+11H}{\sqrt3},
\qquad
Q_{E_c}=16\sqrt6\rho I,
\qquad
[H,Q_{T_{yz}}]=0.
\]

The direct component contributes one nonidentity `E` direction modulo identity; the ring component contributes four `A1+T2` directions.  Their tensor sectors and diagonal/off-diagonal support are separate, proving rank five for every nonzero `rho`, not merely at the two numerical replay values.

## 4. Commutator versus response rank

The independent Hilbert-Schmidt commutator Gram splits exactly into a direct rank-one `E` block and two orthogonal ring directions:

\[
C_E=\rho^2
\begin{pmatrix}
960&-960\sqrt3\\
-960\sqrt3&2880
\end{pmatrix},
\qquad
C_{T_{xy}T_{xy}}=C_{T_{xz}T_{xz}}=25920.
\]

The direct/ring cross Gram vanishes.  Hence `rank(ad_H)=3` for `rho != 0`.  This is an operator rank; it is not the Kubo or ground-state spectral rank.

The exact zero-momentum characteristic polynomial is

\[
x^2(x-2)^2(x^2+4x-4),
\]

with normalized ground vector and energy

\[
g_0=\left(\frac12,\frac1{2\sqrt2},\frac1{2\sqrt2},
\frac1{2\sqrt2},\frac1{2\sqrt2},\frac12\right),
\qquad E_0=-2-2\sqrt2.
\]

The audit rebuilt the exact polynomial projectors onto the degenerate energies zero and two.  Applied to all six source vectors, they give

\[
\begin{aligned}
r_1&=\left(0,\frac\rho{\sqrt2},-\rho\sqrt{\frac32},
-\frac3{\sqrt2},-\frac3{\sqrt2},0\right),\\
r_2&=\left(0,0,0,\frac3{\sqrt2},-\frac3{\sqrt2},0\right),
\end{aligned}
\]

and no centered remainder outside the two projectors.  Three exact scale points certify the degree-two residue polynomial for arbitrary `rho`.  Thus both residues are rank one, their directions are orthogonal, and the combined retarded image has rank two.  A direct exact commutator calculation independently gives

\[
M_0=0,
\qquad
M_1=-2(\Delta_1R_1+\Delta_2R_2),
\qquad
\operatorname{rank}M_1=2.
\]

## 5. Hostile boundary tests

- **Component restriction:** passed.  The selected 180-state winding component excludes its complement partner and is not treated as the complete periodic ice Hilbert space.
- **Identity removal:** passed.  Identity is included as an explicit seventh Gram direction before quotienting; the component source has rank five modulo identity.
- **Normalization and sign:** passed.  FV11 multiplied by `-8/63` gives `-31 I/6+9D_d/2` in `J6` units for all four missing labels; the direct source retains `lambda=-1/2` and the independent ratio `rho=U_d/J6`.
- **Translation and state restriction:** passed.  All twelve direct/ring coordinate pieces commute with cyclic translation, the unique ground lies at zero momentum, and the response reduces exactly to the six-orbit block.
- **Degeneracy artifacts:** excluded.  Residues are formed with exact polynomial projectors onto full two-dimensional eigenspaces; no eigenvector in a degenerate subspace is individually scored.
- **Commutator/Kubo conflation:** excluded.  Operator rank five, commutator rank three, and ground retarded/first-moment rank two are separately recomputed.
- **Conserved `T2`:** verified as `[H,Q_Tyz]=0` on the full selected component, but retained only as a component-specific conservation law—not a Ward identity.
- **Ground-dark direction:** verified.  After expectation subtraction, `Q_dark|0>=0`, while its exact commutator is nonzero.  It is state-dark, not conserved or gauge-generated.
- **Finite `k=0` overgeneralization:** excluded.  The theorem withholds nonzero momentum, a local source, thermodynamic scaling and complete-sector conclusions.
- **Complete-source overclaim:** excluded in the frozen packet.  Generated diagonal derivatives/folds remain uncomputed and can change all complete-response ranks.
- **Ward/helicity/gravity promotion:** excluded.  No finite conserved direction or discrete pole is called a Ward identity, helicity-two particle, tensor pole, RGRL-B, gravity or `G` result.

## 6. Verdict and exact ceiling

**PASS.**  The repaired FW packet is an exact finite, homogeneous, zero-momentum response theorem for `FV-WITNESS`.  Its rank hierarchy and two-pole residue certificate are correct.

The audit does not validate a complete fixed-order FV response, a CTP/Ward packet, a nonzero-momentum mode, a thermodynamic pole, helicity two, gravity emergence or a value of `G`.  The immediate calculation remains the generated `Q_diag^(2,4,6)` and differentiated-fold completion on the same component.
