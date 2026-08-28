# Self-audit: inherited F3/q4 even-channel kernel boundary

**Lane:** `GRA-FM-F3-Q4-IEKB-V001`

**Audit date:** 2026-08-27

**Disposition:**
`PASS_EXACT_FIXED_SUPPORT_ORDER8_CLASSIFICATION_AND_SCOPED_POLE_BOUNDARY__INDEPENDENT_AUDIT_SEALED_SEPARATELY`

## 1. Scope audited

The audit checks only the unchanged `d_*=2`, `E_R=0` single-link-flip F3
parent on supplied finite simple coordination-four bipartite support of girth
at least six. The exact coefficient applies to every alternating simple
octagon. The complete small-torus equality remains subject to the same support
and cycle-set typing already imposed by `CW`; no support is derived here.

The continuum statement is a strict single-insertion spectral boundary. It
does not claim that the TT matching coefficient has been evaluated, that a
fully dressed massless 1PI kernel is analytic, or that a nonperturbatively
solved finite-order Hamiltonian cannot bind or shift a pole.

## 2. Hostile algebra checks

### A. Missing order-eight endpoint

**Attack.** Two cycles, a four-cycle, or another projected endpoint could have
been omitted.

**Result.** Rejected. The endpoint odd-multiplicity set is an even-degree
subgraph. With eight letters and girth six, a nonempty such set has size six
or eight. Size six is one hexagon plus two repeated occurrences; size eight is
one octagon. Two cycles require at least twelve odd links.

### B. Incomplete octagon coefficient

**Attack.** A proper prefix could return to ice, or a Feshbach fold could
alter the coefficient.

**Result.** Rejected. Every proper subset of the edges of a chordless simple
octagon is a forest and violates degree two. All 40320 paths are irreducible.
The verifier obtains thirteen denominator classes and independently reproduces

\[
 J_8={429h^8\over16U_d^7}
\]

by subset recursion. No lower off-diagonal kernel has an octagon endpoint, so
no fold contributes at its first appearance.

### C. Hidden order-eight potential

**Attack.** Four-edge diagonal histories could detect local ice structure and
generate a flippability potential.

**Result.** Rejected under the declared regular bipartite girth-six premises.
Every irreducible diagonal eight-word touches at most four links. Their
colored incidence subgraphs are forests, whose color-preserving counts are
fixed by the per-vertex occupied/empty degrees. Folds are scalar or dress the
existing hexagon transition. The executable hostile `PG(2,3)` coordination-
four example compares two degree-two ice states separated by a hexagon flip,
checks identical colored two/three/four-edge censuses, and obtains the same
exact irreducible diagonal checksum
`2526594309109/13608000`.

### D. Uniform hexagon coefficient silently assumed

**Attack.** Collapse every order-eight hexagon correction to one bulk number.

**Result.** Rejected. The theorem retains an exact dressed-hexagon operator
`R_6^(8)` and permits boundary or allowed local-environment dependence. Only
the endpoint topology is promoted. No numerical hexagon correction is needed
to establish the new octagon coefficient or the pole boundary.

## 3. Hostile spectral checks

### E. Calling order eight the first interaction

**Attack.** Set the two-photon 1PI kernel identically to zero until an
order-eight operator appears.

**Result.** Rejected. The leading hard-core compact hexagon exchange is
already a nonquadratic microscopic many-body Hamiltonian. Gaussian Maxwell is
an infrared fixed point. The first unevaluated TT four-point matching
coefficient therefore belongs to the order-six pure-ice model; order eight
corrects it.

### F. Promoting a TT bare vertex to a pole

**Attack.** A nonzero TT projection of the octagon or dressed-hexagon vertex
is itself the graviton.

**Result.** Rejected at the strict single-insertion level. A bare finite-range
loop vertex is analytic in external lattice momenta, and one insertion does
not contain the geometric Bethe--Salpeter denominator. The fully dressed
massless 1PI function may nevertheless have threshold nonanalyticities, and
the finite-order Hamiltonian may bind or shift a pole if it is solved or
resummed nonperturbatively. The verifier checks only the narrower algebraic
fact that the Bethe--Salpeter denominator is absent from the strict
single-kernel truncation.

### G. Hiding a massless bound state below threshold

**Attack.** Place a stable positive-energy tensor branch below the two-photon
continuum.

**Result.** Not excluded by kinematics alone. At fixed nonzero total momentum
a slower collective branch could lie below `c|k|`, but establishing it
requires a nonperturbative spectral calculation and it fails the common-cone
gate until its velocity is shown to equal `c`. A common-cone pole lies at the
continuum edge, where protection and nonzero thermodynamic residue are
load-bearing. No rank-two constraint or Ward identity is derived at order
eight.

## 4. Derivability boundary

The lattice order-eight operator is derived, not fit. The missing physical
quantity is first the commonly normalized connected TT four-point function of
the existing pure-ice Hamiltonian. External-leg amputation and subtraction of
the pieces reducible in the selected two-photon channel then define the 2PI
matching coefficient `g_TT^match(mu_RG)` and its momentum-dependent kernel.
Two-point flux scaling and the Maxwell photon do not determine any of these.
Exact diagonalization, a sign-free projector calculation, or a controlled
linked-cluster calculation can advance this target without a laboratory.

## 5. Reproduction and ceiling

Run:

```text
python3 LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/verify_inherited_tt_kernel_boundary.py
```

This self-audit is not itself independent; the separate hostile disposition is
in `INDEPENDENT_AUDIT.md`. The packet earns no volume-uniform expansion,
matched four-point coefficient, bound-state pole, helicity-two
representation, Ward identity, universal stress coupling, RGRL-B, Einstein
response, gravity, or `G`.
