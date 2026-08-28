# Self-audit: Gaussian-Maxwell composite pole screen

**Lane:** `GRA-FL-F3-Q4-MCPS-V001`

**Audit date:** 2026-08-27

**Disposition:** `PASS_WITH_EXPLICIT_IMPORTED_MAXWELL_PREMISE`

## 1. Question actually audited

The audit asks only what infrared singularities the already physical FJ
one-link and pair operators can carry after the exact `d_*=2` ice projection
and after importing the leading pure-kinetic model's Gaussian Maxwell
infrared phase.  It does not ask whether six local coordinates can be called
a tensor, and it does not infer continuum spin from an `S4` label.

The exact/internal portion consists of:

- enumeration of the six local ice states;
- the incidence-corrected tetrahedral link-to-flux inverse and shared-edge
  gluing;
- the collapse of the pair module from `A1+E+T2` to constant `A1` plus
  nonconstant `E`;
- complement oddness of one-link observables and evenness of pairs; and
- the finite-group projector and character identities.

The thermodynamic `U(1)` liquid, emergent isotropy, and Gaussian Maxwell
description are imported through the explicit `MAXWELL-IR` premise.  The
public GFMC/ED result is evidence for that premise, not an internal proof of
the complete F3 parent's infinite-volume phase.

## 2. Hostile challenges and results

### Challenge A: carry all six FJ pair directions into ice

**Attack.**  Use the unprojected `A1+E+T2` pair count as six independent
infrared coordinates.

**Result.**  Rejected exactly.  For every ice state, opposite-edge products
are equal and the six-pair sum is `-2`.  The pair `T2` projector is killed,
the `A1` direction is constant in connected response, and only a
two-dimensional `E` variation survives.  The verifier checks these claims on
all six states and with the exact `A1/E/T2` projectors.

### Challenge B: interpret the surviving local `S4` sectors as spin two

**Attack.**  Since a continuum traceless symmetric tensor restricts as
`E+T2`, call the local sectors a graviton.

**Result.**  Rejected.  The same `T2` label is the full-tetrahedral
`T_d \simeq S4` restriction of a polar vector; the proper-rotation subgroup
alone is `A4`.  The actual one-link realization is linear electric flux.  The
verifier distinguishes rotation characters:
`2 cos(theta)` for transverse spin one and `2 cos(2 theta)` for helicity two.
An equivariant coordinate rearrangement cannot alter a pole's little-group
representation.

### Challenge C: hide a one-photon pole in the even pair channel

**Attack.**  Allow renormalization of `s_a s_b` to generate a term linear in
the Maxwell field.

**Result.**  Rejected under the stated zero-flux complement-symmetric
premise.  The microscopic pair is even, the Maxwell field and one-photon
state are odd, and the vacuum is even.  Therefore the one-photon matrix
element vanishes exactly.  A truncated Fock replay independently verifies
that an odd field reaches one particle while its normal-ordered even square
does not.

**Boundary.**  A nonzero background flux or broken complement symmetry can
linearize a bilinear around the background.  That borrows the same vector
photon pole; it does not create a helicity-two particle.  This case is
explicitly outside `MAXWELL-IR`.

### Challenge D: mistake the two-photon light-cone threshold for a pole

**Attack.**  Because both the photon pole and the two-photon threshold occur
at `|omega|=c|k|`, promote the threshold singularity to a new massless
particle.

**Result.**  Rejected.  Wick contraction of a local field-strength bilinear
is a momentum convolution.  The triangle inequality fixes the threshold,
while continuously varying the internal momentum fills the region above it.
This is branch-cut spectral weight, with an operator-dependent onset which
need not have one universal suppression exponent, rather than a delta
function with a TT one-particle residue.  The verifier separately checks the
threshold inequality and continuous two-particle kinematics.

### Challenge E: confuse equal dimension with equal residue type

**Attack.**  Both a transverse vector and a TT tensor have two physical
polarizations, so identify their rank-two residue spaces.

**Result.**  Rejected.  The link residue is the pullback of the vector
projector `P^T_ij`; a helicity-two pole requires the four-index projector
`Pi^TT_ij,kl`.  Their rotations have different characters even though both
projectors have rank two.  The verifier checks symmetry, idempotence,
transversality, trace removal, rank, and the distinct rotation laws.

### Challenge F: promote numerical phase evidence to an F3 theorem

**Attack.**  Treat published QMC at the pure-ice point as proof that every
all-orders thermodynamic realization of F3 lies at the same fixed point.

**Result.**  Rejected.  CROSS-CW earns the exact leading sixth-order model on
supplied plaquette-complete support and imports the public phase result.  A
volume-uniform all-orders expansion and phase-stability transfer remain open.
Every spectral theorem in this packet is labeled conditional on
`MAXWELL-IR`.

### Challenge G: use the negative screen to rule out gravity

**Attack.**  Conclude that no same-parent gravity emergence route exists.

**Result.**  Rejected.  The theorem excludes only the direct Gaussian claim
that the existing local link/pair composite already contains an isolated
helicity-two particle.  Inherited non-Gaussian higher-order dynamics could,
in principle, bind a protected even-channel mode.  A distinct same-parent
rank-two constrained phase is also logically open.  Both alternatives must
pass pole, residue, Ward/constraint, thermodynamic-residue, and universal
stress-coupling tests.

## 3. Dependency and independence audit

The FJ theorem and hostile audit and the CROSS-CW theorem and primary-source
binding are hash-pinned.  The separately owned local ice-sector theorem is
currently only a consistency cross-check; it is not load-bearing because
(FL03)-(FL09) rederive every local identity used here.  Its dependency hash
must not be frozen until its independent hostile audit is final.

No shared MODEL or register file is edited by this lane.  No new microscopic
interaction, fitted attraction, or external scalar is introduced.

## 4. Reproduction and verdict

Run:

```text
python3 LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/verify_maxwell_composite_pole_screen.py
```

Expected result before the post-audit ice-sector hash is pinned:
`Maxwell composite pole screen verification: 108 passed, 0 failed`.

The proof verdict is:

```text
ONE_LINK: ISOLATED MASSLESS SPIN-1 MAXWELL POLE (CONDITIONAL ON MAXWELL-IR)
CENTERED PAIR: TWO-PHOTON CONTINUUM + CONTACTS, NO ONE-PHOTON POLE
HELICITY-2: NO ISOLATED POLE AT THE GAUSSIAN FIXED POINT
NEXT: SAME-PARENT NON-GAUSSIAN TT KERNEL; IF NEGATIVE, DISTINCT RANK-2 WARD ARCHITECTURE
```

This is a substantive negative result: it removes a false short route and
leaves one bounded inherited-interaction calculation as the next direct
composite test.
