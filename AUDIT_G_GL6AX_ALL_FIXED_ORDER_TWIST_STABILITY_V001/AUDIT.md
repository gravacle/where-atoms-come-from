# Distinct hostile audit — GL6AX all-fixed-order twist stability

**Target:** `LANE_CROSS_RFT_GRA_GL6AX_ALL_FIXED_ORDER_TWIST_STABILITY_V001/`  
**Frozen theorem SHA-256:** `bbdeffc955c5dbc2ce1bec00578ab925e531413188f47893e867b99c6ff6c56d`  
**Frozen author-manifest SHA-256:** `e10907affe77c84ae175aa60a6ae60429d3393ea117be32c81ede45e23357f36`  
**Frozen author-seal-file SHA-256:** `619234ef81e5d739aa7c166d567844ac1c2773b5b316bce64edf57c0f87e3263`  
**Disposition:** `PASS__EXACT_ALL_FIXED_ORDER_PORT_CONSERVATION_AND_ANISOTROPIC_TWIST_STABILITY__FINITE_COUPLING_AND_PHYSICAL_PROMOTION_OPEN`

## 1. Independence and custody

All eleven author files are byte-pinned in `AUDITED_TARGETS.sha256`.  The
independent replay imports no author module.  Both author verifiers and both
independent audit verifiers pass in normal and optimized Python modes.  The
eighteen frozen GL6AN, GL6AO, and GL6AW author/audit dependencies resolve,
including every pinned manifest and seal.

The audit began from the frozen author bytes and did not alter the author
packet.

## 2. Affine-moment law

For a finite difference `delta_(x,a)` between two locked configurations, the
parent and child equations are

```text
sum_a delta_(x,a)=0,
sum_a delta_(y-d_a,a)=0.
```

Summing gives `sum_a Delta N_a=0`.  Multiplying the child equation by the
lifted coordinate `y_j`, substituting `y=x+d_a`, and subtracting the
coordinate-weighted parent equation gives

```text
sum_a (d_a)_j Delta N_a=0,  j=0,1,2.
```

The four columns `(1,d_a)` are the vertices of a unimodular affine simplex:

```text
det[(1,d_0) (1,d_1) (1,d_2) (1,d_3)] = -1.
```

Therefore the only solution is

```text
Delta N_0=Delta N_1=Delta N_2=Delta N_3=0.
```

This proof never decomposes the symmetric difference into a single simple
cycle.  It consequently covers degree-four vertices, self-touching even
subgraphs, disconnected unions, and simultaneous multi-loop changes.  The
independent replay additionally enumerates closed alternating parent walks
and constructs a locked degree-four/self-touching periodic example to make
sure no hidden simple-cycle premise entered.

## 3. Periodic seams and the exact `2L_min` threshold

On a torus, using coordinate representatives `0,...,L_j-1`, a changed
port-`j` edge crossing the positive seam contributes `1-L_j` rather than
`1` to the coordinate moment.  The exact identities are therefore

```text
Delta N_j=L_j W_j,                     j=0,1,2,
Delta N_3=-sum_j L_j W_j.
```

The sign changes if the seam orientation convention is reversed, but the
conservation and threshold statements do not.  Zero winding restores the
finite affine law; nonzero winding is the only exception.

Every alternating parent-to-parent step uses two changed links and has
coordinate increment `d_a-d_b`, whose components have magnitude at most
one.  A cycle with nonzero winding in direction `j` therefore has at least
`L_j` parent steps and at least `2L_j` changed links.  Decomposing an even
signed difference into alternating circuits shows that one winding circuit
must exist, so every winding difference has at least `2L_min` links.

The lower bound is attained.  In a uniform locked configuration containing
ports `j` and a spectator port, replacing `j` by `3` around a periodic
`x_j` row preserves both endpoint locks and changes

```text
Delta N_j=-L_j,
Delta N_3=+L_j
```

on exactly `2L_j` links.  The independent replay checks this in all three
directions, checks the seam formula directly, and constructs two crossing
winding cycles whose symmetric difference has a degree-four parent.  The
sharp minimum is thus exactly `2L_min`, not merely an upper or lower bound.

## 4. Fixed-order SW/Kato/Feshbach scope

An order-`r` microscopic word contains exactly `r` one-link flip factors;
diagonal unperturbed resolvents and projectors do not change the final bit
parity.  Hence the Hamming distance between any two locked endpoints of a
nonzero order-`r` matrix element is at most `r`.  If `r<2L_min`, section 3
excludes winding and section 2 forces equality of all four port totals.

This reasoning survives repeated flips, folded terms, intermediate returns
to the locked space, disconnected endpoint differences, complex
coefficients, and the coefficient-by-coefficient expansion of the Feshbach
self-energy.  It proves the displayed charge blocks in the canonical
physical-basis word representation.  An arbitrary extra unitary rotation
inside the locked subspace can hide those blocks by changing the displayed
basis; GL6AX expressly excludes that gauge convention rather than claiming
an invariant matrix form under every locked-subspace rotation.

The argument does not infer a linked interaction decomposition from a raw
reducible Feshbach expression.  Locality or controlled quasi-locality is a
separate premise of the twist theorem.  Nor does fixed-order conservation
imply exact finite-coupling conservation: words of order `2L_min` can wrap.

## 5. Termwise `U(1)` symmetrization and constrained locality

If `H=sum_X Phi(X)` commutes with global `N_0`, replace each interaction term
by

```text
bar Phi(X)=(1/2pi) integral dtheta
           exp(i theta N_0) Phi(X) exp(-i theta N_0).
```

The sum remains `H`.  Since occupation operators outside `X` commute with a
term supported in `X`, the average is generated by `N_(0,X)`.  It therefore
preserves support, translation covariance, Hermiticity, and norm while
making every term commute with its supported charge.  Projection to the
locked configuration span is diagonal and does not invalidate this
support-commutation statement.

For a support of cyclic `x1` diameter below `L1/2`, choose a lift and an
origin.  Changing a lifted coordinate by `L1` has unit twist phase, and
subtracting the origin contributes `y_X N_(0,X)`, which commutes with the
symmetrized term.  Thus the global twist acts on that term exactly through
the author's local generator `A_X`; no seam or origin term was omitted.

## 6. `+/-` twist energy and constants

For

```text
g_X(s)=exp(-is A_X) Phi(X) exp(is A_X),
```

twice differentiating gives an isometric conjugate of
`-[A_X,[A_X,Phi(X)]]`.  The symmetric Taylor remainder therefore obeys

```text
||(g_X(q)+g_X(-q))/2-g_X(0)||
 <= (q^2/2)||[A_X,[A_X,Phi(X)]]||.
```

This operator-norm proof is insensitive to complex amplitudes.  All odd
terms, including a nonzero equilibrium-current expectation, cancel between
the two twist directions.  The independent replay tests arbitrary complex
two-state hopping amplitudes, arbitrary diagonal energies, and both signs
of every position difference.

With `q=2pi/L1`, `V=L0 L1 L2`, and

```text
D_2(L)=(1/V) sum_X ||[A_X,[A_X,Phi(X)]]||,
```

the average excess energy is at most

```text
(q^2/2)V D_2(L)=2pi^2 D_2(L)L0L2/L1.
```

In the centered sector, `N_0=V/2` and odd `L0L2` give the exact translation
character `-1` for both `U psi_0` and `U^-1 psi_0`.  If the centered-sector
ground is unique, both trials are orthogonal to it; otherwise the degenerate
branch already holds.  Each trial energy is at least the sector ground
energy, and at least one is no larger than their average.  Min--max gives

```text
either centered-sector ground degeneracy,
or Delta_L <= 2pi^2 D_2(L)L0L2/L1.
```

The coefficient and factor of two are correct.  Configuration-dependent
diagonal terms commute with the twist and contribute zero.  Local
simultaneous multi-loop terms and complex Hermitian amplitudes are included
through the same double commutator.

## 7. Interaction-tail and order-of-limits attacks

For a wrapping interaction term, the difference between any twisted and
untwisted expectation is at most `2||Phi(X)||`.  Hence

```text
T_L=sum_(X wrapping) 2||Phi(X)||
```

is the correct tail allowance, not a missing factor four.  A uniformly
exponentially decaying interaction with finite second moment makes the
wrapping tail exponentially small up to the polynomial volume count.  A
half-system transfer acquires phase `pi`; algebraic tails with divergent
second moment, all-to-all terms, and nonvanishing wrapping norm are genuine
counterexamples and are expressly excluded.

For every fixed truncation `K`, a translation-covariant linked interaction
with finite range `R_K` and finite `D_(2,K)` obeys the twist theorem once

```text
K<2L_min,  2R_K<L1.
```

Along `(L0,L1,L2)=(m,2m^3,m)` this gives

```text
Delta_L^(<=K) <= pi^2 D_(2,K)/m -> 0
```

on the unique-ground branch.  The order is fixed before the volume limit.
Nothing in the proof controls growth of `R_K` or `D_(2,K)`, convergence as
`K->infinity`, or a finite-coupling dressed locked projector.  Reversing or
interchanging those limits is not licensed.

## 8. Promotion attacks and verdict

The exact theorem is narrower than several tempting readings:

1. It does not select the centered port sector as the physical ground.
2. Degenerate twisted states may enter the full zero-energy projection; a
   compatible selected-state/GNS bridge remains open.
3. The anisotropic Følner obstruction is not an isotropic small-character
   dispersion, pole, or common physical cone.
4. The finite twisted trial is not record-authenticated merely because its
   energy closes.
5. The result is order-by-order, not a uniform finite-`h/U_d` theorem.
6. No photon, graviton, Ricci tensor, Einstein equation, gravity law, or
   numerical `G` is inserted or derived.

Within those exact ceilings, every hostile attack survived.

**Hostile verdict: PASS.**
