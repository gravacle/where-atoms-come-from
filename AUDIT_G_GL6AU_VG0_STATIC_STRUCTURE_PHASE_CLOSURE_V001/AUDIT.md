# Distinct hostile audit — GL6AU `v/g=0` static-structure closure

**Target:** `LANE_CROSS_RFT_GRA_GL6AU_VG0_STATIC_STRUCTURE_PHASE_CLOSURE_V001/`  
**Frozen theorem SHA-256:** `ce2cdf1394053778300c71f3b2d25f79914efe82390b5b66f84fde807e0f2612`  
**Frozen author-manifest SHA-256:** `347bfade725e1cbefbe6d205f596ea059e550d98f08aa452591bb47a45ae9705`  
**Frozen author-seal-file SHA-256:** `3a55ca000cc044305d5daae932a0c935cc73131f4c3f55db716d82753d1f7a4f`  
**Disposition:** `PASS__FIRST_CHARACTER_STATIC_EXPONENT_CLOSURE_EXACT__ICE0_STATIC_AND_GNS_BRIDGE_EXPLICITLY_UNPROVED__NO_PHYSICAL_CONE_STRESS_GRAVITY_OR_G`

## 1. Independence and custody

All twelve frozen author files are pinned in `AUDITED_TARGETS.sha256`.  The
independent replay imports no author module.  Author and independent scripts
were run in normal and optimized Python modes.  The dependency chain resolves
exactly to sealed GL6AO, GL6AR, GL6AS, and the complete author-plus-audit
GL6AT snapshot.

## 2. First-character plane and oscillator coefficient

At `z=(exp(iq),1,1,1)`, the exact transverse plane is

```text
u=(0,u1,u2,u3),  u1+u2+u3=0.
```

Reconstructing all four cycle columns gives

```text
sum_d |u^* C_d(z)|^2
 = 3|1-exp(iq)|^2
 = 12 sin^2(q/2)
```

for every normalized vector in that plane.  A separate two-state expansion
reproduces the per-flip double-commutator factor

```text
(J/2)<T_c>|delta rho|^2.
```

PF positivity gives `0<=<T_c><=1` without orientation symmetry.  With
`q=2pi/L`, these facts yield the exact numerator bound

```text
f_u(q)<=6J sin^2(pi/L).
```

The nontrivial translation character makes the density trial orthogonal to
the unique same-component PF ground vector.  Division occurs only when its
strictly positive spectral weight is nonzero.  Thus the author's

```text
Delta_C(L)<=6J sin^2(pi/L)/S_(u,L)(q)
```

is a valid upper bound on the gap of the same connected component.

## 3. Static exponent and quadrature scaling

For `S_(u,L)(2pi/L)>=sL^-alpha`, the sine bound gives

```text
Delta_C(L)<=(6pi^2J/s)L^(alpha-2).
```

The threshold is strictly `alpha<2`; no conclusion is made at `alpha=2`.
Translation invariance at the nontrivial character gives exactly

```text
Var(F_c)+Var(F_s)=L^3 S_(u,L),
```

so `alpha=1` requires only an `L^2` lower bound in at least one real
quadrature.  The author correctly repairs the earlier sufficient-but-too-
strong extensive-variance target.

## 4. PF transform and the missing measure theorem

For `H_C=-JA_C`, the PF transform

```text
P_C(n,m)=A_C(n,m)psi(m)/(rho_C psi(n)),
pi_C=psi^2
```

is stochastic and reversible, and

```text
Delta_C=J rho_C gap(P_C).
```

The independent replay verifies the corresponding Dirichlet identity on an
unbounded family of nonregular star graphs.  This establishes where the
unproved datum lives: the low-character variance of the nonuniform PF law.
It does not manufacture a lower bound for that variance.

## 5. Shortcut and counterexample attacks

At the RK point the equal-amplitude measure belongs to `J(D-A)`, whereas
the target operator is `-JA`.  Frozen degree-zero configurations and an
active configuration of degree at least `L^3/64` imply

```text
inf_c ||D-cI|| >= L^3/128.
```

This rules out treating the degree potential as a uniformly bounded small
perturbation.  GL6AU correctly stops short of ruling out every possible
special-component comparison theorem.

A sector tower can close while the selected-sector excitation gap remains
positive, so it is not a GNS proof.  Likewise the continuous-link
Wilson/Villain massless theorems have no exact transfer-measure map to this
projected spin-`1/2` adjacency Hamiltonian.

On even cubic zero-port tori, the one-port twist has phase
`exp(pi i L^2)=1`, exactly as stated.  The packet also correctly notes that
odd-transverse-area rectangular tori restore a nontrivial twist but do not
close the elementary bound in the isotropic three-dimensional scaling.

## 6. Post-freeze independent extension, not an AU defect

The hostile replay found a distinct exact route that was not required by any
AU claim.  On rectangular zero-port tori, let the twist length be even and
the transverse area odd.  In a translation-invariant connected PF component
the one-port twisted state has character `-1`, is orthogonal inside the same
component, and only two hexagon orientations per cell acquire the unit twist
phase.  Hence

```text
Delta_C <= 2J V[1-cos(2pi/L_parallel)]
        <= 4pi^2 J A_perp/L_parallel.
```

Choosing side lengths `(m,2m^3,m)` with odd `m` makes all three lengths tend
to infinity and gives `Delta_C<=2pi^2J/m`.  This is an anisotropic Følner
finite-component closure that avoids a structure-factor hypothesis.  It
does not prove isotropic finite-size scaling or selected-GNS gaplessness:
the twist sequence may enter the full zero-energy/topological sector unless
the separate GNS bridge is established.  This extension belongs in a new
lane, not as a silent modification of frozen GL6AU.

## 7. Evidence and promotion ceiling

The packet cleanly separates exact algebra from published finite-size QMC
and fitted Gaussian evidence.  `ICE0-STATIC` and `ICE0-GNS-BRIDGE` are named
hypotheses and are explicitly not counted as theorems.  The evidence may
support the expected `alpha=1` law, but it is not used as an all-`L`
inequality.

No translation character is called calibrated physical momentum.  No pole,
physical cone, stress/Ricci law, gravity, or `G` is derived.

**Hostile verdict: PASS.**
