# GL6AW — anisotropic Følner twist-closure theorem

## Status and claim class

**Status:** author frozen and sealed; distinct independent hostile audit is
required before promotion.

This lane starts from sealed and independently audited `GL6AR`, `GL6AS`, and
`GL6AU`.  It derives an exact finite-size alternative to the unproved
`ICE0-STATIC` structure-factor bound.  On odd-transverse-area rectangular
tori, the centered half-filled port charge gives a nontrivial large-twist
translation character.  The result is an exact dichotomy: translation-
related ground-component degeneracy, or a same-component excitation whose
energy closes along an anisotropic three-dimensional Følner sequence.

This is not an isotropic phase theorem and not an infinite-volume GNS
gaplessness theorem.  It proves no low-character density weight, pole,
dispersion, physical momentum, cone, stress/Ricci law, gravity, or `G`.

## 1. Rectangular pure-loop model

Let

```text
Q_(L0,L1,L2)=Z/L0 x Z/L1 x Z/L2,
V=L0 L1 L2,
```

with all periods at least four.  Reduce the GL6AR incidence modulo these
periods, with

```text
d_0=(1,0,0), d_1=(0,1,0), d_2=(0,0,1), d_3=(0,0,0).
```

The locked configurations obey degree two at every parent and child node.
There are four elementary hexagon orientations per cell, one for each
unordered port triple.  The exact pure-loop Hamiltonian is

```text
H=-J sum_c T_c,
J=(63/8)h^6/U_d^5>0.                                  (AW.1)
```

Every flip preserves each port total

```text
N_a=sum_x n_(x,a).                                     (AW.2)
```

The total locked occupation is `sum_a N_a=2V`.  Hence the centered port-zero
charge is

```text
Q_0=N_0-V/2.                                           (AW.3)
```

This lane restricts only to the nonempty sector `Q_0=0`.

## 2. The centered sector is nonempty

Suppose `L1` is even.  One perfect matching chooses port zero at cells with
even `x1` and port three at cells with odd `x1`.  It is a perfect matching at
the child nodes as well: an incoming port-zero edge has the same `x1`, and
an incoming port-three edge starts at the child cell itself, so exactly one
of the two enters each child.  A second, edge-disjoint perfect matching uses
port one at every cell.

Their union is locked and has

```text
N_0=V/2.                                               (AW.4)
```

Thus the sector used below is not a formal empty sector.  This witness does
not select its controlling flip component or assert that the witness itself
is active.

## 3. Exact large twist and translation character

Let `C` be a connected flip component inside `Q_0=0`, and define

```text
q=2pi/L1,
W_0=sum_x x1 n_(x,0),
U_0=exp(iq W_0).                                       (AW.5)
```

`U_0` is diagonal in the locked configuration basis, so it preserves the
linear span of every connected component even though it does not commute
with `H`.

Let `Y` translate cells by one unit in the `x1` direction.  The wrap-plane
term in translating `W_0` is an integer multiple of `L1`, so its exponential
is one.  Therefore

```text
Y U_0 Y^-1=exp(-2pi i N_0/L1) U_0.                    (AW.6)
```

Take

```text
L0 and L2 odd,  L1 even.                              (AW.7)
```

In `Q_0=0`, equations (AW.3) and (AW.7) give

```text
exp(-2pi i N_0/L1)
 =exp(-pi i L0 L2)=-1.                                (AW.8)
```

This is the nontrivial character that is absent on even cubic half-filled
tori.

## 4. Same-component orthogonality

Assume first that `Y C=C`.  The component Hamiltonian is

```text
H_C=-J A_C.                                            (AW.9)
```

Its PF ground vector `psi_C` is unique and strictly positive.  Translation
permutes its positive entries and commutes with `H_C`, so uniqueness and
positivity force

```text
Y psi_C=psi_C.                                         (AW.10)
```

Equations (AW.6)--(AW.8) then give

```text
Y U_0 psi_C=-U_0 psi_C,
<psi_C,U_0 psi_C>=0.                                  (AW.11)
```

The twisted trial is therefore normalized, lies in the same connected
component, and is orthogonal to its unique ground vector.

## 5. Exact affected-term count and energy

In a hexagon containing port zero, its two port-zero links are separated by
the difference of the other two port displacements.  The four unordered
triples have the following `x1` twist changes:

```text
{0,1,2}: +/-1,
{0,1,3}: +/-1,
{0,2,3}: 0,
{1,2,3}: no port-zero link.                            (AW.12)
```

Modulo a wrap, `+/-(L1-1)` has the same exponential as `-/+1`.  Thus exactly
two orientations per cell, or `2V` terms, acquire phase `exp(+/-iq)`.

Put

```text
t_c=<psi_C,T_c psi_C>.
```

The PF entries and partial-flip matrix give `0<=t_c<=1`; no equal-amplitude
or regular-graph premise is used.  Pairing the two directions of every
configuration edge yields exactly

```text
<U_0 psi_C,H_C U_0 psi_C>-E_0(C)
 =J[1-cos(q)] sum_(c affected) t_c
 <=2JV[1-cos(2pi/L1)].                                (AW.13)
```

By (AW.11), the min--max principle applies inside the same component:

```text
boxed:
Delta_C
 <=2JV[1-cos(2pi/L1)]
 <=4pi^2 J L0 L2/L1.                                  (AW.14)
```

Equation (AW.14) is exact for every translation-stable component satisfying
the displayed arithmetic and centered-charge premises.

## 6. Ground-component dichotomy without assuming component stability

Let `C_*` be a component with minimum energy among all components in the
finite `Q_0=0` sector.  Translation sends `C_*` to a component with the same
energy.

- If `Y C_*` differs from `C_*`, the centered sector has at least two exactly
  degenerate, orthogonal ground components.
- If `Y C_*=C_*`, equations (AW.11)--(AW.14) give a same-component
  orthogonal excitation with the stated energy bound.

> **Finite centered-sector twist dichotomy.** On every rectangular torus
> satisfying (AW.7), the `Q_0=0` sector has either translation-related exact
> ground-component degeneracy or a controlling translation-stable component
> with gap bounded by (AW.14).

This is a sector theorem.  It does not prove that the actual all-sector
ground is centered in port zero.

## 7. An exact three-dimensional closing sequence

For odd integers `m>=5`, choose

```text
(L0,L1,L2)=(m,2m^3,m).                                (AW.15)
```

All three periods tend to infinity and the injectivity radius tends to
infinity.  The corresponding rectangular boxes form a Følner/van Hove
sequence.  Their face-area-to-volume ratio is

```text
2(1/L0+1/L1+1/L2) -> 0.                               (AW.16)
```

For every fixed lattice collar width, the actual boundary-site fraction is
bounded by that width times the same vanishing sum.

For every translation-stable controlling component in the dichotomy,
(AW.14) becomes

```text
Delta_C(m)<=2pi^2 J/m ->0.                             (AW.17)
```

Thus the pure model has an exact anisotropic finite-size obstruction to a
unique uniformly isolated centered-sector ground: along (AW.15), there is
either exact component degeneracy or a same-component closing gap.

## 8. What the theorem does not close

The route bypasses the `ICE0-STATIC` lower bound only for this finite-size
anisotropic conclusion.  It does not prove that

```text
S_T(2pi/L)>=s/L,
```

or any isotropic gap law.  The twist gradient is small, but its total lattice
translation character is `-1`; this is an LSM-like obstruction, not a
derived small-character collective dispersion.

For an infinite-volume GNS conclusion one must still transport the twisted
trials through compatible selected states and prove that their normalized
limits remain outside the full zero-energy projection.  Otherwise they may
become a symmetry-related, topological, or other zero-energy sector while a
selected pure GNS phase retains a positive internal gap.  This is precisely
the finite/GNS ceiling retained in GL6AR and GL6AU.

The theorem also uses the exact pure order-six Hamiltonian.  It proves no
stability under unknown higher-order finite-`h/U_d` terms.

## 9. Final result

```text
odd transverse area + centered port charge
  -> exact twist character -1;
translation-stable controlling component
  -> same-component gap <=4pi^2 J A_perp/L_parallel;
otherwise
  -> exact translation-related ground-component degeneracy.
```

This is an unconditional finite-torus dichotomy within the declared
centered sector and an exact closing result along the declared anisotropic
Følner sequence.  It is not a selected-GNS mode, isotropic phase proof,
physical cone, stress/Ricci response, gravity, or `G`.
