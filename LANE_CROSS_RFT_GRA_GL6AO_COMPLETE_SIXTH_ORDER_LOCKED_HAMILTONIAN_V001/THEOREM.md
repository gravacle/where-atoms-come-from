# GL6AO — COMPLETE SIXTH-ORDER LOCKED HAMILTONIAN THEOREM

## Status and exact scope

This packet starts only from the sealed and independently audited GL6AN
degree-lock result.  It computes the complete canonical effective Hamiltonian
through order `h^6/U_d^5` on GL6AN's declared period-four, girth-at-least-six
quotient `Q_4`.  It then extracts the finite-range linked interaction earned
at that order.

The theorem proves an exact native collective operation.  It does **not**
derive or assume a photon, graviton, gauge constraint, continuum momentum,
pole, phase, physical cone, Ricci tensor, Einstein equation, gravity, or
Newton's constant `G`.

## 1. Locked parent and perturbative convention

On `Q_4`, let

```text
D = sum_v (k_v-2)^2,
W = -sum_e X_e,
H = U_d D + h W + C,
P = 1_{D=0},
Q = 1-P,
R = -Q D^{-1} Q.
```

Here `M=|E(Q_4)|=256`, `|V(Q_4)|=128`, and the spectral gap from `P` to
`Q` at `h=0` is `2U_d`.  The finite-volume canonical Kato parallel-transport
effective Hamiltonian is analytic near `h/U_d=0`.  Write

```text
H_eff = C P
      + (h^2/U_d) K_2
      + (h^4/U_d^3) K_4
      + (h^6/U_d^5) K_6
      + O(h^8/U_d^7).
```

All odd orders through seven vanish because a locked-to-locked flip set is an
even subgraph of the bipartite incidence graph.  There is no nonempty such
set of odd cardinality.

The first nonscalar term occurs only at order six.  Consequently it is
unchanged by analytic basis changes within `P`: commutators with the order-zero,
order-two, and order-four scalar operators vanish.  The order-six operator
below is therefore not a basis-gauge artifact.

## 2. Lower-order scalar data

Define

```text
B   = P W R W P,
A_2 = P W R^2 W P,
A_3 = P W R^3 W P,
T_4 = P W R W R W R W P.
```

The sealed GL6AN census, replayed here, gives

```text
B   = -(M/2) P,
A_2 =  (M/4) P,
A_3 = -(M/8) P,
T_4 = -((3M^2+7M)/24) P.
```

The canonical fourth-order fold is

```text
K_4 = T_4 - (1/2){A_2,B} = -(7M/24)P,
```

while

```text
K_2 = -(M/2)P.
```

Thus there is no configuration-changing operator through fourth order.

## 3. Complete sixth-order formula

Let

```text
T_6 = P W R W R W R W R W R W P,

X_4 = P W R^2 W R W R W P
    + P W R W R^2 W R W P
    + P W R W R W R^2 W P.
```

Because `K_2=bP` and `K_4=dP` are scalar, with

```text
b = -M/2,
d = -7M/24,
```

the sixth-order intermediate-normalized Kato recursion reduces exactly to

```text
K_6 = T_6 - b X_4 + b^2 A_3 - d A_2.                 (AO.1)
```

The sum `X_4` is Hermitian: its first and third summands are adjoints and its
middle summand is self-adjoint.  Equation (AO.1) therefore gives the same
first nonscalar order as the canonical Hermitian Kato transport.

## 4. The only order-six configuration-changing term

Let `s` and `s'` be distinct degree-two locked occupation configurations.
Their symmetric difference has even degree at every constraint node.  Any
nonempty finite even subgraph is a union of cycles.  Since `Q_4` has no
two-cycle or four-cycle, a difference reachable in at most six flips must
contain exactly six distinct links and must be one six-cycle.

At each vertex of that cycle, preserving degree two requires one removed and
one inserted link.  The six occupations therefore alternate.  Conversely,
toggling an alternating six-cycle preserves every constraint.  Hence:

> **AO configuration-change classification.** Every nonzero off-diagonal
> matrix element through order six is an alternating six-cycle toggle, and
> every alternating six-cycle toggle occurs at order six.

All folded terms in (AO.1) contain at most four flips between locked states
and are diagonal.  The off-diagonal coefficient therefore comes only from
`T_6`.  Exact enumeration of all `6!=720` orders gives

```text
<s'|K_6|s> = -63/8                                  (AO.2)
```

when `s xor s'` is one alternating six-cycle, and zero for every other pair.

The declared `Q_4` contains exactly 256 undirected six-cycles, and every link
belongs to exactly six of them.

## 5. Complete diagonal order-six census

A six-flip diagonal word has even multiplicity on every used link.  A
Q-only word cannot use one link alone, so the direct term `T_6` has exactly
two types:

1. one link four times and another twice;
2. three links twice each.

For a two-link parity subset, its excitation energy in units of `U_d` is

```text
p=2  adjacent links with opposite initial occupations,
p=6  adjacent links with equal initial occupations,
p=4  disjoint links.
```

The complete repeated-pair word sums, including both choices of the
four-times link, are

| `p` | direct six-word weight |
|---:|---:|
| 2 | `-1/4` |
| 4 | `-1/16` |
| 6 | `-1/36` |

For three links used twice each, the exact Q-only word sum depends on the
three pair energies and the triple energy `t`:

| graph/occupation class | `(p_12,p_13,p_23;t)` | count | word sum |
|---|---:|---:|---:|
| matching | `(4,4,4;6)` | `C(M,3)-3M^2+19M` | `-9/32` |
| one adjacent, opposite | `(2,4,4;4)` | `2M(M-10)` | `-9/16` |
| one adjacent, equal | `(4,4,6;8)` | `M(M-10)` | `-29/144` |
| three-edge star | `(2,2,6;4)` | `2M` | `-109/144` |
| path, both joins opposite | `(2,2,4;2)` | `4M` | `-41/32` |
| path, one join equal | `(2,4,6;6)` | `4M` | `-337/864` |
| path, both joins equal | `(4,6,6;10)` | `M` | `-209/1440` |

These counts are independent of the chosen locked configuration.  They
follow only from degree four, girth at least six, and two occupied links at
every node.  In particular, every middle link has two opposite-occupation
and one equal-occupation continuation at each endpoint.

Summing the table gives

```text
T_6|diag = [
    -(3/64) M^3
    -(215/576) M^2
    -(893/1080) M
] P.                                                   (AO.3)
```

The three folded four-word sums give

```text
X_4 = [(5/32)M^2 + (173/288)M] P.                    (AO.4)
```

Using `A_3=-(M/8)P`, `b=-M/2`, `d=-7M/24` in (AO.1),

```text
-bX_4   =  (5/64)M^3 + (173/576)M^2,
b^2 A_3 = -(1/32)M^3,
-d A_2  =  (7/96)M^2.
```

The `M^3` and `M^2` terms cancel exactly against (AO.3).  Therefore

```text
<s|K_6|s> = -(893/1080) M                           (AO.5)
```

for **every** locked configuration `s`.

There is no configuration-dependent diagonal term at order six and, in
particular, no order-six diagonal potential proportional to the number of
flippable six-cycles.

## 6. Finite-volume theorem

For each undirected six-cycle `c` in `Q_4`, define the self-adjoint locked
toggle

```text
T_c = P (product_{e in c} X_e) P.
```

It annihilates a nonalternating locked configuration and toggles an
alternating one.  Combining (AO.2) and (AO.5) yields the complete result:

```text
H_eff = C P
  - (M/2)       (h^2/U_d)   P
  - (7M/24)     (h^4/U_d^3) P
  - (h^6/U_d^5) [
        (893M/1080) P
        + (63/8) sum_{c in Hex(Q_4)} T_c
    ]
  + O(h^8/U_d^7).                                    (AO.6)
```

Equation (AO.6) includes every diagonal and off-diagonal connected
contribution through sixth order.  The cancellation of all superextensive
terms is explicit, rather than assumed from a linked-cluster slogan.

## 7. Earned thermodynamic linked interaction

No global locked projector is inserted into the infinite quasi-local
algebra.  For an elementary six-cycle `c` of the infinite `A3` incidence,
let

```text
p_v = 1_{k_v=2},
P_c = product_{v in V(c)} p_v,
tau_c = P_c (product_{e in c} X_e) P_c.
```

This is a bounded, finite-support, self-adjoint operator.  On a globally
locked configuration it is exactly the alternating-cycle toggle.  The
nontrivial formal sixth-order interaction is

```text
Phi_c^(6) = -(63/8)(h^6/U_d^5) tau_c.                (AO.7)
```

Every link belongs to six elementary cycles, so (AO.7) is uniformly
finite-range and has finite interaction norm.  It therefore defines a
well-typed quasi-local linked interaction at the computed order.  The scalar
energy density through order six is

```text
e_scalar/link =
  -(1/2)(h^2/U_d)
  -(7/24)(h^4/U_d^3)
  -(893/1080)(h^6/U_d^5).
```

The scalar density does not affect the finite-order dynamics.

## 8. Exact ceiling and next gate

GL6AO proves the complete first native collective Hamiltonian generated by
the locked parent.  It does **not** prove:

- convergence of the all-orders Schrieffer-Wolff/Kato series uniformly in
  volume;
- selection or existence of a thermodynamic phase;
- a gapless pole, physical momentum, propagation cone, or universal stress
  coupling;
- a gauge photon or graviton;
- Ricci/Einstein form, gravity, or `G`.

Those are separate infrared and physical-identification gates.  GL6AO earns
their microscopic input without presupposing their answer.
