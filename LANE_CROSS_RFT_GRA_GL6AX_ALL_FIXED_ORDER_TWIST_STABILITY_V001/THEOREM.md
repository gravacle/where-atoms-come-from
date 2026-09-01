# GL6AX — all-fixed-order port conservation and twist-stability theorem

## Status and claim class

**Status:** author frozen and sealed; distinct independent hostile audit is
required before promotion.

This lane starts from sealed and independently audited `GL6AN`, `GL6AO`, and
`GL6AW`.  It proves that the `GL6AW` anisotropic finite-size obstruction is
not peculiar to the isolated order-six hexagon Hamiltonian.  Affine
independence of the four native `A3` port displacements forces every finite,
contractible locked-to-locked change to conserve all four port totals.  As a
result, every fixed perturbative order below the wrapping threshold inherits
the port charge used by the large twist.  A `+/-` twist average then gives an
`O(area/length)` variational bound for any Hermitian, translation-invariant,
finite-range locked Hamiltonian with uniformly bounded local second twist
moment.  No reality, Perron--Frobenius, time-reversal, equal-amplitude, or
pure-hexagon premise is used.

The theorem is order-by-order.  It does not prove uniform convergence of the
Schrieffer--Wolff/Kato/Feshbach expansion at finite `h/U_d`, selection of the
centered sector, a selected infinite-volume GNS mode, an isotropic cone,
gravity, or `G`.

## 1. Native incidence and the affine-moment lemma

Write the infinite `A3` parent/child incidence as

```text
edge (x,a): parent x -> child x+d_a,
d_0=(1,0,0), d_1=(0,1,0), d_2=(0,0,1), d_3=(0,0,0).       (AX.1)
```

Let `n` and `n'` be degree-two locked configurations whose difference

```text
delta_(x,a)=n'_(x,a)-n_(x,a)                              (AX.2)
```

has finite support.  The parent and child locks give, respectively,

```text
sum_a delta_(x,a)=0,
sum_a delta_(y-d_a,a)=0.                                  (AX.3)
```

Put `Delta N_a=sum_x delta_(x,a)`.  Summing the first equation and then
subtracting its coordinate-weighted sum from the coordinate-weighted child
equation gives

```text
sum_a Delta N_a=0,
sum_a d_a Delta N_a=0.                                    (AX.4)
```

The four augmented columns `(1,d_a)` are affinely independent; their
`4 x 4` determinant has magnitude one.  Therefore

```text
boxed: Delta N_0=Delta N_1=Delta N_2=Delta N_3=0.          (AX.5)
```

This is stronger than an elementary-loop count.  A finite symmetric
difference can be a self-touching even subgraph, a disjoint union of
alternating cycles, or a simultaneous multi-loop change.  Equation (AX.5)
uses only the two endpoint locks and finite support, so all of those cases
are covered.

The same result holds on a periodic quotient whenever the symmetric
difference has an injective finite lift, or more generally has zero winding.
That is the precise meaning of *contractible* in this theorem.

## 2. Exact wrapping exception and sharp threshold

On

```text
Q_L=Z/L0 x Z/L1 x Z/L2,                                  (AX.6)
```

with every period at least four,

the coordinate moments acquire seam terms.  With `W_j` the signed number of
changed port-`j` links crossing the `x_j` seam, the exact periodic identity is

```text
Delta N_j=L_j W_j,                 j=0,1,2,
Delta N_3=-sum_(j=0)^2 L_j W_j.                           (AX.7)
```

Thus affine independence does not conserve port totals across a genuinely
noncontractible update.

This exception is realized, not merely algebraically possible.  Start from
the uniform locked configuration containing ports `j` and `b` at every
parent, where `b` differs from `j` and `3`.  Along one periodic `x_j` row,
toggle the alternating `j/3` cycle.  Its length is `2L_j` and

```text
Delta N_j=-L_j,
Delta N_3=+L_j.                                           (AX.8)
```

Conversely, each two-edge parent-to-parent step changes any lifted coordinate
by at most one.  Any nonzero winding symmetric difference therefore contains
at least `2L_min` changed links, where

```text
L_min=min(L0,L1,L2).                                      (AX.9)
```

The minimum wrapping Hamming distance is consequently exactly `2L_min`.

## 3. Fixed-order Schrieffer--Wolff/Kato/Feshbach consequence

Let the microscopic perturbation be the native sum of one-link flips.  A
word of order `r` changes the parity of at most `r` links.  Every nonzero
matrix element between two locked endpoint configurations therefore has
symmetric-difference size at most `r`.  Diagonal elements conserve every
port total trivially.  Equations (AX.5)--(AX.9) imply:

> **Fixed-order conservation theorem.** Every order-`r` locked-to-locked
> coefficient in the canonical word-based Hermitian Schrieffer--Wolff/Kato
> representation used by `GL6AO`, and every order-`r` coefficient in the
> expanded Feshbach self-energy, commutes with every `N_a` whenever
> `r<2L_min`.

An arbitrary additional unitary gauge rotation inside the locked subspace
could obscure this block structure by definition; it is not part of the
canonical physical-basis statement.

Folded terms, complex coefficients, repeated flips, intermediate returns to
the locked sector, and disconnected final symmetric differences do not alter
this conclusion: only the locked endpoints and the number of odd final
flips enter the proof.

There are two separate limitations.

1. Conservation of matrix elements does not by itself prove that a chosen
   effective-Hamiltonian representation has a uniformly local linked
   interaction.  A raw energy-dependent Feshbach expression can contain
   reducible or disconnected contributions.  The twist bound below requires
   a finite-range or controlled quasi-local interaction representation.
2. The exact finite-volume resummation contains arbitrarily high orders.
   Terms first appearing at order `2L_min` can wrap and violate `N_a` as in
   (AX.8).  No exact finite-`h/U_d` port symmetry is claimed for the original
   microscopic Hamiltonian.

## 4. General locked Hamiltonian and centered sector

Take the same arithmetic as `GL6AW`:

```text
L0,L2 odd; L1 even; V=L0 L1 L2; N_0=V/2.                 (AX.10)
```

Let `H_L` be any Hermitian, translation-invariant Hamiltonian on the full
locked Hilbert space that commutes with `N_0`.  Assume it has an interaction
decomposition

```text
H_L=sum_X Phi_L(X)                                      (AX.11)
```

of range at most `R<L1/2`.  It is no restriction to take each term to
commute with the port-zero occupation in its own support: average each term
over the global `U(1)` generated by `N_0`.  This preserves the sum, support,
translation covariance, Hermiticity, and does not increase its norm.

For each support `X`, lift its `x1` coordinates into an interval of length at
most `R`, choose any origin `y_X`, and put

```text
A_X=sum_((x,0) in X) (x1_tilde-y_X)n_(x,0).               (AX.12)
```

Changing the lift by `L1` does not change the unit twist, and subtracting the
origin does not change conjugation because `[Phi_L(X),N_(0,X)]=0`.  Define
the local second-twist density

```text
D_2(L)=(1/V) sum_X ||[A_X,[A_X,Phi_L(X)]]||.              (AX.13)
```

For a uniformly bounded finite-range interaction, `sup_L D_2(L)=D_2` is
finite.

## 5. `+/-` twist average without a current assumption

Let

```text
q=2pi/L1,
W_0=sum_x x1 n_(x,0),
U=exp(iqW_0).                                             (AX.14)
```

Translation `Y` along `x1` obeys, exactly as in `GL6AW`,

```text
Y U Y^-1=exp(-2pi iN_0/L1)U=-U.                          (AX.15)
```

The same character `-1` holds for `U^-1`.  If the centered-sector ground is
unique, translation invariance makes it a `Y` eigenvector, so both
`U psi_0` and `U^-1 psi_0` are orthogonal to it.  If the ground is not unique,
the alternative in the theorem is already satisfied.

For a single local interaction term define

```text
g_X(s)=exp(-isA_X) Phi_L(X) exp(isA_X).                   (AX.16)
```

Taylor's formula with integral remainder gives

```text
||(g_X(q)+g_X(-q))/2-g_X(0)||
 <=(q^2/2)||[A_X,[A_X,Phi_L(X)]]||.                      (AX.17)
```

The terms odd in `q`, including any equilibrium-current expectation, cancel
between the two twists.  This uses only Hermiticity.  It does not assume a
real matrix, time reversal, zero current, Perron--Frobenius positivity, or
equal amplitudes.

Both trial energies are at least the centered-sector ground energy.  At
least one is no greater than their average.  Min--max and (AX.13)--(AX.17)
therefore give the exact dichotomy

```text
boxed:
either dim ker(H_L-E_0)>=2,
or Delta_L <= 2pi^2 D_2(L) L0 L2/L1.                    (AX.18)
```

Here `Delta_L` is the gap above the unique centered-sector ground.  The
bound includes arbitrary configuration-dependent diagonal terms, which
commute with `U` and contribute zero, as well as local simultaneous
multi-loop terms and complex Hermitian amplitudes.

## 6. Controlled quasi-local extension

Finite range is not essential.  Split a periodized quasi-local interaction
into supports with cyclic `x1` diameter below `L1/2` and the remaining
wrapping tail.  If the local part has a uniform second-twist density `D_2`
and

```text
T_L=sum_(X wrapping) 2||Phi_L(X)||,                       (AX.19)
```

then the same proof gives

```text
Delta_L <=2pi^2D_2 L0L2/L1+T_L                           (AX.20)
```

on the unique-ground branch.  An exponentially decaying interaction with a
uniform finite second moment has `T_L` exponentially small in `L_min` up to
the polynomial volume factor.  Along the `GL6AW` sequence

```text
(L0,L1,L2)=(m,2m^3,m),  m odd,                            (AX.21)
```

equation (AX.20) closes when `T_L ->0`:

```text
Delta_L <= pi^2D_2/m+T_L ->0.                            (AX.22)
```

The finite second moment is substantive.  A system-spanning hop of port-zero
occupation through distance `L1/2` acquires phase `pi`, not a small phase.
Algebraic tails whose second twist moment diverges, all-to-all terms, or
wrapping terms of nonvanishing total norm can evade (AX.22).

## 7. Higher-order stability actually earned

Fix a truncation order `K`.  Suppose the canonical locked effective
Hamiltonian through order `K` has a translation-covariant linked interaction
of finite range `R_K` and finite `D_(2,K)`.  For tori with

```text
K<2L_min,  2R_K<L1,                                     (AX.23)
```

section 3 supplies exact port conservation and section 5 applies.  Along
(AX.21), after `K` is fixed,

```text
either centered-sector ground degeneracy,
or Delta_L^(<=K) <= pi^2 D_(2,K)/m ->0.                  (AX.24)
```

Thus no finite fixed order of controlled local corrections can turn the
`GL6AW` centered-sector model into a unique ground separated by a uniform
positive gap along this anisotropic Følner family.  This is a genuine
order-by-order stability theorem.

It is not a uniform finite-coupling theorem.  The order of limits is fixed:
first `K`, then volume.  Interchanging `K -> infinity` with the thermodynamic
limit requires a volume-uniform convergent quasi-local block diagonalization,
control of the wrapping tail, and a dressed locked-sector identification.
None has been derived here.

## 8. Remaining physical gates

GL6AX removes one live ambiguity: the anisotropic obstruction is stable under
all controlled fixed-order local locked corrections and is not an artifact
of retaining only the elementary order-six ring.

The following remain open:

1. a volume-uniform finite-`h/U_d` construction of the dressed locked
   dynamics, with a finite second twist moment and vanishing wrapping tail;
2. proof that the actual or selected physical ground lies in the centered
   port sector used by the character;
3. a compatible selected-state/GNS bridge showing that the finite twisted
   vectors remain outside the full zero-energy projection;
4. an isotropic small-character mode, its local dispersion and one common
   physical cone;
5. record authentication of that selected collective mode and its
   constitutive coupling to complete physical stress.

No graviton, gauge photon, Ricci ansatz, Einstein equation, gravity law, or
numerical `G` is inserted or derived by this packet.

## 9. Final theorem

> **GL6AX all-fixed-order twist-stability theorem.**  Affine independence of
> the four native `A3` port displacements makes every finite contractible
> locked-to-locked change conserve every port total.  The only finite-torus
> exception is a winding change, whose exact threshold is `2L_min` link
> flips.  Consequently every fixed perturbative order below that threshold
> conserves the centered port charge.  Any Hermitian translation-invariant
> finite-range locked Hamiltonian with this charge and bounded local second
> twist density obeys the degeneracy-or-`O(L0L2/L1)` gap dichotomy (AX.18),
> without a zero-current or time-reversal assumption.  Controlled quasi-local
> interactions obey (AX.20).  Hence the `GL6AW` anisotropic finite-size
> obstruction is stable at every fixed controlled local order, while finite
> coupling and selected-GNS promotion remain explicitly open.
