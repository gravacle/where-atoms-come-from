# GL6AU — `v/g=0` static-structure phase-closure theorem

## Status and claim class

Author packet status: **author frozen and sealed**.  Promotion requires a
distinct independent hostile audit of these exact bytes.

This lane starts from sealed and independently audited `GL6AO`, `GL6AR`,
`GL6AS`, and the exact `GL6AT` operator/evidence crosswalk.  It derives a
sharper exact finite-size soft-mode
criterion for the pure order-six model, identifies the one missing static
inequality, and determines what the admitted zero-potential quantum-ice
evidence would imply if used as a controlled phase hypothesis.

Three claim classes are kept disjoint:

1. **proved here:** exact finite-component oscillator/structure bound and
   exact shortcut obstructions;
2. **scientific evidence:** published finite-size QMC and its fitted Gaussian
   comparison at the exactly crosswalked `v/g=0` point;
3. **not adopted as theorem:** a stated static-structure phase hypothesis and
   the separate infinite-volume spectral bridge.

Nothing here proves the all-orders finite-`h/U_d` phase, identifies a
translation character with physical momentum, or derives a physical cone,
stress law, Ricci/Einstein response, gravity, or `G`.

## 1. Exact comparison model

After removal of the common scalar, sealed `GL6AO` gives

```text
H_0=-J sum_c T_c,
J=(63/8)h^6/U_d^5>0.                                  (AU.1)
```

`GL6AT` proves that (AU.1), on the degree-two diamond-link Hilbert space, is
exactly the fully-packed-loop quantum-ice comparison Hamiltonian at

```text
v/g=0,                                                 (AU.2)
```

not the Rokhsar--Kivelson point `v/g=1`.  Equations (AU.1)--(AU.2) concern
the displayed order-six interaction.  Unknown order-eight and longer-loop
terms are outside this lane.

## 2. Exact first-character finite-component bound

Take a standard periodic quotient `Q_L`, `L>=4`, with `N=L^3` parent cells.
Let `C_L` be a connected flip component preserved by all three unit cell
translations.  Its Hamiltonian is

```text
H_C=-J A_C.                                             (AU.3)
```

The Perron--Frobenius ground vector `psi_L>0` is unique inside the component.
Every cell translation therefore fixes it with eigenvalue one.

Put

```text
q_L=2pi/L,
z=(exp(iq_L),1,1,1).                                    (AU.4)
```

Choose any real normalized port vector

```text
u=(0,u_1,u_2,u_3),
u_1+u_2+u_3=0,
||u||=1.                                                (AU.5)
```

This is the exact two-dimensional plane transverse to both rows of the
incidence symbol at (AU.4).  Define the normalized density and its equal-time
weight

```text
rho_u(q_L)=N^(-1/2) sum_(x,a)
  exp(iq_L x_0) u_a [n_(x,a)-1/2],
S_(u,L)(q_L)=<psi_L|rho_u(q_L)^* rho_u(q_L)|psi_L>.
                                                               (AU.6)
```

The nontrivial translation character makes `rho_u(q_L)psi_L` orthogonal to
the unique component ground vector.  Since the adjacency Hamiltonian and
`psi_L` are real, the positive spectral weights at `q_L` and `-q_L` agree.
Consequently (AU.6) is entirely strictly positive-frequency weight whenever
it is nonzero, and the double commutator below is its first spectral moment.
For the four elementary hexagon orientations, the exact `GL6AS` cycle
columns satisfy

```text
||C(z)^*u||^2=3|1-exp(iq_L)|^2
             =12 sin^2(pi/L).                           (AU.7)
```

Let `t_d=<psi_L|T_(x,d)|psi_L>` for the orientation whose missing port is
`d`.  Full cell-translation invariance makes it independent of `x`.  The
nonnegative PF entries
of `psi_L` and the partial-flip form of `T_(x,d)` give

```text
0<=t_d<=1.                                              (AU.8)
```

The exact double commutator is therefore

```text
f_(u,L)(q_L)
 =(J/2) sum_d t_d |u^*C_d(z)|^2
 <=6J sin^2(pi/L).                                      (AU.9)
```

This estimate also has an exact Perron--Frobenius/Dirichlet reading.  If
`rho_C` is the top eigenvalue of `A_C`, define on the component

```text
P_C(n,m)=A_C(n,m) psi_L(m)/[rho_C psi_L(n)],
pi_C(n)=psi_L(n)^2.                                     (AU.9a)
```

Then `P_C` is a reversible Markov kernel with stationary law `pi_C`, is
similar to `A_C/rho_C`, and

```text
Delta_C=J rho_C gap(P_C).                               (AU.9b)
```

The numerator in (AU.9) is precisely `J rho_C` times the Dirichlet form of
the diagonal trial function in (AU.6), while `S_(u,L)` is its variance under
`pi_C`.  Thus the unresolved phase datum is not a missing operator: it is
control of the exact, generally nonuniform PF stationary measure.

If `S_(u,L)(q_L)>0`, the positive-frequency single-mode inequality inside
the same component yields

```text
boxed:
Delta_C(L)
 <= 6J sin^2(pi/L)/S_(u,L)(q_L).                       (AU.10)
```

No coherent phase, continuum field, pole, or physical momentum is assumed
in (AU.10).

## 3. Static-exponent closure theorem

Suppose a selected sequence of translation-stable components and one fixed
choice of transverse polarization obeys, for constants `s>0`, `L_0`, and
`alpha<2`,

```text
S_(u,L)(q_L) >= s L^(-alpha),
L>=L_0.                                                 (AU.11)
```

Then (AU.10) and `sin(pi/L)<=pi/L` give

```text
Delta_C(L)
 <= (6pi^2J/s)L^(alpha-2) -> 0.                        (AU.12)
```

> **Static-exponent closure theorem.** Any first-character transverse
> structure exponent strictly below two forces selected finite-component gap
> closure.  A Coulomb/linear-mode equal-time law `S(q_L)>=s/L` needs only
> `alpha=1` and gives `Delta_C(L)=O(J/L)`.

This sharpens the sufficient condition isolated in `GL6AR`.  Define the two
real quadratures

```text
F_c=sum_(x,a) cos(q_L x_0)u_a[n_(x,a)-1/2],
F_s=sum_(x,a) sin(q_L x_0)u_a[n_(x,a)-1/2].             (AU.13)
```

Translation invariance and `2q_L != 0` for `L>=4` give the exact identity

```text
Var(F_c)+Var(F_s)=N S_(u,L)(q_L).                       (AU.14)
```

Thus (AU.11) implies that at least one real slowly varying observable has

```text
Var(F_c) or Var(F_s) >= (s/2)L^(3-alpha).               (AU.15)
```

For the expected `alpha=1` law this is an `L^2` variance, not the stronger
`L^3` variance used in the first `GL6AR` corollary.  Nevertheless it is
already sufficient for the `O(1/L)` closure (AU.12).  Extensive variance was
sufficient, not necessary.

## 4. What first principles still do not supply

The sealed exact results do not prove (AU.11).  They give:

- the degree lock and transverse rank-two cycle image;
- exact port conservation and oscillator strength `O(q^2)`;
- an active finite-periodic energy bound;
- finite-component PF positivity and a variational gap formula.

None lower-bounds the low-character static weight.  Frozen lawful components
have zero variance, and the sealed inequalities do not exclude
`S(q)=O(q^2)` with a nonclosing single-mode quotient.  The exact obstruction is

```text
f(q)=O(q^2) does not imply Delta(q)->0
without S^+(q) decaying more slowly than q^2.            (AU.16)
```

Nor do the sealed inputs select the globally controlling component/flux
sector or prove that it is translation stable.

## 5. Exact shortcut disposition

### 5.1 RK/Markov comparison does not transfer the phase

For a flip graph with adjacency `A` and flippability degree matrix `D`,

```text
H_RK=J(D-A),
H_0=-JA=H_RK-JD.                                       (AU.17)
```

At the RK point the equal-amplitude state is exact within a connected
component.  At zero potential the ground vector is instead the PF vector of
`A`.  In the actual locked space `D` is not constant.  Exact frozen vertices
have degree zero.  For `4|L`, the lifted active construction used in `GL6AR`
contains at least `L^3/64` simultaneously flippable disjoint hexagons, so its
degree is at least `L^3/64`.  Consequently, on the full periodic locked
space,

```text
inf_c ||D-cI|| >= L^3/128.                              (AU.18)
```

Thus `JD` is an extensive, configuration-dependent operator, not a bounded
small perturbation.  This does not rule out a special regular component; it
does show that no uniform Radon--Nikodym, Dirichlet-form, or spectral-gap
comparison between the RK and `v/g=0` ground measures follows from the sealed
data.

### 5.2 Flux-sector towers are not a GNS proof

The pure Hamiltonian preserves exact flux/port sectors.  A vanishing
finite-size energy difference between different sectors would prove only a
global sector tower.  As the repaired `GL6AR` audit emphasizes, such a tower
can coexist with a positive excitation gap inside each selected pure
infinite-volume phase.  A compatible local or energy-form-domain excitation
outside the full zero-energy GNS subspace is still required.

The conserved half-filled port numbers do not supply an isotropic shortcut
either.  In a zero-port sector `N_a=L^3/2`.  The elementary large twist of one
port along coordinate zero changes translation character by

```text
exp(2pi i N_a/L)=exp(pi i L^2)=1                       (AU.18a)
```

on the even cubic tori on which half filling is integral.  Rectangular tori
with odd transverse area can restore a nontrivial LSM character, but the
elementary twist energy scales as transverse area divided by twist length;
it does not give a closing bound in the isotropic three-dimensional limit.
Stronger quasi-adiabatic flux-insertion theorems distinguish a unique gapped
state from degeneracy/topological order; they do not by themselves construct
the local selected-GNS soft excitation required here.  No such theorem is
silently imported.

### 5.3 Rigorous Wilson/Villain compact-`U(1)` theorems do not map exactly

Fröhlich--Spencer-type massless-phase theorems concern continuous compact
link variables in a Euclidean Wilson/Villain weak-coupling action.  The
present model is a spin-`1/2` quantum-link Hilbert space with a strict
degree-two projection and only the discrete hexagon kinetic operator at the
compared order.  No exact transfer-matrix, measure, or coupling map from
(AU.1) to the hypotheses of those theorems has been derived.  Their
conclusion therefore cannot be imported as a proof of (AU.11).

## 6. Controlled scientific evidence at the exact comparison point

The frozen `GL6AT` screen identifies two mutually reinforcing results at
the exact `v/g=0` comparison point:

1. Shannon et al.'s zero-temperature GFMC and finite exact diagonalization
   find zero well inside the numerically inferred liquid region and find the
   flux-sector scaling expected for a three-dimensional Coulomb phase.
2. Benton et al.'s lattice Gaussian theory, fitted to equal-time QMC and
   finite-size ground energies, has two transverse linear character modes.
   Its microscopic link structure factor has low-character weight
   proportional to the mode frequency, hence the fitted law is
   `S_T(q) proportional to |q|` with nonzero fitted coefficient.

This is direct scientific support for the `alpha=1` premise of (AU.11), and
(AU.12) shows exactly what that evidence means for the microscopic
comparison Hamiltonian: its selected finite-component gap should close as
`O(J/L)`.

The evidence is not a certified lower bound valid for every `L`, a proof of
component connectivity/selection, or a direct microscopic real-frequency
pole measurement.  The GL6AT evidence classification has passed its distinct
hostile audit; no external paper byte is made part of this local seal.

## 7. Explicit working hypotheses, not hidden proofs

For later model-conditional work, define:

> **`ICE0-STATIC`.** Standard diamond tori at `v/g=0` admit a selected
> globally controlling, translation-stable component sequence and fixed
> transverse `u` for which `S_(u,L)(2pi/L)>=s/L` for all sufficiently large
> `L`, with `s>0`.

Under `ICE0-STATIC`, (AU.12) rigorously gives

```text
Delta_C(L) <= (6pi^2J/s)L^(-1).                        (AU.19)
```

For an infinite-volume conclusion one must add:

> **`ICE0-GNS-BRIDGE`.** The selected PF states have a compatible weak-*
> ground-state limit, and the finite-character trial transports to normalized
> local or energy-form-domain quasi-local excitation vectors outside the full
> zero-energy projection whose limiting energy form tends to zero.

`ICE0-STATIC + ICE0-GNS-BRIDGE` implies zero GNS spectral gap for that
selected pure-model representation.  Neither hypothesis is proved here.
The published comparison evidence strongly supports the first and gives
effective-theory support for the second; adopting them is a scientific phase
choice, not a mathematical deduction.

## 8. Final disposition

GL6AU closes the **algebra-to-test** part of the phase question.  The missing
quantity is no longer an unspecified extensive variance: it is the precise
first-character lower bound (AU.11), and any exponent below two suffices.
The exact `v/g=0` numerical/effective literature supports the expected
`alpha=1` law, so the direct comparison-model route to a soft collective
sector is scientifically strong and mathematically conditional on one
explicit static inequality.

No result here derives a physical photon, physical momentum or speed, a
common cone, stress conservation, Ricci/Einstein response, gravity, or `G`.
