# GL6AR — LOCKED HEXAGON THERMODYNAMIC-SECTOR THEOREM

## Status and claim boundary

This lane starts only from sealed, independently audited `GL6AN` and
`GL6AO`.  It constructs the pure order-six hexagon-flip dynamics on the
degree-two locked configuration space, gives exact finite boundary and
periodic definitions, derives its native conserved cut fluxes, proves finite
Perron--Frobenius structure and fixed-boundary constrained ground-state
limit-point existence, and reduces selected finite-size component-gap closure
to one explicit variance gate.  It does not identify that closure with
infinite-volume GNS spectral gaplessness.

No gauge theory, photon, graviton, physical momentum, Lorentz cone, Ricci or
Einstein response, gravity, or `G` is assumed or derived.

The scalar terms in `GL6AO` are omitted because they shift every locked
configuration equally.  Put

```text
t := (63/8) h^6/U_d^5 > 0.
```

## 1. Native infinite incidence and locked configurations

Let cells be `x in Z^3`, ports be `a in {0,1,2,3}`, and

```text
d_0=(1,0,0), d_1=(0,1,0), d_2=(0,0,1), d_3=(0,0,0).
```

The active link `e=(x,a)` joins the parent node `P_x` to the child node
`C_{x+d_a}`.  An occupation configuration `n_e in {0,1}` is locked when

```text
sum_{e incident to v} n_e = 2                         (AR.1)
```

at every parent and child node.

For distinct ports `a,b,c`, one elementary hexagon based at `x` has the
ordered links

```text
(x,a),
(x+d_a-d_b,b),
(x+d_a-d_b,c),
(x-d_b+d_c,a),
(x-d_b+d_c,b),
(x,c).                                                (AR.2)
```

Each chosen port occurs twice, at opposite parity in the order.  A hexagon
is flippable precisely when its six occupations alternate.  Toggling all six
then preserves (AR.1).

## 2. Exact finite boundary and periodic models

### 2.1 Fixed exterior boundary

Let `eta` be any infinite locked configuration and `D` a finite link set.
Define

```text
Omega(D,eta)={n locked : n_e=eta_e for e not in D}.    (AR.3)
```

This is finite and nonempty.  For every elementary hexagon `c` whose six
core links lie in `D`, define

```text
T_c |n> = |n xor c>  if c is alternating in n,
          0          otherwise.                       (AR.4)
```

The exact fixed-boundary Hamiltonian is

```text
H_(D,eta) = -t sum_{c subset D} T_c                  (AR.5)
```

on `ell^2(Omega(D,eta))`.  The external links in `eta` determine the degree
test at boundary vertices.  A boundary-independent collared interior keeps
only cycles whose full 18-link projected support lies in `D`.

The sealed GL6AN target gives an exact one-hexagon example: fixing every
other link leaves exactly two alternating locked configurations, and (AR.5)
is `-t` times the adjacency matrix of `K_2`.

### 2.2 Periodic quotient

For `L>=4`, reduce the cell coordinates modulo `L`.  Let `Omega_L` be all
periodic locked configurations and let `Hex_L` contain the four port triples
at every cell:

```text
|E_L|=4L^3,
|V_L|=2L^3,
|Hex_L|=4L^3.                                         (AR.6)
```

The pure periodic Hamiltonian is

```text
H_L = -t sum_{c in Hex_L} T_c.                       (AR.7)
```

## 3. Quasi-local interaction

In the ordinary quasi-local spin algebra put

```text
p_v = 1_{k_v=2},
P_c = product_{v in V(c)} p_v,
tau_c = P_c (product_{e in c} X_e) P_c,
Phi(c) = -t tau_c.                                    (AR.8)
```

`tau_c` is bounded, self-adjoint, and supported on exactly 18 links: the six
cycle links and the other links incident to its six vertices.  Each link
belongs to 18 such projected supports.  Hence

```text
sup_e sum_{c:e in supp(tau_c)} ||Phi(c)|| <= 18t.      (AR.9)
```

The interaction is uniformly finite range and defines a unique quasi-local
dynamics.  Every local degree projector commutes with every term, so the
degree-two sector is invariant.  Equation (AR.8) is a local projected
interaction; no global infinite locked projector is inserted.

## 4. Native divergence and conserved topological sectors

Define the signed link variable

```text
s_(x,a)=2n_(x,a)-1 in {-1,+1}.                        (AR.10)
```

The parent and child forms of (AR.1) are

```text
sum_a s_(x,a)=0,
sum_a s_(x-d_a,a)=0.                                  (AR.11)
```

Subtracting them gives the exact native cell-lattice continuity identity

```text
sum_{i=0}^2 [s_(x,i)-s_(x-e_i,i)] = 0.                (AR.12)
```

Thus the signed occupation is a divergence-free integer link assignment.
This statement is derived from the lock; no external redundancy principle
is assumed.

On `Q_L`, let

```text
N_a(n)=sum_x n_(x,a),
S_a(n)=2N_a(n)-L^3.                                   (AR.13)
```

Summing (AR.12) over a slab proves that

```text
Phi_i(n)=sum_{x:x_i=r} s_(x,i) = S_i(n)/L             (AR.14)
```

is independent of the cut `r`, for `i=0,1,2`.  The dependent fourth
normalized signed port-count value is

```text
Kappa_3:=S_3/L=-Phi_0-Phi_1-Phi_2,                    (AR.15)
```

because `sum_a S_a=0`.  Thus `Phi_0,Phi_1,Phi_2` are three independent
coordinate cut fluxes, while `Kappa_3` is a dependent conserved port-count
value.  No fourth geometric cut has been derived.

Every hexagon contains two occurrences of each of three ports.  Alternation
makes one occurrence leave and the other enter, so every `N_a`, hence every
coordinate `Phi_i` and `Kappa_3`, is exactly conserved by every flip.

> **Sector theorem.** Every connected component of the finite flip graph is
> contained in one fixed tuple `(Phi_0,Phi_1,Phi_2;Kappa_3)`, where the first
> three entries are coordinate cut fluxes and the fourth is fixed by
> (AR.15).

No completeness claim is made: the conserved flux can label components
without proving that all configurations with the same flux are connected.

There are six exact frozen periodic configurations, obtained by occupying
the same two ports at every cell.  Every elementary hexagon uses each of its
ports twice with equal occupation, so none is alternating.  These are
isolated flip-graph vertices at extreme flux.  Active components also exist,
as witnessed by the sealed period-four alternating background.

## 5. Exact finite connected-component Hamiltonian

Let `G_C` be any connected component of a finite boundary or periodic flip
graph, with adjacency matrix `A_C`.  Then

```text
H_C=-t A_C.                                            (AR.16)
```

If `C` has more than one configuration, Perron--Frobenius gives a simple
largest adjacency eigenvalue `rho_C` and a unique normalized vector
`psi_C(n)>0`.  Therefore

```text
E_0(C)=-t rho_C                                        (AR.17)
```

is the unique component ground energy/vector.  If
`lambda_2(C)` is the second adjacency eigenvalue, the exact component gap is

```text
Delta_C=t[rho_C-lambda_2(C)]>0.                        (AR.18)
```

For a one-vertex frozen component there is no internal excitation.  Across
the full finite Hilbert space, ground components are exactly those with
maximal `rho_C`; degeneracy is the number of maximizing components.

Writing `d(n)` for the number of flippable hexagons and `d_max` and `d_bar`
for its maximum and component average gives the elementary energy bounds

```text
d_bar <= rho_C <= d_max,
sqrt(d_max) <= rho_C,
-t d_max <= E_0(C) <= -t max(d_bar,sqrt(d_max)).       (AR.19)
```

The square-root inequality follows from the principal star at a
maximum-degree vertex.

## 6. Exact ground-state transform and variational gap formula

Put `pi_C(n)=psi_C(n)^2`.  For every real function `f` on the component with
`sum_n pi_C(n)f(n)=0`, direct use of the PF eigenvalue equation gives

```text
<f psi_C,(H_C-E_0)f psi_C>
 = t sum_{{n,n'} in Edge(G_C)}
       psi_C(n)psi_C(n')[f(n)-f(n')]^2.               (AR.20)
```

Consequently

```text
Delta_C = inf_{f not constant}
  t sum_{{n,n'}} psi_C(n)psi_C(n')[f(n)-f(n')]^2
  / Var_pi(f).                                         (AR.21)
```

Equations (AR.20)--(AR.21) are exact.  They identify the missing
finite-size spectral question without assuming a continuum field.

## 7. Fixed-boundary exhaustion ground-state limits exist

Fix any infinite locked exterior `eta` and take an increasing locally
complete exhaustion `D_R`.  Choose a ground vector of (AR.5) in the complete
fixed-exterior configuration set, extend its state by `eta` outside `D_R`,
and take a weak-* convergent subnet.  Compactness of the state space supplies
such a subnet.

Every local degree projector has expectation one at every stage and in the
limit.  For any finite degree-preserving variation, the finite-volume ground
inequality applies once its support and interaction collar lie inside
`D_R`.  Passing to the limit proves the local ground inequality.

> **Fixed-boundary exhaustion theorem.** For each infinite locked boundary
> condition `eta`, each increasing locally complete finite-link exhaustion
> `D_R`, and each choice of a finite-volume ground state along that
> exhaustion, there is a weak-* convergent subnet.  Every such limit point is
> an infinite-volume constrained locked ground state for the pure hexagon
> interaction.

This is a ground state of the constrained degree-preserving dynamics.  It is
not a selected physical vacuum and is not promoted to a state of gravity.
The theorem does not prove retention of a nonlocal asymptotic sector in the
limit or that different boundary conditions yield distinct limit states.

## 8. A finite periodic active-energy bound is rigorous

The period-four GL6AN background lifts to every `Q_L` with `4|L`.  Choose one
translated sealed target hexagon per `4 x 4 x 4` block.  Their vertex collars
are disjoint, so all

```text
N_L=(L/4)^3=L^3/64                                    (AR.22)
```

toggles remain independently available.  The equal positive superposition
over their `2^{N_L}` subsets has expectation at most `-tN_L`; all additional
hexagon terms have nonnegative expectation before the overall minus sign.
Since each `||T_c||<=1`,

```text
-4tL^3 <= E_0(Q_L) <= -(t/64)L^3.                     (AR.23)
```

Thus, for every finite periodic quotient with `4|L`, the ground lies in an
active component and obeys

```text
|E_0(Q_L)|/L^3 >= t/64.                               (AR.24)
```

The frozen zero-energy vertices are not global periodic ground states on
those quotients.  This is a finite periodic subsequence bound; convergence
of the energy density or existence of a thermodynamic energy-density limit
is not proved.

## 9. Conditional finite-size component-gap theorem

Choose fixed real port weights `w_a`, and on `Q_L` let

```text
f_L(r)=min(r,L-r)/L,
F_L(n)=sum_{x,a} f_L(x_0) w_a n_(x,a).                 (AR.25)
```

`f_L` is a purely quotient-defined tent function with nearest-step
Lipschitz constant `1/L`; it is not called physical position or momentum.

In a hexagon, the two occurrences of each port carry opposite flip signs,
and their parent cells differ by at most one cyclic step in coordinate zero.
Therefore

```text
|F_L(n xor c)-F_L(n)| <= 3||w||_infinity/L.            (AR.26)
```

Let `psi_L` be the PF ground state of a chosen nontrivial component, center
`F_L` in `psi_L^2`, and put `V_L=Var_(psi_L^2)(F_L)`.  If `V_L>0`, applying
(AR.21), then
`2psi(n)psi(n')<=psi(n)^2+psi(n')^2`, gives the exact upper bound

```text
Delta_L <= 18t ||w||_infinity^2 L / V_L.               (AR.27)
```

Hence:

> **Conditional finite-size component-gap theorem.** If a sequence of
> selected ground components has `V_L >= v L^3` for one fixed nonzero `w`
> and some `v>0`, then
>
> ```text
> Delta_L <= [18t||w||_infinity^2/v] L^{-2} -> 0.       (AR.28)
> ```

This is a rigorous route to closure of selected finite-volume component
gaps using only a native locked occupation fluctuation.  Finite-size
component-gap closure is not, by itself, a theorem that the selected
infinite-volume GNS representation is gapless: a quasi-degenerate tower can
close while each selected pure thermodynamic phase remains spectrally
gapped.

### 9.1 The additional infinite-volume spectral bridge

To promote (AR.28) to GNS spectral gaplessness, one must additionally choose
compatible embeddings and boundary/sector data such that the selected PF
ground states `omega_L` converge weak-* to one constrained locked ground state
`omega`.  Let `(pi_omega,H_omega,Omega_omega)` be its positive-energy GNS
representation and `P_0=1_{\{0\}}(H_omega)` the projection onto the full
zero-energy subspace.  The finite-volume trial must be transported into
local observables, or quasi-local observables in the ground-state
energy-form domain, `A_R` such that

```text
xi_R := (1-P_0) pi_omega(A_R) Omega_omega,
||xi_R||=1,
mathcal E_omega(A_R)=<xi_R,H_omega xi_R> -> 0.         (AR.29)
```

For local `A`, this energy form is

```text
mathcal E_omega(A)
  := sum_c omega(A^* [Phi(c),A]),                      (AR.30)
```

where only finitely many finite-range interaction terms contribute.  A
strictly positive GNS spectral gap would lower-bound (AR.30) by that gap for
every normalized vector outside the full zero-energy subspace, contradicting
(AR.29).  For quasi-local `A`, (AR.30) means the closure of this local
quadratic form.  Thus
(AR.29), or an equivalent limiting spectral-weight construction, is the
needed bridge.  Neither compatible state/excitation transport nor the
vanishing boundary and sector errors required for it is proved here.
Merely requiring `omega(A_R)=0` would be insufficient when the zero-energy
ground space is degenerate.

## 10. Sharp obstruction to unconditional closure

Neither `GL6AN` nor `GL6AO` supplies:

1. selection of the component/flux sector that controls a thermodynamic
   ground-state limit;
2. an extensive lower bound on `V_L` in that selected ground state;
3. a conductance or correlation theorem that could replace that variance
   bound;
4. the compatible infinite-volume state/excitation bridge (AR.29).

The local lock and topology alone cannot supply item 2: exact frozen
components have zero variance for every such mode.  Although (AR.23) proves
the global periodic ground is active, it does not prove its low-character
variance is extensive.

Therefore unconditional selected-component finite-size gap closure is **not
yet proved or refuted**.  The finite-size result has been reduced to the
precise native alternative:

```text
extensive slow port-density variance
  -> selected component gap closes at least as L^{-2};
otherwise a new suppression/sector theorem is required.        (AR.31)
```

Even if that variance hypothesis is established, infinite-volume GNS
spectral gaplessness remains a separate conclusion until (AR.29), or an
equivalent spectral argument, is supplied.

Calling the conditional mode a photon or graviton would add exactly the
physics not yet derived.

## 11. Exact ceiling and next gate

GL6AR proves a well-defined thermodynamic locked hexagon model, three native
coordinate cut fluxes plus one dependent port-count invariant, finite PF
structure, fixed-boundary constrained ground-state limit points, a finite
periodic active-energy bound along `4|L`, and a sharp conditional gap bound.
It does not
prove a unique vacuum, same-sector connectivity, all-orders stability, an
unconditional finite-size component-gap closure, infinite-volume GNS
gaplessness, physical momentum, a propagation cone, stress coupling,
Ricci/Einstein response, gravity, or `G`.
