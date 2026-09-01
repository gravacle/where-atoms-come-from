# Distinct hostile audit — GL6AR locked hexagon thermodynamic sector

**Target:** `LANE_CROSS_RFT_GRA_GL6AR_LOCKED_HEXAGON_THERMODYNAMIC_SECTOR_V001/`  
**Frozen theorem SHA-256:** `960af10c683496ac921c3371a8182db7f018cf788f38dc486ecb0af95c089555`  
**Frozen author-manifest SHA-256:** `cdb1acc430b9c928807d8daabbf85848c40aeb1469da9a8fa4ec052ce34042cc`  
**Frozen author-seal-file SHA-256:** `42860871df4d86af4ce554b35460f928aade1f470de27fd617ebe9fc3dd15466`  
**Disposition:** `PASS__FINITE_AND_INFINITE_LOCKED_DYNAMICS_WELL_DEFINED__FLUX_AND_COMPONENT_SCOPE_EXACT__PF_AND_DIRICHLET_IDENTITIES_EXACT__FINITE_PERIODIC_ACTIVE_ENERGY_BOUND_ONLY__DELTA_LE_18_T_L_OVER_VARIANCE_EXACT_FOR_POSITIVE_VARIANCE__FINITE_SIZE_CLOSURE_CONDITIONAL__FULL_GNS_BRIDGE_STILL_MISSING__NO_MOMENTUM_PHOTON_CONE_GRAVITY_OR_G`

## 1. Custody and independence

The author packet contains exactly eleven files.  All eleven frozen bytes are
pinned in `AUDITED_TARGETS.sha256`; its nine-row manifest and one-row seal
verify exactly.  The twelve transitive dependency rows contain only the
sealed `GL6AN` and `GL6AO` author theorems/manifests/seals and their distinct
hostile-audit theorems/manifests/seals.  No later mutable lane is admitted.

The independent replay imports no author code.  On the frozen bytes, the
author physics replay passes `72161/72161` in normal and optimized modes, the
author packet passes `145/145`, the upstream `GL6AO` author/audit packets pass
`82/82` and `156/156`, and the upstream `GL6AN` author/audit packets pass
`79/79` and `58/58`.

## 2. Incidence, cycles, and the infinite interaction

For cells `x in Z^3`, the link `(x,a)` joins `P_x` to `C_(x+d_a)`, with

```text
d_0=e_0, d_1=e_1, d_2=e_2, d_3=0.
```

Direct traversal of the declared `Q_4` reconstructs `64` cells, `128`
constraint nodes, `256` links, and `256` distinct elementary hexagons.  Each
hexagon has six distinct links and vertices, uses three ports twice at
opposite parity, and each link belongs to six hexagon cores.  No quotient
degeneracy occurs for `L>=4`.

For a hexagon `c`,

```text
tau_c=P_c (product_(e in c) X_e) P_c
```

has precisely the six core links plus the twelve other links incident to its
six vertices.  The replay obtains an 18-link support for every `Q_4` term and
support multiplicity eighteen for every link.  Analytically, a vertex lies
on twelve elementary hexagons, while a link's two endpoints share the six
hexagons containing that link, giving `12+12-6=18`.  Therefore

```text
sup_e sum_(c:e in supp tau_c) ||-t tau_c|| <= 18t.
```

The interaction is uniformly finite range.  Because a degree projector on a
cycle vertex is a factor of `P_c`, and every other overlapping degree
projector meets only diagonal collar factors, every local degree projector
commutes with every `tau_c`.  Thus the ordinary quasi-local dynamics exists
and preserves the locked face without inserting an infinite global locked
projector.

## 3. Native invariants and finite components

Writing `s_(x,a)=2n_(x,a)-1`, the parent and child degree-two constraints give

```text
sum_a s_(x,a)=0,
sum_a s_(x-d_a,a)=0.
```

Their difference cancels the zero-displacement port and yields

```text
sum_(i=0)^2 [s_(x,i)-s_(x-e_i,i)]=0.
```

Slab summation proves three plane-independent integer coordinate cut fluxes.
If `S_a=2N_a-L^3`, then `Phi_i=S_i/L` for `i=0,1,2`, while the conserved
fourth normalized port count is only the dependent value
`Phi_3=-Phi_0-Phi_1-Phi_2`; it is not asserted to be a fourth geometric cut.
Each alternating flip changes one occurrence of each participating port up
and the opposite-parity occurrence down, so all four `N_a` are conserved.
Consequently components refine the invariant tuple; equality of flux is not
claimed to imply connectivity.

The six uniform two-port configurations are locked isolated vertices.  An
independent bipartite-capacity construction finds a period-four locked
configuration with an alternating target hexagon and confirms that toggling
it preserves every degree and port count.  With all exterior links fixed,
the target collar has exactly the two alternating configurations and flip
graph `K_2`.

For any connected finite component `C`, the restricted Hamiltonian is
`H_C=-tA_C`.  Irreducibility of the nonnegative adjacency matrix gives a
simple Perron root `rho_C` and a unique strictly positive component ground
vector for every nontrivial component.  If adjacency eigenvalues are ordered
decreasingly, then

```text
E_0(C)=-t rho_C,
Delta_C=t(rho_C-lambda_2(C)).
```

The bounds `d_bar<=rho_C<=d_max` follow from the constant-vector Rayleigh
quotient and the maximum row sum; `sqrt(d_max)<=rho_C` follows from the
principal star at a maximum-degree vertex.  Across the direct sum of
components, global finite ground components are precisely the maximizers of
`rho_C`.

## 4. Exact ground-state transform and the constant 18

Let `psi>0` be the normalized component PF vector and `pi=psi^2`.  For real
`f` centered in `pi`, the PF eigenvalue equation gives, with unordered graph
edges,

```text
<f psi,(H_C-E_0)f psi>
 = t sum_{{n,n'} in Edge(C)} psi(n)psi(n')[f(n)-f(n')]^2.
```

Since multiplication by `psi` is invertible on a finite component, minimizing
this quotient over nonconstant `f` is exactly the component gap, not merely
an upper bound.

For the quotient tent observable, the two occurrences of a fixed port in a
hexagon have opposite flip sign and parent coordinates separated by at most
one cyclic step.  Summing over the three participating ports gives

```text
|Delta_c F_L| <= 3||w||_infinity/L.
```

There are `4L^3` possible hexagons.  The elementary inequality
`2psi(n)psi(n')<=pi(n)+pi(n')` therefore gives

```text
sum_{{n,n'}} psi(n)psi(n')
 <= (1/2) sum_n pi(n)d(n)
 <= 2L^3.
```

For a nontrivial selected component with `V_L=Var_pi(F_L)>0`, the exact
Rayleigh estimate is consequently

```text
Delta_L
 <= t (9||w||_infinity^2/L^2)(2L^3)/V_L
 = 18t||w||_infinity^2 L/V_L.
```

Thus `V_L>=vL^3` implies only the conditional finite-component statement
`Delta_L<=(18t||w||^2/v)L^-2`.  The repaired theorem correctly excludes
zero-variance and one-vertex components from the displayed division.

## 5. Ground-state limits and active-energy scope

For a locked exterior `eta` and locally complete exhaustion `D_R`, every
finite fixed-boundary Hamiltonian has a ground state.  Extending it by `eta`
gives states on the quasi-local algebra, whose weak-* compact state space has
a convergent subnet.  Degree projectors retain expectation one.  Once the
support and interaction collar of a finite degree-preserving variation lie
inside `D_R`, the finite ground inequality applies and passes to the limit.
Every such limit point is therefore a constrained locked ground state.

This argument does not retain or distinguish a nonlocal asymptotic boundary
sector.  The repaired packet says exactly that and makes no uniqueness or
physical-vacuum claim.

For every `4|L`, disjoint translates of the sealed target produce
`N_L=L^3/64` independent toggles.  Their equal positive hypercube
superposition has expectation at most `-tN_L`; every additional flip matrix
has nonnegative expectation before the overall minus sign.  Together with
the interaction-norm lower bound,

```text
-4tL^3 <= E_0(Q_L) <= -(t/64)L^3.
```

This is a finite periodic subsequence bound.  It neither proves convergence
of `E_0(Q_L)/L^3` nor establishes a thermodynamic energy-density limit.

## 6. Finite-size closure is not a GNS gap theorem

The original hostile review rejected identification of conditional
finite-size component-gap closure with infinite-volume spectral
gaplessness.  The frozen repaired packet now keeps those statements
separate.  A quasi-degenerate finite tower can close while a selected pure
thermodynamic phase stays gapped.

Even compatible weak-* convergence is insufficient by itself.  If
`H_omega` is the positive GNS Hamiltonian and
`P_0=1_{\{0\}}(H_omega)` projects onto its full zero-energy subspace, a valid
spectral bridge needs local, or energy-form-domain quasi-local, observables
whose vectors

```text
xi_R=(1-P_0)pi_omega(A_R)Omega_omega
```

are normalized and satisfy `<xi_R,H_omega xi_R> -> 0`.  Centering only
against `Omega_omega` would not remove other zero-energy vectors when the
kernel is degenerate.  The repaired theorem states the full-`P_0` condition
and explicitly leaves state/excitation transport and boundary/sector-error
control unproved.

Neither the sealed inputs nor GL6AR supplies selection of the maximizing
component/flux sector, an extensive variance lower bound in it, or the GNS
bridge.  Therefore unconditional finite-size component-gap closure and
infinite-volume GNS gaplessness both remain open.

## 7. Promotion attacks and verdict

The quotient tent is not physical position or momentum; the integer cut
invariants are not promoted to a gauge field; PF positivity is not physical
state selection; and the conditional trial is not a photon or graviton.
There is no propagation cone, stress coupling, Ricci/Einstein response,
gravity, or `G`.

**Hostile verdict: PASS.**
