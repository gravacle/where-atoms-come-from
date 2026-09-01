# GL6AY — finite-coupling prethermal locked bridge

## Status and claim class

**Status:** author frozen and sealed; a distinct independent hostile audit is
required before promotion.

**Claim class:** primary-source-backed, volume-uniform, small-but-finite
`h/U_d` prethermal normal-form theorem for the exact native F3 degree-lock
parent; exact port conservation and finite second twist moment for the
`N_def`-conserving effective interaction; exact separation of local charged
leakage from winding locked-to-locked mixing; sharp whole-band proof boundary.

**Not claimed:** a convergent all-orders Schrieffer--Wolff, Kato, or Feshbach
series; an exact invariant locked space for the full microscopic
Hamiltonian; a selected equilibrium phase or GNS state; a physical momentum,
isotropic cone, photon, graviton, Ricci tensor, Einstein equation, gravity,
or `G`.

## 1. Exact F3 parent and finite-volume setting

Use the frozen native `A3` incidence of `GL6AN`.  A coarse cell `x` contains
the four link qubits

```text
e=(x,a),  a=0,1,2,3,                                  (AY.1)
```

so its local Hilbert dimension is `16`.  The link joins parent `P_x` to
child `C_(x+d_a)`, with

```text
d_0=(1,0,0), d_1=(0,1,0), d_2=(0,0,1), d_3=(0,0,0).   (AY.2)
```

For every original incidence vertex put

```text
k_v=sum_(e incident v)n_e,  q_v=k_v-2,
N_def=sum_v q_v^2.                                     (AY.3)
```

On a finite periodic quotient `Lambda_L`, after removal of the common scalar,
the exact Hamiltonian is

```text
H_L=U_d N_def+W_L,
W_L=-h sum_e X_e,  U_d>0, h>0.                         (AY.4)
```

All summands `q_v^2` commute because they are diagonal in link occupation.
Since `k_v` is in `{0,1,2,3,4}`,

```text
spec(q_v^2)={0,1,4} subset Z,
exp(2 pi i N_def)=I.                                   (AY.5)
```

On each finite quotient, the locked projection is exactly

```text
P_L=chi(N_def=0),  Q_L=I-P_L.                          (AY.6)
```

Thus no independent locked projector or gauge constraint has been inserted.
There is no corresponding global spectral projection in the infinite-volume
quasi-local algebra; the local replacement is constructed in section 4.

## 2. Strong support and the exact local norm

The constraint summands in (AY.3) overlap, so the single-site form of the
closed-system prethermal theorem cannot be used literally.  Use instead the
finite-range extension pinned in `PRIMARY_SOURCES.md`.

A potential term `Z_S` is **strongly supported** on the connected cell set
`S` when it is supported on `S` and commutes with every constraint summand
`q_v^2` whose complete support is not contained in `S`.  For a potential
`Z=sum_S Z_S`, define

```text
||Z||^str_kappa
 =sup_x sum_(S connected, x in S) exp(kappa |S|)||Z_S||,
kappa>0.                                                (AY.7)
```

This is a potential norm, not the global operator norm.  Its constants are
uniform in `|Lambda_L|`.

The parent constraint at `P_x` is supported in cell `x`.  The child
constraint at `C_y` is supported on the four cells `y-d_b`.  Therefore a
flip `X_(x,a)` is strongly supported on

```text
S_(x,a)={x+d_a-d_b: b=0,1,2,3},                         (AY.8)
```

a connected set of at most four coarse cells.  A cell belongs to at most
sixteen such strong supports.  Hence every norm in (AY.7) is finite and
volume-independent.

Define the exact `N_def` average

```text
<Z>_N=(1/2pi) integral_0^(2pi)
      exp(i theta N_def)Z exp(-i theta N_def)dtheta,      (AY.9)
D_0=<W_L>_N,  V_0=W_L-D_0.                              (AY.10)
```

The strong-support extension makes (AY.9) support-preserving term by term.
For any fixed `kappa_0>0`, put

```text
d_0=||D_0||^str_(kappa_0),
v_0=||V_0||^str_(kappa_0).                              (AY.11)
```

Contractivity of the average and (AY.8) give the deliberately crude but
explicit sufficient bounds

```text
d_0 <=16 exp(4kappa_0)h,
v_0 <=32 exp(4kappa_0)h.                                (AY.12)
```

The theorem below uses the exact `d_0,v_0`, not the crude envelope.

## 3. Primary closed-system normal-form theorem

Define

```text
nu_0=(54pi/kappa_0^2)(d_0+2v_0),                       (AY.13)
n_* =floor{(U_d/nu_0)/[1+ln(U_d/nu_0)]^3}-2,            (AY.14)
kappa_*=kappa_0/[1+log(n_*+1)].                         (AY.15)
```

Assume the exact primary-source hypotheses

```text
U_d >=9pi v_0/kappa_0,
n_*>=1.                                                 (AY.16)
```

Then the closed-system theorem, with its strong-support extension, supplies
Hermitian potentials `D_hat_L,V_hat_L` and a quasi-local unitary `Y_L` such
that the following finite-volume identity is exact:

```text
Y_L H_L Y_L^*
 =U_d N_def+D_hat_L+V_hat_L,                             (AY.17)
[D_hat_L,N_def]=0,                                      (AY.18)
||V_hat_L||^str_(kappa_*)
 <=(2/3)^n_* v_0,                                      (AY.19)
||Y_L Z Y_L^*-Z||^str_(kappa_*)
 <=C(nu_0/U_d)||Z||^str_(kappa_0).                      (AY.20)
```

The numerical `C` is independent of every model parameter and, most
importantly, of the volume.  The potentials and automorphisms have a direct
infinite-volume quasi-local meaning even when the extensive operators do
not.

Equation (AY.17) is an exact decomposition **with the remainder retained**.
Equations (AY.18)--(AY.19) do not say that the full transformed Hamiltonian
commutes with `N_def`.

The average (AY.9) is not globally trivial.  If a link flip has
`delta=1-2n_e` and endpoints `u,v`, then

```text
Delta N_def=2delta(q_u+q_v)+2.                          (AY.21)
```

An absent-edge addition is resonant when `k_u+k_v=3`; an occupied-edge
removal is resonant when `k_u+k_v=5`.  Those charged-sector processes belong
to `D_0`.  On the locked space, `k_u=k_v=2`, every one-link flip costs
`2U_d`, and only the restricted statement

```text
P_L D_0 P_L=0.                                          (AY.22)
```

holds.

## 4. Exact port conservation of the effective locked interaction

Choose the potential representation in the primary construction termwise:

```text
D_hat_L=sum_S D_hat_L(S),
D_hat_L(S)=<H_hat_L(S)>_N.                              (AY.23)
```

Each term is strongly supported on `S` and commutes with `N_def`.  Do not
compress it with a nonexistent infinite-volume global projector.  Instead
define the contained-constraint collar

```text
N_S=sum_(v: supp(q_v^2) subset S)q_v^2,
P_S^0=chi(N_S=0)
     =product_(v: supp(q_v^2) subset S)chi(q_v^2=0),     (AY.24)
Phi(S)=P_S^0 D_hat(S)P_S^0.                             (AY.24a)
```

Strong support gives

```text
[D_hat(S),q_v^2]=0
 whenever supp(q_v^2) is not contained in S.             (AY.24b)
```

Termwise pinching gives `[D_hat(S),N_def]=0`; subtracting the outside
constraint sum therefore yields `[D_hat(S),N_S]=0`.  Hence `Phi(S)` is a
finite-support operator that agrees with `D_hat(S)` on every globally locked
input and maps that input to another globally locked local configuration:
contained constraints stay in the nonnegative zero eigenspace and every
noncontained constraint is individually unchanged.  The collar projection
preserves strong support and obeys
`||Phi(S)||<=||D_hat(S)||`, so it inherits the same exponential potential
bound.

On the infinite lattice, every finite `S` is contractible.  The sealed and
independently audited `GL6AX` affine theorem consequently gives the exact
locked-endpoint implication

```text
<n'|Phi(S)|n> !=0,  n,n' globally locked
 implies N_a(n')-N_a(n)=0,  a=0,1,2,3.                  (AY.25)
```

Here `N_a` is the native port-`a` occupation.  Equivalently, the interaction
`{Phi(S)}` has the four local port `U(1)` actions when represented on the
locked configuration algebra.  It is not a full-Hilbert-space symmetry of
each collar operator.  Thus

```text
boxed: {Phi(S)} is a well-typed finite-support interaction that
       agrees locally with D_hat on locked configurations and has
       exact termwise port U(1)^4 there.                  (AY.26)
```

No port symmetry has been asserted for the full microscopic Hamiltonian.
In finite volume the legitimate global compression satisfies

```text
P_L D_hat_L P_L
 =P_L[sum_S Phi_L(S)]P_L.                               (AY.26a)
```

On a torus, split (AY.24a) into supports with an injective lift and wrapping
supports.  Every injectively lifted term still obeys (AY.25).  The wrapping
tail is controlled by

```text
T_L=sum_(S wrapping)2||Phi_L(S)||.                       (AY.27)
```

For fixed positive `kappa_*`, connected exponential strong-support decay
makes `T_L` exponentially small in the wrapped period, up to its polynomial
volume count.

## 5. Finite second twist moment

For a nonwrapping support `S`, let `m=|S|` and use the local port-zero twist
generator `A_S` of `GL6AX`.  Bounded coordination gives a geometry constant
`C_geo`, independent of `S,L`, for which

```text
||A_S||<=C_geo m^2.                                     (AY.28)
```

Because `P_S^0` and `A_S` are both diagonal in link occupation,

```text
[A_S,[A_S,Phi_L(S)]]
 =P_S^0[A_S,[A_S,D_hat_L(S)]]P_S^0.                     (AY.29)
```

Projection does not increase norm, and therefore

```text
||[A_S,[A_S,Phi_L(S)]]||
 <=4 C_geo^2 m^4||D_hat_L(S)||.                         (AY.29a)
```

After summing per unit volume,

```text
D_2(L)
 <=4 C_geo^2 sup_(m>=1){m^4 exp(-kappa_*m)}
   ||D_hat_L||^str_(kappa_*)
 <infinity.                                             (AY.30)
```

The bound is uniform in volume at fixed finite `h/U_d` satisfying (AY.16).
Consequently the complete quasi-local `GL6AX` dichotomy applies to the
effective locked interaction:

```text
either the centered-sector effective ground is degenerate,
or Delta_L^eff
   <=2pi^2 D_2(L)L_0L_2/L_1+T_L.                        (AY.31)
```

This remains an effective-sector, anisotropic finite-size statement.  It
does not select the centered sector or an infinite-volume state.

## 6. Contact with the sealed order-six parent

The normal-form interaction has a formal expansion in `h/U_d`.  On each
finite-volume `P_L`, the orders below six are common scalars and the first
configuration-changing matrix element is therefore invariant under harmless
within-`P_L` normal-form gauge changes.  By (AY.26a), the sealed `GL6AO`
coefficient gives

```text
P_L[sum_S Phi_L(S)]P_L
 =C_eff P_L-(63/8)(h^6/U_d^5)sum_c T_c
  +higher normal-form terms.                            (AY.32)
```

Equation (AY.32) identifies the first non-scalar coefficient of the exact
normal form.  It is not a convergent-series assertion and does not replace
the nonperturbative remainder bound (AY.19).

## 7. Strictly conditional finite-observation-horizon corollary

Let `tau_t^H` denote the Heisenberg dynamics of the full microscopic
Hamiltonian.  For spatial dimension `d=3`, any

```text
0<r_1<ln(3/2)/4                                         (AY.33)
```

and any local observable `O`, the primary theorem supplies a constant
`K_3(O)` independent of `U_d` and the volume such that

```text
||tau_t^H(O)-tau_t^(U_d N_def+D_hat)(O)||
 <=K_3(O)/U_d,
t<=t_*:=exp(r_1 n_*).                                   (AY.34)
```

The equation is stated in the theorem's `hbar=1` clock convention.  If an
observation or record-retention protocol separately declares a finite
horizon `t_obs` in the same calibrated units and verifies

```text
t_obs<=t_*,                                             (AY.35)
```

then local F3 dynamics over that window shadows, within the displayed
local-observable error, the exactly `N_def`-conserving effective dynamics.
On a finite quotient, effective data initially in `P_L` remain there; their
local locked matrix elements are represented by `{Phi_L(S)}` and obey the
exact port conservation (AY.25)--(AY.26).

Neither the physical value of `U_d/h` nor conversion of (AY.34) to a physical
clock is supplied by this theorem.  Equations (AY.34)--(AY.35) do not infer a
record lifetime, gravity, or any gravitational constant.

## 8. Local charged excursions versus winding return

By definition of the pinching,

```text
P_L V_hat_L P_L=0.                                      (AY.36)
```

This finite-volume identity is not an invariant-space theorem.  A finite
strongly supported `V_hat_L(S)` may take locked data into `Q_L`, creating
local degree defects.  Such a local `P_L -> Q_L` excursion may change one or
more bare port totals and does **not** require noncontractible support.

Equation (AY.20) controls conjugation of declared local potentials, and the
primary theorem also supplies the corresponding fixed-local-observable
operator-norm control.  Neither statement bounds
`||Y_L^*P_LY_L-P_L||` uniformly in volume, and no global dressed spectral
subspace exists in the infinite quasi-local algebra.  No such projector or
subspace-closeness claim is made.

If a finite-volume sequence later returns to `P_L` and the union of its
support is finite and contractible, its two locked endpoints obey `GL6AX`
and hence

```text
Delta N_a=0 for every a.                                (AY.37)
```

On a torus, a returning `P_L -> P_L` process can change a port total only by
winding, and its symmetric difference contains at least `2L_min` links.
Thus the exact distinction is

```text
local P_L -> Q_L leakage:       no winding required;
port-changing P_L -> P_L return: winding required.       (AY.38)
```

The upper bound (AY.19) does not by itself prove a nonzero remainder; it
proves that no topological argument can remove local charged leakage from
the theorem's allowed remainder class.

## 9. Nonzero first winding coefficient on a fixed torus

The wrapping ceiling is realized perturbatively.  Start from the uniform
locked configuration containing ports `j` and a spectator port.  Toggle the
alternating `j/3` row around direction `j`.  Let `C` be its simple cycle of

```text
r=2L_j                                                  (AY.39)
```

changed links and let `|n'>` be the resulting locked configuration.  No
nonempty proper subset `S` of a simple cycle has locked endpoints at every
cycle vertex.  If `b(S)` is the number of boundary vertices of `S` along
`C`, then

```text
E(S)=U_d b(S),  2<=b(S)<=r.                              (AY.40)
```

The first possible canonical off-diagonal coefficient between `|n>` and
`|n'>` is order `r`.  Every order uses each changed link once.  There is no
intermediate return to `P_L`, so folded terms and lower-order endpoint terms
cannot contribute.  The direct Kato/Feshbach word sum is

```text
A_C^(r)
 =-h^r sum_(pi in S_r) product_(s=1)^(r-1)
   [1/(U_d b(S_(pi,s)))],                                (AY.41)
```

where `S_(pi,s)` is the set of the first `s` links in the order `pi`.  All
denominators are positive and all words have the same sign.  Therefore

```text
r! h^r/(rU_d)^(r-1)
 <=|A_C^(r)|
 <=r! h^r/(2U_d)^(r-1),                                 (AY.42)
A_C^(r)!=0.                                             (AY.43)
```

For each fixed torus, ordinary finite-dimensional analytic perturbation
therefore gives nonzero locked-to-locked port mixing for sufficiently small
nonzero `h`, with (AY.41) its leading coefficient.  Its analytic radius is
not uniform in volume.

## 10. Sharp whole-band proof boundary

The unperturbed gap from `P_L` to `Q_L` is `2U_d`, while the commuting flip sum
has global norm

```text
||W_L||=h|E_L|.                                         (AY.44)
```

A raw Feshbach expansion of

```text
Q_L(E-Q_L H_L Q_L)^(-1)Q_L                             (AY.45)
```

by a global Neumann series, or a global direct-rotation construction for the
entire descendant of `P_L`, therefore obtains only a sufficient condition of
the form

```text
h|E_L|<c U_d,                                           (AY.46)
```

whose radius collapses as the volume grows.  The first non-scalar locked
interaction is itself extensive at order `h^6/U_d^5`, so current inputs do
not supply a volume-uniform spectral contour isolating the entire split
locked band from all charged bands.

This is exactly the many-body boundary identified by the pinned primary
Schrieffer--Wolff source: fixed-order linked coefficients remain controlled,
but the global whole-band transformation and its Taylor series are not
licensed uniformly in volume.  The pinned prethermal source likewise stops
at the optimal `n_*` because local ranges and norms grow and resonant
connections can appear.

Equations (AY.44)--(AY.46) prove failure of the standard whole-band proof
route under the present premises.  They do not prove that every specially
integrable or otherwise algebraically closed model must have a nonzero
remainder.  No such additional structure has been derived for F3.

## 11. Exact classification and next theorem

The result separates three logically different statements.

1. **Exact at finite coupling:** (AY.17) with the remainder retained;
   `[D_hat,N_def]=0`; the well-typed local-collar interaction `{Phi(S)}` with
   exact locked-endpoint port `U(1)^4` and finite `D_2`; the topology
   distinction (AY.38); the fixed-torus winding coefficient (AY.41).
2. **Asymptotic/prethermal:** replacement of the full dynamics by
   `U_d N_def+D_hat`, controlled by (AY.19) and (AY.34) only up to the stated
   quasi-exponential horizon.
3. **Not established:** `V_hat=0`, an exact all-time finite-coupling locked
   phase in the thermodynamic limit, selected-GNS gaplessness, an isotropic
   propagation law, gravity, or `G`.

The most direct next theorem is therefore not another finite perturbative
order.  It is a **record-authenticated finite-horizon application theorem**:
place an independently derived F3 ratio `U_d/h` and physical clock map into
(AY.13)--(AY.16), verify a declared record-retention window against
(AY.34)--(AY.35), and carry the exact port-conserving effective response plus
its explicit remainder into the already typed same-parent observable.  That
would decide whether the proven prethermal window covers the physical record
mission without presuming an exact phase or importing gravity machinery.
