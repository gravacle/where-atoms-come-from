# Minimal FU09b fixed-charge encoding and U(1) current lift

**Lane:** `LANE_GRA_GA_F3_Q4_FU09B_ENCODED_CHARGE_CURRENT_LIFT_V001`  
**Short name:** `ECCL`  
**Date:** 2026-08-28  
**Claim class:** exact finite encoded-reservoir construction; exact inherited-
flip unitary equivalence; exact internal U(1) transfer current and outer-port
balance; exact reservoir-placement nonuniqueness; conditional spatial-source
inheritance and explicit `FV-PURE` failure modes

**Not claimed:** a physical laboratory reservoir, an autonomous reset, a
spatial bond current, a vertex divergence, visible-QED completion, `T^{0j}`,
a metric Ward identity, continuum locality, gravity, or `G`.

## 1. Frozen question

FU proved that the prospective charge solder `q_a=q_* Z_a` makes the bare
inherited flip `X_a` charge nonconserving.  It required

\[
 \widetilde X_a=\sigma_a^+T_{a,-}+\sigma_a^-T_{a,+},
 \qquad [Q_R,T_{a,\pm}]=\pm2q_*T_{a,\pm},
 \qquad[Q_{\rm tot},\widetilde X_a]=0.             \tag{GA01}
\]

FU left open whether a fixed-total-charge encoded subspace can make
`X_tilde_a` exactly equivalent to the inherited `X_a`.  That question matters
to FV--FY: if the charged lift changes the source-off Hamiltonian or first
spatial source, `FV-PURE` and every projected H6 result must be recomputed.

This lane constructs the smallest exact finite witness, derives the current
and port equations it actually owns, and tests the source.  It adds no record
field and does not reinterpret charge current as stress or gravity.

## 2. One-link fixed-total-charge encoder

Let the link basis diagonalize `Z` with eigenvalues `-1,+1`.  Give one paired
reservoir qubit the charge operator

\[
 Q_R=q_*\bigl(|+\rangle\langle+|-|-\rangle\langle-|\bigr),
 \tag{GA02}
\]

and let `T_- = |-><+|`, `T_+=T_-^dagger`.  Define the isometry

\[
 V|-\rangle=|-\rangle_L|+\rangle_R,
 \qquad
 V|+\rangle=|+\rangle_L|-\rangle_R.               \tag{GA03}
\]

Its image `K_enc` is exactly the zero-total-charge eigenspace of
`q_*Z+Q_R`.  Direct multiplication gives

\[
 \boxed{
 (Z\otimes I)V=VZ,
 \qquad \widetilde X V=VX,
 \qquad(q_*Z+Q_R)V=0.}                             \tag{GA04}
\]

The image is invariant under `X_tilde`.  Hence `V` is a unitary map from the
original two-state link Hilbert space onto the fixed-charge encoded subspace.

### Theorem `ECCL-1` -- exact inherited-flip equivalence

The minimal link-plus-reservoir pair proves rather than assumes the encoded
equivalence demanded by FU09b.  Tensoring (GA03) over any finite link set
gives

\[
 \widetilde O\,V_E=V_EO                               \tag{GA05}
\]

for every polynomial `O` in the inherited link `Z_e,X_e` algebra, when each
`X_e` is replaced by its dressed transfer and each `Z_e` by `Z_e tensor I`.
The verifier checks all four q4 link algebras simultaneously on the complete
sixteen-state node.

This is an exact finite encoding, not proof that nature supplies one
independent paired reservoir at every physical link.

## 3. Exact internal U(1) current and boundary balance

For `H_flip=-h X_tilde`, Heisenberg evolution gives

\[
 \begin{split}
 \dot q_L&={2ihq_*\over\hbar}
 (\sigma^+T_- -\sigma^-T_+),\\
 \dot Q_R&=-\dot q_L.
 \end{split}                                         \tag{GA06}
\]

With current oriented from link to paired reservoir,

\[
 I_{L\to R}:=-\dot q_L,
 \qquad
 \boxed{\dot q_L+I_{L\to R}=0,
 \quad\dot Q_R-I_{L\to R}=0.}                     \tag{GA07}
\]

This is a nonzero scalar U(1) exchange current.  It is internal to a
link--reservoir pair; it is not a spatial current along a diamond bond.

If the reservoir couples to an explicit exterior charge through a conserving
exchange block, orient `I_(R->partial)=dot Q_partial`.  Then

\[
 \boxed{
 \dot q_L+I_{L\to R}=0,
 \quad
 \dot Q_R-I_{L\to R}+I_{R\to\partial}=0,
 \quad
 \dot Q_\partial-I_{R\to\partial}=0.}             \tag{GA08}
\]

The verifier constructs the eight-dimensional link/reservoir/exterior model
and checks all three equations exactly.  The same calculation exposes a
boundary condition: a reservoir-only outer exchange changes `Q_R` without
changing `Z`, so it leaks out of `K_enc`.  Exact equivalence to the inherited
closed hold therefore requires that outer reset/transport be off during the
hold or be replaced by a separately proved code-preserving composite process.

### Theorem `ECCL-2` -- current and active-port boundary

The encoded lift owns a complete internal charge-current identity and an
explicit outer-port term.  It does not own a simultaneously active outer port
that preserves the inherited encoded dynamics.  Preparation, response hold,
and reset are distinct dynamical stages unless a new code-preserving port
construction is proved.

## 4. Native placement is not selected

For one declared co-located allocation, place each reservoir with its FY link
midpoint and use the inherited `m=1` phases.  Then

\[
 \rho_m^L={1\over\sqrt{60}}\sum_e e^{ik\cdot r_e}q_{L,e},
 \qquad
 I_m^{L\to R}={1\over\sqrt{60}}\sum_e
 e^{ik\cdot r_e}I_{L\to R,e},                     \tag{GA09}
\]

obey

\[
 \dot\rho_m^L+I_m^{L\to R}=0.                    \tag{GA10}
\]

If the paired reservoir charge is assigned the same phase, the full
link-plus-reservoir `m=1` density is conserved.  If instead all reservoirs are
placed at one common support phase, its time derivative is nonzero and a
connector/boundary current is required.  The global `m=0` total charge remains
conserved in either allocation.

### Theorem `ECCL-3` -- reservoir-placement nonuniqueness

FU09b fixes the algebraic exchange but does not locate the reservoir.  The
co-located construction is one exact native-support witness, not a derived
physical placement.  A displaced placement is equally compatible with the
one-link charge algebra but needs an additional connector current at nonzero
momentum.  No vertex divergence or spatial bond current follows until the
physical reservoir, connector, and endpoint ownership are fixed.

In particular, (GA10) is a port-normal exchange equation.  It must not be
substituted for the stress continuity equation or called `T^{0j}`.

## 5. Exact first-spatial-source test

Freeze the bounded **`GA-CLOSED-FULL-CODE-SCALAR-HOLD`**:

1. one paired reservoir qubit is prepared with every link in `K_enc`;
2. every inherited flip is replaced by (GA01);
3. the outer port is off during the response hold;
4. every added reservoir, self-energy, support, and held-port term restricts
   to one common scalar identity on the **full** encoded Hilbert space,
   including all projected `P` and virtual/off-ice `Q` sectors, at source off;
   and
5. its first derivative with respect to the frozen spatial source `j_ij` is
   absent or one declared common identity/reference source on the full encoded
   Hilbert space.  Source independence alone is insufficient: a constant
   `mu Q_R` is source independent but restricts to `-mu q_* Z` and changes the
   source-off gaps.

On this hold, lift the source family by the same substitution as (GA05).  For
one common scalar `c[j]`, the exact intertwining statement is

\[
 \widetilde H[j]V_E=V_E\bigl(H[j]+c[j]I\bigr),
 \qquad
 \widetilde Q^{ij}V_E=V_E\bigl(Q^{ij}+q_{\rm id}^{ij}I\bigr),
 \quad q_{\rm id}^{ij}=-2\partial_{j_{ij}}c.        \tag{GA11}
\]

The verifier checks separately that the nonidentity parts of the complete
six-coordinate Coulomb pair source and hopping-numerator source intertwine,
and that a common identity source remains an identity.  The nonidentity source
Gram matrix and rank are unchanged.

Let `z` be the Feshbach energy reference.  The common scalar cancels every
virtual gap; after shifting the reference consistently, every invariant
encoded Feshbach projection and resolvent obeys

\[
 \widetilde H_{\rm eff}[j;z+c[j]]V_P
 =V_P\bigl(H_{\rm eff}[j;z]+c[j]I_P\bigr).         \tag{GA12}
\]

Thus all nonidentity H2/H4/H6 coefficients, folds, virtual gaps, commutators,
transition matrix elements, connected responses, ranks, and native phases are
unchanged on the encoded hold.  The literal Hamiltonian and first source may
carry the displayed identity/reference shifts; those are not silently called
zero.

### Theorem `ECCL-4` -- conditional `FV-PURE/FY` preservation

`GA-CLOSED-FULL-CODE-SCALAR-HOLD` preserves the frozen source-off Hamiltonian
and first spatial source **modulo one full-code identity/reference shift**.
It therefore preserves every nonidentity FV--FY H6 result listed after
(GA12).  The charge-conserving flip need not by itself invalidate the H6
source calculation, but literal equality of the uncentered Hamiltonian/source
is claimed only when `c[j]=0`.

The condition is substantive, not automatic.  Exact counterexamples are:

| reservoir/port term | encoded restriction | consequence |
|---|---|---|
| local `Q_R^2` | `q_*^2 I` | identity; compatible |
| local chemical potential `mu Q_R` | `-mu q_* Z` | new nonidentity one-link term |
| shared `Q_(R,1)Q_(R,2)` | `q_*^2 Z_1Z_2` | new cross-link pair term |
| active reservoir-only outer exchange | leaves `K_enc` | source-off parent changes |
| strain-dependent `T_+/-`, support, or port | additional `dH/dj` | source must be re-ranked |

Accordingly the minimal witness proves existence of a preserving charged
lift.  It does not prove that an arbitrary grounded conductor, common
reservoir, reset mechanism, or physical support satisfies `FV-PURE`.

## 6. What has and has not closed

The lane closes the finite algebraic question left by FU09b:

\[
 \boxed{
 \text{bare inherited link dynamics}
 \ \cong\ 
 \text{charge-conserving dressed dynamics on }\mathcal K_{\rm enc}.}
 \tag{GA13}
\]

It also proves the internal current, the outer boundary term, the allocation
ambiguity, and the exact conditions under which the FY source survives.

It does not discharge FU's full physical solder `S3--S4`: reservoir hardware,
placement, common grounding, source work, support/recoil, and a code-preserving
active port remain physical inputs.  It supplies no spatial propagation of
charge, no electromagnetic gauge-field Ward packet, no stress current, and no
gravity result.

## 7. Disposition

`EXACT_MINIMAL_FIXED_TOTAL_CHARGE_ENCODER__DRESSED_FLIP_UNITARILY_EQUIVALENT_TO_INHERITED_X_ON_THE_ENCODED_SUBSPACE__EXACT_NONZERO_LINK_TO_RESERVOIR_U1_CURRENT_AND_OUTER_PORT_BALANCE__ACTIVE_RESERVOIR_ONLY_PORT_BREAKS_THE_HOLD_CODE__RESERVOIR_PLACEMENT_AND_NONZERO_MOMENTUM_CONNECTOR_CURRENT_NOT_SELECTED__CLOSED_FULL_CODE_SCALAR_HOLD_PRESERVES_NONIDENTITY_FV_PURE_AND_FY_RESULTS_MODULO_REFERENCE_IDENTITY__COMMON_SCALAR_CANCELS_VIRTUAL_GAPS_UNDER_CONSISTENT_ENERGY_SHIFT__SOURCE_INDEPENDENCE_ALONE_IS_INSUFFICIENT__BIAS_SHARED_CHARGING_ACTIVE_PORT_OR_NONIDENTITY_SOURCE_DERIVATIVE_INVALIDATES_UNCONDITIONAL_INHERITANCE__NO_SPATIAL_BOND_CURRENT_T0J_METRIC_WARD_GRAVITY_OR_G`
