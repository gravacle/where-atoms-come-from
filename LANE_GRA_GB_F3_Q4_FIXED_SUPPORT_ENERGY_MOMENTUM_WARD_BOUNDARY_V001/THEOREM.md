# Fixed-support energy continuity and momentum-Ward boundary

**Lane:** `LANE_GRA_GB_F3_Q4_FIXED_SUPPORT_ENERGY_MOMENTUM_WARD_BOUNDARY_V001`  
**Short name:** `FSEMWB`  
**Date:** 2026-08-28  
**Claim class:** exact finite local-Hamiltonian energy-continuity theorem;
exact finite-translation theorem; sharp physical-momentum ownership and
stress-Ward dependency theorem

**Not claimed:** a unique local energy density, a physical momentum density,
`T^{0j}`, a stress tensor, a metric Ward identity, a continuum limit, gravity,
or `G`.

## 1. Frozen question

FY supplies a complete first spatial source on the fixed `Z30` diamond
support.  FZ proves that the temporal current, physical discrete divergence,
and contact packet needed for a stress Ward identity are not present.  GA then
constructs an exact U(1) link--reservoir charge current, while explicitly
showing that reservoir placement and spatial connector allocation are not
selected.

This lane asks one narrower physics question on that unchanged parent:

> What exact energy continuity follows from the local H6 Hamiltonian, and
> does the finite translation of the externally fixed support supply the
> missing local physical momentum current?

The answer separates energy from momentum.  The inherited ring Hamiltonian
has an exact local term-current continuity law.  Its `Z30` translation is also
an exact symmetry.  But the latter supplies only a global finite-group
representation label.  It does not select a local momentum density, a
physical support-recoil operator, or `T^{0j}`.

## 2. Exact inherited energy continuity

On the frozen 180-state FO component, write the H6 generator as the inherited
sum of 120 Hermitian hexagon terms,

\[
 H_6=\sum_{a=1}^{120}h_a .                         \tag{GB01}
\]

For the declared ring-term allocation define, with current oriented from
term `a` to term `b`,

\[
 J^E_{a\to b}={i\over\hbar}[h_a,h_b]
              =-J^E_{b\to a}.                     \tag{GB02}
\]

Heisenberg evolution gives the exact operator identity

\[
 \boxed{\dot h_a+\sum_bJ^E_{a\to b}=0.}            \tag{GB03}
\]

For every subset `R` of ring terms, let
`H_R=sum_(a in R) h_a` and `H_Rbar=H_6-H_R`.  Then

\[
 \boxed{
 \dot H_R+J^E_{R\to\bar R}=0,
 \qquad
 J^E_{R\to\bar R}={i\over\hbar}[H_R,H_{\bar R}].} \tag{GB04}
\]

The verifier reconstructs all 120 terms exactly, checks their sum against
the frozen H6 matrix, verifies (GB03) for every ring, exhibits a nonzero pair
current and a nonzero twelve-term boundary current, and verifies total energy
conservation.  It
also checks all 5,700 edge-disjoint ring pairs: their commutators vanish, so
the inherited term-current graph adds no current between disjoint supports.

### Theorem `FSEMWB-1` -- exact term-local energy balance

Equations (GB02)--(GB04) are an exact local Hamiltonian-energy continuity
theorem on the frozen support.  They are not a continuum energy-momentum Ward
identity.  The index `a` labels an inherited ring-energy term, not a derived
spatial vector component.

The result is also allocation-dependent.  For any Hermitian `C`, replacing

\[
 h_a\mapsto h_a+C,
 \qquad h_b\mapsto h_b-C                           \tag{GB05}
\]

leaves `H_6` fixed but can change the attributed pair current.  The inherited
ring decomposition is a lawful and local convention, not a proof that it is
the unique physical energy density.

GA preservation is used only in its repaired form.  A single common scalar
`cI` on the complete encoded `P+Q` Hilbert space is a reference shift and
obeys

\[
 [H_6+cI,h_a]=[H_6,h_a].                            \tag{GB06}
\]

Thus it changes none of (GB02)--(GB04).  A sector-dependent term, a reservoir
bias descending to `Z`, or an active port is not covered by (GB06).

## 3. Finite translation is not a local momentum current

The frozen support has an exact unitary cell translation `U_tau` satisfying

\[
 U_\tau^{30}=I,
 \qquad [H_6,U_\tau]=0.                             \tag{GB07}
\]

This proves conservation of the global `Z30` representation label.  It does
not provide a unique infinitesimal generator.  If `Pi_m` are the exact
character projectors and `theta_m=2 pi m/30`, then every

\[
 K_{\boldsymbol n}=\sum_{m=0}^{29}
       (\theta_m+2\pi n_m)\Pi_m,
 \qquad n_m\in\mathbb Z,                            \tag{GB08}
\]

is Hermitian and exponentiates to the same `U_tau`.  The verifier constructs
the actual decomposition into six length-30 state orbits, constructs two
distinct logarithm branches explicitly on one orbit type, and shows that a
standard spectral branch is nonlocal in the cell-orbit basis.

Equation (GB08) does not prove that no enlarged or continuum completion could
possess a local momentum generator.  It proves the needed derivability
boundary: the supplied finite-group action alone neither selects such a
generator nor decomposes it into local densities and fluxes.

The FO response Hilbert space contains only incidence configurations.  Its
diamond coordinates and bond vectors are fixed numerical support data.  It
contains no canonical support positions and momenta, no recoil factor, and no
autonomous support or clamp Hamiltonian.

### Theorem `FSEMWB-2` -- fixed-support translation boundary

`[H_6,U_tau]=0` is a genuine exact symmetry but is insufficient to derive a
local physical `T^{0j}`.  Calling its global finite-group label physical
momentum would silently add both an infinitesimal translation action and a
local density allocation that the parent does not own.

## 4. Energy conservation is not spatial momentum conservation

For the closed, time-independent H6 block, (GB03) is exact.  If an interaction
with support, reservoir, controller, or boundary factors is activated, the
bulk energy ledger instead has the form

\[
 {dH_{\rm bulk}\over dt}
 = {\partial H_{\rm bulk}\over\partial t}
 +{i\over\hbar}[H_{\rm int},H_{\rm bulk}],          \tag{GB09}
\]

with the opposite exchange assigned to the owned exterior factors in a
closed completion.  The current GA construction owns scalar U(1) charge
transfer.  It does not assign the energy or impulse carried by a spatial
connector, reservoir placement, exterior port, or clamp.  A charge current is
not thereby `T^{0j}`.

Likewise, a completed momentum ledger would require an operator identity of
the schematic form

\[
 \dot{\mathcal P}_{\rm incidence}
 =-\dot{\mathcal P}_{\rm support}
  -\dot{\mathcal P}_{\rm reservoir/field}
  -\dot{\mathcal P}_{\rm controller}
  -\dot{\mathcal P}_{\rm boundary}.                 \tag{GB10}
\]

None of the operators on the right of (GB10), nor a physical incidence
momentum on its left, is defined in the fixed-support response parent.  The
equation is therefore a required ownership ledger, not a result asserted by
this lane.  In a clamped realization, the omitted support/controller is
precisely where unmatched impulse can go.  In an emergent-space realization,
an alternative dynamical sector may own that impulse, but it must be exhibited
rather than inferred from the cell permutation.

### Theorem `FSEMWB-3` -- exact stress-Ward dependency

The frozen parent proves local Hamiltonian-energy continuity but neither
proves nor falsifies local spatial-momentum conservation.  Its missing
`T^{0j}` is a physical ownership gap: support/recoil, reservoir/field,
controller, and boundary momentum exchange have not been made operators in
one closed parent.  This is independent of the positive U(1) current result
in GA.

## 5. Minimum physical completion

A true stress Ward packet requires, without adding any new record mechanism:

1. **A momentum-owning realization of the support.**  Supply dynamical
   positions and conjugate recoil momenta for the diamond realization, or an
   explicit alternative dynamical sector that realizes translations and owns
   the equal-and-opposite impulse.  Fixed coordinate labels are insufficient.
2. **One autonomous closed Hamiltonian.**  Include reservoirs/fields,
   connectors, support or clamp, controller/clock, and boundary ports, with
   their energy and momentum exchanges counted once.
3. **One complete source family.**  Derive `T^{00}`, `T^{0i}`, and `T^{ij}`
   from the same `H[j_00,j_0i,j_ij; ports]`, retaining all second derivatives
   and equal-time contacts.
4. **A native local divergence.**  Derive the finite bond/port divergence
   from those owned currents, prove its boundary law, and only then compare it
   with FY's supplied embedding contraction.
5. **Projection after ownership.**  Carry the complete density/current,
   contact, recoil, and port packet through Feshbach reduction together.

Only after these steps can one test the stress identity identified in FZ,

\[
 {i\over\hbar}[H,T^{0j}_m]+(\Delta_m)_iT^{ij}_m=0, \tag{GB11}
\]

and its retarded contact relation.  Discrete translation symmetry by itself
does not fill any missing slot in (GB11).

## 6. Disposition and ceiling

This lane closes the strongest positive result available on the unchanged
fixed-support parent: exact local ring-energy continuity.  It also closes a
false shortcut: a conserved `Z30` cell-translation label is not already the
missing local momentum current.

The next physics step is a translation-owning support/field/port completion,
or an explicit emergent local generator derived from a larger dynamical F3
parent.  The present lane adds no interaction and makes no gravity promotion.
