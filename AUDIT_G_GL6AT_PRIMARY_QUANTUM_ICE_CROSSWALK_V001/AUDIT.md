# Distinct hostile audit — GL6AT primary quantum-ice crosswalk

**Target:** `LANE_CROSS_RFT_GRA_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001/`  
**Frozen author-manifest SHA-256:** `9acb465079e8dbf5fa15ffc6b08d8c4d34c65ee086f1f822e8276156f484ce61`  
**Frozen author-seal-file SHA-256:** `2a668784df7a4f3702e3e27c076ca8c9b2d15394e3f0e46135d6e440649ee575`  
**Disposition:** `PASS__EXACT_ORDER6_FPL_POINT_V_OVER_G_ZERO__PRIMARY_SOURCE_LABELS_SOUND__NO_SPECTRAL_OR_GRAVITY_PROMOTION`

## 1. Custody and audit independence

The eleven frozen GL6AT author files are pinned in
`AUDITED_TARGETS.sha256`. The author manifest and its one-row seal pass, as do
the twelve frozen GL6AO/GL6AP dependency rows. The author verifier passes
`174/174` in normal and optimized modes. The independent replay imports no
author or upstream verifier.

The four papers are not vendored. I therefore checked their versioned arXiv
or author-hosted published bytes directly, while this local seal authenticates
only the audit and author snapshots. GL6AT states this external-custody
boundary accurately.

## 2. Exact order-six model map

For `A3={x in Z^4: 1^T x=0}`, put parent `P_x` at `x` and child
`C_y` at `y-(1/4)1`. An incidence link `P_x--C_(x+e_a)` has bond vector

\[
 t_a=e_a-\tfrac14\mathbf 1,
 \qquad t_a\cdot t_b=\delta_{ab}-\tfrac14.
\]

The four bonds are a regular tetrahedral star and the two `A3` cosets are the
two FCC sublattices of the diamond net. This is an exact graph embedding up to
an arbitrary overall length; it supplies no physical lattice calibration.

With `n_e in {0,1}` and `S_e^z=n_e-1/2`, the degree-two lock is exactly
two dimers touching every degree-four diamond vertex. Thus occupied links are
fully packed loops, not the degree-one diamond dimer model. The normalization
is

\[
 Z_e=2n_e-1=2S_e^z.
\]

On a simple hexagon, toggling all six bits preserves the two incident cycle
occupancies at every cycle vertex only for the two alternating strings. It
exchanges those strings with unit matrix element and annihilates all other
locked local patterns, exactly matching

\[
 S_1^+S_2^-S_3^+S_4^-S_5^+S_6^-+\mathrm{h.c.}
\]

with no factor of two. The GL6AO coefficient therefore gives

\[
 H^{(6)}-E_{\rm scalar}P=-g\sum_cT_c,
 \qquad g=\frac{63}{8}\frac{h^6}{U_d^5}>0.
\]

Hermele--Fisher--Balents Eqs. (7)--(9) explicitly implement the patterned
sitewise `pi` rotation that changes the ring sign to the negative sign-free
form. Shannon et al. Eq. (2) sets the kinetic coefficient to one and calls the
dimensionless potential `mu`; restoring `g` gives exactly `mu=v/g`. Because
GL6AO's entire order-six diagonal census is a common scalar, `v=0` and hence
`v/g=0`. The RK point is instead `v/g=1`. No unknown order-eight term is
silently discarded into `v`.

## 3. Primary-source evidence replay

- [Hermele, Fisher, and Balents](https://doi.org/10.1103/PhysRevB.69.064404)
  map pyrochlore sites to diamond links with two dimers per site, give the
  alternating hexagon operator and sign change in Eqs. (6)--(9), and place the
  exact RK state at `V/J_ring=1` in Eqs. (10)--(11). Their conclusion calls the
  pure easy-axis point a speculation; GL6AT does not turn that into a theorem.
- [Shannon et al.](https://doi.org/10.1103/PhysRevLett.108.067204) study the
  same fully-packed-loop family by zero-temperature GFMC and finite ED. Their
  crossing near `mu=-0.50+-0.03`, the 320/640/1280-site GFMC data, and the
  80-site ED comparison support the liquid assignment at `mu=0` numerically,
  not rigorously.
- [Benton, Sikora, and Shannon](https://doi.org/10.1103/PhysRevB.86.075154)
  derive two transverse modes in their Gaussian lattice theory, a generic
  linear branch (quadratic at RK), Eq. (105) one-spin weight proportional to
  `omega(k)`, and the QMC-calibrated estimate
  `c=(0.6+-0.1)g a_0/hbar` at zero potential. GL6AT correctly labels the
  real-frequency atom and dispersion as effective-theory results, not direct
  microscopic real-time QMC.
- [Fu et al.](https://doi.org/10.1103/PhysRevB.96.035136) give the active
  `T2g` differences in published Eqs. (45a)--(45c), the electric bilinear in
  published Eq. (55c)/arXiv-v1 Eq. (57c), and a two-opposite-momentum-photon
  mechanism. Their published Eq. (70c)/arXiv-v1 Eq. (72c) contains
  `epsilon_p^2` times a nonzero small-`p` form factor. In three dimensions a
  linear `epsilon_p` therefore gives `Omega^4`, whereas their prose after
  published Eq. (73) describes `Omega^2`, the unweighted two-photon density
  power. GL6AT correctly imports no exponent and records the discrepancy.

All four DOI identities and versioned public links resolve to the named
primary research papers.

## 4. Operator and representation selection

The six unordered port pairs form `A1+E+T2` under `S4`. The three centered
opposite-pair sums form `E`, while

\[
 (M_{01}-M_{23},M_{02}-M_{13},M_{03}-M_{12})
\]

form `T2`. Direct enumeration of all six two-in/two-out assignments gives
zero for every local `T2` difference and a generally nonzero centered `E`
sum. Fu's three cubic `T2g` complementary-pair differences therefore cannot
be relabeled as the sealed local pair `E` observable. In addition, Fu use an
enlarged slave-particle Hilbert space and an uncontrolled leading gauge-field
expansion. Their result is evidence for a composite two-photon mechanism,
not an exact nonzero `N_E` matrix element, residue, threshold, exponent, or
pole in the projected `v/g=0` model.

The repaired single-link statement is also exact: `Z=2S^z`, so the Benton
spin structure factor and the `Z` structure factor differ by an overall
factor four. That changes normalization, not the vanishing `omega(k)` weight
or its explicitly conditional use.

## 5. Scope and hostile verdict

The fixed period-four quotient has no infrared sequence. Thermodynamic paper
results apply only after choosing the formal infinite interaction or a
compatible growing diamond-torus family, a flux/ground-state limit, and the
pure order-six operator. The crosswalk cannot be promoted through unknown
higher F3 orders.

No cited source proves the `v/g=0` thermodynamic phase, a microscopic gap or
gapless pole, a nonzero pair-`E` spectral measure, or survival at finite
`h/U_d`. No uncalibrated `A3` character is called physical momentum; no
physical speed, common cone, Lorentz invariance, graviton, Ricci/Einstein
dynamics, gravity, or Newton's constant `G` is derived.

**Hostile verdict: PASS.**
