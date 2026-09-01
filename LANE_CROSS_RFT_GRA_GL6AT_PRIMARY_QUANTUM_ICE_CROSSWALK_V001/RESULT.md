# GL6AT result — exact effective-model crosswalk and evidence ceiling

## 1. Sealed input and scalar removal

The sealed and independently audited GL6AO result is

\[
 H_{\rm eff}=E_{\rm scalar}P
 -{63\over8}{h^6\over U_d^5}\sum_{c\in\operatorname{Hex}}T_c
 +O(h^8/U_d^7),                                           \tag{AT01}
\]

where `E_scalar` includes the common order-two, order-four, and order-six
diagonal shifts.  Define

\[
 g={63\over8}{h^6\over U_d^5}>0.                           \tag{AT02}
\]

Removing `E_scalar P` changes absolute energies but not eigenvectors,
energy differences, or locked-sector dynamics.  The literature crosswalk is
to the displayed order-six operator in (AT01), not to the uncomputed
`O(h^8/U_d^7)` completion.

## 2. Exact graph isomorphism

GL6AN/GL6AO use parent vertices `P_x`, `x in A3`, and child vertices
`C_(x+e_a)`, with one edge from `P_x` to each `C_(x+e_a)`, `a=1,2,3,4`, where

\[
 A_3=\{x\in\mathbb Z^4:\mathbf1^Tx=0\}.
\]

Embed the parent label `x` at `x` in the three-plane `1^T r=0`, and embed a
child label `y`, `1^T y=1`, at `y-1/4`.  The four bond vectors are

\[
 t_a=e_a-{1\over4}\mathbf1,
 \qquad t_a\mathbin\cdot t_b=\delta_{ab}-{1\over4}.        \tag{AT03}
\]

They are the four vertices of a regular tetrahedral star.  Thus the infinite
parent/shared-child incidence is abstractly, and after this optional
embedding geometrically, the standard diamond net: two `A3`/FCC cosets joined
by four tetrahedral bonds.  The graph isomorphism is exact.  Its physical
length scale is not inherited; multiplying (AT03) by a lattice spacing is an
extra calibration.

## 3. Exact Hilbert-space and operator map

For each diamond link set

\[
 n_e=1\ \longleftrightarrow\ \hbox{dimer present}
 \ \longleftrightarrow\ S_e^z=+1/2,
 \qquad
 n_e=0\ \longleftrightarrow\ S_e^z=-1/2.                  \tag{AT04}
\]

The GL6AO lock `sum_(e incident v) n_e=2` becomes two dimers touching every
diamond vertex, equivalently the two-in/two-out ice manifold after choosing
the usual bipartite orientation.  Occupied links form a fully-packed loop
configuration.  This is not the standard degree-one diamond quantum dimer
model.

For an elementary six-cycle `c`,

\[
 T_c=P\left(\prod_{e\in c}X_e\right)P                    \tag{AT05}
\]

has unit matrix element between the two alternating occupations and vanishes
on a nonalternating locked configuration.  This is exactly the ring operator

\[
 S_1^+S_2^-S_3^+S_4^-S_5^+S_6^-+\mathrm{h.c.}            \tag{AT06}
\]

used in the quantum-ice papers.  Hermele et al.'s explicit sitewise `pi`
rotation changes its sign, so the sign convention `-g sum_c T_c`, `g>0`, is
the sign-free convention used by Shannon et al.

## 4. Exact comparison parameter

Write the comparison family only as a coordinate system,

\[
 H(v)=-g\sum_cT_c+v\sum_cF_c,                              \tag{AT07}
\]

where `F_c` projects onto the two flippable alternating states of `c`.  Since
GL6AO proves that no configuration-dependent diagonal term survives at order
six,

\[
 \boxed{v/g=0\quad\hbox{for GL6AO}.}                       \tag{AT08}
\]

The exactly soluble RK point is `v/g=1`.  Equation (AT07) does not add `v` to
the F3 model; it merely locates the already derived `v=0` Hamiltonian in the
published one-parameter family.

## 5. Conditions on importing thermodynamic results

The paper-level phase crosswalk is valid only if all of the following are
declared:

1. use the formal infinite order-six interaction of GL6AO, or a growing
   sequence of standard diamond tori whose local hexagons agree with it;
2. retain the exact degree-two Hilbert space and only the elementary
   alternating-hexagon kinetic operator at the order being compared;
3. remove only the common scalar and set the flippability potential exactly
   to zero;
4. choose a thermodynamic flux/sector and stationary ground-state limit
   comparable to the one used in the numerical papers; and
5. do not promote the comparison through the unknown higher-order F3 terms.

The fixed `Q_4` theorem by itself has no literal infrared limit.  Conversely,
the published thermodynamic calculations do not establish the finite-`h/U_d`
phase of the parent F3 model.

## 6. Phase, gap, and dispersion disposition

- **Exact at `v/g=0`:** only the Hamiltonian crosswalk (AT01)--(AT08).  No
  exact ground state, thermodynamic phase, or spectral gap is known from the
  admitted papers at this point.
- **Exact at a different point:** at `v/g=1`, the RK Hamiltonian has exact
  equal-weight ground states within connected sectors.  Its special
  long-wavelength theory has quadratic rather than the generic linear photon
  dispersion.  Neither fact is an exact statement about `v/g=0`.
- **Numerical at `v/g=0`:** zero-temperature GFMC, supplemented by small-cluster
  exact diagonalization, finds flux-energy scaling consistent with a `U(1)`
  liquid.  Shannon et al. place the level crossing near `v/g=-0.50+-0.03` and
  report unambiguous liquid evidence for `v/g>-0.3`, so zero lies well inside
  their numerically inferred liquid region.
- **Effective/QMC-calibrated at `v/g=0`:** Benton et al. infer two transverse
  modes with `omega(k)~c|k|` and estimate

  \[
  c=(0.6\mathbin{+-}0.1)g a_0/\hbar.                       \tag{AT09}
  \]

  Their `a_0` and wave vector belong to a chosen physical diamond embedding.
  GL6AO supplies neither that calibration nor physical momentum.
- **Gap ceiling:** the numerical and effective results support a gapless
  collective branch in the thermodynamic comparison model, but do not prove a
  microscopic pole or a zero many-body spectral gap.  Degree-violating
  electric charges are outside the projected GL6AO Hilbert space.  Gapped
  spinon/charge or monopole statements in enlarged spin-ice descriptions are
  therefore not imported as spectra of the pure locked operator.

## 7. Operator spectral-overlap crosswalk

### Single-link observable

Under (AT04), the microscopic link observable `Z_e=2n_e-1=2S_e^z` is twice,
and hence proportional to, the ice spin component.  Benton et al.'s Gaussian
lattice field theory predicts near its long-wavelength origin

\[
 S^{yy,zz}(k,\omega)\ \propto\
 \omega(k)\,\delta(\omega-\omega(k)).                     \tag{AT10}
\]

Thus their one-photon atom has weight tending to zero with frequency.  This is
a useful controlled warning that a soft pole need not have nonvanishing
visibility.  Equation (AT10) is an effective-theory prediction whose
equal-time integral was compared with QMC; it is not a direct microscopic
real-frequency computation.

### Pair/tensor composites

Fu et al.'s Raman vertex contains local bilinears
`S^z_mu S^z_nu`, represented in their gauge description by
`E_mu E_nu`.  At leading Gaussian order this operator creates a pair of
opposite-momentum photons in the `T2g` Raman channel.  No low-energy power is
imported: the paper's prose and its displayed three-dimensional integral have
an unresolved power-counting discrepancy.

This does **not** close the sealed pair-`E` question.  Their three active
components are complementary-pair differences,

\[
 R_{01}-R_{23},\quad R_{02}-R_{13},\quad R_{03}-R_{12},    \tag{AT11}
\]

which form the three-dimensional cubic/local-pair `T2` sector.  The sealed
GL6AP pair `E` sector is instead the centered plane of complementary-pair
sums

\[
 M_{01}+M_{23},\quad M_{02}+M_{13},\quad M_{03}+M_{12}.    \tag{AT12}
\]

Moreover, Fu et al. use an enlarged XXZ slave-particle/gauge treatment and
state approximation limitations.  Their result is mechanism-level evidence
that a tensor bilinear can overlap a two-photon continuum, not an exact
nonzero overlap, infrared exponent, residue, threshold, or pole theorem for
(AT12).

Hermele et al. also obtain, in Gaussian Coulomb theory, a scalar loop kinetic
energy density whose connected equal-time correlation falls as `1/R^8`.
That is another composite-field mechanism, not a dynamic pair-`E` spectral
calculation.

## 8. Final ceiling

The literature raises the pure order-six locked interaction from an unnamed
loop Hamiltonian to an exact published quantum-ice model point, and provides
strong but numerical/effective evidence for its `U(1)` liquid behavior.  It
does not supply a rigorous `v/g=0` phase or gap theorem, an exact local
pair-`E` overlap, a state-independent pole, physical momentum without an
embedding, a physical cone, Ricci/Einstein dynamics, gravity, or `G`.
