# Finite-family native-source and common-cone boundary

**Lane:** `LANE_GRA_GC_F3_Q4_FINITE_FAMILY_COMMON_CONE_BOUNDARY_V001`  
**Short name:** `FFCCB`  
**Date:** 2026-08-28  
**Claim class:** exact covering-family H6 inventory and local-coefficient
inheritance; exact reciprocal-volume, small-momentum, refinement, and TT
projector scaling; exact response-identifiability boundary and finite-size
acquisition target

**Not claimed:** a complete native source on every `G_L` momentum, nested
many-body Hilbert spaces, gap or residue scaling, a Ward-derived tensor pole,
a common physical cone, gravity, or `G`.

## 1. Frozen question

Five existing results now meet at one finite-family boundary:

- FS owns the covering-matched periodic diamond family
  `G_L`, `L=5*2^r`, but explicitly not inclusions of many-body Hilbert spaces.
- FD owns the exact `A3` affine volume and the conditional dispersion/common-
  cone formulas once a massless collective action is supplied.
- FV owns the environment-independent local H6 ring-source coefficients and
  six exact off-shell source witnesses on every covering member.
- FY owns the complete native-support H6 source and sampled response at one
  nonzero momentum of one 30-cell quotient.
- FL proves that the direct even composite has no isolated helicity-two pole
  at the Gaussian-Maxwell fixed point; its baseline is a multiphoton
  continuum.

This lane asks:

> What volume, momentum, projector, and local-source scaling already follows
> exactly on `G_L`, and what response scaling is still required before a
> common massless tensor cone can be claimed?

No interaction is added, and no thermodynamic ice Hilbert space is enumerated.

## 2. Exact family and H6 coefficient inheritance

For cells `x in (Z/LZ)^3`, FS defines two vertices `A_x,B_x` and the four
links

\[
 e_{x,a}:A_x\longrightarrow B_{x+s_a},
 \qquad
 (s_0,s_1,s_2,s_3)=(0,e_1,e_2,e_3).               \tag{GC01}
\]

Hence

\[
 |V_L|=2L^3,
 \qquad |E_L|=4L^3,
 \qquad |\mathcal S_L|=6L^3,                       \tag{GC02}
\]

where `S_L` contains the two vertex and four link-midpoint support types per
cell.

For each missing q4 label `d`, the other three labels `(a,b,c)` generate the
elementary hexagon

\[
\begin{aligned}
A_x&\to B_{x+s_a}\to A_{x+s_a-s_b}
\to B_{x+s_a-s_b+s_c}\\
&\to A_{x-s_b+s_c}\to B_{x+s_c}\to A_x .          \tag{GC03}
\end{aligned}
\]

There are exactly `L^3` distinct translates for each `d`, hence `4L^3`
elementary H6 rings.  Reduction modulo `L` sends the rings of `G_(2L)` to
those of `G_L` with fiber eight.

At the Coulomb slope `lambda=-1/2`, FV's Hermitian local ring tensor is

\[
 \overline T_d={651\over16}I-{567\over16}D_d,       \tag{GC04}
\]

independently of the external local ice assignment.  Together with FV's two
direct-pair `E` witnesses, the four missing-label rows have

\[
 \operatorname{rank}(E_1,E_2,\overline T_0,\ldots,
 \overline T_3)=6,
 \qquad
 \det W=-{4678629417\over256}.                     \tag{GC05}
\]

Because (GC04) is a bounded local coefficient indexed only by the missing
label, every marked local witness pulls back unchanged under the graph cover.
Thus the source has an exact volume-uniform **local coefficient ancestry** and
an extensive inventory: two vertex/four link support slots and four H6 ring
orientations per cell.

### Theorem `FFCCB-1` -- local source survives every cover

The exact direct-pair and H6 coefficient witnesses do not dilute or change
under `G_(2L)->G_L`.  The local off-shell tensor list remains rank six, and
the number of translated H6 opportunities scales as `4L^3`.

This is not spectral inheritance.  A graph cover is not an isometry or an
inclusion of the ice Hilbert spaces.  It does not map ground states,
excitation gaps, or residues between volumes.

There is a second strict ceiling.  FY assigned every pair, virtual gap,
numerator, fold, and ring derivative to a physical support on its 30-cell
cyclic quotient.  The complete diagonal proportionalities and within-cell
support phases were proved only at `m=0,1` there.  FV and FS license (GC04)
and the local support inventory on `G_L`; they do not by themselves license
copying FY's complete `m=1` diagonal ledger to every `G_L` momentum.  That
bounded native-ledger lift remains the first calculation in the acquisition
target below.

### Exact FS-to-FD affine join

The graph-to-affine identification is not inferred from the volume count.
Identify GC's zero-shift q4 label `0` with label `4`, write FD's tetrahedral
vectors as `n_a`, and put `alpha_(i4)=n_4-n_i`.  On the universal cover of
the periodic graph define

\[
 X(A_x)=a_*\sum_{i=1}^3x_i\alpha_{i4},\qquad
 X(B_x)=X(A_x)-a_*n_4 .                              \tag{GC05a}
\]

Then the zero-shift edge has displacement `-a_*n_4`, while for `i=1,2,3`

\[
 X(B_{x+e_i})-X(A_x)
 =a_*\alpha_{i4}-a_*n_4=-a_*n_i .                  \tag{GC05b}
\]

Thus the four FS edge labels are exactly the four FD diamond directions and
the FS translation steps are exactly
`a_*alpha_(14),a_*alpha_(24),a_*alpha_(34)`.  This constructs the affine
join used below.  Periodic reduction makes an affine torus; turning its free
scale and cell map into physical rods and volumes remains FD's independent
scale/cell-binding premise.

## 3. Exact momentum and affine-volume scaling

Use FD's primitive `A3` basis

\[
 A=a_*[\alpha_{14}\ \alpha_{24}\ \alpha_{34}].     \tag{GC06}
\]

Its dimensionless Gram matrix and inverse are

\[
 G={4\over3}
 \begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix},
 \qquad
 G^{-1}={3\over16}
 \begin{pmatrix}3&-1&-1\\-1&3&-1\\-1&-1&3\end{pmatrix}.    \tag{GC07}
\]

The periodic FD-affine spatial volume is therefore

\[
 \boxed{V_L^{\rm aff}=L^3v_3
 ={16\over3\sqrt{3}}(La_*)^3.}                    \tag{GC08}
\]

This becomes a physical volume only after FD's independent scale/cell
binding.  The graph and source theorem fixes the dependence on `L` and the
free affine scale `a_*`; it does not calculate an absolute length.

The `L^3` exact translation characters are indexed by classes
`[n] in (Z/LZ)^3`.  Choose for each class a minimum-norm Brillouin-zone
representative `hat n`; with the primitive basis (GC06), its reciprocal
vector is

\[
 k_{L,[n]}={2\pi\over L}A^{-\mathsf T}\widehat n,
 \qquad
 |k_{L,[n]}|^2=\left({2\pi\over La_*}\right)^2
 \widehat n^{\mathsf T}G^{-1}\widehat n .          \tag{GC09}
\]

For every integer representative `n=(n_1,n_2,n_3)`, put

\[
 q(n)=4\sum_i n_i^2-\left(\sum_i n_i\right)^2
 =\sum_i n_i^2+\sum_{i<j}(n_i-n_j)^2,
 \qquad n^{\mathsf T}G^{-1}n={3\over16}q(n).        \tag{GC10}
\]

The sum-of-squares form is positive, and
`q(n)=-[sum_i n_i]^2 mod 4`, so a positive value is `0` or `3 mod 4` and is
at least three.  Equality `q=3` forces `sum_i n_i^2<=3` and gives exactly
`+/-e_1,+/-e_2,+/-e_3,+/-(1,1,1)`, eight vectors.  Equality `q=4` gives exactly
the six same-sign two-coordinate vectors.  These representatives are
distinct modulo every `L>=5`.  Hence the shortest nonzero Brillouin-zone
value is `9/16`, and the next is `3/4`.  Consequently

\[
 \boxed{|k_{\min}(L)|={3\pi\over2La_*}.}            \tag{GC11}
\]

The cover and infrared limits are different operations.  A character `n` on
`G_L` pulls back to `2n` on `G_(2L)`.  The new odd-index modes of `G_(2L)`,
including its lowest nonzero mode, have no base ancestor.  Local coefficient
inheritance therefore cannot determine their response.

### Theorem `FFCCB-2` -- exact IR/refinement trichotomy

Under one cover step:

1. fixed `a_*` gives `V_(2L)^aff=8V_L^aff` and
   `k_min(2L)=k_min(L)/2`: this is an infrared/thermodynamic step;
2. simultaneous `a_*->a_*/2` leaves both the affine volume and `k_min`
   fixed: this is geometric refinement, not an infrared limit; and
3. for `L_r=5*2^r`, `a_r=a_0 2^{-sr}`, both `a_r->0` and
   `L_ra_r->infinity` hold exactly when `0<s<1`.

Thus an infrared response test requires `L_ra_r->infinity`, equivalently
`k_min->0`.  If one also claims simultaneous mathematical refinement of the
affine support, the stronger joint conditions are

\[
 \boxed{a_r\to0,
 \qquad L_ra_r\to\infty.}                          \tag{GC12}
\]

FD's fixed-affine-volume refinement (`s=1`) cannot by itself supply new
momenta approaching zero.  Conversely, fixed `a_*` is sufficient for a
thermodynamic long-wavelength test even though it is not an `a_*->0`
mathematical refinement.

## 4. Exact kinematic TT scaling

For every nonzero affine reciprocal vector `k` (physical only after the same
scale binding used for (GC08)), define

\[
 P_{ij}(k)=\delta_{ij}-{k_ik_j\over|k|^2},
\]

\[
 \Pi^{TT}_{ij,kl}(k)
 ={1\over2}\left(P_{ik}P_{jl}+P_{il}P_{jk}-P_{ij}P_{kl}\right). \tag{GC13}
\]

The verifier constructs (GC13) in exact rational arithmetic on a shortest
reciprocal ray.  It is idempotent and has rank two on the six-dimensional
symmetric-tensor coordinate space.  Since

\[
 P(\rho k)=P(k),
 \qquad
 \Pi^{TT}(\rho k)=\Pi^{TT}(k)\quad(\rho\ne0),       \tag{GC14}
\]

the same exact two-dimensional kinematic quotient persists as fixed `n` is
sent toward `k=0` by growing `L a_*`.

The six `A3` roots also obey

\[
 \sum_{a<b}\alpha_{ab}\alpha_{ab}^{\mathsf T}
 ={16\over3}I.                                      \tag{GC15}
\]

This is enough for FD's isotropic quadratic principal symbol **if** a healthy
massless collective action is supplied.  Equations (GC13)--(GC15) do not
derive that action, a pole, a Ward identity, or a residue.  Kinematic TT rank
two at every volume is compatible with either a massive response or no
source-visible excitation at all.

## 5. Exact response target for a common massless tensor cone

After the complete native `G_L` source is assembled, choose exact plus/cross
basis tensors in the image of (GC13).  If `q_x(n)` denotes the completed
six-support within-cell block and `x` runs over the `L^3` primitive cells,
FY's vertex-count convention is

\[
 Q_{L,\mathrm{cnt}}^{ij}(n)={1\over\sqrt{2L^3}}
 \sum_x e^{2\pi i n\cdot x/L}q_x^{ij}(n).          \tag{GC16a}
\]

This is count-normalized, not yet physically volume- or wavefunction-
normalized.  At fixed `a_*`, vertex count is proportional to affine volume,
but along `a_*->0` the conversion is scale dependent.  More generally a
physical or canonical operator must be fixed prospectively as

\[
 Q_{L,\mathrm{phys}}(n)=Z_Q(a_*,k)Q_{L,\mathrm{cnt}}(n),
 \qquad R_{\mathrm{phys}}=|Z_Q|^2R_{\mathrm{cnt}},  \tag{GC16b}
\]

where `Z_Q` must follow from a physical cell measure, a source/operator map,
or a Ward/canonical normalization before the response is inspected.  It may
also contain a known momentum power when the native observable is a
derivative or field-strength source.  Post-result reweighting is not a
normalization theorem.

For either prospectively fixed convention, define the positive-time
imaginary-time connected correlator

\[
 C^{TT}_{AB,L}(\tau,n)
 =\langle Q^{TT}_{A,L}(\tau,n)
          Q^{TT\dagger}_{B,L}(0,n)\rangle_c
 =\sum_s R^{(s)}_{AB,L}(n)e^{-\Delta_{s,L}(n)\tau}. \tag{GC16}
\]

Here `s` labels a full spectral projector at one energy, not necessarily one
eigenvector, so a degenerate projector can have residue rank two.  Pin the
quadratic-denominator convention by the even Euclidean cosine transform:

\[
 G^{\cos}_{E,AB}(\omega,n)
 =2\int_0^\infty d\tau\,\cos(\omega\tau)C_{AB}(\tau,n)
 =\sum_s{2\Delta_{s,L}R^{(s)}_{AB,L}
             \over\omega^2+\Delta_{s,L}^2}.        \tag{GC17}
\]

Thus `2 Delta R` is the exact numerator only in this declared even-Euclidean
quadratic-denominator convention.  It need not tend to a constant for a
genuine massless pole.  FL's exact one-link electric source is the in-parent
counterexample: its positive-frequency weight is proportional to `|k|`, so
`2 Delta R` is proportional to `k^2` although the photon pole is massless.

The native source establishes a **source-visible tensor cone** only if a
degenerate spectral projector or a controlled cluster `S_L` of split finite-
volume polarizations satisfies all of the following on an infrared sequence
`L_ra_r->infinity`.  If simultaneous mathematical refinement is also
claimed, (GC12) is required in addition:

1. **Linear coalescing gaps:** for several inequivalent fixed reciprocal rays
   and every `s in S_L`, `Delta_(s,L)(n)/|k_(L,n)| -> c_T` with one finite
   `c_T>0`, while the cluster splitting is `o(|k|)`.  Exact degeneracy is a
   special case.
2. **Isotropy:** the same `c_T` is obtained on all directional sequences;
   lattice anisotropy vanishes in the principal term.
3. **Raw pole visibility and scaling:** the isolated/projected cluster has a
   nonzero rank-two summed residue at generic finite `k`, and the two
   eigenvalue scalings of

   \[
   \mathcal N^{\mathrm{raw}}_{AB,L}(S_L)
   =\sum_{s\in S_L}2\Delta_{s,L}R^{(s)}_{AB,L}     \tag{GC18}
   \]

   are reported rather than forced to approach a constant.  A vanishing
   power law is allowed; an identically absent or persistently rank-deficient
   complete-source cluster fails source-visible rank two.
4. **Canonical amplitude, if claimed:** a finite positive isotropic limit of
   (GC18) is required only after a separately derived, prospectively fixed
   nonderivative canonical/source normalization (GC16b).  Without that map,
   the raw scaling exponent is the result and absolute pole normalization
   remains open.
5. **Atomic stability:** the candidate does not dissolve into an increasing
   set of finite-volume multiparticle levels.  A below-threshold isolated
   cluster or a proved symmetry-protected threshold delta function is
   required.
6. **Physical Ward ancestry:** the temporal source, native divergence,
   recoil/port terms, and contacts must give the stress identity.  Applying
   (GC13) by hand is not this proof.
7. **Common cone:** using the same clock and length calibration, `c_T` must
   equal the admitted Maxwell/matter probe speed.  In FD's factorized class,
   this is the condition `m_X=0` and common positive
   `kappa_X/chi_X` for every retained sector.
8. **Parent stability:** the pole and residue must survive the controlled
   higher inherited terms or a volume-uniform remainder bound.  A leading-H6
   finite-size pattern alone is not an all-orders phase theorem.

Passing items 1--8 would establish the missing source-visible, Ward-qualified
common massless tensor cone within the declared parent and scaling domain.
It would not by itself establish universal nonlinear gravitational coupling.

The route is refuted in this source channel if the lowest TT gap stays finite,
has `z!=1`, or has direction-dependent limiting velocity; if no isolated or
protected pole cluster survives; if its speed differs from the common probe
speed; or if the physical Ward/contact test fails.  A raw numerator that
vanishes with a controlled momentum power is **not** a refutation.  An
identically zero or persistently rank-deficient overlap refutes only this
complete source channel, not a wholly invisible tensor sector.  A failed
amplitude limit refutes a canonical claim only when (GC16b) was derived and
frozen before seeing the result.

## 6. Theorem-level identifiability boundary

Current inputs cannot decide the tests above:

- FY supplies one 30-cell quotient, not a member of the `L^3`-cell FS family,
  and four strictly positive finite gaps at two source weights.  It supplies
  no volume sequence.
- The graph cover transports local coefficients but neither Hilbert spaces
  nor the new lowest momentum character.
- FD leaves the collective mass, inertia, and stiffness as antecedents.  The
  same family and source kinematics are compatible with both
  `omega^2=r lambda_A3(k)` and
  `omega^2=m^2/chi+r lambda_A3(k)`.
- FL's Gaussian-Maxwell baseline has a spin-one photon and an even composite
  continuum, not an isolated helicity-two atom.  A positive tensor result
  must therefore come from inherited non-Gaussian binding/protection or a
  distinct same-parent collective rank-two sector.
- The public Shannon-et-al. results currently in custody establish flux
  scaling and a static small-wavevector structure factor.  They do not
  provide the volume-resolved even-source TT correlator and Lehmann residue
  matrices required by (GC16)--(GC18).

### Theorem `FFCCB-3` -- gap/residue nonidentifiability

Exact local source ancestry, family small-momentum support, A3 isotropy, and a
rank-two kinematic projector do not determine mass, velocity, atomic residue,
or common-cone status.  No theorem from FS+FD+FV+FY+FL can replace the missing
matched-family response sequence.  This is a physics/data boundary, not a
request for additional record machinery.

## 7. Smallest executable no-lab/public-data spectral screen

The bounded spectral screen has two parts.

1. **Complete the native family source.**  Lift FY's support-resolved diagonal
   histories and folds to the translation-covariant `G_L` cell, using bounded
   local H6 histories rather than a many-body enumeration.  Verify recovery
   of FY on its quotient and FV/FX at zero momentum.  This fixes all
   within-cell form factors for (GC16).
2. **Measure one three-size response packet.**  On `G_5,G_10,G_20`, use a
   sign-free projector/GFMC or equivalent method for the inherited
   pure-kinetic ice Hamiltonian.  Measure the full `2x2` TT correlator and the
   one-link transverse photon correlator at

   \[
   n^{(1)}=(1,0,0),\quad
   n^{(2)}=(1,1,0),\quad
   n^{(3)}=(1,-1,0),                               \tag{GC19}
   \]

   plus conjugate momenta.  These have distinct reciprocal norms and supply
   the first bounded directional/dispersion screen.  Report raw imaginary-
   time covariance, normalization, flux sector, source form factors, fitted
   gaps, residue matrices, continuum spectral weight, and common photon
   speed.

Three sizes cannot prove an asymptotic theorem, but they are the minimum
packet that can expose or falsify a spectral trend.  They are not a decisive
closure packet.  Additional doubled sizes are required if corrections remain
unresolved.  A public release can replace a new simulation only if it
contains (GC16), or raw worldline/projector data with sufficient operator-
estimator information to reconstruct the same insertions and covariance.
Published static aggregates alone cannot do so.

The spectral screen also cannot supply the independent momentum-owning
theory dependency exposed by GB.  Even a perfect linear pole cluster leaves
the stress Ward gate open until the physical parent assigns local spatial
momentum and its exchange with support/recoil/controller/port degrees of
freedom.  That is a separate physical-ownership derivation, not an additional
finite-size observable and not permission to identify the scalar GA charge
current with `T^{0j}`.

## 8. Disposition

The finite-family lane has advanced from one finite momentum to an exact
family-level kinematic theorem:

\[
 \boxed{
 \begin{gathered}
 4L^3\text{ inherited H6 rings and rank-six local source coefficients},\\
 V_L^{\rm aff}\propto(La_*)^3,
 \quad k_{\min}={3\pi\over2La_*},
 \quad \operatorname{rank}\Pi^{TT}=2,\\
 \text{but no inherited gap/residue relation across the covers.}
 \end{gathered}}                                    \tag{GC20}
\]

The next executable work is the native-ledger lift and finite-size spectral
screen; the separate momentum-ownership/Ward dependency remains required.
Neither is a new interaction or a claim of gravity from geometric kinematics.
