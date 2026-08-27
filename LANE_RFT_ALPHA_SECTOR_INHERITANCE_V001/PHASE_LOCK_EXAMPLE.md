# Exact phase-lock example: one parent, two alpha sectors

**Purpose:** demonstrate one mathematically explicit mechanism by which a parent
can admit several alpha phases, while every record inside one realized phase
inherits one common value

**Status:** exact classical two-derivative EFT construction and static planar
wall theorem; not asserted to be nature's microscopic or quantum parent

## 1. Parent action

In natural units, take one \(U(1)\) connection and an order parameter \(\phi\):

\[
 {\cal L}=
 -{1\over2}(\partial\phi)^2
 -{\lambda\over4}(\phi^2-v^2)^2
 -{1\over4}Z(\phi)F_{\mu\nu}F^{\mu\nu}
 +\sum_a\bar\psi_a
 \left[i\gamma^\mu(\partial_\mu+i n_aqA_\mu)-m_a\right]\psi_a,
                                                               \tag{1}
\]

where \(\lambda,v,Z_0>0\), \(q\ne0\), and the fixed charge lattice is
\(q_a=n_aq\).  Define

\[
 Z(\phi)=Z_0e^{\beta h(\phi/v)},
 \qquad
 h(x)={15x-10x^3+3x^5\over8}.                           \tag{2}
\]

The polynomial obeys

\[
 h(\pm1)=\pm1,
 \qquad h'(\pm1)=h''(\pm1)=0,                           \tag{3}
\]

and the exponential ensures \(Z(\phi)>0\).

## 2. Exact classical homogeneous phases

At \(F_{\mu\nu}=0\) and \(\psi_a=0\), the written classical action has two
degenerate homogeneous vacua

\[
 \phi_\pm=\pm v,
 \qquad Z_\pm=Z_0e^{\pm\beta}.                          \tag{4}
\]

Canonical normalization in either vacuum gives

\[
 \alpha_\pm={q^2\over4\pi Z_\pm}
 =\alpha_0e^{\mp\beta},
 \qquad \alpha_0={q^2\over4\pi Z_0}.                   \tag{5}
\]

Any prescribed positive pair \((a_+,a_-)\) can be represented algebraically by

\[
 \alpha_0=\sqrt{a_+a_-},
 \qquad
 \beta={1\over2}\log{a_-\over a_+}.                    \tag{6}
\]

For the displayed sign convention the scalar equation is

\[
 \Box\phi-V'(\phi)-{1\over4}Z'(\phi)F_{\mu\nu}F^{\mu\nu}=0. \tag{7}
\]

Because \(V'(\pm v)=Z'(\pm v)=0\), \(\phi=\pm v\) solves this scalar equation
for any finite electromagnetic field.  It is a solution of the full coupled
system only when \(A_\mu\) and \(\psi_a\) also obey their equations.

Writing \(\phi=\pm v+\eta\), the scalar quadratic block is

\[
 {\cal L}^{(2)}_\eta
 =-{1\over2}(\partial\eta)^2
 -{1\over2}(2\lambda v^2)\eta^2.                        \tag{8}
\]

Since \(Z'(\pm v)=Z''(\pm v)=0\), there is no
\(F_{\rm bg}\eta\), \(F_{\rm bg}\eta^2\), or
\(\eta\,\delta F\) term from the gauge-kinetic function.  This proves local
scalar stability in the written classical action, not nonlinear/global
stability of an arbitrary electromagnetic background.

Every charged subsystem constituted as a restriction of one homogeneous phase
therefore uses its phase's common classical base coupling.  A common quantum RG
trajectory additionally requires SAI4–SAI6, including gauge-preserving
quantization and threshold matching.  The model contains no record-specific
alpha coordinate.

## 3. An exact interface cost

For a static planar configuration with vanishing background \(F_{\mu\nu}\) that
connects \(-v\) and \(+v\), the minimum wall is

\[
 \phi_{\rm wall}(z)=
 v\tanh\!\left[v\sqrt{\lambda\over2}(z-z_0)\right].      \tag{9}
\]

Its tension is

\[
 \begin{aligned}
 \sigma
 &=\int_{-\infty}^{\infty}
 \left[{1\over2}(\phi')^2+V(\phi)\right]dz\\
 &=\int_{-v}^{v}\sqrt{2V(\phi)}\,d\phi
 ={2\sqrt{2\lambda}\over3}v^3>0.                       \tag{10}
 \end{aligned}
\]

The lower bound follows from completing the square:

\[
 {1\over2}(\phi')^2+V
 ={1\over2}\left(\phi'-\sqrt{2V}\right)^2
 +\phi'\sqrt{2V},                                      \tag{11}
\]

for a monotone \(-v\) to \(+v\) transition.  The kink (9) saturates it.
For a one-dimensional planar path with host vacuum at both
\(z\to\pm\infty\) that actually reaches the opposite vacuum, the total-variation
form

\[
 {1\over2}(\phi')^2+V(\phi)
 \ge |\phi'|\sqrt{2V(\phi)}                              \tag{12}
\]

gives a lower bound \(2\sigma\) per unit transverse area for the two complete
traversals.  A generic finite kink–antikink excursion need not reach the opposite
vacuum and can have a smaller interacting energy.  No universal area-only or
nucleation bound is claimed for arbitrary curved inclusions.

Thus this parent makes the sharp alternatives explicit:

- within the specified classical planar class, a conserved scalar/interface
  energy budget below \(2\sigma\) per unit area excludes a complete excursion to
  the opposite vacuum and back;
- a foreign-alpha region can exist only as the other \(\phi\) phase together
  with a physical interpolation region; or
- a separate gauge connection can be added, which is another sector rather than
  an ordinary host subsystem.

## 4. What this example proves and does not prove

It proves by construction that:

1. one parent action can admit multiple internally consistent alpha phases;
2. selecting one phase makes one alpha common to its same-sector records; and
3. a complete static planar transition between the two exact vacua has the
   positive classical tension (10).

It does not prove that nature uses (1), that either phase contains realistic
atoms or records for arbitrary parameters, or that record formation initially
selects the phase.  Charged-field and photon corrections can shift the quantum
effective potential and gauge function; exact quantum vacua would require
renormalization conditions or a symmetry preserving
\(V'_{\rm eff}(\pm v)=Z'_{\rm eff}(\pm v)=0\).

Finite quantum systems can tunnel or nucleate interior bubbles.  A fixed exterior
boundary does not forbid that.  Exact superselection requires an imposed sector
restriction, an exact conserved label, or an appropriate infinite-volume
disjoint-representation limit.  Conditional on a selected homogeneous phase,
the phase fixes the common base coefficient; no selection or permanent lock is
derived.

The quintic \(h\) is unbounded, so for \(\beta\ne0\) the gauge coefficient tends
toward zero in one large-field direction and infinity in the other.  This is
harmless on the kink interval \([-v,v]\), but the example does not establish a
globally weakly coupled or UV-complete parent.  A bounded smooth extension with
the same endpoint values and first two endpoint derivatives would be needed for
that stronger claim.
