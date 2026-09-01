# Native hexagon collective-response and composite-channel theorem

**Short name:** `GL6AS V001`  
**Date:** 2026-08-31  
**Status:** author frozen and sealed; exact replay passes; independent hostile
audit required before promotion  
**Inputs:** frozen, independently audited `GL6AO`, `GL6AP`, and `GL6AQ`
snapshots  
**Claim class:** exact conserved-variable and local-continuity theorem for the
pure degree-two hexagon Hamiltonian; exact representation and observable
selection rules; conditional harmonic and Gaussian spectral formulas; no-go
for an unconditional mode, cone, or tensor/gravity promotion

**Not claimed:** convergence of the microscopic all-orders expansion; a
selected thermodynamic phase; a conventional gauge description; spontaneous
symmetry breaking; a gapless excitation, isolated pole, photon, graviton,
physical momentum or speed, common cone, stress tensor, Ricci/Einstein
response, gravity, or `G`.

## 1. Exact pure-loop model and scope

Remove the scalar terms from sealed GL6AO and declare the pure degree-two
hexagon Hamiltonian

\[
 H_{\rm hex}=-J\sum_c\tau_c,\qquad
 J={63\over8}{h^6\over U_d^5}>0.                            \tag{AS01}
\]

On the declared `Q_4`, this is the complete nonscalar operator through order
six.  On the infinite incidence it is the formal uniformly finite-range
linked interaction of GL6AO.  Treating (AS01) as an exact pure-loop model is a
new declaration for this lane; it is not a claim that the all-orders
Schrieffer--Wolff series truncates.

The finite `Q_4` character set has no literal infrared limit.  Every
statement involving \(\chi\to1\) below is conditional on a
translation-covariant growing-quotient or thermodynamic completion of the
local rule (AS01).

## 2. Exact native collective variables

Write the four parent-cell link occupations as

\[
 e_a(x)=n_{(x,a)}-{1\over2},\qquad a=1,2,3,4.              \tag{AS02}
\]

The degree-two constraint is the incidence equation.  In the GL6AP character
gauge,

\[
 B(\chi)=
 \begin{pmatrix}1&1&1&1\\z_1&z_2&z_3&z_4\end{pmatrix},
 \qquad z_a=\chi(d_a),                                    \tag{AS03}
\]

and every locked fluctuation obeys \(B(\chi)e(\chi)=0\).
At \(\chi=1\), the centered four-port space is

\[
 {\cal T}=\{u\in{\mathbb R}^4:\mathbf1^Tu=0\}\cong T_2.   \tag{AS04}
\]

For \(|\sum_a z_a|<4\), \(\ker B(\chi)\) has rank two.  Near the trivial
character, \(z_a=e^{i\theta_a}\), \(\sum_a\theta_a=0\), and it is the plane

\[
 {\cal T}_{\perp\theta}
 =\{u:\mathbf1^Tu=0,\ \theta^Tu=0\}+O(|\theta|).           \tag{AS05}
\]

Thus the native finite-character variable has two transverse polarizations
inside the moving `T2` bundle.  Its dimension two does not make it the local
pair `E`.

There is a stronger dynamical statement specific to (AS01).  Define the four
port totals on a periodic quotient,

\[
 N_a=\sum_x n_{(x,a)},\qquad
 \widetilde N_a=N_a-\frac14\sum_bN_b.                     \tag{AS06}
\]

Every elementary hexagon uses three distinct port labels, each once on each
alternating half of the cycle.  Its toggle removes and inserts exactly one
link of each used port.  Hence

\[
 [H_{\rm hex},N_a]=0\quad(a=1,\ldots,4),\qquad
 \sum_aN_a=2|{\cal X}|.                                   \tag{AS07}
\]

The fixed `A1` total leaves three independent centered charges transforming
as `T2`.  This conservation is exact for the pure hexagon Hamiltonian.  It
must not be promoted to the unprojected finite-\(h\) Hamiltonian, whose
one-link flips do not conserve (AS06).

The degree-four complement automorphism

\[
 {\cal C}:n_e\mapsto1-n_e                                \tag{AS08}
\]

also preserves the locked space and commutes with every \(\tau_c\).  It sends
\(e\mapsto-e\), while every pair read \(M_{ab}=Z_aZ_b\) and every loop toggle
is even.

## 3. Exact cycle symbol and local continuity

For a hexagon using distinct ports \(a,b,c\), choose either orientation and
discard an irrelevant unit character phase.  Its Fourier boundary column is

\[
 \begin{aligned}
 C^{abc}_a(\chi)&=z_b-z_c,\\
 C^{abc}_b(\chi)&=z_c-z_a,\\
 C^{abc}_c(\chi)&=z_a-z_b,\\
 C^{abc}_d(\chi)&=0\quad(d\notin\{a,b,c\}).                \tag{AS09}
 \end{aligned}
\]

Let \(C(\chi)\) be the four-by-four matrix of the four unordered port
triples.  Direct algebra gives

\[
 \mathbf1^TC(\chi)=0,\qquad z^TC(\chi)=0,\qquad
 \operatorname{im}C(\chi)=\ker B(\chi)                    \tag{AS10}
\]

whenever the character is nontrivial.  The last equality follows because
both spaces have rank two.  At the trivial character \(C(1)=0\), which is the
Fourier statement of (AS07).

Writing \(C(\chi)=iC_1(\theta)+O(|\theta|^2)\) and
\(P_{\cal T}=I-\mathbf1\mathbf1^T/4\), the four exact leading columns obey

\[
 \boxed{
 C_1(\theta)C_1(\theta)^T
 =4\left[I_2(\theta)P_{\cal T}-\theta\theta^T\right],
 \qquad I_2(\theta)=\sum_a\theta_a^2.}                    \tag{AS11}
\]

Therefore, for \(u\in{\cal T}_{\perp\theta}\),

\[
 u^\dagger C_1C_1^Tu=4I_2(\theta)\,u^\dagger u.           \tag{AS12}
\]

Splitting each toggle into its two oriented alternating transitions defines
a Hermitian local cycle-current operator.  The Heisenberg equation then has
the exact Fourier continuity form

\[
 \dot e(\chi)=C(\chi)\,{\cal I}(\chi).                    \tag{AS13}
\]

Equations (AS09)--(AS13) use only incidence and local moves.  They introduce
no gauge variable or gauge phase.

## 4. What conservation proves, and what it does not

Take a translation- and `S4`-invariant ground state of a thermodynamic
completion, normalize a transverse density
\(\rho_u(\chi)=u^\dagger e(\chi)/\sqrt{|{\cal X}|}\), and put
\(t_{\rm hex}=\langle\tau_c\rangle\), common to all elementary hexagons.
The exact ground-state oscillator strength is

\[
 \begin{aligned}
 f_u(\chi)
 &={1\over2}
 \left\langle[\rho_u(\chi)^\dagger,
 [H_{\rm hex},\rho_u(\chi)]]\right\rangle\\
 &={Jt_{\rm hex}\over2}\,
 u^\dagger C(\chi)C(\chi)^\dagger u.                     \tag{AS14}
 \end{aligned}
\]

For a normalized leading transverse \(u\),

\[
 f_u(\chi)=2Jt_{\rm hex} I_2(\theta)+O(|\theta|^3).       \tag{AS15}
\]

Let \(S_u^+(\chi)\) be the strictly positive-frequency density weight and
\(\Delta_{T_2}(\chi;u)\) the least positive support point.  The single-mode
bound gives

\[
 \Delta_{T_2}(\chi;u)
 \le {f_u(\chi)\over S_u^+(\chi)},\qquad S_u^+(\chi)>0.   \tag{AS16}
\]

This is a conditional diagnostic, not a dispersion theorem:

- if \(S_u^+\) approaches a positive constant, (AS16) forces an
  \(O(|\theta|^2)\) or softer excitation;
- if \(S_u^+\sim s_1|\theta|\), it permits the linear bound
  \(O(|\theta|)\);
- if \(S_u^+=O(|\theta|^2)\), it does not force a closing gap;
- no case forces an isolated pole rather than a continuum.

The exact native conclusion within the nontrivial link/pair sectors being
compared is therefore: `T2` is the only charge/continuity-supported soft
**candidate**, while no soft mode is proved.  The local pair `E` has no
corresponding charge, as independently proved by GL6AP.  This statement does
not classify the always-conserved `A1` energy density or possible nonlocal
homology labels.

## 5. Conditional harmonic continuum

A harmonic mode requires extra phase information absent from the sealed
inputs.  Make the following explicit additional premises:

1. a selected translation-, `S4`-, and complement-invariant coherent
   thermodynamic phase;
2. continuous transverse variables \(e\) and conjugates \(a\) with the
   projected canonical bracket;
3. a positive analytic density Hessian \(K_0\) and a positive loop-phase
   Hessian \(G_0\) on the image of \(C\).

The most direct harmonic completion has

\[
 H_{\rm harm}={1\over2}e^\dagger K_0e
 +{1\over2}a^\dagger C(\chi)G_0C(\chi)^\dagger a.         \tag{AS17}
\]

Its two nonzero transverse frequencies satisfy

\[
 \omega_\lambda(\chi)^2
 =\operatorname{eig}_\lambda\!\left[
 C(\chi)G_0C(\chi)^\dagger K_0\right],
 \qquad\lambda=1,2.                                      \tag{AS18}
\]

If the additional isotropic Hessian conditions
\(K_0=\kappa P_{\cal T}\) and \(G_0=gI_4\) hold at the trivial character,
then (AS11) gives

\[
 \boxed{\omega_1^2=\omega_2^2
 =4g\kappa I_2(\theta)+O(|\theta|^3).}                   \tag{AS19}
\]

Thus that specific coherent harmonic completion has two degenerate,
linearly soft **character** modes.  More general positive Hessians can split
their directional character velocities.  A null Hessian, nonanalytic
response, damping, or absence of a coherent conjugate invalidates (AS19).
The coefficient \(g\) includes the phase normalization and any
state-dependent flippability curvature; it cannot be set equal to the bare
\(J\) without an additional rotor/coherent-state construction.

GL6AO supplies no configuration-dependent diagonal term through order six.
It therefore supplies neither the density curvature \(\kappa\) nor a
selected coherent saddle.  These may arise from the state-dependent sector
energy or higher orders, but they cannot be calculated from the sealed pure
operator alone.  Consequently neither (AS19) nor any cone is earned.

Even when (AS19) holds, \(\theta\) is a dimensionless translation character.
A physical velocity or cone additionally needs an embedding, length and time
calibration, a stable scaling limit, and agreement with independently
sourced sectors.

## 6. Exact pair-read overlap and one- versus two-mode channels

For the authenticated GL6AQ pair source, let \(c\in{\cal E}=\ker R\), where
\(R\) is the four-port/six-pair incidence.  Since \(Z_a=-2e_a\),

\[
 \boxed{
 O_x(c)=\sum_{a<b}c_{ab}M_{x,ab}
 =4\sum_{a<b}c_{ab}e_a(x)e_b(x).}                         \tag{AS20}
\]

The authenticated local `E` read is therefore an exact quadratic composite
of the native link variable.  Its state-independent local locked overlap is

\[
 {1\over6}\sum_{k(x)=2}O_x(c)^2={8\over3}\,c^Tc.          \tag{AS21}
\]

Equation (AS21) is an integrated operator norm, not a pole residue.

At the trivial character, `S4` gives
\(\operatorname{Hom}_{S_4}(T_2,E)=0\).  In a complement-invariant state there
is also an all-character parity rule: a one-density `T2` excitation is odd
under \({\cal C}\), while \(O(c)\) is even.  Hence

\[
 \langle {\rm one}\text{-}T_2|O(c)|0\rangle=0             \tag{AS22}
\]

whenever the ground/state vector has definite complement parity.  Without
that state symmetry, the zero-character matrix element still vanishes by
`S4`, but a small-character amplitude proportional to the unique
\(E\)-projection of \(\theta\odot T_2\) is symmetry allowed.

An elementary even `E` excitation or an `E` bound state is not forbidden.
Its existence, energy, and residue are state-dependent.

Under the stronger Gaussian premise that the only elementary soft quanta are
the two transverse `T2` modes, write their one-mode density form factors as
\(r_{a\lambda}(q)\).  Here \(q\) and \(k\) are additive character coordinates,
not calibrated physical momenta.  The exact two-mode form factor inherited
from (AS20) is

\[
 F_c^{\lambda\mu}(q,k-q)
 =4\sum_{a<b}c_{ab}
 \left[
 r_{a\lambda}(q)r_{b\mu}(k-q)
 +r_{b\lambda}(q)r_{a\mu}(k-q)
 \right].                                                 \tag{AS23}
\]

Up to the standard identical-particle counting convention, the Gaussian
positive-frequency measure is

\[
 d\mu_{E,c}^{(2)}(k,\omega)
 ={1\over2}\sum_{\lambda,\mu}\int dq\,
 |F_c^{\lambda\mu}(q,k-q)|^2
 \delta\!\left(\omega-\omega_\lambda(q)
                    -\omega_\mu(k-q)\right)d\omega.       \tag{AS24}
\]

The exact representation reason this channel exists is

\[
 \operatorname{Sym}^2(T_2)=A_1\oplus E\oplus T_2.        \tag{AS25}
\]

For isotropic linear constituents \(\omega=v|q|\), the kinematic two-mode
lower edge is \(v|k|\); for quadratic constituents \(\omega=Dq^2\), it is
\(Dk^2/2\).  These values are actual `E` support thresholds only if the form
factor (AS23) is nonzero arbitrarily near the minimizing momenta; otherwise
the threshold is higher or the channel is absent.  They are continuum
thresholds, not isolated poles.  Interactions can bind an `E` state, but
neither its existence nor its dispersion follows from (AS20)--(AS25).

There is one exact local check that must not be confused with a bulk mode.
For an isolated flippable-hexagon doublet \(|i\rangle,|f\rangle\), choose the
GL6AQ read \(c=\delta M_x\).  Then \(O_f-O_i=16\).  In the isolated
two-state block \(|\pm\rangle=(|i\rangle\pm|f\rangle)/\sqrt2\),

\[
 |\langle+|O_x(\delta M_x)|-\rangle|=8,\qquad
 E_--E_+=2J.                                               \tag{AS26}
\]

The full lattice does not preserve this two-state block, so (AS26) is not a
collective pole or dispersion.

## 7. Retained-lineage overlap

The exact GL6AQ obstructions remain:

\[
 P_QX_eP_Q=0,\qquad
 \operatorname{tr}_0(O_x(c)X_e)=0,                        \tag{AS27}
\]

and a one-cell four-port source is \(A_1\oplus T_2\), with no direct
`E` projection.  At sixth order the retained coefficients enter only
through

\[
 J_c(\kappa)=J\prod_{e\in c}\kappa_e.                     \tag{AS28}
\]

For the multilinear extension about homogeneous \(\kappa=1\), a uniform
port perturbation \(r_a\) changes the four hexagon orientations, indexed by
their missing port \(d\), by

\[
 \delta J_d/J=2\sum_{a\ne d}r_a.                          \tag{AS29}
\]

The map (AS29) has eigenvalue \(6\) on `A1` and \(-2\) on centered `T2`;
it contains no `E`.  The four unoriented hexagon-source components also
transform as \(A_1\oplus T_2\).  Therefore every
translation- and `S4`-invariant state's zero-character cross kernel obeys

\[
 K_{E\leftarrow\mathrm{loop}}(\omega,1)=0.                \tag{AS30}
\]

The scalar source at exactly zero character is even more restricted:
\(\sum_c\tau_c=-H_{\rm hex}/J\).  It has no finite-frequency off-diagonal
matrix elements in the exact energy basis.

Let \(P_E(\theta\odot\ell_{T_2})\) denote the unique `E` Clebsch projection,
and let \(Q_E(\theta)\) be GL6AP's unique quadratic `E` harmonic.  Analytic
covariance permits the leading nonzero small-character structures

\[
 K_{E\leftarrow\mathrm{loop}}(\omega,\theta)\ell
 =\alpha(\omega)P_E(\theta\odot\ell_{T_2})
 +\beta(\omega)Q_E(\theta)\ell_{A_1}
 +O(|\theta|^2\|\ell_{T_2}\|+|\theta|^3|\ell_{A_1}|).     \tag{AS31}
\]

Thus a centered retained-loop source can enter an `E` cross channel at
linear character order, while its scalar part first enters at quadratic
order.  Symmetry allows these coefficients but does not make them nonzero.
In a complement-invariant phase, the even retained-loop source cannot excite
the odd one-density `T2` channel; it can still overlap even composite or
bound channels.

In the conditional coherent harmonic expansion, a complement-even loop
toggle begins with a constant plus a quadratic form in the conjugate cycle
coordinate.  Its leading Fock overlap is therefore a two-mode channel.  The
coefficient is not fixed by AO because no rotor normalization, coherent
saddle, or flippability Hessian was sealed.  A retained-source one-mode pole
would require a separate even elementary excitation or broken complement
selection and remains state-dependent.

In the isolated doublet of (AS26), \(\tau_c\) is diagonal in the
\(|\pm\rangle\) energy basis, so modulation of its own coefficient has zero
transition matrix element between those two states.  Other noncommuting
loops can change that conclusion in the full lattice; their spectral overlap
is state-dependent.

## 8. Composite tensor test

Equation (AS25) gives the exact algebraic opening:

\[
 \operatorname{Sym}^2_0(T_2)=E\oplus T_2.                 \tag{AS32}
\]

This five-component sum is the same finite-group decomposition carried by a
traceless symmetric rank-two object after restriction to the tetrahedral
`S4` symmetry.  Therefore a tensor-like **composite response channel is
algebraically possible** even though no bare `E` charge exists.

The authenticated one-node pair read supplies only its `E` part.  In the
strict degree-two node, the local pair `T2` part vanishes identically and
the `A1` part is fixed.  A nonzero composite `T2` requires a bilocal,
derivative, or coarse-grained product not authenticated by that one-node
read.

The representation coincidence (AS32) proves none of the following:

1. an elementary or isolated composite pole rather than the continuum
   (AS24);
2. degeneracy, common residues, or a common velocity between the `E` and
   composite-`T2` pieces;
3. enhancement from finite `S4` to physical rotations;
4. a physical metric, transverse-traceless helicity reduction, stress
   conservation, or a diffeomorphism Ward identity;
5. a Ricci/Einstein response, gravity, or `G`.

In particular, a linear two-particle threshold inherited from conditional
linear `T2` constituents is not a graviton pole and is not a physical cone.

Finally, GL6AQ's normalized product trace remains an exact universal
counterexample: it can have positive pair-`E` correlation norm while every
retarded commutator vanishes.  None of the exact overlaps above guarantees a
nonzero stationary response in every lawful state.

## 9. Exact conclusion and strict ceiling

\[
 \boxed{\begin{gathered}
 \text{native conserved/continuity candidate}=T_2
 \text{ with two finite-character transverse polarizations};\\
 \text{authenticated local pair read}=E
 \text{ quadratic composite};\\
 \text{no selected state/Hessian}\Longrightarrow
 \text{no derived dispersion or pole};\\
 \operatorname{Sym}^2_0(T_2)=E\oplus T_2
 \text{ allows only an algebraic tensor-like composite channel}.
 \end{gathered}}                                             \tag{AS33}
\]

The strongest earned answer is therefore conditional.  Among the nontrivial
link/pair sectors in scope, the pure native move identifies `T2`, not `E`,
as the only charge-supported soft candidate.
If a coherent positive harmonic completion is independently established, it
has two linearly soft character modes with (AS18), reducing to (AS19) under
isotropic Hessians.  The sealed inputs alone permit a gapped response,
gapless continuum, diffusive/nonanalytic response, isolated bound pole, or no
visible pair-read weight.

Nothing here assumes a conventional gauge phase or derives a photon,
graviton, physical momentum/speed/cone, stress or Ricci/Einstein law, gravity,
or `G`.
