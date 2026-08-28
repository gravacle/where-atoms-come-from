# Gaussian-Maxwell infrared pole screen for q4/F3 link-pair observables

**Lane ID:** `GRA-FL-F3-Q4-MCPS-V001`

**Short name:** `MCPS`

**Date:** 2026-08-27

**Claim class:** exact ice-projected parity and finite-group typing; exact
conditional Gaussian-Maxwell spectral theorem; exact exclusion of a direct
one-particle helicity-two pole from the FJ link-pair composite at the Maxwell
fixed point; bounded next-lawful-route theorem

**Status:**
`ONE_LINK_ODD_FLUX_HAS_SPIN1_MAXWELL_POLE__CENTERED_PAIR_EVEN_COMPOSITE_HAS_TWO_PHOTON_CONTINUUM_AND_CONTACTS__PAIR_A1_CONSTANT_AND_PAIR_T2_KILLED_BY_ICE__NO_ISOLATED_HELICITY2_POLE__S4_TENSOR_COUNTING_NOT_CONTINUUM_SPIN__DIRECT_COMPOSITE_ROUTE_FAILS_AT_GAUSSIAN_FIXED_POINT__INHERITED_NON_GAUSSIAN_KERNEL_OR_DISTINCT_RANK2_WARD_ROUTE_NEXT`

**Not claimed:** that the imported `U(1)` phase is an all-orders phase of the
complete F3 parent; that the emergent photon is visible electromagnetism; that
FJ's response-state preparation or PMMDC solder is complete; that every
microscopic spectral contribution is Gaussian; that an independently owned
scalar mode is present; a protected tensor pole, RGRL-B, universal stress
coupling, Einstein dynamics, gravity, or numerical `G`.

## 1. Exact internal input and imported infrared input

FJ supplies four physical active link observables at an authenticated,
programmed q4/F3 diamond site,

\[
 s_a=Z_{e_a},\qquad j_{ab}=s_as_b,\qquad 1\le a<b\le4,       \tag{FL01}
\]

and proves that the unprojected six `j_ab` span the tetrahedral edge module
`A1+E+T2`.  CROSS-CW independently proves that, on supplied
plaquette-complete coordination-four diamond support, the leading controlled
`d_*=2` F3 Hamiltonian is exactly the pure-kinetic quantum-ice model

\[
 H_{\rm ice}^{(6)}=E_0-J_6\sum_C B_C,
 \qquad J_6={63h^6\over8U_d^5}>0,                           \tag{FL02}
\]

up to `O(h^8)` on each fixed graph.

The thermodynamic input is deliberately separate:

**`MAXWELL-IR`.**  The zero-flux, complement-symmetric thermodynamic phase of
the pure-kinetic model in (FL02) is governed in the infrared by a deconfined
Gaussian compact-`U(1)` Maxwell fixed point with nonzero electric and magnetic
stiffnesses and one linearly dispersing transverse spin-one photon.

This input is supported by the zero-temperature GFMC/ED flux scaling and
small-wave-vector structure factor reported by Shannon et al.,
“Quantum Ice: a quantum Monte Carlo study,” *Physical Review Letters* 108,
067204 (2012), [arXiv:1105.4196](https://arxiv.org/abs/1105.4196),
[DOI](https://doi.org/10.1103/PhysRevLett.108.067204).  It is not an internal
finite-F3 theorem.  The primary-source custody is arXiv v3: Eq. (1) gives the
Maxwell action, Eq. (2) gives the pure-ice Hamiltonian, Fig. 5(b) tests
`mu=0` flux scaling, and Fig. 6 tests the small-wave-vector structure factor.
The exact statements below are conditional consequences of combining the
internal operator map with `MAXWELL-IR`.

The question is narrower than whether the local operator list can be arranged
as six symmetric-tensor coordinates:

> Which one-particle poles, multiparticle continua, and contact terms do the
> actual one-link and pair observables carry at the Gaussian Maxwell fixed
> point?

The load-bearing dependency custody is frozen to the following exact files:

| dependency | SHA-256 |
|---|---|
| `LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md` | `05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769` |
| `LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/INDEPENDENT_AUDIT.md` | `44690ad431c85af7a4947a431a4c57ad4cd8b19a346e3d26720b341f77256f90` |
| `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md` | `5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe` |
| `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/PRIMARY_SOURCES.md` | `4ee84b4f9b78003cdc5ce80a86cba6cbab618feb1fcd78d25903cb5e97c42a62` |

The separately owned local ice-sector theorem is a consistency cross-check,
not a logical premise of the pole classification, and its dependency hash
will be pinned only after its hostile audit is final.  Equations
(FL03)-(FL09) rederive exactly the minimal identities this screen needs, so
the two packets do not form a circular proof chain.

## 2. Minimal exact ice projection and complement parity

For every vertex `v`, let `eta_v=+1` on `V_+` and `eta_v=-1` on
`V_-`, and let the four vectors `e_{v,a}` point outward from `v`.  At a
reference `V_+` vertex choose the tetrahedral bond frame

\[
 \begin{aligned}
 \mathbf e_1&=(1,1,1)/\sqrt3,&
 \mathbf e_2&=(1,-1,-1)/\sqrt3,\\
 \mathbf e_3&=(-1,1,-1)/\sqrt3,&
 \mathbf e_4&=(-1,-1,1)/\sqrt3.
 \end{aligned}                                                \tag{FL03}
\]

The outward frame at either sublattice obeys

\[
 \sum_a\mathbf e_a=0,
 \qquad \mathbf e_a\!\cdot\!\mathbf e_b=
 \begin{cases}1&a=b,\\-1/3&a\ne b,\end{cases}
 \qquad \sum_a e_a^ie_a^j={4\over3}\delta^{ij}.              \tag{FL04}
\]

The FJ variable `s_a(v)=Z_{e_a(v)}` has the same sign at the two endpoints
of a shared link, so it is not by itself an outward-oriented flux.  Define
the incidence-corrected variable

\[
 \varepsilon_{v,a}:=\eta_v s_a(v).
\]

It changes sign between the endpoints of a shared link.  On the exact local
ice fiber,

\[
 \Omega_{2,v}=\{\varepsilon\in\{-1,+1\}^4:
                 \sum_a\varepsilon_{v,a}=0\},               \tag{FL05}
\]

define the local flux vector

\[
 \mathcal E_v={3\over4}\sum_a\mathbf e_{v,a}
                              \varepsilon_{v,a}.              \tag{FL06}
\]

Equations (FL04)-(FL05) give the exact inverse

\[
 \boxed{\varepsilon_{v,a}=\mathbf e_{v,a}\cdot\mathcal E_v.}\tag{FL07}
\]

This is the microscopic lattice-electric-flux map.  On a shared edge
`e_{w,a'}=-e_{v,a}` and `varepsilon_{w,a'}=-varepsilon_{v,a}`, so the vector
`varepsilon e` agrees at both endpoints.  Equation (FL05) is the zero-charge
lattice Gauss law with correct incidence signs.  Since `eta_v` is fixed,
`s_a` and `varepsilon_{v,a}` have the same local representation and complement
parity, while `j_ab=s_as_b=varepsilon_{v,a}varepsilon_{v,b}`.  Passing from
this exact lattice identity to a smooth Gaussian field is the imported
infrared step.  Calling this divergence-free link variable electric flux is
the lattice-Hamiltonian convention used here; Shannon et al. call the same
ice-arrow field magnetic flux.  Source-free Gaussian Maxwell duality exchanges
those names without changing any pole, parity, or continuum classification
below.  This convention choice is not an identification with visible
electromagnetism.

The ice constraint also gives

\[
 j_{ab}=j_{cd}\quad
 (\{a,b,c,d\}=\{1,2,3,4\}),
 \qquad \sum_{a<b}j_{ab}=-2I.                                \tag{FL08}
\]

Therefore the pair `T2` opposite-edge differences vanish, the pair `A1` is a
fixed constant, and the centered pair-state variations span only `E`.
This is the minimum local projection needed for the pole screen; the fuller
hybrid tensor-rank construction is a separate lane.

At symmetric detuning the exact global complement

\[
 \mathcal C=\prod_eX_e,
 \qquad \mathcal Cs_a\mathcal C^\dagger=-s_a,
 \qquad \mathcal Cj_{ab}\mathcal C^\dagger=j_{ab}             \tag{FL09}
\]

commutes with the pure-kinetic ice Hamiltonian.  In its zero-flux symmetric
ground state, one-link operators are complement odd and pair operators are
complement even.

### Theorem `MCPS-1` - exact microscopic sector boundary

Inside `d_*=2` ice:

1. the nonconstant one-link sector is a three-dimensional `T2` lattice-vector
   module and is complement odd;
2. the nonconstant connected pair sector is the two-dimensional `E` module and is
   complement even;
3. pair `A1` contributes no connected response because it is constant; and
4. pair `T2` is identically zero after projection.

No local six-pair pole analysis which carries the unprojected `A1+E+T2`
independence through the ice projection is correctly typed.

## 3. The one-link observable carries the Maxwell photon pole

Under `MAXWELL-IR`, choose transverse gauge and write the quadratic action as

\[
 S_M={1\over2}\int dt\,d^3x\,
 \left[\chi\,\dot{\mathbf A}_T^2
       -\kappa(\nabla\!\times\!\mathbf A_T)^2\right],
 \qquad c^2={\kappa\over\chi},
 \qquad \nabla\cdot\mathbf A_T=0.                            \tag{FL10}
\]

The coarse link flux has the operator expansion

\[
 s_a(v)=\eta_v Z_1\,e_{v,a}^iE_i(\mathbf x_v)
 +\hbox{derivative, higher-odd, and noninfrared terms},       \tag{FL11}
\]

with nonzero overlap `Z_1` fixed up to the continuum field normalization.
Let

\[
 P^T_{ij}(\mathbf k)=\delta_{ij}-{k_ik_j\over|\mathbf k|^2}. \tag{FL12}
\]

Put `d_{v,a}=eta_v e_{v,a}`, the fixed globally oriented link form factor.
Separating equal-time/UV contact terms from the propagating part, the leading
retarded response is

\[
 \boxed{
 \chi^R_{s_{v,a}s_{w,b}}(\omega,\mathbf k)
 =Z_E\,d_{v,a}^iP^T_{ij}(\mathbf k)d_{w,b}^j
 {c^2|\mathbf k|^2
  \over(\omega+i0)^2-c^2|\mathbf k|^2}
 +C_{ab}(\omega,\mathbf k)+\chi^R_{\rm higher},}             \tag{FL13}
\]

where `Z_E>0`, `C` is analytic/contact at the light cone, and the higher odd
terms begin with multiphoton or irrelevant contributions.  The equivalent
form with `omega^2` in the numerator differs by an analytic contact term.

For every nonzero momentum, the residue matrix

\[
 d_{v,a}^iP^T_{ij}d_{w,b}^j                              \tag{FL14}
\]

has rank two on the four-link space.  It is the two-polarization transverse
**vector** residue.  Individual link form factors can vanish for special
orientations, but the complete one-link `T2` sector contains the two photon
polarizations.  The poles are

\[
 \omega=\pm c|\mathbf k|.                                  \tag{FL15}
\]

### Theorem `MCPS-2` - one-link pole classification

Conditional on `MAXWELL-IR`, the FJ one-link observables contain an isolated
massless spin-one Maxwell photon pole with helicities `+1,-1` in the emergent
isotropic continuum.  Equation (FL14) is the normalized polarization/form-
factor matrix; the positive-frequency complex pole residue of (FL13) carries
the additional scalar `Z_E c|k|/2`.  That gauge-invariant electric-field pole
weight therefore vanishes linearly as `|k| -> 0`, but is nonzero for generic
nonzero momentum.  It is not a helicity-two pole.

## 4. The centered pair observable has a two-photon continuum

Let `q_A` be either independent centered pair direction in the ice-surviving
`E` sector.  Complement parity gives the exact selection rule

\[
 \boxed{\langle0|q_A|\gamma(\mathbf k,\lambda)\rangle=0}       \tag{FL16}
\]

for every one-photon state in the complement-symmetric zero-flux phase:
the vacuum and `q_A` are even, while a one-photon state created by the odd
Maxwell field is odd.  Thus microscopic renormalization cannot generate a
single-photon pole in this channel without breaking the stated symmetry or
expanding about a nonzero flux background.

At the Gaussian fixed point the leading nonconstant local expansion is a sum
of normal-ordered gauge-invariant bilinears,

\[
 q_A=C_A^{ij}:E_iE_j:
     +D_A^{ij}:B_iB_j:
     +F_A^{ij}:E_iB_j:
     +\hbox{derivatives and higher-even operators}.           \tag{FL17}
\]

Some coefficients can vanish by additional microscopic symmetries; that only
removes spectral weight and cannot create a one-particle pole.  Wick's theorem
gives the nonanalytic connected spectral density in the form

\[
 \rho_{AB}(\omega,\mathbf k)
 =\sum_{\lambda\lambda'}\int{d^3p\over(2\pi)^3}
 {\mathcal M_A^{\lambda\lambda'}
  \mathcal M_B^{\lambda\lambda' *}
  \over4\omega_{\mathbf p}\omega_{\mathbf k-\mathbf p}}
 \delta\!\left(\omega-\omega_{\mathbf p}
                       -\omega_{\mathbf k-\mathbf p}\right)
 -(\omega\to-\omega),                                     \tag{FL18}
\]

with `omega_p=c|p|`.  Here `mathcal M` is the unnormalized
field-strength/polarization vertex; operator coefficients, the
identical-particle convention, and the overall spectral-density convention
are absorbed into it, while the displayed `1/(4 omega_p omega_q)` is the
canonical two-leg normalization.  By the triangle inequality,

\[
 \omega_{\mathbf p}+\omega_{\mathbf k-\mathbf p}
 \ge c|\mathbf k|,                                          \tag{FL19}
\]

and varying `p` supplies continuous support above that threshold.  The
light cone is therefore the edge of a two-particle branch cut, not a delta
function.  The threshold onset and its coefficient depend on the chosen
field-strength bilinear, tensor component, and normalization; no universal
threshold exponent is asserted or needed for the no-pole result.

Local composite renormalization and equal-time commutators add polynomials or
other analytic contact terms in `(omega,k)`.  Microscopic gapped modes can add
higher thresholds.  Neither contribution is a massless propagating pole.

### Theorem `MCPS-3` - pair pole classification

At the complement-symmetric Gaussian Maxwell fixed point, the connected FJ
pair response contains:

- no one-photon pole, by (FL16);
- a two-photon continuum beginning at `|omega|=c|k|`, when the allowed
  bilinear overlap in (FL17) is nonzero;
- analytic/contact terms and possible noninfrared gapped weight; and
- no isolated `1/[(omega+i0)^2-c_2^2k^2]` helicity-two pole.

In a supplied nonzero-flux background, an even bilinear can acquire a term
linear in the fluctuation and hence borrow the same spin-one photon pole.
That would be a background-induced vector pole, not a new tensor particle,
and is outside the zero-flux `MAXWELL-IR` premise.

## 5. Why the local `S4` tensor count does not change the spectrum

For polar `O(3)` tensors restricted to the full tetrahedral point group
`T_d \simeq S_4`, the convention called `T2` occurs both as the
three-dimensional polar-vector representation and as one part of a restricted
even rank-two representation:

\[
 \ell=1\downarrow S_4=T_2,
 \qquad
 \ell=2\downarrow S_4=E\oplus T_2.                          \tag{FL20}
\]

The proper-rotation subgroup alone is `A_4`; the `T2` name in (FL20) uses the
full `T_d` polar action, including its improper classes.  Finite-group labels
therefore do not determine continuum spin.  In this
specific operator realization, the leading infrared `T2` is linear electric
flux and carries the spin-one pole, whereas the even `E` pair is bilinear and
carries the two-photon continuum.  An `S4`-equivariant local linear map which
places these five directions into entries of a symmetric `3 x 3` array is a
kinematic coordinate isomorphism.  It does not change complement parity,
scaling dimension, particle number, Ward identity, or pole denominator.

An isolated helicity-two excitation would require, for a symmetric-tensor
probe `Q_ij`, a term of the form

\[
 \boxed{
 \chi^R_{ij,kl}(\omega,\mathbf k)
 \supset {Z_2\,\Pi^{TT}_{ij,kl}(\mathbf k)
 \over(\omega+i0)^2-c_2^2|\mathbf k|^2},
 \qquad Z_2>0,}                                             \tag{FL21}
\]

where

\[
 \Pi^{TT}_{ij,kl}
 ={1\over2}\left(P^T_{ik}P^T_{jl}+P^T_{il}P^T_{jk}
                  -P^T_{ij}P^T_{kl}\right)                 \tag{FL22}
\]

is the transverse-traceless rank-two projector.  The vector residue (FL14)
is not (FL22), and the convolution (FL18) has no denominator of the form
(FL21).

### Theorem `MCPS-4` - direct composite helicity-two no-go

Under the exact ice map, exact complement symmetry, and `MAXWELL-IR`, no
linear combination of the FJ one-link and pair observables obtains an
isolated massless helicity-two pole merely because its local `S4` directions
can be arranged as `A1+E+T2`.  Odd/even symmetry block-diagonalizes the
one-photon and pair-composite channels.  The earned massless pole is spin one;
the even tensor-looking spectral weight is a two-photon continuum.

This is a no-go for the **direct Gaussian composite route**, not a no-go for
all possible collective rank-two phases of the complete F3 parent.

## 6. Exact versus imported custody

| statement | custody | status |
|---|---|---|
| physical `s_a,j_ab` factors on supplied authenticated q4/F3 support | FJ finite theorem | exact internal |
| `d_*=2` ice constraint and arrow/Gauss map | CROSS-CW plus (FL03)-(FL08) | exact internal on supplied support |
| microscopic complement parity and Hamiltonian commutation | CROSS-CW symmetry plus FJ operators | exact internal at symmetric detuning |
| one-photon parity selection rule | (FL09), (FL16), and complement-symmetric `MAXWELL-IR` state | exact conditional deduction |
| leading pure-kinetic ring Hamiltonian and `J_6` | CROSS-CW | exact internal to sixth order on each admitted fixed graph |
| thermodynamic `mu=0` point is in a `U(1)` liquid | Shannon et al. GFMC/ED | imported numerical phase evidence |
| Gaussian Maxwell action, emergent isotropy, and photon description | `MAXWELL-IR` | imported infrared phase identification |
| vector pole and two-photon convolution once `MAXWELL-IR` is assumed | Gaussian spectral calculus | exact conditional deduction |
| complete F3 all-orders phase and volume-uniform stability | none | open |
| visible-EM identity, helicity two, gravity, or `G` | none | open |

The distinction is load-bearing.  The internal model map determines which
microscopic operator enters which symmetry channel.  Public QMC supports the
phase in which Gaussian spectral calculus is applicable.  Neither alone
proves the complete F3 parent has that thermodynamic limit.

## 7. Next lawful route

The direct pair-composite route has failed at the Gaussian fixed point.  The
next no-lab calculation must test dynamics rather than add more local tensor
labels:

1. derive the inherited `O(h^8)` and higher quasi-local operators of the same
   F3 parent on plaquette-complete diamond support, with a volume-uniform
   remainder/stability bound;
2. project their **even-channel** one-particle-irreducible kernel onto the
   continuum spin-two/TT source and solve the corresponding spectral or
   Bethe-Salpeter problem;
3. require either an isolated branch below the two-photon threshold or a
   symmetry-protected delta function at threshold, with nonzero thermodynamic
   residue, the projector (FL22), a linear common cone, and the necessary
   tensor Ward/constraint identities; and
4. reject the composite route if finite-size residue dissolves into the
   two-photon continuum or the TT/Ward tests fail.

At the Gaussian point the irreducible binding kernel is zero, so this route
starts with a negative baseline.  No `j-j` attraction or fitted tensor
interaction may be inserted to force a pole.  Only interactions generated by
the already declared F3 parent may be used.

If that inherited non-Gaussian calculation does not produce the protected
pole, the lawful gravity route is a distinct same-parent collective
rank-two architecture: derive a symmetric gauge variable, its rank-two
Gauss/diffeomorphism-type constraint and Ward identity, and its universal
stress vertex from record-lineage relational dynamics.  The Maxwell photon
and pair continuum may remain constituents or environmental channels, but
they may not be renamed as the graviton.

## 8. Disposition

The strongest earned infrared chain is

\[
 \boxed{
 \begin{gathered}
 \text{authenticated q4/F3 links}+d_*=2
 \longrightarrow\text{odd lattice electric flux }s_a
 \longrightarrow\text{spin-one Maxwell photon pole},\\
 \text{centered even pair }q_A
 \longrightarrow\text{two-photon continuum}+\text{contacts},\\
 \text{local }S_4\text{ tensor count}
 \not\Longrightarrow\text{isolated helicity-two pole}.
 \end{gathered}}                                             \tag{FL23}
\]

This advances the gravity lane by eliminating a tempting false closure and
by converting the remaining question into one precise spectral calculation
on inherited non-Gaussian dynamics.  It does not justify more record machinery
or a new microscopic rescue interaction.
