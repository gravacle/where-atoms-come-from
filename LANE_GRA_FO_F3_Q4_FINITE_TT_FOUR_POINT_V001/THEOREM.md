# Finite-volume inherited F3/q4 TT composite-cumulant screen

**Lane ID:** `GRA-FO-F3-Q4-FTT4-V001`

**Short name:** `FTT4`

**Date:** 2026-08-27

**Claim class:** exact plaquette-complete periodic diamond quotient; exact
finite ice-sector construction; deterministic complete finite-matrix
diagonalization of the inherited sixth-order ring Hamiltonian; exact
complement and spatial-tensor typing; finite-volume static connected two- and
four-**composite** cumulants; scalar composite-source Legendre transform;
finite-pole and two-one-link proxy screen; sharp four-one-link/amputation and
thermodynamic boundary

**Status:**
`INHERITED_H6_ONLY__180_STATE_TRANSLATION_CLOSED_RING_SECTOR__STATIC_TWO_Q_AND_FOUR_Q_CONNECTED_CUMULANTS_NONZERO__SCALAR_COMPOSITE_LEGENDRE_QUARTIC_POSITIVE_IN_BOTH_SUSCEPTIBILITY_EIGENCHANNELS__SELECTED_PLUS_COS_HAS_FOUR_FINITE_POLES__MOMENTUM1_POLE_ABOVE_FINITE_TWO_LINK_THRESHOLD_PROXY__NO_BELOW_PROXY_OR_ENERGY_EXCLUSIVE_TENSOR_CANDIDATE_IDENTIFIED__CONNECTED_FOUR_ONE_LINK_CHANNEL_2PI_AND_THERMODYNAMIC_HELICITY2_OPEN`

**Not claimed:** that the selected ring component is the complete periodic
ice Hilbert space or its global ground sector; that finite discrete levels
are a continuum; a photon-amputated four-link vertex; an attractive TT
kernel; a bound state; a massless pole; helicity two; a rank-two Ward
identity; universal stress coupling; RGRL-B; gravity; or `G`.

## 1. Exact question and dependency custody

`FM` identified the normalized connected four-one-link response of the
already inherited pure-ice Hamiltonian as the next no-laboratory
discriminator.  This packet computes a bounded composite-source precursor on
a lawful finite periodic diamond quotient without adding any interaction; it
does not complete the four-one-link calculation.  The microscopic Hamiltonian
is exactly

\[
 H_6=-J_6\sum_C B_C,
 \qquad J_6={63h^6\over8U_d^5}>0,                 \tag{FO01}
\]

where the sum is over every simple six-cycle of the supplied
plaquette-complete quotient and `B_C` flips its two alternating ice
orientations.  There is no pair attraction, tensor projector, fitted
coupling, diagonal flippability term, or order-eight rescue term.

The final load-bearing dependency bytes are:

| dependency | SHA-256 |
|---|---|
| `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md` | `5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe` |
| `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md` | `a91caa20d16b0a1194333f9b51d96546a4ea24d55e23bf1f04c7d249641af8db` |
| `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md` | `cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98` |
| `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md` | `c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4` |
| `LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md` | `98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452` |
| `LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/INDEPENDENT_AUDIT.md` | `327bf6a4476c4c6382757dc156a96c6032233d34c25c1f7935e2582acf6c607a` |
| `LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md` | `78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25` |
| `LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md` | `53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9` |

The final `FM` hostile correction is therefore inherited: a bare finite-range
vertex may be analytic, but a dressed massless 1PI function need not be, and a
finite-order Hamiltonian may bind or shift a pole when solved or resummed
nonperturbatively.

## 2. A small plaquette-complete periodic diamond quotient

Let the Bravais-cell group be `Z_30`.  Put one `A_x` and one `B_x` diamond
vertex in every cell and connect

\[
 e=(x,a): A_x\longrightarrow B_{x+s_a},
 \qquad (s_0,s_1,s_2,s_3)=(0,1,5,19)\pmod {30}.  \tag{FO02}
\]

The four shifts are the q4 append labels.  The twenty unordered sums of
three shifts with repetition are all distinct modulo 30.  This `B_3`
property implies the lower `B_2` property and prevents a nonlocal quotient
identification from making a new four- or six-cycle.  Direct enumeration then
gives

\[
 |V|=60,qquad |E|=120,qquad |\mathcal C_6|=120. \tag{FO03}
\]

Every simple six-cycle uses exactly three q4 labels twice and is one of the
four elementary diamond hexagon types per cell.  Thus this quotient satisfies
the `CW` elementary-plaquette-completeness condition.  It is a compact
certified witness, not a claim that it is the globally smallest possible
diamond supercell under every boundary convention.

Begin with the frozen degree-two state occupying labels `0,1` in every cell
and reverse the alternating noncontractible octagon

\[
 (84,11,9,114,112,39,37,87).                     \tag{FO04}
\]

Closure under all inherited hexagon flips produces one exact connected ring
sector with

\[
 |\Omega_{\rm sec}|=180,qquad
 N_{\rm transitions}=420,qquad
 4\le d_{\rm ring}\le6.                          \tag{FO05}
\]

Every state obeys two-in/two-out ice.  Translation `x -> x+1` preserves the
sector and commutes with (FO01).  Global occupation complement maps it to a
disjoint isomorphic 180-state sector.  That separation is retained rather
than hidden: the 180-state sector is **not** the full periodic ice Hilbert
space.

The translation action is free on this component: it decomposes into exactly
six orbits of length 30.  Fourier reduction therefore gives thirty Hermitian
`6 x 6` blocks whose spectra reproduce the full `180 x 180` spectrum.  The
zero-momentum block has exact spectrum

\[
 \{-2-2\sqrt2,0,0,-2+2\sqrt2,2,2\},
\]

and the `m=5` block has exact spectrum

\[
 \{-1-2\sqrt3,-1,1,1,1,-1+2\sqrt3\}.
\]

The conjugate `m=25` block is identical in spectrum.  A complete scan of the
remaining blocks places their lowest eigenvalue at or above
`-4.410987667370205`, separated from `-1-2 sqrt(3)` by more than `0.053`.
This independently certifies which algebraic block values supply the sector
ground energy and first gap; it does not make them the global ground and gap
of all periodic ice sectors.

In units `J_6=1`, deterministic complete diagonalization of the exact
`180 x 180` matrix gives a unique sector ground state

\[
 {E_0\over J_6}=-2(1+\sqrt2),                    \tag{FO06}
\]

and its lowest sector gap is

\[
 {\Delta_{\rm sec}\over J_6}
 =1+2\sqrt2-2\sqrt3=0.364325509608438\ldots .   \tag{FO07}
\]

These algebraic values are checksums of the finite matrix, not continuum
dispersion claims.

## 3. Correct complement-even `E` tensor and TT source

Use the tetrahedral q4 bond vectors

\[
 u_0={1\over\sqrt3}(1,1,1),\quad
 u_1={1\over\sqrt3}(1,-1,-1),\quad
 u_2={1\over\sqrt3}(-1,1,-1),\quad
 u_3={1\over\sqrt3}(-1,-1,1).                   \tag{FO08}
\]

At either sublattice vertex let `s_(v,a)=1-2n_(v,a)`.  Reversing all local
bond orientations on the `B` sublattice changes the sign of the vector below
but not its tensor.  Define

\[
 F_v=\sum_as_{v,a}u_a,
 \qquad
 Q_{ij}(v)=F_iF_j-{\delta_{ij}\over3}F^2.         \tag{FO09}
\]

Because `sum_a u_a=0` and `s_a^2=1`, the one-link square contribution to
(FO09) is isotropic and is removed by the trace subtraction.  Hence

\[
 Q(v)=\operatorname{STF}\!\left[
 \sum_{a<b}s_{v,a}s_{v,b}
       (u_au_b^T+u_bu_a^T)\right],               \tag{FO10}
\]

which is exactly an even pair observable.  On two-in/two-out ice it lies in
the surviving local `E` sector of `FK`.  It is symmetric traceless and
unchanged by global complement `s -> -s`.

The cyclic quotient is the image of the primitive-cell map
`(x_1,x_2,x_3) -> x_1+5x_2+19x_3 mod 30`.  Select cyclic momentum `m=1` and
the shortest reciprocal alias

\[
 q={1\over30}(1,5,-11),
 \qquad k=2\pi A^{-T}q,                           \tag{FO11}
\]

where the columns of `A` are `u_1-u_0,u_2-u_0,u_3-u_0`.  With
`P^T=I-\hat k\hat k^T`, choose orthonormal transverse vectors and the usual
two tensors

\[
 \epsilon_+={e_1e_1^T-e_2e_2^T\over\sqrt2},
 \qquad
 \epsilon_\times={e_1e_2^T+e_2e_1^T\over\sqrt2}. \tag{FO12}
\]

The basis offset of `B` is `c=(-1/4,-1/4,-1/4)` in primitive coordinates.
The complex source is

\[
 \mathcal O_\lambda(m)={1\over\sqrt{|V|}}
 \left[
 \sum_xe^{2\pi imx/30}\epsilon_\lambda:Q(A_x)
 +e^{2\pi iq\cdot c}
 \sum_xe^{2\pi imx/30}\epsilon_\lambda:Q(B_x)
 \right].                                       \tag{FO13}
\]

It transforms with momentum `+1`, is complement even, and is explicitly TT.
Its Hermitian cosine and sine parts for `lambda=+,x` give four real source
channels.  The polarization tensors have unit Frobenius norm and are mutually
orthogonal.  The `1/sqrt(|V|)` Fourier normalization and the `sqrt(2)` real-
quadrature normalization are fixed throughout this packet.  They are a finite
quotient convention, not a derived continuum field residue, cell-volume
normalization, or common photon normalization.  Here “TT” means the exact
kinematic transverse-traceless projection at the selected nonzero reciprocal
vector; it does not assign a helicity-two particle representation.

## 4. Exact finite static connected two-composite matrix

Let `R=Q_0(H_6-E_0)^(-1)Q_0`, where `Q_0` removes the unique sector ground
state.  In the ordered basis

\[
 (O_{+,c},O_{+,s},O_{\times,c},O_{\times,s}),
\]

the zero-frequency connected response

\[
 W^{(2)}_{ab}=2\langle0|O_aRO_b|0\rangle         \tag{FO14}
\]

is the static, time-integrated connected susceptibility of two **composite**
`Q` insertions.  Since each `Q` is bilinear in one-link variables, this is a
four-one-link expectation before the one-link disconnected/Wick channels are
subtracted.  It is not yet the connected four-one-link function or its
photon-amputated vertex.

Restoring `J_6`, its matrix is

\[
 W^{(2)}={1\over J_6}
 \begin{pmatrix}
 1.130847135996&0&-0.037434360320&0\\
 0&1.130847135996&0&-0.037434360320\\
 -0.037434360320&0&0.114433012322&0\\
 0&-0.037434360320&0&0.114433012322
 \end{pmatrix}.                                  \tag{FO15}
\]

Its eigenvalues are two translation-quadrature doublets,

\[
 {0.113056176225351\over J_6}\quad(2\times),
 \qquad
 {1.132223972092660\over J_6}\quad(2\times).   \tag{FO16}
\]

Thus the selected finite even/TT composite response is nonzero and has full rank in
the two TT polarizations and their two real momentum quadratures.  The small
polarization mixing is a property of this skew finite quotient, not a
continuum anisotropy theorem.

## 5. Static connected four-composite cumulant and scalar Legendre vertex

For one real composite source `O`, define

\[
 W(\lambda)=-E_{\rm g}(H_6-\lambda O)
 =W(0)+{W_2\over2!}\lambda^2
       +{W_3\over3!}\lambda^3
       +{W_4\over4!}\lambda^4+\cdots .           \tag{FO17}
\]

These derivatives are zero-frequency, time-integrated connected cumulants of
the composite observable, not equal-time moments.  For the plus-cosine
channel, momentum conservation gives `W_3=0`.  Rayleigh--Schrodinger recursion
in the complete 180-state matrix, independently replayed with both spectral
resolvents and an augmented linear solve, gives

\[
 \boxed{
 W_2={1.130847135995723\over J_6},\qquad
 W_4={-0.136825085605100\over J_6^3}.}            \tag{FO18}
\]

The connected subtraction is explicit.  With

\[
 A=\langle ORO\rangle,
 \quad C=\langle OR^2O\rangle,
 \quad B=\langle ORORORO\rangle,
\]

one finds `E_4=-B+AC=0.005701045233546/J_6^3` and
`W_4=-24E_4`.  The scalar Legendre transform with respect to the **composite
source** therefore has

\[
 \boxed{
 \Gamma^{(4)}_{\rm comp}=-{W_4\over W_2^4}
 =0.083666214307836J_6>0.}                        \tag{FO19}
\]

In the two cosine-polarization eigenchannels of (FO15), the corresponding
checks are

| channel | `W_2 J_6` | `W_4 J_6^3` | `Gamma_comp^(4)/J_6` |
|---|---:|---:|---:|
| low susceptibility | 0.113056176225351 | -0.003036590242056 | 18.586988116578 |
| high susceptibility | 1.132223972092660 | -0.148796476796671 | 0.090544748133189 |

Thus the inherited ring model has a nonzero connected finite-volume
**four-`Q` composite cumulant** without an added interaction.  In this static
one-coordinate Legendre convention the two susceptibility eigenchannels
sampled have a positive composite quartic coefficient.  This sign is not
promoted to an all-direction, all-frequency attractive or repulsive photon
kernel.

The distinction is physical, not terminological.  A quadratic composite of
an exactly Gaussian underlying field can itself have nonzero higher connected
cumulants and a nonquadratic composite Legendre action.  Therefore the nonzero
`W_4` and `Gamma_comp^(4)` do not by themselves isolate the inherited compact
interaction from Gaussian composite kinematics, prove a nonzero four-one-link
1PI kernel, or diagnose binding.

The distinction in the subscript is load-bearing: (FO19) is **not the
connected or photon-amputated four-one-link vertex** in `FM17`.  A scalar
Legendre transform of a bilinear composite source amputates only that scalar
composite susceptibility.  A four-photon vertex instead requires four
separately normalized one-link external legs, a relative-momentum assignment
inside each pair, subtraction of the one-link disconnected terms followed by
the channels reducible in the selected two-particle channel, and a common
finite-volume residue convention.  “Composite 1PI” below refers only to this
one-source Legendre coordinate and is never used as a synonym for photon-leg
1PI or channel-2PI.

## 6. Finite spectral poles and the below-proxy diagnostic

Above the stated numerical residue floor, the selected plus-cosine spectral
measure has four pole groups; an independent four-step Lanczos closure
reproduces the same support:

| `Delta/J_6` | spectral weight |
|---:|---:|
| 3.194109035554332 | 0.005026104004432 |
| 3.490165912028476 | 1.965864248197576 |
| 6.166688337463908 | 0.003649484732984 |
| 9.139267639373482 | 0.000000908886868 |

They reconstruct (FO18) through
`W_2=2 sum_n |<n|O|0>|^2/Delta_n`.  These are ordinary discrete finite-sector
poles, not a massless branch or continuum.

As a conservative comparison, define a finite two-one-link threshold proxy
by extracting the lowest transverse one-link pole at every nonzero cyclic
momentum and minimizing

\[
 \Delta_{2\gamma}^{\rm proxy}(m)
 =\min_{p\ne0,m}\left[\Delta_\gamma(p)
 +\Delta_\gamma(m-p)\right].                    \tag{FO20}
\]

For total momentum `m=1`, the minimum is attained at `p=9,22` and equals

\[
 \Delta_{2\gamma}^{\rm proxy}(1)
 =2.059674505691458J_6.                           \tag{FO21}
\]

The lowest selected TT pole is

\[
 3.194109035554332J_6>
 2.059674505691458J_6.                            \tag{FO22}
\]

Moreover, that same energy occurs with nonzero transverse one-link residue.
After adjoining the disjoint complement sector, complement-even and
complement-odd combinations are exactly degenerate because the two blocks do
not mix.  The finite spectrum therefore supplies no energy-exclusive level
that this diagnostic can identify uniquely as a tensor composite.  The
correct result is:

\[
 \boxed{\text{no below-proxy or energy-exclusive tensor candidate is identified at }m=1.} \tag{FO23}
\]

Equation (FO20) is a finite noninteracting threshold proxy, not a
thermodynamic two-photon/ring continuum.  Being above it fails this declared
finite screen, but it is not a no-bound-state theorem: interactions,
finite-volume shifts, other ring sectors, and the thermodynamic threshold are
not determined by the proxy.  Energy degeneracy with one-link response is
also not a spin decomposition; after complement-sector doubling the even and
odd combinations remain distinct degenerate parity states.

## 7. Sharp remaining four-one-link datum

The finite calculation proves three narrower pieces:

1. the unchanged inherited ring Hamiltonian has nonzero static connected
   two-`Q` and four-`Q` composite cumulants;
2. its scalar composite-source Legendre quartic is nonzero, with positive
   coefficients along both sampled `W_2` polarization eigenchannels; and
3. this selected diagnostic has only finite poles and no identified
   energy-exclusive candidate below its declared two-one-link proxy.

These are precursors to, not evaluations of, the four-one-link quantity `FM`
left open.  The sharp next lawful datum is the frequency- and
relative-momentum-resolved connected function of four independently sourced
transverse one-link observables,

\[
 \boxed{
 W^{(4),c}_{\lambda_1\lambda_2\lambda_3\lambda_4}
 (\omega_1,p_1;\ldots;\omega_4,p_4),
 \quad \sum_i(\omega_i,p_i)=0,}                  \tag{FO24}
\]

computed from a four-source generating functional with every one-link
disconnected contraction removed.  Its matching prescription must freeze

\[
 \boxed{((\omega_i,p_i,\lambda_i)_{i=1}^4;\
 G^{(2)}_\gamma;\ Z_\gamma;\ \mathcal S_{\rm 2PR};\
 V_{\rm cell},\mathcal N_\gamma)},               \tag{FO25}
\]

where `G_gamma^(2)` and `Z_gamma` own the external one-link propagator and
residue conventions, `S_2PR` declares the subtraction of contributions
reducible in the selected two-particle channel, and the final entries own the
cell/field normalization.  External-leg amputation of (FO24), followed by the
channel-2PI subtraction, yields the finite kernel that can enter a
Bethe--Salpeter calculation.  A total-momentum bilinear source contains none
of the required relative-momentum information.  Equations (FO24)--(FO25) are
a matching specification, not a new interaction, and can be evaluated in a
successor calculation without a laboratory.

Even after that calculation, a thermodynamic helicity-two claim requires
sector- and volume-scaling, an isolated pole and residue that survive at fixed
physical momentum, a common linear cone, a helicity-two transformation test,
and a protecting rank-two Ward identity.  None follows from one 180-state
sector.

## 8. Disposition

The strongest earned chain is

\[
 \boxed{
 \begin{gathered}
 \text{unchanged F3/q4 parent}
 \longrightarrow H_6=-J_6\sum_CB_C,\\
 \text{exact periodic ice sector}+\text{TT pair source}
 \longrightarrow W_{2Q}\ne0,\ W_{4Q}\ne0,\\
 \text{scalar composite Legendre transform}
 \longrightarrow \Gamma^{(4)}_{\rm comp}>0
 \text{ in both susceptibility eigenchannels},\\
 \text{finite spectrum}
 \longrightarrow \text{no below-proxy or energy-exclusive candidate at }m=1,\\
 \text{finite-volume diagnostic}
 \not\Longrightarrow \text{four-one-link 1PI/2PI, thermodynamic helicity-two, or gravity}.
 \end{gathered}}                                  \tag{FO26}
\]

This advances the physics lane by supplying an exact finite-sector composite
response and a negative below-proxy screen from unchanged inherited ice
dynamics.  It does not authenticate the ice occupation as a record, isolate a
fundamental non-Gaussian four-one-link interaction from composite kinematics,
or demonstrate an isolated tensor carrier.
