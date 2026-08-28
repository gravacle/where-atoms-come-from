# Inherited F3/q4 even-channel kernel and order-eight pole boundary

**Lane ID:** `GRA-FM-F3-Q4-IEKB-V001`

**Short name:** `IEKB`

**Date:** 2026-08-27

**Claim class:** exact fixed-support Feshbach operator classification through
order eight at symmetric detuning; exact alternating-octagon coefficient;
exact order-eight diagonal-scalar theorem; exact distinction between the
leading compact interaction and its Gaussian infrared fixed point; exact
strictly-truncated response boundary for the bare finite-range order-eight
insertion; sharp continuum matching target

**Status:**
`ORDER8_OPERATOR_TOPOLOGIES_CLASSIFIED__EXACT_J8_429_OVER_16__COMPLETE_ORDER8_DIAGONAL_SCALAR_AND_V8_ZERO__HEXAGON_DRESSING_PLUS_NEW_OCTAGON_RING_ONLY__LEADING_H6_COMPACT_RING_ALREADY_INTERACTING__ORDER8_NOT_FIRST_NON_GAUSSIANITY__BARE_TT_VERTEX_ALLOWED_BUT_NORMALIZED_CONNECTED_FOUR_POINT_UNEVALUATED__STRICT_SINGLE_INSERTION_DOES_NOT_ESTABLISH_NEW_TENSOR_POLE__NONPERTURBATIVE_RESUMMATION_AND_WARD_RESIDUE_TESTS_REMAIN`

**Not claimed:** a volume-uniform all-orders expansion; a numerical reduction
of every dressed-hexagon coefficient; a matched continuum four-photon
amplitude; a tensor bound state; a rank-two Ward identity; a protected
helicity-two pole; support selection; visible electromagnetism; RGRL-B;
gravity; or Newton's constant.

## 1. Exact question and dependency custody

`FL` proved that the physical one-link field carries the spin-one Maxwell pole
while the local centered even pair carries a two-photon continuum at the
Gaussian Maxwell fixed point. It correctly nominated the inherited
non-Gaussian even-channel kernel as the next direct composite test. This
packet asks the first order-by-order question without adding an interaction:

> What does the unchanged F3 parent generate at order `h^8`, which part can
> enter a continuum even/TT two-photon kernel, and can that finite-order term
> itself produce an isolated tensor pole?

One terminology correction is load-bearing. Order eight is the first
microscopic correction **beyond the leading projected hexagon Hamiltonian**.
It is **not the first inherited non-Gaussian interaction**. The exact
order-six hard-core, flippability-projected compact ring Hamiltonian is already
an interacting microscopic many-body model. `MAXWELL-IR` is its Gaussian
long-distance fixed-point description, not an assertion that the complete
order-six lattice Hamiltonian is quadratic at every scale.

The load-bearing files are:

| dependency | SHA-256 |
|---|---|
| `LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md` | `00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec` |
| `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md` | `5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe` |
| `LANE_CROSS_ALPHA_GRA_F3_DIAMOND_SIXTH_ORDER_V001/THEOREM.md` | `211b1aa61917c98dccae278129a8016a1a14f73587908bfeceeba090a808536c` |
| `LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md` | `05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769` |
| `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md` | `cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98` |
| `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md` | `c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4` |
| `LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md` | `98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452` |

No statement below modifies those parents. No `j-j` attraction, tensor
projector, continuum graviton field, fitted coefficient, or rescue term is
inserted.

## 2. Frozen parent and irreducible word convention

On a supplied finite simple coordination-four bipartite support `G=(V,E)` of
girth at least six, freeze the `d_*=2`, `E_R=0` slice of the declared F3
parent,

\[
 H=H_0+V_X,
 \qquad H_0=U_d\sum_v(d_v-2)^2,
 \qquad V_X=-h\sum_{e\in E}X_e,
 \qquad U_d>0.                                      \tag{FM01}
\]

Let `P_2` project onto the ice manifold `Omega_2(G)`. For one ice state `n`
and a set `S` of links, put

\[
 \sigma_e(n)=1-2n_e,
 \qquad
 \Delta_n(S)=U_d\sum_v
 \left(\sum_{\substack{e\in S\\e\ni v}}\sigma_e(n)\right)^2. \tag{FM02}
\]

For a flip word `w=(e_1,...,e_m)`, its prefix parity is

\[
 S_r(w)=\{e_1\}\triangle\cdots\triangle\{e_r\}.   \tag{FM03}
\]

The irreducible Feshbach contribution retains exactly the words for which
every proper prefix lies in `Q_2=1-P_2`. At symmetric detuning every such
even-length word contributes

\[
 -h^m\prod_{r=1}^{m-1}{1\over\Delta_n(S_r(w))}.    \tag{FM04}
\]

The sign follows from `m` factors of `-h` and `m-1` negative resolvents. This
is the same self-consistent Feshbach convention used by the sixth-order
parent. The global occupation-parity unitary makes the projected expansion
even in `h`.

## 3. Complete order-eight transition topology

The symmetric difference of two degree-two configurations is an even-degree
subgraph and therefore a disjoint union of alternating even cycles. Girth six
then gives an exact order-eight classification.

### Theorem `IEKB-1` -- word and endpoint classification

Every nonzero projected order-eight matrix element has exactly one of the
following endpoint types:

1. the same ice state (diagonal);
2. an ice state differing by one alternating length-six cycle, with exactly
   one link flipped two additional times; or
3. an ice state differing by one alternating simple length-eight cycle, with
   all eight cycle links flipped once.

There is no other endpoint topology. In particular, two disjoint cycles need
at least twelve odd link occurrences, while a two- or four-link cycle is
forbidden.

#### Proof

The number of odd multiplicities in an eight-letter word is the size of its
endpoint symmetric difference. An off-diagonal projected endpoint requires
that odd set to be a union of alternating cycles. With girth six and only
eight letters, its size is six or eight. Six odd multiplicities use six
letters, so the remaining two occurrences repeat one link; that link may be
on or off the hexagon. Eight odd multiplicities force eight distinct links,
which form one simple octagon. The diagonal partitions are precisely

\[
 (8),\quad(6,2),\quad(4,4),\quad(4,2,2),\quad(2,2,2,2). \tag{FM05}
\]

The `(8)` class necessarily revisits `P_2` after two flips and is folded.
QED.

Thus the complete new order-eight operator content consists of:

- scalar diagonal energy;
- dressed length-six transitions; and
- new alternating length-eight ring transitions.

The phrase “dressed length-six” permits the exact finite Feshbach coefficient
to depend on the flippable hexagon and its allowed local environment. It is
not silently replaced by one translation-invariant number at a boundary or
on a nonhomogeneous supplied graph.

## 4. Exact alternating-octagon coefficient

Let `C_8=(e_1,...,e_8)` be one simple alternating octagon in an ice state.
Every nonempty proper subset of its edges is a forest, hence violates the ice
constraint. All `8!=40320` flip orders therefore remain in `Q_2` at every
proper prefix. Define

\[
 J_8=h^8\sum_{\pi\in S_8}\prod_{r=1}^{7}
 {1\over\Delta_n(S_r(\pi))}.                       \tag{FM06}
\]

At `E_R=0`, the denominator coefficients in units of `U_d` have thirteen
classes:

| sorted coefficients | multiplicity |
|---|---:|
| `(2,2,2,2,2,2,2)` | 512 |
| `(2,2,2,2,2,2,4)` | 1280 |
| `(2,2,2,2,2,4,4)` | 2816 |
| `(2,2,2,2,4,4,4)` | 4672 |
| `(2,2,2,2,4,4,6)` | 1152 |
| `(2,2,2,4,4,4,4)` | 5632 |
| `(2,2,2,4,4,4,6)` | 3456 |
| `(2,2,2,4,4,6,6)` | 2304 |
| `(2,2,4,4,4,4,4)` | 4096 |
| `(2,2,4,4,4,4,6)` | 4608 |
| `(2,2,4,4,4,6,6)` | 5184 |
| `(2,2,4,4,6,6,6)` | 3456 |
| `(2,2,4,4,6,6,8)` | 1152 |

Their multiplicities sum to `40320`, and exact summation gives

\[
 \boxed{J_8={429h^8\over16U_d^7}>0.}               \tag{FM07}
\]

No lower-order fold can have an octagon endpoint: all lower off-diagonal
kernels have a hexagon endpoint, while the lower diagonal kernels are scalar.
Equation (FM07) is therefore the complete coefficient at the first appearance
of this operator, not merely an irreducible fragment. If the supplied graph
contains no length-eight cycle, its sum is empty. A literal periodic diamond
witness with an alternating octagon is checked executably; the theorem does
not assume every boundary completion contains the same octagon count.

With

\[
 B_C=P_C\prod_{e\in C}X_e,                         \tag{FM08}
\]

the new topology is

\[
 H_{8,\mathrm{oct}}=-{429h^8\over16U_d^7}
 \sum_{C\in\mathcal C_8(G)}B_C.                   \tag{FM09}
\]

## 5. Complete order-eight diagonal remains scalar

An irreducible diagonal eight-word touches at most four distinct links. Its
denominators depend only on:

1. the occupied/empty color of each selected link; and
2. which selected link pairs share a vertex.

No untoggled link enters (FM02).

### Lemma `IEKB-L2` -- colored four-edge census

On a finite simple `z`-regular bipartite graph of girth at least six, every
fixed-degree-`d` edge configuration has the same census of occupied/empty
colored incidence types on at most four selected links.

#### Proof

The one-, two-, and three-edge statements are the frozen fixed-degree census
of `CW`. For four selected links, the absence of triangles and squares means
their original-graph incidence subgraph is a forest: a four-star, a fork, a
path, or a disjoint union of smaller trees. A color-preserving embedding of
each connected tree is counted by choosing its first colored edge and then
adding leaves. At every used endpoint the number of available edges of a
nominated color is `d` or `z-d`, minus the already used edges of that color.
These numbers are configuration independent. A collision during a four-edge
tree embedding would make a cycle of length at most four and is forbidden.
Counts of the disconnected incidence types follow successively by subtracting
the already fixed types from the fixed total number of four-link color
combinations. Thus every exact colored incidence type has a fixed count. QED.

### Theorem `IEKB-2` -- no order-eight potential

The complete order-eight diagonal is scalar on `Omega_2(G)`:

\[
 \boxed{H^{(8)}_{\rm diag}=\epsilon_8(G,U_d)P_2,
 \qquad V_8=0.}                                    \tag{FM10}
\]

#### Proof

The irreducible `(6,2)`, `(4,4)`, `(4,2,2)`, and `(2,2,2,2)` sums are
functions only of the colored four-edge census and are scalar by `IEKB-L2`.
The `(8)` histories are folds. Through fourth order every lower kernel and
its energy derivatives are scalar. At sixth order the only non-scalar lower
kernel is the alternating-hexagon kinetic operator. Consequently every
order-eight self-consistency fold is either scalar or a renormalization of an
existing hexagon transition; none is a configuration-dependent diagonal.
QED.

`V_8=0` means that no standard diagonal flippability potential or other
configuration-dependent diagonal operator is generated through order eight.
It does not say that the off-diagonal hexagon dressing is zero.

In the fixed Feshbach basis the projected operator can therefore be written
without inventing a term as

\[
 \boxed{
 H_{\rm eff}=E_{\rm scalar}^{(\le8)}(h)P_2
 -J_6\sum_{C_6}B_{C_6}
 +h^8R^{(8)}_6
 -{429h^8\over16U_d^7}\sum_{C_8}B_{C_8}
 +O_G(h^{10}),}                                    \tag{FM11}
\]

where `R_6^(8)` has matrix elements only between configurations already
related by one alternating hexagon. Its entries are exact finite sums of
eight-word denominators plus the stated folds. Equation (FM11) classifies
the complete operator support without claiming a volume-uniform remainder or
silently imposing translation invariance on `R_6^(8)`.

## 6. Where the first non-Gaussian even kernel actually enters

The leading projected Hamiltonian

\[
 H^{(6)}_{\rm ice}=-J_6\sum_{C_6}B_{C_6},
 \qquad J_6={63h^6\over8U_d^5},                    \tag{FM12}
\]

is a constrained compact ring-exchange model. It is not a free oscillator
Hamiltonian. A harmonic long-wavelength expansion and renormalization can
flow to `MAXWELL-IR`, but the omitted compactness, hard-core constraint, and
higher-field terms already supply the first allowed even multiphoton
interaction. Thus the continuum amputated even-channel kernel has the
schematic matching structure

\[
 {\cal K}_{TT}^{\rm match}
 =J_6\,{\cal C}^{(6)}_{TT}
 +h^8\left({\cal C}^{(8,6)}_{TT}
 +{429\over16U_d^7}{\cal C}^{(8,8)}_{TT}\right)
 +O_G(h^{10}),                                     \tag{FM13}
\]

after the powers of `h` already included in `J_6` are accounted for. The
symbols `C_TT` are continuum matching functionals of the actual loop
operators, field residue, geometry, and source normalization; they are not
new interactions or fitted numbers.

Equation (FM13) makes the scope precise:

- the order-six pure-ice model owns the leading possible non-Gaussian kernel;
- order eight adds a dressed-hexagon correction and an exactly normalized new
  octagon operator;
- `MAXWELL-IR` alone fixes neither matching functional; and
- two-point flux scaling and a photon pole do not determine a connected
  amputated four-point function.

The exact local ice projection still applies. The on-site centered pair source
contains only the even `E` sector; its local pair `T2` is killed. A continuum
TT source must therefore be built from the correctly glued nonlocal
field-strength bilinear and then projected with

\[
 \Pi^{TT}_{ij,kl}(\mathbf k)
 ={1\over2}\left(P^T_{ik}P^T_{jl}+P^T_{il}P^T_{jk}
                 -P^T_{ij}P^T_{kl}\right).        \tag{FM14}
\]

The allowed order-eight loop vertices can have a nonzero even/TT projection;
symmetry does not force it to vanish. Because each **bare** loop insertion has
finite lattice support, its tree-level Fourier vertex is analytic in external
lattice momenta. This analyticity statement is deliberately not promoted to
the fully dressed massless 1PI four-point function: photon loops may produce
threshold nonanalyticities or logarithms. Nor is a bare vertex by itself a
new propagator denominator.

## 7. Exact strict-truncation tensor-pole boundary

Let `Pi_2` be the Gaussian two-photon bubble in the even/TT channel and let
`K_TT` be the amputated kernel that is two-particle irreducible in that chosen
two-photon channel. Perturbatively,

\[
 \chi_{TT}=\Pi_2+\Pi_2{\cal K}_{TT}\Pi_2
 +O({\cal K}_{TT}^2).                              \tag{FM15}
\]

An isolated composite pole would instead require a zero of the resummed
Bethe--Salpeter denominator,

\[
 \chi_{TT}^{\rm resum}
 =\Pi_2\left(1-{\cal K}_{TT}\Pi_2\right)^{-1},
 \qquad
 \det(1-{\cal K}_{TT}\Pi_2)=0.                    \tag{FM16}
\]

### Theorem `IEKB-3` -- a strict single insertion does not establish a new pole

The inherited order-eight operators can change Maxwell stiffnesses and supply
an analytic bare even/TT two-photon vertex. The susceptibility truncated
strictly at the single-insertion expression (FM15) **does not establish a new
isolated tensor pole**: it inherits the singular support of its constituent
bubbles and contains no Bethe--Salpeter denominator of the form (FM16).

This is not a universal assertion that a Hamiltonian term carrying a finite
power of `h` can never create or move a pole. If the order-six-plus-order-eight
Hamiltonian is diagonalized nonperturbatively, or if its kernel/self-energy is
iterated through a Bethe--Salpeter/Dyson equation, it can in principle bind a
state or shift an already existing pole. A finite-order term could also expose
an already present exchange pole in a different model. None of those
resummations or pre-existing tensor denominators is contained in the strict
single-insertion calculation here. The exact result is therefore a boundary
on what (FM15) proves, not an all-settings no-pole theorem.

For two massless photons, the continuum at fixed total momentum begins at
`omega=c|k|`; at zero total momentum its lower edge is zero. A hypothetical
slower collective branch could lie below the continuum for nonzero `k`, but
its existence would already be a nonperturbative spectral result and it would
not satisfy the required common cone until its velocity was shown to equal
`c`. A genuinely massless common-cone pole lies at the continuum edge and
therefore requires a protecting constraint or Ward identity and a nonzero
thermodynamic residue. No rank-two gauge constraint or Ward identity is
generated by the order-eight operator classification.

This is stronger than saying that the coefficient is small, but narrower than
a nonperturbative no-pole theorem. An actual pole test requires the full
Bethe--Salpeter resummation or a nonperturbative finite-volume spectral
calculation followed by residue scaling.

## 8. Sharp derivability and identifiability boundary

The lattice-side order-eight boundary is not blocked:

- the new ring topology is derived;
- its coefficient (FM07) is exact;
- the full diagonal is proved scalar; and
- all remaining order-eight matrix elements are exactly typed as hexagon
  dressing.

What is not contained in the owned data is the continuum, commonly normalized
connected TT four-point function of the interacting pure-ice model. The first
required observable is

\[
 G^{(4),c}_{TT}:=
 \Pi^{TT}\,\langle {\cal O}{\cal O}{\cal O}{\cal O}\rangle_c\,\Pi^{TT},
                                                               \tag{FM17}
\]

with the cell volume, q4 field/source normalization, external-state
convention, and matching scale frozen. After external photon legs are
amputated and the pieces reducible in the selected two-photon channel are
separated, the corresponding first 2PI matching coefficient may be written

\[
 \boxed{
 g_{TT}^{\rm match}(\mu_{\rm RG})
 :=\mathcal N_{TT}^{-1}
 \Pi^{TT}{\cal K}^{\mathrm{2PI}}_{\rm ice,TT}({p_i};\mu)
 \Pi^{TT}\Big|_{\text{declared symmetric low-momentum point}}.} \tag{FM18}
\]

with `mu_RG` denoting the common matching scale (not Shannon's QDM potential
parameter), the external photon residues amputated, and the q4
cell/field/source normalization held fixed. Its leading contribution belongs to the order-six
pure-ice Hamiltonian; order eight supplies a calculable correction once the
same matching prescription is frozen. Existing public phase data establish
two-point Maxwell behavior, not (FM17), and none of the finite local response
packets evaluates it.

The connected function (FM17), its amputated vertex, and the two-particle-
irreducible kernel are distinct objects; the subtraction needed to pass from
one to the next may not be skipped. One scalar value of (FM18) is only the
first discriminator. If it is nonzero,
the pole question needs the momentum-dependent kernel
`K_TT(omega,k;p,p')`, the resummation (FM16), a TT residue, finite-volume
survival, a common linear cone, and a rank-two Ward/constraint test. No number
is inserted or fit here.

This target is accessible without a laboratory: exact diagonalization,
sign-free projector methods where applicable, or a controlled linked-cluster
calculation can evaluate the connected four-point function and TT finite-size
level shifts of the already fixed pure-ice Hamiltonian. That is physics work,
not model expansion.

## 9. Disposition

The strongest earned chain is

\[
 \boxed{
 \begin{gathered}
 \text{unchanged F3 single-link parent}+d_*=2
 \longrightarrow -J_6\sum B_{C_6}\quad\text{(already interacting)},\\
 O(h^8)\longrightarrow
 \text{hexagon dressing}
 -{429h^8\over16U_d^7}\sum B_{C_8},
 \qquad V_8=0,\\
 \text{even/TT projection}
 \longrightarrow\text{analytic bare vertex; dressed kernel unevaluated},\\
 \text{strict single-insertion truncation}
 \not\Longrightarrow\text{new isolated helicity-two pole}.
 \end{gathered}}                                   \tag{FM19}
\]

No new microscopic interaction has been added. The next direct no-lab physics
calculation is the matched connected four-point/TT finite-volume response of
the existing order-six pure-ice model, with the exact order-eight octagon and
hexagon-dressing corrections added only after that baseline is measured. A
negative pole/residue result would close the direct two-photon composite route
and leave the distinct same-parent rank-two Ward architecture as the lawful
gravity route.
