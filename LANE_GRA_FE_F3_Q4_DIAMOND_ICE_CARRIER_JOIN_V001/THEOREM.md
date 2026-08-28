# q4 append-incidence to diamond-ice carrier join theorem

**Lane ID:** `GRA-FE-F3-Q4-DICJ-V001`

**Short name:** `Q4DICJ`

**Date:** 2026-08-27

**Claim class:** exact infinite-graph identification; exact finite-slab local
exhaustion theorem; exact conditional q4-carrier/F3-edge support join; exact
inheritance ledger for the existing leading diamond-ice Hamiltonian and
spin-one U(1) comparator

**Builder status:** `READY_FOR_INDEPENDENT_HOSTILE_AUDIT`

**Not claimed:** that mutually exclusive q4 front labels are already
coexisting physical sites; that append incidence is already a material F3
link; an unbounded physical q4 stream; a support-binding or support-stability
law; a physically selected periodic completion; a volume-uniform all-orders
F3 phase theorem; visible electromagnetism, alpha, a tensor mode, gravity, or
`G`

## 1. Frozen inputs and the narrow question

The bounded q4 record-stream theorem supplies, below its mission cap, the
fronts

\[
 S_N=\{m\in\mathbb N_0^4:\mathbf1^{\mathsf T}m=N\}
 \tag{FE01}
\]

and the append incidences

\[
 m\longrightarrow m+e_a,
 \qquad a\in\{1,2,3,4\}.
 \tag{FE02}
\]

`CCMAC` has already proved the exact finite incidence identity
`B_N^dagger B_N=4I+A_N` and isolated `Q4-CARRIER-LIFT` as a physical, not
combinatorial, antecedent.  Independently, `CROSS-CW-F3-PKU1S` proves that the
F3 parent on supplied coordination-four diamond support at `d_*=2` has the
leading pure diamond quantum-ice Hamiltonian, while `CROSS-ALPHA-GRA-F3-D6KA`
and `CROSS-ALPHA-GRA-F3-CST` supply its exact short-loop and stiffness
audits.

The present question is exactly:

> Is the support graph between two consecutive q4 count fronts already the
> diamond net in its deep interior, and if so, which parts of the existing F3
> diamond-ice result can be inherited without silently turning alternative
> record histories into coexisting matter?

The graph answer is exact and positive.  The physical promotion remains an
explicit antecedent.

## 2. The infinite two-front completion

For one integer `N`, define the two affine integer hyperplanes

\[
 \widetilde S_N=\{m\in\mathbb Z^4:\mathbf1^{\mathsf T}m=N\},
 \qquad
 \widetilde S_{N+1}=\{c\in\mathbb Z^4:\mathbf1^{\mathsf T}c=N+1\}.
 \tag{FE03}
\]

Let `D_infinity(N)` be the bipartite graph with these two sets as its vertex
classes and with

\[
 \{m,c\}\in E(D_\infty(N))
 \quad\Longleftrightarrow\quad
 c=m+e_a\text{ for exactly one }a\in\{1,2,3,4\}.
 \tag{FE04}
\]

The word “completion” is mathematical: negative count coordinates are used
only to expose the translation-invariant bulk graph.  They are not physical
negative record counts.  Section 4 proves that arbitrarily large rooted balls
of this graph occur unchanged inside the nonnegative fronts (FE01).

Put

\[
 n_1={1\over\sqrt3}(1,1,1),\quad
 n_2={1\over\sqrt3}(1,-1,-1),\quad
 n_3={1\over\sqrt3}(-1,1,-1),\quad
 n_4={1\over\sqrt3}(-1,-1,1).
 \tag{FE05}
\]

Then

\[
 \sum_an_a=0,
 \qquad
 n_a\cdot n_b=
 \begin{cases}1,&a=b,\\-1/3,&a\ne b.\end{cases}
 \tag{FE06}
\]

For a declared bond scale `a_*>0`, embed both parts by

\[
 X(m)=a_*\sum_am_an_a.
 \tag{FE07}
\]

The four append bonds from every `N`-part vertex are exactly the four regular
tetrahedral vectors `a_* n_a`.

### Theorem `Q4DICJ-1` -- exact standard-diamond identification

`D_infinity(N)` is, up to a rigid motion and the overall bond scale `a_*`, the
standard three-dimensional diamond net.  More explicitly, its two vertex
classes are

\[
 X(\widetilde S_N)=x_0+\Lambda_{\rm FCC},
 \qquad
 X(\widetilde S_{N+1})=x_0+a_*n_1+\Lambda_{\rm FCC},
 \tag{FE08}
\]

where

\[
 \begin{aligned}
 \Lambda_{\rm FCC}
 &=a_*\operatorname{span}_{\mathbb Z}\{n_1-n_4,n_2-n_4,n_3-n_4\}\\
 &={2a_*\over\sqrt3}
 \{(i,j,k)\in\mathbb Z^3:i+j+k\text{ is even}\}.
 \end{aligned}
 \tag{FE09}
\]

Thus the Bravais lattice is `A_3`, equivalently FCC.  In conventional cubic
coordinates its lattice constant is

\[
 a_{\rm cubic}={4a_*\over\sqrt3},
 \tag{FE10}
\]

and the second basis point is displaced by
`(a_cubic/4)(1,1,1)`, the standard diamond basis displacement.  Both
bipartition classes have degree four.

#### Proof

Fix `m_0 in widetilde S_N`.  Every other element of `widetilde S_N` is
`m_0+z` with `z in A_3={z in Z^4:1^Tz=0}`.  Every element of the next
hyperplane is `m_0+e_1+z`.  Equation (FE07) therefore gives the two cosets in
(FE08).  The three displayed root generators become, after relabelling,
`(2a_*/sqrt3)(0,1,1)`, `(2a_*/sqrt3)(1,0,1)`, and
`(2a_*/sqrt3)(1,1,0)`.  Their integer span is exactly the even-coordinate-sum
lattice (FE09).  The offset is `(a_*/sqrt3)(1,1,1)`, which gives (FE10).
Finally, (FE04) and (FE07) make every edge one of the four tetrahedral bonds,
and every vertex has the four inverse incidences on the opposite part.  This
is the conventional FCC-plus-two-point-basis construction of diamond. QED.

This proves an isomorphism of an abstract graph together with its regular
tetrahedral embedding.  It does not yet prove that the graph is physical
space.

## 3. Girth six and complete local ring classification

### Theorem `Q4DICJ-2` -- the first cycles are diamond hexagons

`D_infinity(N)` is simple, connected, bipartite, and has girth six.  Every
simple six-cycle uses exactly three distinct append labels, each once in the
forward and once in the reverse direction.  Consequently its simple
six-cycles are exactly the elementary diamond rings.

#### Proof

Bipartiteness is built into (FE03)--(FE04).  Two distinct `N`-part vertices
share at most one child: if

\[
 m+e_a=m'+e_b,
 \tag{FE11}
\]

then `m'-m=e_a-e_b`, and the ordered root `e_a-e_b` uniquely fixes `(a,b)`.
Hence no two vertices share two children and no simple four-cycle exists.
The two-step moves `e_a-e_b` generate all of `A_3`, so one bipartition class
is connected through length-two paths; every vertex of the other class is
incident to it.  The whole graph is therefore connected.

A six-cycle exists explicitly:

\[
\begin{aligned}
m&\to m+e_1\to m+e_1-e_2
\to m+e_1-e_2+e_3\\
 &\to m-e_2+e_3\to m+e_3\to m.
\end{aligned}
\tag{FE12}
\]

For classification, write any closed six-step path as three forward labels
`p_1,p_2,p_3` and three reverse labels `q_1,q_2,q_3`.  Closure in `Z^4`
forces equality of the two label multisets.  A simple path cannot reverse the
same edge immediately.  If only two labels occurred, avoiding reversal at
each opposite-part vertex would exchange the two labels, forcing their
multiplicities to be equal; three steps cannot split equally between two
labels.  One label is still more clearly impossible.  Therefore three
distinct labels occur, each once in each direction.  Avoiding reversal at
the intervening same-part vertices fixes the usual alternating diamond
hexagon. QED.

If all simple six-cycles are admitted as plaquettes, the infinite graph is
elementary-plaquette complete.  A finite quotient requires the separate
nonwrapping condition in Section 6.

## 4. The physical nonnegative slabs locally exhaust diamond

Let `D_N^+` be the actual finite incidence graph induced by
`S_N sqcup S_(N+1)`.  Every `m in S_N` has four children.  A child
`c in S_(N+1)` has

\[
 \deg_{D_N^+}(c)=|\{a:c_a>0\}|,
 \tag{FE13}
\]

so the finite graph has a simplex boundary and is not globally
coordination-four.

### Theorem `Q4DICJ-3` -- exact rooted local exhaustion

For every graph radius `r>=1`, choose

\[
 N=4r,
 \qquad m_r=(r,r,r,r)\in S_N.
 \tag{FE14}
\]

The rooted radius-`r` ball around `m_r` in `D_N^+` is exactly the rooted
radius-`r` ball around the same vertex in `D_infinity(N)`.  More generally the
same holds whenever `min_a m_a>=r`.

#### Proof

One incidence step changes one coordinate by one.  A path of length at most
`r` from a base point with every coordinate at least `r` therefore never
leaves the nonnegative orthant.  At every vertex of distance strictly less
than `r`, all neighbors which could remain inside the rooted ball are also
nonnegative.  The finite and infinite adjacency rules are otherwise
identical. QED.

Thus the family of q4 slabs locally exhausts the diamond net.  For the
bounded q4 witness, realizing the displayed ball also requires a cap
`R>=N+1`.  One fixed finite cap does not prove an unbounded physical diamond
support.

## 5. The exact physical interface: `Q4-CARRIER/EDGE-LIFT`

The q4 theorem makes `|m>` a count-front state.  Different `m` values can be
mutually exclusive alternatives in one Hilbert space.  A state-transition
graph is not automatically a material lattice whose vertices and links all
exist at once.  The following interface is therefore explicit and
load-bearing.

**`Q4-CARRIER`.**  In one complete F3 parent and one fixed adjacent-depth
block, the front labels in `S_N sqcup S_(N+1)` are realized as coexisting
carrier modes.  The authenticated order/provenance factor and every source,
work, controller, boundary, failure, quarantine, and reference port are
retained.  The carrier action is identity on the retained/future-blind
history factor.

**`EDGE-LIFT`.**  Each keyed append incidence `(m,a)` is realized as one
distinct physical, undirected binary F3 link `e={m,m+e_a}`, oriented from the
`N` part to the `N+1` part only for bookkeeping.  The same label-symmetric
 parameters `(E_R,U_d,h)` act on all lifted links through

\[
 H_{\rm F3}=E_R\sum_en_e+U_d\sum_v(d_v-2)^2-h\sum_eX_e.
 \tag{FE15}
\]

For the perturbative inheritance below, retain the frozen parent domain
`U_d>0` and `|E_R|<2U_d`.

The operator `X_e` toggles the carrier-link occupation; it does **not** undo
an authenticated append, erase a record, or run record formation backward.

The conjunction is called **`Q4-CARRIER/EDGE-LIFT`**.  It is a precise
construction target, not a theorem of URFT, BQ4RSW, or F3.  Under this
antecedent, Theorems Q4DICJ-1--3 remove an otherwise independent choice of a
diamond graph: the lifted bulk support is forced to be diamond by q4
incidence.  They do not prove the antecedent or the stability of the support.

## 6. What the existing F3 theorems inherit

Take a coordination-four bulk block of the lifted graph and freeze `d_*=2`.
Orient each link from the `N` part to the `N+1` part when `n_e=1` and reverse
it when `n_e=0`.  At either bipartition class,

\[
 d_v=2
 \quad\Longleftrightarrow\quad
 \text{two arrows enter and two arrows leave}.
 \tag{FE16}
\]

Thus the F3 low manifold is exactly diamond ice, equivalently a pair of
complementary degree-two fully packed loop coverings.  This is the exact
`CROSS-CW` Hilbert-space bijection on support now identified by q4 incidence.

Because the support has girth six, no non-scalar off-diagonal locked-sector
transition occurs below sixth order.  On every alternating diamond hexagon,
`D6KA` and `CROSS-CW` give

\[
 \langle n\triangle C|H_{\rm eff}^{(6)}|n\rangle=-J_6,
 \qquad
 J_6=h^6\sum_{\pi\in S_6}\prod_{j=1}^{5}
 {1\over\Delta(S_j(\pi))}>0,
 \tag{FE17}
\]

and at symmetric detuning

\[
 \boxed{J_6={63h^6\over8U_d^5}.}
 \tag{FE18}
\]

On a finite, coordination-four, plaquette-complete realization, the complete
sixth-order diagonal and all Feshbach folds are scalar, so

\[
 \boxed{
 H_{\rm eff}^{(6)}=E_{\rm scalar}P_2
 -J_6\sum_{C\in\mathcal P_6}B_C,
 \qquad V_6=0.}
 \tag{FE19}
\]

At `E_R=0`, global occupation complement is an exact symmetry.  These are
inherited theorems, not new coefficient assumptions.

The raw nonnegative slab is not regular at its boundary, so (FE19) is not a
complete finite-slab Hamiltonian theorem.  It applies exactly to each bulk
linked hexagon, and globally only after one of the following additional
premises:

1. an infinite-support quasi-local definition with the required control; or
2. a finite coordination-four periodic completion which introduces no
   wrapping six-cycle absent from the admitted plaquette set.

For reference, a purely combinatorial periodic quotient exists.  Quotient
the FCC translation lattice by the sublattice generated by
`L(e_i-e_4)`, `i=1,2,3`.  For `L>=4` the shortest nontrivial same-part
translation needs at least `2L>6` graph steps, so the quotient has `2L^3`
vertices, degree four, girth six, and no extra wrapping six-cycle.  Selecting
and physically identifying this quotient is nevertheless additional to the
q4 record stream; distinct authenticated front labels are not periodically
identified by BQ4RSW.

There is also a sharp internal comparator.  On the same q4-derived diamond
support, `d_*=1` gives the one-dimer model whose leading pure-kinetic point is
on the ordered R-state side of the published phase diagram.  The positive
U(1) route is therefore not a consequence of the support alone.  It uses the
`d_*=2` ice sector.

## 7. Exact phase and Coulomb-tangent inheritance ledger

On a plaquette-complete periodic completion, dividing (FE19) by `J_6` gives
exactly Shannon et al.'s diamond quantum-ice Hamiltonian

\[
 H_\mu=-\sum_CB_C+\mu\sum_CP_C
 \tag{FE20}
\]

at

\[
 \boxed{\mu=0.}
 \tag{FE21}
\]

The already audited public GFMC and exact-diagonalization evidence imported
by `CROSS-CW` places this pure-kinetic point inside the reported spin-one U(1)
liquid and observes the Maxwell `Phi^2/L` flux scaling.  Therefore the join
supports the following precise conditional statement:

\[
\boxed{
\begin{gathered}
\text{q4 append incidence}
 +\text{`Q4-CARRIER/EDGE-LIFT'}
 +d_*=2\\
\Longrightarrow
\text{the leading controlled F3 bulk Hamiltonian is the published}
\\
\text{pure diamond-ice model with positive numerical U(1)-phase evidence.}
\end{gathered}}
\tag{FE22}
\]

This is the first direct support-shape join between the q4 record-front
construction and the existing F3 U(1) result.  “Positive numerical evidence”
is not upgraded to an exact all-orders F3 phase theorem.

The `Coulomb-tangent` finite-volume targets also acquire a canonical
translation support after a periodic completion.  In particular the matched
electric and magnetic curvatures retain the exact form

\[
 {L\over\Phi_0^2}
 [E_L(+\Phi_0)+E_L(-\Phi_0)-2E_L(0)]
 \tag{FE23}
\]

in the corresponding electric-flux and twisted magnetic-flux sectors.  The
q4/FCC identification supplies the three primitive translation directions;
it does not calculate either infrared stiffness, suppress all-orders flux
mixing, bind `a_*` to an absolute physical length, identify a compact visible
charge, or determine alpha.

## 8. Exact closure boundary

The theorem closes the **combinatorial support-shape gap**:

\[
 \text{two adjacent q4 count-front cosets}
 \Longleftrightarrow
 \text{FCC plus a two-point basis with tetrahedral bonds}
 =\text{diamond}.
 \tag{FE24}
\]

It also proves that every finite-radius diamond neighborhood is already
present in a sufficiently deep bounded nonnegative q4 slab.  The prior
diamond-ice calculation is therefore not geometrically unrelated to the q4
record construction.

The following remain supplied or open, and none is hidden in (FE24):

1. **Physical coexistence.**  q4 front labels must become simultaneous
   carrier modes rather than alternative histories.
2. **Edge and port binding.**  append incidences must carry the binary F3
   links with complete lineage custody; link flips must not erase records.
3. **Support stability.**  no restoring free-energy basin for diamond has
   been proved.  The existing support-force lane explicitly leaves this
   defect response open.
4. **Unbounded/periodic realization.**  one finite cap gives one finite
   simplex slab.  A periodic quotient or a controlled infinite physical
   family is additional.
5. **Thermodynamic and all-orders control.**  the exact result is leading
   order on fixed finite graphs/local linked cells.  A volume-uniform
   expansion and stability of the U(1) basin under every generated operator
   remain open.
6. **Visible electromagnetism.**  the phase is a compact spin-one U(1)
   comparator.  No visible-current map, charge normalization, RG trajectory,
   or alpha calculation follows.
7. **Tensor gravity.**  no symmetric rank-two Gauss law, helicity-two pole,
   universal stress vertex, nonlinear back-reaction, RGRL-B law, gravity, or
   `G` follows.

Accordingly this lane is a support/phase join, not gravity closure.

## 9. Disposition

`INFINITE_TWO_FRONT_Q4_APPEND_GRAPH_IS_EXACTLY_STANDARD_DIAMOND__A3_FCC_TRANSLATION_LATTICE_TWO_COSETS_DEGREE4_TETRAHEDRAL_BONDS_GIRTH6__ALL_SIMPLE_SIX_CYCLES_ARE_DIAMOND_HEXAGONS__NONNEGATIVE_Q4_SLABS_LOCALLY_EXHAUST_THE_NET__Q4_CARRIER_EDGE_LIFT_EXPLICIT_AND_UNPROVED__D2_F3_ICE_BIJECTION_J6_63_OVER8_V6_ZERO_INHERITED_CONDITIONALLY__PURE_MU0_SPIN1_U1_NUMERICAL_COMPARATOR_INHERITED__PHYSICAL_COEXISTENCE_BINDING_STABILITY_PERIODICITY_ALL_ORDERS_VISIBLE_EM_ALPHA_TENSOR_GRAVITY_OPEN`
