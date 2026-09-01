# One-bridge `E2` no-go and retained-FPSS active-block escape theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6AB_E2_MULTI_CONNECTOR_TRIANGLE_V001`  
**Short name:** `GL6AB V001`  
**Date:** 2026-08-31  
**Status:** author frozen after independent hostile pre-freeze review  
**Claim class:** exact all-orders symmetry theorem for a genuinely one-port
bridge; exact shortest simple inherited connector-only FPSS graph motif,
without external port-resolving dressing, that removes that obstruction; exact
complete-`N=1` sixteen-link active-block F3 mixed-jet witness with nonzero
`E2` projection

**Not claimed:** that an arbitrary full-parent dressing preserves the
one-bridge no-go; a full-rank propagated `E2` block at the displayed order;
an autonomous support law; an authenticated physical atlas, spatial direction,
length, momentum, or cone; an all-time Schur kernel; Ward/Bianchi closure; an
infrared operator; Ricci or Einstein form; gravity; or `G`.

## 1. Question left by GL6Z

Let

\[
 W:=\mathbb R^{{\cal E}_4},\qquad
 P_{A,a}:=\mathbf 1_{\{a\in A\}},                 \tag{AB01}
\]

where `A` runs over the six unordered pairs of four q4 half-ports.  The
audited GL6Z first cross-cell coefficient factors through `P` and therefore
has

\[
 E:=\ker P^T\cong E_2                               \tag{AB02}
\]

as an exact two-dimensional null space.  The present packet asks two sharply
separated questions.

1. Is that null forced at every order by one physical shared-child bridge?
2. If so, what is the smallest inherited FPSS structure that can escape it,
   and does the complete finite F3 active block actually have a nonzero
   coefficient?

The answers are: **yes, under an explicit one-port equivariance premise; and,
within the inherited connector-only graph with no external port-resolving
dressing, the first escape is a distinct-child sibling triangle.**  The complete
sixteen-link `N=1` active block gives a nonzero exact witness while all other
raw factors remain retained spectators.  Thus the GL6Z null is not an
all-orders no-go for F3.

## 2. Theorem `GL6AB-1` -- exact all-orders one-port null

Fix a port `a`.  Its stabilizer

\[
 H_a:=\{\sigma\in S_4:\sigma(a)=a\}\cong S_3       \tag{AB03}
\]

has exactly two orbits on the pair set: the three pairs containing `a` and
the three pairs not containing `a`.  Hence

\[
 W^{H_a}=\operatorname{span}\{u_a^{\rm in},u_a^{\rm out}\},             \tag{AB04}
\]

where

\[
 (u_a^{\rm in})_A=P_{A,a},\qquad
 (u_a^{\rm out})_A=1-P_{A,a}.                       \tag{AB05}
\]

Both vectors belong to `im P`, because

\[
 \mathbf1_W={1\over2}\sum_{c=1}^4P_{\cdot,c}.      \tag{AB06}
\]

Equivalently, if an `H_a`-invariant vector has value `alpha` on incident
pairs and `beta` on nonincident pairs, then `P^Tx=0` implies
`3 alpha=0` at `a` and `alpha+2 beta=0` at every other port.  Therefore

\[
 \boxed{W^{H_a}\cap E=\{0\}.}                     \tag{AB07}
\]

Now join port `a` of cell `m` to port `b` of cell `n` by one inherited
shared-child interaction.  Let `K_nm` be any complete dressed **linear**
pair-response kernel, at any two times or frequencies, which satisfies the
independent endpoint equivariance

\[
 \rho_n(h_b)K_{nm}\rho_m(h_a)^{-1}=K_{nm}
 \quad\forall(h_b,h_a)\in H_b\times H_a.           \tag{AB08}
\]

This premise permits arbitrary repetitions of the same bridge and arbitrary
local state, contact, source, read, and self-energy dressing only when those
objects preserve the displayed endpoint stabilizers.  Let

\[
 R_a={1\over|H_a|}\sum_{h\in H_a}\rho(h)           \tag{AB09}
\]

be the Reynolds projector.  Equation (AB08) gives

\[
 K_{nm}=R_bK_{nm}R_a.                              \tag{AB10}
\]

By (AB07), `R_a Pi_E=0=Pi_E R_a`.  Thus

\[
 \boxed{K_{nm}\Pi_E=0=\Pi_EK_{nm}}                \tag{AB11}
\]

exactly, at every response-time order, every frequency, and every
perturbative order in repeated use of that one bridge.

Equation (AB11) is deliberately conditional on (AB08).  It applies to an
isolated or exactly refocused one-bridge block.  It does **not** apply merely
because two cells have one direct common child.  Other incident FPSS bridges,
port-addressed controls, unequal boundary degrees, or a state which resolves
a second port can break `H_a x H_b`; the complete parent must then be
calculated.

## 3. Under connector-only endpoint equivariance, two ports are necessary

For two distinct ports `a,c`, define

\[
 w^{ac}_A=P_{A,a}P_{A,c}=\delta_{A,\{a,c\}}.       \tag{AB12}
\]

Since

\[
 P^TP=2I_4+\mathbf1\mathbf1^T,qquad
 (P^TP)^{-1}={1\over2}I_4-{1\over12}\mathbf1\mathbf1^T,                \tag{AB13}
\]

the exact orthogonal projection obeys

\[
 \boxed{\|\Pi_Ew^{ac}\|^2={1\over3}>0.}           \tag{AB14}
\]

Thus two distinct port insertions are the minimum representation-theoretic
requirement for an `E2` component within an inherited connector-only block
which has no external port-resolving state or apparatus.  A bare two-edge
wedge `m-l-n` does not
yet transmit `E_m` to `E_n`: each external endpoint still exposes only one
port.  The middle cell has two ports, but the endpoint Reynolds projections
continue to kill `E2` under their one-port symmetry premise.

The shortest simple inherited connector-only FPSS motif with two distinct
ports at **both** endpoints is a three-edge triangle with a direct bridge and a two-edge alternative
path.  Let `x in S_(N-1)` and let `a,b,c` be distinct.  Put

\[
 m=x+e_b,\qquad n=x+e_a,\qquad \ell=x+e_c.         \tag{AB15}
\]

Then:

\[
\begin{array}{c|c|c}
\text{cell edge}&\text{shared child}&\text{endpoint ports}\\ \hline
m-n&x+e_a+e_b&(a\text{ at }m,\ b\text{ at }n)\\
m-\ell&x+e_b+e_c&(c\text{ at }m,\ b\text{ at }\ell)\\
\ell-n&x+e_c+e_a&(a\text{ at }\ell,\ c\text{ at }n).
\end{array}                                        \tag{AB16}
\]

The endpoints therefore see `{a,c}` and `{b,c}`.  This is the inherited
`A3` sibling triangle

\[
 (e_a-e_b)+(e_c-e_a)+(e_b-e_c)=0.                 \tag{AB17}
\]

It must not be confused with three parents lying around one common child;
that other orientation can reuse the same port at an endpoint and does not
by itself remove the one-port stabilizer.  Representation theory only
**permits** the sibling-triangle response.  A nonzero coefficient still has
to be earned from the physical parent.

## 4. Complete `N=1` active block inside the full physical parent

Use GL6Y's complete `N=1` FPSS member.  Its four cells are indexed by
`i=0,1,2,3`, and its sixteen active physical link registers are

\[
 e_{ia}:=\ell_{e_i,a},\qquad i,a=0,1,2,3.          \tag{AB18}
\]

No active link is deleted.  The interaction graph has all six parent-clique edges
inside each cell and all six shared-child bridges:

\[
 {\cal L}_1=
 \{\{e_{ia},e_{ic}\}:i=0,\ldots,3,\ a<c\}
 \cup
 \{\{e_{ij},e_{ji}\}:i<j\}.                      \tag{AB19}
\]

Thus `|L_1|=4*6+6=30`.  The full raw FPSS device also retains its eighty-four
blank nonedges, six parent guards, formation, route, controller, clock, work,
failure, quarantine, boundary, reference, and source/read port factors.
GL6Y proves that the raw nonedges remain blank/invariant and that the retained
non-active factors have the inherited factorized/common future during this
source-off interval.  Hence every displayed commutator reduces exactly to
the following **complete active-link block**; the other raw factors have not
been deleted.

Expanding the inherited degree law and discarding only its common scalar gives

\[
 H_1=-h\sum_{e=1}^{16}X_e+\delta\sum_{e=1}^{16}n_e
 +2U_d\sum_{\{e,f\}\in{\cal L}_1}n_en_f,          \tag{AB20}
\]

where

\[
 \delta=\Delta+2(1-2d_\star)U_d.                 \tag{AB21}
\]

Choose the admissible one-parameter exact witness, with arbitrary energy
scale `E_0>0`,

\[
 h=U_d=E_0,\qquad d_\star=2,\qquad\Delta=6E_0,
 \qquad\delta=0.                                  \tag{AB22}
\]

Then `H_1=E_0 \widehat H_1`, where the dimensionless integer operator
`widehat H_1` is the one replayed below.  Setting `E_0=1` is used only inside
that dimensionless arithmetic and is not a physical unit choice.

The pair order is

\[
 (01,02,03,12,13,23),                              \tag{AB23}
\]

and choose an exact basis of `E=ker P^T`,

\[
 e^{(1)}=(1,-1,0,0,-1,1),\qquad
 e^{(2)}=(1,0,-1,-1,0,1).                         \tag{AB24}
\]

For cell `i` and pair `A={a,b}`, let

\[
 M_{i,A}:=Z_{ia}Z_{ib},                            \tag{AB24a}
\]

and let

\[
 M_{i,e}:=\sum_{A\in{\cal E}_4}e_AM_{i,A}.        \tag{AB25}
\]

The GL6Z finite-family source compiler realizes `A` by applying the six
commuting GL6V pair phases with the signed amplitudes in `e^(2)`; every source
ancilla is uncomputed to blank.  `B` is a signed deterministic function of
the complete terminal read of all sixteen active links.  The six pair values
are not substituted for that complete read.  Thus (AB25) selects one linear
source/read direction inside the already normalized apparatus; it adds no
autonomous force or pair interaction.

Inside the exact all-formation/KEEP branch, the source-off prewait state is

\[
 \rho_\tau=e^{-iH_1\tau/\hbar}|0^{16}\rangle\langle0^{16}|
 e^{+iH_1\tau/\hbar}.                             \tag{AB25a}
\]

This is the same branch-normalized entrance construction as GL6Z; no
stationarity or time-translation invariance is assumed.

Take the source direction `A=M_(0,e^(2))`, the read direction
`B=M_(1,e^(1))`, and the common blank-prewait entrance state.  Define the
dimensionless exact mixed commutator integers

\[
 \widehat Q_{p,r}:=
 \langle0^{16}|\operatorname{ad}_{\widehat H_1}^{p}
 ([\operatorname{ad}_{\widehat H_1}^{r}B,A])|0^{16}\rangle,
 \qquad
 Q_{p,r}=E_0^{p+r}\widehat Q_{p,r}.                \tag{AB26}
\]

## 5. Theorem `GL6AB-2` -- exact full-active-block `E2` escape

Exhaustive integer recursion on all `2^16=65536` computational states gives

\[
 \boxed{\widehat Q_{p,r}=0\quad\text{for every }p+r\le18,}            \tag{AB27}
\]

and, at the first nonzero total order,

\[
\boxed{\begin{array}{c|r}
(p,r)&\widehat Q_{p,r}\\ \hline
(16,3)&-67530641899520\\
(15,4)&-135061283799040\\
(14,5)&-165289008267264\\
(13,6)&-158213815304192\\
(12,7)&-132193339604992\\
(11,8)&-105585215864832\\
(10,9)&-87422498291712\\
(9,10)&-77413660606464\\
(8,11)&-70744951750656\\
(7,12)&-62882751594496\\
(6,13)&-51947031625728\\
(5,14)&-38655920619520\\
(4,15)&-25051091828736\\
(3,16)&-13213688233984\\
(2,17)&-4594932678656.
\end{array}}                                      \tag{AB28}
\]

The other total-order-nineteen distributions vanish exactly:
`(p,r)=(18,1),(17,2),(1,18),(0,19)`; the `r=0` entry also vanishes trivially
because the two cell observables commute at the response entrance.

The `r=3,p=16` coefficient was independently evaluated for all four entries
of the nonorthogonal basis (AB24).  The exact `E2 x E2` bilinear block is

\[
 \boxed{\widehat Q^{E}_{16,3}=-67530641899520
 \begin{pmatrix}1&1\\1&1\end{pmatrix}.}           \tag{AB29}
\]

It has rank one.  Hence the inherited multi-connector dressing propagates
one genuine `E2` shear combination at this order; it does not yet produce a
full-rank two-shear block.  The `S4` orbit of other port-pair orientations is
the shortest next finite-family assembly to test for `E2` rank two.

### Corollary -- exact algebraic orbit span

The basis (AB24) has Gram matrix

\[
 G_E=\begin{pmatrix}4&2\\2&4\end{pmatrix}.         \tag{AB29a}
\]

Therefore the dimensionless linear operator represented by the bilinear
matrix in (AB29) is `G_E^(-1) widehat Q^E_(16,3)`.  Its nonzero operator
eigenvalue is

\[
 \boxed{\widehat Q_{16,3}/3}                       \tag{AB29b}
\]

on the line generated by

\[
 v_{01}=e^{(1)}+e^{(2)}=(2,-1,-1,-1,-1,2).         \tag{AB29c}
\]

The three distinct `S4` orbit lines may be represented by

\[
\begin{aligned}
v_{01}&=(2,-1,-1,-1,-1,2),\\
v_{02}&=(-1,2,-1,-1,2,-1),\\
v_{03}&=(-1,-1,2,2,-1,-1).
\end{aligned}                                      \tag{AB29d}
\]

They obey

\[
 v_{01}+v_{02}+v_{03}=0,qquad
 \|v_{0i}\|^2=12,qquad
 v_{0i}\cdot v_{0j}=-6\ (i\ne j).                 \tag{AB29e}
\]

Consequently any two lines span `E2`.  If three physically covariant copies
of the rank-one operator could be composed with equal weight, their
**algebraic average** would have bilinear matrix and operator eigenvalue

\[
 \overline Q_E={\widehat Q_{16,3}\over6}G_E,\qquad
 \lambda(\overline Q_E)={\widehat Q_{16,3}\over6},                 \tag{AB29f}
\]

while their unaveraged sum would be `(widehat Q_(16,3)/2)G_E`.  This is an exact
orbit-span sufficiency result, not yet a physical composition theorem.  The
current packet has not proved that all three oriented sources/reads coexist
with equal coefficient in one authenticated, glued parent kernel.  That is
the finite target for `GL6AC`, not a conclusion imported into GL6AB.

### Exact replay identity

At (AB22), write a computational word as `s in {0,1}^16` and put

\[
 D(s)=2\sum_{\{e,f\}\in{\cal L}_1}s_es_f.          \tag{AB30}
\]

Then

\[
 (\widehat H_1\psi)(s)=D(s)\psi(s)-\sum_{q=1}^{16}
 \psi(s\mathbin{\mathsf{xor}}2^{q-1}).             \tag{AB31}
\]

Let `v_a=widehat H_1^a|0>`, and define the exact integer table

\[
 T_{abc}:=\langle v_a|B \widehat H_1^b A|v_c\rangle.                 \tag{AB32}
\]

Binomial expansion of both nested commutators gives

\[
\boxed{
\widehat Q_{p,r}=\sum_{l=0}^{p}\sum_{j=0}^{r}
(-1)^{l+j}{p\choose l}{r\choose j}
\left(T_{p-l+r-j,j,l}-T_{l+j,r-j,p-l}\right).}    \tag{AB33}
\]

Equations (AB19), (AB23)--(AB25), and (AB30)--(AB33) are a complete exact
replay specification.  The packet verifier evaluates the representation and
graph identities; the supplied arbitrary-precision replay evaluates every
nontrivial `r>=1` entry through total order twenty.  The `r=0` entries vanish
identically because the two response-entrance cell observables commute.  No
fitted continuum operator enters.

## 6. Physical normalized CTP coefficient

With GL6Z's normalized source convention,

\[
 {\cal G}^{R}_{BA}(t,0)=
 {iE_\star^2\over2\hbar}\Theta(t)
 \langle[B(t),A]\rangle_{\rho_\tau}.              \tag{AB34}
\]

Equations (AB26)--(AB28) imply the exact selected bivariate coefficient

\[
\boxed{
[\tau^{16}t^3]_{t\to0^+}{\cal G}^{R}_{BA}
={\widehat Q_{16,3}E_\star^2E_0^{19}\over12\,16!\,\hbar^{20}}
=-{103043582\over383107725}{E_\star^2E_0^{19}\over\hbar^{20}}\ne0.} \tag{AB35}
\]

For fixed `d_star=2`, every finite nested-commutator coefficient is a finite
polynomial in `(h,U_d,Delta)`.  In particular,

\[
 Q_{16,3}(h,U_d,\Delta;d_\star=2)
 \in\mathbb Z[h,U_d,\Delta]                       \tag{AB35a}
\]

after the displayed integer observable normalization.  Since this polynomial
is nonzero at `(E_0,E_0,6E_0)` for every `E_0>0`, continuity gives an open
parameter neighborhood on which it remains nonzero.  Thus `delta=0` is an
exact evaluation slice, not a claim that the multi-connector effect exists
only at a tuned point.  No characterization of the polynomial zero locus is
needed here.

This is evaluated inside the same all-formation/KEEP branch and complete
sixteen-link read used by GL6Y/GL6Z.  In matched BREAK the transverse terms
are absent, the Hamiltonian and pair queries are diagonal, and every such
response coefficient vanishes.

The physical interpretation is narrow but important.  The first cubic
shared-child operator jet did not need replacement.  Its leading
`tau^4` scalar dressing was too short to distinguish two-port shear.  At a
longer exact accumulation horizon, coherent dressing by the other inherited
shared-child connectors feeds back into that same jet and opens `E2`.

## 7. Program frame versus physical atlas

Every register, formation lineage, degree interaction, source, and complete
read in (AB18)--(AB35) is physical in the fixed FPSS program.  Therefore
(AB35) is a physical finite, program-conditional, record-qualified response.

The labels `i,a` and the triangle organization enter this calculation through
the supplied FPSS program.  Audited GL6AA separately proves a record-
authenticated endpoint-identity/port-label query, inverse consistency, and
cocycle for the selected relational atlas.  This packet does not formally
compose its order-nineteen response apparatus with that atlas packet, nor
does either result assign physical length or establish a common cone.  Until
that exact composition is sealed, the present result remains a **multi-cell
`E2` response**, not spatial shear propagation or curvature.

## 8. Result and shortest next gate

The exact earned chain is

\[
\boxed{\begin{gathered}
\text{one equivariant physical bridge}
\Rightarrow E_2\text{ null at all orders};\\
\text{under connector-only endpoint equivariance, two ports are required};\\
\text{the distinct-child FPSS sibling triangle is the minimal simple motif};\\
\text{complete }N=1\text{ F3 active block has a nonzero rank-one }E_2
\text{ coefficient at total order }19.
\end{gathered}}                                    \tag{AB36}
\]

The next response calculation is not more one-bridge machinery.  It is to
construct the covariant same-parent source/read composition whose algebraic
target is (AB29f), and test whether the rotated rank-one images coexist and
fill both `E2` directions while preserving the already derived `A1/T2`
principal scaling.
Audited GL6AA has independently turned the selected port/cell organization
into an earned relational chart.  Its formal composition with this
order-nineteen response remains to be proved.  Only after the results compose
should the complete six-channel kernel, contacts,
Schur quotient, common cone, conservation, and infrared operator be tested.

`PASS__ONE_PORT_S3xS3_EQUIVARIANT_BRIDGE_E2_NULL_ALL_ORDERS__CONNECTOR_ONLY_NO_EXTERNAL_PORT_RESOLVER_TWO_DISTINCT_ENDPOINT_PORTS_MINIMAL__TWO_EDGE_WEDGE_INSUFFICIENT__DISTINCT_CHILD_A3_SIBLING_TRIANGLE_SHORTEST_SIMPLE_MOTIF__FULL_N1_SIXTEEN_LINK_THIRTY_INTERACTION_F3_ACTIVE_BLOCK_RETAINED__EXACT_MIXED_E2_COEFFICIENTS_ZERO_TOTAL_ORDER_LE18__FIRST_NONZERO_TOTAL_ORDER19__TAU16_T3_CTP_COEFFICIENT_EXACT_NONZERO_OPEN_PARAMETER_NEIGHBORHOOD__DISPLAYED_E2_BLOCK_RANK1__MATCHED_BREAK_ZERO__PROGRAM_CONDITIONAL_PHYSICAL_RESPONSE_BUT_GL6AA_COMPOSITION_LENGTH_CONE_WARD_IR_RICCI_GRAVITY_G_OPEN`
