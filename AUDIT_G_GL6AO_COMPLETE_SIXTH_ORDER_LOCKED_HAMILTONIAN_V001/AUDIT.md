# Distinct hostile audit — GL6AO complete sixth-order locked Hamiltonian

**Target:** `LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/`  
**Frozen theorem SHA-256:** `f75edcb115c3f7c86c6598f4597366b36e363df2d03ad919cc607b57dfb6b20c`  
**Frozen author-manifest SHA-256:** `c690665043fbbb277aae307a4308e8d30a41f0fbf87be8b1501d0ba86874a494`  
**Frozen author-seal-file SHA-256:** `9df82d1cdc53822bb88b1d419f67db367caf872e2168a9bbeb69d1a6acc9f0ae`  
**Disposition:** `PASS__CANONICAL_KATO_ORDER6_FORMULA_EXACT__DIRECT_AND_FOLDED_WORD_CENSUS_COMPLETE__M3_M2_CANCEL__COMMON_DIAGONAL_MINUS893_OVER1080_M__ONLY_ALTERNATING_HEXAGON_OFFDIAGONAL_MINUS63_OVER8__Q4_AND_FORMAL_LINKED_SCOPE_SOUND__NO_PHASE_POLE_MOMENTUM_CONE_GRAVITY_OR_G`

## 1. Custody and independence

The author froze and sealed eleven GL6AO files before this replay.  Their
exact hashes are pinned in `AUDITED_TARGETS.sha256`; the author manifest and
seal pass.  The six direct dependencies are confined to the sealed GL6AN
author snapshot and its distinct sealed hostile audit.  GL6AL and every later
mutable lane are absent.

The independent replay imports no author code.  In addition to its own exact
checks, the frozen author physics replay passes `120304/120304`, the author
packet passes `82/82`, the upstream GL6AN packet passes `79/79`, the GL6AN
audit packet passes `58/58`, and the independent GL6AN replay passes
`3939/3939`.

## 2. Canonical sixth-order recursion

Put `H=H_0+lambda W`, `H_0P=0`, and

\[
 R=-QH_0^{-1}Q.
\]

For the intermediate-normalized Bloch wave operator `Omega=P+chi`, the exact
equations are

\[
 \chi=\lambda RW(P+\chi)-R\chi K,
 \qquad K=\lambda PW\chi.                                      \tag{A01}
\]

Expanding (A01), using the parity-forced vanishing of odd locked-to-locked
blocks, gives

\[
 K_2=PWRWP=:B,
 \qquad K_4=T_4-A_2B.                                         \tag{A02}
\]

On the sealed homogeneous quotient, `B=bP` and `K_4=dP`, with

\[
 b=-{M\over2},\qquad d=-{7M\over24}.
\]

The fifth wave-operator coefficient then gives exactly

\[
 \boxed{K_6=T_6-bX_4+b^2A_3-dA_2,}                            \tag{A03}
\]

where `X_4` is the sum over the three possible positions of the squared
resolvent.  The first and third positions are adjoints and the middle is
self-adjoint.  Since every lower-order block is scalar on `P`, canonical
Hermitian Kato transport cannot modify this first nonscalar order by a basis
commutator.  A separate exact finite-matrix eigenseries in the audit replay
checks (A03), including all signs, without using the author's recursion.

## 3. Independent quotient and configuration-change classification

The replay rebuilds the period-four incidence from its parent/child rule.  It
has 64 cells, 128 constraint nodes, 256 links, degree four, is connected and
simple, and has girth six.  An independent simple-cycle traversal finds
exactly 256 undirected six-cycles and six cycle cores per link.

No quotient-only six-cycle is hidden in that count.  A six-step lift contains
only three parent-to-parent port differences, so each integer coordinate
displacement lies in `[-3,3]`; closure modulo four therefore forces exact
closure.  The cycles are native local hexagons rather than wrapped artifacts.

For two distinct locked configurations, the symmetric difference is an even
subgraph.  A nonempty even subgraph with at most six edges on a simple graph
of girth six can only be one six-cycle.  Degree two is preserved at each
cycle vertex exactly when the removed and inserted occupations alternate.
Hence the author classification is exhaustive: every order-six off-diagonal
entry is one alternating hexagon, and there are no others.

Dynamic programming over proper toggle subsets, rather than the author's
ordering loop, independently gives

\[
 \sum_{\pi\in S_6}\prod_{j=1}^{5}{-1\over E(S_{\pi,j})}
 =-{63\over8}.                                                  \tag{A04}
\]

Folded terms contain at most four explicit flips and are diagonal because
the quotient has no two- or four-cycle.  Thus (A04) is the complete
off-diagonal coefficient.

## 4. Complete diagonal census

A diagonal six-flip word has even multiplicity on every used link.  After
the `Q` resolvents remove early returns, the only surviving multisets are
`4+2` on two links and `2+2+2` on three links.  The independent multiset
recursion reproduces the repeated-pair weights

\[
 p=2,4,6:\qquad -{1\over4},\ -{1\over16},\ -{1\over36},        \tag{A05}
\]

and the seven three-link weights

\[
 -{9\over32},\ -{9\over16},\ -{29\over144},\ -{109\over144},
 -{41\over32},\ -{337\over864},\ -{209\over1440}.             \tag{A06}
\]

The graph shapes are matching, one adjacent pair plus an isolated edge,
three-edge star, and three-edge path.  Exact enumeration on `Q_4` and a
separate local degree-two proof agree on the author counts.  In particular,
each marked link has two opposite-occupation and one equal-occupation
continuation at each endpoint, so the counts do not depend on the locked
configuration.

Summing (A05)--(A06) gives

\[
 T_6\big|_{\rm diag}
 =-{3\over64}M^3-{215\over576}M^2-{893\over1080}M.             \tag{A07}
\]

The independently reconstructed squared-resolvent fold is

\[
 X_4={5\over32}M^2+{173\over288}M.                             \tag{A08}
\]

Using `A_3=-(M/8)P` and `A_2=(M/4)P`, the three folded pieces are

\[
 -bX_4={5\over64}M^3+{173\over576}M^2,
 \quad b^2A_3=-{1\over32}M^3,
 \quad -dA_2={7\over96}M^2.                                   \tag{A09}
\]

Equations (A07)--(A09) cancel `M^3` and `M^2` exactly and leave

\[
 \boxed{\langle s|K_6|s\rangle=-{893\over1080}M}              \tag{A10}
\]

for every locked `s`.  No order-six flippable-hexagon diagonal potential
survives.

## 5. Finite-volume result and linked-interaction boundary

The audited finite result is therefore

\[
 H_{\rm eff}^{(6)}=-{h^6\over U_d^5}
 \left[{893M\over1080}P+{63\over8}\sum_cT_c\right].           \tag{A11}
\]

This is complete through order six on the declared `Q_4`; its validity is
not assigned to a smaller quotient with wrapped short cycles or to a generic
finite-open boundary.

The infinite expression remains correctly formal.  Each local extension

\[
 \tau_c=P_c\left(\prod_{e\in c}X_e\right)P_c
\]

uses the six degree projectors on the cycle vertices.  The audit closes the
support count explicitly: every term has an 18-link collar, and every link
lies in the support of exactly 18 cycle terms (six as a cycle edge and twelve
through endpoint collars).  Hence the displayed order-six interaction has a
uniform finite interaction norm.  This does not construct an infinite global
locked projector or establish convergence of the all-orders Kato series.

## 6. Promotion attacks and verdict

The audit found no surviving promotion to a thermodynamic phase, a gap or
pole, a physical translation momentum or cone, a photon or graviton, a
stress/Ricci/Einstein law, gravity, or `G`.  The theorem is a finite-order
microscopic collective Hamiltonian plus one formal local linked interaction,
and says no more.

**Hostile verdict: PASS.**
