# GL6CN — COMPLETE DIAGONAL SIXTH-ORDER PURE-`T2` FIRST-SOURCE THEOREM

## Status and exact scope

This packet computes the part left open by `GL6CH`: the complete **diagonal**
pure-`T2` first-source vertex at order `h^6`.  It uses the canonical `GL6AO`
Kato/Feshbach convention, inserts the microscopic pair source before
elimination, differentiates every reduced resolvent, and folds every diagonal
six-flip return history exactly once.

The theorem is pointwise on every locked basis state of the inherited
degree-four, girth-at-least-six parent.  The combinatorial proof applies to any
finite simple degree-four bipartite incidence graph in that domain; a literal
all-history replay is also performed on `Q4`.  The result is a zero theorem:
the complete diagonal order-`h6` pure-`T2` first vertex vanishes.  Combined
with the sealed lower-order result and the sealed `GL6CH` off-diagonal writer,
this closes the complete pure-`T2` first-source operator through sixth order.

It does **not** compute source-second contacts, higher orders, a stationary
phase, record authentication, a bulk or continuum limit, Ricci response,
gravity, or `G`.

## 1. Parent, source, and normalization

Let

\[
 D=\sum_v(k_v-2)^2,\qquad W=-\sum_eX_e,
\]

and insert the six-component local pair source before eliminating the defect
sector:

\[
 H(j)=U_dD+hW+\sum_vj_v^T\mathbf M_v,\qquad
 (\mathbf M_v)_{ab}=Z_{v,a}Z_{v,b}.                    \tag{CN01}
\]

Put `eta=j/U_d` and `r=h/U_d`.  At a locked degree-two node, in pair order
`(01,02,03,12,13,23)`, use the pure tensor basis

\[
 t_1=(1,0,0,0,0,-1),\quad
 t_2=(0,1,0,0,-1,0),\quad
 t_3=(0,0,1,-1,0,0).                                  \tag{CN02}
\]

Because `z_0z_1z_2z_3=1` in every locked word,

\[
 M_{01}=M_{23},\quad M_{02}=M_{13},\quad M_{03}=M_{12},
 \qquad P_T\mathbf M_v=0.                              \tag{CN03}
\]

Thus `P M_T P=0` pointwise.  There is no endpoint source insertion for a
pure-`T2` first derivative; every possible contribution comes from the
source dependence of an intermediate `Q` resolvent or a canonical fold.
The dimensionless order-six coefficient derived below multiplies
`h^6/U_d^6` in `dH_eff/dj`.

## 2. Complete differentiated Kato expression

`GL6AO` proves, with `M=|E|`,

\[
 K_6=T_6-bX_4+b^2A_3-dA_2,\qquad
 b=-{M\over2},\qquad d=-{7M\over24}.                   \tag{CN04}
\]

Here `T6` contains five first-power reduced resolvents; `X4` is the sum of
the three four-flip products with resolvent powers `(2,1,1)`, `(1,2,1)`, and
`(1,1,2)`; and `A_p=PWR^pWP`.

For a parity subset `S`, a source score `m(S)`, and dimensionless defect
energy `E(S)`,

\[
 R_S(\eta)=-{1\over E(S)+\eta m(S)},\qquad
 R'_S(0)=+{m(S)\over E(S)^2}.                           \tag{CN05}
\]

Equivalently, differentiating a factor `R_S^p` multiplies its source-free
value by `-p m(S)/E(S)`.  The bare, order-`h2`, and order-`h4` pure-`T2`
first vertices vanish.  Hence `b'_T=d'_T=0`, while the unsourced `b,d` in
(CN04) are scalars, and the exact first derivative reduces to

\[
 \boxed{K'_{6,T}=T'_{6,T}-bX'_{4,T}+b^2A'_{3,T}-dA'_{2,T}}. \tag{CN06}
\]

No commutativity of a general source-dependent Kato product is assumed in
(CN06).  It follows specifically from the scalar unsourced lower operators
and their already-zero pure-`T2` first derivatives.

The `GL6AO` diagonal classification is exhaustive.  A six-flip diagonal word
has either one link four times and another twice, or three distinct links
twice each.  Proper intermediate returns to `P` are excluded in each `Q`-only
word.  The replay differentiates all retained words: `4+4` for each repeated
pair, four for each of the three `X4` folds, and 60 for each labelled
three-link signature.

## 3. Exact differentiated history kernels

For a repeated pair label the mask columns `1,2,3` denote the source scores
after toggling the first link, the second link, or both.  The complete
coefficients are

| pair energy `p` | `T6':1` | `T6':2` | `T6':3` | `X4':1` | `X4':2` | `X4':3` |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | `3/16` | `3/16` | `1/4` | `-1/2` | `-1/2` | `-1/2` |
| 4 | `3/64` | `3/64` | `1/32` | `-7/32` | `-7/32` | `-3/32` |
| 6 | `1/48` | `1/48` | `1/108` | `-5/36` | `-5/36` | `-1/27` |

For a three-link word, masks `(1,2,3,4,5,6,7)` mean
`(m_0,m_1,m_01,m_2,m_02,m_12,m_012)`.  Up to relabelling, the seven complete
60-word kernels are

| `(p01,p02,p12;t)` | `1` | `2` | `3` | `4` | `5` | `6` | `7` |
|---|---:|---:|---:|---:|---:|---:|---:|
| `(4,4,4;6)` | `1/8` | `1/8` | `3/64` | `1/8` | `3/64` | `3/64` | `1/64` |
| `(2,4,4;4)` | `17/64` | `17/64` | `1/4` | `3/16` | `5/64` | `5/64` | `1/16` |
| `(4,4,6;8)` | `5/48` | `49/576` | `7/192` | `49/576` | `7/192` | `1/54` | `1/144` |
| `(2,2,6;4)` | `7/16` | `19/72` | `5/16` | `19/72` | `5/16` | `19/432` | `49/576` |
| `(2,2,4;2)` | `5/8` | `29/64` | `1/2` | `29/64` | `1/2` | `9/64` | `25/64` |
| `(2,4,6;6)` | `41/192` | `19/108` | `13/72` | `217/1728` | `35/576` | `19/648` | `121/5184` |
| `(4,6,6;10)` | `41/576` | `41/576` | `9/320` | `7/120` | `2/135` | `2/135` | `49/14400` |

The executable also constructs every needed permutation of these signatures,
for 22 labelled cases.  A distinct dual-number implementation inserts
arbitrary rational intermediate scores into the denominators and reproduces
every kernel.  It independently obtains

\[
 (A_2)'[m]=-{1\over4}m,\qquad
 (A_3)'[m]=+{3\over16}m,                                \tag{CN07}
\]

and verifies all four signs in (CN06).

## 4. Universal rooted census

Fix a source node `v`, one incident root link `a`, and a locked word at `v`.
The root has six adjacent links: four with opposite occupation and pair
energy 2, two with equal occupation and pair energy 6.  The remaining
`M-7` links are disjoint from it and have pair energy 4.

For triples containing the root, degree four, girth at least six, and
degree-two locking give the following exhaustive owner-once refinement of
`C(M-1,2)`:

| rooted class | count |
|---|---:|
| matching | `(M^2-21M+116)/2` |
| one adjacent/opposite, root adjacent | `4(M-10)` |
| one adjacent/opposite, root disjoint | `2(M-10)` |
| one adjacent/equal, root adjacent | `2(M-10)` |
| one adjacent/equal, root disjoint | `M-10` |
| star, root minority / majority | `2 / 4` |
| path, root middle: both opposite / one equal / both equal | `4 / 4 / 1` |
| path, root end: both opposite / first opposite / first equal / both equal | `8 / 4 / 4 / 2` |

The counts sum exactly to `(M-1)(M-2)/2`.  One useful independent count is
that the `M-7` root-disjoint links contain `3M-30` adjacent pairs, leaving
`C(M-7,2)-(3M-30)=(M^2-21M+116)/2` matchings.  Of the root-disjoint adjacent
pairs, `2(M-10)` are opposite and `M-10` are equal.  The finite star and path
counts follow by choosing the two endpoints and their locked same/opposite
continuations.  Girth at least six prevents double ownership by a square.

All corrections in which one or both of the other selected links also touch
the source node are then enumerated over the six possible degree-two central
words.  For **each** of those six words, and for **each** of its four ports,
the resulting local odd-parity coefficients are identical:

\[
\begin{aligned}
 c^{T_6}_1&={15\over128}M^2+{3049\over3456}M+{8653\over4800},
&c^{T_6}_3&={49\over576},\\
 c^{X_4}_1&=-{5\over16}M-{487\over432},
&c^{A_3}_1&={3\over16},\qquad c^{A_2}_1=-{1\over4}.     \tag{CN08}
\end{aligned}
\]

Here subscript 1 means a singleton set of toggled ports and subscript 3 its
three-port complement.  The `X4`, `A3`, and `A2` terms have no odd
three-port pattern.  After the exact Kato folds,

\[
 c^{K_6}_1={1\over128}M^2+{283\over1152}M+{8653\over4800},
 \qquad c^{K_6}_3={49\over576}.                         \tag{CN09}
\]

This is a coverage proof, not an extrapolation from selected `Q4` nodes:
the seven energy classes are the complete `GL6AO` diagonal-history classes;
the rooted table partitions every possible second and third link; and all six
central locked words are evaluated exactly.

## 5. Pointwise cancellation theorem

For a locked local word `s`, define the six-pair word after one local port
flip by

\[
 w_a=\mathbf M_v(s\mathbin\triangle\{a\}),\qquad a=0,1,2,3.
\]

For every pair component, two of the four single-port flips change its sign
and two preserve it.  Therefore the stronger full-vector identity holds:

\[
 \boxed{\sum_{a=0}^3w_a=0_6}.                          \tag{CN10}
\]

Flipping the complementary three ports gives the full complement of the
one-port-flipped word.  Pair products are invariant under full complement,
so

\[
 \mathbf M_v(s\mathbin\triangle(\{0,1,2,3\}\setminus\{a\}))=w_a. \tag{CN11}
\]

Every even local flip pattern remains `T2`-dark because it preserves
`z_0z_1z_2z_3=1`.  Equations (CN08), (CN10), and (CN11) now show separately
that

\[
 P_TT'_{6,\mathrm{diag}}=P_TX'_{4,\mathrm{diag}}
 =P_TA'_3=P_TA'_2=0                                    \tag{CN12}
\]

at every node and on every locked basis state.  Inserting (CN12) into
(CN06) yields the main theorem:

\[
 \boxed{
 \langle s|K'_{6,T}|s\rangle=0
 }
 \quad\text{for every locked }s.                       \tag{CN13}
\]

Because a first derivative is linear and (CN13) is pointwise in the source
node and each of the three `T2` directions, the theorem includes every
spatially nonuniform pure-`T2` first-source profile.  It does not include a
second derivative with two source insertions; mixed two-site or nonuniform
source-source contacts remain open.

## 6. Literal `Q4` all-history replay

The proof above establishes the declared domain.  As an independent finite
parent check, the executable constructs `Q4` with 128 nodes, 256 links, and
768 owner-once adjacent-link pairs.  It verifies the locked `GL6CC`
background, locates representatives of all six local locked words, and for
each representative includes every unordered pair and triple meeting the
source support:

\[
 1014={256\choose2}-{252\choose2},\qquad
 128020={256\choose3}-{252\choose3}.                    \tag{CN14}
\]

It also performs legal incident-ring toggles around one fixed source node
and checks three changed environments/local words.  Thus nine literal cases
are evaluated.  In every case, `T6'`, `X4'`, `A3'`, `A2'`, and their complete
Kato combination have zero contraction with each `T2` basis vector.  The
literal coefficients match (CN08) port by port.  The replay passes
`10775/10775` exact rational checks.

The nine cases validate the implementation against the finite `Q4` parent;
they are not used as a substitute for the universal rooted proof.

## 7. Integrated first-source corollary through `h6`

On the shared degree-four, girth-at-least-six `Q4`/linked-parent convention,
the sealed results now give:

\[
 PM_TP=V_T^{(2)}=V_T^{(4)}=V_{T,\mathrm{diag}}^{(6)}=0. \tag{CN15}
\]

`GL6CH` supplies the complete off-diagonal sixth-order term.  Therefore, for
an arbitrary pure-`T2` source profile `j_T`, the complete **first-source**
effective operator through sixth order is

\[
\boxed{
 \left.{\partial\over\partial\epsilon}
 H_{\rm eff}(\epsilon j_T)\right|_{\epsilon=0}^{[0,2,4,6]}
 ={105\over16}{h^6\over U_d^6}
  \sum_cT_c\sum_{v\in c}j_{T,v}^T\Theta_{v,c}
 +O\!\left({h^8\|j_T\|\over U_d^8}\right)
}.                                                       \tag{CN16}
\]

Here `T_c` is the owner-once alternating six-cycle toggle and
`Theta_{v,c}=e_{ab}-e_{cd}` is the complementary-pair tensor at a cycle
vertex.  Equation (CN16) means the first pure tensor source does not create a
diagonal potential through this order: its sole surviving action is the
configuration-changing `GL6CH` writer.

## 8. Exact ceiling

GL6CN closes only the first-source operator through sixth order.  It leaves
open:

- source-second contacts, including mixed two-site/nonuniform
  source-source derivatives;
- pure-`T2` first-source vertices at order `h8` and above;
- the complete contact-plus-spectral Hessian in a stationary accumulated
  state;
- derivation of a candidate field from authenticated record lineage;
- thermodynamic/refinement limits, a common cone, Ricci or other continuum
  tensor form, gravity, and `G`.

No phase, record, metric, Ricci, gravity, graviton, or numerical-`G` claim is
made here.

`PASS__GL6CN_COMPLETE_DIAGONAL_H6_T2_FIRST_SOURCE__POINTWISE_ZERO__ALL_6_LOCKED_WORDS__ROOTED_OWNER_ONCE_CENSUS__ALL_REPEATED_PAIR_AND_TRIPLE_RETURN_HISTORIES__DUAL_RESOLVENT_SIGN_REPLAY__Q4_9_CASE_LITERAL_REPLAY__COMPLETE_FIRST_T2_THROUGH_H6_EQUALS_GL6CH_WRITER__NONUNIFORM_FIRST_SOURCE_CLOSED__SOURCE_SECOND_HIGHER_PHASE_RECORD_BULK_RICCI_GRAVITY_G_OPEN`
