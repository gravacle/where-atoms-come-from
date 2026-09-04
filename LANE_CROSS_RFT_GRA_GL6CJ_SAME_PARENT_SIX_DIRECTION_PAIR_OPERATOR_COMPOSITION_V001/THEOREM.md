# GL6CJ — SAME-PARENT SIX-DIRECTION PAIR-OPERATOR COMPOSITION THEOREM

## Status and scope

This packet composes two exact derivatives of one microscopic source,

\[
 H(j)=H(0)+\sum_v j_v^TM_v,                              \tag{CJ01}
\]

where the six-pair source is inserted **before** Feshbach/Kato elimination.
Its locked diagonal evaluation reads `A1+E`, while its first global
configuration-changing tensor vertex at order `h^6` writes `T2`.  The two
maps have complementary kernels and exact reconstruction inverses.

This is a perturbative finite-range operator-jet/source-access theorem on the
declared `Q4` parent.  It does **not** prove that `j` is autonomously generated
by records, select a stationary state or phase, compute a response function,
identify a metric, establish `RGRL-B`, derive Ricci or Einstein dynamics,
prove gravity, or calculate `G`.

## 1. Six-pair decomposition

Use pair order

\[
 {cal P}=(01,02,03,12,13,23),                           \tag{CJ02}
\]

and the orthogonal decomposition

\[
 \mathbb R^6=A_1\oplus E\oplus T_2,
 \qquad I_6=P_A+P_E+P_T,                                 \tag{CJ03}
\]

where `A=(1,1,1,1,1,1)`, `P_A=AA^T/6`, and

\[
 T_2=\operatorname{span}\{e_{01}-e_{23},
 e_{02}-e_{13},e_{03}-e_{12}\}.                         \tag{CJ04}
\]

The exact replay constructs all three projectors and verifies ranks
`1,2,3`, mutual orthogonality, and (CJ03).

## 2. Locked diagonal read: exact rank three on `A1+E`

At a degree-two locked node let

\[
 M(s)=(z_az_b)_{a<b},\qquad z_a=1-2s_a,qquad \sum_as_a=2.
                                                               \tag{CJ05}
\]

There are six locked words.  Complementary words `s` and `1-s` give the same
pair vector, so there are three distinct vectors, each twice.  Every one
obeys `P_TM(s)=0`, and the six-by-six evaluation map

\[
 {cal D}:j\longmapsto d_s=j^TM(s)                       \tag{CJ06}
\]

has rank three and kernel exactly `T2`.  Its normal operator is

\[
 \boxed{{\cal D}^*{\cal D}=\sum_{s:\,|s|=2}M(s)M(s)^T
       =4P_A+16P_E.}                                      \tag{CJ07}
\]

Consequently the compatible locked read reconstructs the `A1+E` component.
Writing the displayed operator as `R_D`, its exact generalized-inverse
identities are `R_D D=P_A+P_E` and `D R_D=Pi_im(D)`:

\[
 \boxed{
 j_{A+E}={1\over4}P_A\sum_sd_sM(s)
         +{1\over16}P_E\sum_sd_sM(s).}                   \tag{CJ08}
\]

This is not merely a bare-source accident.  Define the exact displayed
diagonal first-source jet through fourth order, with `r=h/U_d`, by

\[
 \begin{aligned}
 V_{v,\rm diag}^{[0,2,4]}(j;s)
 &=j^TM(s)+r^2j^T[-M(s)]\\
 &\quad+r^4j^T[-\tfrac49{\bf1}_6-\tfrac{37}{12}M(s)] .  \tag{CJ09}
 \end{aligned}
\]

so its formal leading rank remains three and its source support remains
`A1+E` through the displayed diagonal orders.  Equation (CJ09) is a
truncation identity, not an asymptotic equality for the complete vertex; it
does not claim that any diagonal order-six term has been classified.

## 3. Global sixth-order tensor writer at one node

For every elementary hexagon `c` incident on node `v`, let `{a,b}` be its two
local cycle ports and `{c,d}` the complementary pair.  Define

\[
 \Theta_{v,c}=e_{ab}-e_{cd}.                             \tag{CJ10}
\]

`GL6CH` proves that the pure-`T2` part of the same pre-Feshbach source changes
the order-six cycle-toggle coefficient by

\[
 \lambda_T\,j_v^T\Theta_{v,c},\qquad
 \lambda_T={105\over16}{h^6\over U_d^6}.                \tag{CJ11}
\]

Define the unscaled local writer map

\[
 {cal W}_v:j\longmapsto
 w_c=j^T\Theta_{v,c},\qquad c\ni v.                     \tag{CJ12}
\]

An independent reconstruction of all 256 `Q4` hexagons proves, at every one
of its 128 constraint nodes:

- exactly twelve elementary hexagons meet the node;
- each of the six unordered local port pairs occurs exactly twice;
- all `Theta_{v,c}` lie in `T2`; and

\[
 \boxed{{\cal W}_v^*{\cal W}_v
 =\sum_{c\ni v}\Theta_{v,c}\Theta_{v,c}^T=8P_T.}       \tag{CJ13}
\]

Thus `W_v` has rank three, kernel `A1+E`, and the exact reconstruction is

\[
\boxed{j_T={1\over8}\sum_{c\ni v}w_c\Theta_{v,c}.}    \tag{CJ14}
\]

Equivalently `R_W W=P_T` and `W R_W=Pi_im(W)`.  These are left inverses of
the forward maps on their source sectors (and right inverses of the
corresponding synthesis maps), rather than nonexistent inverses onto all
arbitrary six- or twelve-component target data.

If the actual dressed coefficient change is
`delta a_c=lambda_T w_c`, then

\[
 \boxed{j_T={2\over105}{U_d^6\over h^6}
              \sum_{c\ni v}\delta a_c\Theta_{v,c}.}    \tag{CJ15}
\]

The complete canonical order-six pair gradient also has scalar and `E`
pieces before irrep projection.  Equations (CJ11)--(CJ15) isolate the
non-scalar tensor writer proved by `GL6CH`; they do not assert that other
off-diagonal source projections are absent.

## 4. One source, complementary operator slots, full rank

For every locked word and every incident cycle,

\[
 M(s)^T\Theta_{v,c}=0.                                  \tag{CJ16}
\]

The diagonal map and tensor-writer map therefore occupy orthogonal source
sectors.  They also occupy disjoint operator support: (CJ06) consists of
locked-basis diagonal matrix elements, whereas a six-link toggle has no
diagonal matrix element.  They cannot cancel one another.

Stacking the two maps gives

\[
 {\cal C}_v=\begin{pmatrix}{\cal D}\\{\cal W}_v\end{pmatrix},
 \qquad
 \boxed{{\cal C}_v^*{\cal C}_v=4P_A+16P_E+8P_T},        \tag{CJ17}
\]

which is positive definite and has rank six.  Combining (CJ08) and (CJ14)
gives an explicit left inverse `R_C C_v=I_6` on the full pair-source space:

\[
 \boxed{
 j={1\over4}P_A\sum_sd_sM(s)
  +{1\over16}P_E\sum_sd_sM(s)
  +{1\over8}\sum_{c\ni v}w_c\Theta_{v,c}.}             \tag{CJ18}
\]

The replay checks (CJ18) on every standard pair coordinate at all 128 nodes,
not only at one symmetry representative.

At the effective-operator level, define brackets with subscripts to mean the
indicated independently classified projection and perturbative orders.  The
two exact identities are

\[
 \boxed{
 \left[\partial_\epsilon H_{\rm eff}(\epsilon j)|_{\epsilon=0}
 \right]_{\rm diag}^{[0,2,4]}
 =V_{\rm diag}^{[0,2,4]}(j_{A+E}),}                     \tag{CJ19a}
\]

\[
 \boxed{
 \left[\partial_\epsilon H_{\rm eff}(\epsilon j)|_{\epsilon=0}
 \right]_{{\rm off},T_2}^{(6)}
 ={105\over16}{h^6\over U_d^6}
   \sum_cT_c\sum_{v\in c}j_{T,v}^T\Theta_{v,c}.}       \tag{CJ19b}
\]

No equality for the complete first-source vertex at order six is asserted:
the diagonal order-six vertex and the `A1/E` off-diagonal order-six pieces
remain unclassified.  The important composition fact is not an equality
between a source and a read.  Both displayed projections are derivatives of
the same microscopic Hamiltonian with respect to the same `j_v`, separated
only after the calculation by exact orthogonal projectors and
diagonal/off-diagonal operator support.

## 5. Relation to the earlier EW/AV typed split

The mutable `EW` model showed that six pair coordinates can span all of
`Sym^2(V)` as deformations of a Fisher tensor in a separate exponential
family.  `GL6CJ` neither imports that mutable theorem nor identifies its
Fisher tensor with space.  Instead, it establishes an exact six-direction
operator realization inside the actual F3 perturbative parent.

The sealed `GL6AV` theorem correctly refused to identify its formal
`A1+T2` orientation coefficient with its authenticated `E` pair read.  Those
were differently typed variables, and no common source Jacobian or rank-six
right inverse existed.  Equations (CJ06)--(CJ19b) change that precise status:

\[
 \boxed{
 \text{one pre-Feshbach pair source }j_v
 \longrightarrow
 \begin{cases}
 \text{locked diagonal }A_1+E\text{ read},\\
 \text{order-six off-diagonal }T_2\text{ writer},
 \end{cases}
 \quad\operatorname{rank}=3+3=6.}                        \tag{CJ20}
\]

Thus the former **algebraic source/read type split is closed at the
operator-jet level**: no post-hoc identification of two unrelated source
coordinates is required.  The stronger `AV-CONSTITUTIVE` and `AV-UPDATE`
gates remain open.  A laboratory query/control source is not thereby an
autonomous bulk field; a rank-six operator derivative is not a reciprocal
stationary response; and neither is yet a metric.

## 6. Exact disposition

The degree-two F3 parent now has one six-coordinate microscopic pair source
whose locked diagonal and sixth-order off-diagonal projections are jointly
faithful.  This closes the local source-access architecture that was missing
from the earlier typed atlas.  The next physics gate is to make the same
field endogenous to retained formation and evaluate both operator slots in
one stationary accumulated state.  Only after that can common-cone and
Ricci-form tests be meaningful.

`PASS__GL6CJ_ONE_PREFESHBACH_SIX_PAIR_SOURCE__SIX_LOCKED_WORDS_THREE_COMPLEMENT_CLASSES__DIAGONAL_A1_PLUS_E_RANK3_NORMAL_4PA_PLUS16PE__EVERY_Q4_NODE_TWELVE_HEXAGONS_EACH_LOCAL_PAIR_TWICE__H6_T2_WRITER_RANK3_GRAM_8PT__DIAGONAL_OFFDIAGONAL_ORTHOGONAL_OPERATOR_SUPPORT__COMBINED_RANK6_NORMAL_4PA_PLUS16PE_PLUS8PT__EXACT_RECONSTRUCTION_INVERSES_WITH_TYPED_IDENTITIES__EW_AV_TYPED_SPLIT_CLOSED_ONLY_AT_OPERATOR_JET_LEVEL__NO_AUTONOMOUS_SELF_GENERATION_STATIONARY_RESPONSE_METRIC_RGRLB_RICCI_GRAVITY_OR_G`
