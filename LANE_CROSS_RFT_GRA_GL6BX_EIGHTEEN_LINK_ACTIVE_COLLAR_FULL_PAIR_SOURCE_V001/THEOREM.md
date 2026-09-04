# Eighteen-link active-collar source-first effective-response theorem

**Packet:** `GL6BX`

**Claim class:** exact finite-collar geometry; exact source-before-projection
canonical Hermitian effective Hamiltonian and source jets through `h^6`; exact
retained-doublet Kubo response and eliminated-`Q` contact; exact response of
an explicitly selected equal one-per-orientation-class finite composition;
exact contact-polynomial and truncated-static root exclusions.

**Not claimed:** a parent-selected bulk phase; a full-parent or all-frequency
CTP response; a connected-cluster accumulation; a Ward identity; Ricci or
Einstein dynamics; a graviton; gravity; a cell threshold; or `G`.

## 1. Result

The prospective eighteen-link collar requested by `GL6BV/GL6BW` is a
well-typed finite regulator.  Its exact locked subspace has dimension two and
is precisely the two alternating configurations of the canonical `Q4`
hexagon.  With the full six-pair source inserted on all six core constraint
nodes before elimination, the canonical source-free effective Hamiltonian is

\[
 {H_{\rm eff}\over U_d}=
 \left(-9r^2-{11\over4}r^4-{169\over90}r^6\right)I
 -{63\over8}r^6\sigma_x+O(r^8),\qquad r={h\over U_d}.       \tag{BX01}
\]

The order-six off-diagonal coefficient is the exact `GL6AO` alternating-
hexagon coefficient.  Thus the normalized stationary branches are
`|+>` and `|->`, with `|+>` the ground branch and retained gap

\[
 \Delta={63\over4}U_dr^6                                      \tag{BX02}
\]

in the displayed `h<=6` functional.

A nonzero low-energy `T2` first vertex does appear.  It starts at `h^4`, but
it is proportional only to `I` and `sigma_x`; hence it is diagonal in the
actual stationary basis and has no retained-doublet transition matrix
element.  The `E2` source supplies the noncommuting transition.

For the prospectively imposed equal one-per-class composition, the exact
ground-branch **energy-Hessian contact** at spatial order `a_*^2 k^2` is, in
units `1/U_d`,

\[
\begin{aligned}
 h_A^{\rm ct}&=4r^2+{662\over27}r^4
                 +{235186909\over364500}r^6,\\
 h_E^{\rm ct}&=4r^2+{62\over27}r^4
                 -{1902787001\over364500}r^6,\\
 h_T^{\rm ct}&=4r^2-14r^4+{1627679\over4860}r^6.          \tag{BX03}
\end{aligned}
\]

The connected CTP local-in-time contact is `-C delta(t)`, so its three signs
are the opposites of (BX03).  The solder-normalized contact defect is

\[
 D_{\rm ct}(r):={h_E^{\rm ct}\over2}-h_T^{\rm ct}
 =r^2\left[-2+{409\over27}r^2
 -{2146938851\over729000}r^4\right].                       \tag{BX04}
\]

It is strictly negative for every `r>0`.  Therefore even this selected
finite collar cannot be tuned onto `h_E/2=h_T` within the complete contact
polynomial through `h^6`.

The retained-doublet connected Kubo piece is

\[
 K^R_{\rm dyn}(t,k)={\theta(t)\over\hbar}
 R(k)\sin\!\left({\Delta t\over\hbar}\right).             \tag{BX05}
\]

This uses the inherited packet convention

\[
 K^R(t)={i\over2\hbar}\theta(t)
 \langle[V(t),V(0)]\rangle,                                \tag{BX05a}
\]

which is one-half the standard raw susceptibility for a source Hamiltonian
`H+jV`.  In this convention the stationary relation is
`F''=-2K^R(0)`; the same normalization fixes the sign and weight of the
local contact.

At `k=0`, `R` vanishes in all three pair sectors.  At `a_*^2k^2`,

\[
 R_E(r)=128-256r^2-{3424\over9}r^4
              -{1621288\over2025}r^6+O_{\rm trunc}(r^8),
 \qquad R_A=R_T=0.                                         \tag{BX06}
\]

Only the leading small-time physical coefficient is complete at the current
Hamiltonian order:

\[
 K^R_{E}(t,k)={2016\,U_dr^6t\over\hbar^2}
 a_*^2k^2P_E+O(r^8t,t^3),
 \qquad D^R_{\rm dyn}(t,k)={1008\,U_dr^6t\over\hbar^2}
 a_*^2k^2.                                                  \tag{BX07}
\]

The energy/contact and dynamical pieces have different temporal support and
must not be equated coefficientwise.  The zero-momentum Kubo term vanishes,
whereas the order-`h^2` energy contacts are

\[
 (h_A,h_E,h_T)_{k=0}=(-24,-72,-56)r^2/U_d.                 \tag{BX08}
\]

Hence the declared low-energy CTP assembly does **not** cancel its `k^0`
contact.  No Ward identity is present.

## 2. Exact collar and source custody

Use the `L=4` `Q4` incidence witness and canonical links

\[
 (x,0),(x+e_0-e_1,1),(x+e_0-e_1,2),
 (x+e_2-e_1,0),(x+e_2-e_1,1),(x,2).                       \tag{BX09}
\]

The active microscopic set is the union of every link incident on a core
hexagon node.  Exact enumeration gives:

* six cycle links and twelve spokes, hence eighteen active links;
* six core constraint nodes, each with all four incident flips active;
* twelve boundary constraint nodes, each with one active link and three
  exterior links frozen;
* six boundary nodes with one frozen occupied link and six with two;
* exactly two fully locked active words, whose symmetric difference is the
  six-link cycle.

All other `Q4` links are frozen to the authenticated deterministic witness.
This is a selected finite collar/regulator, not a state or phase chosen by
the full parent.

At every core node use pair order

\[
 {\cal P}=(01,02,03,12,13,23),\qquad M_{ab}=Z_aZ_b.        \tag{BX10}
\]

The microscopic source rule is

\[
 H(j)=U_d\sum_{v\in V_{\rm affected}}q_v^2
      -h\sum_{e\in E_{18}}X_e
      +\sum_{v\in V_{\rm core}}j_v^TM_v.                  \tag{BX11}
\]

Thus there are thirty-six independent core source coordinates.  The twelve
boundary constraints remain in `H_0`, but their source is explicitly fixed
to zero.  `GL6BV`'s four-incident-edge formula is used only at the six core
nodes; it is never attached to a one-active-edge boundary node.

## 3. Canonical Hermitian construction

Let `P` project onto the two exact locked words and `Q=1-P`.  In the
`P+Q` block decomposition write

\[
 H=\begin{pmatrix}A&B\\C&D\end{pmatrix},\qquad
 \chi: P\longrightarrow Q.                                \tag{BX12}
\]

The exact Bloch wave operator is solved order by order from

\[
 C+D\chi=\chi(A+B\chi).                                    \tag{BX13}
\]

The non-Hermitian Bloch representative and its metric are

\[
 H_B=A+B\chi,\qquad G=I+\chi^\dagger\chi.                  \tag{BX14}
\]

The canonical des-Cloizeaux representative used everywhere in this packet
is

\[
 \boxed{H_C=G^{1/2}H_BG^{-1/2}}.                            \tag{BX15}
\]

Differentiating (BX13)--(BX15) once and twice with respect to the microscopic
sources includes every active direct history, reducible/folded term,
source-induced `P`-space energy, and normalization term through `h^6`.
All arithmetic is rational.  Hermiticity is checked at every retained order.

At `h^2`, a pure-`T2` onsite source gives

\[
 D^2H_C[j,j]=-2\|j\|^2r^2/U_d,                             \tag{BX16}
\]

exactly reproducing `GL6BV` at a core node.  The restriction `P_Tj=j` in
that predecessor is respected: (BX16) is not incorrectly applied to `A1`
or `E2`, whose direct locked-space energy is retained by (BX11)--(BX15).

## 4. First vertices and stationary basis

Use the rational orthogonal pair basis

\[
\begin{aligned}
 A&=(1,1,1,1,1,1),\\
 E_a&=(1,1,-2,-2,1,1),\qquad
 E_b=(1,-1,0,0,-1,1),\\
 T_1&=(1,0,0,0,0,-1),\quad
 T_2=(0,1,0,0,-1,0),\quad
 T_3=(0,0,1,-1,0,0).                                     \tag{BX17}
\end{aligned}
\]

In the sorted core-node order

\[
 C_{001},C_{100},C_{1,-1,1},P_{000},P_{0,-1,1},P_{1,-1,0},
\]

the only nonzero `T2` directions are

\[
 (-T_3,+T_1,+T_2,+T_2,+T_1,-T_3).                         \tag{BX18}
\]

For sign `s=+-1`, each corresponding first vertex is

\[
 V_T=s\left[{79\over36}r^4I+r^6
 \left({33307\over6480}I+{105\over8}\sigma_x\right)\right]. \tag{BX19}
\]

All other local `T` directions vanish through this order.  Equation (BX19)
is nonzero but commutes with (BX01).  It changes stationary energies and
contacts without producing a `|+> <-> |->` pole.

The spatially uniform selected-orientation first vertices are

\[
\begin{aligned}
 V_A&=(-12+12r^2+47r^4/6+86369r^6/5400)I
       +(315/4)r^6\sigma_x,\\
 V_{E_a}&=V_{E_b}=0,\\
 V_{T_1}&=V_{T_2}={79\over18}r^4I
       +r^6\left({33307\over3240}I+{105\over4}\sigma_x\right),\\
 V_{T_3}&=-V_{T_1}.                                      \tag{BX20}
\end{aligned}
\]

## 5. Contacts, `E/T` mixing, and selected orientation completion

Before orientation completion the selected collar has a nonzero `E/T`
contact.  In bases `(E_a,E_b)` and `(T_1,T_2,T_3)`, its ground-branch
bilinear block is a scalar times

\[
 B_{ET}=\begin{pmatrix}1&1&2\\1&-1&0\end{pmatrix}.          \tag{BX21}
\]

At `k=0` the scalar through orders `(r^2,r^4,r^6)` is

\[
 \left(8,{937\over27},{21054127\over48600}\right),         \tag{BX22}
\]

and at normalized `a_*^2k^2` it is

\[
 \left(-4,-{1003\over27},-{31370461\over72900}\right).     \tag{BX23}
\]

The script constructs the literal `S4` orbit of the six links.  It finds
four translation-inequivalent orientation classes with stabilizer six.  At
each of `k^0,k^2` and `h^2,h^4,h^6`, (BX21) is invariant under the canonical
six-member stabilizer.  Therefore one literal representative from each of
the four classes is well defined, and their transformed equal-weight sum
kills the complete `E/T` block.  The choice of one copy from each class with
equal multiplicity is imposed prospectively for this selected finite
composition: it is not generated or weighted by the parent and is not a
homogeneous accumulation.  No abstract Schur cancellation is assumed in
place of the explicit orbit check.

The spatial normalization ledger is

\[
 R_{E,P}=8,\qquad R_{E,C}=8,\qquad
 \nabla A_P=\nabla A_C,\qquad
 R_{E,P+C}=32,\qquad R_{E,4\,selected\ classes}=128.        \tag{BX24}
\]

The factor from coherent parent-plus-child amplitudes is therefore distinct
from the imposed multiplicity four of the selected orientation classes.

For completeness, the selected equal-weight four-class `k=0` energy contacts
are

\[
\begin{aligned}
 h_A^{(0)}&=-24r^2-{812\over27}r^4-{90246287\over121500}r^6,\\
 h_E^{(0)}&=-72r^2-{11168\over27}r^4-{40866533\over12150}r^6,\\
 h_T^{(0)}&=-56r^2-{284\over9}r^4-{3735793\over4860}r^6.   \tag{BX25}
\end{aligned}
\]

Equations (BX03) and (BX25) include the expectation of the second canonical
source derivative in `|+>`.  They are energy-Hessian contacts; the connected
CTP delta-contact has the opposite sign.

## 6. Branchwise CTP and static energy

For a real symmetric first vertex `V`,

\[
 V_{+-}={V_{00}-V_{11}\over2}.                              \tag{BX26}
\]

Only `E2` has nonzero (BX26).  Its selected equal-weight four-class Fourier
residue gives (BX06).  The terms displayed above `r^6` by the executable are
only algebraic products of first vertices truncated at `r^6`; they are not
complete higher-order parent coefficients.

The eliminated-`Q` second derivative is local in time in the low-energy
functional.  Restoring the high-gap microscopic time dependence would
replace/resum this contact; it must not be added as a second independent
effect.  Consequently (BX05) plus `-C delta(t)` is a low-energy effective
CTP assembly, not an exact all-frequency response.

The stationary ground-energy Hessian of the same truncated two-level
functional is

\[
 H^{\rm stat}_{AB}=C_{AB}-{2V^A_{+-}V^B_{-+}\over\Delta}.   \tag{BX27}
\]

At `a_*^2k^2`, its solder defect is the Laurent polynomial

\[
\begin{aligned}
 D_6(r)={}&-{512\over63}r^{-6}+{1024\over63}r^{-4}
 +{13696\over567}r^{-2}+{6485152\over127575}\\
 &-{4982834\over42525}r^2
 -{172145767\over1148175}r^4
 -{6525772091869\over2066715000}r^6.                      \tag{BX28}
\end{aligned}
\]

Only the leading `-512 r^-6/63` spectral coefficient is protected from
unknown `gamma_8+` and gap corrections.  The remaining Laurent coefficients
are exact for the **summed `h<=6` truncated two-level functional**, not for
the untruncated parent.

## 7. Exact root exclusions

For the contact polynomial, put `x=r^2`.  The bracket in (BX04) has

\[
 \operatorname{disc}=-{236225414\over10125}<0              \tag{BX29}
\]

and negative leading coefficient.  Its maximum occurs at

\[
 x={5521500\over2146938851},\qquad
 D_{\rm ct}(r)/r^2=-{4252057452\over2146938851}<0.          \tag{BX30}
\]

Thus (BX04) has no positive root.

Multiplying (BX28) by `r^6` gives `Q_6(x)` with primitive ascending integer
coefficients

\[
 (-16796160000,33592320000,49921920000,105059462400,
 -242165732400,-309862380600,-6525772091869).               \tag{BX31}
\]

The exact Sturm endpoint sign lists are

\[
 (-,+,+,-,+,+,+)\quad(x=0^+),\qquad
 (-,-,+,+,-,-,+)\quad(x=+\infty).                          \tag{BX32}
\]

Both have three variations.  Hence `Q_6` has zero positive roots, and since
`Q_6(0)<0`,

\[
 \boxed{D_6(r)<0\quad\hbox{for every }r>0}                  \tag{BX33}
\]

inside the truncated functional.  A conservative Weyl bound
`||sum_e X_e||<=18` keeps the doublet separated from the source-off
`2U_d` defect threshold for `0<r<1/18`.  There is no root in that rigorously
controlled domain because there is no positive root at all.  Equation
(BX33) is not a convergence or Ricci theorem for the untruncated parent.

## 8. Conclusion and next route

The eighteen-link regulator closes the finite-coupling representation gap:
virtual defects generate a lawful static/contact `T2` first vertex in the
same canonical effective functional that carries the hexagon `E2`
transition.  They do not generate a `T2` doublet commutator through `h^6`.
The selected equal-weight four-class completion is diagonal in `A1/E2/T2`,
but its like-coefficient solder test is off ray, its `k^0` contact does not
cancel, and neither the contact polynomial nor the summed truncated
stationary surrogate has a positive tuning threshold.

Therefore single-collar tuning is closed at this order.  The next physical
route is connected accumulation of authenticated collars/owners under one
common source and boundary rule—not a fitted isotropic weight, an inserted
Ricci term, or an inference of gravity.

`PASS__GL6BX_SELECTED_EIGHTEEN_LINK_ACTIVE_COLLAR__EXACT_TWO_WORD_LOCK__FULL_SIX_PAIR_CORE_SOURCE_BEFORE_PROJECTION__CANONICAL_HERMITIAN_H2_H4_H6__AO_MINUS_63_OVER_8__NONZERO_COMMUTING_T2_FIRST_VERTEX_AT_H4__ALL_SECOND_CONTACTS_AND_SELECTED_ET_MIXING__EXPLICIT_STABILIZER_INVARIANT_EQUAL_ONE_PER_CLASS_SUM__PARENT_CHILD_8_8_COHERENT_32_SELECTED_CLASS_MULTIPLICITY_128__NOT_PARENT_WEIGHTED_OR_HOMOGENEOUS_ACCUMULATION__K0_CONTACT_NOT_CANCELLED__CONTACT_SOLDER_DEFECT_STRICTLY_NEGATIVE__SUMMED_H_LE_6_STATIC_TWO_LEVEL_SURROGATE_STURM_ZERO_POSITIVE_ROOTS__NO_FULL_PARENT_CTP_RICCI_GRAVITY_G_OR_THRESHOLD_CLAIM`
