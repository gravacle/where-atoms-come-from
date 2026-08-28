# Coulomb-DPAR projected source-rank closure on the F3/q4 diamond-ice family

**Lane:** `LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001`  
**Date:** 2026-08-27  
**Status:** exact conditional theorem on the frozen FU physical completion,
the additional `FV-PURE` source premise, and the covering-matched q4
diamond-ice family  
**Claim class:** exact direct ice-pair source; exact differentiated H6
Hermitian ring source from all 720 histories in every orientation; exact
off-shell projected operator rank; formal through-order-eight consequence

## 1. Question and dependencies

FU derived, subject to its complete physical-solder conditions `S1`--`S9`, a
grounded, charge-conserving, noncircular pair-resolved Coulomb contribution to
the q4 degree square.  It fixed the ideal Coulomb slope

\[
 \lambda=-\frac12
\]

and required the physical source to be inserted before the fixed incidence
Feshbach reduction while the FS one-edge affine source remains unchanged.
FU proved microscopic rank six but correctly left the projected rank open.
Its `S1`--`S9` conditions still permit residual local mutual kernels and
ranked covariant remainders.  They do not by themselves identify the complete
source derivative with the pure ideal-Coulomb law used below.
FT established that the direct pair source restricts to `A1+E` on one ice
fiber, with its `A1` part scalar and `T2` null.  CW supplied the exact 720-path
H6 coefficient, FS fixed the one-edge source, and FM classified all projected
operator endpoints through H8.

This lane asks the remaining calculation:

> Does the same complete source, inserted before Feshbach, have six
> independent **nonidentity projected operators**, rather than merely six
> microscopic coefficient tensors?

The answer is yes for the FU Coulomb contribution under the additional
`FV-PURE` premise below.  This is a conditional theorem because it inherits
FU's physical-solder premises and adds a complete-source identification; it
is not an independent derivation of either premise set.

## 2. Frozen source and conventions

In addition to FU `S1`--`S9`, impose the following explicit premise.

**`S10 / FV-PURE` (complete pure-source identification).** On the frozen
scored subspace, the complete nonidentity first derivative of the
pre-Feshbach source is exactly the pair-resolved ideal-Coulomb derivative in
(FV01) plus the unchanged FS one-edge flip derivative.  Every residual local
pair term either vanishes or coalesces into the same normalized radial law;
cross-node, boundary, controller, and other remainder derivatives are absent
or identity/reference terms; and the dressed flip source is exactly the
inherited FS source.  No additional nonidentity source operator is present.

`FV-PURE` is stronger than FU `S1`--`S9`.  It is introduced here rather than
silently inferred from FU.  If an allowed FU remainder has a different strain
law, the exact coefficients, cancellation slopes, and rank theorem below
must be recomputed for that completed source.

Let

\[
 D_a=n_an_a^{\mathsf T},\qquad
 \widehat R_{ab}=
 \frac{(n_b-n_a)(n_b-n_a)^{\mathsf T}}{|n_b-n_a|^2},
 \qquad P_{v,ab}=Z_{v,a}Z_{v,b}.
\]

The complete linear source used here is

\[
\begin{split}
 H[j]={}&H_C[j]
 -h\sum_e\left(1-\frac12j:D_{a(e)}\right)\widetilde X_e
 +O(j^2),\\
 H_C[j]={}&H_C[0]-\frac{U_d\lambda}{2}
 \sum_v\sum_{a<b}(j:\widehat R_{ab})P_{v,ab}
 +H_{\rm id}[j]+O(j^2),                         \tag{FV01}
\end{split}
\]

where `H_C[0]=U_d sum_v(d_v-2)^2`, the dressed flips are exactly equivalent to
the inherited flips on FU's frozen scored subspace, and

\[
 Q^{ij}:=-2\frac{\partial H}{\partial j_{ij}}\bigg|_{j=0}.    \tag{FV02}
\]

Self-elastance and reference terms in `H_id` are retained.  They are identity
operators.  They cancel from every virtual energy difference and identity
shifts do not count toward nonidentity rank.

All symmetric-tensor coordinates below use
`(xx,yy,zz,2xy,2xz,2yz)`.  The physical tetrahedral normalization is retained:
`D_a` contains the factor `1/3`, while every normalized root dyad contains the
factor `1/8` in the unnormalized sign-vector coordinates.

## 3. Direct projected pair source

The order-zero projected pair conjugate is

\[
 Q_{\rm pair}^{(0),ij}
 =U_d\lambda\sum_v\sum_{a<b}\widehat R_{ab}^{ij}P_{v,ab}
 +Q_{\rm id}^{(0),ij}.                            \tag{FV03}
\]

On the six local ice states, its `A1` value is common and hence an identity;
its centered nonidentity image is exactly the two-dimensional `E` irrep, and
its `T2` image vanishes.  At `lambda=-1/2`, normalized diagonal expectation
differences between the actual direction-pair ice coverings `(01)`, `(02)`,
and `(03)` give

\[
 e_1=(-1,1,0,0,0,0),\qquad
 e_2=(-1,0,1,0,0,0).                              \tag{FV04}
\]

These are matrix-element differences, so any identity contribution cancels.
They prove nonidentity **operator** rank two in `E`, not merely coefficient
rank two.

## 4. Differentiating every H6 history

Fix a flippable elementary hexagon `C` in an ice state `n`.  Its three link
labels each occur twice; let `d` be the missing tetrahedral label.  For a
proper prefix subset `S` of toggled cycle links, use the endpoint-referenced
gap

\[
 \Delta_n(S;j)=E_{n\triangle S}(j)-E_n(j).        \tag{FV05}
\]

It is essential to subtract the derivative of the appropriate endpoint ice
energy.  At source off, write `Delta_n(S;0)=U_d delta_n(S)`.  The pair part of
the exact gap derivative is

\[
 \frac{\partial\Delta_n(S;j)}{\partial j}\bigg|_0
 =-\frac{U_d\lambda}{2}\sum_{v,a<b}\widehat R_{ab}
 \left[Z^{n\triangle S}_{v,a}Z^{n\triangle S}_{v,b}
       -Z^n_{v,a}Z^n_{v,b}\right].                \tag{FV06}
\]

All self/reference identities cancel in (FV06).  For an order
`pi in S_6`, let `S_r(pi)` be its first `r` links and

\[
 w_\pi=\prod_{r=1}^{5}\frac1{\delta_n(S_r(\pi))}.
\]

The differentiated forward matrix element, with the common factor
`-h^6/U_d^5` removed, has tensor

\[
 T^\rightarrow_{C,n}(\lambda)=
 \sum_{\pi\in S_6}w_\pi\left[
 \sum_{e\in C}D_{a(e)}
 +2\sum_{r=1}^{5}
 \frac{\partial_j\Delta_n(S_r(\pi);j)|_0}
      {U_d\,\delta_n(S_r(\pi))}
 \right].                                         \tag{FV07}
\]

The first term differentiates all six flip numerators.  The second
differentiates all `720 x 5=3600` resolvents.  Exact enumeration reproduces
the CW gap classes with multiplicities `96,144,216,192,72` and

\[
 \sum_{\pi\in S_6}w_\pi=\frac{63}{8}.             \tag{FV08}
\]

A forward Bloch entry alone need not be Hermitian once the source splits the
ice endpoints.  The projected conjugate operator is therefore defined by the
standard Hermitian forward/reverse average

\[
 \overline T_{C,n}(\lambda)=\frac12
 \left(T^\rightarrow_{C,n}(\lambda)
      +T^\leftarrow_{C,n\triangle C}(\lambda)\right).          \tag{FV09}
\]

The reverse term uses (FV05) with `n triangle C` as its initial endpoint; it
is not obtained by reusing the forward initial energy.

### Theorem `CPSR-1` -- exact Hermitian H6 tensor

For every elementary hexagon orientation and every one of the `2^6` locally
allowed assignments of the six external occupied links,

\[
 \boxed{
 \overline T_{C,n}(\lambda)=
 \frac{21}{8}(8-15\lambda)I
 -\frac{63}{8}(2-5\lambda)D_d .}                 \tag{FV10}
\]

#### Proof

Every proper prefix has positive defect gap, so all 720 orders are
irreducible.  Evaluate (FV06)--(FV09) exactly.  The forward tensor equals the
right side of (FV10) plus an endpoint/environment-dependent diagonal
traceless `E` tensor.  Complementing the six cycle occupations maps it to the
reverse history with the same denominator weight and the negative `E`
residue.  Thus the Hermitian average cancels that residue.  The remaining
trace and off-diagonal coordinates are independent of the external choices.
Exact exhaustive evaluation of all `4 x 64` labeled local environments gives
the displayed rational coefficients.  Because the derivative is affine in
`lambda`, evaluation of its constant and linear parts proves (FV10) for every
real `lambda`.  QED.

No environment-averaging assumption is used.  The cancellation is between
the two endpoint-referenced orientations of the same Hermitian operator
matrix element.

At the ideal Coulomb slope,

\[
 \overline T_{C,n}\!\left(-\frac12\right)
 =\frac{651}{16}I-\frac{567}{16}D_d.              \tag{FV11}
\]

For missing labels `d=0,1,2,3`, its four coordinate rows are

\[
\begin{array}{c|rrrrrr}
d&xx&yy&zz&2xy&2xz&2yz\\ \hline
0&231/8&231/8&231/8&-189/8&-189/8&-189/8\\
1&231/8&231/8&231/8& 189/8& 189/8&-189/8\\
2&231/8&231/8&231/8& 189/8&-189/8& 189/8\\
3&231/8&231/8&231/8&-189/8& 189/8& 189/8
\end{array}                                        \tag{FV12}
\]

These rows have exact rank four and tetrahedral content `A1+T2`.  Their `A1`
is not an identity: it multiplies off-diagonal ring transitions
`|n triangle C><n|+h.c.`.  Their `T2` distinguishes the four missing-label
orientations.

The exact special slopes are also exposed.  `lambda=2/5` cancels `T2` and
leaves ring rank one; `lambda=3/5` cancels ring `A1` and leaves rank three;
`lambda=0` retains ring rank four but removes the direct pair `E`.  Coulomb
`lambda=-1/2` is none of these values.

## 5. Folds, diagonal terms, and the complete leading-H6 source

Every order-two and order-four projected endpoint is diagonal because the
covering-matched `G_5` graph has no two- or four-link cycle.  At order six,
an off-diagonal ice endpoint must be one elementary hexagon with each link
flipped once.  Every proper prefix is in `Q_2`, so no Feshbach fold contributes
to that matrix element.  Every H6 fold is made from lower diagonal kernels
and their energy derivatives and remains diagonal.  Differentiating it can
make it a nonidentity diagonal source, but cannot change any ring entry.

Thus the exact rank-relevant fixed-`P_2` source through leading H6 is

\[
\begin{split}
 Q_{\rm eff}^{(\le6)}={}&Q_{\rm pair}^{(0)}
 +Q_{\rm diag}^{(2)}+Q_{\rm diag}^{(4)}+Q_{\rm diag}^{(6)}
 +Q_{\rm id}^{(\le6)}\\
 &-\frac{h^6}{U_d^5}\sum_{n,C\ {\rm flippable}}
 \overline T_{C,n}^{ij}
 \left(|n\triangle C\rangle\langle n|+|n\rangle\langle n\triangle C|\right)
 +O(h^8).                                          \tag{FV13}
\end{split}
\]

Here each `Q_diag^(m)` is the complete derivative of the irreducible diagonal
words and all self-consistency folds at that order; none is dropped.
`Q_id^(<=6)` contains the complete projected scalar/reference derivatives.
The diagonal and identity terms cannot contribute to the four selected
off-diagonal ring functionals.  Conversely, their higher formal powers of
`h` cannot cancel the order-zero diagonal `E` witnesses (FV04).  This is why
their exact values are unnecessary for the rank, while their operator support
and custody remain explicit.

## 6. Exact projected operator rank

On `G_5`, an exact integer-capacity bipartite completion constructs one global
degree-two ice state with a flippable hexagon of each missing-label type.
Switching each hexagon remains in the projected Hilbert space.  Therefore the
four rows (FV12) are four distinct off-diagonal **operator matrix elements**.
Pulling these states and marked transitions through the covering maps frozen
by FS gives the same four local witnesses on every `G_(5*2^r)`; the local
720-path coefficients are unchanged.  The direction-pair diagonal coverings
likewise pull back with degree two.
Together with the two normalized diagonal expectation differences (FV04),
the six evaluation functionals give

\[
 \det W=-\frac{4678629417}{256}\ne0.              \tag{FV14}
\]

Identity operators vanish on the diagonal differences and on every
off-diagonal functional.  Since a symmetric source has only six coordinates,
the nonzero determinant is both a lower and upper bound.

### Theorem `CPSR-2` -- Coulomb projected rank closure

Under the frozen FU `S1`--`S9` completion **and `S10 / FV-PURE`**, on the
covering-matched q4 diamond-ice family with `h != 0` as a formal expansion
parameter,

\[
 \boxed{\operatorname{rank}_{\rm nonid}
 D_jH_{\rm eff}^{(\le6)}\big|_{j=0}
 =2_E+1_{A_1}+3_{T_2}=6.}                        \tag{FV15}
\]

The `E` lower bound comes from the direct diagonal pair potential.  The
`A1+T2` lower bound comes from nonidentity Hermitian H6 ring operators.  No
identity shift is counted.

## 7. Formal H8 consequence and strict ceiling

FM proves that H8 contains only diagonal terms, dressed-H6 transitions, and
new octagon transitions.  Source-before-Feshbach differentiation preserves
that endpoint classification.  In the six witness functionals, the two
diagonal rows begin at order `h^0` and the four ring rows begin at order
`h^6`.  Hence the normalized determinant has the formal structure

\[
 \det W(h)=
 -\frac{4678629417}{256}h^{24}+O(h^{26})           \tag{FV16}
\]

after common nonzero `U_d`, volume, and sign factors are removed.  An H8 term
cannot cancel a distinct lower formal order.  Therefore the formal
through-order-eight rank remains six.

Equation (FV16) is not an unrestricted numerical finite-`h` statement.  It
does not provide a calibrated convergence radius and does not exclude a
tuned algebraic zero at larger finite `h` after higher orders are included.

Finally, (FV15) is an off-shell projected operator rank.  It does not
establish a retarded or CTP rank, a state-dependent spectral residue, or a
continuum tensor response.  It does not prove a Ward identity, tensor pole,
gravity, common-metric coupling, `RGRL-B`, or a value of `G`.  Those require
the separate dynamical/continuum tests already identified by FQ and FM.
