# Degree-pair geometric-strain underdetermination and one sufficient closure law

**Lane ID:** `GRA-FT-F3-Q4-DPGSB-V001`

**Short name:** `DPGSB`

**Date:** 2026-08-27

**Claim class:** exact same-source-off-source nonuniqueness theorem; exact
degree-pair and ice-pair representation audit; exact conditional microscopic
rank-six root-pair source construction; exact inherited order-six local
commutator boundary; one compact sufficient constitutive-law identification

**Status:**
`CURRENT_REDUCED_PARENT_DOES_NOT_DERIVE_E_STRAIN_COUPLING__FJ_FK_PAIR_OBSERVABLES_SUPPLY_AN_AVAILABLE_E_QUERY_NOT_AN_INHERITED_GEOMETRIC_SOURCE__FQ17A_FS_ADDITIVE_SOURCE_REMAINS_A1_PLUS_T2_RANK4__SAME_SOURCE_OFF_HAMILTONIAN_ADMITS_RANK4_AND_RANK6_MICROSCOPIC_SOURCE_EXTENSIONS__ROOT_PAIR_EXTENSION_IS_CONDITIONAL_AND_HAS_AN_INHERITED_LOCAL_H6_E_COMMUTATOR_CHANNEL__DPAR_IS_ONE_SUFFICIENT_UNADOPTED_CONSTITUTIVE_LAW_NOT_A_UNIQUENESS_THEOREM`

**Not claimed:** that an arbitrary external pair query is physical strain;
that the root/cross-dyad source is already inherited by F3; that equality at
`j=0` fixes a source derivative; that the full state-dependent CTP kernel is
known; that the conditional source has the RGRL-B Ward packet, a tensor pole,
gravity, or `G`.

## 1. Exact question and custody

`FR` and `FS` proved that the complete linear source of the selected periodic
CW/FM pure-incidence branch has exact rank four, `A1+T2`, with the two
diagonal-traceless `E` strain directions in its kernel.  `FJ` and `FK` had
already proved two facts which could appear to evade that obstruction:

1. the degree interaction contains physical link-pair operators; and
2. after ice projection their normalized pair variations are exactly an
   `E` representation with nonzero inherited ring commutators.

This lane asks whether those existing observables *derive* the missing `E`
part of the already frozen block-strain source while preserving all of the
following:

- the FR/FS one-edge source remains `A1+T2`;
- source deformation occurs before fixed Feshbach reduction;
- the source-off Hamiltonian is exactly unchanged; and
- no source-off interaction or post-result fitted tensor weight is inserted.

The dependency hashes are frozen in `DEPENDENCIES.sha256`.  The answer is a
sharp distinction: an `E` pair **query** is available, but the present parent
does not choose it as a geometric strain derivative.  A single additional
physical constitutive premise would do so without changing the source-off
Hamiltonian.

## 2. Three objects that must not be conflated

### 2.1 External pair query

For any already existing Hermitian pair observable `P_ab=Z_a Z_b`, one may
mathematically introduce a probe

\[
 H_{\rm probe}[j]=H-\frac12\sum_{a<b}(j:C_{ab})P_{ab}. \tag{FT01}
\]

This is a lawful generating-functional query when `C_ab`, its normalization,
the state, ports, and contacts are frozen prospectively.  It vanishes from
the source-off dynamics.  But `H_probe[0]=H` alone does not prove that `j` is
physical strain or that `C_ab` is the derivative of the F3 parent under a
geometric deformation.  It also does not make the resulting source a new
source-off interaction: a source query and a source-off Hamiltonian are
different objects.

### 2.2 Parent-derived geometric strain

A geometric strain source requires a physical one-parameter family
`F -> H(F)` whose value at `F=I` is the frozen parent and whose derivative is
fixed by how the *existing* interaction depends on lengths, angles, areas, or
other earned geometric data.  The BS06 incidence Hamiltonian explicitly
contains no distance or dimension.  Its degree penalty depends on occupation
and incidence only.  The q4 coframe supplies candidate geometric tensors, but
does not supply this missing constitutive dependence.

### 2.3 Unearned cross-dyad assignment

Writing a root or cross-dyad next to a pair operator after noticing that it
fills `E` is not a derivation.  It changes the source map even if it leaves
`H[j=0]` unchanged.  In particular, the existence and uniqueness up to scale
of an equivariant `E` intertwiner in `FK` do not determine its physical scale,
sign, or identification with the derivative of the degree energy.

## 3. Exact local algebra

Use the tetrahedral coframe

\[
 n_1={1\over\sqrt3}(1,1,1),\quad
 n_2={1\over\sqrt3}(1,-1,-1),\quad
 n_3={1\over\sqrt3}(-1,1,-1),\quad
 n_4={1\over\sqrt3}(-1,-1,1),                    \tag{FT02}
\]

and define

\[
 D_a=n_an_a^{\mathsf T},\qquad
 R_{ab}=(n_b-n_a)(n_b-n_a)^{\mathsf T}.           \tag{FT03}
\]

`FR/FS` give

\[
 \operatorname{span}\{D_a\}=A_1\oplus T_2,
 \quad\operatorname{rank}=4,
 \quad\ker\{j\mapsto(j:D_a)_a\}=E.              \tag{FT04}
\]

Every FQ17a additive weight on a two-link term is `D_a+D_b` and remains in
that four-dimensional span.  By contrast, the six root dyads obey

\[
 R_{ab}=D_a+D_b-n_an_b^{\mathsf T}-n_bn_a^{\mathsf T}           \tag{FT05}
\]

and span all of `Sym^2(R^3)` with rank six.  The last two cross terms are
exactly what the additive rule does not contain.

The degree penalty at one q4 vertex has the exact operator identity

\[
 U_d(d_v-2)^2
 =U_d I+{U_d\over2}\sum_{a<b}P_{ab},
 \qquad P_{ab}=Z_aZ_b.                             \tag{FT06}
\]

On the local ice fiber, complementary pairs agree and
`sum_(a<b) P_ab=-2I`.  Thus the centered pair space has dimension two and is
exactly `E`.  Contracting the six `R_ab` with either diagonal-traceless source
basis vector gives two independent centered functions on the six ice states;
contracting with any `T2` off-diagonal source gives zero after ice
restriction.  The root-pair map therefore has type `A1+E` in ice, precisely
in the sense that its nonconstant `E` part supplies the two directions absent
from the edge-weight span.  Its `A1` part overlaps the edge/degree scalar
sector and is not a second independent scalar.

### Theorem `DPGSB-1` -- observables are available but coupling is not derived

The selected parent already contains pair operators whose ice-restricted
nonconstant span is the missing `E`.  FJ/FK also prove that those operators
are dynamically nontrivial.  Neither fact changes the derivative of the
frozen FQ17a/FS source.  Operator availability, response availability, and a
parent-derived geometric source are three different statements.

## 4. Same-source-off nonuniqueness

Keep the FR/FS one-edge deformation unchanged.  There are at least two
prospective microscopic extensions of the degree term which agree exactly at
source off.

The frozen FS extension keeps the degree square unsplit and assigns it the
scalar support sum `W_v=sum_a D_a=(4/3)I`,

\[
 H_{\rm deg}^{\rm FS}[j]
 =\left[1-{1\over2}j:W_v\right]U_d(d_v-2)^2+O(j^2). \tag{FT07}
\]

Its union with the edge weights has rank four and both `E` contractions
vanish.  Even if the degree square is prospectively expanded while retaining
FQ17a additive custody, each pair receives only `D_a+D_b`; that alternative
bookkeeping remains in the same rank-four span.

A root-pair external query instead uses the normalized sibling dyad
`Rhat_ab=R_ab/|n_b-n_a|^2`:

\[
 H_{\rm deg}^{\rm root}[j]
 =U_d I+{U_d\over2}\sum_{a<b}
 \left[1-\lambda j:\widehat R_{ab}\right]P_{ab}+O(j^2),
 \qquad\lambda\ne0.                               \tag{FT08}
\]

The edge terms are untouched.  Equations (FT07) and (FT08) satisfy

\[
 H_{\rm deg}^{\rm FS}[0]
 =H_{\rm deg}^{\rm root}[0]
 =U_d(d_v-2)^2                                   \tag{FT09}
\]

exactly, yet their microscopic linear source ranks are four and six
respectively.  Equation (FT08) **replaces** the FS scalar deformation of the
degree term in (FT07); it is not added on top of it.  The unchanged one-edge
terms remain present.  Thus the `A1` component of (FT08) is an alternative
degree-sector response, not a double-counted scalar.  In the full microscopic
operator algebra the six pair Pauli strings are independent and the six root
dyads already give exact rank six.

Consequently `DPAR` is not an extension *inside* the unchanged FQ17a/FS
degree-source rule.  If adopted, it prospectively revises that rule for the
degree-pair sector while leaving the one-edge source and `H[0]` fixed.  This
is exactly why it is additional physical input rather than a consequence of
the additive source custody.

After ice projection the direct root-pair image has only `A1+E`, hence rank
three as a diagonal operator map.  The current chain proves only an upper
bound for the projected edge/Feshbach image, not exact projected `T2` rank.
Accordingly the exact rank-six statement in this lane is microscopic and
source-before-Feshbach; it is not a claim that the complete ice-projected or
CTP source already has rank six.

Both deformations can be inserted before the same fixed Feshbach map and both
leave every source-off H6/H8 coefficient unchanged.  Feshbach ordering
therefore cannot decide between them; the first derivative must be fixed at
the microscopic parent.

### Theorem `DPGSB-2` -- sharp underdetermination

The source-off reduced Hamiltonian, the degree-term identity, the FJ link-pair
operators, the FK ice-pair representation, and the FM H6/H8 histories do not
uniquely determine a geometric strain derivative.  Two local, Hermitian,
same-parent source families with identical `H[0]` have different exact linear
ranks.  Consequently the current inherited data do **not** lawfully promote
(FT08) from an external query to physical strain.

This is stronger than saying that a calculation is missing: `H[0]` is
mathematically insufficient to select `D_jH|_0`.

## 5. One explicit sufficient physical premise

One compact sufficient premise is:

> **Degree-pair affine-response law (`DPAR`).** The existing pair-resolved
> degree energy is a real differentiable physical function of the affinely
> strained sibling separation, with one common prospectively fixed law and a
> nonzero, independently derived or calibrated derivative at the tetrahedral
> point.

One explicit normal form is

\[
 {U_d\over2}P_{ab}
 \longmapsto
 {U_d\over2}
 g\!\left({|F(n_b-n_a)|^2\over|n_b-n_a|^2}\right)P_{ab},
 \quad g(1)=1,\quad g'(1)=\lambda\ne0,             \tag{FT10}
\]

with real symmetric `j` and `F(j)=I-j/2+O(j^2)`.  Since

\[
 { |F(j)r_{ab}|^2\over |r_{ab}|^2}
 =1-j:\widehat R_{ab}+O(j^2),                    \tag{FT10a}
\]

differentiating (FT10), rather than assigning a tensor by inspection, yields
(FT08).  With the established convention
`Q^(ij)=-2 partial H/partial j_(ij)|_0`, it gives

\[
 Q^{ij}_{\rm pair}=U_d\lambda\sum_{a<b}
 \widehat R_{ab}^{ij}P_{ab}.                     \tag{FT10b}
\]

This checks the factor of two and fixes the sign through `g'(1)`.  Real `g`
and real `j` preserve Hermiticity.  The common law is `S4` covariant because
permuting tetrahedral labels orthogonally permutes the six normalized roots.
The equality `g(1)=1` proves that no source-off interaction is added.

`DPAR` is one compact law, not six fitted weights.  Once the coframe, affine
map, exact `S4` covariance, and source convention are held fixed, its linear
content is one common nonzero slope.  But this lane does **not** prove that
`DPAR` is the unique or logically necessary closure: another prospectively
derived covariant constitutive family could also supply `E`.  `DPAR` must be
derived from a physical realization of the degree energy, adopted explicitly,
or calibrated independently before the response is inspected.  It is neither
inherited nor adopted here; the current BS/FJ/FK/FM/FQ/FR/FS chain does not
contain it.

### Theorem `DPGSB-3` -- conditional closure under `DPAR`

If `DPAR` is separately adopted, derived, or independently calibrated with
`lambda != 0`, the source-before-Feshbach microscopic block has exact
off-shell rank six while retaining the FR/FS one-edge `A1+T2` component and
exactly the same source-off Hamiltonian.  Relative to the previous four-rank
microscopic source, its two new directions are carried on the already
existing centered ice-pair `E` operators; no new source-off term is required.

This conditional theorem closes source **rank**, not the CTP pole, Ward, or
gravity tests.

## 6. Exact source-off commutator boundary on the frozen periodic family

Use the FS covering-matched periodic family `G_L`, `L=5*2^r`, and the inherited
leading ice Hamiltonian

\[
 H_6=E_{\rm scalar}P_2-J_6\sum_C B_C,
 \qquad J_6={63h^6\over8U_d^5}>0.                 \tag{FT11}
\]

Under `DPAR`, an `E` source direction `e` localized at one vertex block `v`
has the leading projected microscopic source-off conjugate

\[
 Q_{E,v}(e)=U_d\lambda\sum_{a<b}
 (e:\widehat R_{ab})P_{ab}(v).                   \tag{FT12}
\]

For every ring `C`, a pair is odd precisely when the ring crosses exactly one
of its two incident links.  The FK parity identity gives the exact operator
formula

\[
 [H_6,Q_{E,v}(e)]
 =2J_6U_d\lambda\sum_{C\ni v}\sum_{a<b:\,P_{ab}(v)\ \mathrm{odd\ on}\ C}
 (e:\widehat R_{ab})P_{ab}(v)B_C.                \tag{FT13}
\]

At every vertex the periodic q4 family contains hexagons realizing all six
local ring-edge pair types.  Their odd-projection coefficient map has rank
two on `E`.  For every nonzero `e in E`, at least one such local ring type
therefore has a nonzero formal commutator.  To turn that local coefficient
into a nonzero matrix element of the flippability-projected global operator,
one must also supply a compatible global ice state in which the selected ring
is alternating, exactly as required by FK.  The independent hostile audit
constructs such a state on `G_5` for both alternating orientations of every
hexagon through the reference vertex.  Distinct ring symmetric differences
then prevent cross-ring cancellation by FK23.

Thus, **conditional on DPAR**, the new block-local linear `E` conjugates are
not static labels; the inherited source-off H6 dynamics already acts on them
at formal leading order.  The full Feshbach source conjugate may also contain
generated higher-order corrections; (FT12)--(FT13) isolate its nonzero
order-six coefficient.  FM's order-eight terms remain in the parent and are
source-differentiated from the same microscopic law; they are not assigned
post-Feshbach weights and cannot erase a distinct formal order-six
coefficient.

The retarded kernel is

\[
 \chi^R_{EE'}(v,w;t)=-{i\over\hbar}\Theta(t)
 \langle[Q_{E,v}(e,t),Q_{E,w}(e',0)]\rangle.      \tag{FT14}
\]

Equation (FT13) supplies an exact source-off block-local formal operator
moment, and the compatible-state construction supplies nonzero finite matrix
elements.  A uniform sum over all vertices can have additional conservation
or cancellation, so no `k=0` response is inferred from the local result.  A
numerical or state-specific rank for (FT14) is not claimed because FS did not
freeze a stationary density family and `DPAR` is neither inherited nor
adopted.
Operator nonconservation is necessary but not sufficient for a nondegenerate
state-dependent CTP matrix.

## 7. Disposition

The current lawful conclusion is

\[
 \boxed{
 \begin{gathered}
 \text{existing degree/link/ice pair observables}
 \Longrightarrow \text{available dynamical }E\text{ query},\\
 \text{current FQ17a/FS geometric source}
 \Longrightarrow A_1+T_2\text{ rank }4\text{ and }E\text{ null }2,\\
 H[0]\text{ plus those observables}
 \not\Longrightarrow D_EH|_0,\\
 \text{DPAR}
 \Longrightarrow \text{conditional rank }6+\text{nonzero inherited H6 }E
 \text{ commutator}.
 \end{gathered}}                                  \tag{FT15}
\]

Therefore the reduced model has not yet derived the two missing **geometric
strain** directions.  It has localized the entire gap to one physical
question: does the existing degree-pair energy respond to affine sibling
geometry according to a prospectively fixed nonzero law?  Until that premise
is derived or adopted, using the root/cross-dyad weights would be an external
query, not a proof of gravity emergence.
