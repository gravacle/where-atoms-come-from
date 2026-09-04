# GL6CO — INVARIANT CYCLE-RESPONSE AND TENSOR-EXTENSION MATCHING THEOREM

## Status and scope

Assume that a translation-invariant analytic stationary cycle-response
symbol exists for the four elementary-hexagon orientations of the inherited
tetrahedral coordinate parent.  This packet classifies every real,
reciprocal, inversion-even, `S4`-invariant symbol through order `k^2`, pulls
it back through the complete `GL6CL` tensor writer, and derives the exact
condition under which the resulting `T2-T2` block can be extended to an
`SO(3)`-covariant symmetric-rank-two quadratic operator.

The result is a matching/obstruction theorem.  It proves that the required
matching is possible and has codimension one, but symmetry does not enforce
it and the microscopic response coefficients have not yet been calculated.
The `GL6BV` order-`h^2` contact is added only as a separately typed,
same-state conditional block.

This packet does **not** establish the assumed infinite stationary symbol,
select a phase, add contact and spectral pieces from different states,
construct the missing `E2` blocks, perform response-to-1PI inversion, derive
spacetime or a metric, establish a Ricci/Einstein law, prove gravity, or
calculate `G`.

## 1. Four cycle orientations as `A1+T2`

Let the four cycle outputs be indexed by their missing port `d=0,1,2,3`, and
use the inherited tetrahedral coordinate vectors

\[
 T_0=(1,1,1),\quad T_1=(1,-1,-1),\quad
 T_2=(-1,1,-1),\quad T_3=(-1,-1,1).                  \tag{CO01}
\]

Define the orthogonal cycle solder

\[
 u={1\over2}(1,1,1,1)^T,\qquad
 Q_{di}={1\over2}(T_d)_i,\qquad S=(u\;Q).             \tag{CO02}
\]

Then `S^T S=I_4`; `u` carries `A1`, while the three columns of `Q` carry the
tetrahedral three-dimensional irrep `T2`.  The exact replay constructs all
24 port permutations, obtains their signed-permutation coordinate matrices
`R_g`, and verifies

\[
 S^TP_gS=1\oplus R_g.                                  \tag{CO03}
\]

The coordinate `k` in this theorem is the character coordinate of the
inherited embedding.  It is not independently calibrated physical momentum.

## 2. Complete invariant symbol through quadratic order

Put

\[
 r^2=|k|^2,\qquad
 D(k)=\operatorname{diag}(k_x^2,k_y^2,k_z^2),\qquad
 O(k)=kk^T-D(k),                                       \tag{CO04}
\]

and let

\[
 q(k)=(k_yk_z,k_zk_x,k_xk_y)^T.                       \tag{CO05}
\]

In the centered gauge, reciprocity plus inversion makes the static symbol
real symmetric and even.  Before imposing the stationary common-amplitude
null, there are two constant invariant matrices.  `GL6CM` proves that a
uniform change of every cycle amplitude is a stationary spectral null, so
the constant `A1` coefficient vanishes.

The normalization used here is deliberately **bare**.  With `R` the reduced
positive resolvent of the same stationary cycle Hamiltonian, define

\[
 K^{\rm bare}_{cc'}(k)
 =2\operatorname{Re}\langle0|T_cR T_{c'}|0\rangle.    \tag{CO05a}
\]

It contains no pair-source writer coefficient.  Thus the `GL6CM` matrix
(CM10) is `lambda_T^2 K^bare`, whereas the actual `GL6CL` common tensor
source-to-cycle-amplitude derivative used below is `mu B_T`.  This prevents
double counting the writer.  We henceforth abbreviate `K^bare` as
`K_cyc`.  Its coefficients have inverse-energy units and, in the isolated
spectral regime with cycle scale `J`, are `O(J^{-1})`.

The most general remaining symbol is

\[
\boxed{
 S^TK_{\rm cyc}(k)S=
 \begin{pmatrix}
  \alpha r^2 & \eta q(k)^T\\
  \eta q(k) & \kappa I+b r^2I+cD(k)+dO(k)
 \end{pmatrix}+O(r^4).}                                \tag{CO06}
\]

Thus the stationary-null form has one constant coefficient `kappa` and five
independent quadratic coefficients `alpha,eta,b,c,d`.  These are functions
of whichever static or fixed-frequency response coefficient is being
classified; symmetry does not calculate them.

The nonzero `k=0` tensor coefficient `kappa` is a local susceptibility or
background-potential datum.  This packet retains it while classifying the
quadratic jet; it does not prove the zero-derivative/background-stationarity
or masslessness condition required by an Einstein-Hilbert endpoint.

Completeness is exact.  The character replay gives

\[
 \dim\operatorname{Hom}_{S_4}
 [\operatorname{Sym}^2(T_2),
  \operatorname{Sym}^2(A_1\oplus T_2)]=5,             \tag{CO07}
\]

and separately verifies covariance of all five displayed quadratic basis
terms under every one of the 24 group elements.  There is no omitted
quadratic invariant.  If the response is positive semidefinite and
`kappa>0`, then `kappa>=0` and a nonnegative analytic scalar branch requires
`alpha>=0`; complete higher-order Schur constraints are outside this
quadratic classification.

## 3. Pullback through the complete tensor writer

Let `B_T(k)=B_+(k)P_T` be the unscaled complete common tensor writer of the
repaired `GL6CL` packet, and put

\[
 \mu={105\over8}{h^6\over U_d^6}.                     \tag{CO08}
\]

The actual common tensor source-to-cycle-amplitude derivative in the
`GL6CL` common coordinate `j_+` is `mu B_T(k)`.  Recall that `GL6CL` uses

\[
 j_P=j_++j_-,\qquad j_C=j_+-j_-,
\]

so its common coordinate is `j_+=(j_P+j_C)/2`, not the orthonormal
parent/child coordinate.  Therefore the cycle-response contribution to the
pair-source Hessian in that `GL6CL` coordinate is

\[
 \mathcal H^{H6}_T(k)=\mu^2 B_T(k)^*K_{\rm cyc}(k)B_T(k). \tag{CO09}
\]

For comparison with the normalized `GL6BV` common projection, put

\[
 \widehat j_+={j_P+j_C\over\sqrt2}=\sqrt2j_+.
\]

Then `delta a=(mu/sqrt2)B_T widehat j_+`, and the same cycle Hessian in the
orthonormal common coordinate is

\[
 \boxed{\widehat{\mathcal H}^{H6}_T(k)
 ={\mu^2\over2}B_T(k)^*K_{\rm cyc}(k)B_T(k).}          \tag{CO09a}
\]

Use the orthonormal pair basis

\[
 {1\over\sqrt2}(e_{01}-e_{23}),\quad
 {1\over\sqrt2}(e_{02}-e_{13}),\quad
 {1\over\sqrt2}(e_{03}-e_{12}).                       \tag{CO10}
\]

Exact composition through quadratic order gives

\[
\boxed{
 {1\over\mu^2}\mathcal H^{H6}_T(k)
 =8\kappa I+{\cal A}r^2I+{\cal B}D(k)+{\cal C}O(k)+O(r^4),} \tag{CO11}
\]

with

\[
\boxed{
 {\cal A}=-2\kappa+8b,\qquad
 {\cal B}=-16\kappa+8c,\qquad
 {\cal C}=12\kappa+8d.}                               \tag{CO12}
\]

The scalar-cycle coefficients `alpha` and `eta` do not enter (CO11) at this
order.  At `k=0` the tensor writer lands entirely in the cycle `T2` sector;
its first cycle-`A1` component is order `k^2`, so an `A1-A1` or `A1-T2`
quadratic kernel can contribute only at order `k^4` or above.

The replay obtains (CO11)--(CO12) by composing the full four-by-three writer
rows with the four-by-four cycle matrices.  It does not infer them from a
continuum ansatz.

## 4. The exact `SO(3)` tensor-extension condition

If `x=(x_1,x_2,x_3)` are the coordinates in the orthonormal pair basis
(CO10), the inherited solder gives

\[
 (h_{yz},h_{zx},h_{xy})=-{1\over\sqrt2}(x_1,x_2,x_3). \tag{CO13}
\]

The common sign and normalization in (CO13) change a bilinear only by an
overall congruence factor, so they do not change the extension conditions
below.  These coordinates are **not** an ordinary three-vector
representation under `SO(3)`.  A general parity-even, `SO(3)`-covariant
quadratic operator on a symmetric tensor contains, among other full-tensor
terms,

\[
 a r^2 h_{ij}
 +v[k_i k^lh_{lj}+k_j k^lh_{li}].                     \tag{CO14}
\]

Restricting both input and output to (CO13) yields the most general
extendible `T2-T2` form

\[
 u r^2I+v[(r^2I-D)+O].                                 \tag{CO15}
\]

Conversely every matrix (CO15) is the restriction of (CO14), so this is a
necessary-and-sufficient extension test.  In the coefficient convention of
(CO11), it is precisely

\[
 \boxed{{\cal B}+{\cal C}=0.}                          \tag{CO16}
\]

Substitution of (CO12) gives the microscopic cycle-kernel condition

\[
 \boxed{c+d={\kappa\over2}.}                           \tag{CO17}
\]

This answers the matching question exactly:

- it is **not automatic**, because `c` and `d` are independent `S4`
  invariants;
- it is **not impossible**, because (CO17) has positive analytic witnesses;
- it is **one independent linear condition**, not a demand that the `T2`
  block look vector-isotropic.

For example,

\[
 \kappa=1,\quad\alpha=1,\quad\eta=0,\quad
 b=c=d={1\over4}                                      \tag{CO18}
\]

gives

\[
 K_{VV}=I+{1\over4}(r^2I+kk^T),                       \tag{CO19}
\]

which is positive near zero and satisfies (CO17).  Thus positivity does not
obstruct the matching condition.

The constant cycle response alone has `b=c=d=0`, hence
`\({\cal B}+{\cal C}=-4\kappa\)`; it fails (CO16) whenever `kappa` is
nonzero.  The
momentum dependence of the stationary cycle response is therefore essential,
not an optional refinement of the constant response.

## 5. Stronger reference-shape diagnostic is not a Ricci proof

In the same solder convention, the static Euclidean Fierz--Pauli/linearized-
Einstein **reference** `T2-T2` bilinear is, up to overall sign and scale,

\[
 M_{TT}^{\rm ref}(k)=2D-kk^T=D-O.                      \tag{CO20}
\]

Algebraic proportionality of (CO11)'s quadratic term to (CO20) requires
both (CO17) and

\[
 \boxed{b={\kappa\over4},}                             \tag{CO21}
\]

because (CO21) sets `\({\cal A}=0\)`.  The witness (CO18) then gives
`-14(D-O)`.

Equations (CO20)--(CO21) are only a held-out reference-shape diagnostic.  A
stationary **response** Hessian is not automatically the 1PI/inverse kernel;
the full `E2-T2` and `E2-E2` blocks, their relative solder normalization, the
contact, and the response-to-1PI inversion must be supplied before any
Ricci comparison is lawful.  No such promotion is made here.

The replay enforces that ceiling with a full-solder kill test.  In the raw
pair-irrep basis

\[
 A=(1,1,1,1,1,1),\quad
 E_1=(1,-1,0,0,-1,1),\quad
 E_2=(1,1,-2,-2,1,1),                                \tag{CO22}
\]

followed by the three unnormalized `T2` vectors of (CO10), the inherited
solder maps `A` to `-I` and maps the `T2` basis to minus the three symmetric
off-diagonal tensors.  At the generic test momentum `k=(2,3,5)`, direct
application of the static linearized-Einstein reference symbol gives

\[
\begin{pmatrix}
-38&5&37&15&10&6\\
5&100&20&-30&20&0\\
37&20&4&-30&-20&24\\
15&-30&-30&4&-6&-10\\
10&20&-20&-6&9&-15\\
6&0&24&-10&-15&25
\end{pmatrix}.                                       \tag{CO23}
\]

This matrix has rank three, has `T2-T2` block `D-O`, and has nonzero
`A1/E2-T2` cross blocks.  Therefore a later Ricci claim must reproduce the
full cross-block pattern and gauge quotient, after the lawful
response-to-1PI operation; passing only (CO20) is insufficient.

## 6. Separately typed same-state order-`h^2` contact

`GL6BV`, under independent hostile audit, gives the source-before-Feshbach
contact for any translation- and `S4`-invariant locked density.  Let

\[
 p=\Pr[\sigma_{z_v}(a)=\sigma_{z_w}(a)],\qquad
 \zeta=4p-{4\over3},\qquad
 \gamma_c={8\over3}(4p-1).                            \tag{CO24}
\]

With `theta_a=k.T_a`, exact tetrahedral algebra gives

\[
 \sum_a\theta_a^2=4r^2,
 \qquad
 \sum_a\theta_a^2Q_a={4\over3}r^2I+{8\over3}O(k).    \tag{CO25}
\]

The fixed-common-sublattice, inversion-even contact block in the orthonormal
`T2` basis is therefore

\[
\begin{aligned}
 G_{+,T}^{\rm ct}(k)
 ={}&(8+\gamma_c)I\\
 &+{4\over3}(1-4p)r^2I
 +{8\over3}(2p-1)O(k)+O(r^4).                         \tag{CO26}
\end{aligned}
\]

There is no `D` term in (CO26).  Its standalone tensor-extension mismatch is

\[
 \boxed{\Delta_{\rm ct}={8\over3}(2p-1),}             \tag{CO27}
\]

so this block alone passes (CO16) only at `p=1/2`.  The bounded `Q4` orbit
witness of `GL6BV` has `p=109/128` and `Delta_ct=15/8`, but that orbit mixture
has not been proved stationary under the order-six Hamiltonian and is not
substituted for the required same state.

In the connected-functional sign convention, the contact scale is

\[
 g_{\rm ct}={h^2\over4U_d^3}.                          \tag{CO28}
\]

Only if one completed source-first functional proves that the same
stationary state owns both (CO09) and (CO26) may their quadratic coefficients
be added.  In that conditional future calculation, the single tensor-
extension equation would be

\[
\boxed{
 {\mu^2\over2}[-4\kappa+8(c+d)]
 +{h^2\over4U_d^3}{8\over3}(2p-1)=0.}                 \tag{CO29}
\]

The additional reference-shape equation would be

\[
 {\mu^2\over2}[-2\kappa+8b]
 +{h^2\over4U_d^3}{4\over3}(1-4p)=0.                 \tag{CO30}
\]

Equations (CO29)--(CO30) do not assert a cancellation.  They expose exactly
what a same-state calculation would have to determine.  The order-`h^2`
contact and the spectral order-six cycle response may not be fitted against
one another or imported from different states.

There is also a useful strong-lock power-counting consequence.  In the
`GL6CM` spectral regime,

\[
 J={63\over8}{h^6\over U_d^5},\qquad
 \kappa,b,c,d=O(J^{-1}),                              \tag{CO31}
\]

so the cycle term in (CO29)--(CO30) is
`O(h^6/U_d^7)`, while the contact is `O(h^2/U_d^3)`.  Therefore repetition
of the isolated strong-lock cycle block cannot generically cancel a nonzero
leading contact mismatch order by order as `h/U_d -> 0`.  If the full match
is realized, it must involve a finite-ratio or collective stationary regime,
an independently vanishing leading relation, or an additional same-order
block.  This is power counting, not a phase-transition claim.

## 7. Exact disposition

What is proved:

1. the general invariant stationary cycle symbol has five and only five
   quadratic coefficients after its one constant tensor coefficient;
2. only the three cycle-`T2` quadratic coefficients enter the tensor-writer
   pullback through order `k^2`;
3. the pulled-back coefficients are exactly (CO12);
4. `SO(3)` symmetric-tensor extendibility is possible and requires exactly
   the one relation (CO17), not vector isotropy;
5. the stronger held-out reference shape adds (CO21); and
6. the full reference kill test has rank three and indispensable `A1/E2-T2`
   cross blocks, so `T2` matching alone cannot close a Ricci claim; and
7. the order-`h^2` contact has the separately derived mismatch (CO27), with
   the lawful same-state total test given by (CO29).

What remains open:

- existence and calculation of the infinite-volume analytic stationary
  coefficients `kappa,b,c,d` from the actual parent;
- the complete diagonal order-six first-source contribution now being
  treated separately in `GL6CN`;
- a common stationary state for the spectral and contact blocks;
- full `E2-T2/E2-E2` completion and the inherited `h_E/2=h_T` solder test;
- response-to-1PI inversion, time dependence, record authentication,
  refinement, and causal continuation.

Accordingly this theorem narrows the next physics calculation to a single
coefficient relation on the cycle response, while preserving every later
gravity gate as open.
