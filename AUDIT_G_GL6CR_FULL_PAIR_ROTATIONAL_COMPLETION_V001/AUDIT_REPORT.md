# Independent hostile audit — GL6CR full-pair rotational completion

**Target:** `LANE_CROSS_RFT_GRA_GL6CR_FULL_PAIR_ROTATIONAL_COMPLETION_V001`  
**Disposition:** **PASS**  
**Target edits by auditor:** none

## 1. Verdict

The GL6CR algebraic theorem is correct on its declared surface.  An
independent standard-library rational reconstruction, which neither imports
nor executes the author calculation, enumerates the full raw space and finds

\[
 \dim {\cal K}^{(2)}_{S_4}=9,\qquad
 \dim {\cal K}^{(2)}_{SO(3)}=4.
\]

The five stored residual covectors are independent, annihilate the complete
four-dimensional rotational space, and are therefore necessary and
sufficient.  The `T2` condition `B+C=0` is exactly one nonzero projection of
those five.

Most importantly, the independent direct calculation produces a `180 x 9`
polynomial Ward matrix of rank eight.  Its one-dimensional nullspace is

\[
 (0,2,-2,0,-2,-2,4,-1,1),
\]

exactly sixteen times the Reynolds-basis coordinates of the static
Einstein/Fierz--Pauli bilinear.  Thus the shortcut is genuine as an
algebraic theorem: within the complete declared cubic family, the declared
longitudinal null leaves no non-Einstein alternative.

It remains a conditional classifier, not a physical F3 gravity derivation.
The target correctly states that F3 has not yet been shown to generate this
Ward null in a complete same-state 1PI/quotient kernel.

## 2. Independent construction

The replay starts only with the four tetrahedral vectors, the six unordered
port pairs, and the solder definition.  It uses exact fractions throughout.
It:

- constructs every one of the 24 orthogonal tetrahedral actions directly as
  `R_g=(1/4) sum_a |T_g(a)><T_a|`;
- constructs the induced six-pair permutation matrices;
- Reynolds-projects all `21 x 6 = 126` raw symmetric-pair/quadratic-momentum
  seeds, rather than taking the author's nine-element basis as complete;
- separately obtains the dimension from characters;
- reconstructs all four rotational tensor bilinears through the solder;
- derives the five-dimensional annihilator and the `T2` restriction;
- inverts the solder exactly and forms every coefficient of
  `K^(2)(k)D_C^{-1}(k odot xi)`; and
- compares the unique Ward-null result with a separately built
  Einstein/Fierz--Pauli bilinear.

## 3. Tetrahedral action and solder — PASS

With tensor-coordinate order `(xx,yy,zz,xy,xz,yz)`, the audit obtains the
declared solder and the exact inverse

\[
D_C^{-1}=\begin{pmatrix}
0&-1/2&-1/2&0&0&-1\\
-1/2&0&-1/2&0&-1&0\\
-1/2&-1/2&0&-1&0&0\\
-1/2&-1/2&0&1&0&0\\
-1/2&0&-1/2&0&1&0\\
0&-1/2&-1/2&0&0&1
\end{pmatrix}.
\]

Both products with `D_C` are the six-dimensional identity.  For every group
element and every pair basis vector, the replay verifies

\[
 D_CP_g=\rho_{\rm sym}(R_g)D_C.
\]

This checks the direction of the pair permutation and the tensor action,
rather than accepting equivariance from a group average.

The off-diagonal tensor coordinates store the actual entries `h_xy`,
`h_xz`, and `h_yz`, not doubled entries.  Accordingly the Frobenius form in
coordinate space has weights `(1,1,1,2,2,2)`.  The independent construction
uses those weights before pulling a tensor bilinear back with `D_C`; this is
the convention under which all target coefficients agree.

## 4. Complete cubic response space — PASS

For each of the 126 raw seeds, the audit applies the exact group projector

\[
 F(k)\longmapsto\sum_{g\in S_4}P_g^{\mathsf T}F(R_gk)P_g.
\]

Their combined rational rank is nine.  A greedy basis produces the same
nine seed labels stored by GL6CR, and each basis element independently
passes the declared covariance equation

\[
 K(R_gk)=P_gK(k)P_g^{\mathsf T}
\]

for all 24 group elements.

The separate character computation gives the same result.  It also checks
the detailed decomposition: `A1`, `E`, and `T2` each occur three times in
the symmetric bilinears of the pair field, while
`Sym^2(T2)=A1+E+T2`.  Hence the invariant dimension is `3+3+3=9`, not merely
the numerical rank of one proposed parametrization.  The constant cubic
space has dimension three.

## 5. Rotational subspace and five completion conditions — PASS

The audit directly builds the four parity-even self-adjoint rotational
bilinears

\[
 |k|^2h:g,\quad |k|^2\operatorname{tr}h\operatorname{tr}g,\quad
 (hk)\cdot(gk),\quad
 \operatorname{tr}h(k^Tgk)+\operatorname{tr}g(k^Thk).
\]

After exact solder pullback, their rank is four and their coordinates in the
nine-element Reynolds basis agree entry by entry with the target ledger.
All five stored residual covectors annihilate these four rows, have rank
five, and complement them to dimension nine.  This proves both directions
of the completion criterion.

Using the raw pair directions

\[
 e_{01}-e_{23},\qquad e_{02}-e_{13},\qquad e_{03}-e_{12},
\]

the audit independently fits every basis element to

\[
 A|k|^2I+B\,\operatorname{diag}(k_i^2)
 +C\,[kk^T-\operatorname{diag}(k_i^2)].
\]

All 27 stored `A,B,C` coefficients agree exactly.  Their `B+C` row has rank
one.  It is therefore one projection of the full five-condition test, not a
surrogate for full rotational completion.

## 6. Direct Ward-null shortcut — PASS

For each Cartesian basis vector `xi`, the audit constructs the linear pair
polynomial

\[
 x(k,\xi)=D_C^{-1}(k\odot\xi)
\]

and symbolically verifies `D_C x=k odot xi` coefficient by coefficient.
For each of the nine independently constructed cubic basis tensors it then
forms all six output components.  Three choices of `xi` times six outputs
times ten degree-three monomials give the complete `180 x 9` exact
constraint matrix.  No momentum sampling or floating-point rank decision is
used.

The result is

\[
 \operatorname{rank}A_{\rm Ward}=8,\qquad
 \dim\ker A_{\rm Ward}=1.
\]

Separately, the replay builds

\[
 B_E={1\over2}\left(|k|^2h:g-|k|^2\operatorname{tr}h\operatorname{tr}g
 -2(hk)\cdot(gk)+\operatorname{tr}h(k^Tgk)
 +\operatorname{tr}g(k^Thk)\right).
\]

Its Reynolds coordinates are

\[
(0,1/8,-1/8,0,-1/8,-1/8,1/4,-1/16,1/16),
\]

and multiplication by sixteen gives the unique direct Ward ray.  Direct
polynomial multiplication also verifies `K_E x=0` identically.

Within the four-parameter rotational family, the independently derived Ward
rowspace has rank three and is exactly the rowspace of

\[
2a+c=0,\qquad c+2d=0,\qquad b+d=0.
\]

At `k=(2,3,5)`, the pair-coordinate Einstein matrix has rank three, its
three solder-inverted longitudinal vectors are independent, and all three
are killed.  Transforming that matrix to the raw `A1,E1,E2,T1,T2,T3` basis
reproduces the complete GL6CO held-out matrix, providing an independent
cross-lane check of the solder and off-diagonal normalization.

## 7. Scientific typing and ceiling — PASS

The target does not claim that the physical parent already satisfies the
Ward identity.  It distinguishes the currently available connected source
response from the required physical 1PI/quotient kernel and leaves open the
complete same-state `A1+E2+T2` response, contacts, causal continuation,
common cone, record authentication, normalization, gravity, and `G`.

The logical result is precisely

\[
 \left[S_4\text{-covariant quadratic kernel}\right]
 +\left[K^{(2)}(k)D_C^{-1}(k\odot\xi)=0\right]
 \Longrightarrow K^{(2)}\propto K_E^{(2)}.
\]

It does **not** establish the second bracket from F3.  Moreover, a future
physical use must derive the same declared leading longitudinal generator;
inserting it by hand would insert the target condition.  The author packet
states this boundary explicitly and makes no Ricci, gravity, or numerical
`G` promotion.

## 8. Hostile disposition

The audit attacked completeness of the nine-dimensional space, the
Reynolds covariance direction, symmetric off-diagonal factors, solder
inversion, `T2` normalization, residual sufficiency, use of all cubic
monomials, exact rank/nullity, and possible physical overpromotion.  No
material defect remains on the declared algebraic surface.

`PASS__GL6CR_INDEPENDENT_HOSTILE_AUDIT__S4_DIMENSION_9__SO3_DIMENSION_4__FIVE_RELATIONS_NECESSARY_AND_SUFFICIENT__T2_TEST_ONE_PROJECTION__DIRECT_WARD_180_BY_9_RANK_8_NULLITY_1__UNIQUE_EINSTEIN_RAY__PHYSICAL_F3_WARD_1PI_GRAVITY_G_OPEN`
