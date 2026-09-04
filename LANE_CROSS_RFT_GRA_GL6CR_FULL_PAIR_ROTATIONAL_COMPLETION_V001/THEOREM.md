# GL6CR — FULL SIX-PAIR ROTATIONAL-COMPLETION THEOREM

## Status and scope

This packet gives the exact algebraic test that the complete long-wave
six-pair response must pass before it can be interpreted as a rotationally
covariant symmetric-tensor field.  It uses the inherited tetrahedral
four-port action and the already derived pair-memory solder.  It does not
assume a continuum, select a stationary phase, perform the microscopic
response calculation, invert a connected response into a 1PI kernel, or
prove Ricci/Einstein dynamics, gravity, or `G`.

The central result is quantitative.  At order `k^2`, tetrahedral/cubic
symmetry permits nine independent coefficients.  A rotationally covariant
self-adjoint symmetric-tensor operator permits only four.  Therefore the
complete accumulated microscopic kernel must obey five independent
relations.  The `GL6CO` tensor-only condition is exactly one projection of
those five; it cannot replace the other four.

There is also a shorter exact route.  On the full nine-dimensional cubic
family, impose the longitudinal Ward null directly, without first assuming
rotational completion.  The resulting exact constraint matrix has rank
eight, and its one-dimensional nullspace is precisely the
Einstein/Fierz--Pauli ray.  Thus a physically derived Ward identity would
force rotational completion and the Einstein tensor form simultaneously.
The Ward identity is a target condition in this packet, not yet a derived
property of the complete F3 1PI kernel.

## 1. Inherited six-pair field and metric solder

Use pair order

\[
 {cal P}=(01,02,03,12,13,23)                           \tag{CR01}
\]

and tetrahedral vectors

\[
 T_0=(1,1,1),\ T_1=(1,-1,-1),\
 T_2=(-1,1,-1),\ T_3=(-1,-1,1),\qquad v_a=T_a/2.       \tag{CR02}
\]

The inherited `EW/GJ` pair-memory tangent is

\[
 \boxed{D_C(e_{ab})=v_a\odot v_b
 =v_av_b^{\mathsf T}+v_bv_a^{\mathsf T}.}              \tag{CR03}
\]

In coordinate order `(xx,yy,zz,xy,xz,yz)`, its matrix is

\[
 D_C={1\over2}
 \begin{pmatrix}
  1&-1&-1&-1&-1& 1\\
 -1& 1&-1&-1& 1&-1\\
 -1&-1& 1& 1&-1&-1\\
  0& 0&-1& 1& 0& 0\\
  0&-1& 0& 0& 1& 0\\
 -1& 0& 0& 0& 0& 1
 \end{pmatrix},\qquad \operatorname{rank}D_C=6.         \tag{CR04}
\]

For all 24 port permutations, the exact replay verifies

\[
 D_C P_g=\rho_{\rm sym}(R_g)D_C,                       \tag{CR05}
\]

where `P_g` permutes the six pairs and
`rho_sym(R_g)h=R_ghR_g^T`.  Thus (CR03) is not a fitted continuum map: it is
an exact equivariant isomorphism between the inherited pair-memory chart and
`Sym^2(R^3)`.

## 2. Exact dimension of the cubic response space

Let `K^(2)(k)` be a real symmetric six-by-six matrix homogeneous of degree
two in the inherited character coordinate `k`, and require

\[
 K^{(2)}(R_gk)=P_gK^{(2)}(k)P_g^{\mathsf T}.            \tag{CR06}
\]

The exact character count and an independent Reynolds-basis construction
both give

\[
 \boxed{\dim {cal K}^{(2)}_{S_4}=9.}                  \tag{CR07}
\]

Representation-theoretically, the domain is
`Sym^2(T2)=A1+E2+T2`.  Each of those three irreducible sectors occurs three
times in the symmetric bilinears of the six-pair field
`A1+E2+T2`, giving `3+3+3=9` invariant maps.

At `k=0`, cubic symmetry separately permits one scalar coefficient on each
of `A1`, `E2`, and `T2`, hence

\[
 \dim {cal K}^{(0)}_{S_4}=3.                           \tag{CR08}
\]

## 3. Necessary and sufficient rotational form

Write `h=D_Cx` and `g=D_Cy`.  Every parity-even, self-adjoint,
`SO(3)`-covariant bilinear that is quadratic in `k` is uniquely of the form

\[
\boxed{\begin{aligned}
 B_{a b c d}(h,g;k)={}&a|k|^2 h:g
 +b|k|^2\operatorname{tr}h\operatorname{tr}g\\
 &+c(hk)\cdot(gk)\\
 &+d\{\operatorname{tr}h\,(k^{\mathsf T}gk)
       +\operatorname{tr}g\,(k^{\mathsf T}hk)\}.
\end{aligned}}                                        \tag{CR09}
\]

Consequently

\[
 \boxed{\dim {cal K}^{(2)}_{SO(3)}=4,\qquad
        \operatorname{codim}_{S_4}{\cal K}^{(2)}_{SO(3)}=5.}  \tag{CR10}
\]

The replay constructs an exact nine-element Reynolds basis
`I_1,...,I_9`.  If

\[
 K^{(2)}=\sum_{r=1}^9x_rI_r,                            \tag{CR11}
\]

then rotational completion is necessary and sufficient precisely when

\[
 \boxed{\ell_s\cdot x=0,\qquad s=1,\ldots,5,}          \tag{CR12}
\]

for the five exact independent rational covectors stored in
`EXACT_LEDGER.json`.  This is a directly executable test for the future
same-state microscopic kernel; it is not a qualitative isotropy judgment.

At constant order, rotational covariance permits only `h:g` and
`tr(h)tr(g)`, so

\[
 \dim {cal K}^{(0)}_{SO(3)}=2.                         \tag{CR13}
\]

The single missing constant condition removes the cubic splitting between
the `E2` and `T2` parts of one traceless tensor.  It is separate from
background stationarity or masslessness, which require a dynamical
statement rather than a symmetry count.

## 4. Exact relation to the `GL6CO` tensor test

Restricting the complete cubic kernel to the three `T2` pair directions
always gives

\[
 K_{TT}^{(2)}(k)=A|k|^2I+B\,\operatorname{diag}(k_i^2)
 +C\,[kk^{\mathsf T}-\operatorname{diag}(k_i^2)].       \tag{CR14}
\]

Restriction of the four-dimensional rotational family (CR09) has dimension
only two and obeys

\[
 \boxed{B+C=0.}                                        \tag{CR15}
\]

This recovers the content of `GL6CO`, including its
`c_cycle+d_cycle=kappa/2` condition after the H6 writer pullback.  But (CR15)
is only one linear projection of (CR12).  Four further conditions live in
the `A1/E2`, `A1/E2-T2`, and relative-normalization data that a `T2`-only
calculation cannot see.  This proves algebraically why full `E2+T2`
completion is mandatory.

## 5. Direct Ward-null shortcut to the Einstein ray

Within (CR09), requiring invariance under the longitudinal transformation

\[
 h_{ij}\mapsto h_{ij}+k_i\xi_j+k_j\xi_i               \tag{CR16}
\]

gives three independent equations,

\[
 2a+c=0,\qquad c+2d=0,\qquad b+d=0.                   \tag{CR17}
\]

Their unique nonzero ray is

\[
 \boxed{(a,b,c,d)\propto(1,-1,-2,1).}                 \tag{CR18}
\]

Equivalently, with an overall factor `1/2`, the held-out static
linearized-Einstein/Fierz--Pauli bilinear is

\[
\begin{aligned}
 B_{\rm E}(h,g;k)={1\over2}\{&|k|^2h:g
 -|k|^2\operatorname{tr}h\operatorname{tr}g
 -2(hk)\cdot(gk)\\
 &+\operatorname{tr}h(k^{\mathsf T}gk)
 +\operatorname{tr}g(k^{\mathsf T}hk)\}.
\end{aligned}                                         \tag{CR19}
\]

At the exact generic test momentum `(2,3,5)`, its pair-coordinate matrix has
rank three and kills the three longitudinal modes (CR16).  Thus the route
has an exact nested target:

\[
 \boxed{9\ \text{cubic coefficients}
 \xrightarrow{\;5\ \text{relations}\;}4\ \text{rotational coefficients}
 \xrightarrow{\;3\ \text{gauge relations}\;}1\ \text{Einstein ray}.} \tag{CR20}
\]

Equation (CR20) is a target classifier, not a derivation of the target from
the parent.  The microscopic object currently available is a connected
source response.  It must be assembled in one stationary state across all
six channels, and then lawfully converted to the appropriate 1PI/quotient
kernel before (CR17)--(CR19) can be promoted from a held-out comparison to a
gravity result.

The five-plus-three decomposition in (CR20) is useful diagnostically, but it
is not the shortest algebraic route.  Since (CR04) is invertible, define the
unique pair coordinate `x(k,xi)` by

\[
 D_Cx(k,\xi)=k\odot\xi.                                \tag{CR21}
\]

Now impose the Ward null directly on the complete cubic family:

\[
 \boxed{K^{(2)}(k)x(k,\xi)=0
 \quad\hbox{for every }k\hbox{ and }\xi.}              \tag{CR22}
\]

Each component of (CR22) is a homogeneous cubic polynomial in `k`.
Equating all ten cubic-monomial coefficients for all six outputs and three
independent `xi` directions gives an exact rational constraint matrix on
the nine coordinates in (CR11).  The replay proves

\[
 \boxed{\operatorname{rank}A_{\rm Ward}=8,
 \qquad \dim\ker A_{\rm Ward}=1.}                      \tag{CR23}
\]

In the frozen Reynolds basis, its sole ray is

\[
 (0,2,-2,0,-2,-2,4,-1,1),                             \tag{CR24}
\]

which is exactly sixteen times the invariant-basis coordinates of (CR19).
Therefore

\[
 \boxed{S_4\hbox{-covariance}+(\mathrm{CR22})
 \Longrightarrow K^{(2)}\propto K^{(2)}_{\rm E}.}      \tag{CR25}
\]

In particular, rotational completion follows from the direct Ward null; it
does not have to be supplied independently by orientation averaging.  This
removes an algebraic gate, not a physical one: the complete same-state F3
1PI/quotient kernel must still be shown to possess (CR22) from its own
relational redundancy.  Imposing (CR22) by hand would merely insert the
desired answer.

## 6. Disposition

What is proved:

1. the inherited pair-memory solder is a rank-six `S4`-equivariant map to a
   symmetric spatial tensor;
2. the complete cubic quadratic response has exactly nine coefficients;
3. the complete rotational quadratic response has exactly four;
4. five explicit independent relations are necessary and sufficient for
   rotational completion;
5. the `GL6CO` condition is one, and only one, projection of those five; and
6. the gauge-null Einstein reference is the unique one-dimensional ray
   inside the four-dimensional rotational family; and
7. more strongly, the direct Ward null on the full nine-dimensional cubic
   family already has rank eight and forces that same Einstein ray, so
   rotational averaging is not an independent algebraic prerequisite.

What remains physical work:

- calculate the complete same-state `A1+E2+T2` connected response, including
  owner-once contacts and H6 writer contributions;
- derive the longitudinal Ward null for the complete same-state F3
  1PI/quotient kernel from the parent relational redundancy; if successful,
  this enforces all five rotational relations automatically;
- establish the lawful response-to-1PI/quotient step, causal continuation,
  common cone, and record authentication; and
- calculate the surviving coefficient and its physical normalization, from
  which the model value of `G` can be read.

`PASS__FULL_SIX_PAIR_S4_QUADRATIC_DIMENSION_9__SO3_SELF_ADJOINT_DIMENSION_4__FIVE_EXACT_ROTATIONAL_MATCHING_RELATIONS__T2_B_PLUS_C_IS_ONLY_ONE_PROJECTION__DIRECT_CUBIC_WARD_NULL_RANK_8__UNIQUE_EINSTEIN_RAY__SO3_FOLLOWS__MICROSCOPIC_F3_WARD_1PI_CAUSAL_REFINEMENT_GRAVITY_G_OPEN`
