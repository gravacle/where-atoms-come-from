# q4/F3 ice `T2` Fisher-solder boundary theorem

**Lane ID:** `GRA-FN-F3-Q4-ITFSB-V001`

**Short name:** `ITFSB`

**Date:** 2026-08-27

**Claim class:** exact six-state ice-axis realization; exact all-orders
one-link Fisher/QFI covariance; exact lowest nonzero odd-source Taylor
tensor; exact `S4` and `O(3)` classification; exact complement-preserving
no-go; exact generic complement-broken rank theorem conditional on an
independently owned scalar

**Status:**
`FIRST_ODD_DERIVATIVE_ZERO__SECOND_JET_FULL_S4_SYM2__SECOND_JET_NOT_LINEAR_TANGENT__T2_TERM_IS_VECTOR_MEAN_DYAD__COMPLEMENT_PRESERVING_FISHER_TANGENT_E_ONLY__GENERIC_NONZERO_FLUX_PLUS_INDEPENDENT_SCALAR_HAS_RANK6__THAT_RANK6_IS_VECTOR_BACKGROUND_RESPONSE_NOT_TENSOR_MODE__PHYSICAL_METRIC_SOLDER_REMAINS_OPEN`

**Not claimed:** a physical preparation of the saturated local family; a
same-parent owner for the scalar `A1`; an autonomous nonsymmetric phase; an
`O(3)`-covariant metric response; an isolated tensor pole; helicity two; a
massless graviton; RGRL-B; universal stress coupling; Einstein dynamics;
gravity; or numerical `G`

## 1. Exact question and frozen input

`FK` proves that the diagonal functions on one degree-four `d_*=2` ice fiber
decompose as

\[
 \mathbb R^{\Omega_2}=A_1\oplus E\oplus T_2,
 \qquad \operatorname{Sym}^2(V)=A_1\oplus E\oplus T_2,       \tag{FN01}
\]

and that the inherited one-link statistics `s_a` carry `T2`, while the
normalized pair statistics carry `E`.  It also proves that the most direct
one-link Fisher query has zero first derivative in every odd `T2` source at
the complement-symmetric point.  The independent audit accepts that result
and leaves the physical `T2` metric solder open.

This packet answers the next finite question without adding an interaction:

1. What is the exact lowest nonzero odd-source response?
2. Does that response give a six-dimensional **linear** metric tangent?
3. Can any complement-preserving state or another equivariant local vector
   query recover the missing `T2` tangent?
4. What is the minimal background on which a full linear rank can occur,
   and what physical mode does that construction actually probe?

The load-bearing frozen inputs are:

| role | dependency | SHA-256 |
|---|---|---|
| exact finite input | `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md` | `cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98` |
| exact finite input | `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md` | `c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4` |
| conditional spectral interpretation only | `LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md` | `98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452` |
| conditional spectral interpretation only | `LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/INDEPENDENT_AUDIT.md` | `327bf6a4476c4c6382757dc156a96c6032233d34c25c1f7935e2582acf6c607a` |

All statements through Section 6 are finite exact consequences of the six
ice states and use only the `FK` inputs.  The pole interpretation in Section
7 invokes the separately audited `FL` classification and remains conditional
on its named infrared premise; it is not an additional finite theorem.

## 2. The six ice states are three antipodal orthogonal axes

Let

\[
 V=\mathbf1^\perp\subset\mathbb R^4,
 \qquad
 \Omega_2=\{s\in\{-1,+1\}^4:\mathbf1^{\mathsf T}s=0\}.     \tag{FN02}
\]

Define the orthonormal basis of `V`

\[
 \begin{aligned}
 u_1&={1\over2}(1,1,-1,-1),\\
 u_2&={1\over2}(1,-1,1,-1),\\
 u_3&={1\over2}(1,-1,-1,1).
 \end{aligned}                                               \tag{FN03}
\]

Then the local ice fiber is exactly

\[
 \boxed{\Omega_2=\{+2u_1,-2u_1,+2u_2,-2u_2,+2u_3,-2u_3\}.} \tag{FN04}
\]

Thus a local ice state occupies one of three orthogonal axes with one of two
complement signs.  In the `u_i` basis write its coordinates as
`z_i=u_i dot s`.  Every state has exactly one nonzero coordinate, equal to
`+2` or `-2`.

This realization is useful because it separates two otherwise easy-to-blur
structures:

- the three **axes** are the even matching/`E` information; and
- the sign on an occupied axis is the odd one-link/`T2` information.

The full `S4` action becomes a signed permutation action on the three axes.
In this basis the restriction of a continuum symmetric tensor is

\[
 \operatorname{Sym}^2(V)
 =\underbrace{\mathbb RI}_{A_1}
 \oplus\underbrace{\{\text{diagonal trace-free}\}}_E
 \oplus\underbrace{\{\text{off-diagonal symmetric}\}}_{T_2}. \tag{FN05}
\]

## 3. Exact all-orders Fisher/QFI covariance

Turn on only an odd one-link source `t in V`, expressed in the basis (FN03).
The inherited exponential query is

\[
 p_t(z)={e^{t\cdot z}\over Z(t)},
 \qquad z\in\{\pm2e_1,\pm2e_2,\pm2e_3\}.                    \tag{FN06}
\]

Put

\[
 C(t)=\sum_{i=1}^3\cosh(2t_i),\qquad
 w_i(t)={\cosh(2t_i)\over C(t)},\qquad
 m_i(t)={2\sinh(2t_i)\over C(t)}.                          \tag{FN07}
\]

The partition function is `Z=2C`, the mean is `m=E_t[z]`, and the second
moment is `4 diag(w)`.  Therefore the exact covariance is

\[
 \boxed{
 F(t)=\operatorname{Cov}_{p_t}(z)
 =4\operatorname{diag}(w(t))-m(t)m(t)^{\mathsf T}.}          \tag{FN08}
\]

Because (FN06) is a positive commuting exponential family, this covariance
is both its classical Fisher information matrix and its commuting-state SLD
quantum Fisher information matrix for the one-link source parameters.  No
noncommuting QFI identification is being made.

Equation (FN08) is all-orders and already identifies the only source of its
off-diagonal `T2` sector:

\[
 \boxed{F_{ij}(t)=-m_i(t)m_j(t)\quad(i\ne j).}               \tag{FN09}
\]

The even second moment is diagonal.  Every `T2` entry comes from subtracting
the dyad of a nonzero odd vector mean.

## 4. Lowest nonzero odd-source tensor

Global complement sends `t -> -t` and `z -> -z`, so `F(-t)=F(t)`.  All odd
Taylor coefficients vanish.  Expanding (FN08) at the uniform point gives

\[
 \boxed{
 F(t)={4\over3}I+Q(t)+O(\|t\|^4),}                          \tag{FN10}
\]

where

\[
 \boxed{
 Q(t)={8\over3}\operatorname{diag}(t_1^2,t_2^2,t_3^2)
       -{8\over9}\|t\|^2I
       -{16\over9}tt^{\mathsf T}.}                          \tag{FN11}
\]

Equivalently, for two source directions `x,y in V`,

\[
 \boxed{
 \begin{aligned}
 D^2F_0[x,y]
  ={}&{16\over3}\operatorname{diag}(x_1y_1,x_2y_2,x_3y_3)
       -{16\over9}(x\cdot y)I\\
     &-{16\over9}(xy^{\mathsf T}+yx^{\mathsf T}).
 \end{aligned}}                                             \tag{FN12}
\]

This is the lowest nonzero coefficient for every nonzero one-direction
source because

\[
 \operatorname{tr}Q(t)=-{16\over9}\|t\|^2<0\quad(t\ne0).   \tag{FN13}
\]

### Theorem `ITFSB-1` -- exact second-jet classification

For `S=(xy^T+yx^T)/2`, define the polarized Hessian operator

\[
 \mathcal H(S)
 ={16\over3}\operatorname{diag}(S)
  -{16\over9}\operatorname{tr}(S)I
  -{32\over9}S.                                             \tag{FN14}
\]

On the `S4` decomposition (FN05), it acts by the exact nonzero scalars

\[
 \boxed{
 \mathcal H\big|_{A_1}=-{32\over9},\qquad
 \mathcal H\big|_E={16\over9},\qquad
 \mathcal H\big|_{T_2}=-{32\over9}.}                       \tag{FN15}
\]

Hence the **polarized second jet** is an `S4`-equivariant isomorphism

\[
 \operatorname{Sym}^2(T_2)\longrightarrow
 \operatorname{Sym}^2(V).                                  \tag{FN16}
\]

This is the strongest positive result in the packet.  It is not a
six-dimensional linear state tangent.  One physical source vector `t` has
only three components, and `t tensor t` lies on the rank-one quadratic cone.
The six-dimensional span in (FN16) is obtained by polarizing two independent
interventions, or equivalently by taking linear combinations of second
derivatives.  It is a nonlinear susceptibility.

Under `O(3)`, `Sym^2(V)=ell=0 plus ell=2`, while under the tetrahedral `S4`
subgroup `ell=2` restricts to `E plus T2`.  An `O(3)`-equivariant Hessian
would act with one scalar on all five `ell=2` directions.  Equation (FN15)
has different eigenvalues on `E` and `T2`.  Therefore

\[
 \boxed{\mathcal H\text{ is }S4\text{-equivariant but not }O(3)
 \text{-equivariant}.}                                      \tag{FN17}
\]

The isotropic value `F(0)=4I/3` does not by itself provide rotationally
covariant nonlinear metric dynamics; the second response exposes the
underlying tetrahedral axes.  In plain language, the Hessian is not `O(3)`-equivariant.

## 5. Exact complement-preserving no-go for every local state

The preceding result is not restricted to the exponential curve.  Write an
arbitrary interior probability distribution on the six states as

\[
 p_i^+=p(+2u_i),\qquad p_i^-=p(-2u_i),qquad
 w_i=p_i^++p_i^-,\qquad
 m_i=2(p_i^+-p_i^-),                                        \tag{FN18}
\]

with `sum_i w_i=1`.  Direct summation gives the exact universal formula

\[
 \boxed{F(p)=4\operatorname{diag}(w)-mm^{\mathsf T}.}        \tag{FN19}
\]

Here `F(p)` means the covariance of the fixed one-link statistic `z` at
background state `p`.  Equivalently, it is the Fisher/commuting-SLD-QFI
matrix at `theta=0` of the local tilted query

\[
 p^{(p)}_\theta(z)
 ={p(z)e^{\theta\cdot z}\over
   \sum_{z'}p(z')e^{\theta\cdot z'}}.                       \tag{FN19a}
\]

It is not being identified with the intrinsic five-dimensional Fisher
metric of the whole probability simplex.

Complement symmetry is exactly `p_i^+=p_i^-`, or `m=0`.  Hence every
complement-symmetric local state satisfies

\[
 F(p)=4\operatorname{diag}(w),                               \tag{FN20}
\]

and every complement-preserving normalized tangent satisfies

\[
 dF=4\operatorname{diag}(dw),\qquad \sum_i dw_i=0.          \tag{FN21}
\]

Its image is precisely the two-dimensional diagonal trace-free `E` sector.
Adding an independently controlled common scalar `rho I` raises this only to
`A1 plus E`, of rank three.  It still supplies no `T2` tangent.

### Theorem `ITFSB-2` -- sharp symmetric-family boundary

No open complement-preserving probability family on the exact six-state ice
fiber gives a first-order `T2` tangent in the one-link Fisher/QFI covariance.
This holds at symmetric and nonsymmetric axis weights `w`; it is not merely a
zero derivative at the uniform point.

There is also no unused equivariant local vector query hidden in the same
diagonal observable algebra.  The six-state function module contains `T2`
with multiplicity one:

\[
 \mathbb R^{\Omega_2}=A_1\oplus E\oplus T_2.                \tag{FN22}
\]

Consequently every `S4`-equivariant `V`-valued diagonal statistic is a scalar
multiple of the inherited one-link vector `z`, up to the choice of basis in
`V`.  Its covariance obeys the same boundary.  Mapping the inequivalent
pair `E` directly into a vector, or selecting a different set of axes,
requires an extra frame or symmetry-breaking rule not owned by this fiber.

This does not exclude every imaginable nonlocal or noncommuting query.  It
is the sharp no-go for the exact local diagonal state/query domain inherited
from `FK`.

## 6. The minimal full-rank background and its cost

Differentiate (FN19):

\[
 dF=4\operatorname{diag}(dw)
     -(m\,dm^{\mathsf T}+dm\,m^{\mathsf T}),
 \qquad \sum_i dw_i=0.                                      \tag{FN23}
\]

The off-diagonal projection of the `dm` map, in the component order
`(12),(13),(23)`, has coefficient matrix

\[
 -\begin{pmatrix}
 m_2&m_1&0\\
 m_3&0&m_1\\
 0&m_3&m_2
 \end{pmatrix},
 \qquad
 \det=2m_1m_2m_3.                                          \tag{FN24}
\]

Therefore all three off-diagonal `T2` tangent directions occur exactly when

\[
 m_1m_2m_3\ne0.                                             \tag{FN25}
\]

The five-parameter normalized probability family then has metric-map rank
five: three odd/vector directions plus two even/axis-weight directions.  If
one additionally owns an independent common scalar and defines the augmented
candidate

\[
 G(p,\rho)=F(p)+\rho I,                                     \tag{FN26}
\]

then the two `dw` directions plus `d rho` span all diagonal matrices, while
(FN24) spans all off-diagonal matrices.  Thus:

### Theorem `ITFSB-3` -- exact generic broken-background rank

On the open interior set `m_1m_2m_3 != 0`, the augmented map (FN26) has full
linear rank six.  Without the scalar it has rank five.  At `m=0` its ranks
are three with the scalar and two without it.  If any component of `m`
vanishes, even the augmented map has rank below six.

These are ranks of the map from five background-probability controls, plus
the optional sixth scalar control, into the six components of the
three-by-three one-link query covariance.  They are not the rank of a
six-parameter Fisher metric and do not supply six propagating degrees of
freedom.

The saturated five-parameter exponential family already written in `FK`
mathematically reaches every interior distribution, so such open points
exist without adding a local operator.  But their physical costs are
load-bearing:

1. `m != 0` is a complement-breaking one-link/vector background;
2. full rank requires nonzero components on all three orthogonal ice axes;
3. preparation and stabilization of that background by the F3 parent have
   not been derived;
4. the `A1` scalar remains independently supplied and has not been soldered
   to the same parent, clock, or calibration; and
5. the `T2` response is exactly the linearization of `-m m^T`:

\[
 \delta F_{T_2}
 =-\operatorname{offdiag}(m\,\delta m^{\mathsf T}
                          +\delta m\,m^{\mathsf T}).         \tag{FN27}
\]

It is therefore a vector-background response, not an independent symmetric-
tensor degree of freedom.

There is no nonzero threshold hidden here.  An arbitrarily small generic
`m` gives algebraic rank six after the scalar is added, but the determinant
of the off-diagonal block scales as `m_1m_2m_3`; the construction becomes
singular as the complement-symmetric point is approached.

## 7. Symmetrization does not evade the boundary

Let `C p` be the complemented distribution, so its mean is `-m` and its
covariance is the same `F(p)`.  The unlabelled symmetric mixture

\[
 \bar p={1\over2}(p+Cp)                                     \tag{FN28}
\]

has zero mean and obeys the exact total-covariance identity

\[
 \boxed{
 F(\bar p)=F(p)+mm^{\mathsf T}=4\operatorname{diag}(w).}     \tag{FN29}
\]

Thus randomized `+m/-m` control has two different meanings:

- averaging the **conditional** covariances while retaining the control-sign
  label keeps `F(p)` and its `T2` dyad, but the answer is conditioned on an
  external classical record;
- discarding that label to form one genuinely complement-symmetric physical
  state adds the between-branch covariance `mm^T` and cancels the `T2` dyad.

No complement-preserving randomized-query loophole exists within the owned
single-copy Fisher covariance.

The dynamical interpretation is now precise.  At the symmetric point, the
full `S4` second jet (FN15) is a two-vector nonlinear susceptibility.  Under
the audited `FL` premise `MAXWELL-IR`, that response is in the even two-vector
channel and does not by itself create a one-particle tensor pole.  At a
broken background, (FN27) makes the response linear in the same vector
fluctuation.  It therefore inherits the spin-one photon pole already carried
by that vector under `MAXWELL-IR`, not an independently derived helicity-two
pole.  This pole statement is conditional on `MAXWELL-IR` and is not needed
for the finite rank/no-go proofs above.

## 8. Disposition and promotion ceiling

The exact result is

\[
 \boxed{
 \begin{gathered}
 D_tF(0)=0,\\
 D_t^2F(0):\operatorname{Sym}^2(T_2)
      \overset{S4}{\cong}A_1\oplus E\oplus T_2,\\
 \text{but the }T_2\text{ term is }-mm^T\text{ and is not a
 linear tensor mode},\\
 m=0\Longrightarrow\text{first-order Fisher image }=E,\\
 m_1m_2m_3\ne0+\text{independent }A_1
 \Longrightarrow\text{rank six vector-background response}.
 \end{gathered}}                                             \tag{FN30}
\]

This closes the lowest-order local Fisher calculation and the proposed
complement-breaking workaround.  It does **not** close physical metric
solder.  Promotion beyond this packet requires at least one of:

1. an inherited even rank-two collective variable with its own linear
   response and tensor Ward identity;
2. a same-parent non-Gaussian kernel producing an isolated tensor branch
   rather than a two-vector continuum; or
3. a derived nonlocal/noncommuting information query whose `T2` metric
   tangent survives complement symmetry and whose continuum transformation
   law is proved.

No new fitted interaction, unowned scalar, programmed vector background, or
conditional control label may be renamed as that closure.
