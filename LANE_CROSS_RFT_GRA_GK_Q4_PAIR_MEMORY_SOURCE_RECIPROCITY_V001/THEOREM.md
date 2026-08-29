# q4 pair-memory expectation / physical-source reciprocity theorem

**Lane ID:** `CROSS-RFT-GRA-GK-Q4-PAIR-MEMORY-SOURCE-RECIPROCITY-V001`

**Official short name:** `PMSR`

**Date:** 2026-08-29

**Status:** `MUTABLE_BUILDER_VERIFIED__INDEPENDENT_HOSTILE_AUDIT_PENDING`

**Claim class:** exact finite commuting same-parent theorem; exact separation
of natural controls from observable pair-memory expectation coordinates;
exact q4 expectation-to-Fisher identity; exact conditional composition with
the FU/FT physical strain source; exact mixed-derivative reciprocity; exact
finite-temperature rank and pure-ice rank boundary

**Not claimed:** that a control knob is a record; that an ensemble moment is
a qualified record without an authenticated retained sufficient statistic;
that the commuting Gibbs carrier is the complete noncommuting F3 parent;
that FY carries the EW localization query or fidelity metric through its
Feshbach reduction; physical-space soldering, gluing, a causal cone, a Ward
identity, a massless pole, Einstein dynamics, gravity, or `G`

## 1. Exact question and dependency boundary

EW proved that six q4 pair natural parameters can deform all six components
of one rank-three localization Fisher tensor.  FT and FU separately proved,
conditionally on a complete physical solder, that the same six pair
operators can be conjugate to an affine physical strain.  FV and FY then
carried the resulting direct pair source into an ice projection and a finite
H6 response calculation.

Those results do not by themselves prove a common physical mechanism.  In
particular, the EW natural parameters `J_ab` are controls of an exponential
family.  A retained record is not a value of a control knob.  This lane asks
for the smallest exact common-parent statement that uses observable memory
coordinates instead, and asks exactly how far that statement survives toward
the projected FY calculation.

The dependency hashes are frozen in `DEPENDENCIES.sha256`.  The lawful answer
is:

1. observable pair expectations deform the EW Fisher tensor by an exact
   affine isomorphism;
2. in one commuting Gibbs realization, the FU pair-sector physical source and that
   Fisher deformation are mixed derivatives of the same generating
   functional;
3. the resulting tangent has rank six at every finite positive temperature
   when the **complete** pair-sector tangent has the full DPAR form with a
   nonzero common slope; but
4. exact projection to the ice fiber collapses the pair-memory susceptibility
   to rank two, and FY has not yet transported the EW query/state family.

Thus the theorem is an exact microscopic bridge, not a gravity endpoint.

## 2. q4 notation and the control/record distinction

Let

\[
 \mathbf 1=(1,1,1,1)^{\mathsf T},\qquad
 P=I_4-\frac14\mathbf1\mathbf1^{\mathsf T},\qquad
 V=\mathbf1^\perp,\qquad v_a=Pe_a.                 \tag{GK01}
\]

For `s in {-1,+1}^4`, define

\[
 X(s)=\sum_av_as_a=Ps,\qquad
 Y_{ab}(s)=s_as_b,\qquad a<b.                      \tag{GK02}
\]

Write `e={a,b}` for one of the six unordered edges and

\[
 B_e=v_a\odot v_b:=v_av_b^{\mathsf T}+v_bv_a^{\mathsf T}.
                                                               \tag{GK03}
\]

The EW family is

\[
 p_{\theta,J}(s)=\exp\!\left[\theta\!\cdot\!X(s)
       +\sum_eJ_eY_e(s)-W(\theta,J)\right].        \tag{GK04}
\]

There are two categorically different six-vectors:

\[
 \boxed{J_e=\hbox{natural control},\qquad
 C_e:=\langle Y_e\rangle=\partial_{J_e}W
       =\hbox{observable pair expectation}.}       \tag{GK05}
\]

Only `C`, or an authenticated retained sufficient statistic that realizes
it, is eligible to be a collective memory coordinate.  `J` is not renamed as
memory in this lane.

At `theta=0`, global spin flip gives `<X>=0`.  The pair covariance is

\[
 G_{ef}:=\partial_{J_f}C_e
        =\operatorname{Cov}(Y_e,Y_f).              \tag{GK06}
\]

The six Walsh pair characters are linearly independent modulo constants.
Every finite `J` gives full support, so `G` is positive definite: for a
nonzero vector `z`, `z^T G z` is the variance of the nonconstant function
`sum_e z_eY_e`.  Consequently `J -> C` is locally invertible throughout the
finite family.  This local invertibility licenses a change of coordinates;
it does not turn `J` into a record.

## 3. Theorem `PMSR-1` -- exact expectation-memory/Fisher identity

The following identity holds state by state:

\[
 \boxed{
 X(s)X(s)^{\mathsf T}
 =P+\sum_{a<b}B_{ab}Y_{ab}(s).}                   \tag{GK07}
\]

Indeed, expanding `XX^T` separates the four diagonal terms from the six
unordered cross terms, and
`sum_a v_av_a^T=P`.  At `theta=0`, `<X>=0`, so the localization Fisher tensor
is therefore

\[
 \boxed{
 {\cal F}_\theta(C)
 =\operatorname{Cov}(X)
 =I_V+\sum_e C_eB_e.}                             \tag{GK08}
\]

Here `P|_V=I_V`.  Define

\[
 {\cal B}:\mathbb R^6\longrightarrow\operatorname{Sym}^2(V),
 \qquad {\cal B}(c)=\sum_ec_eB_e.                 \tag{GK09}
\]

EW proved that the six `B_e` form a basis of
`Sym^2(V)`.  Hence

\[
 \boxed{D_C{\cal F}_\theta={\cal B}\quad
        \hbox{is an exact isomorphism}.}           \tag{GK10}
\]

This is stronger and cleaner than calling `J` pair memory.  The metric is an
exact affine function of observable pair expectations, while the response of
those expectations to a control is separately governed by `G`:

\[
 D_J{\cal F}_\theta={\cal B}G.                    \tag{GK11}
\]

For finite full support, both factors are invertible.

Equation (GK08) is a collective-state statement.  One individual four-bit
outcome supplies only a discrete constrained vector `Y(s)`.  An open
six-dimensional expectation neighborhood requires a joint query law or
retained sufficient counts across an actual qualified collection.

## 4. Conditional physical same-parent solder

Introduce a real symmetric physical strain `j` with six coordinates `j_A`.
Let

\[
 \widehat R_{ab}=
 { (n_b-n_a)(n_b-n_a)^{\mathsf T}\over|n_b-n_a|^2},
 \qquad M_{eA}=\widehat R_e:E_A,                  \tag{GK12}
\]

where the `n_a` are the regular tetrahedral directions and the `E_A` form a
basis of real symmetric three-tensors.  The six normalized root dyads form a
basis, so the strain-to-edge map `M` is invertible.

This lane uses the following explicit, load-bearing same-parent packet.

- **`GK-S1 -- EW carrier`:** the scored local carrier is the full-support
  commuting family (GK04), with the complete diagonal query used by EW.
- **`GK-S2 -- physical Gibbs solder`:** the same pair coefficients that
  define the physical strain source define the pair natural controls,

  \[
   p_{\theta,j}(s)\propto
   \exp[\theta\cdot X(s)-\beta H_C(j;s)],
   \qquad 0<\beta<\infty.                         \tag{GK13}
  \]

- **`GK-S3 -- FU completion`:** FU's complete physical conditions `S1`--`S9`
  hold, including the grounded field, conserved dressed transfer, source
  ordering, source-off match, and noncircular physical length/charge/kernel.
- **`GK-S4 -- complete DPAR form and noncancellation`:** after every
  pair-sector remainder is included, the complete nonidentity tangent is

  \[
   \partial_{j_A}K_e\big|_0
   =-{U_d\lambda_{\rm pair}^{\rm net}\over2}M_{eA},
   \qquad \lambda_{\rm pair}^{\rm net}\ne0.       \tag{GK-S4}
  \]

  Thus the remainders vanish, are identity-only, or share the same common
  normalized radial law.  FU's weaker conclusion
  `lambda_E^net != 0` proves a noncancelled `E` component, but it does not by
  itself prove that the whole six-direction pair tangent equals one scalar
  multiple of `M`.  With a general complete tangent `L` in place of
  `lambda_pair^net M`, equations (GK17)--(GK19) must use `L`, and rank six
  additionally requires `det L != 0`.
- **`GK-S5 -- fixed localization statistic`:** `X`, its material-coordinate
  coframe, and the complete query have no omitted explicit `j` dependence in
  the scored derivative.  If they do, their additional derivative terms must
  be retained and (GK17) below is not the complete source response.

In the commuting scored sector write

\[
 H_C(j;s)=h_{\rm id}(j)+\sum_eK_e(j)Y_e(s),        \tag{GK14}
\]

with the complete FU/DPAR tangent

\[
 K_e(0)={U_d\over2},\qquad
 \partial_{j_A}K_e\big|_0
 =-{U_d\lambda_{\rm pair}^{\rm net}\over2}M_{eA}.
                                                               \tag{GK15}
\]

Since `J_e(j)=-beta K_e(j)`, the sign and factor are fixed:

\[
 \boxed{
 \partial_{j_A}J_e\big|_0
 ={\beta U_d\lambda_{\rm pair}^{\rm net}\over2}M_{eA}.}
                                                               \tag{GK16}
\]

Combining (GK10), (GK06), and (GK16) gives

\[
 \boxed{
 D_j{\cal F}_\theta
 ={\beta U_d\lambda_{\rm pair}^{\rm net}\over2}\,{\cal B}GM.}
                                                               \tag{GK17}
\]

More generally, if the audited complete pair tangent is
`partial_j J=(beta U_d/2)L`, with no common-DPAR reduction, the exact formula
is

\[
 D_j{\cal F}_\theta={\beta U_d\over2}{\cal B}GL.   \tag{GK17g}
\]

At finite full support, `B` and `G` are invertible, so this generalized
response has rank six if and only if `L` has rank six.

### Theorem `PMSR-2` -- finite-temperature rank-six composition

Under `GK-S1`--`GK-S5`,
`U_d lambda_pair^net != 0`, and finite positive `beta`, all three maps in
(GK17) are invertible.  Therefore

\[
 \boxed{\operatorname{rank}D_j{\cal F}_\theta=6.} \tag{GK18}
\]

The exact pair-sector physical conjugate is

\[
 \boxed{
 Q_A:=-2\partial_{j_A}H_C\big|_0
 =U_d\lambda_{\rm pair}^{\rm net}\sum_eM_{eA}Y_e
  +q_{{\rm id},A}I.}                             \tag{GK19}
\]

Equations (GK15), (GK16), and (GK19) fix the factor of two and the sign.  An
identity source is retained in the work ledger but cannot change a normalized
state or the localization Fisher tensor.

This theorem has a deliberately stronger whole-pair-source premise than
FU's `E`-only noncancellation result.  Merely proving a nonzero `E`
coefficient is insufficient for the six-direction formula (GK17).  The
source-off operator `H_C(0)` alone does not determine its strain derivative;
FT's same-source-off nonuniqueness theorem still applies if `GK-S3` or
`GK-S4` is absent.

## 5. Theorem `PMSR-3` -- generating-functional reciprocity

Let

\[
 W(\theta,j)=\log\sum_s
 \exp[\theta\cdot X(s)-\beta H_C(j;s)].            \tag{GK20}
\]

Analyticity is automatic for this finite full-support family.  From
`Q_A=-2 partial_{j_A}H_C`,

\[
 \partial_{j_A}W={\beta\over2}\langle Q_A\rangle.\tag{GK21}
\]

Since
`F_{theta,mn}=partial_{theta_m}partial_{theta_n}W`, commuting mixed
derivatives gives the exact reciprocity law

\[
 \boxed{
 \partial_{j_A}{\cal F}_{\theta,mn}
 ={\beta\over2}
 \partial_{\theta_m}\partial_{\theta_n}
 \langle Q_A\rangle.}                            \tag{GK22}
\]

For the EW squared complete-query fidelity,

\[
 \gamma(\theta,\theta+d\theta;j)
 =1-{1\over4}d\theta^{\mathsf T}{\cal F}_\theta(j)d\theta
   +O(|d\theta|^3),                               \tag{GK23}
\]

and hence

\[
 \boxed{
 \partial_{j_A}[-\log\gamma]
 ={\beta\over8}d\theta^{\mathsf T}
 \left[
  \partial_\theta\partial_\theta\langle Q_A\rangle
 \right]d\theta+O(|d\theta|^3).}                \tag{GK24}
\]

This is the exact common-generating-functional bridge.  It says that the
pair-sector physical-source response and the change in record
distinguishability share one parent derivative.  It does not say that `gamma`
is a force, energy, or curvature, and it does not supply the omitted complete
F3 source sectors.

## 6. Exact uniform thermal spectrum and singular limits

At source off, put

\[
 u=-{\beta U_d\over2},\qquad
 a=e^{2u}=e^{-\beta U_d}>0,\qquad
 D=a^4+4a+3.                                      \tag{GK25}
\]

Multiplying all Boltzmann weights by the same factor, the magnetization
sectors `|sum s_a|=4,2,0` have respective weights `a^4,a,1`.  For every edge
and for the four-spin character,

\[
 c:=\langle Y_e\rangle={a^4-1\over D},\qquad
 q:=\langle s_1s_2s_3s_4\rangle
 ={a^4-4a+3\over D}.                              \tag{GK26}
\]

The pair covariance has diagonal, adjacent-edge, and opposite-edge entries

\[
 1-c^2,\qquad c-c^2,\qquad q-c^2.                 \tag{GK27}
\]

On the edge representation
`A1 direct-sum E direct-sum T2`, its eigenvalues are

\[
 \boxed{
 \begin{aligned}
 g_{A_1}&={8a(3a^4+4a^3+1)\over D^2},&&\text{multiplicity }1,\\
 g_E&={8\over D},&&\text{multiplicity }2,\\
 g_{T_2}&={8a\over D},&&\text{multiplicity }3.
 \end{aligned}}                                   \tag{GK28}
\]

All are strictly positive for finite `a>0`, independently confirming the
finite-temperature rank claim.

Three exact controls show why operator availability is not enough.

1. **Infinite-temperature control.**  At `beta=0`, `a=1` and `G=I_6`, but
   `D_jJ=(beta U_d lambda_pair^net/2)M=0`.  A nonzero source operator need not deform
   the state or its Fisher metric.
2. **Zero-slope control.**  If the complete `lambda_pair^net` vanishes, both (GK19)'s
   nonidentity pair source and (GK17) vanish.
3. **Pure-ice control.**  In the exact `a -> 0` limit,

   \[
    g_{A_1}\to0,\qquad g_E\to{8\over3},\qquad
    g_{T_2}\to0,                                  \tag{GK29}
   \]

   so `G` has rank two.  This matches FV's result that the direct pair source
   has only a nonidentity `E` image after ice restriction.  The unchanged
   edge source can supply `A1+T2` operator directions in the larger FT/FV
   construction, but this lane has not identified those directions with the
   same EW expectation-memory metric.

## 7. Exact FY interface and next physical calculation

FV and FY use the same `P_ab=Z_aZ_b`, the same normalized root dyads, and the
same source convention `Q=-2 partial_j H`.  FY therefore retains exact
**operator lineage** from the direct FU pair source into its complete H6
finite-graph source.

That is the full lawful FY conclusion here.  FY's selected 180-state ice
component is a noncommuting, Feshbach-reduced finite quantum system.  It did
not insert the EW localization source `theta`, construct the complete EW
query on that parent, or prove that its ground-state/CTP information metric
equals (GK08).  For a general noncommuting Gibbs parent, the Hessian of
`log Z` is the Bogoliubov--Kubo--Mori metric, not automatically the SLD/QFI
metric obtained from squared state fidelity.  Therefore FY is not promoted
to a gamma theorem.

The next exact calculation is now sharply typed: insert both `theta` and `j`
into the same FU09b-dressed, port-complete F3 parent before Feshbach reduction;
retain the temporal, current, work, boundary, and contact terms; carry both
derivatives through the reduction; and then compare the resulting physical
fidelity or CTP metric with the EW complete-query metric.  The product
`BGM` must be re-earned there rather than imported from the commuting model.

## 8. Record binding and final disposition

To call `C` retained physical memory, supply EW's authenticated same-parent
record bind: formation margin, retention, complete-query census,
reference-stable descent, lineage custody, and actual joint law.  For `N`
qualified copies, the retained sufficient statistics `sum_r Y_e(s^(r))`
make the empirical pair moments operational.  Pair preparation, controller,
reservoir, heat, work, recoil, boundary, failure, and quarantine variables
remain part of any KEEP/BREAK attribution test.

The exact disposition is

\[
 \boxed{
 \begin{gathered}
 C=\langle Y\rangle\text{, not }J,
 \qquad {\cal F}_\theta=I_V+{\cal B}C,\\
 \text{commuting EW carrier}+\text{complete FU solder}
 \Longrightarrow
 D_j{\cal F}_\theta=
 {\beta U_d\lambda_{\rm pair}^{\rm net}\over2}{\cal B}GM,\\
 0<\beta<\infty,\ \lambda_{\rm pair}^{\rm net}\ne0
 \Longrightarrow \operatorname{rank}D_j{\cal F}_\theta=6,\\
 \partial_j{\cal F}_\theta
 ={\beta\over2}\partial_\theta^2\langle Q\rangle,\\
 \text{exact ice projection}\Longrightarrow
 \operatorname{rank}G=2,\\
 \text{FY}\Longrightarrow\text{direct-source operator custody only},\\
 \text{full F3 fidelity solder, gluing, Ward response, gravity, and }G
 \text{ remain open.}
 \end{gathered}}                                  \tag{GK30}
\]
