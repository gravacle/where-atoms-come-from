# Complete homogeneous H6 Coulomb-DPAR source response on the FO component

**Lane:** `LANE_GRA_FX_F3_Q4_COULOMB_COMPLETE_H6_SOURCE_RESPONSE_V001`  
**Date:** 2026-08-28  
**Status:** exact conditional finite-component theorem  
**Claim class:** complete homogeneous (`k=0`) nonidentity source modulo
`H_id`, on the selected 180-state FO winding component, fixed through H6

## 1. Scope and inherited premises

This theorem is conditional on FV `S1`--`S10`, including `FV-PURE`, and on
the frozen FO 180-state translation-closed H6 ice component.  FW already
proved the exact finite response of the displayed pair-plus-irreducible-ring
source, `FV-WITNESS`.  FV13 left one physics term explicitly unevaluated:

`Q_diag^(2,4,6)`, the derivative of every diagonal closed word and every
Brillouin-Wigner/Feshbach fold through sixth order.

This lane computes that term and composes it with FW.  Here **complete** means
only the homogeneous source on this selected finite winding component,
through H6, and modulo Hilbert-space identities arising from `H_id`.  It does
not mean the local or nonzero-momentum source, every ice component, H8 and
higher orders, a closed-time-path (CTP) response, a Ward theorem, RGRL-B, a
thermodynamic massless tensor, gravity, or Newton's constant.

Set

\[
 x=h/U_d,\qquad J_6={63\over8}U_dx^6,\qquad
 \rho={U_d\over J_6}={8\over63x^6}.
 \tag{FX01}
\]

The strain source is inserted before reduction, the endpoint reference is
`E_R=0`, and

\[
 Q^{ij}=-2{\partial H\over\partial j_{ij}}\bigg|_{j=0}.
 \tag{FX02}
\]

The six stored source covector coordinates are
`(xx,yy,zz,2xy,2xz,2yz)`: the factor two is the one required by
`j:R` for symmetric tensors.

## 2. Exhaustive diagonal-word theorem

For a length-`m` flip word to return to the same bit state, every used edge
must occur an even number of times.  Therefore the complete diagonal
partitions through order six are

\[
 m=2:(2),\qquad m=4:(4),(2,2),\qquad
 m=6:(6),(4,2),(2,2,2).                           \tag{FX03}
\]

The one-edge `(4)` and `(6)` words return to `P` after their second flip and
have no irreducible order.  For the remaining multisets the exact unique
order counts before removal of intermediate `P` returns are respectively
`1,6,15,15,90` for `(2)`, `(2,2)`, `(4,2)`, `(2,4)`, and `(2,2,2)`.

This census is complete.  The FO graph has no two- or four-link ice cycle,
so a proper nonzero parity prefix below length six cannot be a second ice
state.  At length six, the only nonzero ice endpoint is an elementary
hexagon with six distinct edges; that is precisely the off-diagonal
irreducible ring already computed by FV/FW.  Every diagonal word is
therefore in (FX03), and every intermediate `P` return is a fold rather than
an omitted irreducible history.

For selected edges and a proper prefix `S`, use the endpoint-referenced gap

\[
 \Delta_n(S;j)=E_{n\triangle S}(j)-E_n(j).
 \tag{FX04}
\]

At ideal Coulomb slope `lambda=-1/2`, `U_d=1`, the exact gap derivative is
obtained from

\[
 {\partial E_C\over\partial j}\bigg|_0
 ={1\over4}\sum_{v,a<b}\widehat R_{ab}Z_{v,a}Z_{v,b},
 \tag{FX05}
\]

while differentiating all hopping numerators gives

\[
 {\partial\over\partial j}\prod_{r=1}^{m}
 \left(1-{1\over2}j:D_{a_r}\right)\bigg|_0
 =-{1\over2}\sum_{r=1}^{m}D_{a_r}.               \tag{FX06}
\]

The verifier evaluates every gap and derivative with exact rational
arithmetic.  It independently compares the closed resolvent-derivative
formula with sequential dual-polynomial multiplication for every
multiplicity family.

## 3. Complete BW/Feshbach folds

Write the endpoint-referenced diagonal kernel as

\[
 \delta=h^2k_2(\delta)+h^4k_4(\delta)+h^6k_6(\delta)+O(h^8),
 \qquad k_m(\delta)=k_{m0}+k_{m1}\delta+k_{m2}\delta^2+\cdots .
 \tag{FX07}
\]

Putting `delta=h^2 a_2+h^4 a_4+h^6 a_6+...` restores every fold through H6:

\[
\begin{aligned}
 a_2&=k_{20},\\
 a_4&=k_{40}+a_2k_{21},\\
 a_6&=k_{60}+a_2k_{41}+a_4k_{21}+a_2^2k_{22}.    \tag{FX08}
\end{aligned}
\]

Differentiation gives

\[
\begin{aligned}
 da_2={}&dk_{20},\\
 da_4={}&dk_{40}+da_2k_{21}+a_2dk_{21},\\
 da_6={}&dk_{60}+da_2k_{41}+a_2dk_{41}
 +da_4k_{21}+a_4dk_{21}\\
 &+2a_2da_2k_{22}+a_2^2dk_{22}.                 \tag{FX09}
\end{aligned}
\]

Equations (FX03), (FX08), and (FX09) jointly cover all irreducible diagonal
words and all folds; no static coefficient rank is substituted for a
retarded rank.

## 4. Exact orbit result

FO translation preserves each tetrahedral edge label, graph incidence, and
the local labeled occupations.  Its six orbits all have length thirty, so a
homogeneous diagonal operator is determined by one representative per orbit
with normalized zero-momentum weight `1/sqrt(30)`.

The source-off folded energy is scalar on all six orbits:

\[
 \boxed{a_2=-60,\qquad a_4=-35,\qquad a_6=-{893\over9}.}       \tag{FX10}
\]

In orbit order `0,...,5`, the exact derivative rows are

\[
\begin{array}{c|c|c|c}
o&da_2&da_4&da_6\\ \hline
0&(46,22,22)&(541/6,97/6,97/6)&(635503/1350,50611/1350,50611/1350)\\
1&(46,22,22)&(541/6,97/6,97/6)&(635503/1350,50611/1350,50611/1350)\\
2&(46,21,23)&(541/6,157/12,77/4)&(635503/1350,52481/2700,149963/2700)\\
3&(46,23,21)&(541/6,77/4,157/12)&(635503/1350,149963/2700,52481/2700)\\
4&(46,22,22)&(541/6,97/6,97/6)&(635503/1350,50611/1350,50611/1350)\\
5&(46,22,22)&(541/6,97/6,97/6)&(635503/1350,50611/1350,50611/1350)
\end{array}                                                     \tag{FX11}
\]

The three omitted off-diagonal columns are exactly zero.  The exhaustive
replay contains 13,725 distinct signature/family weights; the pair/triple
signature counts on the six orbits are

\[
 (325,5444),(362,6021),(275,4793),(277,4823),(360,5968),(316,5220).
 \tag{FX12}
\]

## 5. Exact complete-source reduction

Let

\[
 D^{ij}:={Q_{\rm pair}^{(0),ij}\over U_d},\qquad
 \mathbb I^{ij}:=\delta^{ij}I_{\cal H}.
\]

Applying (FX02) to (FX11) gives the exact 180-state operator identities

\[
\begin{aligned}
 Q_{\rm diag}^{(2)}&=U_dx^2\left(-D-40\mathbb I\right),\\
 Q_{\rm diag}^{(4)}&=U_dx^4\left(-{37\over12}D-20\mathbb I\right),\\
 Q_{\rm diag}^{(6)}&=U_dx^6\left(-{16247\over900}D-{374\over135}\mathbb I\right).
                                                               \tag{FX13}
\end{aligned}
\]

Thus the generated diagonal source creates no new nonidentity direction.  It
renormalizes the direct `E` source and adds only a Hilbert identity.  With
`R=Q_ring/J_6`, the complete homogeneous source through H6 is

\[
 {Q_{\rm eff}^{(\le6)}\over J_6}
 =\rho_E(x)D+R+\rho s(x)\mathbb I,                \tag{FX14}
\]

where

\[
\boxed{
 f_E(x)=1-x^2-{37\over12}x^4-{16247\over900}x^6,
 \qquad \rho_E(x)=\rho f_E(x),}
                                                               \tag{FX15}
\]

and `s(x)=-40x^2-20x^4-(374/135)x^6`.

## 6. Dynamic response theorem

### Theorem `CH6SR-1`

As a formal through-H6 series, `f_E(x)` has constant term one and is a
**formal power-series unit**.  Therefore the direct `E` witness cannot be
cancelled order by order.  Exact composition of (FX14) with the audited FW
algebra proves, generically,

\[
 \boxed{5\ \longrightarrow\ 3\ \longrightarrow\ 2\ \longrightarrow\ 2}
                                                               \tag{FX16}
\]

for operator rank modulo identity, source-to-commutator rank, unique-ground
retarded/spectral rank, and first nonzero commutator-moment rank.

The two poles remain

\[
 \Delta_1=2+2\sqrt2,\qquad \Delta_2=4+2\sqrt2,                 \tag{FX17}
\]

with rank-one residue vectors

\[
 r_1=\left(0,{\rho_E\over\sqrt2},
 -\rho_E\sqrt{3\over2},-{3\over\sqrt2},-{3\over\sqrt2},0\right),
 \quad
 r_2=\left(0,0,0,{3\over\sqrt2},-{3\over\sqrt2},0\right).
                                                               \tag{FX18}
\]

There is no third ground-state pole through H6.  This exact response proof is
the exact pair-plus-identity reduction (FX13)--(FX15) composed with FW's
audited algebraic projectors and residues.  The independent 180-state NumPy
commutator/eigensystem/Lehmann calculations are replay checksums, not the
foundation of the exact claim.

## 7. The finite polynomial root is not a threshold

If the finite polynomial (FX15) is evaluated as though it were an ordinary
function, it has one positive root.  Put `y=x^2`:

\[
 p(y)=1-y-{37\over12}y^2-{16247\over900}y^3,
\quad
 p'(y)=-1-{37\over6}y-{16247\over300}y^2<0\quad(y\ge0).         \tag{FX19}
\]

Moreover,

\[
 p(1/4)={15853\over57600}>0,
\qquad
 p(729/2500)={-2157513587\over1562500000000}<0.                \tag{FX20}
\]

Hence there is exactly one root with
`1/2<x_0<27/50`; numerically `x_0=0.5398271903...`.  At the algebraic
**finite through-H6 cancellation** stratum, the ranks are `4 -> 2 -> 2 -> 2`
because the two ring-pole directions survive.

This is **not a physical threshold** or critical point.  No convergence
radius has been established near that value, and H8 or higher terms can
shift or remove the zero.  It cannot be used as a gravity-onset prediction.

## 8. Scientific meaning and ceiling

The missing FV diagonal terms are dynamically real, but on this exact
component they do not enlarge the response space: they dress the same
record-compatible `E` channel already supplied by the Coulomb pair source.
The finite rank-two ground response therefore survives completion of every
diagonal source term and fold through H6.

This advances the no-lab gravity program by closing a previously explicit
microscopic source-completion gap.  It does not establish nonzero-momentum
locality, a CTP/Ward identity, a thermodynamic massless tensor, RGRL-B,
gravity, a value of `G`, or Newton's constant.
