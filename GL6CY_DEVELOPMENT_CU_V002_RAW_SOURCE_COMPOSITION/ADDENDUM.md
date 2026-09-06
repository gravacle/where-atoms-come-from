# Exact raw source composition addendum

## 1. Declared path

Let `m in Z4^3`, let `chi_m(v)=i^(m dot cell(v))`, and put
`sigma(P)=+1`, `sigma(C)=-1`.  The local child frame reverses all half-port
vectors on `C`, so its pair dyads are unchanged while its center `T2` source
changes sign.

For tensor row `s=(xx,yy,zz,xy,xz,yz)`, choose the symmetric half-basis
`A_s` such that `h_s=A_s+A_s^T`, and let `q_s=D_C^* h_s`.  For center
direction `r=(x,y,z)`, the declared probe-plus-literal-anchor path is

\[
 F_v=I+u_s\chi_m(v)^*A_s,
 \qquad b_v={d_v\over a_*}=u_r\chi_m(v)e_r,
\]

\[
 \eta_A(u_s,u_r)
 =u_s\chi_m(v)^*q_{A,s}
  +\eta_A^{\rm phys}(F_v,b_v).
\]

Thus, explicitly,

\[
 A_r=0,\qquad b_s=0,\qquad b_r=\chi_m(v)e_r,
 \qquad b_{sr}=0.
\]

The last equality is a property of this declared bilinear-free anchor path,
not a general assumption.  CY13--CY14 then give

\[
 \eta_{A,r}=-2\sigma(v)\chi_m(v)t_{A,r},
 \qquad
 \boxed{\eta_{A,sr}=2\sigma(v)t_A(A_se_r)\ne0\ \text{generically}.}
\]

The calculation therefore includes both CY02 terms.  This path is an exact
algebraic diagnostic on the inherited finite Q4 incidence structure; it is
not asserted to be an authenticated physical source law.

In pair order `(01,02,03,12,13,23)`, the columns of `D_C^*` are

\[
\begin{pmatrix}
 1/2&-1/2&-1/2&0&0&-1\\
-1/2& 1/2&-1/2&0&-1&0\\
-1/2&-1/2& 1/2&-1&0&0\\
-1/2&-1/2& 1/2& 1&0&0\\
-1/2& 1/2&-1/2&0& 1&0\\
 1/2&-1/2&-1/2&0&0& 1
\end{pmatrix}.
\]

## 2. Exact k=0 matrices

Rows below are `(xx,yy,zz,xy,xz,yz)` and columns are `(x,y,z)`.  Define

\[
 S=\begin{pmatrix}
0&0&0\\
0&0&0\\
0&0&0\\
0&0&-1\\
0&1&0\\
0&0&0
\end{pmatrix}.
\]

The exact physical `-2` center normalization gives

\[
 {\cal R}^{h0}_{sr}=0,\qquad
 {\cal R}^{h2,\rm diag}_{sr}=0,
\]

\[
 \boxed{{\cal R}^{h4,\rm diag}_{sr}={32128\over27}S,}
\]

\[
 \boxed{{\cal R}^{h6,\rm diag}_{sr}
 ={20042848\over6075}S,}
\qquad
 \boxed{{\cal R}^{h6,\rm cyc}_{sr}=-2112S,}
\]

and hence

\[
 \boxed{{\cal R}^{h6,\rm diag+cyc}_{sr}
 ={7212448\over6075}S.}
\]

Every displayed matrix is the sum of the two CY02 terms.  On this declared
path the separately evaluated result is

\[
 \boxed{H_{,A}\eta_{A,sr}=0}
\]

for the bare, `h2`, `h4`, `h6` diagonal, every active `h6` cycle, and every
Q4 character.  This is not because `eta_sr` was set to zero: it is nonzero.
The contraction vanishes by the inherited locked/staggered `T2` structure
and exact alternating-cycle cancellation.  Therefore the displayed nonzero
blocks are supplied entirely by `H_,AB eta_A,s eta_B,r` for this path.

As a regression in the earlier raw-pair/unit-center convention, the exact
matrices are

\[
 H^{h2}_{\rm diag}=0,\quad
 H^{h4}_{\rm diag}={8032\over27}P,\quad
 H^{h6}_{\rm diag}={5010712\over6075}P,\quad
 H^{h6}_{\rm cyc}=-528P,
\]

where

\[
P=\begin{pmatrix}
0&0&0\\0&1&0\\0&0&-1\\0&0&1\\0&-1&0\\0&0&0
\end{pmatrix}.
\]

## 3. Exact finite-range export

For raw coordinate `A=(v,p)` and `B=(w,q)`, the executable aggregates the
complete CU Hessian once by `Delta=cell(w)-cell(v) mod 4`:

\[
\begin{split}
 K^{(n)}_{sr}(\Delta)=\sum_{A,B:\,w-v=\Delta}
 H^{(n)}_{,AB}\,[D_C^*]_{p s}
 [-2\sigma(w)t_{q r}] + \text{the symmetric orientation}.
\end{split}
\]

It then evaluates the exact Laurent/character function

\[
 {\cal R}^{(n)}_{sr}(m)
 =\sum_\Delta K^{(n)}_{sr}(\Delta)i^{m\cdot\Delta}
  +H^{(n)}_{,A}\eta_{A,sr}.
\]

The `--json` interface emits the `h2`, `h4`, `h6` diagonal stencils, the
aggregate 64-cycle stencil, their `h6` sum, the 64 exact character-function
classifications, and per-cycle function/stencil hashes.  These are local
raw-Hamiltonian inputs for overlap gluing.  Pairwise/cyclic seam
cancellation, stationary state/measure, connected response, and the unique
1PI/Schur quotient remain downstream owners.

