# GL6CL — GLOBAL FOURIER PAIR-WRITER AND SMOOTH-ACCESS THEOREM

## Status and scope

This packet derives the exact translation-covariant Fourier symbol of the
complete `GL6CH` order-six `T2` pair-source-to-ring-writer map on the infinite `Q4`
(diamond-incidence) parent.  Parent and child sublattices are retained until
the final common/relative decomposition.  The result proves full common-field
pair-source access at zero momentum and throughout an explicit neighborhood
of it, once the locked diagonal read is included.  It also finds two genuine
ceilings: independent parent/child fields are underdetermined, and the common
tensor writer has finite-momentum rank loss.  The full canonical direct pair
gradient is retained as bookkeeping, but its arbitrary-profile `A1/E` part
is not promoted to a complete effective writer.

This is a global linear Fourier/operator-jet theorem.  It does **not** prove
autonomous generation of the source, a stationary response, a phase, a
continuum limit, spacetime, a metric, `RGRL-B`, Ricci or Einstein dynamics,
gravity, or `G`.

## 1. Inherited tetrahedral coordinate embedding

Use the inherited tetrahedral coordinate link vectors

\[
 T_0=(1,1,1),\quad T_1=(1,-1,-1),\quad
 T_2=(-1,1,-1),\quad T_3=(-1,-1,1),                 \tag{CL01}
\]

so that `T_a^2=3` and `T_a.T_b=-1` for `a != b`.  For
`x in Z^3`, with `e_3=0`, place the two sublattices at

\[
 r(P_x)=\sum_{i=0}^2x_i(T_i-T_3),\qquad
 r(C_y)=r(P_y)+T_3.                                    \tag{CL02}
\]

The link `(x,a)` then joins `P_x` to `C_{x+e_a}` with coordinate
displacement `T_a`.  This embedding has not independently derived physical
distance or spacetime.

Let `a<b<c` be the three ports of a canonical elementary hexagon, and let
`d` be its missing port.  With anchor `x=0`, its six node/pair entries are

\[
\begin{array}{c|c|c}
\text{node}&\text{coordinate position}&\text{cycle-port pair}\\ \hline
P&0&ac\\
C&T_a&ab\\
P&T_a-T_b&bc\\
C&T_a-T_b+T_c&ac\\
P&T_c-T_b&ab\\
C&T_c&bc.
\end{array}                                             \tag{CL03}
\]

Its coordinate center is

\[
 R_{abc}={1\over2}(T_a+T_c-T_b).                       \tag{CL04}
\]

For each of its three pairs, the parent and child sites are exactly opposite
about `R`.  Writing the child displacement from the center as `rho`,

\[
 \rho_{d,ab}={T_a+T_b-T_c\over2},\quad
 \rho_{d,ac}={T_a-T_b+T_c\over2},\quad
 \rho_{d,bc}={-T_a+T_b+T_c\over2}.                    \tag{CL05}
\]

All twelve orientation-pair offsets obey `|rho|^2=11/4`.  The replay derives
(CL03)--(CL05) from node incidence rather than installing the centered
positions by hand.

## 2. Exact parent/child Fourier rows

Use pair order

\[
 {\cal P}=(01,02,03,12,13,23),                         \tag{CL06}
\]

and let `e_p` denote its standard pair vector.  `GL6CH`, under independent
hostile audit, gives the canonical **direct-history** local gradient

\[
 \mu e_p,\qquad \mu={105\over8}{h^6\over U_d^6}.       \tag{CL07}
\]

Its `T2` projection is the complete first tensor-source term at this order.
The arbitrary-profile `A1/E` off-diagonal completion is not classified by
`GL6CH`; only the spatially uniform `A1` direction is separately protected by
the exact denominator identity in Section 8.

After factoring the phase at the ring center, the exact canonical-direct row
for orientation `d` is

\[
 \boxed{
 B_d^P(k)=\sum_{p\subset\bar d}e^{-ik\cdot\rho_{d,p}}e_p^T,
 \qquad
 B_d^C(k)=\sum_{p\subset\bar d}e^{+ik\cdot\rho_{d,p}}e_p^T.}
                                                               \tag{CL08}
\]

Here `p subset bar d` runs over the three pairs among the ports other than
`d`.  Define Fourier amplitudes by

\[
 j_P=j_++j_-,\qquad j_C=j_+-j_-.                       \tag{CL09}
\]

Then the exact common and relative rows are

\[
 \boxed{
 B_d^+(k)=2\sum_{p\subset\bar d}\cos(k\cdot\rho_{d,p})e_p^T,
 \qquad
 B_d^-(k)=-2i\sum_{p\subset\bar d}\sin(k\cdot\rho_{d,p})e_p^T,}
                                                               \tag{CL10}
\]

and the **complete tensor** sourced ring coefficient is

\[
 \boxed{\delta a_d^T(k)=\mu\,[B_d^+(k)P_Tj_+(k)
                         +B_d^-(k)P_Tj_-(k)].}          \tag{CL11}
\]

Thus common tensor source profiles are even in `k`, while relative
parent/child tensor source profiles are odd.  Equation (CL11) is an externally charted
source derivative; it is not an autonomous equation of motion.  Equations
(CL08)--(CL10) without `P_T` are canonical-direct bookkeeping and must not be
read as a complete arbitrary-profile `A1/E` writer.

## 3. Exact zero-mode rank and reconstruction

Let

\[
 A=(1,1,1,1,1,1),\qquad
 t_1=e_{01}-e_{23},\quad t_2=e_{02}-e_{13},\quad
 t_3=e_{03}-e_{12},                                    \tag{CL12}
\]

and use the orthogonal pair decomposition
`R^6=A1 direct-sum E direct-sum T2`, with projectors `P_A,P_E,P_T`.
At `k=0`,

\[
 B_d^+(0)=2\sum_{p\subset\bar d}e_p^T,qquad
 \boxed{B_+(0)^*B_+(0)=24P_A+8P_T.}                   \tag{CL13}
\]

Thus the canonical-direct row has rank four on `A1+T2`.  The complete tensor
writer established at arbitrary profile is its tensor projection,

\[
 W_T(k)=B_+(k)P_T,
 \qquad \boxed{W_T(0)^*W_T(0)=8P_T,}                  \tag{CL13a}
\]

which has rank three and null sector `A1+E`.  For compatible four-component
tensor-writer data `w`,

\[
 \boxed{j_T={1\over8}P_TB_+(0)^*w.}                   \tag{CL14}
\]

Independently re-deriving the six locked local words gives the diagonal read
`D` and

\[
 D^*D=4P_A+16P_E.                                      \tag{CL15}
\]

For the common smooth field, stack the locked read and complete tensor writer
as `C(k)=(D;B_+(k)P_T)`.  At zero momentum,

\[
 \boxed{C(0)^*C(0)=4P_A+16P_E+8P_T,}                  \tag{CL16}
\]

which has rank six and determinant `524,288`.  Its explicit zero-mode left
inverse is

\[
 \boxed{j=\left({1\over4}P_A+{1\over16}P_E+{1\over8}P_T\right)
                 [D^*d+P_TB_+(0)^*w].}                \tag{CL17}
\]

This is the same-parent six-direction access result promoted from a single
node to the translation-invariant zero mode without using unclassified
off-diagonal `A1/E` pieces.

## 4. Rigorous smooth-field neighborhood

The exact even expansion is

\[
 B_d^+(k)=2\sum_pe_p^T
 -\sum_p(k\cdot\rho_{d,p})^2e_p^T
 +{1\over12}\sum_p(k\cdot\rho_{d,p})^4e_p^T+O(|k|^6).
                                                               \tag{CL18}
\]

Since there are twelve entries, every `|rho|^2=11/4`, and
`|2(cos x-1)|<=x^2`,

\[
 \|B_+(k)-B_+(0)\|_2
 \le \|B_+(k)-B_+(0)\|_F
 \le {11\sqrt3\over2}|k|^2.                            \tag{CL19}
\]

The smallest nonzero tensor singular value in (CL13a) is `sqrt(8)`.  Hence
both the common tensor-writer rank three and the locked-read-plus-tensor-
writer rank six are rigorously preserved whenever

\[
 \boxed{|k|^4<{32\over363}.}                           \tag{CL20}
\]

On this open ball the exact analytic left inverse is

\[
 \boxed{L(k)=[C(k)^*C(k)]^{-1}C(k)^*.}                 \tag{CL21}
\]

Thus six-direction access does persist for sufficiently slowly varying
common pair-source profiles: `A1+E` through the locked read and `T2` through
the complete writer.  No continuum or stationary-response limit is needed
for that rank statement.

## 5. Exact cubic `T2` block; rotational diagnosis remains open

The twelve centered offsets have an isotropic aggregate second moment,

\[
 \sum_{d,p}\rho_{d,p}\rho_{d,p}^T=11I_3,              \tag{CL22}
\]

but this scalar identity does not make the tensor writer rotationally
invariant.  In the orthonormal `T2` basis `t_i/sqrt(2)`, its normal operator
through quadratic order is

\[
 \boxed{
 N_T(k)=8I-2|k|^2I+12kk^T
 -28\,\operatorname{diag}(k_x^2,k_y^2,k_z^2)+O(|k|^4).}
                                                               \tag{CL23}
\]

The last term displays the cubic structure of the `T2` block.  But `T2` is
the three-dimensional off-diagonal part of the five-dimensional traceless
symmetric tensor space `E2+T2`; it is not closed under `SO(3)`.  A
`diag(k_i^2)` term can occur in the `T2-T2` restriction of an `SO(3)`-
covariant rank-two-tensor operator, with `E2-T2` and `E2-E2` blocks completing
the transformation law.  Equation (CL23) therefore does **not**, by itself,
prove physical rotational anisotropy or isotropy.  It is an exact diagnostic
that must be compared with a consistently normalized full `E2+T2`
completion.

There is a useful but weaker canonical-direct scalar diagnostic.  In the
unnormalized `(A,t_1,t_2,t_3)` basis,

\[
\begin{aligned}
 \det B_+(k)|_{A+T}
 &=768-1408|k|^2+1072|k|^4\\
 &\quad-{416\over3}(k_x^4+k_y^4+k_z^4)+O(|k|^6).
                                                               \tag{CL24}
\end{aligned}
\]

Its quadratic scalar is rotationally invariant as a polynomial, while a
cubic invariant appears at fourth order.  This canonical-direct determinant
includes the separately typed `A1` bookkeeping row and is not a completed
five-shear response.  The raw fourth offset contraction is

\[
 \sum_{d,p}(k\cdot\rho_{d,p})^4
 ={57\over4}|k|^4+{13\over2}\sum_i k_i^4.              \tag{CL25}
\]

Equations (CL23)--(CL25) prevent either an isotropic or anisotropic physical
conclusion from being inferred from rank or from the `T2` block alone.

## 6. Relative parent/child sector and the soldering obstruction

At zero momentum `B_-(0)=0`: a ring cannot see a uniform source that changes
sign between the two sublattices.  The leading `T2` map is

\[
 B_-(k)|_{T_2}=-2i\,{\cal L}(k)+O(|k|^3),              \tag{CL26}
\]

where, after removing the displayed factor, the four rows and three columns
of `L` are the dot products with

\[
{1\over2}\begin{pmatrix}
(3,-1,-1)&(-1,3,-1)&(-1,-1,3)\\
(3,1,1)&(1,3,-1)&(1,-1,3)\\
(3,1,-1)&(1,3,1)&(-1,1,3)\\
(3,-1,1)&(-1,3,1)&(1,1,3)
\end{pmatrix}.                                         \tag{CL27}
\]

The exact Cauchy--Binet sum of squares of its four `3x3` minors is

\[
\begin{aligned}
 S(k)&=9(k_x^6+k_y^6+k_z^6)
 -9\sum_{i\ne j}k_i^4k_j^2+58k_x^2k_y^2k_z^2.        \tag{CL28}
\end{aligned}
\]

A Hadamard recombination of those four minors gives, up to irrelevant
overall signs,

\[
 6k_y(k_x^2+k_z^2-k_y^2),\quad
 6k_z(k_x^2+k_y^2-k_z^2),\quad
 4k_xk_yk_z,\quad
 6k_x(k_y^2+k_z^2-k_x^2).                              \tag{CL29}
\]

Therefore the leading relative `T2` map has rank three generically, rank two
on the six nonzero Cartesian face-diagonal directions

\[
 k\parallel(1,\pm1,0),\ (1,0,\pm1),\ (0,1,\pm1),      \tag{CL30}
\]

and rank zero at `k=0`.  The rank-two face-diagonal dependence is not merely
leading-order: the exact sine symbol has one dependent pair of `T2` columns
for every momentum along each line.  For example, along `k=q(1,1,0)`, its
first two columns coincide and the null direction is `t_1-t_2`.

There is also a dimension obstruction independent of this special-direction
calculation.  If `j_P` and `j_C` are unrelated six-vectors, separate locked
reads supply at most `3+3` outputs and the ring writer supplies four.  The
combined map from twelve source directions has rank at most ten at every
momentum.  At zero it has rank exactly nine, with all three relative `T2`
directions in the kernel.  Thus an arbitrary-support inverse requires either
a physical common-field/soldering law or additional independent writer
channels.  Incidence alone does not supply that law.

## 7. Exact finite-momentum rank loss

Let reciprocal coordinates be `q_i=k.(T_i-T_3)`.  At the Brillouin-boundary
point

\[
 q=(\pi,0,0),\qquad k={\pi\over4}(1,1,1),              \tag{CL31}
\]

the three common `T2` columns are identical.  In terms of
`C=cos(pi/8)` and `S=cos(3pi/8)`, their four rows are proportional to

\[
 (-2C,-2C,-2C),\quad(2S,2S,2S),\quad
 (2S,2S,2S),\quad(2S,2S,2S).                           \tag{CL32}
\]

The complete common tensor writer therefore has rank one there.  This exact example
rules out a full-Brillouin-zone six-direction inverse.  It does not affect
the rigorous smooth-field ball (CL20); it precisely limits the result to
long-wavelength/common-field gluing.

## 8. Uniform scalar source links storage energy to the writer

For a local word of occupancy `n`,

\[
 \boxed{\sum_{a<b}Z_aZ_b=2(n-2)^2-2.}                 \tag{CL33}
\]

Consequently a uniform pair source `qA` changes the constraint coefficient
as \(U_d \to U_d+2q\), plus the scalar `-2q`.  Differentiating the source-free
ring amplitude gives

\[
 {d\over dq}\left[-{63\over8}{h^6\over(U_d+2q)^5}\right]_{q=0}
 ={315\over4}{h^6\over U_d^6}.                         \tag{CL34}
\]

This equals the sum of the six canonical `GL6CH` scalar vertices,

\[
 6\left({105\over8}{h^6\over U_d^6}\right)
 ={315\over4}{h^6\over U_d^6}.                        \tag{CL35}
\]

Thus the same microscopic pair operator that measures local storage-energy
cost also changes a later ring-transition amplitude.  This is an exact
same-parent storage-energy/future-writer linkage.  It still does not show
that records autonomously generate the source or that the resulting writer
forms a reciprocal stationary bulk response.

## 9. What is established and what remains

Established here:

1. the exact inherited-coordinate parent/child Fourier symbol of the audited `GL6CH`
   complete `T2` writer, with the larger canonical-direct row separately
   typed;
2. rank-three common `T2` writer access at zero momentum;
3. rank-six pair access after composing that writer with the locked `A1+E`
   read;
4. an explicit analytic inverse on a rigorous nonzero smooth-field ball;
5. the exact cubic form of the common `T2` block, the need for a full
   `E2+T2` rotational diagnostic, and finite-momentum rank loss;
6. the relative-sublattice rank structure and a precise soldering
   obstruction; and
7. an exact uniform-source consistency between storage energy and future
   ring writing.

Still open are arbitrary-profile order-six `A1/E` off-diagonal completion,
the physical law that makes the source autonomous, the law that relates or
dynamically solders parent and child fields, the consistently normalized
`E2+T2` completion needed to interpret the cubic `T2` block, a stationary
reciprocal response, refinement to a common bulk limit, and every later
identification with geometry or gravity.
