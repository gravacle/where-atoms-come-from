# S4 lineage-to-metric response-kernel reduction theorem

**Theorem ID:** `RGRL-LMRK-V001`

**Date:** 2026-08-27

**Claim class:** exact finite-dimensional representation theorem for the q4
pair-memory tangent; exact point-local/zero-spatial-momentum response-kernel
count at an `S4`-fixed background; exact local scalar-gamma and record-count
underdetermination theorem; conditional causal, passivity, and Ward ceilings

**Status:**
`EXACT_S4_1_PLUS_2_PLUS_3_PAIR_MEMORY_TO_TRACE_AND_SHEAR_DECOMPOSITION__THREE_POINT_LOCAL_RESPONSE_FORM_FACTORS_OR_TWO_UNDER_FULL_O3_ISOTROPY__SCALAR_GAMMA_AND_RECORD_COUNT_DO_NOT_FIX_THE_KERNEL__PHYSICAL_SOLDERING_ABSOLUTE_SCALE_SPATIAL_DISPERSION_AND_NUMERICAL_G_REMAIN_INDEPENDENT`

**Not claimed:** that the EW Fisher tensor is physical space; that the six
pair coordinates are physical fields without RGRL/GSGB soldering; that three
form factors exhaust a nonzero-momentum dispersive kernel; that `S4` implies
full rotational invariance; that causality, passivity, or the Ward identity is
supplied by Fisher positivity; a numerical value of any form factor, length
scale, Einstein--Hilbert stiffness, or Newton's constant.

## 1. Exact scope and notation

Let

\[
 V=\mathbf 1^\perp\subset\mathbb R^4,\qquad \dim V=3,
 \qquad v_a=\left(I_4-\frac14\mathbf1\mathbf1^{\mathsf T}\right)e_a,
 \tag{LM01}
\]

and let

\[
 \mathscr E:=\mathbb R^{\mathcal E_4},\qquad
 \mathcal E_4=\{12,13,14,23,24,34\}                 \tag{LM02}
\]

be the six-dimensional unordered-pair space.  The exact EW kinematic
Jacobian is

\[
 D:\mathscr E\longrightarrow\operatorname{Sym}^2(V),
 \qquad D(e_{ab})=v_a\odot v_b,
 \qquad u\odot v:=uv^{\mathsf T}+vu^{\mathsf T}.   \tag{LM03}
\]

It is `S4`-equivariant and invertible.  At EW's symmetric point
`(theta,J)=(0,0)`, it is the exact tangent derivative of the localization
Fisher tensor with respect to the pair-memory coordinate:

\[
 \delta F_X\big|_{J=0}=D\,\delta J.                \tag{LM04}
\]

Equation (LM04) is information-coordinate kinematics.  It becomes a physical
spatial-metric variation only after an independently typed soldering and
scale join.  This distinction is load-bearing below.

The response count in section 4 concerns the derivative-zero constitutive
term, equivalently a point-local kernel or its zero-spatial-momentum value
`k=0`, at an `S4`-fixed background.  A spatial wavevector transforming in
`V` is not a spectator and can permit additional tensor structures.

## 2. Theorem LMRK-1 -- exact `1+2+3` decomposition

Let `O` be the opposite-edge involution

\[
 12\leftrightarrow34,\qquad
 13\leftrightarrow24,\qquad
 14\leftrightarrow23,                              \tag{LM05}
\]

and let `\mathbb J_6` denote the six-by-six all-ones matrix.  Define

\[
 P_A={\mathbb J_6\over6},\qquad
 P_T={I-O\over2},\qquad
 P_E={I+O\over2}-P_A.                              \tag{LM06}
\]

Then these are mutually orthogonal projectors and

\[
 \boxed{
 \mathscr E=A_1\oplus E_2\oplus T_2,
 \qquad \dim(A_1,E_2,T_2)=(1,2,3).}               \tag{LM07}
\]

The target has the matching multiplicity-free decomposition

\[
 \boxed{
 \operatorname{Sym}^2(V)=A_1\oplus E_2\oplus T_2
 =\mathbb R I_V\oplus\operatorname{STF}(V).}       \tag{LM08}
\]

Here `E_2` and `T_2` are two inequivalent `S4` irreducibles, while together
they make the five-dimensional symmetric-trace-free space after a full
rotational structure is supplied.  More exactly,

\[
 DP_A=P_{\rm tr}D,
 \qquad D(P_E+P_T)=P_{\rm STF}D.                   \tag{LM09}
\]

An explicit tetrahedral coframe makes the content transparent.  Choose an
orthogonal identification `V\simeq\mathbb R^3` for which

\[
 2v_1=(1,1,1),\quad 2v_2=(1,-1,-1),\quad
 2v_3=(-1,1,-1),\quad 2v_4=(-1,-1,1).             \tag{LM10}
\]

For `x\in\mathscr E`, set

\[
\begin{aligned}
 p&=(x_{12}+x_{34},\ x_{13}+x_{24},\ x_{14}+x_{23}),\\
 d&=(x_{12}-x_{34},\ x_{13}-x_{24},\ x_{14}-x_{23}),\\
 s&=p_1+p_2+p_3=\sum_{a<b}x_{ab},
 \qquad p'_i=p_i-{s\over3}.
\end{aligned}                                      \tag{LM11}
\]

In that orthonormal coframe,

\[
 \boxed{D(P_Ax)=-{s\over6}I_3,}                   \tag{LM12}
\]

\[
 \boxed{D(P_Ex)=\operatorname{diag}(p'_1,p'_2,p'_3),} \tag{LM13}
\]

and

\[
 \boxed{
 [D(P_Tx)]_{yz}=-{d_1\over2},\quad
 [D(P_Tx)]_{xz}=-{d_2\over2},\quad
 [D(P_Tx)]_{xy}=-{d_3\over2},}
 \tag{LM14}
\]

with zero diagonal in (LM14).  Thus the six pair channels resolve exactly
into one scale/trace mode, two tetrahedral diagonal-shear modes, and three
tetrahedral off-diagonal-shear modes.

### Proof

The `+1` eigenspace of `O` consists of the three opposite-pair sums and the
`-1` eigenspace consists of the three opposite-pair differences.  Removing
the uniform vector from the former gives dimensions `1+2+3`, and direct
substitution proves (LM06).  The standard `S4` character decomposition of
the edge permutation representation is `A1+E2+T2`; the symmetric square of
the standard three-dimensional contrast representation has the same
decomposition.  Equivariance and invertibility of `D` identify the matching
summands.  Substitution of (LM10) into (LM03) gives (LM12)--(LM14).  QED.

## 3. Theorem LMRK-2 -- exact pullback Gram and sector conditioning

Give `\mathscr E` its standard Euclidean inner product and
`\operatorname{Sym}^2(V)` its Frobenius inner product.  Then

\[
 G_D:=D^*D
 ={5\over4}I-{1\over4}A+{1\over4}O,              \tag{LM15}
\]

where `A` is the adjacency matrix of the line graph of the tetrahedron: two
edges are adjacent exactly when they share one endpoint.  Equivalently, the
Gram entry is

\[
 \langle D e,D f\rangle_F=
 \begin{cases}
  5/4,&e=f,\\
 -1/4,&e\ne f\text{ share a vertex},\\
  1/4,&e\text{ and }f\text{ are opposite}.
 \end{cases}                                       \tag{LM16}
\]

The three exact sector eigenvalues are

\[
 \boxed{
 G_D={1\over2}P_A+2P_E+P_T.}                      \tag{LM17}
\]

Hence the singular values of `D` are

\[
 \boxed{2^{-1/2}\ (A_1),\qquad 2^{1/2}\ (E_2),
 \qquad 1\ (T_2).}                               \tag{LM18}
\]

This proves invertibility while also proving that the raw equal-coordinate
norm on the six pair channels is not the Frobenius metric norm: the three
sectors have exact relative weights `1/2:2:1`.

### Proof

For `D(e_{ab})=v_a\odot v_b`,

\[
 \langle u\odot v,x\odot y\rangle_F
 =2[(u\cdot x)(v\cdot y)+(u\cdot y)(v\cdot x)].   \tag{LM19}
\]

Using `v_a\cdot v_b=\delta_{ab}-1/4` gives (LM16), hence (LM15).
The line-graph adjacency eigenvalues on `(A1,E2,T2)` are `(4,-2,0)`
and the opposite involution eigenvalues are `(+1,+1,-1)`.  Substitution
gives (LM17)--(LM18).  QED.

## 4. Theorem LMRK-3 -- the dynamical kernel has three `S4` form factors

Let `\delta L` denote a complete six-channel lineage intervention tangent,
carrying the same natural edge action of `S4` as `\mathscr E`, and let `H^R`
be the independently physical, retarded lineage-to-pair response:

\[
 \delta L\ \xrightarrow{\ H^R\ }\ \delta J
 \ \xrightarrow{\ D\ }\ \delta F_X.              \tag{LM20}
\]

The exact factorization is therefore

\[
 \boxed{
 \delta F_X=D H^R\delta L.}                       \tag{LM21}
\]

After a valid physical soldering in the QFI convention, and at fixed
independently calibrated `ell_F`, it becomes

\[
 \boxed{
 \delta s_{\rm sp}=\ell_F^2D H^R\delta L.}        \tag{LM22}
\]

Neither `H^R` nor `\ell_F` is part of the algebraic Jacobian `D`.

At an `S4`-fixed background, suppose the derivative-zero or `k=0` kernel is
`S4`-equivariant.  Real Schur reduction of the multiplicity-free decomposition
(LM07) gives exactly

\[
 \boxed{
 H^R(\omega,0)
 =h_A^R(\omega)P_A+h_E^R(\omega)P_E+h_T^R(\omega)P_T.} \tag{LM23}
\]

Consequently the direct lineage-to-information-metric kernel has exactly the
same three-function freedom:

\[
 \boxed{
 K^R(\omega,0)
 =D\,[h_A^RP_A+h_E^RP_E+h_T^RP_T].}               \tag{LM24}
\]

Equivalently, every `S4`-equivariant edge kernel has matrix entries `a` on the
same edge, `b` on adjacent edges, and `c` on opposite edges, with sector
eigenvalues

\[
 \boxed{
 h_A=a+4b+c,\qquad h_E=a-2b+c,\qquad h_T=a-c.}    \tag{LM25}
\]

If an independently established physical `O(3)` action makes the response
isotropic in the coframe of (LM10), it must act compatibly on both sides of
the kernel.  Explicitly, require the lineage/source tangent to carry the
pulled-back action

\[
 \rho_{\mathscr E}(R)x
 :=D^{-1}\!\left[R(Dx)R^{\mathsf T}\right],
 \qquad R\in O(3),                                 \tag{LM26a}
\]

and require `K^R` (equivalently `H^R` under `D`) to intertwine that action
with the ordinary tensor action on the output.  The five-dimensional
trace-free tensor is then one irreducible `l=2` sector.  In the `D`-relative
normalization of (LM23)--(LM24), isotropy requires

\[
 \boxed{h_E^R=h_T^R=:h_{\rm sh}^R,}               \tag{LM26}
\]

and only two point-local functions remain: a trace response `h_tr=h_A` and a
shear response `h_sh`.  Equation (LM26) is an added rotational conclusion
with the compatible-source premise (LM26a); it does not follow from
tetrahedral `S4` symmetry or from an `O(3)` action on the output alone.

At nonzero spatial momentum, tensors made from `k_i` can split longitudinal
and transverse responses and generate additional allowed structures.  Thus
(LM23) is not a three-function classification of the full dispersive kernel.
Gradient propagation must be supplied and classified independently, for
example at the EX/Einstein endpoint.

### Proof

Every equivariant endomorphism preserves the three inequivalent irreducible
summands of (LM07), and acts as one scalar convolution kernel on each (one
generally complex scalar value per sector after frequency transformation).
Since `D` is an
equivariant isomorphism, the same statement applies to maps from
`\mathscr E` to `Sym^2(V)` after factoring out `D`.  Counting the three edge
orbits of a stabilizer gives the equivalent `(a,b,c)` form and (LM25).
Under the compatible actions in (LM26a),
`Sym^2(V)=l=0\oplus l=2` on both sides of the response.  Schur reduction then
joins `E2` and `T2`, proving (LM26).  QED.

## 5. Theorem LMRK-4 -- scalar gamma and record count cannot fix the kernel

At the symmetric point `J=0` of the complete-query EW pair family, the six
pair characters are orthonormal in covariance, so

\[
 F_J(0)=I_6.                                       \tag{LM27}
\]

For EW's squared-fidelity complete-query gamma convention,

\[
 \boxed{
 -\log\gamma_Q(0,\delta J)
 ={1\over4}\|\delta J\|^2+O(\|\delta J\|^3),
 \qquad
 \gamma_Q(0,\delta J)
 =1-{1\over4}\|\delta J\|^2+O(\|\delta J\|^3).} \tag{LM28}
\]

Thus the quadratic scalar gamma datum leaves a five-sphere of pair-memory
directions.  Unit vectors chosen separately in `A1`, `E2`, and `T2` have the
same gamma decrement through second order, but (LM12)--(LM14) send them to a
scale deformation, diagonal shear, and off-diagonal shear, respectively.
Their squared metric-output norms are moreover weighted by `1/2`, `2`, and
`1` by (LM17).

More strongly, even complete knowledge of this scalar gamma family and its
record count does not determine the dynamical law `H^R`.  Let `f^R` be any
nonzero retarded scalar response already admitted by the chosen causal and
passive class.  The two laws, sector-nondegenerate at every frequency where
`f^R` is nonzero, are

\[
 H_1^R=f^R(P_A+P_E+P_T),                           \tag{LM29}
\]

\[
 H_2^R=f^R(2P_A+3P_E+4P_T)                        \tag{LM30}
\]

Both are `S4`-equivariant and preserve retardedness; with the work convention
in section 6, positive sector multipliers preserve passivity.  They have
different metric-response kernels, but they can be attached to exactly the
same EW state/query family, the same complete-query gamma, and the same
record count.  Under full `O(3)` one may instead compare `H_1` with
`f^R[2P_A+3(P_E+P_T)]`; the two-function underdetermination remains.

No KMS/fluctuation--dissipation same-observable axiom is assumed in this
counterexample: the EW state/query family does not specify a Hamiltonian,
physical source, or response observable.  Supplying those additional data can
constrain a dissipative spectrum, but that is precisely information beyond
scalar gamma and count.

For `N` actual independent equal carriers,

\[
 F_J^{(N)}=NF_J^{(1)},\qquad
 \gamma_Q^{(N)}=(\gamma_Q^{(1)})^N.               \tag{LM31}
\]

The count fixes this accumulation factor, not the direction in pair space,
the sector response functions, or their physical normalization.  Adopted
RGRL-C supplies its declared nonzero/full-rank lineage ancestry statement;
it does not supply the amplitudes in (LM23).

Equations (LM28)--(LM31) prove the exact ceiling:

\[
 \boxed{
 \text{scalar gamma or record count}
 \ \not\Longrightarrow\ H^R, K^R,
 \text{ curvature, or }G.}                        \tag{LM32}
\]

The first nonimplication is already exact by the two same-gamma response
models (LM29)--(LM30); the local five-sphere calculation identifies the
missing directional information explicitly.

## 6. What causality, passivity, and Ward custody do constrain

### 6.1 Retarded causality

Once an operational time and causal-response class have been independently
identified, retardedness requires, sector by sector,

\[
 h_r^R(t)=0\quad(t<0),\qquad r\in\{A,E,T\}.        \tag{LM33}
\]

In a standard tempered/stable causal-response class, with Fourier--Laplace convention
`h(\omega)=\int_0^\infty e^{i\omega t}h(t)dt`, each form factor is analytic
for `Im omega>0`, obeys the appropriate reality reflection, and satisfies a
Kramers--Kronig dispersion relation once its subtraction/contact data are
fixed.  Causality does not determine the spectral density or contact terms.
Before a physical spacetime cone is earned, (LM33) means mission-time
retardedness only; it cannot be promoted to spacetime-cone support by words.

### 6.2 Passivity

Fisher positivity is a positive information metric, not dynamical passivity.
If, additionally, `\delta L` is established as a conjugate generalized force
and the work convention is

\[
 \dot W=\langle\delta L,\partial_t\delta J\rangle_{G_D},
 \qquad
 \langle x,y\rangle_{G_D}:=\langle Dx,Dy\rangle_F, \tag{LM34}
\]

then harmonic passivity gives

\[
 \overline{\dot W}
 ={1\over2}\sum_{r=A,E,T}
 \omega g_r\,\operatorname{Im}h_r^R(\omega)
 \|P_r\delta L\|^2\ge0,                          \tag{LM35}
\]

where `(g_A,g_E,g_T)=(1/2,2,1)`.  For independently drivable sectors and
positive frequency this implies `Im h_r^R(omega)>=0` in the stated sign
convention.  It constrains dissipative signs and spectral residues, not their
magnitudes, their equality, or the absolute calibration.  A different
physical work pairing must be stated and rederived rather than imported from
(LM34).

### 6.3 Diffeomorphism Ward identity

At the fully soldered continuum endpoint, a complete diffeomorphism-invariant
effective action and on-shell nonmetric sectors give a covariant residual
identity of the EX form

\[
 \nabla_\mu\mathcal E^{\mu\nu}=0.                 \tag{LM36}
\]

Together with six complete spatial variations and prospectively certified
zero initial constraints, this propagates the four normal constraints and
closes the ten metric equations.  It does not by itself identify `D` as a
physical metric, determine `(h_A,h_E,h_T)`, impose (LM26), select an
Einstein--Hilbert principal operator, or normalize its coefficient.  In an
Einstein endpoint the gauge and constraint packet can remove an independently
propagating trace scalar and leave the healthy spin-two pole, but that is
additional response physics, not a consequence of the six-dimensional
kinematic decomposition alone.

## 7. Absolute-calibration ceiling

The information-to-length convention must remain explicit.  In the GSGB/EW
QFI convention,

\[
 \boxed{s=\ell_F^2F.}                              \tag{LM37}
\]

In the ET Bures-line convention,

\[
 \boxed{q={\ell_B^2\over4}Q,\qquad \ell_B=2\ell_F} \tag{LM38}
\]

when `F` and `Q` denote the same Fisher/QFI form.  Therefore the physical
version of the EW Jacobian is (LM22), not `D` alone.  Setting `ell_F=1` would
erase an unproved physical calibration.

The following are exact, dimensionless outputs of the q4 construction:

1. the `1+2+3` sector decomposition;
2. the invertibility of `D`;
3. the Gram ratios `1/2:2:1`; and
4. the three- or, with added isotropy, two-form-factor count.

The construction does **not** fix:

1. `ell_F` or an absolute spatial ruler;
2. the physical normalization of the lineage intervention `delta L`;
3. the amplitudes, poles, contacts, or spatial dispersion of `H^R`;
4. the proper-time normalization or common transport cone;
5. the Einstein--Hilbert stiffness or the stress normalization needed to
   read off numerical `G`.

Indeed `ell_F -> c ell_F` changes the physical metric scale in (LM22) while
leaving the dimensionless EW family, gamma, `D`, and its representation
theory unchanged.  A source-coordinate recalibration can likewise be
absorbed into `H^R`.  Numerical `G` therefore requires the independent
RGRL/GSGB physical join plus the complete response/stress packet; `ell_F` is
not `G`.

## 8. Exact scientific advance and remaining target

The result reduces the unknown local gravity-formation response from an
arbitrary six-by-six tensor kernel to

\[
 \boxed{
 \text{three }S4\text{ response functions}
 \quad\text{or}\quad
 \text{trace plus shear under established }O(3)\text{ isotropy}.}   \tag{LM39}
\]

It also proves why accumulated scalar distinguishability alone cannot finish
the gravity calculation: gamma supplies a positive record seed and local
information norm, while the response needs sector direction, retardation,
physical soldering, and calibrated amplitude.

The next bounded physics target is consequently not more pair-memory
machinery.  It is to identify or measure the retarded lineage-response
functions in (LM23), test the isotropy relation (LM26), classify the leading
nonzero-momentum terms, and join them to the independently controlled EX/EH
response packet.

## 9. Dependency and custody ledger

This theorem uses, without changing, the following canonical sources:

- `LANE_CROSS_RFT_GRA_EW_Q4_PAIR_MEMORY_METRIC_DEFORMATION_CLOSURE_V001/THEOREM.md`,
  SHA-256 `495e4e99171f4e3e5809f24e5a9a5b68116e996f4f29115bd22f780127d4714e`;
- `LANE_CROSS_RFT_GRA_EX_SPATIAL_DEFORMATION_CONSTRAINT_PROPAGATION_V001/THEOREM.md`,
  SHA-256 `59487e2ce0585e291ecb032215ed3a9d23883e418df2994692448ced4cf5a1f2`;
- `LANE_CROSS_RFT_MGFT_INDUCED_EH_BACKREACTION_V001/THEOREM.md`,
  SHA-256 `e73c077d1402ebab9a6061b060b9b24e9ac323b25809cccc706f4365da7a5e2f`;
- `LANE_CROSS_RFT_GRA_ES_GAMMA_SOLDERING_GRAVITY_BRIDGE_V001/THEOREM.md`,
  SHA-256 `8721183a72f4b864d06f79f6e68405a5393e37d1a6f4d24d5f4c3c2b79a81075`;
- `ET_ES_GAMMA_LENGTH_CONVENTION_CROSSWALK_V001.md`, SHA-256
  `e3037d5fcc0b449b8c46414de94365075e372f8da8ca1b6a87df85c7ef85359b`;
- `GRAVITY_RGRL_POST_ADOPTION_STRUCTURAL_THEOREM_V001.md`, SHA-256
  `733b18ecaa29c7acd755db6947b790a9ae37240a3c74d199752d5e278280783d`.

## 10. Disposition

`Q4_PAIR_MEMORY_KINEMATICS_DECOMPOSES_EXACTLY_AS_ONE_TRACE_PLUS_TWO_AND_THREE_TETRAHEDRAL_SHEARS__THE_POINT_LOCAL_S4_LINEAGE_RESPONSE_HAS_EXACTLY_THREE_FORM_FACTORS_AND_FULL_O3_ISOTROPY_REDUCES_THEM_TO_TRACE_PLUS_SHEAR__SCALAR_GAMMA_AND_RECORD_COUNT_PROVABLY_DO_NOT_DETERMINE_THOSE_FUNCTIONS_OR_ABSOLUTE_SCALE__PHYSICAL_SOLDERING_DISPERSION_RESPONSE_CALIBRATION_AND_NUMERICAL_G_REMAIN_SEPARATE_TARGETS`
