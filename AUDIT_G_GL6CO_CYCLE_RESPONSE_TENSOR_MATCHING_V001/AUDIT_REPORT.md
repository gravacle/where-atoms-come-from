# Independent hostile audit — GL6CO cycle-response tensor matching

**Target:** `LANE_CROSS_RFT_GRA_GL6CO_CYCLE_RESPONSE_TENSOR_MATCHING_V001`  
**Disposition:** **PASS_AFTER_AUTHOR_REPAIR**  
**Target edits by auditor:** none

## 1. Verdict

The repaired GL6CO theorem is correct on its declared symmetry-classification
and matching surface.  An independent standard-library rational replay,
which neither imports nor executes the author derivation, contains
`277/277` mathematical checks; the frozen-result replay adds its byte check
and passes `278/278`.  It reconstructs the full `S4` invariant space, the
centered `GL6CL` writer jet, the exact tensor solder, the necessary-and-
sufficient `SO(3)` extension plane, the stronger full-reference diagnostic,
and the `GL6BV` contact including its sign and common-sublattice
normalization.

Two defects were found and reported before this verdict:

1. a packet guard had been weakened to the unsigned substring `ordinary
   three-vector`, so it would also have accepted the opposite physical
   claim; and
2. the conditional contact-plus-cycle equations mixed the unnormalized
   `GL6CL` common source coordinate with the normalized `GL6BV` common
   projection, missing a relative factor of two.

The author repaired the first guard to require the explicit negation and
repaired the second issue by using the normalized common source throughout:
the cycle block is now `mu^2 B_+^* K B_+/2`, while the displayed contact
block remains `g_ct(8+Re C)`.  The core cycle-only coefficients and matching
conditions did not change.

This is a matching theorem, not the calculation of the unknown stationary
coefficients.  It does not establish the assumed infinite stationary symbol,
perform the response-to-1PI inversion, complete the `E2` blocks, derive a
metric/Ricci operator, prove gravity, or calculate `G`.

## 2. Independent method

The independent replay starts only from the declared tetrahedral vectors,
the centered six-cycle incidence, the microscopic pair-to-tensor solder,
and the source-before-Feshbach defect vectors.  It uses exact fractions and
degree-two polynomial dictionaries.  It does not import the GL6CO module or
its ledger.

The replay:

- enumerates all 24 port permutations and constructs their signed-permutation
  coordinate actions;
- evaluates the relevant representation characters and separately checks
  covariance and linear independence of the five displayed quadratic
  covariants;
- reconstructs all twelve centered cycle/pair offsets and composes the full
  four-by-three writer jet with each cycle-kernel basis term;
- derives the symmetric-tensor restriction directly, rather than treating
  the three `T2` coordinates as an ordinary vector;
- constructs the full six-component pair solder and the generic-momentum
  linearized-Einstein reference matrix independently; and
- rebuilds the four defect frames and their common contact expansion,
  including the physical common/relative coordinate transformation.

## 3. `S4` classification — PASS

For port permutation `g`, the replay obtains the exact coordinate action

\[
 R_g={1\over4}\sum_{d=0}^3T_{g(d)}T_d^T,
\]

and verifies that it is an orthogonal signed permutation and maps
`T_d` to `T_{g(d)}`.  With

\[
 u={1\over2}(1,1,1,1)^T,\qquad Q_{di}={1\over2}(T_d)_i,
 \qquad S=(u\;Q),
\]

the audit finds `S^T S=I` and `S^T P_g S=1+R_g` for every group element.

The exact character average is

\[
 {1\over24}\sum_g
 \chi_{\operatorname{Sym}^2(T_2)}(g)
 \chi_{\operatorname{Sym}^2(A_1+T_2)}(g)=5.            \tag{A-CO01}
\]

The constant symmetric invariant space has dimension two.  Its projectors
are `uu^T` and `QQ^T`; the stationary common-amplitude null kills the former,
leaving only `kappa QQ^T`.

The five quadratic matrices

\[
 \alpha r^2|A\rangle\langle A|,\quad
 \eta(Aq^T+qA^T),\quad
 b r^2I_{T},\quad cD,\quad dO                       \tag{A-CO02}
\]

are each covariant under all 24 elements and are exactly linearly
independent.  Dimension (A-CO01) therefore proves completeness: there is no
sixth quadratic invariant hidden by the displayed parametrization.

## 4. Independent writer pullback — PASS

From the centered incidence, the audit reconstructs twelve offsets
`rho_(d,p)`, each with `|rho|^2=11/4`, and verifies

\[
 \sum_{d,p}\rho_{d,p}\rho_{d,p}^T=11I.
\]

In the raw pair basis `t_i`, the common writer has

\[
 B_{\rm raw}(0)=-2T,
\]

where the four rows of `T` are the tetrahedral vectors.  Passing to the
orthonormal pair basis `t_i/sqrt(2)` supplies the required factor
`1/sqrt(2)`.  Composing the full quadratic cosine jet with the bare cycle
kernel gives, independently,

\[
 {1\over\mu^2}{\cal H}^{H6}_T
 =8\kappa I+(-2\kappa+8b)r^2I
 +(-16\kappa+8c)D+(12\kappa+8d)O+O(r^4).              \tag{A-CO03}
\]

The audit obtains the six parameter contributions separately:

```text
kappa : 8 I - 2 r2 I - 16 D + 12 O
b     : 8 r2 I
c     : 8 D
d     : 8 O
alpha : 0 through k2
eta   : 0 through k2
```

Thus `alpha` and `eta` do not enter by accidental cancellation.  The
zero-mode writer lies wholly in cycle `T2`; inserting its first cycle-`A1`
component costs another two powers of momentum.

## 5. Tensor solder and exact extension test — PASS

The pair-to-symmetric-tensor map is independently reconstructed as

\[
 {cal S}(j)={1\over4}\sum_{a<b}j_{ab}
 (T_aT_b^T+T_bT_a^T).                                  \tag{A-CO04}
\]

It sends the pair trace `A` to `-I` and sends the three raw tensor vectors
to minus the `yz`, `zx`, and `xy` symmetric off-diagonal matrices.  Hence in
the orthonormal pair basis

\[
 (h_{yz},h_{zx},h_{xy})=-{1\over\sqrt2}(x_1,x_2,x_3). \tag{A-CO05}
\]

These coordinates are not an ordinary three-vector.  Restricting a general
parity-even `SO(3)`-covariant quadratic operator on a symmetric tensor to
off-diagonal input and output leaves exactly

\[
 u r^2I+v[(r^2I-D)+O].                                  \tag{A-CO06}
\]

Therefore a cubic form `A r^2I+B D+C O` extends if and only if

\[
 \boxed{B+C=0}.                                         \tag{A-CO07}
\]

Using (A-CO03), this is exactly

\[
 \boxed{c+d={\kappa\over2}}.                           \tag{A-CO08}
\]

It is one condition, not ordinary-vector isotropy.  The positive analytic
witness `kappa=alpha=1`, `eta=0`, `b=c=d=1/4` passes and pulls back to
`-14(D-O)` at quadratic order.

The stronger held-out Fierz--Pauli/linearized-Einstein `T2-T2` reference is
`D-O`.  Proportionality additionally requires

\[
 \boxed{b={\kappa\over4}}.                              \tag{A-CO09}
\]

At `k=(2,3,5)`, the independently soldered full six-by-six reference matrix
is exactly

\[
\begin{pmatrix}
-38&5&37&15&10&6\\
5&100&20&-30&20&0\\
37&20&4&-30&-20&24\\
15&-30&-30&4&-6&-10\\
10&20&-20&-6&9&-15\\
6&0&24&-10&-15&25
\end{pmatrix}.                                         \tag{A-CO10}
\]

It has rank three, its `T2-T2` block is `D-O`, and its `A1/E2-T2` blocks are
nonzero.  This confirms the target's kill test: passing (A-CO09) in a
three-dimensional subblock is not a Ricci proof.

## 6. `GL6BV` contact and sign — PASS AFTER NORMALIZATION REPAIR

The four one-defect pair words obey

\[
 \tau_a^T\tau_b=8\delta_{ab}-2,
 \qquad\sum_a\tau_a=0,
 \qquad\sum_a\tau_a\tau_a^T=8P_T.                     \tag{A-CO11}
\]

With `Q_a=tau_a tau_a^T/6`, `theta_a=k.T_a`, and
`p=Pr[sigma_v(a)=sigma_w(a)]`, the replay derives

\[
 \sum_a\theta_a^2=4r^2,
 \qquad\sum_a\theta_a^2Q_a={4\over3}r^2I+{8\over3}O. \tag{A-CO12}
\]

It then obtains the normalized common contact block

\[
 G_{+,T}^{\rm ct}={16\over3}(1+2p)I
 +{4\over3}(1-4p)r^2I
 +{8\over3}(2p-1)O+O(r^4),                            \tag{A-CO13}
\]

with no `D` term.  Its mismatch is therefore

\[
 \Delta_{\rm ct}={8\over3}(2p-1).                     \tag{A-CO14}
\]

At the bounded `Q4` witness `p=109/128`, this is `15/8`.  The fixed-branch
energy Hessian has sign `-g_ct G`, while the connected functional has sign
`+g_ct G`, where `g_ct=h^2/(4U_d^3)`.  The target's repaired sign is correct.

The hostile normalization test is essential.  `GL6CL` defines

\[
 j_P=j_++j_-,\qquad j_C=j_+-j_-,\qquad B_+=B_P+B_C.
\]

Thus `B_+` differentiates along the unnormalized physical vector `(j_+,j_+)`.
The normalized common vector is `(j_+,j_+)/sqrt(2)` and its writer is
`B_+/sqrt(2)`.  In the normalized convention consistently adopted by the
repaired target,

\[
 {\cal H}^{H6}_{+,T}={\mu^2\over2}B_+^*K_{\rm cyc}B_+,
 \qquad {\cal H}^{\rm ct}_{+,T}=g_{\rm ct}G_{+,T}^{\rm ct}. \tag{A-CO15}
\]

Equivalently, one may retain the unnormalized cycle block `mu^2 B_+^*KB_+`
only if every contact coefficient is doubled.  The first target snapshot
mixed these conventions.  The repaired conditional tests are

\[
 {\mu^2\over2}[-4\kappa+8(c+d)]
 +{h^2\over4U_d^3}{8\over3}(2p-1)=0,                  \tag{A-CO16}
\]

and

\[
 {\mu^2\over2}[-2\kappa+8b]
 +{h^2\over4U_d^3}{4\over3}(1-4p)=0.                 \tag{A-CO17}
\]

No match is asserted by displaying these equations.

## 7. State, 1PI, and scaling guards — PASS

The contact parameter `p` and the cycle coefficients
`kappa,b,c,d` must be evaluated in one stationary state and one complete
source-first functional before (A-CO16)--(A-CO17) may be used.  The target
does not substitute the nonstationary `Q4` orbit witness into that total.

The classified object is a response/susceptibility Hessian.  It is not an
ordinary 1PI quadratic kernel.  The full `E2-T2` and `E2-E2` blocks, contact,
same-state completion, and lawful response-to-1PI inversion remain mandatory
before a Ricci comparison.

Units and strong-lock power counting also pass:

\[
 \mu={105\over8}{h^6\over U_d^6},\qquad
 J={63\over8}{h^6\over U_d^5},\qquad
 {\mu^2\over J}=O(h^6/U_d^7),\qquad
 g_{\rm ct}=O(h^2/U_d^3).                              \tag{A-CO18}
\]

The factor `1/2` in the normalized common projection does not change this
order separation.  Isolated strong-lock repetition therefore cannot
generically cancel a nonzero leading contact mismatch term by term.  The
target correctly labels the alternatives—finite-ratio/collective response,
an independently vanishing leading relation, or another same-order block—as
power counting, not as a proved phase transition.

## 8. Claim boundaries and contemporaneous GL6CN result

The target assumes rather than derives an analytic infinite stationary
cycle symbol.  Symmetry classifies `kappa,b,c,d`; it does not calculate them.
Physical calibration of the inherited character coordinate is also absent.

The target froze while the diagonal sixth-order first-source term was still
assigned to `GL6CN`.  That separate author packet has now proved the term is
pointwise zero.  This removes a possible same-order diagonal first vertex
but does not calculate GL6CO's stationary coefficients, create the missing
`E2` response, perform the 1PI inversion, or promote the matching test to
gravity.

No source-second contact is double-counted as a spectral two-first-vertex
term.  No stationary bulk state, massless pole, physical cone, metric,
Ricci/Einstein law, gravity, graviton, or value of `G` is claimed.

## 9. Custody

The repaired author packet passes `217/217` science checks and `175/175`
packet checks.  `TARGET.sha256` pins all twelve repaired target bytes.  The
independent mathematical result contains `277/277` checks and its frozen
replay passes `278/278`; `INDEPENDENT_RESULT.json` freezes the exact output.
This audit has its own manifest and one-line seal.

**Final disposition: PASS_AFTER_AUTHOR_REPAIR.**

`PASS__GL6CO_INDEPENDENT_HOSTILE_AUDIT__S4_DIMENSION_5__CL_PULLBACK_EXACT__PAIR_TENSOR_SOLDER_EXACT__SO3_EXTENSION_IFF_B_PLUS_C_ZERO__REFERENCE_ADDS_A_ZERO__FULL_SOLDER_RANK3_CROSS_BLOCK_GUARD__BV_CONTACT_AND_SIGNS_EXACT__COMMON_SUBLATTICE_FACTOR_REPAIRED__SAME_STATE_AND_1PI_GUARDS__MATCHING_ONLY__NO_RICCI_GRAVITY_OR_G`
