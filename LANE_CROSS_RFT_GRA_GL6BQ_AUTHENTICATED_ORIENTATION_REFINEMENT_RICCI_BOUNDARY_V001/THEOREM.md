# Authenticated orientation-refinement Ricci boundary

**Lane:** `LANE_CROSS_RFT_GRA_GL6BQ_AUTHENTICATED_ORIENTATION_REFINEMENT_RICCI_BOUNDARY_V001`  
**Short name:** `GL6BQ V001`  
**Date:** 2026-09-02  
**Status:** exact narrow author result; exact replay passed  
**Claim class:** exact joint-frame orbit average of the audited `GL6BO`
four-cell response; exact tetrahedral-versus-Haar moment separation; exact
metric-normalized `E2/T2` condition; conditional rotational-design theorem;
exact parent-ownership boundary

**Not claimed:** that the F3 parent prepares a rotational ensemble; that an
`SO(3)` average may be inserted as a new law; an all-cluster convergence
theorem; Ward/Bianchi closure; a common cone; a full Ricci operator; gravity;
a graviton; or `G`.

## 1. Result first

`GL6AA` authenticates complete four-port frames whose transitions lie in
`S4`. In the tetrahedral realization these give exactly 24 orthogonal frame
maps

\[
 \Gamma=\{R_\sigma:\sigma\in S_4\},\qquad
 R_\sigma n_a=n_{\sigma(a)}.                         \tag{BQ01}
\]

The `A3` translations of `GL6AK` change the anchor but not this orientation
set. Every product of authenticated frame transitions remains in the same
finite group \(\Gamma\).

This exhausts the complete four-cell star embeddings, not merely their
labels. For any interior shared child \(c\), its four parents are
\(m_a=c-e_a\). Under the inherited affine tetrahedral chart their positions
relative to their common centroid are \(-a_*\nu_a\). Changing \(c\) is only
an \(A_3\) translation; changing complete port frames applies one
\(R_\sigma\); the dual incoming/outgoing choice applies central inversion.
The latter is invisible to the parity-even \(k^2h^2\) coefficient. Thus all
authenticated embeddings of this complete motif occupy one orientation
coset. Boundary fragments are smaller clusters, not additional embeddings
of the completed four-cell response.

Let \(F_p(k,h)\) be the genuine-four-cell `GL6BO` TT quadratic form at one
of the retained coefficients \(\tau^p t^3 k^2\), \(p=4,6,\ldots,16\), after
the causal Schur quotient and Boolean connected subtraction. The exact
replay proves

\[
 F_p(R_\sigma k,R_\sigma hR_\sigma^T)=F_p(k,h)
 \quad(\sigma\in S_4).                              \tag{BQ02}
\]

Consequently the average over **all orientations authenticated by the
present atlas** is the identity operation:

\[
 \boxed{{1\over24}\sum_{\sigma\in S_4}
 F_p(R_\sigma k,R_\sigma hR_\sigma^T)=F_p(k,h).}     \tag{BQ03}
\]

It is a covariance orbit of one crystallographic tetrahedral frame, not an
ensemble of 24 orientations in \(SO(3)/\Gamma\). In particular, at the
leading coefficient the authenticated-orbit six-TT row remains

\[
 a_*^2\left(0,{5411224\over243675},{87831613\over23880150},
 {74936332\over11940075},{130214356\over21492135},
 {130214356\over21492135}\right),                  \tag{BQ04}
\]

with the unchanged nonzero Ricci-ray residual norm

\[
 \boxed{\|\rho\|^2=
 {81124093179109711373\over277147120114935000}a_*^4>0.} \tag{BQ05}
\]

Thus the physically available finite frame average does not repair the
four-cell non-Ricci result.

## 2. Tetrahedral/cubic invariance is not rotational refinement

The 24 full-tetrahedral maps include improper maps. Because the scored
response is parity even (\(k^2h^2\)), multiplying those maps by central
inversion gives the equivalent proper-cubic action. This changes neither
the response nor the following conclusion: the point group is finite and is
not rotationally isotropic at the degree needed here.

For one coordinate axis the exact orbit moments are

\[
 {1\over24}\sum_{R\in\Gamma}(e_3^TRe_3)^4={1\over3}\ne{1\over5},
 \qquad
 {1\over24}\sum_{R\in\Gamma}(e_3^TRe_3)^6={1\over3}\ne{1\over7}. \tag{BQ06}
\]

The right sides \(1/5\) and \(1/7\) are the corresponding Haar moments.
Hence the authenticated point group is not even a degree-four spherical
design, and in particular is not a degree-six rotational design.

An aligned exhaustion, a translation Følner average, or a scale change
\(a_*\to a_L\) cannot alter a dimensionless normalized residual. Both the
desired coefficient and (BQ04)'s residual scale by \(a_L^2\); after the
nontrivial \(k^2\) coefficient is normalized, their ratio is unchanged. If
they are not normalized, both vanish together. Therefore no refinement flow
presently derived in `GL6AA`, `GL6AK`, or `GL6BL` replaces a missing
orientation law.

## 3. Exact \(E_2/T_2\) equality condition

In pair order \((01,02,03,12,13,23)\), let

\[
 P_A={J_6\over6},\qquad P_T={I-O\over2},\qquad
 P_E={I+O\over2}-P_A.                              \tag{BQ07}
\]

After a legitimate \(S_4\) twirl, a reciprocal pair block has the form

\[
 H_\Gamma=h_A P_A+h_E P_E+h_TP_T.                  \tag{BQ08}
\]

For the inherited tetrahedral pair-to-\(\operatorname{Sym}^2\) solder \(D\),
the exact Gram operator is

\[
 D^*D={1\over2}P_A+2P_E+P_T.                       \tag{BQ09}
\]

It follows that spatial equality of the two traceless sectors is not the raw
condition \(h_E=h_T\). It is

\[
 \boxed{{h_E\over2}=h_T\quad\Longleftrightarrow\quad h_E=2h_T.} \tag{BQ10}
\]

Equation (BQ10) must hold separately at every independent retained
time/frequency and derivative coefficient of the **total finite-family or
converged accumulated response** to which an analytic Ricci identity is
claimed. It need not hold cluster by cluster: connected-cluster `E/T`
anisotropies may cancel in the owner-correct accumulated density. Apply
(BQ10) to an individual cluster only when making an explicitly clusterwise
Ricci claim. The scalar \(A_1\) coefficient remains independent.

## 4. Conditional orientation-measure theorem

Suppose a future same-parent construction supplies physical cluster
orientations \(Q_C\in SO(3)/\Gamma\) and prospectively owned nonnegative
motif frequencies \(w_{C,L}\). Define their normalized empirical measure

\[
 \mu_L={1\over W_L}\sum_Cw_{C,L}\delta_{Q_C},\qquad
 W_L=\sum_Cw_{C,L}.                                  \tag{BQ11}
\]

For one fixed coefficient, let \(r(Q)\in\mathbb R^6\) be its six-TT row after
jointly rotating both the spatial support and the tensor solder. The exact
condition for that weighted response to lie on the diagnostic Ricci ray is

\[
 \boxed{\left(I_6-{\mathbf1\mathbf1^T\over6}\right)
        \int r(Q)\,d\mu_L(Q)=0.}                    \tag{BQ12}
\]

For the whole class of \(k^2h^2\) kernels, a checkable sufficient finite rule
is that \(\mu_L\) be a weighted \(SO(3)\) design through rotation degree six:

\[
 \int\prod_{s=1}^{m}Q_{i_sj_s}\,d\mu_L(Q)
 =\int_{SO(3)}\prod_{s=1}^{m}Q_{i_sj_s}\,dQ,
 \quad0\le m\le6.                                  \tag{BQ13}
\]

Moment convergence in (BQ13), together with the uniform response bounds
required by `GL6BL`, is sufficient in a refinement sequence. Dense support
alone is not: the weights must become Haar-equidistributed through degree
six. Conversely, if a rule is to isotropize every kernel in this polynomial
class, equality of these moments is required.

The exact unit-quaternion Haar replay gives, conditionally on (BQ13), the
leading genuine-four-cell common TT value

\[
 \boxed{c_{\rm Haar}^{(4)}={2927638768\over417902625}a_*^2,\qquad
 \bar r_{\rm Haar}^{(4)}=c_{\rm Haar}^{(4)}\mathbf1_6.} \tag{BQ14}
\]

For the full \(N=1\) four-cell response rather than its connected Boolean
weight, the corresponding conditional leading value is \(512a_*^2/675\).
These are consequences of a hypothetical Haar/design ensemble; they are not
coefficients of the current physical parent. A common six-TT row is also
only the held-out Ricci classifier. Full promotion still requires the
owner/contact, Ward/Bianchi, temporal, cone, and refinement gates.

## 5. Ownership verdict

The present parent does not supply (BQ11)--(BQ13):

1. `GL6AA` authenticates selected IDs, incidences, and \(S_4\) port frames,
   but explicitly does not autonomously select its program, support, or
   frames.
2. `GL6AK` proves \(A_3\) translation and finite \(S_4\) covariance. Its
   averaged stationary state is an existence construction, not a unique or
   dynamically selected orientation distribution, and it contains no
   \(SO(3)/\Gamma\) orientation variable.
3. `GL6BL` gives conditional refinement/convergence criteria; it does not
   generate new rotated embeddings or their frequencies.

Therefore assigning equal weights to rotated copies outside \(\Gamma\), a
Haar measure, random crystallite orientations, or a degree-six design at
this stage would be an inserted state/preparation law. The narrow conclusion
is

\[
 \boxed{\text{current accumulated parent: finite }S_4\text{ covariance only;}
 \quad\text{SO(3)-isotropic/Ricci TT averaging is not required.}} \tag{BQ15}
\]

`PASS__ALL_GL6AA_ORIENTATIONS_ENUMERATED__JOINT_S4_ORBIT_AVERAGE_EQUALS_GL6BO_KERNEL__NONRICCI_RESIDUAL_SURVIVES__FINITE_TETRAHEDRAL_CUBIC_GROUP_FAILS_HAAR_MOMENTS__METRIC_NORMALIZED_ISOTROPY_HE_OVER_2_EQUALS_HT__SO3_DEGREE6_DESIGN_CONDITION_EXACT__DENSE_SUPPORT_INSUFFICIENT__NO_CURRENT_PARENT_GENERATED_OR_REQUIRED_ORIENTATION_WEIGHTS__HAAR_RESULT_CONDITIONAL_ONLY__NO_RICCI_GRAVITY_OR_G_CLAIM`
