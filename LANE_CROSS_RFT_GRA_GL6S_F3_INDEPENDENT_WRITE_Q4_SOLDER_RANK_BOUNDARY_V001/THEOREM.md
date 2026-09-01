# F3 independent-write q4 solder-rank boundary theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6S_F3_INDEPENDENT_WRITE_Q4_SOLDER_RANK_BOUNDARY_V001`  
**Short name:** `GL6S V001`  
**Date:** 2026-08-31  
**Status:** author packet frozen pending independent post-freeze custody audit  
**Claim class:** exact four-write retained-pair moment law; exact global
four-parameter image ceiling; exact `S4` tangent decomposition and missing
`E`-shear theorem

**Not claimed:** a no-go for F3 outside the frozen independent-write slice; a
no-go for sources which physically change writer pair correlations; a physical
q4 solder, common cone, metric, Ricci response, gravity, or `G`.

## 1. The shortcut being tested

The direct repetition shortcut is to place four local F3 writer-retained pairs
at the q4 labels `a=1,...,4`, drive the four existing BS10/BS08 blank-target
copy/write operations independently, and hope that their six retained pair
correlations already furnish the six EW deformation coordinates.

This is the strongest amplitude-only version of that shortcut.  Four local
writers with perfectly correlated supplied contents and one common amplitude
are a more restricted special case.  A literal writer coupled across four
separated vertices is not attributed to local BS10.  The theorem below keeps
the complete four-writer content law arbitrary but fixed and asks only what
the four physical copy amplitudes can change.

## 2. Frozen four-write mission and exact pair moments

Let the four local writer contents be `s_a in {-1,+1}` with one prospectively
fixed diagonal joint law `pi(s)`.  Here `+1,-1` relabel the two F3 content
states, while the third retained state remains the unique blank.  The writer
law may contain arbitrary correlations, but it is independent of the four
later copy envelopes in this theorem.  Define

\[
 C_{ab}:=\mathbb E_\pi[s_as_b],\qquad a<b.          \tag{S01}
\]

Each retained qutrit starts blank and is driven by its local BS10 instance

\[
 U_a=\exp[-i\Phi_aK_{w_a\to r_a}],\qquad
 \Phi_a={1\over\hbar}\int_{I_a}j_a(t)dt,qquad
 p_a:=\sin^2\Phi_a.                                \tag{S02}
\]

All noncommuting terms are off or echoed on the same differentiable
open-neighborhood terms-off slice as GL6Q.  Define the retained content
observable

\[
 Z_a^r:=|+1\rangle\!\langle+1|_{r_a}
       -|-1\rangle\!\langle-1|_{r_a},              \tag{S03}
\]

which annihilates the blank.  Because the four local copy operations act on
disjoint writer-target factors, their sector weights factor at fixed writer
content.  The exact later retained pair moment is

\[
 \boxed{
 m_{ab}:=\langle Z_a^rZ_b^r\rangle
 =C_{ab}p_ap_b.}                                   \tag{S04}
\]

Equation (S04) is a coherent block-diagonal expectation, not an assumed
classical collapse.  It says nothing about whether the supplied writer
correlation is itself a record lineage.

## 3. Theorem `GL6S-1` -- exact global and tangent rank ceiling

For fixed `C=(C_ab)`, the complete amplitude-only image is the map

\[
 f_C:[0,1]^4\longrightarrow\mathbb R^6,qquad
 [f_C(p)]_{ab}=C_{ab}p_ap_b.                       \tag{S05}
\]

Its differential is

\[
 [D_pf_C(\delta p)]_{ab}
 =C_{ab}(p_b\delta p_a+p_a\delta p_b).             \tag{S06}
\]

Therefore

\[
 \boxed{\operatorname{rank}D_pf_C\le4}             \tag{S07}
\]

at every point.  The map (S05) is polynomial and hence Lipschitz on the compact
cube; its image has Hausdorff dimension at most four and therefore has empty
interior in the six-dimensional pair-coordinate space.  Equivalently, Sard's
theorem applies because every point has rank below six.
Time-shaped variations inside a terms-off write slot do not evade (S07),
because each envelope enters only through the one integrated coordinate
`Phi_a` in (S02).

### Proof

Equation (S06) is a linear map from a four-dimensional domain, proving (S07).
The polynomial/Lipschitz dimension bound above upgrades the local rank result
to the global empty-interior statement.  Allowing all four amplitudes is
already more general than the common-amplitude local-writer special case.
QED.

## 4. Theorem `GL6S-2` -- exact missing `E` shear at the `S4` point

At an `S4`-fixed background put

\[
 0<p_a=p<1,\qquad C_{ab}=C.                         \tag{S08}
\]

Let `A:R^4 -> R^E4` be the unsigned vertex-edge incidence map

\[
 (Au)_{ab}=u_a+u_b.                                \tag{S09}
\]

Then

\[
 D_pf_C=CpA,\qquad A^{\mathsf T}A=2I_4+\mathbb J_4. \tag{S10}
\]

Thus `A` is injective, with squared singular values `6` on the uniform mode
and `2` on the three-dimensional vertex-contrast mode.  Equivariance gives

\[
 \boxed{
 \operatorname{im}A=A_1\oplus T_2,\qquad
 (\operatorname{im}A)^\perp=E_2.}                 \tag{S11}
\]

Explicitly, the missing `E_2` directions obey

\[
 x_{12}=x_{34}=q_1,\qquad
 x_{13}=x_{24}=q_2,\qquad
 x_{14}=x_{23}=q_3,\qquad q_1+q_2+q_3=0.         \tag{S12}
\]

Every vertex is incident to one edge from each opposite pair, so (S12) is
orthogonal to every vector `(u_a+u_b)`.  This proves (S11) directly.

For the physical Jacobian `CpA`, the exact branch statement is

\[
 C\ne0\Longrightarrow\operatorname{im}D_pf_C=A_1\oplus T_2,
 \qquad
 C=0\Longrightarrow\operatorname{im}D_pf_C=\{0\}. \tag{S12a}
\]

Universally, `im D_p f_C` is a subset of `A1+T2`, so the `E_2` sector is
absent on both branches.

The exact EW/LMRK deformation map carries `E_2` to the two tetrahedral
diagonal-shear directions of `Sym^2(V)`.  Hence independent local write
amplitudes cannot close the complete six-direction q4 deformation tangent.
At EW's uniform independent symmetric point `(theta,J)=(0,0)`, one has
`C_ab=0`, and the copy-amplitude derivative (S06) vanishes altogether:
changing only the copy probabilities cannot create pair content correlations
which are absent from the writer law.  This does not null EW's different
derivative `D_J F_theta(0)`, which remains the exact six-dimensional
isomorphism EW20--EW21.

## 5. Functional-source form and relation to GL6R

At a common write angle `Phi`,

\[
 {\delta p_a\over\delta j_c(t')}
 ={\sin(2\Phi)\over\hbar}\delta_{ac}\mathbf1_{I_c}(t'). \tag{S13}
\]

Multiplication of (S10) by this nonzero scalar does not change its rank or
representation content.  GL6R's downstream nonlinear retained-occupation
coefficient shows that two simultaneous writes can change one local causal
block.  Its two envelope derivatives are scalar multiples of the same block,
consistent with rather than evading (S07).

This theorem does not constrain a parent in which the physical source changes
the joint writer law itself, `delta C_ab != 0`, or in which an autonomous
interaction generates pair correlations after the write.  Those are exactly
the missing dynamics.

## 6. Smallest missing response datum on the current EW/q4 route

The necessary datum on this route is a same-parent physical response with a
nonzero `E_2` projection.  One sufficient canonical realization is a
pair-correlation source or autonomous pair-history response

\[
 R^R_{ef}={\delta J_e\over\delta j_f^{\rm pair}}
 \quad\hbox{or}\quad
 R^R_{ef}={\delta m_e\over\delta j_f^{\rm pair}},  \tag{S14}
\]

with lifecycle/FORMATION-SHAM/BREAK custody, complete CTP contacts, and a
nonzero `E_2` sector.  At an `S4`-fixed point, LMRK reduces a square pair
response to

\[
 h_A=a+4b+c,\qquad h_E=a-2b+c,\qquad h_T=a-c.     \tag{S15}
\]

Full q4 tangent rank requires all three factors to be nonzero; in particular
`h_E != 0` is the direction independent writes cannot supply.  The missing
projection need not be driven by a source literally named `j_pair`: it may
come from source-dependent writer correlations, a one-body source propagated
through an owned interaction, an owned cross-target operation, or autonomous
post-write dynamics.  Existing F3, URFT, gamma/QFI accumulation, and
E-EMERGENTSPACE do not by themselves assign any of those completed response
laws.  The response must be derived from a complete existing same-parent
interaction; only after that search fails would a constitutive theory decision
be required.

No incidence coordinate is renamed `J`, no identity solder is inserted, and
no material modulus, Ricci preterm, observed `G`, or graviton appears here.

`PASS__EXACT_FOUR_LOCAL_BS08_WRITE_PAIR_MOMENTS__GLOBAL_AMPLITUDE_ONLY_IMAGE_DIMENSION_AT_MOST_FOUR__S4_NONZERO_C_IMAGE_A1_PLUS_T2__ZERO_C_IMAGE_ZERO__E2_DIAGONAL_SHEAR_MISSING__INDEPENDENT_EW_POINT_RESPONSE_ZERO__GL6R_COLLINEARITY_CONSISTENT__PAIR_CORRELATION_RESPONSE_SHARPLY_IDENTIFIED_AS_NEXT_LAW__NO_Q4_SOLDER_METRIC_RICCI_GRAVITY_OR_G_PROMOTION`
