# Independent hostile audit — GL6W normalized operator-neutral q4 CTP kernel V001

**Target:** `LANE_CROSS_RFT_GRA_GL6W_Q4_NORMALIZED_OPERATOR_NEUTRAL_CTP_KERNEL_V001/`  
**Frozen theorem SHA-256:** `c416e9fe329e30e868e0332f27f3c3ef065d7c13530400f536e6f5ad152d7aaf`  
**Frozen MANIFEST-file SHA-256:** `2282c8b52652aeee54c9a123f5438d4ce9483d861228ba1d186feb3ac373a819`  
**Disposition:** `PASS_AT_FORMED_KEEP_ACTIVE_SYSTEM_PULSE_CTP_SCOPE__FULL_APPARATUS_SCHUR_GLUE_IR_RICCI_GRAVITY_AND_G_OPEN`

## Custody and independence

The audit pins every byte named by the frozen author manifest, separately pins
the manifest file, and rechecks every antecedent named in the author dependency
ledger.  The supplied theorem and manifest-file hashes match the frozen bytes.

The independent replay does not import or execute either author verifier.  It
rebuilds the four-qubit Pauli parent, all six literal pair operators, the CTP
pulse functional, the line-graph association scheme, and the interacting total-
spin blocks directly.  Three factorized and three interacting parameter/time
witnesses are used rather than the author's single numerical witness.

## Negative half-source sign and CTP factors

The frozen GL6V source is

\[
 H_{\rm src}(t)=-\frac12\sum_AJ_A(t)Q_A,
 \qquad Q_A=E_\star M_A,
 \qquad J_A(t)=j_AT_\star f_A(t),
 \qquad E_\star=\frac{\hbar}{T_\star}.
\]

The audit independently checks both the clock-coordinate identity

\[
 -\frac12J_AQ_A=-\frac{\hbar}{2}j_Af_A M_A
\]

and the actual pulse unitary

\[
 P_A(\eta)=\exp\!\left(+\frac{i\eta Q_A}{2\hbar}\right).
\]

Therefore, with the frozen GL6T convention

\[
 \chi^R_{MM}=-\frac{i}{\hbar}\Theta(t-s)
 \langle[M_B(t),M_A(s)]\rangle,
\]

the physical full-query response has the opposite sign:

\[
 \boxed{{\cal G}^R=-\frac{E_\star^2}{2}\chi^R_{MM}.}
\]

The audit differentiates the complete two-branch two-pulse functional by
independent central differences.  It recovers

\[
 2W_{\Delta c}={\cal G}^R,
 \qquad 2W_{c\Delta}={\cal G}^A,
 \qquad W_{\Delta\Delta}=\frac{i}{4\hbar}{\cal N},
 \qquad W_{cc}=0.
\]

The equal-source trace is one for unrelated common pulse values.  This replay
directly catches the source sign, the two factors of one-half, the noise factor,
and the distinction between equal-source unitarity and retarded response.

## Closed factorized kernel

For `U_d=0`, the independent sixteen-state evolution reproduces the frozen
one-link correlator

\[
 \langle Z(t)Z(s)\rangle=u+iv
\]

and the complete pair kernels

\[
 \boxed{{\cal G}^R=-\frac{E_\star^2}{\hbar}\Theta(t-s)
 \left(2uvI_6+vz_tz_sA_{L(K_4)}\right),}
\]

\[
 \boxed{{\cal N}=E_\star^2\left[
 (u^2-v^2-z_t^2z_s^2)I_6+z_tz_s(u-z_tz_s)A_{L(K_4)}
 \right].}
\]

All thirty-six correlator entries and the `A1`, `E2`, and `T2` spectra agree at
each witness.  The exact equal-time noise eigenvalues are nonnegative.  A
separate direct double-commutator calculation and short-time response give

\[
 \left.\partial_t{\cal G}^R(t,0)\right|_{0^+}
 =-\frac{E_\star^2}{2\hbar^2}D(\tau),
\]

including the frozen physical sign and normalization.

## Interacting compiler and `S4` typing

The audit independently diagonalizes the complete interacting sixteen-state
Hamiltonian and obtains

\[
 \operatorname{spec}(H_K-E_{\rm blank})
 =\operatorname{spec}(\mathsf H_{J=2})
 \cup3\operatorname{spec}(\mathsf H_{J=1})
 \cup2\operatorname{spec}(\mathsf H_{J=0}).
\]

The multiplicities are exactly `5 + 3*3 + 2*1 = 16`.  The independently
constructed pair projectors have ranks one, two, and three and resolve the
identity.  For every interacting witness, the complex correlator, retarded
kernel, and connected noise kernel each have one common diagonal, adjacent,
and opposite entry.  Thus they have the exact association-scheme form

\[
 K=k_dI_6+k_aA_L+k_bA_O
\]

with sector eigenvalues

\[
 K_{A_1}=k_d+4k_a+k_b,\qquad
 K_{E_2}=k_d-2k_a+k_b,\qquad
 K_{T_2}=k_d-k_b.
\]

The direct interacting double commutator agrees with

\[
 D=-8hxI_6-4hyA_L,
\]

and its short-time physical response again has the frozen normalization.  The
matched BREAK Hamiltonian gives identically zero retarded and connected-noise
blocks at every interacting witness.

## Contact, apparatus, and physics scope

The theorem correctly claims only the direct affine **system** seagull as
zero.  It does not infer that controller, clock, work, switching-boundary,
reference, source-setting, or full formation-instrument blocks vanish.  The
formed-KEEP state is a normalized conditional response sector rather than an
unconditional formation average.

The scheduled GL6V source remains a pulse-insertion interface with parent
evolution off or refocused during each sandwich.  It is not promoted to
arbitrary simultaneous smooth source evolution.  Active overlap calculations
may proceed in parallel, but complete gluing or infrared promotion is
explicitly withheld until the missing apparatus ledger and causal Schur
quotient are closed.

No pair coordinate is renamed a metric component.  No common-query
localization, gluing, causal cone, Ward/Bianchi identity, infrared operator,
Ricci response, gravity, or `G` is claimed.

**Audit verdict: PASS.**

