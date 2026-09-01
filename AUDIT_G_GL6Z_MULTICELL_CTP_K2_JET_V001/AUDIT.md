# Independent hostile audit — GL6Z multicell CTP `k^2` jet V001

**Target:** `LANE_CROSS_RFT_GRA_GL6Z_MULTICELL_CTP_K2_JET_V001/`  
**Frozen theorem SHA-256:** `7d938b97c3d7f6a2c4ea53f39f092beb44bb6cd6cc842de5dbf3384b9d071a07`  
**Frozen MANIFEST-file SHA-256:** `a433bae80357b5a4b29c45920c0d5e4262b06f5ee0e2d16e3fefec1347fe0c23`  
**Frozen SEAL-file SHA-256:** `b1189a8f369212d37467c27030e461276bc31dcabe948f8b196a751a19c7fd63`  
**Disposition:** `PASS__BRANCH_NORMALIZED_ENTRANCE_TIME_MULTICELL_JET_AND_A3_PROGRAM_FRAME_SYMBOL_EXACT__RANK4_E2_NULL_AND_A1_K2_EXACT__PHYSICAL_ATLAS_RICCI_GRAVITY_AND_G_OPEN`

## Frozen custody

The audit pins every author byte named by the frozen manifest, the manifest
file, the seal file, and every author antecedent.  The explicit theorem,
manifest-file, and seal-file digests equal the supplied post-freeze custody
hashes.  The seal content pins the frozen author manifest.  The dependency
chain reaches the audited FPSS/FJ response parent, EO/FD program-frame
algebra, and independently audited GL6T through GL6Y chain without drift.

The independent replay uses only the Python standard library.  It neither
imports nor executes the author physics verifier.

A second independent post-freeze hostile reader separately reconstructed the
mixed coefficient, Fourier orientation, rank/sector result, scalar `k^2`
coefficient, conditional source-time scope, and every custody row.  Its
verdict was `CLEAN` with no material theorem or custody defect.

Custody note: the frozen `VERIFICATION.txt` transcribes the earlier packet
count `61/61`, while the frozen `verify_packet.py` now executes `81/81` on
the frozen manifest and seal.  The audit pins both bytes and reports both
numbers.  This is a stale transcript count, not a failed custody or physics
check.

## Conditional CTP and source/read scope

The finite source/read construction is exact on each fixed FPSS member.
GL6V's six commuting pair sources tensor over disjoint parent-cell link
factors, use independent ancillas and one common registered source clock,
and uncompute every ancilla.  The complete read is the computational-basis
read of all `4a_N` links; the `6a_N` pair values are deterministic functions
of that complete outcome rather than a purported complete pair-only read.

The frozen theorem correctly defines

\[
 {cal G}^{R,N}_{nB,mA}(t,s)
 ={iE_\star^2\over2\hbar}\Theta(t-s)
 \langle[M_{n,B}(t),M_{m,A}(s)]\rangle_{\rho_{\tau,N}}
\]

as a branch-normalized conditional active-system kernel on the exact
all-formation/KEEP branch.  Retaining the complete `F/S` instrument supplies
custody but does not turn this normalized branch into an unconditional
full-instrument response.  A selected-apparatus lift remains conditional on
the exact GL6X factorization premises with the whole interacting link block
inside the system.

Sections 3--5 explicitly specialize to source time `s=0`.  This is
load-bearing because the prewait state is not stationary.  For general
source time,

\[
 \langle[M_B(t),M_A(s)]\rangle_{\rho_\tau}
 =\langle[M_B(t-s),M_A(0)]\rangle_{\rho_{\tau+s}},
\]

so the entrance moment shifts from `tau` to `tau+s`.  The theorem does not
silently assume time-translation invariance.

## Complete first cross-cell jet

For an adjacent cell move `n=m+e_a-e_b`, the audit independently obtains

\[
 [M_{n,B},M_{m,A}]
 =[{\rm ad}_HM_{n,B},M_{m,A}]
 =[{\rm ad}_H^2M_{n,B},M_{m,A}]=0,
\]

\[
 \boxed{
 [{\rm ad}_H^3M_{n,B},M_{m,A}]
 =-8h^2U_dP_{B,b}P_{A,a}
 X_{m,a}X_{n,b}Z_{m,A\setminus a}Z_{n,B\setminus b}.}
\]

The unique minimal route contains one transverse flip on each half-port and
the inherited shared-child `2U_dn_en_f` connector.  One-link fields and
external incident links cannot enter the cubic operator jet without
exceeding its insertion count.

A separate sparse computational-basis Krylov replay tests the fourth
blank-prewait moment on the core path and on two distinct externally dressed
interaction graphs.  Writing

\[
 \delta=\Delta+2(1-2d_\star)U_d,
\]

all external dependence cancels and the exact result is

\[
 8h^2(\delta+U_d)(3\delta+2U_d)
 =8h^2A_\star B_\star.
\]

The lower blank moments vanish.  Combining the fourth-moment factor `1/4!`,
the Heisenberg cubic factor `i^3/3!`, the GL6T diagnostic sign, and GL6W's
physical half-source normalization gives

\[
 \boxed{
 [\tau^4t^3]{\cal G}^{R,N}_{nB,mA}(t,0)
 =-{2E_\star^2h^4U_dA_\star B_\star\over9\hbar^8}
 P_{B,b}P_{A,a}.}
\]

The sign, factor `2/9`, eight powers of `hbar`, and physical response units
are exact.  Matched BREAK has no transverse connector-half-port flips and is
response-silent.  The outer `U_d`, not the entrance moment by itself, owns
the cross-cell effect.

## Program-frame Fourier symbol

The audit reconstructs the six-by-four pair/port incidence matrix `P` and
finds

\[
 P^TP=2I_4+\mathbf1\mathbf1^T,
 \qquad \operatorname{rank}P=4.
\]

For the inherited coordinate

\[
 X(m)=a_\star\sum_am_an_a,
 \qquad z_a=e^{-ia_\star k\cdot n_a},
\]

the frozen positive-exponent convention

\[
 \widetilde F(k)=\sum_me^{+ik\cdot X(m)}F_m
\]

assigns the move `e_a-e_b` the Laurent phase `z_b/z_a`.  Entry-by-entry
reconstruction then gives precisely

\[
 \boxed{{\cal K}(k)=P[zz^\dagger-I_4]P^T.}
\]

With the negative-exponent convention the symbol is its complex
conjugate/transpose.  The theorem correctly does not Fourier-diagonalize the
finite simplex.  It uses the infinite `A3` translation completion, or the
identical local principal twelve-neighbor stencil at an all-positive
`N>=4` interior cell.

## Rank, sectors, and scalar quadratic seed

Since `zz^dagger-I` has eigenvalues `(3,-1,-1,-1)`, it is invertible.  Full
column rank of `P` then yields

\[
 \operatorname{rank}{\cal K}(k)=4,
 \qquad \ker{\cal K}(k)=\ker P^T\cong E_2.
\]

At zero program-frame wavevector the exact sectors are

\[
 {\cal K}_{A_1}=18,
 \qquad {\cal K}_{E_2}=0,
 \qquad {\cal K}_{T_2}=-2.
\]

For normalized `u=1_6/sqrt(6)`, every half-port belongs to three pair
coordinates, so

\[
 u^T{\cal K}(k)u
 ={3\over2}\left(\left|\sum_az_a\right|^2-4\right).
\]

The tetrahedral root tight frame supplies the quadratic change
`-8a_star^2|k|^2` after the `A1` contraction.  Multiplication by the negative
mixed coefficient therefore gives the exact `+16/9` coefficient displayed
in the frozen theorem.

This is a physical, record-qualified response organized on a supplied
program frame.  It is not yet a physical spatial momentum kernel.  The
shared-child address/incidence query, inverse atlas, cocycle, and physical
step calibration remain explicit gates.  The rank-four result also leaves a
real two-channel `E2` propagation deficit rather than a complete metric
stiffness.

## Scope verdict

GL6Z closes the branch-normalized entrance-time multicell connector jet, its
exact infinite/local-principal `A3` program-frame symbol, a fixed propagated
`E2` null, and an isotropic `A1` quadratic program-frame seed.

It does not close an unconditional formation average, full selected-apparatus
Schur kernel, physically authenticated spatial atlas, common physical cone,
complete six-channel stiffness, frequency kernel, Ward/Bianchi structure,
infrared operator, Ricci or Einstein form, gravity, or `G`.

**Audit verdict: PASS.**
