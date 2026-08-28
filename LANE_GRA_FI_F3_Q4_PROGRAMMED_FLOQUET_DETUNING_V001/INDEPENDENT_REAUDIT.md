# Independent hostile re-audit -- dual-flip-free programmed Floquet detuning

**Lane:** `GRA-FI-F3-Q4-PFCD-V001`

**Date:** 2026-08-27

## Verdict

`ACCEPT_AFTER_EXPLICIT_PESC_FLIP_AND_STROBOSCOPIC_SUPPORT_REPAIR__EXACT_FLOQUET_RESULT_SURVIVES__NO_REMAINING_FI_MATERIAL_DEFECT`

The exact finite Floquet spectrum, uniform principal-branch gap,
dressed-parent functional calculus, and operator kernel bound remain correct
after the physical slice is amended.  The amendment is load-bearing: both the
raw BS06 incidence flip and the PESC `P^KX_n` actuator must be exactly zero in
the joint generator throughout both carrier pulses.  Returning the incidence
word only at the end of an echo is not sufficient because BS09 carrier hopping
depends on the instantaneous `n` word.

The historical `INDEPENDENT_AUDIT.md` is preserved unchanged.  Its generic
requirement that "incidence flips" be isolated is now made explicit and
operator-complete by the canonical theorem and this re-audit.

## 1. Exact support-slice audit

The old `V_0--V_1` q4 incidence word is saturated and the new `V_1--V_2`
incidence word is blank.  `K_e` and `n_e` are distinct physical factors.  In
the local basis `|K,n>`, the two incidence-changing generators are

\[
 X_{\rm raw}=I_K\otimes X_n,
 \qquad
 X_{\rm PESC}=P_1^K\otimes X_n.
\]

On an old support edge, `|K=1,n=1>`, both operators move the state out of the
saturated block.  On a new blank edge, `|K=0,n=0>`, the PESC actuator vanishes
but the raw operator still leaves the blank block.  Therefore the exact common
condition for both slabs is

\[
 X_{\rm raw}=0,qquad X_{\rm PESC}=0
\]

in the pulse generator, either by switching them off or by continuous exact
cancellation.  Every remaining displayed incidence energy is diagonal in
`n`, BS09 and BS11 contain `n` only as a multiplicative control, and all
formation/copy terms are off.  The old saturated and new blank incidence
words are consequently invariant throughout both pulse intervals.

This condition does not collapse `K` into `n`.  The authenticated support word
is passively retained on `K`; the active word `n` is separately held fixed by
the qualified controller slice.  Concrete control matrices, cancellation
work, timing, recoil, failure, and calibration remain supplied physical
antecedents.

## 2. Why an incidence echo is insufficient

The distinction between continuous support invariance and Floquet sampling is
exact.  Consider one active link and a two-state carrier with hopping operator
`T`.  The desired saturated-link interval is

\[
 U_{\rm desired}=e^{i\eta T}.
\]

Now flip `n=1` to `n=0`, wait under the incidence-gated Hamiltonian
`n\otimes T`, and flip `n` back.  The incidence bit returns exactly, but the
carrier sees

\[
 U_{\rm echo}=I\ne e^{i\eta T}
\]

for nonzero generic `eta`.  The strengthened verifier constructs this joint
unitary explicitly.  Thus a stroboscopic echo which is identity only on `n`
cannot be counted as a carrier-identity switching interval in `T_F`.  Only an
interval proved identity on the joint `n`-plus-carrier sector, or continuous
exact cancellation that leaves the target generator, can be inserted without
recomputing the Floquet unitary.

The word "stroboscopic" in PFCD remains correct in a different sense: it
describes the carrier quasienergy law sampled after complete two-pulse cycles.
It does not license support excursions within either pulse.

## 3. Exact physical reductions under the corrected slice

With the old incidence word invariant, the BS09 hop interval restricts
exactly to

\[
 H_H=-tX_B,qquad
 X_B=\begin{pmatrix}0&B_N^\dagger\\B_N&0\end{pmatrix}.
\]

The uniform old-slab onsite term is a scalar in the one-carrier sector and is
factored out.  During the next-slab interval, old-slab transfer is off, the
new incidence word remains blank, and `V_2` remains carrier-blank.  The same
inherited onsite generator then restricts exactly to

\[
 H_D=\epsilon_\psi\Pi_C,qquad
 U_D=\Pi_P-i\Pi_C
\]

for `tau_D=pi hbar/(2 epsilon_psi)`.  No staggered bulk coefficient, BS06 link
detuning, `K_eT_e` term, new field, or symbolic-port rescue is used.

## 4. Floquet spectrum recomputation

For a singular value `sigma>0` of the injective q4 incidence matrix, put
`a=eta sigma`.  The exact two-pulse block remains

\[
 U_F=
 \begin{pmatrix}
  \cos a&i\sin a\\
  \sin a&-i\cos a
 \end{pmatrix}.
\]

Its eigenvalues are

\[
 \lambda_P=e^{i(\omega-\pi/4)},\qquad
 \lambda_C=e^{-i(\omega+\pi/4)},\qquad
 \omega=\arccos(\cos a/\sqrt2).
\]

Every q4 column has four entries and every row at most four, so `||B_N||<=4`.
The ceiling `0<|eta|<=pi/16` gives `|a|<=pi/4` and
`omega in [pi/4,pi/3]`.  Parent phases lie in `[0,pi/12]`; coupled and dark
child phases lie at or below `-pi/2`.  The principal-branch separation is
therefore at least `pi hbar/(2T_F)` exactly as claimed.  The incidence
amendment changes the physical antecedent for this block, not its spectrum.

## 5. Dressed-parent kernel

The parent-connected eigenphase depends only on
`K_N=B_N^dagger B_N=4I+A_N`, so the pulled-back principal Floquet logarithm
remains

\[
 H_P^F=-{\hbar\over T_F}
 \left[\arccos\!\left({\cos(\eta\sqrt{K_N})\over\sqrt2}\right)
       -{\pi\over4}I\right].
\]

The previously audited scalar inequality on
`0<=z<=pi^2/16`,

\[
 {z\over2}-{z^2\over6}
 \le
 \arccos(\cos\sqrt z/\sqrt2)-{\pi\over4}
 \le {z\over2},
\]

again gives

\[
 0\preceq H_P^F+{\hbar\eta^2\over2T_F}K_N
 \preceq {\hbar\eta^4\over6T_F}K_N^2.
\]

The leading nonconstant term is the controlled q4 sibling kernel.  No
massless phase, continuum time binding, tensor response, gravity law, or
numerical `G` follows from this finite programmed result.

## 6. Instantaneous-versus-stroboscopic dependency search

Every canonical FI use of the saturated or blank incidence word now points to
the dual-flip-free slice.  PFCD-1, PFCD-2, and PFCD-3 all inherit that premise.
The cycle-time definition admits only intervals that are identity on the joint
incidence-and-carrier sector.  No remaining FI theorem statement derives an
exact carrier block from an incidence word that is allowed to leave and later
return.

The control architecture remains supplied rather than autonomous.  If a
future apparatus uses a nontrivial echo instead of continuous cancellation,
its complete joint `n`/carrier/control unitary must be calculated afresh; the
current closed-form Floquet result cannot simply be imported.

## 7. Verification and final ceiling

The strengthened verifier reports `57/57 checks passed`.  It retains all q4
incidence, exact spectrum, phase-window, dark-mode, dressed-parent, and kernel
checks.  It additionally verifies old-slab leakage under both flips,
next-slab leakage under the raw flip, vanishing of the gated flip at `K=0`,
dual-flip-free conservation of `n`, distinct `K/n` typing, and the explicit
incidence-return/carrier-unitary counterexample.

The final result is an exact finite programmed Floquet substitute for the
positive-detuning role on a supplied dual-flip-free incidence slice.  Static
source-off detuning, autonomous scheduling, finite-error robustness,
collective-phase stability, physical clock/refinement binding, visible
electromagnetism, tensor gravity, and `G` remain open.
