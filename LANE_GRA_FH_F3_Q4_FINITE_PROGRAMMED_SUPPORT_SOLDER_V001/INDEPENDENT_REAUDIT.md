# Independent hostile re-audit -- amended finite programmed q4/F3 support solder

**Lane:** `GRA-FH-F3-Q4-FPSS-V001`

**Date:** 2026-08-27

## Verdict

`ACCEPT_AFTER_INSTANTANEOUS_VS_STROBOSCOPIC_AND_FD_SLICE_REPAIR__NO_REMAINING_FH_MATERIAL_DEFECT`

The amended packet now proves the advertised finite programmed preparation
and quarantine with the required qualification: the raw ungated BS06 flip is
absent from, or continuously cancelled in, every Hamiltonian interval for
which instantaneous blank-nonedge invariance is claimed.  The hostile re-audit
found and repaired two residual defects.  First, a merely stroboscopic echo
was still grouped with switching off even though it can leave and later return
to the blank sector.  Second, the FD saturated-incidence carrier block needed
the PESC `K`-gated flip off as well as the raw flip.  Both ceilings are now
explicit and executable checks pass.

The first `INDEPENDENT_AUDIT.md` is preserved unchanged as historical custody.
Its unqualified sentence that nonedges are invariant after formation is
superseded by the canonical amended theorem and this re-audit.

## 1. Exact finite preparation

For every fixed supplied `N`, the stars-and-bars and padding census replays
exactly:

\[
 |S_N|={N+3\choose3},\qquad
 |S_{N+1}|={N+4\choose3},\qquad
 |E_N|=4|S_N|.
\]

The append edges are distinct, every active parent has degree four, and a
child's degree equals its number of positive coordinates.  Embedding both
layers at width `M=|S_(N+1)|` leaves exactly `{N+3 choose 2}` parent guards and
`M^2-4|S_N|` explicit nonedges.

The programmed formation uses distinct endpoint-owned writer slots and
fresh relation targets on every edge, so its factors commute.  Each inherited
FPMH dilation is reversible, the parallel KEEP route is a reversible
`L_e <-> K_e` move with all route history retained, and no register is reset,
traced, or identified.  The optional pulse

\[
 U_{KX}=\prod_e\exp(i\pi P_e^KX_{n_e}/2)
\]

is identity on `K_e=0` and equals `iX_(n_e)` on `K_e=1`.  It therefore maps
the blank incidence word exactly to the supplied q4 edge word when the raw
flip and every other noncommuting term are isolated for that pulse.  The
address map, edge list, cap, source tokens, hardware, switching, clock, work,
and physical port matrices remain supplied program data; no autonomous graph
selection is inferred.

## 2. `K`/`n` typing and qualified hold

The literal FPMH link factor is reused as the F3/PESC active incidence factor
`n_e=a_e`; it is not duplicated.  The support memory `K_e` remains a separate
binary factor.  In basis order `|K,n>`, the relevant qualified local hold is

\[
 H_e^{\rm qual}=-hP_1^K\otimes X_n
                 +\Delta I_K\otimes n
                 +\text{terms diagonal or multiplicatively controlled in }n.
\]

Let `Q_00=|0,0\rangle\!\langle0,0|`.  Direct matrix replay gives

\[
 (I-Q_{00})H_e^{\rm qual}Q_{00}=0,
 \qquad [H_e^{\rm qual},P_1^K\otimes I_n]=0.
\]

Thus the blank nonedge is reducing and the authenticated `K` word is
conserved.  In contrast,

\[
 [H_e^{\rm qual},I_K\otimes n]\ne0
\]

when the `K`-gated actuator is on: active `n` may evolve even though `K` does
not.  This exactly confirms the repaired `K`/`n` distinction.

Adding the raw BS06 term `-h_N I_K\otimes X_n` gives

\[
 (I-Q_{00})(I_K\otimes X_n)Q_{00}\ne0,
\]

so raw evolution immediately leaks out of the blank block.  FH17--FH19 now
require this term to be exactly zero in the hold generator, either switched
off or continuously cancelled.  A finite echo which returns the state only
at selected times proves a stroboscopic return map, not continuous invariance,
and is explicitly excluded from the Hamiltonian claim.

The corrected full nonedge sector is the tensor product of the one-dimensional
`|K_e=0,n_e=0>` blocks over all nonedges, tensored with the unrestricted
remaining factors.  Diagonal degree/incidence energies preserve it, and
BS09/BS11 terms vanish on a nonedge because they carry `n_e`.  Blank guard
carriers also remain blank: every incident link is a quarantined nonedge,
their onsite term is diagonal in occupation, and formation/copy couplings are
off.

## 3. History factorization

The projector `Pi_(p_N,E_N)` fixes the one orthogonal program and the retained
`K` word, but not the active `n` word.  The repaired theorem now states

\[
 [H_{\rm hold},\Pi_{p_N,E_N}]=0
\]

before writing the fixed-block factorization.  This follows only for the
qualified programmed hold: all BQ4/compiler history-writing couplings have
ended, the program sector is not coherently superposed, and remaining
controller/port evolution is independent of formation history.  Independent
port Hamiltonians may evolve their own factors but act identically on support
and retained history.

This is not a theorem about an arbitrary source-off F3 Hamiltonian.  If
different support programs were coherently superposed, the controlled F3
evolution could entangle with program/support and the single history-blind
block would fail.  ASSC conservation is also only passive: it preserves every
wrong support word as exactly as the intended q4 word.

## 4. FD and FE downstream statements

The optional saturation pulse prepares `n_e=1` on every eligible q4 edge.
That word is not invariant while the PESC `K`-gated actuator remains on.  The
FD claim has therefore been repaired to require both the raw ungated and
`K`-gated incidence flips exactly zero during the carrier comparator.  All
remaining incidence terms are diagonal in `n`, so the saturated word is then
an invariant block and BS09 restricts exactly to

\[
 -t\begin{pmatrix}0&B_N^\dagger\\B_N&0\end{pmatrix}.
\]

This closes only the finite physical support/off-diagonal block.  It does not
supply the positive child/parent detuning, autonomous maintenance, or a
collective phase.

The FE boundary no-go is unchanged.  The extreme child `(N+1,0,0,0)` has only
one eligible parent, so no subgraph of the raw slab can have degree two at
every active vertex.  A regular or periodic completion is supplied boundary
physics.  FD saturation has parent degree four and the FE ice sector requires
degree two; their same-`n` intersection remains empty.

## 5. Search for residual raw-flip assumptions

The canonical FH theorem and amended self-audit no longer assume that the raw
ungated flip preserves quarantine.  The historical first audit does contain
the obsolete unqualified invariant sentence; it remains frozen for custody
and is expressly superseded here.

The current FJ response theorem explicitly turns the raw flip off and uses
only the gated actuator.  FF explicitly states that the raw flip destroys
nonedge invariance.  The programmed-Floquet FI theorem ends/refocuses raw
incidence-flip pulses before freezing its support.  No current dependent
theorem located by repository-wide search requires the raw ungated flip to
preserve a blank nonedge.

One separate dependent wording caution remains: FI's saturated carrier pulse
must also include the PESC `K`-gated actuator among its refocused incidence
terms.  FH now makes that requirement exact for every saturated FD
comparator.  This does not defect FH, but FI should inherit that explicit
wording when its own packet is next amended.

The shared experiment registers still display the historical `37/37` FH
count.  They were not edited under this task's no-shared-file boundary and
should be refreshed during the next authorized register update.

## 6. Verification and final ceiling

The strengthened verifier reports `49/49 checks passed`.  In addition to the
finite combinatorics and preparation pulse, it now constructs the exact local
`K\otimes n` matrices, verifies blank-block invariance, verifies `K`
conservation and possible `n` evolution, proves raw-`X` leakage, and verifies
that the saturated FD slice conserves `n` only with both flip actuators off.
It also guards the instantaneous-versus-stroboscopic, history-projector, and
claim-ceiling language.

The accepted theorem remains a finite, programmed, ideal-model compiler and
qualified hold.  It does not derive the physical array, addresses, edge list,
controls, calibration, autonomous/scalable support selection, stability
basin, positive detuning, global FE phase, visible electromagnetism, tensor
gravity, or `G`.
