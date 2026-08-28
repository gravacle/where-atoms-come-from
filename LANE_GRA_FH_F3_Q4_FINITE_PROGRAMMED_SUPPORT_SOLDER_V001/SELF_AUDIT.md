# Self-audit -- finite programmed q4-to-F3 support solder

**Lane:** `GRA-FH-F3-Q4-FPSS-V001`

**Date:** 2026-08-27

## Verdict

`PASS_49_OF_49__RAW_FLIP_AND_K_N_TYPING_CORRECTION_EXACT__READY_FOR_INDEPENDENT_REAUDIT`

The packet proves a narrow positive result: for every supplied finite q4 slab,
existing BQ4 label/append structure, F3 equal-layer site/link allocation, and
FPMH/PESC formation/KEEP/controlled-link gates suffice to prepare an exact
physical support encoding.  Padding and nonedges remain in explicit invariant
blank sectors under the qualified raw-flip-free hold, and that fixed-program
hold is blind to retained q4 and compiler histories.

## Load-bearing checks

1. **Census.** `|S_N|=binom(N+3,3)`, `|S_(N+1)|=binom(N+4,3)`, so an equal F3
   slab with `M=|S_(N+1)|` needs `binom(N+3,2)` parent guards.  There are
   exactly `4|S_N|` distinct eligible edges and `M^2-4|S_N|` nonedges.
2. **No factor conflation.** BQ4 count/order factors, complete F3 node
   allocations,
   PESC `L/K/G`, F3/PESC `n`, endpoint writer slots, and all controller/port
   factors remain distinct.  The theorem does not reinterpret one BQ4
   `Q_N` factor as many sites.  The FPMH link `a_e` is reused literally as
   F3/PESC `n_e`; it is not duplicated.  Any within-layer generic FPMH pair
   factors remain blank spectators and are not called F3 links.
3. **Reversibility.** Formation uses inherited FPMH unitary dilations, KEEP is
   an inherited reversible route, and the optional saturation step is a
   finite unitary `exp(i pi P_K X_n/2)`.  No register is erased or traced.
4. **Program status.** The address map, edge list, cap, hardware allocation,
   source tokens, schedule, and ports are fixed supplied inputs in one
   orthogonal program state.  Coherent support programs are excluded from the
   history-blind theorem.
5. **Quarantine.** During the exact controlled pulse **and throughout any
   claimed blank-nonedge hold**, the raw ungated BS06 flip is exactly zero in
   the hold generator, either switched off or continuously cancelled, with
   switching/work ports retained.  A merely stroboscopic echo does not prove
   instantaneous invariance.  Any active incidence flip in the hold is the
   admitted PESC `P^K X` actuator.  On `K=0,n=0` it vanishes, while all
   remaining displayed diagonal/link/carrier controls preserve blankness.
   The exact operator replay confirms both invariance and the raw-`X`
   counterexample.  Guard carriers have no eligible incident links and their
   formation couplings are off.
6. **History factorization.** The qualified fixed-program hold preserves its
   program/`K` projector and turns every history-writing coupling off.  The
   factorization is a declared hold premise with independent controller/port
   evolution, not a theorem about an arbitrary source-off F3 Hamiltonian.
7. **Retention ceiling.** ASSC proves all support sectors invariant, including
   wrong ones.  This is passive memory, not autonomous selection, correction,
   or a stability basin.
8. **FD scope.** Saturating the eligible links gives the exact finite q4
   append-incidence off-diagonal BS09 block when both raw and `K`-gated
   incidence flips are zero during the carrier comparator, so the saturated
   `n` word is invariant.  The positive uniform child/parent detuning, its
   ports, and all phase/scaling claims remain open.
9. **FE scope.** The extreme child `(N+1,0,0,0)` has eligible degree one, so
   the raw finite slab has no global `d_*=2` subgraph.  Only finite physical
   edge binding and deep-interior local inheritance are earned.  Periodic or
   boundary completion remains supplied.
10. **Same-`n` scope.** FD saturation has degree four at every active parent;
   FE requires degree two.  The compiler does not make those exact sectors
   coexist and introduces no second field or `K_eT_e` law.
11. **No gravity promotion.** No collective phase, visible EM, tensor mode,
    universal stress law, gravity, or numerical `G` is claimed.

## Correction incorporated after the first audit

The first independent audit correctly checked isolation of the preparation
pulse but repeated an overbroad hold statement: it omitted the raw ungated
BS06 flip from the nonedge invariant-block census.  A raw `X_e` does not
preserve `n_e=0`.  The canonical theorem and verifier now require that term
to be absent or continuously cancelled in the qualified hold generator and
explicitly test the counterexample.  They also distinguish instantaneous
invariance from a merely stroboscopic return.

## Residual risks deliberately left visible

- A fully calibrated apparatus must freeze physical port matrices, energies,
  timing errors, recoil, and work.  Logical/custody completeness is not a
  measurement of those quantities.
- Repeating the compiler for each finite `N` is not an autonomous `N`-uniform
  preparation law and does not establish thermodynamic support stability.
- A separately programmed periodic diamond graph uses the same hardware
  types but is not the raw BQ4 slab and cannot be promoted as BQ4-derived.
- Optional diagonal F3 terms are nonuniform at the raw child boundary.  The
  FD statement is exact for the off-diagonal q4 transfer block; the detuned
  Schur parent remains conditional on a genuinely supplied positive offset.

## Verification

The dependency-free deterministic verifier passes `49/49`.  Its captured
stdout is frozen in `VERIFICATION.txt`.  The manifest hashes the theorem,
verifier, verification capture, self-audit, historical first audit, and fresh
re-audit; it does not silently include or modify shared model files.
