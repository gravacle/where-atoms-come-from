# Self-audit -- finite programmed q4-to-F3 support solder

**Lane:** `GRA-FH-F3-Q4-FPSS-V001`

**Date:** 2026-08-27

## Verdict

`PASS_37_OF_37__SCOPE_CEILINGS_RETAINED__READY_FOR_MANIFEST_FREEZE`

The packet proves a narrow positive result: for every supplied finite q4 slab,
existing BQ4 label/append structure, F3 equal-layer site/link allocation, and
FPMH/PESC formation/KEEP/controlled-link gates suffice to prepare an exact
physical support encoding.  Padding and nonedges remain in explicit invariant
blank sectors, and the source-off fixed-program hold is blind to retained q4
and compiler histories.

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
5. **Quarantine.** During the exact controlled pulse, the raw ungated BS06
   flip and every noncommuting incidence/carrier term are off or refocused,
   with switching/work ports retained.  On `K=0,n=0`, the PESC
   flip vanishes and all remaining displayed link/carrier controls preserve
   blankness.  Guard carriers have no active incident links and their
   formation couplings are off.
6. **Retention ceiling.** ASSC proves all support sectors invariant, including
   wrong ones.  This is passive memory, not autonomous selection, correction,
   or a stability basin.
7. **FD scope.** Saturating the eligible links gives the exact finite q4
   append-incidence off-diagonal BS09 block.  The positive uniform child/
   parent detuning, its ports, and all phase/scaling claims remain open.
8. **FE scope.** The extreme child `(N+1,0,0,0)` has eligible degree one, so
   the raw finite slab has no global `d_*=2` subgraph.  Only finite physical
   edge binding and deep-interior local inheritance are earned.  Periodic or
   boundary completion remains supplied.
9. **Same-`n` scope.** FD saturation has degree four at every active parent;
   FE requires degree two.  The compiler does not make those exact sectors
   coexist and introduces no second field or `K_eT_e` law.
10. **No gravity promotion.** No collective phase, visible EM, tensor mode,
    universal stress law, gravity, or numerical `G` is claimed.

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

The dependency-free deterministic verifier passes `37/37`.  Its captured
stdout is frozen in `VERIFICATION.txt`.  The manifest hashes the theorem,
verifier, verification capture, and this self-audit; it does not silently
include or modify shared model files.
