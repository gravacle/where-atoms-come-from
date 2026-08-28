# Self-audit -- programmed F3/q4 Floquet carrier detuning

**Lane:** `GRA-FI-F3-Q4-PFCD-V001`

**Date:** 2026-08-27

## Verdict

`SUBSTANTIVE_EXACT_PROGRAMMED_RESULT__DUAL_FLIP_FREE_INCIDENCE_SLICE_EXPLICIT__STATIC_AND_AUTONOMOUS_CEILINGS_RETAINED`

The packet closes the finite **programmed** detuning role, not the static or
autonomous detuning problem.  Its positive result uses only the already
declared F3 carrier onsite and hopping generators, the finite q4 support
compiler, one additional admitted blank F3 layer, and an owned orthogonal
pulse schedule.

## Exact checks and ownership

1. **Static boundary retained.**  The source-off regular F3 block still has no
   child-only onsite term.  The theorem does not alter or repeal `CLDNG-3` and
   does not call the static `CCMAC` Hamiltonian derived.
2. **Why the scheduled pulse works.**  In the next active slab, the physical
   child layer is present and the earlier parent layer is absent.  With the
   next layer carrier-blank, the existing uniform onsite term restricts to
   `epsilon_psi Pi_C` exactly.  No logical label or q4 count is assigned an
   energy.
3. **No hidden bulk interaction.**  The construction does not use a BS06 link
   detuning as a carrier energy, does not add `K_e T_e`, and does not fill the
   symbolic BS12 slot with an unlisted stagger.  It switches only already
   admitted bulk generators.  Exact controller couplings, matrices, and
   calibration are still supplied physical antecedents within the inherited
   controller/clock/work port types; the symbolic BS12 slot does not prove
   them.
4. **Positive-parameter ceiling.**  `epsilon_psi>0` and `t!=0` are required.
   The parent family permits zero values, so the theorem is an admitted-domain
   construction, not a statement about every parameter point.
5. **Schedule ceiling.**  Repetition, exact isolation/refocusing, and the
   pulse durations are supplied by the controller program.  The result does
   not derive a natural autonomous oscillator or zero-error finite ramp.
   During both carrier pulses, the raw BS06 incidence flip and the PESC
   `P^KX_n` actuator are exactly zero in the joint generator, either switched
   off or continuously cancelled.  A merely stroboscopic incidence echo does
   not preserve the exact carrier unitary because hopping is controlled by
   the instantaneous `n` word.
6. **Work is owned, not erased.**  Switching, timing, heat, recoil, reset, and
   failure remain assigned to the inherited physical ports.  Their concrete
   matrices, work values, and calibration are not calculated here.
7. **No history conflation.**  BQ4 and compiler histories remain retained.
   On the fixed orthogonal program, the carrier unitary factors from those
   histories; no trace converts a formation append into hopping.
   `K` support and active `n` remain distinct: `K` is passively retained,
   while old `n` is held saturated and next-slab `n` blank by the dual-flip-
   free pulse slice.
8. **Exact spectrum.**  The SVD of `B_N` reduces the cycle to independent
   two-level blocks plus `dim(C)-dim(P)` dark child modes.  Under
   `|eta|<=pi/16`, every parent-connected phase is separated from all child
   phases by at least `pi/2` per cycle.
9. **Branch ownership.**  `H_P^F` is the principal Floquet logarithm pulled
   back through the dressed parent spectral isometry.  It is not the raw
   compression of `log U_F` onto a subspace which is not invariant.
10. **Kernel bound.**  The exact function of `K_N=B_N^dagger B_N` obeys the
    stated positive operator remainder inequality.  The leading sibling
    term is controlled, while all higher common-child returns remain present.
11. **Physical ceiling.**  No collective coordinate, masslessness, stable
    thermodynamic phase, physical clock/probe binding, pair-field dynamics,
    tensor constraint, universal stress vertex, RGRL-B, gravity, or `G` is
    promoted.

## Known non-defects

- The additional blank layer is not a new field type; it is another layer of
  the already declared finite F3 allocation.  Its blank state and boundary
  realization are nevertheless supplied antecedents, not emergent results.
- The child/parent contrast is schedule-relative.  That is permitted because
  F3 already declares layer order and controller timing as physical program
  data; it is not advertised as an observer-independent static bulk mass.
- Quasienergy is defined modulo `2pi hbar/T_F`.  The pulse ceiling places all
  displayed phases in one explicit principal interval, so the separation is
  unambiguous.
- The result is stroboscopic in the carrier Floquet sense only.  Continuous
  incidence-block invariance is a premise inside each pulse; a return of `n`
  at the end of an echo is not substituted for it.
- The leading coefficient depends on calibrated `t`, `tau_H`, and `T_F`.
  The theorem derives its form and sign, not a numerical gravitational
  constant or a parameter-free scale.

## Amendment and verification

The first independent audit correctly required isolation of incidence flips
but did not name the PESC `P^KX_n` actuator separately from the raw BS06
flip.  The amended theorem now requires both generators exactly zero during
both Floquet pulses and distinguishes continuous incidence-block invariance
from a merely stroboscopic return.  The historical first audit is preserved;
a fresh independent re-audit owns this amendment.

The strengthened deterministic verifier passes `57/57`.  It replays the
exact q4/Floquet spectrum and kernel bounds and now also constructs the local
`K\otimes n` operators, checks both old-slab and next-slab flip leakage,
checks the dual-flip-free invariant slice, and gives an explicit echo in which
`n` returns while the carrier hop unitary changes.  Its stdout is retained in
`VERIFICATION.txt`.

## Disposition

The finite detuning obstruction has split cleanly.  The static source-off
stagger remains absent, while an exact finite programmed Floquet substitute
is constructible with existing owned interactions and ports on the qualified
dual-flip-free incidence slice.  The next physics question is whether this
schedule is realized autonomously and stably in a collective phase; further
record machinery is not justified by this result.
