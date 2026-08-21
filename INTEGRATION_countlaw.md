# INTEGRATION_countlaw — the count-law layer (C-86) folded into the URM (T-54, family countlaw)

Builder: countlaw family lane, 2026-08-21, under the principal's directive: **the URM is the world
model — the framework new observations are added INTO.** This family homes C-86 (the surviving-
record count law) as an extensible layer: a new device census enters through the provenance gate,
the law enters as a layer method, and the sealed numbers are validator gates.

## Files (all new; nothing outside them was touched)

- `model/countlaw.py` — the machinery, ported from the sealed lanes (LANE_T47_A_WIDTH,
  LANE_T47_B_STAIRCASE; registered via LANE_T47_D_REGISTER). Importable standalone; `python3
  model/countlaw.py` prints the sealed six-record census as a demo. No file I/O, returns are data.
- `model/checks_countlaw.py` — the check block, validate_geometry.py idiom:
  `run_countlaw_checks(check)`. Standalone: `python3 model/checks_countlaw.py` → **40 PASS, 0
  FAIL, ~2 s** (far under the 120 s budget).
- this document.

The corner law `k = min_E v2(m_E)` stays homed at `model/count_law.py` — **referenced, not
duplicated**; countlaw's private `_v2` exists only so the T-31 control can exhibit the corner
instrument's reading beside the record-mode census.

## 1. ProjectModel methods to add (the registrar integrates; signatures + docstring text)

Both delegate to `model/countlaw.py` exactly as the GEOMETRY layer delegates to
`model/geometry.py`. Suggested layer header comment:

```python
# ---------------------------------------------------------------- COUNT LAW (T-54; C-86)
# Delegates to model/countlaw.py -- machinery ported from LANE_T47_A_WIDTH and
# LANE_T47_B_STAIRCASE (sealed), registered via LANE_T47_D_REGISTER. The corner law
# k = min_E v2(m_E) remains homed at model/count_law.py (C-14), referenced not duplicated.
```

```python
    def census(self, surfaces, t_m):
        """C-86: THE SURVIVING-RECORD COUNT LAW k(t_m) -- the URM's wholly-owned
           falsifiable census (LANE_T47_A_WIDTH sec. 6, LANE_T47_B_STAIRCASE secs.
           C/D/G, sealed; registered via LANE_T47_D_REGISTER).  surfaces: a LIST of
           RecordSurface objects, world-tier ones built through URM.surface() so D-25
           provenance rides on every row; t_m: the retention spec in s.  Counts a
           record only while BOTH its values are durable -- clause (ii') on the
           record's OWN Liouvillian mode, |lambda_record| <= 1/t_m -- so each record
           dies at its SHALLOWER value's escape, at the parameter-free drop time
           t*_i = f0^-1 exp((B_i-dE_i)/kT)/(1+e^{-dE_i/kT}).  Returns k by the
           instrument AND k_formula by the derived width delta_pop(t_m) =
           kT ln(expm1(B/kT - ln f0 t_m)) so the agreement is CHECKED, never assumed;
           plus the dated drop schedule, delta_coh = hbar/t_m, the departure term
           sum_over_dead tanh(dE_i/2kT) (remanence persists while records die), and
           the declined list (non-thermal surfaces are declared, never silently
           counted).  THE SIGNATURE IS THE C-76 GATE: (surfaces, t_m) and nothing
           else -- no width, tolerance, or clustering parameter exists on this path;
           checks_countlaw.py gates the unreachability.  Owners per C-86 (ownership
           PARTIAL): Neel/Street-Woolley/Sharrock own the activation window as a
           remanence-decay device; Charap-Lu-He and Weller-Moser the dE = 0 corner;
           Preisach-Neel rate-level asymmetry; ours the derived width, the two-face
           unification, the margin-free integer census, the departure term."""
        import countlaw as CLW
        return CLW.census(surfaces, t_m)

    def count_widths(self, s, t_m):
        """C-86: both faces of clause (ii') on one RecordSurface at one retention
           spec -- delta_pop from the record's own constants (B = E_b + dE, the exact
           convention map; None when t_m is beyond the symmetric bound exp(B/kT)/(2 f0),
           which IS the no-crossing condition), delta_coh = hbar/t_m, and the record's
           own instrument drop time beside its closed form (the check, never the
           source).  LANE_T47_A_WIDTH, sealed."""
        import countlaw as CLW
        B = s.E_b + s.dE
        return dict(delta_pop=CLW.delta_pop(B, s.T, s.f0, t_m),
                    delta_coh=CLW.delta_coh(t_m),
                    t_star=CLW.drop_time(s),
                    t_star_formula=CLW.drop_time_formula(s))
```

Deliberately **not** made ProjectModel methods: the T-31 multi-shell machinery (`t31_carrier`,
`t31_basin`, `t31_staircase`) — control instruments on an abstract carrier in carrier units, kept
module-side like the lanes' control batteries, cited from the validator. Note for the register
(carried from T-47 note 3): `grounded.clause_ii` misreads multi-shell carriers; on them the
record's own mode is the basin-lumped slow sector, which is exactly what `t31_basin` implements.

## 2. Where the checks chain in

`model/checks_countlaw.py` exposes `run_countlaw_checks(check)` taking the house
`check(name, cond, detail="")` callable. The builder supplied two equivalent hook-ups:

1. **Into the T-54 validator** (recommended once more families land): the umbrella validator
   defines its `check`, then
   ```python
   from checks_countlaw import run_countlaw_checks
   run_countlaw_checks(check)
   ```
   and keeps the validate_geometry pattern of chaining `validate_project.py` for the conjunction.
2. **Interim**: append the same two lines to `validate_geometry.py`'s gate section (it already
   owns a `check` of the right shape); the summary counts then absorb the 40 gates.

**Registrar disposition:** option 1 landed through `model/validate_urm.py`. COUNTLAW is
kept separately countable at 40 gates, including the consumption-side provenance-bypass
refusal, before the umbrella chains geometry and project/D-25. Standalone remains
available (`python3 model/checks_countlaw.py`, exit 0 iff all pass).

What the gates hold (40 checks; D-8: literals only as sealed anchors, stated as such, with a
computed comparison beside every one; D-15: every zero paired with a positive control):

- **A/B — width**: sealed instrument rates (worst rel 3.3e-13 vs the 12-digit sealed values);
  `delta_pop` == sealed `2.535960987846706e-20 J` at the (1.2 eV, 300 K, 1e9, 10 y) anchor;
  the naive-minus-derived-correction identity; the bisected instrument crossing landing on the
  closed form (1.2e-15); sealed 350 K anchors to 12 decimals; the sealed no-crossing row (None
  both routes, beside a crossing); `delta_coh` sealed anchor + exact identity; the one-modulus
  unification against the Liouvillian's own eigenvalue (sealed 7.211102550928); the coherence-face
  bisection == hbar/t_m with the sealed 4/2 slow-mode counts.
- **C — census**: the sealed six-record staircase [6,6,5,4,2,1,0] with its six sealed t*_i; both
  routes equal at all 1201 grid points on the sealed UNIFORM N=12 and LOGNORMAL N=10 ensembles
  (every sealed tau matched); the symmetric D-15 control flat {8,0} at the sealed threshold
  1.254112e+05 s == e^(B/kT)/(2f0) (the derived factor-2 corner); the departure carrier: k=4,
  dead grains exactly [4..9], departure == sealed 5.4339, zero on the symmetric control.
- **D — T-31**: basin taus and spectrum-read dE at every sealed eps row; Q a classical Markov
  generator (worst 3.73e-11 of the slow scale, the sealed value); the staircase [2,1,0] at
  eps=0.16; **the C-76 kill as control** — exact-multiplicity v2 = 0 at every eps>0 while the
  record-mode census reads 2 (v2 = 2 = k at the corner, count_law.py's reading); the sealed
  splitting/delta_coh ratio class 2.99e36..e41.
- **E — the C-76 gate**: `census` signature exactly `(surfaces, t_m)`; width/tol/delta/
  cluster_width/margin all TypeError beside the succeeding call; module-wide scan finds no public
  callable with a width-like parameter.
- **F — observation entry**: registry-provenance NAND surface through `URM.surface()` → census
  (alive at 1e3 s beside dead at 1 y, both routes agreeing); D-25 refusal of an unregistered
  surface beside the acceptance; consumption-side refusal of both a provenance-bypassed world
  surface and a falsely declared corner; `thermal=False` declared in `declined`, never silently
  counted.
- **G — beyond the gated range** (definition-not-shortcut): the step located at t*_3 by ±0.1%
  probes; a fresh N=5 ensemble with MIXED per-record f0 (no sealed lane ever mixed f0) — routes
  agree at 301 t_m; T-31 at eps=0.20 (beyond the sealed grid); a 500 K instrument crossing
  (beyond the sealed 200–400 K grid) landing on the closed form at 4.9e-16.

## 3. The observation-entry story (what the Saira/Woodside grounding lanes call)

A real device census — MFM grain-by-grain media data, MLC NAND Vt-vs-bake levels, Saira-class
single-electron boxes, Woodside-class DNA-hairpin two-state landscapes — enters in three steps,
all existing machinery:

1. **Provenance gate (D-25)**: each measured record becomes
   `URM.surface(name, mechanism, dE, E_b, T, f0, provenance="<pinned source>")` — measured
   constants, not fitted to the law; the gate REFUSES a world-tier surface without a pinned
   source (gated at CL-F2). The census rechecks this invariant and refuses bypassed objects
   (CL-F3); construction alone is not treated as a security boundary. Non-thermal surfaces may
   enter; the census will decline them by name rather than mis-count them (CL-F4).
2. **The law**: `model.census(surfaces, t_m)` returns the predicted integer count k, the dated
   drop schedule (each record's t*, first-to-die first), both widths, and the departure term —
   with `k` and `k_formula` returned together so the instrument/width agreement is checked on
   the user's own data, not assumed.
3. **The gated comparison**: the measured survivor count (or level-retirement times) lands beside
   the predicted k and ln t*_i as a validator check with a stated extraction tolerance — the
   F1/F2 falsifier classes of LANE_T47_D_REGISTER (JEDEC-class MLC bake data; exchange-biased
   media censuses run as survivor COUNTS, not drop-time asymmetry alone). The departure term is
   the discriminator against Sharrock-style remanence counting: remanence persists while records
   die (sealed 5.4339-of-10 anchor, gated).

No adjustable content anywhere on the path: the only inputs are the carrier's own measured
constants and the retention spec. **There is no way to hand the census a width** — that absence
is itself a gated property (section E), which is C-76's kill kept as a fence.

## 4. Notes for the registrar

- `census` treats a surface whose Arrhenius rates underflow floats as DECLINED (open_system's
  sealed convention), where LANE_T47_A's sweep shim treated rate 0.0 as durable; the divergence
  is declared in the `record_rate` docstring and follows the model's registered decline behavior.
- The two T-31 numerical floors (1e-9 Bohr-grouping, 1e-12 matrix-element zero) are the sealed
  lane's own declared instrument floors, module constants, not caller-adjustable.
- T-55 hygiene: neither new file constructs a raw `RecordSurface` — every surface in module demo
  and checks goes through `URM.surface()` (corner grids self-declare DEF-A).
- Named next step (no route closes without one): the F1 external-data first pass (T47_D's N1 —
  Cai-et-al.-class Vt-vs-bake extraction) now has its entry point: build the level surfaces via
  `URM.surface()` with pinned provenance and gate the measured retirement staircase against
  `census`'s schedule; that run is what C-86's PROVED bar awaits.

## 5. Surprises found while building (for the record, none blocking)

- `PROVENANCE` auto-lookup in `URM.surface` means a *registered* name with `provenance=None`
  is accepted (registry fills it); only unregistered names are refused. The checks gate exactly
  that behavior (CL-F1/F2) — worth keeping in mind when reading T-55's "substring test" item.
- The basin-lumped T-31 instrument reproduces the sealed worst generator deviation 3.73e-11
  bit-for-bit at the printed precision, and holds to eps = 0.20 beyond the sealed grid.
