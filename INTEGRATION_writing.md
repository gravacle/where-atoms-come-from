# INTEGRATION_writing — folding the WRITING TIER (C-91) into the Universal Record Model

T-54/T-55, family: **writing**.  Built under the principal's directive of 2026-08-21: the
URM is the world model — the framework new observations are added INTO.  The writing tier
is homed as an extensible object: a new record surface enters through the provenance gate,
the writing laws enter as layer methods with validator gates, and a new external number
enters as a gated comparison.  This document is the registrar's integration spec; the
builder edited **no existing file**.

## Files delivered (the only writes)

| file | what it is |
|---|---|
| `model/writing.py` | the module: ported machinery from LANE_T48_A_DERIVATION / LANE_T48_B_CORNER / LANE_T48_C_WORLD (T-46 fold-in pattern; every function's docstring names claim row C-91, sealed source, owners) |
| `model/checks_writing.py` | the check block: `run_writing_checks(check)` in the validate_geometry.py idiom; standalone `python3 model/checks_writing.py` → **52 PASS, 0 FAIL, exit 0, ~4.5 s** |
| `INTEGRATION_writing.md` | this spec |

Sealed numbers reproduced exactly (gated as the sealed anchors they are, stated as such):
the four CTRL-LEAK determinants (269059/4000000; 653188856401/4096000000000;
196750721299/2560000000000; 886217035591602121/16384000000000000000), the energy census
1024/2048 and 2304/4608, the pair split 4/18 + 4/18 + 10/18 with conditional 1/4, the
bath column sums [1,1]/[3/2,1/2]/[2,0], the corner invariant tuple (1,2,2,2,1) on
42/48/50 links, mu_c located in-lane (1/4 corner, 1/2 chain, singular/nonsingular-beside),
the coset stratum (w_min, N_min) = (1,1), E1's induced-mu set {1/6} and partial sum 13,
E2's split (1/7, 2/7) at b = 1/2, E3a's mu table {1/7, 9/64, 3/22, 1/8, 1/10, 1/16} and
partial-sum anchor 39867016537742941/4096000000000000, the T44-B comparison row
mu(1/2) = 1/8, deg_NB = 5.  API-fidelity probes run BEYOND the gated range: fresh world
venue 6^3 at fresh (u, b); fresh kernel venue C_12; fresh lazy row c = 2/7; fresh corner
venue (4,5) including a fresh exhaustive coset scan; the field-kernel entry point at a
declared anisotropic field.

## A. ProjectModel methods to add (paste-ready; delegation style of the GEOMETRY layer)

Add under a new banner after the GEOMETRY section:

```python
# ---------------------------------------------------------------- WRITING (T-54/T-55; C-91)
# Every method below delegates to model/writing.py, where the sealed-lane machinery
# lives (exact Fraction kernels, venues from carrier supports, the coset instrument
# reused from geometry.py -- ported, not reimplemented).  Each writing function's
# docstring names its claim row, sealed source, and owners.
def writing_kernel_verdict(self, venue, c=0):
    """C-91 kernel tier: the invariant lazy-family kernel c*I + ((1-c)/deg)*A on a
       kernel venue ('C8' | 'T3' | 'T4' | 'Z27'), with conservation (double
       stochasticity), criticality (exact det(I - K)), link amplitudes, and the
       CTRL-LEAK control beside it.  The sealed identity: conserving <=> critical,
       per-crossing amplitude 1/deg identically for every trivial-writer share
       (LANE_T48_A_DERIVATION, sealed).  Owners: unital channels / Birkhoff
       (re-verified on explicit operators in-lane), Perron, Gershgorin."""
    import writing as WW
    from fractions import Fraction as Fr
    adj = {"C8": WW.ring_venue(8)["adj"], "T3": WW.plaquette_venue(3)["adj"],
           "T4": WW.plaquette_venue(4)["adj"], "Z27": WW.grain_venue(3)["adj"]}[venue]
    K = WW.kernel_uniform(adj, c)
    deg = WW.venue_degree(adj)
    return dict(deg=deg, conserving=WW.is_doubly_stochastic(K),
                critical=WW.crit_det(K) == 0,
                link_amplitudes=WW.link_amplitudes(adj, K),
                per_crossing=(1 - Fr(c)) / deg / (1 - Fr(c)),
                leak_det=WW.crit_det(WW.leak_kernel(adj, Fr(9, 10))))

def writing_uniformity(self, Lx, Ly):
    """C-91 corner tier: uniformity EARNED from the writer algebra on the (Lx, Ly)
       plaquette venue rebuilt from carrier supports -- the elementary writer's
       invariant tuple on every link, mu_c located in-lane (Perron row sums + exact
       resolvent singular at 1/deg, nonsingular beside), and the unique conserving
       member t* = 1/deg = mu_c (LANE_T48_B_CORNER, sealed at (4,6), (3,7), (5,5)).
       Owners: carrier Kitaev quant-ph/9707021; coset instrument o54c lineage
       (geometry.py); Perron/Gershgorin/Feller as in-lane."""
    import writing as WW
    from fractions import Fraction as Fr
    cv = WW.corner_venue(Lx, Ly)
    inv = WW.writer_invariants(cv)
    return dict(invariant_tuple=inv[0], identical=all(v == inv[0] for v in inv),
                n_links=len(inv),
                mu_c=WW.mu_c_locate(cv["rows"], 0, beside=(Fr(1, 8), Fr(23, 100))),
                conserving_member=WW.conserving_member(cv["rows"]))

def writing_transport(self, n, a):
    """C-91 world tier, E1 TRANSPORT (iv' literal): the energy-conserving writer on the
       n^3 census venue at per-attempt amplitude a -- conserving AND critical
       (mu = 1/deg = mu_c) at EVERY dE and barrier; verdicts computed, never narrated
       (LANE_T48_C_WORLD V1, sealed).  Owners: Stinespring (why row sums 1 is the
       dilation's structural property; the sums themselves computed)."""
    import writing as WW
    _c, _i, nbr = WW.torus3(n)
    return WW.transport_verdict(nbr, WW.ensemble_transport(nbr, a))

def writing_trail_retreat(self, n, u, b):
    """C-91 world tier, E2 TRAIL WITH RETREAT: the raw two-rate writer with
       backtracking kept -- conserving at every dE, uniform (1/deg) exactly at dE = 0,
       split (b, 1)/(5b + 1) beside (LANE_T48_C_WORLD V2, sealed).  Owners:
       Goldstein 1951 / Kac 1974 named for the persistence remark, comparison-only."""
    import writing as WW
    _c, _i, nbr = WW.torus3(n)
    return WW.retreat_verdict(nbr, WW.ensemble_trail_retreat(nbr, u, b))

def writing_trail_decay(self, n, u, b, counting="H1"):
    """C-91 world tier, E3 TRAIL WITH DECAY (the model's own erase channel): NEVER
       critical -- mu = b/(deg*b + 1), f0 and E_b dropping out exactly; the COMPUTED
       mass ratio mu_c/mu is returned; the closed form ln(1 + e^{dE/kT}/l) is checked
       against it by the validator, never sourced from it (LANE_T48_C_WORLD V3/V4,
       sealed; the DONE_WHEN control).  counting='NB' uses the venue's own
       directed-edge criticality reference (Hashimoto 1989, earned by row sums)."""
    import writing as WW
    _c, _i, nbr = WW.torus3(n)
    W = WW.ensemble_trail_decay(nbr, u, b, counting)
    return (WW.decay_verdict(nbr, W) if counting == "H1"
            else WW.decay_verdict_nb(nbr, W))

def writing_gap(self, s, n=4, den=10 ** 9):
    """C-91 + D-25, THE OBSERVATION ENTRY: a real record surface's written-trail mass
       gap.  REFUSES without provenance (build s through URM.surface); declines
       non-thermal surfaces.  The surface's b = exp(-dE/kT) is bracketed by exact
       rationals, the E3 kernel is built and measured exactly at both brackets (u
       independence re-computed at entry), the closed form is checked against every
       computed ratio, and the float gap is certified INSIDE the computed bracket."""
    import writing as WW
    return WW.surface_gap(s, n=n, den=den)
```

(`writing_field_kernel` is deliberately NOT added as a model method yet — the C-93
computation it serves is commissioned separately; the entry point lives in the module,
see C below.)

## B. Where the checks chain in — the decision, as option cards

**Option 1 (recommended): a new `model/validate_writing.py`** mirroring
validate_geometry.py — header naming the sealed sources, local `check()`, call
`run_writing_checks(check)`, then chain `validate_geometry.py` (which itself chains
validate_project.py), exit 0 iff conjunction.
  - keeps validate_geometry.py byte-identical (its 31-gate count is cited in VERIFY_T46);
  - the chain becomes validate_writing → validate_geometry → validate_project, one entry
    point for the whole model.

**Option 2: append to `model/validate_geometry.py`** — after the C-71/C-72 formation
section, before the summary: `from checks_writing import run_writing_checks;`
`run_writing_checks(check)`.
  - one validator file, but the T-46 gate count changes (31 → 83) and VERIFY_T46's
    recorded counts no longer describe the file.

Either way, `reproduce.sh` / model.sha256 sidecars are the registrar's to extend.  The
check block is import-clean from any cwd (it inserts model/ on sys.path itself).

## C. The observation-entry story — how NEW observations of this family's kind enter

1. **A new world surface** (a physicist's own written medium: dE, E_b, T, f0 in SI):
   built through `URM.surface(...)` (D-25 — refused without provenance), then
   `writing_gap(s)`.  Its Boltzmann dial b is bracketed by exact rationals and the E3
   kernel is RUN at both brackets: the surface's mass gap arrives certified inside a
   computed bracket, with the closed form checked against the computed ratios at entry —
   never read off the formula alone.  Non-thermal surfaces decline (None), as the model's
   laws do.

2. **C-93's named next step** (commissioned separately — NOT run here): *does one written
   record shift dE for an adjacent write?* — the responsive-venue computation, in C-91's
   own currency (bias is mass).  ITS ENTRY POINT IS `writing.kernel_pos_field(nbr,
   diag_of, link_of)`: a written pattern's dE field enters as the link map
   b(x, d) = exp(-dE(x, d)/kT), rationally declared; the SAME instruments then read the
   result — conservation by `srow_sums` (a computed row-sum deficit IS a computed mass),
   induced amplitudes and uniformity by `extract_pos` (a dE shift appears as computed
   link anisotropy), criticality by the exact resolvent/determinant kit.  The check block
   already gates this entry point's fidelity (constant field == kernel_pos entrywise;
   a declared anisotropic field computes conserving-but-non-uniform).  Whatever the
   computation returns registers as the surface's own law — accumulation is whatever it
   proves to be.

3. **A new ensemble construction or venue** (the O-58 N3 design point, the N5 anisotropy
   frontier): a new honest construction enters as a constructor beside E1/E2/E3 returning
   dict rows; the verdict instruments and the D-15 controls (CTRL-LEAK, CTRL-BIAS-LINK)
   apply unchanged, and its numbers enter through a new check() gate beside the sealed
   ones.  A new venue enters through `ring_venue` / `plaquette_venue` / `grain_venue` /
   `corner_venue` / `torus3` — all built from carrier supports with deg computed, never
   declared.

## D. Registrar notes (scope, and two named items)

- **Ported vs. left in the lanes**: LANE_T48_B's K2 symmetry-orbit family and its
  exhaustive `all_automorphisms` machinery are NOT ported.  The registered finding they
  carry (conservation + symmetry leaves an (R-1)-parameter freedom on two-scale venues,
  every member critical) is homed at kernel tier by the ported CTRL-BIAS-LINK
  counterexample with its invariance witness; the anisotropy frontier itself is O-58 N5,
  open, lane-scoped.  Porting the orbit machinery belongs to whichever lane computes N5.
- **Gate-count note for the register**: the C-91 ledger row's evidence text says "World
  tier (LANE_T48_C_WORLD, 42 gates)"; the sealed OUT, RESULT.json, and the T48_D judge's
  own re-count all say **36** (G01–G36; the five V-verdicts are conjunctions of gates,
  not gates).  126/126 and 112/112 check out as written.  Registrar's call whether the
  row's narrative wants an erratum; no artifact disagrees with any number, only the row's
  count of them.
- **Runtime**: full check block ~4.5 s (two exhaustive coset scans included), well under
  the 120 s budget; module import has no side effects.
- **Collision safety**: no existing file edited; no git run; the three files above are
  the entire write set.
