# INTEGRATION_classes -- folding the reachable-classes family (C-87, C-90) into the URM

T-54 build-out, family: **classes**.  Builder session 2026-08-21.  The T-46 pattern:
ported machinery in one new module, sealed numbers gated in one new check block, this
document tells the registrar exactly what to wire and nothing here edits anything.

## Files delivered (new; nothing existing was touched)

| file | what it is |
|---|---|
| `model/classes.py` | the classes layer: venue entry (D-25 gate), mu_c located by the resolvent route (exact rationals), the Gamma-priced series with proven tails, the certified critical kernel, the cross-dimension discriminator, the class-verdict boolean triple.  Ports from LANE_T44_B_WORLD/t44b_lib.py + t44b_world.py, LANE_T44_A_CORNER/t44a_lib.py, and the sealed verifier's deepening loop (VERIFY/adv_verify.py E-section).  Imports nothing heavier than `fractions`/`math` at import time (numpy only if `corner_venue` is called, via the model's own `geometry.Torus`). |
| `model/checks_classes.py` | `run_classes_checks(check)` -- 47 checks in the validate_geometry.py idiom; standalone `python3 model/checks_classes.py` exits 0 iff all pass.  Measured runtime: **98.6 s** (budget ~120 s; the cost is the M=2800 kernel deepening that the C-90 coefficient bracket requires). |

Standalone run 2026-08-21: **47 PASS, 0 FAIL, 98.6 s**, exit 0.

## ProjectModel methods to add (the CLASSES layer)

Append after the GEOMETRY layer block in `model/project_model.py`, same delegation
pattern (lazy `import classes as CC` inside each method; model dir is already on
sys.path).  Signatures and docstring text, ready to paste:

```python
# ---------------------------------------------------------------- CLASSES (T-54; C-87/C-90)
# Every method below delegates to model/classes.py, where the sealed-lane machinery
# lives (exact ints/Fractions on the measurement path; certified tails; computed
# booleans -- ported, not reimplemented).  Each docstring names its claim row, sealed
# source, and owners.
def coupling_venue(self, name, adj, provenance=None, tier="world", sector=None):
    """C-87 observation entry: a NEW venue graph (a record surface's access geometry)
       enters as an adjacency structure through the D-25 provenance gate -- world tier
       REFUSED without a pinned source, corner tier must self-declare 'DEF-A'
       (LANE_T44_B_WORLD S0 / LANE_T44_A_CORNER S0, sealed).  adj: list over nodes of
       [(neighbor, multiplicity), ...]; sector optionally declares the venue limit
       ("Z3"/"Z2"/"Z1") for the evidence instruments."""
    import classes as CC
    return CC.venue(name, adj, provenance=provenance, tier=tier, sector=sector)

def world_coupling_venue(self, n):
    """C-90: the world venue (earned D=3) -- the n^3 census-grain torus with the PINNED
       provenance (GR1 grains, face adjacency, T42_C/T43_B lineage; one walk step = one
       grain-boundary crossing = one unit of writer weight, the Gamma price C-80/O-54).
       Returns (VenueGraph, cells, idx).  (LANE_T44_B_WORLD S0, sealed.)"""
    import classes as CC
    return CC.world_venue(n)

def corner_coupling_venue(self, Lx, Ly):
    """C-87: the corner venue (earned D=2) -- the dual lattice of the model's OWN ported
       carrier: plaquettes of geometry.Torus adjacent iff their supports share a carrier
       edge, computed from the plaquette masks alone, multiplicity kept; DEF-A
       self-declared.  (LANE_T44_A_CORNER/t44a_lib.py plaquette_adjacency,
       verbatim-in-substance; sealed row sums exactly 4 on (4,6) and (3,7).)"""
    import classes as CC
    return CC.corner_venue(Lx, Ly)

def chain_coupling_venue(self, L):
    """C-87: the D=1 control venue (cycle C_L, DEF-A).  (t44a_lib.py cycle_adjacency,
       verbatim; sealed: row sums exactly 2, mu_c = 1/2 the venue's own.)"""
    import classes as CC
    return CC.chain_venue(L)

def critical_price(self, venue):
    """C-87/C-90: the venue's OWN critical coupling price mu_c, LOCATED by computation,
       never a literal (D-8): Perron row-sum candidate 1/deg, (I - mu_c A) annihilates
       the constant vector, exact-rational resolvent SINGULAR at mu_c (the D-15 zero)
       and SOLVABLE at mu_c*(19/20), mu_c*(21/20) beside it (the positive controls; on
       the sealed deg-6 venue these are exactly the sealed 19/120 and 7/40).  Sealed
       values: 1/6 world D=3, 1/4 corner D=2, 1/2 chain D=1.  DECLINES on a venue that
       is not degree-regular.  Owners: Perron-Frobenius/Gershgorin standard.
       (LANE_T44_B_WORLD S0/S2; LANE_T44_A_CORNER S2, sealed.)"""
    import classes as CC
    return CC.mu_c_of(venue)

def reachable_class(self, venue, mu, evidence=False):
    """C-87: THE REACHABLE-CLASS VERDICT -- the computed boolean triple (exponential,
       critical, divergent) for a declared coupling price mu, by exact rational
       comparison against the venue's own computed mu_c; exactly one True whenever mu_c
       is located.  evidence=True additionally runs the class's own instrument on the
       declared sector: subcritical booleans (ratios <= 1 - 1/20, Cauchy, power
       exclusion), the critical 1/d kernel signature (INV window, deepening-stabilized),
       or the divergence witness.  (LANE_T44_B_WORLD taxonomy + S2/S3/S4, sealed.)"""
    import classes as CC
    return CC.class_verdict(venue, mu, evidence=evidence)

def coupling(self, mu, target, K):
    """C-87 class (1): the Gamma-priced coupling G_mu(d) = sum over admissible strings
       of mu^weight on the D=3 sector -- exact partial sum (Fraction) with the exact
       geometric tail (6mu)^{K+1}/(1-6mu); leading term N_min mu^d with w_min = d the
       confinement cost (C-80/O-54 standing).  Owners: walk generating functions,
       Spitzer/Lawler; comparison tier only.  (t44b_lib.py series_3d, verbatim.)"""
    import classes as CC
    return CC.series_3d(mu, target, K)

def critical_kernel(self, targets, M):
    """C-90: the regularized critical kernel a_M(x) = sum (N_2m(0)-N_2m(x))/36^m at the
       COMPUTED mu_c = 1/6, exact rationals, with certified tails beside it
       (classes.diff_tail_bound / abs_tail_bound; honest for M >~ 1000).  The D=3
       critical member is the 1/d POWER LAW: exponent bracket contains 1, coefficient
       bracket [0.476369, 0.487321] onto owner 3/(2 pi) (comparison).  Owners: Polya
       1921 transience, Watson 1939 G(0), Spitzer P26.1 coefficient -- comparisons only.
       (LANE_T44_B_WORLD S4 + register row C-90, sealed.)"""
    import classes as CC
    return CC.crit_kernel_3d(targets, M)

def class_discriminator(self, K2=6000, K1=80000):
    """C-87: the cross-dimension discriminator -- ONE instrument (doubling-increment
       ratios of the regularized critical kernel), three venues, three pairwise-DISJOINT
       declared windows: D=1 LINEAR [9/5,11/5], D=2 LOG [4/5,5/4], D=3 INV [2/5,3/5].
       The critical class is the earned dimension's own; mu = 1/6 is critical on D=3 and
       subcritical on D=2.  (LANE_T44_B_WORLD S5, sealed.)"""
    import classes as CC
    return CC.discriminator(K2=K2, K1=K1)
```

## Where the checks chain in

`model/checks_classes.py` exposes `run_classes_checks(check)` exactly like the C-72
check block: it issues `check(name, cond, detail)` calls and owns no counter.  Because
the block costs ~99 s (the C-90 M=2800 deepening), the recommendation is a SEPARATE
validator stage so the fast validators stay fast:

```python
# model/validate_classes.py  (registrar-created)
import os, sys, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checks_classes import run_classes_checks
n_pass = n_fail = 0
def check(name, cond, detail=""):
    global n_pass, n_fail
    ...  # the validate_geometry.py idiom, verbatim
run_classes_checks(check)
# then chain validate_geometry.py (which chains validate_project.py), conjunction on
# exit codes: sys.exit(0 if (n_fail == 0 and chain_ok) else 1)
```

Alternative (registrar's call): append `from checks_classes import run_classes_checks;
run_classes_checks(check)` inside validate_geometry.py before the chain -- one edit,
but it makes the geometry validator ~99 s slower.  The separate stage is the
recommendation.

## The observation-entry story

The URM directive (the principal, 2026-08-21): the model is the framework new
observations are added INTO.  For this family:

1. **A new record surface's access geometry** (a new venue) enters as an adjacency
   structure with declared provenance through `coupling_venue()` -- the D-25 gate
   REFUSES a world-tier venue without a pinned source and requires DEF-A
   self-declaration from idealisations.  Checked end-to-end in the check block on a
   venue on no sealed grid (n=3): the universal-cover wrap identity, BFS == earned
   separation, and the venue's OWN mu_c re-earned by the resolvent route.
2. **A new law** would enter as a further layer method whose validator gate reproduces
   its lane's sealed numbers -- the pattern this fold-in itself follows.
3. **A new external number** (a measured coupling price mu for a real medium; a measured
   falloff rate) enters as a declared exact rational and receives the computed verdict
   triple against the venue's own mu_c, plus evidence booleans; external measured values
   compare against computed intervals ONLY in labeled comparison gates with a stated
   tolerance (COMP_TOL = 1/25), the way Watson/Spitzer/OZ anchors are gated now.
   No class label is ever asserted; every label is a computed boolean.

## What the check block gates (sealed sources)

- **mu_c = 1/6 exact three ways** (Perron row-sum; resolvent singular at 1/6, solvable
  at 19/120 and 7/40; sector sandwich) -- LANE_T44_B_WORLD S0/S2.
- **The exponent bracket contains 1**: axis pairs at M=1400 == sealed
  [0.453369, 0.458705] and [0.469869, 0.507905] (contains 1/2, INV window) --
  t44b_world.OUT.txt S4.
- **G in [1.503919, 1.554391]** with Watson inside (comparison), increments
  0.007287 -> 0.005159 -- t44b_world.OUT.txt S4.
- **d*G(d) onto [0.476369, 0.487321]**: the M=2800 deepened coefficient == register row
  C-90 (LANE_T44_B_WORLD/VERIFY/adv_verify.OUT.txt E7), owner 3/(2 pi) inside
  (comparison), strictly tighter than the sealed M=1400 bracket [0.467230, 0.500095].
- **The cross-dimension discriminator**: D=2 [0.967771, 0.986029], [0.953910, 1.028887]
  in LOG; D=1 [1.940862, 2.047476] in LIN; windows pairwise disjoint; mu=1/6 subcritical
  on D=2 at [0.063870, 0.063870] -- t44b_world.OUT.txt S5.
- **Subcritical rows** mu=1/12 and mu=1/8 == RESULT.json strings (ratios, power
  exclusion, G examples to 10dp) with computed class booleans; OZ owner comparisons
  labeled as such.  **Supercritical witnesses** at 13/72, 1/5, 1/4 == sealed term and
  ratio strings.
- **Corner/chain venues**: row sums 4 / 2; resolvent singular at 1/4, solvable at the
  sealed 1/8 and 23/100; chain mu_c = 1/2 -- t44a_corner.OUT.txt S0/S2.
- **API-fidelity probes beyond every gated range**: the resolvent identity EXACT at the
  unswept mu = 1/9; the untested kernel target (10,0,0) strictly interior; the untested
  venue n=3 end-to-end; kernel_pass == crit_kernel_3d at spot depth (one instrument,
  two drivers).
- **D-15 throughout**: every zero (leading-stratum, parity, resolvent-singular) gated
  with a positive control beside it.  **D-8**: no literal on any decision path; sealed
  anchors gated as the sealed anchors they are, stated as such in the check names.

## Deviations and findings the registrar should know

1. **The lane's certified kernel tails are only honest for M >~ 1000**: the edge term
   `EDGE_C * B5 * (M-2)^{-3/2} * RHO^{M+1}/(1-RHO)` (RHO = 199/200) only decays past
   M ~ 1000; at M = 350 it dominates the increments and certified intervals go negative.
   Not a defect in the sealed lane (which used M = 1400/2800) -- but it means a CHEAP
   certified-interval evidence tier is impossible with the ported tail machinery.  The
   evidence tier (`critical_evidence_3d`, used by `class_verdict(evidence=True)`)
   therefore uses POINT kernel values with a computed deepening-stabilization gate
   (M/2 -> M movement <= 1/8 of the smallest increment), honestly scoped in its
   docstring: the certified statements live at the sealed depths in the check block.
   The evidence tier also only claims the INV-window signature on the (2,4)/(4,8)
   increment pair -- the exponent-contains-1 claim needs the deepest pair and stays at
   the sealed-depth gates.
2. **`kernel_pass` is an additive extension** of the ported `crit_kernel_3d`
   (per-target stop depths + snapshots, the sealed verifier's E-section loop) so one
   ~90 s pass serves the M=1400 and M=2800 gates simultaneously; it is gated equal to
   the verbatim port at spot depth in the check block.
3. **The corner venue is built from the model's own ported carrier** (geometry.Torus
   plaquette masks) rather than by importing the sealed lane's o54c_lib -- same
   construction rule (shared carrier edges, multiplicity kept), gated against the
   sealed row sums and resolvent numbers.
4. The generic solvable-beside probe points mu_c*(19/20) and mu_c*(21/20) reproduce the
   sealed 19/120 and 7/40 exactly on the deg-6 venue (gated); the corner lane's own
   sealed probe points 1/8 and 23/100 are gated separately with the same instrument.
5. Venue provenance lives in `classes.WORLD_VENUE_PROVENANCE` (module-level, pinned);
   it is a venue-graph provenance, not a RecordSurface provenance, so it was NOT added
   to `project_model.PROVENANCE` (that registry is for surfaces).  If T-55 wants one
   provenance registry for both, that is a registrar decision.

## Next step (no route closes without one)

The named next piece is unchanged from the sealed lane's own: what earns criticality --
the record contains no Gamma-internal reason for mu to sit AT the venue's mu_c
(O-58 N2/N3, the E1/E2/E3 writing-tier ensembles, C-91).  When that lands, its
machinery enters this same layer as a venue-side method with a validator gate; the
`class_verdict` triple is the object it will decide occupancy OF.
