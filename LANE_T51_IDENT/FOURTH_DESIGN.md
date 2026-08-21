# T-51 FOURTH DESIGN — GEOMETRY RESPONSE — **THE RIGIDITY THEOREM** — V001 — 2026-08-21

**Status: this completes the owed fourth design of T-51's instrument round.** The original
fourth design (geometry-response angle) failed in flight and never reached the judge; the
register records the debt in C-92's row and in `FIELD_INSTRUMENT_V001.md` §4. This document
delivers the angle after the judgment, so it does not re-litigate the ranking: the Second Lump
runs first, and nothing here displaces the judged sequence. What this design was commissioned
to ask, it asks — and the exploration answered decisively enough that, per the commission's
own branch instruction, the deliverable is **the rigidity theorem with its proof mechanisms
and computations, plus the statement of which venue could be responsive**, with the round's
seven design elements carried in §5.

**Computations:** `FOURTH_rigidity.py` (this lane), output `FOURTH_rigidity_OUT.txt` —
**22/22 gates pass**, every quantity by or against the sealed machinery in `model/geometry.py`
(C-78 `dW_class_matrix`, C-80 `Torus`/`coset_min_np` conventions, C-81 `_subgroup_in`).
Gates marked CERTIFICATE are algebraic identities, labeled per the round's own critique
discipline and never counted as controls; every other gate could have come out the other way
(D-8), and two of them reproduce or extend the auditor's counterexample with independent
machinery.

**Discipline:** D-1 absolute — no classical form is required or tested anywhere below; every
shape is an output. D-24 — separations are earned quantities (`w_enc`-connected carries;
`d_gen` descriptive only, per the judgment; lattice coordinates appear only as construction
labels). D-15 — two-way controls, with certificates explicitly demoted. The principal's
directive quoted in C-92 binds: the mechanism for accumulation — and here, for response — is
whatever it proves to be; a computed no-response theorem is a result, not a failure.

---

## 0. THE QUESTION, AND THE COMPUTED ANSWER

**The question the angle owns:** C-78 earned `d_W(A,B)` as minimal admissible-writer
boundary-crossing cost. C-77's shape says a field is the aggregate of boundary-shaping terms
appearing as geometry. So measure the geometry directly: does WRITING content change
`d_W(A,B)` for a pair straddling or neighboring it? Is the earned metric itself responsive —
or provably rigid?

**The computed answer, in one sentence:**

> **On stabiliser venues every earned-geometry quantity — `d_W`, the full writer-cost
> landscape (`w_min`, `w_enc`, degeneracies, whole weight histograms), cut-rank, region
> entropy at the exact tier — is EXACTLY invariant under every content write, at every size,
> by three algebraic mechanisms; and the response that exists instead lives one tier below
> the integers, in the real-valued amplitude functional over the metric's own rigid class
> decomposition — which is precisely the quantity the Second Lump reads.**

The earned metric has no field side on the exact surface. The field side the surface does
have is *the measure over the rigid geometry's classes*. This is the geometry-currency mirror
of C-80's zero-beyond-contact and of Cert-Shadow's static arm (CERT = cut-rank,
configuration-independent — the theorem the judge already ordered registered): three
rigidity statements, one family, three currencies (information, certification, geometry).

---

## 1. THE RIGIDITY THEOREM — THREE MECHANISMS, EACH PROVED AND EXHIBITED

### R-I — Label tier: `d_W` cannot see content, because writers act by translation.

**Statement.** For every stabiliser venue and every content write `c`:
`d_W(s ⊕ c, s' ⊕ c) = d_W(s, s')` — the earned metric is a function of the label difference
alone.

**Proof sketch (3 lines).** An admissible writer `w` acts on the configuration labels by
`s ↦ s ⊕ a(w)` where `a(w)_i = sp(w, Z̄_i)` is a symplectic pairing against fixed operators —
a function of `w` alone, never of the state. Hence the set of writers implementing `s → s'`
is `{w : a(w) = s ⊕ s'}`, the same set for every translate `(s ⊕ c, s' ⊕ c)`. A content write
IS a translation of the base point; the coset, its weights, and its minimum are unchanged. ∎

**Computed exhibit (CERTIFICATE — cannot fail given R-I, labeled so):** per-pair
recomputation of the full `d_W` matrix at `L = 2, 3, 4` *without* the difference-class cache,
equal to the sealed class matrix (`class_min = {L, L, 2L}` at each `L`), and
translation-invariant over all 64 (s, s', c) triples at each size.

### R-II — Writer-class tier: the cost landscape around content is state-blind.

**Statement.** The probe's label-exact writer class
`{w : sp(w, A_p) = 1 (both probe stars), sp(w, q) = 0 (every remaining star),
sp(w, p) = 0 (EVERY plaquette — no source co-write)}`
is defined by symplectic pairings against *operators*. The written value `b` lives in the
*state*. No term of the enumeration reads `b`; therefore the entire weight histogram —
`w_min`, `w_enc` (old-spec and connected), `N_min`, everything — is identically equal in the
written and unwritten sectors. Moreover the same constraint list defines the class whether
the source plaquettes sit in `H` (carrier absent) or are removed (carrier present): the
constraint `sp(w, B_src) = 0` is an admissibility condition in one venue and the no-co-write
label condition in the other, so **carrier presence does not move the probe's landscape
either** — only co-writing operators are added by carrier presence, and those implement a
different label action.

**Computed exhibits.** Full pure-Z-span weight histograms at two 3×2 placements and two 3×3
placements; mixed-span (Y-dressed, closed-X-loop group added) enumeration **could have
produced lower minima and did not** (gated, could fail, passed); the auditor's counterexample
**independently reproduced from separate machinery**: old-spec `w_enc = 3`, connected
`w_enc = 4` at the auditor's exact 3×2 placement.

### R-III — Region tier: where the state responds, the integers still do not.

**Statement.** (a) Cut-rank (C-81's `CERT`) of any region is computed from the stabiliser
group and the region alone — state-free, so content-independent (this is Cert-Shadow's static
arm, confirmed here on the hole venue: cut-ranks 7 and 5 for the two declared regions).
(b) At the exact tier (`λ = 0`), a region's *reduced state* responds to a distant content
write **iff the region supports an enclosure detector** — a Z-loop separating the source
holes; equivalently, every admissible representative of the content's writer crosses the
region oddly. Where no detector fits, the write is transported by an operator supported
entirely off the region and `ρ_R` is literally unchanged. Where a detector fits, the region
*reads the bit* and the response is maximal. (c) Even there, the write is Pauli transport —
a unitary conjugation — so the spectrum of `ρ_R` is conserved: **the region's entropy, the
integer tier, is exactly unmoved even where the state responds at trace distance 1.**

**Computed exhibits (both directions could fail; both gated):**

| λ = 0, source ((1,1),(2,1)) | trace-dist(ρ_written, ρ_unwritten) | ΔS (bits) |
|---|---|---|
| R_away (no detector fits) | **2.9e−16** | 0.0e+00 |
| R_cross (encloses one source hole) | **1.000** — the region reads the bit | **0.0e+00 exactly** |

---

## 2. WHERE THE RESPONSE ACTUALLY LIVES — COMPUTED BESIDE THE THEOREM

The rigidity theorem says what the field side is *not*. The same table computes what it *is*:
switch on the program's own priced mediator (`V = λ Σ_e Z_e`, content-blind,
`[V, B_src] = 0`, written/unwritten two exact sectors of one Hamiltonian) and the
**real-valued amplitude functional over the rigid classes moves with content**, in exactly
the structure the Second Lump's design names:

| quantity | tier | 3×2 src ((1,1),(2,1)) | 3×2 src ((0,1),(2,1)) | 3×3 src ((1,1),(2,1)) | 3×3 src ((1,2),(2,2)) |
|---|---|---|---|---|---|
| `w_direct` | integer — **rigid** | 1 | 1 | 1 | 1 |
| `w_enc` old-spec | integer — **rigid** | 3 | 3 | 4 | 4 |
| `w_enc` CONNECTED | integer — **rigid** | 4 | 3 | **5** | 4 |
| Δ onset order (both sectors) | integer — **rigid** | 0.99 / 0.99 | 0.99 / 0.99 | 0.99 | 0.99 |
| `F` onset order (first rung) | dynamical shadow of the integer tier | **3.995** | **2.980** | **4.988** | **3.980** |
| `F` at λ = 0.05–0.064 | amplitude — **responds** | +1.13e−4 | −8.29e−4 | −2.6e−5 | +1.06e−4 |

Three findings in this table beyond the theorem:

1. **The onset law follows the CONNECTED `w_enc` of each placement** — 4 where connected = 4,
   3 where connected = 3, 5 where connected = 5, ~4 where connected = 4 — an independent,
   placement-resolved confirmation of the auditor's connectivity gate, from machinery that
   shares nothing with the critics' scratchpad rebuilds.

2. **A NEW instance of the auditor's counterexample at the commissioned venue size.** At 3×3
   with probe stars (0,0),(1,0) and source plaquettes ((1,1),(2,1)), the old-spec coset
   minimum is 4 but the CONNECTED minimum is 5 — and the measured onset is 4.99. The
   disconnected-representative cancellation is **not a 3×2 accident**; the pre-commission
   `w_enc`-connectivity gate is live at the exact venue LANE_T51_A will run on, and V3's
   bracket must be scored against the connected value there.

3. **The amplitude-tier response also exists in the information currency:** at λ = 0.05 the
   region-entropy differential between sectors is **+1.80e−3 bits on the crossing region
   against +4.6e−6 bits on the away region** — a content-attributable, geometry-concentrated
   entanglement response (an output; no shape required of it), sitting exactly where R-III
   says the exact tier forbids any integer to move.

Winding annex (V5's discipline, computed): at λ = 0.05, `F` across the four Z-winding
sectors is +1.13e−4, −1.86e−5, −1.23e−4, +1.86e−5 — sign is a winding-sector quantity on
small tori, exactly as the judged repair mandates before any sign attribution.

Probe doublet witness `|⟨0|A_probe|1⟩|` = 0.990–1.000 across every cell (test-object
identity certified per the Second Lump's convention).

---

## 3. WHICH VENUE COULD BE RESPONSIVE — THE OWED STATEMENT

All three rigidity mechanisms use the exactness of the stabiliser algebra and nothing else:
R-I needs writers to act by exact label translation; R-II needs the class constraints to be
pairings against operators with integer weight; R-III needs writes to be exact Pauli
transport. **None of the three survives real-valued costs.** Therefore:

- **Within stabiliser venues** the only field side available is the one computed in §2: a
  real-valued functional over the rigid class structure (the Second Lump's `F`; the Tier-2
  signed kernel; the entanglement differential). The door on the integers is closed at every
  size — by proof, not by the smallness of the venues reached.

- **The metastable/world tier is the venue where the metric itself could respond.** There the
  writer cost is not an integer weight but a free energy — clause (iv′)'s floor
  `kT ln2 + ΔE_config` — and `ΔE_config` reads the configuration, i.e. written content,
  whenever the surface's energy function couples regions. C-91's own vocabulary already
  prices the consequence: a content-induced shift of the per-step bias `dE` moves the writer
  kernel off or onto criticality by the closed-form gap `ln(μ_c/μ) = ln(1 + e^{dE/kT}/ℓ)` —
  **bias is mass** — so a content-responsive world metric would be measurable in the
  surface's own criticality currency, not in an imported one.

**Named first computation for the responsive venue** (no route closes without a next step):
in the world model's writing tier (C-91's E1/E2/E3 ensembles), compute the per-step writer
cost landscape around a written vs unwritten neighboring record — *does one written record
shift `dE` for an adjacent write?* If yes, the shift feeds the closed-form gap and the
world-tier geometry response is a number in earned units; if no, the rigidity theorem
extends to the world tier at the computed range and that too registers. Motivates-only
flag, zero evidential weight (D-25 posture): the census's own record surfaces are engineered
*against* exactly this response — cell-to-cell program disturb in NAND, adjacent-track
interference in HDD are content-dependent write-cost shifts fought as defects — which is why
the world tier is the right place to look and why finding either sign of the answer is a
finding about real surfaces' operating margin, not a toy artifact. Empirical contact
remains zero (X-4).

---

## 4. WHAT THIS DOES AND DOES NOT SETTLE

**Settled (computed here):** the exact-tier earned metric is not responsive to content — not
`d_W`, not the writer-cost landscape, not cut-rank, not region entropy at λ = 0; the response
exists and is confined to the amplitude functional over the rigid classes; the connectivity
gate binds at 3×3 with connected `w_enc = 5` at one commissioned-size placement.

**Not settled, carried honestly:** (1) MEDIATION — the amplitude-tier response is read
through the declared, priced mediator; whose field it is, C-77 must absorb (C-80's standing
sentence; no control can decide it). (2) The world-tier computation of §3 has not run;
nothing here measures a responsive metric, only where one could live. (3) The R-III "iff" is
proved as mechanism and exhibited at one venue with one region pair per direction — the
class-avoidance argument is general, but a swept-region version has not been enumerated.
(4) Contact-scale caveats of small tori carry verbatim from the Second Lump's brief; shape
claims at scale belong to the unbuilt signed connected kernel. (5) Per standing discipline
the rigidity result reads two ways — the exact tier has no metric response, or "metric
response" is a differently-shaped concept on this surface — and both readings travel with
the register entry.

---

## 5. THE SEVEN ELEMENTS (the round's standard), AS THE RIGIDITY INSTRUMENT

**Name: THE RIGID RULER.** Angle: start from the earned geometry itself and condition it on
content — the field as a response of the metric, or the theorem that there is none.

1. **PROBE** — the earned-geometry battery itself is the probe: `d_W` per-pair (C-78's
   instrument), the label-exact probe-writer class with its full weight histogram (C-80's
   conventions), cut-rank and region reduced states (C-81's tier); plus the minimal hole-pair
   doublet as the dynamical shadow-probe (the Second Lump's own test record, reused so the
   two designs' columns sit in one table).

2. **OBSERVABLE** — the conditioned differentials of every earned integer
   (`δd_W`, `δw_min`, `δw_enc`, `δN_min`, `δcut-rank`, `δS` at λ = 0) — **all computed
   identically zero**, three mechanisms, 22/22 gates — and beside them the amplitude
   differentials (`F`, onset ladders, `δS` at λ > 0) — **computed nonzero**, in energy units
   and bits, at places where the content is not.

3. **RESPONSE CLAIM** — existence-only and two-sided (D-8): either some earned integer moves
   under a content write (any nonzero cell falsifies the theorem at that size), or no integer
   moves and the response is confined to the amplitude tier. No shape, exponent, or
   classical form is required of either branch (D-1). The claim that survived: rigidity at
   the integer tier, response at the amplitude tier, onset locked to connected `w_enc`.

4. **CONTROLS** (D-15, two-way, certificates demoted) — live members: the mixed-span
   (Y-dressing) minimum could have beaten pure-Z and did not; the onset battery could have
   tracked old-spec `w_enc` (3, 4) instead of connected (4, 5) and did not; the region pair
   runs both directions (away: zero required; crossing: response required); Δ-onset equality
   across sectors could have split. Certificates, labeled and never counted: R-I translation
   invariance; R-II histogram equality (no `b` input exists); cut-rank state-freeness;
   unwritten ≡ source-absent portage (carried from the Second Lump's table, same
   construction).

5. **VENUE** — exact: toric 3×2 (32-dim sectors, calibration) and the commissioned 3×3
   (256-dim sectors), sector-exact, no sampling; the theorem itself venue-unbounded by
   mechanism (the proofs are size-free; the computations are exhibits at reached sizes).
   Responsive continuation: the world/metastable tier, §3, with its first computation named.

6. **MEANING** — POSITIVE READING (obtained): the record surface's exact tier carries a
   rigidity theorem — geometry is not where its field side lives, and this is a computed
   division of labor (Γ's integers rigid; the coupling's amplitudes responsive), the
   C-80-shaped structure appearing in the geometry currency; it strengthens the Second Lump
   by proving its reading is the *only* place a field can live on this surface. NEGATIVE
   READING (still open at the world tier): if the §3 computation also returns no response,
   the rigidity family closes over both tiers and T-17's ledger carries that sentence beside
   the mediation clause; per C-92's directive neither outcome is a failure against an
   imported standard.

7. **BORROWED**, owners named, program variable first (relevance test) — hole-pair records →
   C-80's construction (Bravyi–Kitaev quant-ph/9811052); rigidity-of-code-distance under
   logical action → standard stabiliser formalism (Gottesman; Calderbank–Shor–Steane
   lineage), comparison only after computation; Pauli-transport invariance of reduced
   spectra → textbook unitary invariance (owner: none needed; stated for completeness);
   entanglement entropy of stabiliser states = cut-rank → Hamma–Ionicioiu–Zanardi,
   Fattal et al., the instrument C-81 already names; enclosure/crossing detector →
   Aharonov–Bohm / toric mutual statistics, already named by the claim; linked-cluster
   cancellation of disconnected representatives → Goldstone linked-cluster / Kato–Bloch
   degenerate PT, the owners the judgment already bound to the `w_enc` gate. Dropped for
   failing the relevance test: metric perturbation, lensing/deflection, elasticity of space —
   no program variable maps to them.

**Honest risks:** (1) the theorem's reach is exactly its scope — DEF-A/stabiliser exactness;
its world-tier fate is uncomputed and could go either way; (2) the crossing-detector
convention rides the sealed `dual_path_x` route and a declared winding-class convention —
winding annex computed, sign attribution stays behind V5's sweep; (3) the R-III iff is
mechanism-plus-instance, not a swept enumeration; (4) first-rung onset estimators drift with
λ (ladders reported in full; the drift direction is the known finite-λ curvature); (5) a
judge may rule that "geometry response" should have meant the world tier all along — in
which case this design's exact-tier theorem is the reason, stated with proofs, that the
world tier is where the question lives.

---

## 6. HANDOFF NOTES FOR THE REGISTRAR

- This artifact clears the fourth-design debt named in C-92's row and
  `FIELD_INSTRUMENT_V001.md` §4 ("the fourth design angle... is owed").
- The rigidity theorem here and Cert-Shadow's static arm (which the judge ordered registered
  regardless of ranking) are one family and can share a register entry or sit adjacent:
  static certifiability has no field side (certification currency); the earned metric has no
  field side (geometry currency); both mirror C-80's zero-beyond-contact (information
  currency).
- **For LANE_T51_A before it fires:** at 3×3, probe stars (0,0),(1,0), source plaquettes
  ((1,1),(2,1)): connected `w_enc = 5` against old-spec 4, measured onset 4.99 — the
  connectivity gate is confirmed live at the commissioned venue and V3 must score against
  the connected value per placement.
- Named next steps: (a) the world-tier responsive-venue computation of §3 (does a written
  record shift `dE` for an adjacent write — the geometry response in C-91's own currency);
  (b) the swept-region enumeration closing R-III's iff; both are lane-sized, neither blocks
  the judged sequence.

*Files: `FOURTH_DESIGN.md` (this document), `FOURTH_rigidity.py` (computation, sealed
machinery read-only), `FOURTH_rigidity_OUT.txt` (output, 22/22 gates).*
