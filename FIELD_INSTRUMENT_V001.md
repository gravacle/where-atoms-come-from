# THE FIELD-SIDE INSTRUMENT — V001 — 2026-08-20 (T-51 identification)

**Status: IDENTIFIED, NOT MEASURED.** This document records the judged outcome of T-51's design
round. It defines what a *field* can mean on the record surface, names the instrument commissioned
first, the sequence behind it, and the pre-registered decision rule. **Nothing here is a registered
number.** The first commissioned computation may return the negative closure, and that path is
genuinely reachable.

**Provenance.** Three instrument designs were produced from deliberately different starting points
(a fourth, the geometry-response angle, failed in flight and is owed — see the register entry); each
was attacked by a D-1/D-24 importer-hunter and a measurability skeptic, all six critiques returning
`VIABLE_WITH_REPAIRS` with no unrepaired D-1 violation found on any path; a judge who added no design
ranked them as a **sequence, not rivals** — the three probe different channels. The full designs,
critiques, and judgment are archived in `LANE_T51_IDENT/`.

---

## 0. WHAT A FIELD IS ALLOWED TO MEAN HERE

O-39's standing verdict: every lane has measured a *source*; none has measured a *field*. A field is
what a **test object feels at a place where the source is not**. Under this program's discipline that
means: a probe placed at **earned** separation (D-24 — `d_W`-class quantities computed by sealed
machinery, never a lattice coordinate), read by the surface's **own laws**, against the **unwritten
surface as positive control** (D-15), with the shape of any separation dependence an **output, never
an input** (D-1 — no classical form is required or tested).

**The binding structural fact, and it is not a defect:** C-80 proves no pure-Γ falloff exists — reach
belongs to the coupling. So *every* field instrument on this surface is a **mediated-response
instrument**. Whether a field read only through a declared, priced mediator is the surface's field or
the mediator's is a judgment C-77 must eventually absorb; no control can decide it, and there is
provably no probe-independent field value at a place on this surface.

---

## 1. THE INSTRUMENT COMMISSIONED FIRST — **THE SECOND LUMP**

*(test-record angle; the only design whose first computation was independently rebuilt by both of its
critics from the repository's own conventions, reproducing every exploration number to three
significant figures — the strongest pre-commissioning replication status any instrument in this
program has had.)*

**The probe** is the smallest, weakest record the surface admits: a minimal hole-pair record — two
adjacent removed star stabilizers whose shared bit is durable and whose admissible connector writer
has weight 1 (C-80's own Bravyi–Kitaev defect construction, sealed in `geometry.py`). It sits at
earned separation from a written hole-pair record of the **dual** type, stated in two earned
quantities computed by the sealed coset-minimum machinery: `d_gen` (the `w_min = d_gen` law) and
`w_enc`, the minimal admissible connector weight in the *enclosing* path class. No lattice coordinate
enters any claim.

**The mediator** is the program's own alpha-type generic local coupling `V = λ Σ_e Z_e` — declared,
priced (C-87), content-blind.

**The observable** — the probe's own laws, read by existing instruments:
- `Δ(b; D)`: the probe doublet's tunneling splitting in source sector `b`, energy units,
  sector-exact (the source bit is an exact quantum number under the conserving-quadrature mediator);
- `τ_P(b; D)`: the probe's lifetime from its own Liouvillian mode (`grounded.py`, C-69's instrument)
  — scoped to the two-level effective GKSL tier (the dense Liouvillian is ~68 GB at 256-dim sectors
  and does not run as first written);
- **the field reading**: `F(D) = Δ(b=−1; D) − Δ(b=+1; D)` — **the written-value differential of the
  probe's own law, in energy units, at a place where the written content is not.** This is literally
  O-39's named minimum: a two-lump induced energy as a function of earned separation.

**The mechanism, from the design exploration (sector-exact, toric 3×2; not yet a registered
number):** the source's record value multiplies the enclosure class of the probe's admissible writer
sum — holonomy, the program's own claim vocabulary. `F(contact) = −1.99e−3` and
`F(far) = +1.94e−4` at `λ = 0.05`; onset order in `λ` measured at 3.0 and 3.9 against enclosure
costs 3 and 4; unwritten control identical to source-absent at machine precision. **Artifact
structure is the cleanest achievable**: written and unwritten are two exact sectors of *one*
Hamiltonian, so the mediator's own furniture cancels in the differential by construction, and a
nonzero `F` requires an enclosure-class amplitude — content-attributable by algebra, not by control.

**One earning debt, found by the audit and gated before scoring:** the sign-constrained coset minimum
admits **disconnected** representatives that the dynamics cancels (computed counterexample at 3×2:
coset minimum 3, measured onset 4.0). `w_enc` must be redefined over **connected** enclosing strings
and gated before the onset verdict is scored.

### The first computation (LANE_T51_A, to be commissioned)

Venue: toric **3×3** torus (18 edges; sector-exact blocks ≲ 2^11), with the exercised 3×2 as
calibration. Probe star-holes at adjacent vertices; source plaquette-hole pair at the two distinct
earned separations 3×3 affords, each stated as `(d_gen, w_enc-connected)`. Mediator
`λ ∈ {0.02, 0.05, 0.10}`.

**Controls (D-15, all two-way):** (1) `F` written vs unwritten; (2) Γ-equivalent placement swap at
equal earned separation; (3) **winding-sector sweep — mandatory before any sign of `F` is attributed
to content** (on small tori the sign decomposes through the winding sector until this sweep says
otherwise); (4) onset-order bracket against connected `w_enc` at both placements; (5) back-action:
source-sector energies with/without probe below a declared tolerance. *Construction certificates* —
unwritten-equals-absent (algebraic identity certifying the port), X-quadrature probe quiescence
(commutation identity) — are reported as certificates and never counted as controls.

### The pre-registered decision rule (D-8: every verdict can come out the other way)

- **V1** `|F|` exceeds the same-table measured control floor beyond contact at the larger separation.
- **V2** the reading follows earned geometry under the placement swap.
- **V3** the onset-order bracket contains **connected** `w_enc` at both placements (scored only
  after the connectivity gate).
- **V4** back-action below the declared tolerance.
- **V5** sign attribution licensed only if the winding-sector sweep separates content sign from
  winding sign.

**ALL PASS →** the record surface has a field side *in the only sense it has earned*: a test record's
own law is modified at a place where the written content is not, by an amount set by the content's
value and earned separation, through a declared, priced, content-blind mediator — C-80's computed
division of labor, not a defect. Then commission, in order: **(a)** the two-source composition at 4×4
— **the mod-2 saturation kill, first, because it is the fastest way the angle dies**; **(b)** the
Tier-2 signed kernel over the connected/dynamical string class, classifying `F` vs `D` against
C-87/C-90's three computed classes per swept `μ`, criticality never assumed (O-58's bit is not
assumable); **(c)** the Tier-3 grounded `τ` readout. The surviving operator combination (source
`Z̄` × probe writer amplitude) feeds **O-53** its field-side candidate.

**V1 FALSE →** sweep the priced mediator family (quadratures, staggered signs, weight-2 local
terms); still false → **register the negative closure**: C-80's pure-Γ zero propagates to the
dynamical probe, action at earned separation does not exist on this surface without an inserted
mediator's own structure, **the emergence claim is exhibitable-not-provable, and T-17 carries that
sentence.** Per standing discipline the null reads **two ways**: no field present, or field-side
concept differently shaped.

---

## 2. THE SECOND LANE — **THE HANGING CLOCK** (bath-mediated)

A two-level probe with nonzero splitting — **a clock; a degenerate probe is exactly gauged away,
contrast ~1e−14, computed** — hanging off a structured bath, coupled only to bath operators, never
the carrier. Observable: the exact four-term Casimir-style subtraction
`E_ind = F(both) − F(source) − F(probe) + F(bath)` — O-39's named minimum through the thermodynamic
channel. Feasibility run: record contrast resolvable at every separation tested (1.3e−1 at contact to
3.3e−3 at bath distance 6) against an unwritten control decaying to ~1e−7. Commissioned **after** the
Second Lump's first computation returns, with its critics' six repairs binding (recorded in
`LANE_T51_IDENT/`), chief among them: `d_sep` on the chain medium is **DECLARED/UNEARNED** and that
scope travels mechanically with every separation claim; and the odd/even-in-`s` split is
pre-registered as two named readings — **PRESENCE** (even part, O-39's quantity) and **VALUE** (odd
part, which screens). Its non-duplicated deliverable: whether an object at separation comes to *hold
record bits* or only an energy response — that split, if found, is new registrable structure.

## 3. THE THIRD LANE — **CERT-SHADOW** (certifiability), plus one theorem to take now

`SHADOW(A←C) = CERT_dyn(A | C written) − CERT_dyn(A | C unwritten)` in **bits** — the currency this
surface has earned. Runs third, renamed **FORMATION SHADOW** until a same-table bridge to C-81's
`CERT` is computed, with its three repairs binding. **Independent of ranking, its static arm is a
cheap, provable rigidity theorem on sealed T-43 machinery and should be registered on its own:**
on stabiliser venues `CERT = cut-rank` is configuration-independent — *static certifiability has no
field side* — the certification-currency mirror of C-80's zero-beyond-contact. A computed no-response
theorem is a result, not a failure.

---

## 4. WHAT THIS DOES NOT SETTLE (carried verbatim into T-17's ledger of limits)

1. **MEDIATION** — every field instrument here is mediated; whose field it is, C-77 must absorb; no
   control can decide it.
2. **COMPOSITION** — both energy-channel instruments read **VALUE**, not accumulating **AMOUNT**: the
   Second Lump's response is Z₂-valued and may saturate mod 2; the Clock's screens. The anchor's
   source standard is sign-definite accumulation. **If the field side cannot compose, the field side
   and the source side do not meet the way emergence at scale needs.** The composition computations
   are scheduled first for exactly this reason.
3. **SHAPE** — all shape-at-scale claims wait on the unbuilt signed connected kernel; below `μ_c` the
   reading is exponentially confined, and whether the surface sits at criticality is O-58's bit.
4. **O-39's FULL ASK** — metric, connection, deflection remain unconnected to anything computed;
   nothing dynamical moves; no tidal concept has been earned.
5. **EMPIRICAL CONTACT: ZERO** (X-4) — every venue here is a model venue; the PROOF bar (two
   structurally different real record surfaces) has not been approached by any of this.

Also carried: the fourth design angle (geometry response — is the earned metric itself responsive to
content, or provably rigid?) never reached the judge and is owed; O-53 remains uncomputed regardless
of which lane fires.
