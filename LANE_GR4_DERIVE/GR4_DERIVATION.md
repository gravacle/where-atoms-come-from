# LANE_GR4_DERIVE — DERIVE THE DEFINITION, DO NOT STIPULATE IT (H-3 head-on)

Companion numbers: `gr4_numbers.py` / `gr4_numbers_output.txt`. Every figure below is
computed there or carries an inline citation.

---

## 1. WHAT A RECORD MUST DO — stated without assuming an answer

A record is a physical configuration **now** whose value licenses an inference about an event
**earlier**. Precisely: there is a variable `X` (the event, at time 0), an observable `R` of some
subsystem (the record, read at time `t_m`), a fidelity `1−δ`, and:

  **(REQ)**  `I(X : R(t_m)) ≥ (1−δ) H(X)` — the mutual information between the event and the
  present observable survives the intervening dynamics — **and** `R(t_m)` is accessible to
  observers whose measurement does not destroy the correlation for the next observer.

This is deliberately close to what the decoherence literature already formalises: quantum
Darwinism states the accessibility half (Zurek, Nat. Phys. 5, 181 (2009)); decoherent
histories states the correlation half (Gell-Mann & Hartle, PRD 47, 3345 (1993), "records"
`R_α` perfectly correlated with history branches). Nothing in (REQ) presumes operator
algebra, degeneracy, admissibility, or topology.

### What follows NECESSARILY from (REQ)

Write the storage interval's dynamics as an open system and coarse-grain nothing yet. Then:

- **N1 — DISTINGUISHABILITY.** `H(X) ≥ 1` bit requires ≥ 2 distinguishable sectors of the
  record's state space: a decomposition `P₊ ⊕ P₋` (or finer) with near-orthogonal
  coarse-grained states. *This is clause (i) in kind* — but only on the record manifold, not
  on the full Hilbert space (see §3-i).
- **N2 — A WRITE CHANNEL EXISTED.** The correlation had to be created: at some earlier time
  the dynamics mapped `X`-differences into `R`-differences. *Clause (iv) in channel form.*
- **N3 — PERSISTENCE.** The flip/leak rate `Γ` must satisfy `Γ·t_m ≲ δ`. *Clause (ii) in
  lifetime form* — exactly O-5's relaxation, with tolerance `τ ≥ t_m/δ`.
- **N4 — PROTECTION FROM TYPICAL NOISE.** N3 must hold not only for the modelled `{L_k}` but
  for the environment actually present; the events able to rewrite `R` must be atypical at
  temperature `T` over `t_m`. *Clause (v) in rate form.*
- **N5 — NON-REDUNDANT WITH THE BACKGROUND.** If `R`'s value is computable from what the
  observer already knows (the macrostate: energy, conserved charges, the Hamiltonian), then
  `I(X:R)` adds nothing. *Clause (iii) in informational form.*
- **N6 — READABILITY / REDUNDANCY.** (REQ)'s second half. The value must be copied into many
  fragments so that reading is non-destructive and many observers agree. **No clause of the
  five states this.** It is quantum Darwinism's entire content, and for real records the
  numbers are enormous (a 1 mm² ink mark in sunlight broadcasts ~1e15 photon-copies/s;
  redundancy ~1e8 for a micron grain after ~1 μs — Riedel & Zurek, PRL 105, 020404 (2010)).
- **N7 — AN ARROW.** Inferring the past from the present at all requires the low-entropy
  initial condition (the past hypothesis); records point backward because entropy was lower.
  The clause set is silent; decoherent histories puts this in `ρ`. Noted, not chargeable to
  the five clauses alone — but it is where the program's third term may actually live (§6).

**And one thing that does NOT follow:** that `R` be a *fixed operator*. (REQ) constrains the
*information*, and the world's longest-lived records survive only by migrating substrate
(DNA replicated ~every generation; texts recopied; digital data scrubbed and rewritten).
Over those histories **no** operator approximately commutes with the dynamics; what persists
is a decodable content — the QEC notion (information survives iff a decoder exists), not the
commutant notion. The five clauses hard-code the fixed-operator special case.

---

## 2. THE EXISTING FORMALISATIONS, CLAUSE BY CLAUSE

| formalisation | its definition of a record | vs the five clauses |
|---|---|---|
| **Quantum Darwinism** (Zurek 2003/2009; Blume-Kohout & Zurek PRA 73, 062310 (2006)) | a pointer observable of S whose value is imprinted redundantly: `I(S:F) ≥ (1−δ)H(S)` for many disjoint environment fragments `F`; redundancy `R_δ = N/N_δ` | Approximate by construction; state- and dynamics-dependent; **redundancy central — absent from the five clauses**. Has no writability or protection clause (formation is the decoherence process itself). The program **adds** (iv)/(v) — a genuine addition — and **subtracts** redundancy — a genuine omission. |
| **Einselection / pointer states** (Zurek RMP 75, 715 (2003)) | pointer observables minimise entropy production; exact limit `[Λ, H_int] = 0` | Clause (ii) **is** exact einselection extended to `[H,R]=0` and all `L_k`. Known, in its zero-tolerance limit. The world uses the sieve (approximate), never the exact limit. |
| **Decoherent histories** (Gell-Mann & Hartle PRD 47, 3345 (1993); Halliwell PRD 60, 105031 (1999)) | projections `R_α` with `Tr(R_α C_β ρ C†_γ) ≈ δ_αβ δ_αγ p_α` — records correlated with branches, **relative to ρ** | Their records are **state-dependent**; the five clauses quantify over no state. State-independence is exactly what forces exact degeneracy on the program (only algebra can hold a bit if no ρ is special). Halliwell's finding that the natural quasiclassical records are **local densities of conserved quantities** (slow hydrodynamic modes) anticipates §4's amendment: slowness, not commutation. |
| **QEC logical operator** | `Z̄` commutes with the stabiliser group, is not in it; `X̄` anticommutes with `Z̄`; distance = smallest support reaching the code | The five clauses are, almost verbatim, "**R is a logical Pauli of a passive (self-correcting) quantum memory** whose stabiliser Hamiltonian is `H` and whose noise respects the code, with (v) = distance exceeds any contractible region." The relevant no-gos then apply: no self-correction in 2D (Bravyi & Terhal, NJP 11, 043029 (2009)); 3D open (Haah's code: barrier ~log L, still finite τ at fixed T); genuine self-correction known only in 4D (Dennis et al., JMP 43, 4452 (2002); Alicki et al. 2010). **At 300 K in ≤3 dimensions, objects satisfying the five clauses with growing protection are not known to exist.** The world's records are instead *classical* self-correcting memories (a 3D ferromagnet is one), which nature has in abundance. |

**Verdict on priority:** clause (ii) = exact einselection (known); clauses (i)+(iv)+(v) =
logical-Pauli pair + distance (known, QEC); clause (iii) is the program's own sharpening
(information beyond the energy) and is a real contribution in its exact form. The genuinely
novel move — making *writability* a defining clause and exhibiting its tension with
durability — is the program's. The genuinely missing move — redundancy/readability — is
Darwinism's, and the definition omits it.

---

## 3. THE FIVE CLAUSES AGAINST THE WORLD'S RECORDS

Test corpus (numbers in `gr4_numbers_output.txt`): **HDD grain** (FePt, 5 nm, KuV = 199 kT,
τ ~ 1e-9·e^199 s, write ~1e12× Landauer at device level, τ/t_w ~ 1e26), **DNA base**
(depurination k = 3e-11 /s at 37 °C → τ ~ 1e3 yr unrepaired, Ea ≈ 50 kT; polymerase write
~24 kT ≈ 35× Landauer), **NAND flash** (3.1 eV = 120 kT barrier; retention set by oxide
defects, not the named barrier), **ink on paper** (chemically inert C; substrate hydrolysis,
centuries–millennia), **K–Ar geochronometer** (written by weak decay, stored under a ~72 kT
Ar-diffusion barrier, read *destructively*), **photograph** (4-atom Ag cluster, developed
×1e9), **baryon number / B−L** (sphaleron barrier 9 TeV = 3.5e14 kT).

| clause (exact, as written) | HDD | DNA | flash | ink | K–Ar | B−L | verdict |
|---|---|---|---|---|---|---|---|
| (i) `R²=I` on the full space | ✗ | ✗ | ✗ | ✗ | ✗ | ✓ | **fails**; holds on a coarse-grained record manifold `P_rec ⊂ H` (the two metastable wells are a sliver of a 2^(10^5)-dim space, and "melted/blank/demagnetised" states are in neither eigenvalue) |
| (ii) `[H,R]=0, [L_k,R]=0` | ✗ | ✗ | ✗ | ✗ | ✗ | **✓ (exact!)** | **fails exactly, holds as O-5's width bound** — with the essential caveat of §4: the small quantity is the *off-block* part of `H, L_k`, while the full commutator norm is O(1) (spin waves roam inside each well) |
| (iii) non-constant on an eigenspace (exact degeneracy) | ✗ (Earth-field split 3.4e-3 kT ≫ many-body level spacing ~e^(−10^4) kT) | ✗ (different bases are different molecules, eV apart) | ✗ | ✗ | ✗ | ≈✓ (p→e⁺π⁰ final states: same energy, different B) | **fails categorically.** No real record stores its bit in an exact eigenspace degeneracy. What real records have instead is **broken ergodicity**: two macrostates each metastable. Exact degeneracy is the T=0, closed-system shadow of that. |
| (iv) admissible writer, `[U,H]=0` (⇒ dE=0, C-60) | ✗ (write = 2 T + 700 K pulse; 1e12× Landauer) | ✗ (24 kT/base, dissipative, ATP-driven) | ✗ (1e7× Landauer) | ✗ | ✗ (writer is the weak interaction) | ✗ (no SM writer at all) | **fails for every record on Earth, and must**: overwriting an unknown bit costs ≥ kT ln 2 (Landauer 1961), so a zero-work writer is thermodynamically forbidden for any record an agent actually uses. O-4/O-44/C-60 already contain this; GR4 confirms it from the world's side. |
| (v) no single contractible region flips it | ✗ (the 5 nm grain IS one contractible region; the write head flips it locally) | ✗ (one UV photon or •OH radical flips a base) | ✗ (one cosmic-ray strike) | ✗ (a match) | ✗ (reheat above closure T) | ✓ (no local SM operator changes B−L) | **fails in the interesting direction: every terrestrial record IS locally flippable.** Real protection is a *rate* statement — the local events able to flip it are exponentially atypical (ν t_m e^(−ΔF/kT) ≤ δ) — an improbability, never an impossibility. |
| (vi — absent) redundancy/readability | needed | needed | needed | needed | needed | needed | the world's records all carry it; the definition never asks |

**Score: of the five exact clauses, the number satisfied by a typical terrestrial record is
zero.** Every one holds in a natural approximate/coarse-grained form — but the approximate
forms of (iii), (iv), (v) are not small deformations of the written clauses; they are
different mechanisms (metastability for degeneracy; paid work for admissibility; large-
deviation suppression for regional impossibility).

**The one exact-clause record found in nature is B−L**: anomaly-free in the SM, so
`[H_SM, B−L] = 0` exactly; non-trivial (same-energy states of different B−L exist);
protected (no local SM process reaches it). And it **fails (iv)**: nothing in the SM writes
it. This is not an accident — it is the theorem the corpus already proved from the other
side (C-60 + Landauer): *within one fixed dynamics, exact durability and writability
coexist only on the measure-zero degenerate locus.* The world's records escape because the
write epoch and the storage epoch have **different effective dynamics** (the field is on,
then off; the universe is at 10^15 K, then 2.7 K; the polymerase is bound, then released).
The five clauses hand one stationary `(H, {L_k})` to both epochs — that single modelling
decision is the root from which the failures of (iii), (iv), (v) all grow.

### The writable/durable tension, resolved the way the world resolves it

The program spent months on (iv)-vs-(v) and resolved it with admissibility + regions. The
world resolves it with **one inequality and a fluctuation theorem**:

- PROTECTED: `ν t_m e^(−ΔF‡/kT) ≤ δ` ⟺ `ΔF‡ ≥ kT ln(ν t_m/δ)` — 40 kT for 10 years at
  ν = 1e9 /s, 45 kT at δ = 1e-2, 59 kT for a Gyr; the disk industry's "60 kT rule" is this
  inequality with margin.
- WRITABLE: a *driven* protocol transiently supplies ΔF‡ (heat pulse, write field, enzyme)
  and dissipates ≥ kT ln 2.
- The bath can only do what the writer does by paying the same work spontaneously, at
  probability e^(−W/kT) (Crooks, PRE 60, 2721 (1999)). **The asymmetry between the writer
  and the noise is a work budget, not an operator support.** The agent has a free-energy
  source; the bath has only kT.

---

## 4. THE AMENDED DEFINITION — the one the world's records satisfy

Fix ambient temperature `T`, mission time `t_m`, failure tolerance `δ`, write time `t_w`.
Let the storage-epoch dynamics be a Lindbladian `L` (generator on density matrices).

**R is a (T, t_m, δ)-RECORD if:**

- **(i′) BIT.** There is a projection `P_rec = P₊ ⊕ P₋` onto a *record manifold* with
  `R = P₊ − P₋` (so `R² = P_rec`, not `I`), and the coarse-grained states on `P₊, P₋` are
  `(1−δ)`-distinguishable.
- **(ii′) DURABLE — the slow-manifold clause.** `P₊, P₋` span a **metastable subspace of
  `L`**: every eigenvalue `λ` of `L` supported on the record manifold has
  `|Re λ| ≤ δ/t_m`, while the manifold is reached from relevant initial states within
  `t_w`. Equivalently (O-5's form): the dressed record's spectral width obeys
  `width ≤ ħ δ/t_m` **on the record manifold** — the tolerance is the inverse lifetime,
  exactly as O-5 registered; GR4 adds *where* the bound lives: the off-block parts
  `P∓(H, L_k)P±` must be small; the within-block parts may be O(1).
  As `t_m → ∞` this contracts to `ker L` = the commutant (Frigerio; Baumgartner–Narnhofer),
  i.e. to clause (ii) as written. **The exact clause is the τ→∞ limit of this one, and the
  metastable-manifold theory needed already exists** (Macieszczak, Guţă, Lesanovsky,
  Garrahan, PRL 116, 240404 (2016)).
- **(iii′) NON-TRIVIAL.** Conditional on the observer's background knowledge (the
  macrostate: energy density, conserved charges, `H` itself), both values retain support:
  the bit is **not computable from the thermodynamic description**. Tolerance: the
  free-energy split `ΔF₊₋` between the two values is unconstrained *except* through (ii′)
  (both lifetimes must clear `t_m/δ`). Exact degeneracy is the special case
  `ΔF₊₋ = 0, T = 0`.
- **(iv′) WRITABLE.** There exists a **physical channel** (time-dependent control and/or
  bath coupling — O-4's untested second disjunct, now the primary text) `W_b` with
  `W_b(P_rec-supported states) → P_b` at fidelity `1−δ` in time `t_w`, with
  `t_w ≤ δ·t_m` (timescale separation; real ratios `t_m/t_w` run 1e13 (DNA) to 1e26 (HDD)),
  at work `W ≥ kT ln 2` (Landauer floor — the definition should *state* the floor, not
  contradict it).
- **(v′) PROTECTED.** No **typical** environmental process at temperature `T` performs
  (iv′) uninstructed: `ν t_m e^(−ΔF‡/kT) ≤ δ`, where `ΔF‡` is the smallest free-energy
  barrier any local process must cross to flip `R`, and `ν` the attempt rate. (iv′)+(v′)
  together say: **the record is flippable for W and unflippable for free** — the
  writer/noise asymmetry is the work budget.
- **(vi′) REDUNDANT.** The value is imprinted in `N_r ≫ 1` disjoint fragments
  (`I(R:F_i) ≥ (1−δ)` bit each), so reading is non-destructive and inter-observer
  consistent. Everyday records: `N_r` from 2 (DNA strands) to 1e15/s (photon environment of
  an ink mark).

Tolerance parameters the definition now carries: `T, t_m, δ, t_w, ΔF‡, ν, N_r` — seven, all
measurable, all with real values in `gr4_numbers_output.txt`. Every exact clause is
recovered in the limit `T → 0, t_m → ∞, W → 0, N_r` dropped.

**Computability.** The model changes by one move: `m.records()` currently returns
projections in the commutant of `alg{I,H,L_k,L_k†}` (C-9/C-12). Amended: return the
eigenmodes of the Lindbladian with `|Re λ| ≤ 1/t_m` — the commutant is exactly the
`1/t_m → 0` limit of this computation, so C-12's machinery deforms continuously instead of
being discarded. This is implementable in the existing `record_model.py` (the generator is
already constructed for `formation()`).

---

## 5. COST ACCOUNTING — row by row under (i′)–(vi′)

Principle: pure mathematics about the exact clauses **survives as mathematics, scoped to the
exact sector**; any row whose *physical* force ran through exactness NARROWS to the τ→∞
limit; any **negative** proved at zero tolerance is UNDECIDED at finite tolerance and must
be re-run (H-13 said precisely this).

| row | status under the amendment | why |
|---|---|---|
| **C-9, C-10, C-11** (commutant lemmas, trace-balance) | **SURVIVE** (math) and remain the `t_m→∞` limit of the amended construction | the commutant = ker L; the lemmas anchor the limit |
| **C-12** (existence ⟺ commutant projection, non-trivial + trace-balanced) | **NARROWS** | still the exact-sector existence theorem; as *the* record-existence criterion it must be re-proved for the slow manifold: candidate statement "a (T,t_m,δ)-record exists iff the Lindbladian has ≥2 eigenvalues with |Re λ| ≤ δ/t_m whose eigenmodes are (1−δ)-distinguishable" — the metastability literature gives the tools; the trace-balance hypothesis (H2) should *disappear* in the amended form (it was the free-writability shadow) |
| **C-60** (admissible ⇒ dE = 0) | **SURVIVES as math; its exclusionary use FALLS** | it is now the proof that DEF-A admissibility is unphysical (Landauer): a definition of writing that forces dE=0 defines writing out of the world. The "entire energy route" it closed **reopens** under (iv′) — writes cost work, records carry energy budgets |
| **C-52** (Z[i] quantisation — no slot for a gravity-strength residual) | **UNDECIDED — must re-run** | proved with trace ratios *exactly* 0/fourth-roots on exact stabiliser records; under (ii′) the algebra acquires corrections of order the commutation defect `~ħ/(t_m·gap)` — and a gravity-strength residual is precisely the size that could hide inside a finite tolerance. The exclusion may return, but it is not inherited |
| **C-53** (symplectic-form determination) | **UNDECIDED — must re-run** | same character: exhaustive over exact Pauli structures; approximate records are not signed Paulis |
| **C-63** (power-law record–record falloff on a gapless mediator) | **SURVIVES, lightly narrowed** | a dynamics computation; its records enter only through their operators' existence, and metastable records reproduce them to `O(δ)`; re-check is cheap |
| **G-5** (H₁ + intersection pairing satisfies R1–R3) | **SURVIVES** (math), scoped to the exact sector | nothing false in it; it just no longer describes the *only* way R1–R3 are met |
| **G-14** (Γ = non-triviality of the homology of EM's complex) | **SURVIVES scoped**; its generality FALLS | under (iii′)/(v′) the record space of a real record is a metastable manifold, not a homology group; homology is one realisation |
| **A-EM** (record IS a holonomy of EM's own complex) | **NARROWS — and should be re-aimed upward** | as registered, scoped to exact-sector carriers. The world-side statement that replaces it is *stronger and true*: every barrier ΔF‡ and every write channel in the test corpus (anisotropy = spin-orbit EM; covalent bonds; oxide barriers) **is electromagnetic**. EM's real role: it supplies the barriers and the channels — all of (ii′), (iv′), (v′) run on EM. The refuter in the ledger ("a record that is not a holonomy") is now instantiated by every record on this page |
| **A-GR** (Γ/topology supplies record space, writer, protection — records *require* genus) | **the necessity direction FALLS; the possibility direction SURVIVES scoped** | under (v′), protection is available from free-energy barriers with **no topology** (a 3D ferromagnet's droplet barrier grows as the cross-sectional *area* — faster than the 2D toric code's O(1) anyon barrier: at 300 K a magnetic domain outlives a 2D toric logical bit by an exponential factor). The fragility dichotomy (W-60/W-61/Thm D) is correct *in the exact game*, where any splitting kills; the amended game tolerates splitting up to ~ΔF‡ (tens of kT ≈ 0.1–1 eV — a tolerance ~1e6× wider than the 1e-6 exact-sector kill threshold, measured against the same perturbation scale). "Records require genus" becomes "**τ=∞ records at fixed size under generic perturbations require genus**" — true, and no longer about the world's records |
| **A-AL** (alpha splits at order n*) | **NARROWS** | intact as perturbation theory on degenerate codespaces; real records have no exact degeneracy to split. Its amended descendant: alpha (generic local-coupling strength, T-22) sets the *attempt-rate prefactor and barrier attack* in (v′) — a re-derivation target, not an inheritance |
| **A-GR2** (Γ supplies the channel; classical gravity is what Γ becomes at scale) | **FALLS with A-GR's necessity direction** | its evidence chain (O-20, intersection-pairing selection of forming couplings) lives in the exact sector; under (iv′) channels are selected by work budgets and barrier geometry, and nothing yet connects those to an intersection pairing |
| **H-13's negatives (route 1 closure via C-26/C-52/C-53 lineage, O-42's exclusion machinery)** | **UNDECIDED — re-run under stated tolerances** | H-13 registered exactly this contingency. A negative proved at tolerance 0 does not survive tolerance e^(−ΔF‡/kT) without a new argument; the corrections in an approximate record algebra sit at the same order a gravity-strength residual would. This is the *good-news* cost: the program's most disappointing negatives are not inherited by the amended definition |
| **F-13, F-16** (formation needs weight ≥ d, non-commuting with writer; Knill–Laflamme suppression) | **SURVIVE in the exact sector; re-derivation target in the amended one** | these are bath-coupling results, the closest thing the corpus has to (iv′) physics already |
| **O-5** | **CONFIRMED and PROMOTED** | its width-bound is the germ of (ii′); GR4's check: it suffices for real records **only** once the bound is placed on the off-block parts on the record manifold (a full-norm bound fails — `‖[H,R]‖` is O(1) for every real record while the record is fine) |
| **O-4 / DEF-A** | **REPLACED** | the untested disjunct ("or a physical channel") becomes the primary text of (iv′); DEF-A survives as the τ→∞, W→0 corner |
| **P-1, P-2, P-3 / C-14 count law** | **NARROW to the exact sector** | the amended count of independent bits is the count of metastable phases (number of slow Lindblad modes), not `v₂(m_E)`; on real carriers these differ (a ferromagnet has 2 phases and no exact degeneracy at all) |

**Net:** the theorems all survive as theorems; what falls is the bridge sentence "and
therefore this is what a record is." The exact sector the corpus proved things about is the
`T→0, t_m→∞, W→0` corner of the amended definition — a corner containing B−L and
engineered quantum memories, and not containing one terrestrial record.

---

## 6. WHAT THE PROGRAM'S THREE TERMS LOOK LIKE FROM THE WORLD'S SIDE

Stated as re-derivation targets, not results:

- **EM**: supplies every barrier and every channel in the test corpus. (ii′), (iv′), (v′)
  are electromagnetic physics at 300 K. This is A-EM's honest successor and is *larger*
  than "supplies the carrier."
- **ALPHA** (generic coupling strength, per T-22): sets the attempt frequency ν and the
  perturbative attack on the barrier — the rate side of (v′).
- **GRAVITY**: absent from every lab-scale clause — and present at the root of N7: the
  low-entropy past that makes inference-from-records possible at all, and the free-energy
  sources (gravitational condensation, stars) that pay every Landauer bill in (iv′). If the
  program wants gravity in the story of the world's records, this is the door that is
  actually open. (Flagged: motivational direction, not a registered claim.)

And the title's own record: **that atoms exist is a record (baryon asymmetry), written when
the effective dynamics permitted it and protected now by a topological barrier of the gauge
vacuum (sphaleron, 9 TeV = 3.5e14 kT).** The program's topological instinct is not wrong
about nature's deepest records; it is wrong as a *requirement* on records in general.

---

## 7. NUMBERS I AM UNSURE OF (bounds given)

- Attempt frequencies ν: quoted 1e9–1e13 /s; conclusions move by < 25% in barrier units
  (logarithmic).
- FePt Ku: literature spread 4.6–10 ×10^6 J/m³; 199 kT could be 140–300 kT.
- Device-level HDD write energy 3e-9 J/bit: drive-power/throughput estimate; true media
  dissipation per bit is smaller, bounded below by ~KuV ~ 2e-19 J. The "×Landauer" factor
  is 1e2 (physics floor) to 1e12 (device).
- DNA write dissipation 24 kT/base: 2-ATP estimate; with proofreading, 10–100 kT
  (Bennett's analyses).
- Riedel–Zurek redundancy 1e8: order of magnitude from their model, not re-derived here.
- Ar diffusion Ea ~180 kJ/mol: mineral-dependent (150–250 kJ/mol range).
- Flash: thermal barrier is not the operative retention channel (oxide defects are); the
  120 kT figure is an upper story, stated as such.
- Sphaleron 9 TeV: Klinkhamer–Manton scale; lattice value ~9.1 TeV, ±few %.

## 8. NEXT STEP (no route closes without one)

Implement the one-move amendment in the model: `m.records(t_m)` = slow eigenmodes of the
already-constructed generator with `|Re λ| ≤ 1/t_m`; verify it returns the exact records at
`t_m → ∞` on the three registered carriers, then run C-52's trace-ratio computation on a
metastable (non-exact) record to decide whether the gravity-exclusion survives finite
tolerance. That single computation adjudicates the largest UNDECIDED above.
