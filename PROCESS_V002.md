# THE PROCESS — as it stands, 2026-08-20

**Supersedes the 2026-08-19 statement of this document, `PROCESS_DESCRIPTION_V003`, and the process
paragraphs of `THE_CLAIM_V001`.** Every sentence cites a ledger row; rows cited in the corner
sections are corner-tier (`FORMAL`) results and say nothing about the world by themselves —
`CORE_FRAMEWORK_V001` states what `PROVED` means and the tool enforces it. **No withdrawn claim
appears here.** Terms are as mapped in `GLOSSARY.md`; the model is the **Universal Record Model
(URM)**, `model/project_model.py` (`C-85`).

---

## 1. WHAT A RECORD IS — the amended definition, grounded

A record lives on **carrier ⊗ local bath** and satisfies the amended clauses **(i′)–(v′)** with
the declared tolerance tuple the anchor enumerates, `(T, t_m, δ, t_w, Δ, W, N_r)`: a commuting
spectral family `(i′)`, **durable** — read throughout this program as `|λ_record| ≤ 1/t_m` on the
record's own Liouvillian mode (`C-75`, `C-86`) — `(ii′)`, non-constant on the energy shell `(iii′)`,
**writable by an energy-conserving dilation unitary** at free-energy floor `kT ln2 + dE_config`
`(iv′)`, and **protected** with every single-region flip costing at least `E_b`, `E_b/kT ≫ 1` `(v′)`
— `O-51`, **adopted** per the register's O-51 entry on the principal's directive. The amendment was
forced by measurement: the census set's **thirteen** real records were tested against the exact
clauses and **zero** passed — ~20 across the four adversarially checked censuses — and all pass
amended (`O-51`; `H-3` stands `PARTIAL`, converted to a standing acceptance criterion by the goal's
grounding clause, with "why THIS definition" still its open question).

**DEF-A — `T → 0`, `t_m → ∞`, write work `W = 0`, `E_b → ∞`, with admissibility `[U,H] = 0` at the
corner (this `W` is the anchor's write work; §6's certification window is a different `W`) — is the
VERIFIED corner of this definition,
not a second definition**: with no dissipation the slow modes at `t_m → ∞` are exactly the commutant
on every carrier tested, 96/1536/24, and `(ii′)` reduces to `(ii)` at rate `0.00e+00`. Durability
requires `|λ| ≤ 1/t_m`, not `|Re λ|` — a rotating mode does not decay and is not durable (`C-75`).

## 2. WHEN A RECORD EXISTS — the corner, carrier-free

`[L,R]=0 ⟺ [L†,R]=0` for Hermitian `R`, so **`R` lies in the commutant of
`alg{I,H,L_k,L_k†}`** (`C-9`). **A record satisfying (i)–(iv) exists iff that commutant contains a
projection non-trivial on some eigenspace of `H` and trace-balanced** — `C-12`, 210/210 exhaustively.
A record exists on `E_λ` iff `P_λ A P_λ` is a **proper** subalgebra (`C-10`); writable iff `Tr R = 0`
(`C-11`), and under DEF-A iff `Tr(P_E R) = 0` on every eigenspace (`C-16`). The naive converse is false — degeneracy alone gives no record; the converse's
failure is history (the failed O-1, recorded in the register) and `C-12` is the criterion that
survived it. **Clause (ii) is a condition on the ENVIRONMENT**: generic single-site noise leaves
the commutant scalar and admits no record at all — `C-17`, three carriers.

## 3. WHAT CARRIES IT

**EM supplies both boundary maps of the chain complex** — `∂₁` is the Gauss law, `∂₂` the plaquette
term (`A-EM`); the record space is the homology of EM's own complex, `Z₁/B₁`, match rate `1.0`
(`G-9`). **The record, its writer, and the minimal coupling that can form it are all EM holonomies on
non-contractible cycles** (`A-EM2`, `A-EM3`), and the object is neither local nor the whole system:
extent `d = L → ∞` with local content exactly zero (`A-EM4`). The EM role's sentence in the
unified roles statement is `T-15`, BLOCKED on the carrier audit.

**Γ's candidate realisation is `(H₁, ⟨·,·⟩)`, defined only by the function it performs — the
CANDIDATE scoping is the row's own** (`G-5`); homology is not necessary —
records exist in dimension 6 with no qubits and no cells (`C-22`) — and the necessary condition is
`C-12` and nothing narrower. **The world protects records by BARRIERS**: all thirteen census records are protected by
metastable barriers, `E_b/kT` from 50 to 780 (`C-66`); the corpus's kinematic protection is the
`E_b → ∞` limit (`O-51`); and on the torus the homology realisation is exact — 1480 regions,
contractible regions carrying a logical: zero (`C-74`). **Clause (v) also holds on the non-abelian
carrier, by an exact integer obstruction with a genuinely non-abelian mechanism** — `D(D₄)` on the
1×2 torus, `t_E = {30,−48,12,0,6}` for `T = r²` central, with the abelian control failing on the same
lattice; scoped as the row scopes it: a boundary-constrained proxy at `Ly = 1`, the stated mechanism
being what a larger lattice must be tested against (`C-89`).

## 4. HOW A RECORD FORMS

> **A coupling opens a channel to a record iff its compression onto the code space has a non-zero
> component along that record** — `G-16`, 31/31 cycles, three carriers.

The coupling makes the environment's evolution conditional on the record; the divergence of the two
conditional environment states **is** the copy, and the record's value never changes (`F-20`). Two
requirements, both necessary: weight `≥ d` (`F-13`, `F-16`), not commuting with the writer (`F-15`);
an ordinary thermal environment at strong coupling suffices (`F-13`); a single local term gives
exactly zero at every order while a sum forms at order `λ^{2d}`, slope `3.8988` against the
predicted `2d = 4` (`F-23`). Gauge invariance plus a locality bound below `d`
forbids formation exactly (`C-18`). Records compose (`C-19`).

**In the world this becomes the two PROVED encoding-level laws.** Formation creates an accumulated
quantity the unwritten surface lacks — extensive, additive, not a count — and **the write mechanism
fixes its sign**: occupancy-encoded surfaces (NAND net trapped charge) accumulate one-signed for
every data pattern; orientation-encoded surfaces (CoCrPt remanent moment) accumulate under
DC-saturation and screen on real data, measured ratio `0.00096` at `N = 1e5` — `C-71`, `C-72`,
`RECORDS VERIFIED` on two structurally different mechanisms, standard device physics owned by name.

## 5. HOW LONG IT LASTS — the laws and the count

**The lifetime law**: `τ = exp(E_b/kT)/(2 f₀ cosh(dE/2kT))`, the record's own Liouvillian mode read
back by the model and confirmed by direct time evolution at `8.4e-15` — the CoCrPt grain is a record
under the amended clauses and is not one under the exact clauses (`C-69`). **The steady-state law**:
`⟨R⟩_ss = tanh(dE/2kT)`, an exact identity of any two-state detailed-balance GKSL record (`C-70`).
External measured anchors are pinned with sources in `LANE_T41_EXTERNAL/CITATIONS.md` (`T-41`).

**The count law**: clause `(ii′)` applied to the record's own mode yields both widths with zero
adjustable content — `δ_pop(t_m) = kT·ln(expm1(B/kT − ln f₀t_m))` and `δ_coh = ħ/t_m`, two corners of
one modulus — and `k(t_m)` is a dated staircase with drop times
`t*_i = f₀⁻¹ exp((B_i−dE_i)/kT)/(1+e^{−dE_i/kT})`; C-14's clustered count survives only as a corner
proxy, exact at `dE = 0`, broken by chain merging (`C-86`, `C-76`).

## 6. WHY IT IS OF THE PAST — and what the outside can know

The environment holds Holevo information about the record and acquiring it requires coupling over a
region at least `d` — `χ = 0.00000000` for all 24 weight-1 observables (`F-17`); a weight-1 coupling
still entangles, `I(S:B) = 0.045`, transferring zero bits about the record (`F-18`); no system-only
operation removes it, invariance `3.686e-14` (`F-19`); redundancy evens out exactly when fragments
are equivalent (`F-28`).

**Certifiability is the interface**: externally certifiable content equals the stabiliser cut-rank
identically on every region, `CERT = 8s − 10` at the corner (`C-81`); the world tier splits exactly —
the boundary bounds the certification **rate** `6n² − 12n + 8` per epoch, an exact min-cut, while
the whole volume certifies cumulatively at `E_min ≤ 3` epochs on the computed range — the general
bound is `⌈n³/(6n²−12n+8)⌉ ~ n/6`, growing without bound (`C-82`). **The surface owns the clock**: a
record must be certified within its own lifetime, `W = τ/t_epoch`, and
`CERT_W(n) = min(n³, W·(6n²−12n+8))` — certifiable content is boundary-bounded at scale for any
finite-lifetime surface, and `W = ∞` is exactly the DEF-A immortal-record corner (`C-83`).

## 7. WHAT IT COSTS

Alpha is a sum of local terms and splits the record space at order `n*` — the fewest available local
terms whose product reaches the code space (`A-AL`, `C-8`, three carriers). `n* = d` only when the
perturbation class is unrestricted; the two separate on `[[5,1,3]]`, `n* = 5` against `d = 3`
(`T-14`'s restatement). In the record laws the cost appears as `E_b`: whether **alpha's value** enters
through the barrier is `T-38`, open, awaiting its derivation.

## 8. THE EMERGENCE CLAIM — shared origin, increments computed

The claim's shape is fixed: **shared origin, boundary-shaping terms, each classical concept earning
its place, with record-level measures permitted to behave unlike their analogues elsewhere**
(`C-77`, `A-GR3`, `A-GR4`, `D-24`). The increments to date: **distance is earned** — `d_W` is
boundary-crossing cost, a true metric, exhaustively closed (`C-78`); **the boundary law** —
`deg IR = deg C − 1` in both tiers, mechanism exact Euler bulk cancellation (`C-79`); **the division
of labor** — no pure-Γ falloff exists; Γ is the boundary law, the coupling carries the falloff
(`C-80`); **the reachable classes** — a Γ-priced coupling reaches exactly three: exponential below
`μ_c`, the earned dimension's own point-source potential at `μ_c` (the critical resolvent identity IS
the venue's discrete Poisson equation), divergent above (`C-87`); **the Newton composition profile
exists at exactly one computed design point**, transparent admissibility plus occupancy encoding,
36/36 gates (`C-88`); **the D=3 critical member is COMPUTED as `1/d`** — exponent bracket containing
1 on all rays, `G` finite, coefficient tightening onto `3/(2π)`, adversarially confirmed (`C-90`).
**The Newton verdict is MATCHES AT MEMBER LEVEL** — the pre-registered upgrade rule fired on its own
stated condition, superseding NEITHER, with the conditional structure load-bearing and kept: unique
critical member computed; occupancy of criticality **not** earned (masslessness, the
measure-conservation hypothesis `μ = 1/deg = μ_c`); design-point enforcement **not** earned; norm
selection open — never a claim that gravity is derived (`C-77`, `O-58`); also open: the universal
coefficient (`O-55`), the thermodynamic pair (`O-56`), the collapse analogue (`O-57`).

## 9. WHAT IS OURS TO FALSIFY

The distinguishing audit (successor to the failed PF-6 search, recorded in the register) returned **non-empty**: three survivors, each attacked
by three rival advocates under the fairness rule that a rival wins only if its account already made
the statement — **S-1** the durability-clocked certification window `CERT_W = min(volume, W·min-cut)`;
**S-2** flash alpha-channel discrimination; **S-3** the thickness profile `1,2,1,0` (`C-84`). The
count law `k(t_m)` is the wholly-owned falsifiable count law, ownership stated with fragments
conceded by name; its external-data run is the pending step to promotion (`C-86`).

---

## WHAT IS **NOT** ESTABLISHED

- **The Newton verdict is MATCHES at member level only, on conditions not yet earned.** What earns
  criticality (masslessness), whether the surface enforces the design point, and which norm the
  record's physics selects are all open (`O-58`, `T-48`); nothing here asserts classical gravity has
  been recovered.
- **Alpha's value has not been connected to the record laws.** The route through `E_b` is stated as
  a task, not a result (`T-38`; `H-2`'s closure is scoped to "alpha names a generic coupling" at the
  corner).
- **The count law has not met external census data** (`C-86`, promotion pending), and **S-1's
  literature check is owed** (`C-84`).
- **Carrier breadth is being marked, not presumed**: the per-row two-carrier audit is in progress and
  rows it marks SINGLE-CARRIER carry that mark into the proof (`T-9`).
- **The field side of the emergence claim has no instrument yet**: every measured quantity is a
  source; what would measure the field side at the world tier is open (`O-39`, `O-53`).
- **Objectivity needs equivalent environment fragments and the account does not supply them** (`O-21`).
