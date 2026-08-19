# THE PROCESS — as it stands, 2026-08-19

**Supersedes `PROCESS_DESCRIPTION_V003` and the process paragraphs of `THE_CLAIM_V001`, both of which
predate the O-, F- and G-series.** Every sentence cites a ledger row. **No withdrawn claim appears
here** — the withdrawn ones are in `REGISTER_V001.md` and `THE_ROAD_V001.md`, kept so they are not
re-made. Terms are as mapped in `GLOSSARY.md`: **BORROWED**, **RENAMING**, or **OURS**.

---

## 1. WHAT A RECORD IS

A record is an observable `R` on a system `(H, {L_k})` that is a **bit** `(i)`, **durable** `(ii)`,
**non-trivial** `(iii)`, **writable** `(iv)` and **protected** `(v)` — `CORE_FRAMEWORK_V001`.

**`admissible` means `[U,H] = 0`** — `O-4`, and **provisional**: tested on one carrier family with no
fallback. **The word appears in BOTH (iv) and (v)**, and clause (v) is **false of the toric code**
without it — a single-edge `X` flips the record, `1158` such operators at `L=3` against `0` admissible
ones (`O-12`).

## 2. WHEN A RECORD EXISTS — **carrier-free**

Clause (ii) is stronger than it reads: `[L,R]=0 ⟺ [L†,R]=0` for Hermitian `R`, so **`R` lies in the
commutant of the `*`-ALGEBRA** `alg{I,H,L_k,L_k†}` (`C-9`).

> **A record satisfying (i)–(iv) exists iff that commutant contains a projection that is non-trivial
> on some eigenspace of `H` and trace-balanced** — `C-12`, verified 210/210 exhaustively.

A record exists on `E_λ` iff `P_λ A P_λ` is a **proper** subalgebra (`C-10`); it is writable iff
`Tr(P_E R) = 0` on every eigenspace (`C-11`), from which **no record is writable in odd dimension**.
The **naive converse is FALSE**: degeneracy alone does not give a record (`O-1`).

**Clause (ii) is a condition on the ENVIRONMENT, not on the record.** Generic single-site noise leaves
the commutant scalar and admits **no record at all** — `C-17`, three carriers.

## 3. WHAT CARRIES IT

**EM supplies both boundary maps of the chain complex** — `∂₁` **is** the Gauss law, `∂₂` **is** the
plaquette term (`A-EM`). **The record space is the homology of EM's own complex**, `Z₁/B₁`, match rate
`1.0`, with numerator and denominator both growing with area while the quotient does not (`G-9`).

**The record is an EM holonomy on a non-contractible cycle; so is its writer; so is the minimal
coupling that can form it** — three separate measurements, `8 of 8` for the third (`A-EM3`, `A-EM2`).
The object is **neither local nor the whole system**: extent `d = L → ∞`, share of EM's content
`2/(L²+1) → 0`, and **local content exactly zero** (`A-EM4`, `D-5`).

**Γ is `(H₁, ⟨·,·⟩)`** — defined **only** by the function it performs, R1–R3 (`THIRD_TERM_V001`,
`G-5`). It is **a requirement on EM's complex, not an ingredient beside it** (`G-14`). **`H₁` is not
unique**: good qLDPC codes exceed Delfosse's compact-2-manifold ceiling by `Θ(n²/log²n)` (`G-6`), and
**homology is not even necessary** — records exist in dimension 6 with no qubits, no stabiliser group
and no cells (`C-22`). **The necessary condition is C-12 and nothing narrower.**

The **count** is an index, `dim H₀ + dim H₂ − χ`, holding on 11 carriers including one with no
geometric realisation (`G-7`); from the clauses alone it is `min_E v₂(m_E)` **when the noise is empty**
(`C-14`, and false with noise present).

## 4. HOW A RECORD FORMS

> **A coupling opens a channel to a record iff its compression onto the code space has a non-zero
> component along that record** — `G-16`, 31 of 31 cycles, three carriers.

**The mechanism:** the coupling makes the environment's evolution **conditional** on the record, so the
two record values drive the environment to **different states**, and that divergence **is** the copy
(`F-20`). **The record's value never changes** — it is being **read**, not written.

**Two requirements, both necessary:** a coupling of **weight ≥ `d`** (`F-13`, `F-16`) that **does not
commute with the writer** (`F-15`). **Nothing about equilibrium** — an ordinary thermal environment at
strong coupling suffices (`O-11`). A **single** local term gives **exactly zero at every order**; a
**sum** of them forms a record at order `λ^{2n*}` (`F-23`).

**Gauge invariance plus a locality bound below `d` forbids formation exactly** (`C-18`, three
carriers). What distinguishes couplings is **homological, not gauge-theoretic** (`G-16`).

**Records compose**: each forms with `χ = 0.908` while the other stays at `0.00000000` and its value
moves by `~1e-16` (`C-19`, two carriers).

## 5. WHY IT IS OF THE PAST

**The environment holds Holevo information about the record**, and **acquiring it requires coupling
over a region at least as large as `d`** — `χ` is `0.00000000` for all 24 weight-1 observables and
non-zero at weight `d` (`F-17`). **This is not ambient decoherence**: a weight-1 coupling still
entangles environment with system, `I(S:B) = 0.045`, while transferring **zero bits about the record**
(`F-18`).

**No system-only operation removes it.** `I(S:B)` is invariant under every system-only unitary,
`3.686e-14`; the system can **move** which observable the environment knows about, not destroy the
correlation (`F-19`). **The arrow is relative, not absolute** — the joint evolution is exactly
reversible (`F-20`).

**Fragments hold it redundantly**, and **redundancy evens out exactly when the fragments are
equivalent** — spread `0.0000` for identical energies, `1.247→1.391` relative spread for random ones
(`F-28`).

## 6. WHAT IT COSTS

**Alpha is a sum of local terms and splits the record space at order `n*`** — the fewest **available**
local terms whose product reaches the code space (`A-AL`, `C-8`, three carriers). **`n* = d` only when
the perturbation class is unrestricted.** Measured: generic perturbations give `2.0211 / 3.0036 /
3.0046`, tracking `d`; record-commuting ones give `2.0389 / 4.9836 / 3.0118`, tracking `n*` — and the
two **separate on `[[5,1,3]]`, where `n* = 5` against `d = 3`**. **`d` is the special case, not the
law.**

---

## WHAT IS **NOT** ESTABLISHED

- **No distinguishing prediction, and no known route to one.** Six rival accounts were examined and
  every measured result is already theirs (`PF-6`). `X-4`, `T-VI.3`, `T-VI.4` all **BLOCKED**.
- **Γ has not been shown to be where gravity comes from.** Γ's record-level content is built and
  gravity was kept out of its definition; **the recovery of classical gravity from accumulated Γ is
  untouched** (`H-1`, plan `T-20`).
- **Alpha has no content beyond a generic coupling** (`H-2`).
- **The five clauses are written, not derived** (`H-3`).
- **Nothing connects these carriers to macroscopic scale** (`H-6`).
- **Objectivity needs equivalent environment fragments and the account does not supply them** (`O-21`).
