# PROOF — RECORD FORMATION AND THE ROLES OF EM, GRAVITY AND ALPHA — V001 — 2026-08-18

Numerical verification of every statement: `LANE_W63_DERIVATION/`, `LANE_W62_CREATION/`,
`LANE_W61_TOPOLOGICAL/`, `LANE_W60_WHICH_HAMILTONIANS/`.

## 0. SETUP

Let `Σ` be a closed orientable surface with a cell decomposition: vertices `V`, edges `E`,
faces `F`, genus `g`. A `Z_2` gauge field assigns `s_e ∈ Z_2` to each edge. Write
`C_1 = Z_2^E` for 1-chains, `∂: C_2 → C_1` the boundary map, `Z_1 = ker ∂_1` the cycles,
`B_1 = im ∂_2` the boundaries, `H_1 = Z_1/B_1`. Cochains, `Z^1`, `B^1`, `H^1` dually.

**The physical Hilbert space** is `span{ |s⟩ : s ∈ Z_1 }` — the Gauss law is `∂_1 s = 0`.
Define `M(z)|s⟩ = |s+z⟩` for `z ∈ Z_1`, and `Z(c)|s⟩ = (−1)^⟨c,s⟩|s⟩` for `c ∈ C^1`.
The Hamiltonian is `H = −Σ_{p∈F} M(∂p)`.

---

## THEOREM A — THE RECORD SPACE HAS DIMENSION `|H_1| = 2^{2g}`

**Proof.** Ground states satisfy `M(∂p)|ψ⟩ = |ψ⟩` for every face `p`, so the amplitude of `|ψ⟩` is
constant on each coset of `B_1` in `Z_1`. Conversely each coset gives one such state. Hence
`dim(ground space) = |Z_1/B_1| = |H_1(Σ;Z_2)| = 2^{2g}`. ∎

**Verified:** `2×2` and `2×3` tori — derived `2^{dim Z_1 − dim B_1} = 4`, measured `4`.
**Disk:** `H_1 = 0`, so dimension `1` — **no record**. Measured: `0` non-contractible cycles.

---

## THEOREM B — THE LOGICALS ARE `H_1` AND `H^1`, AND A WRITER ALWAYS EXISTS

**(i)** `Z(c) M(z) = (−1)^{⟨c,z⟩} M(z) Z(c)`.
*Proof.* `Z(c)M(z)|s⟩ = (−1)^{⟨c,s+z⟩}|s+z⟩` and `M(z)Z(c)|s⟩ = (−1)^{⟨c,s⟩}|s+z⟩`. ∎

**(ii)** `M(z)` preserves the ground space iff `z ∈ Z_1`, and acts as the identity on it iff `z ∈ B_1`.
`Z(c)` preserves it iff `c ∈ Z^1`, and acts as the identity iff `c ∈ B^1`.
*Proof.* Immediate from Theorem A's coset description and from `Z(c)M(∂p) = (−1)^{⟨c,∂p⟩}M(∂p)Z(c)`. ∎

**(iii)** **A WRITER EXISTS FOR EVERY RECORD.** By Poincaré duality the intersection pairing
`H^1 × H_1 → Z_2` is **non-degenerate**. So for every non-trivial `[z] ∈ H_1` there is `[c] ∈ H^1`
with `⟨c,z⟩ = 1`, and by (i) they **anticommute**: `Z(c)` labels the record and `M(z)` flips it. ∎

**Verified:** `‖{M,Z}‖ = 1.26e-15` (anticommute), `‖[M,Z]‖ = 4.00`, record labels `[−1,−1,+1,+1]`.

---

## THEOREM C — NO OPERATOR SUPPORTED ON A CONTRACTIBLE REGION ACTS

Let `O` be supported on `S ⊆ E` containing no non-contractible cycle, and `P` the ground-space
projector. Then `P O P ∝ P`.

**Proof.** Expand `O` in the basis `Z(c)M(z)` with `supp(c), supp(z) ⊆ S`.
If `c ∉ Z^1` then `Z(c)` anticommutes with some `M(∂p)`, so it maps the ground space into an
orthogonal eigenspace and `P Z(c) M(z) P = 0`.
If `c ∈ Z^1` with support in `S`, then since `S` carries no non-contractible cycle, `c ∈ B^1`, so
`Z(c)` acts as the identity by B(ii). Likewise `z ∈ Z_1` supported in `S` forces `z ∈ B_1`.
Every surviving term is therefore proportional to `P`. ∎

**Corollary — durability.** No local noise process can change the record.
**Verified:** worst single-link commutator with the record label `1.59e-16`; contractible loops
`3.01e-17`.

---

## THEOREM D — THE SPLITTING IS `O(ε^d)`, `d` THE MINIMAL NON-CONTRACTIBLE CYCLE

Let `H(ε) = H + ε V` with `V` a sum of local terms. Then the ground-space splitting is `O(ε^d)`,
where `d = min{ |z| : z ∈ Z_1 \ B_1 }`.

**Proof.** In degenerate perturbation theory the ground space is split first at the order `k` at
which some product of `k` terms of `V` acts non-trivially on it. Such a product is supported on the
union of `k` local supports; by Theorem C it acts non-trivially only if that union contains a
non-contractible cycle, which needs `k ≥ d`. Hence no splitting occurs below order `d`. ∎

**Verified — and this is the prediction that could have failed.** Derived `d = 2` for the `2×2`
torus. Measured local slope `d ln(splitting)/d ln ε` = **`2.000, 2.000, 2.000, 2.000`** across two
decades. At `ε = 1e-06` the prediction `ε² = 1e-12` matches W-61's measured `4.9e-13`.
**Symmetry-induced degeneracy has no such protection (`d = 1`): measured `2.0e-06`, linear.**

---

## COROLLARY — THE THREE ROLES

For a `Z_2` gauge theory on a closed orientable surface of genus `g`:

| term | role | by |
|---|---|---|
| **EM** | supplies the field, and the holonomy `Z(c)` **is** the record; its labels are `±1` | setup, B(ii) |
| **GRAVITY (genus)** | supplies the record space, `dim = 2^{2g}`; supplies the **writer**, via non-degeneracy of the intersection pairing; supplies the **protection**, order `d` | **A, B(iii), C, D** |
| **ALPHA** | the electric term is a sum of local operators, so by D it splits the record space at order `d` and destroys it | **D**, and W-60/W-62 |

> **At `g = 0` there is no record: `H_1 = 0`, `dim = 1`, no logicals, nothing to write.
> At `g ≥ 1` the record exists, is writable by exactly the non-local operators, is protected against
> everything local, and is destroyed by alpha at order `d`.**

**This resolves the program's oldest obstruction.** W-29/W-30 proved writable and durable conjugate.
Theorems C and D show why: the conjugacy holds for **local** operations. Theorem B(iii) supplies a
**non-local** writer, which Theorem C shows noise cannot mimic. **Writable and durable simultaneously.**

---

## ATTRIBUTION AND SCOPE — WHAT IS AND IS NOT NEW

**Theorems A–D are standard.** They are Kitaev's toric code: ground-space degeneracy `4^g`, logical
operators as homology classes, code distance, and perturbative stability. **This program did not
discover them and does not claim to.**

**What is this program's:** the chain establishing that record formation *requires* them.
W-60 proves a record exists **iff** `H` has a degenerate eigenvalue; W-61 shows symmetry-sourced
degeneracy is destroyed at `1e-06` while topology-sourced degeneracy survives by four million times;
so **the degeneracy records need must be topological**, and Theorems A–D then say what follows.

**NOT PROVED HERE:** that the relevant topology is spacetime's rather than a lattice's; that it
responds to matter content; the outcome problem (which unravelling is physical); and there is **no
empirical contact anywhere in this program.**
