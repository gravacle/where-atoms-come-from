# THE MODEL — records from `(H, {L_k})` and nothing else

`model/record_model.py`

```python
from record_model import RecordModel
m = RecordModel(H, Ls)          # a Hamiltonian and Lindblad operators. THAT IS THE ENTIRE INPUT.
recs = m.records()              # every record the pair admits
fam, commuting, writable = m.independence(recs)   # the multi-record structure
```

**No lattice. No gauge group. No temperature. No coupling constant. No code. No geometry.**
If a value is not derivable from the pair, the model does not have it and does not invent it.

## THE CONSTRUCTION — each step is a registered theorem

| step | what it computes | theorem |
|---|---|---|
| `star_algebra` | `A = alg{I, H, L_k, L_k†}` | **C-9** — `[L,R]=0 ⟺ [L†,R]=0` for Hermitian `R`, so clause (ii) lands in the `*`-algebra's commutant, not the set's |
| `commutant` | `A'`, by nullspace of the commutator system | **C-9** |
| `minimal_projections` | eigenprojections of a **generic** Hermitian element of `A'` | **C-10** — a record exists on `E` ⟺ `P_E A P_E` is a **proper** subalgebra |
| `clause_iii` | non-trivial on some eigenspace | anchor (iii) |
| `clause_iv` | `Tr(P_E R) = 0` on every eigenspace | **C-11 / O-4** — an *admissible* flipper exists ⟺ trace-balanced |
| `build_writer` | explicit admissible `U`, `[U,H]=0`, `U†RU = −R` | **C-11** |
| `commuting_family` | the independent record **bits** | **C-14** |

## THE COUNT LAW — **C-14**

> ### `k = min over eigenspaces of v₂(m_E)`

Each independent record must **halve every eigenspace**, so the family size is how many times every
multiplicity can be halved. **22 spectra, 22 PASS.** The naive `floor(log₂ min m_E)` fails on
`[3,3]`, `[6,6]`, `[5,5]` — measured `0, 1, 0`.

**On the toric code the ground multiplicity is 4, so `k = v₂(4) = 2 = 2g`.**
**The topological formula is a consequence on that carrier, not the source of the count.**

## THE BOUNDARY THE MODEL ENFORCES — **C-15**

**Clauses (i)–(iv) are carrier-free.** **Clause (v) is not.** Protection needs a locality structure,
which is carrier data. `m.protection()` **raises** rather than supplying a default, so any claim
resting on (v) visibly inherits a carrier.

## WHAT THE MODEL CANNOT DO

- **Dimension — FIXED (T-3).** The model reaches **dim 256 in seconds**. It never builds the `*`-algebra
  (O-19) and never builds a basis of the commutant: `minimal_projections` needs **one generic Hermitian
  element of `A′`**, obtained by projection — twisted averaging over unitary generators, exact block
  projection for Hermitian ones. A commutant *basis* is still `O(n²)` and stays a lazy property, unused
  on the critical path.
- **Clause (v).** By design, above.
- **Formation.** The model constructs records; forming one needs an environment (see `F-13`).

## FORMATION — the model computes the process, not only the object

```python
env = Environment(nq=3)                      # a qubit bath at inverse temperature beta
m.formation(record, coupling, env, lam, t)   # chi(record : bath) after unitary evolution
m.formation(..., fragment=[0])               # what ONE fragment holds (redundancy)
m.channel(record, coupling)                  # does this coupling open a channel at all?
m.channel_map(family, couplings)             # the dependency structure
m.formation_independence(family, couplings, env)   # can one form without disturbing another?
m.evolve(coupling, env, lam, t)              # the joint state, so many readouts share ONE eigh
m.redundancy(record, coupling, env, lam, t)  # what the whole bath and each fragment hold
```

**`evolve()` exists because `formation()` redid the eigendecomposition per fragment**, which put a
7-qubit bath out of reach. Anything needing several readouts of the same evolution should use it.

A coupling may be a **product** `A ⊗ probe`, a **distributed** list of `(A_i, j)` pairs, or a full
interaction operator. **The distinction is not cosmetic** — the lanes distribute each system term to
a specific bath qubit, and a product ansatz silently returns a different number.

**`channel()` is the criterion (G-16, corrected):** a coupling opens a channel iff **its compression
onto the code space has a non-zero component along the record**. Anticommuting with the writer is
implied by this and does not imply it.

## CARRIER INDEPENDENCE

Every carrier-dependent result is verified on **three** carriers — `[[8,2,2]]` toric, `[[8,1,2]]`
**non-manifold**, and `[[4,2,2]]` which is **not a lattice at all**: `LANE_T9_CARRIERINDEP`,
**32 PASS, 0 FAIL**. Note **D-12**: a complex can be redrawn as non-manifold without becoming a
different quantum system — only a change to the **stabiliser group** changes the carrier.

## VALIDATION

```bash
cd model && python3 validate_model.py       # existence half — expect 12 PASS, 0 FAIL
cd model && python3 validate_formation.py   # formation half — expect 17 PASS, 0 FAIL
cd model && python3 count_law.py            # the count law — expect 22 PASS, 0 FAIL
cd LANE_T9_CARRIERINDEP && python3 t9_sweep.py   # carrier independence — expect 32 PASS, 0 FAIL
```
