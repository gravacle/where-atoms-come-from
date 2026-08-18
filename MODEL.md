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

- **Dimension.** The commutant is found by an `n²`-dimensional nullspace, so exact use is small-`n`.
  Large carriers are handled in their own lanes and cross-checked against the model where they overlap.
- **Clause (v).** By design, above.
- **Formation.** The model constructs records; forming one needs an environment (see `F-13`).

## VALIDATION

```bash
cd model && python3 validate_model.py    # expect 12 PASS, 0 FAIL
cd model && python3 count_law.py         # expect 22 PASS, 0 FAIL
```
