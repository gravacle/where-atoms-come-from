# THE UNIVERSAL RECORD MODEL — the program's executable representation

`model/project_model.py`

```python
from project_model import URM

m = URM()
s = m.surface(name, mechanism, dE, E_b, T, f0,
              provenance="<pinned source>")
```

The URM is where this program adds observations and laws. `model/record_model.py` remains its
first-principles corner engine: from `(H, {L_k})` it constructs every record admitted by the five
clauses. `model/project_model.py` carries that engine into the world tier with declared surfaces,
provenance refusals, layer methods, and validator gates.

Exactly one ledger row is `PROVED`: `C-71`. The URM does not upgrade a row's status; it makes the
row's computation, entry conditions, scope, and tests inspectable.

## THE THREE ENTRY DOORS

Every landed feature must use one of these paths:

1. **A new surface or venue** enters through a URM provenance gate. World-tier inputs require a
   pinned source; exact corners must self-declare `provenance="DEF-A"`. The consuming layer rechecks
   the declaration where bypass would otherwise be possible.
2. **A new law** enters as a layer method and a validator gate with a failing branch and a positive
   control. A sealed lane by itself is evidence, not a URM feature.
3. **An external number** enters as a pinned comparison beside computed output, with its units,
   tolerance, extraction semantics, and a power control stated. A source substring is provenance,
   not empirical confirmation.

Python metadata is enforceable at these program boundaries; it is not a cryptographic custody
mechanism. Sealed historical lanes may still contain raw constructors, but no new observation can
silently use them as the public path.

## THE LAYERS

| layer | public URM surface | what it carries |
|---|---|---|
| definition/laws | `clauses`, `lifetime`, `steady_value` | the five clauses, rates, and values |
| formation | `configuration`, `formation_occupancy`, `formation_orientation` | how a record is written and read |
| corner | `corner` | the DEF-A exact idealisation, with explicit self-declaration on `URM` |
| geometry/roles | geometry and role delegates | located record geometry and the three-role ledger |
| arrow | `arrow_threshold`, `arrow_ledger`, `arrow_invariance`, `arrow_history`, `arrow_redundancy`, `arrow_observation` | record-copy threshold, history, and fragments |
| count law | `census`, `count_widths` | the surviving-record staircase and both durability widths |
| classes | `coupling_venue`, `reachable_class`, `critical_kernel`, and related delegates | subcritical, critical, and supercritical coupling classes |
| writing | `writing_kernel_verdict`, `writing_uniformity`, `writing_transport`, `writing_trail_*`, `writing_gap` | conservation, criticality, transport, and trail diagnostics |

The field-instrument family is not listed: T-51 is still independently unverified and nothing from
that lane is registered or folded into the URM.

## THE CORNER ENGINE

```python
from record_model import RecordModel

r = RecordModel(H, Ls)        # Hamiltonian and Lindblad operators are the entire corner input
records = r.records()         # every record the pair admits
family = r.independence(records)
```

No lattice, gauge group, temperature, coupling constant, code, or geometry is invented on this
path. If a value is not derivable from `(H, {L_k})`, the corner engine does not have it.

| step | what it computes | theorem |
|---|---|---|
| `star_algebra` | `A = alg{I,H,L_k,L_k†}` | C-9 |
| `commutant` | `A'` | C-9 |
| `minimal_projections` | a maximal splitting allowed by `A'` | C-10 |
| `clause_iii` | non-triviality on an eigenspace | anchor (iii) |
| `clause_iv` | trace balance on every eigenspace | C-11 / O-4 |
| `build_writer` | an admissible `U`, with `[U,H]=0` and `U†RU=-R` | C-11 |
| `commuting_family` | independent record bits | C-14 |

The exact corner count is `k = min_E v2(m_E)`: every independent record must halve every
eigenspace. The naive `floor(log2 min m_E)` control fails on `[3,3]`, `[6,6]`, and `[5,5]`.
Clause (v) remains carrier data; the model raises rather than inventing locality.

## FORMATION AND OBSERVATION

`RecordModel.evolve` shares one eigendecomposition across many readouts. `formation`, `redundancy`,
and the arrow layer then score what the environment and its fragments hold. A coupling may be a
product operator, a distributed list of system-term/bath-site pairs, or a full interaction operator;
those are physically different inputs and are not silently interchanged.

`channel()` uses the corrected G-16 criterion: the coupling's compression onto the code space must
have a non-zero component along the record. Anticommuting with the writer is necessary on the gated
venues and is not sufficient in general.

## VALIDATION

From the repository root:

```bash
python3 model/validate_urm.py          # all four landed families, then geometry/project chain
python3 model/validate_project.py      # URM surface/delegate/D-25 gates
python3 model/validate_geometry.py     # geometry gates, then project chain
python3 model/validate_formation.py    # formation engine — 17 checks
python3 model/validate_model.py        # corner existence engine — 12 checks
python3 model/count_law.py             # exact corner count — 22 checks
```

The family validator is the landing gate for new URM features. Its printed totals are computed by
the runners; documentation describes the expected composition but the exit status is authoritative.
