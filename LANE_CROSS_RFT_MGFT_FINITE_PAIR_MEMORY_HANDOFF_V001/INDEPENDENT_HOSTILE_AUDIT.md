# Independent hostile audit -- finite pair-memory handoff and composition

**Audit target:** `CROSS-RFT-MGFT-FPMH-V001`  
**Disposition:** `ACCEPT_AFTER_NARROW_MINIMALITY_REPAIR__NO_REMAINING_MATERIAL_DEFECT`

## Base handoff witness

The handoff unitary was independently rederived as a direct sum of two swaps
and identities.  It is Hermitian, unitary, and exactly reversible.  More
strongly than the displayed bit-flip test alone, the two exchange blocks
intertwine the same arbitrary unitary on the two-dimensional content subspace;
therefore the writer transports an arbitrary mixed or coherent content state
without learning its value.  The formation and sham arms have exactly the same
post-write carrier density operator, while only the actual custody handoff
writes the fresh relation bit.  There is no cloned content at the endpoint.

`KEEP=SWAP(L,K)` and `BREAK=SWAP(L,G)` are related by the `K <-> G`
permutation, preserve one total relation excitation, and are individually
reversible.  The quarantine is retained, not traced out.  With the common
pair-gated pulse, the exact query table is

| Arm | KEEP `(K,a)` | BREAK `(K,a)` |
|---|---:|---:|
| formation | `(1,1)` | `(0,0)` |
| sham | `(0,0)` | `(0,0)` |

so the formation/sham total-variation contrast is one under KEEP and zero
under BREAK.  This is an exact state-sector support intervention in the finite
model.  It is not yet an experimental claim that controller work, recoil,
electromagnetic, thermal, and stress ports are all matched.

## Composition and quantum-information boundaries

Induction on the unchanged handoff map gives the serial `P_(M+1)` custody
lineage with `M` retained relations, `O(M)` physical factors and operations,
and an extensive reversible BREAK.  This is linear physical-resource scaling,
not polynomial Hilbert-space dimension.

For the two-slot branch/rejoin construction, the two branch maps commute, the
two rejoin maps commute, and linearity transports an arbitrary joint -- even
entangled -- two-content state to the two recipient slots.  The resulting
support is exactly `C4`, with four vertices, four edges, and cycle rank one.
It is only a finite relational cycle, not a spatial plaquette.

The no-cloning obstruction is valid: preservation of the overlap of two
nonorthogonal inputs is incompatible with two perfect copies.  The one-slot
merge obstruction is also valid: a reversible garbage-free map cannot embed
the four-dimensional two-carrier content subspace into one qutrit.  The packet
correctly uses two source and two recipient slots and claims neither copying
nor quantum-state fusion.

The hostile audit found and repaired one narrow overstatement.  `C4` is not
absolutely the smallest lawful two-token branch/rejoin cycle: one direct route
and one two-hop route give a lawful asymmetric triangle.  The theorem now
correctly calls its `C4` the minimal **symmetric two-stage, equal-depth** cell.
No other result depends on the original wording.

## REC, FHBC, and universal-scope audit

The record is a derivative event record: pre-write carrier location differs,
and the common handoff writes that distinction into a new physical relation
coordinate.  Prior absence, positive source-off retention, complete
label-blind query, noncreation, physical lineage, and explicit quarantine are
all supplied inside the exact constructed episode.

FHBC is not inferred from contrast.  It is an explicit conditional premise of
each fixed finite model: finite device, exact outer factorization, one joint
root, arm-common later dynamics for a fixed prospective route context,
complete outcomes, and a finite source-to-query schedule.  Therefore the
imports are typed exactly as

\[
 \operatorname{FHBC}(r)\Rightarrow DCL_{\rm phys}(r),
 \qquad
 \operatorname{REC}(r)\land\operatorname{FHBC}(r)
 \Rightarrow\operatorname{COV}_{\cup}(r).
\]

This establishes finite model witnesses only.  It neither derives universal
FHBC/DCL membership nor authenticates a laboratory realization.

## Verification and ceiling

The two verifiers and all dependency hashes were independently rerun before
the audit file was added:

```text
Finite handoff:       113/113 PASS
Serial/composition:  128/128 PASS
Dependencies:        all OK
Manifest:            all pre-audit entries OK
```

The manifest must be refreshed to include this audit and the minimality wording
repair.  Acceptance does not establish robust autonomous pair-memory formation,
an open thermodynamic network phase, generic graph selection, dimension above
one, a continuum cone, tensor response, gravity, or `G`.
