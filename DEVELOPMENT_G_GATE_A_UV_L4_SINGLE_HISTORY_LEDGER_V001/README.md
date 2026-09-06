# Gate-A UV L=4 single-history ledger witness

**Status:** development-only, exact model-conditional witness; independent
hostile audit required before any promotion.

This bounded packet instantiates one finite `L=4` F3 blank-target write in a
periodic 64-cell family and evaluates a discrete retained-lineage ledger.  It
does not alter the frozen `m001` response result: that result is a stationary
response sampler and contains no cell-history ledger variables.

The physical source map used here is declared rather than derived from bare
F3.  It is a single admissible member of the adopted F3-MDC working-law
domain, chosen so that the audited raw source-Jet is given one explicit
Hamiltonian write attachment.  The witness therefore proves an exact balance
for its stated history only.  It does not prove that F3 uniquely supplies this
attachment, that every physical history has the same ledger, a physical Ward
identity, accumulation, gravity, or `G`.

Run:

```text
python3 -B verify_single_history_ledger.py
```
