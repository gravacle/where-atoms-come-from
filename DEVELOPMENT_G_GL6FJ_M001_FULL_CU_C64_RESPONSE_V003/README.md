# GL6FJ — fixed-m001 complete CU/C64 response flight

This launch-binding packet carries the shortest response flight selected by
GL6FI: the complete 300-ray real symmetric
response of the audited GL6FA
stationary engine at the non-self-conjugate `L=4`, `m=(0,0,1)` character.
The existing GL6FE flight supplies `m123`; this packet does not duplicate it.

The 300 canonical directions are dimension-minimal for `Sym(24)`.  Every
direction uses the same source-zero random history and retains all six
owner-once terms.  Reconstruction preserves the normal and anomalous P/C
blocks and the real AA/AO/OA/OO blocks.  Translation symmetry predicts a
zero anomalous expectation at non-self-conjugate `m001`, but the finite-run
estimate is retained exactly and is never projected to zero.

## Authorization boundary

Launch is not authorized by this author packet.  GL6FG V003, GL6FH V004, and
GL6FI V003 are now independently audited and load bearing, and their exact
target/audit chains are bound by the 76-row dependency ledger.  This V003
target must itself pass a distinct hostile audit before its controller can
authorize a production flight.

V003 changes admission, version, dependency, and authorization surfaces only.
The C++ measurement engine, its 300-ray physics flight, and the V002
reconstructor are preserved byte for byte.

No range inverse, optical Schur complement, m123 composition, physical Ward
rank, Einstein/Fierz–Pauli identification, gravity, `C_R`, or `G` follows
from this target alone.
