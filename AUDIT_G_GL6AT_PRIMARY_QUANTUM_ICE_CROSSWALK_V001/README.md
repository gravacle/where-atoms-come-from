# Independent hostile audit — GL6AT V001

**Target:** `LANE_CROSS_RFT_GRA_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001`  
**Audit date:** 2026-08-31  
**Disposition:** PASS at the exact snapshot in `AUDITED_TARGETS.sha256`

This audit was performed independently of the GL6AT author. It replays the
exact graph, Hilbert-space, ring-operator, parameter, representation, and
small-frequency power-counting checks without importing the author verifier.
It also checks the four primary papers at equation level and keeps numerical,
Gaussian-field, and slave-particle evidence below exact-theorem status.

The earned result is an order-six operator crosswalk only: after a common
scalar is removed, GL6AO is the two-dimers-per-diamond-site fully-packed-loop
Hamiltonian at `v/g=0`, with `g=(63/8)h^6/U_d^5>0`. The audit does not promote
that identification to all orders, a rigorous thermodynamic phase, an exact
pole or pair-`E` overlap, physical momentum or speed, a cone, gravity, or `G`.

Run:

```text
python3 AUDIT_G_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001/independent_gl6at_replay.py
python3 AUDIT_G_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001/verify_audit_packet.py
```
