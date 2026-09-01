# Distinct hostile audit — GL6AY finite-coupling prethermal locked bridge

This directory is an independent audit of the frozen author packet
`LANE_CROSS_RFT_GRA_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_V001/`.
No author byte was changed.

**Disposition:**
`FAIL__REPAIR_REQUIRED__GLOBAL_LOCK_PROJECTOR_AND_DRESSED_SUBSPACE_SCOPE`.

The primary-source normal form, F3 strong-support mapping, finite second
moment estimate, local-observable horizon, and first winding coefficient all
survive independent replay.  Two related infinite-volume/global-projector
claims do not:

1. `P=chi(N_def=0)` is used as though `P D_hat P` were an infinite-volume
   quasi-local interaction, although no such global projector belongs to the
   quasi-local algebra.
2. The local/potential conjugation estimate is used to claim control of the
   global dressed subspace `Y^*P`; the cited theorem supplies no uniform
   global projector or subspace-norm estimate.

Both defects have a direct local repair, described in `AUDIT.md`.  Promotion
requires author repair and a fresh independent audit.
