# Independent hostile audit — GL6AZ record-authenticated prethermal mission

This packet is a distinct post-repair hostile audit of the frozen author
packet
`LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001`.
It pins every author byte, replays the mathematics without importing the
author verifier, checks dependency custody, and fails closed on scope
promotion.

The audit found the original source-domain omission before sealing this
packet.  The author repair restores `R>=bar(nu)_0`, excludes the spurious
low-ratio logarithmic branch, distinguishes all three sufficient-domain
floors, and gives an explicit inside/outside-domain continuation.  The
repaired packet passes this audit without a remaining material defect.

Run:

```text
python3 AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/independent_gl6az_replay.py
python3 AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/verify_audit_packet.py
```

No author byte is modified by this audit.  No graviton, Ricci template,
gravity identification, or numerical `G` is introduced.
