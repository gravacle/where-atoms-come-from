# GL6AK distinct post-freeze audit

This is the authoritative post-freeze custody and physics replay for GL6AK.
The older `README.md`, `AUDIT.md`, `VERIFICATION.txt`, and mutable replay are
preserved byte-for-byte because the frozen author packet pins them as its
pre-freeze review.  They are historical inputs, not the terminal audit status.

Run `python3 -B verify_postfreeze_audit.py`; after sealing, run
`python3 -B verify_audit_packet.py`.
