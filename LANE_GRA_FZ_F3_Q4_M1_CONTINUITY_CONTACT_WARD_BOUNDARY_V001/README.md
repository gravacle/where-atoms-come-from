# GRA-FZ: F3/q4 `m=1` continuity/contact/Ward boundary

This packet asks whether FY's exact native-support nonzero-momentum spatial
source already closes a continuity/contact Ward identity.  It proves the
sharp finite boundary and the exact TT quotient without adding an interaction
or modifying the record theory.

Run the fast exact replay:

```sh
python3 verify_continuity_contact_ward_boundary.py
```

The optional full replay rebuilds FY's expensive H6 native ledgers and applies
the finite Liouvillian range diagnostic:

```sh
python3 verify_continuity_contact_ward_boundary.py --full-liouvillian
```

An algebraic `ad_H` range pass is not a physical-current or Ward pass.  The
claim ceiling remains the finite Z30 parent through H6 under `FV-PURE`.
