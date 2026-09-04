# GL6CN — Complete diagonal `h6` `T2` first source

This sealed author packet closes the diagonal complement to the `GL6CH`
six-cycle tensor writer.  Its exact result is

```text
complete diagonal h6 pure-T2 first-source vertex = 0 pointwise.
```

Start with `THEOREM.md` for the proof and exact scope, `RESULT.md` for the
short physics result, and `EXACT_LEDGER.json` plus
`derive_complete_diagonal_h6_t2_source.py` for the rational replay.

Run:

```text
python3 derive_complete_diagonal_h6_t2_source.py
python3 verify_packet.py
```

The first-source conclusion is complete through `h6` only.  This packet does
not turn the tensor writer into a phase, record, bulk geometry, Ricci law,
gravity theorem, or calculation of `G`; it also does not evaluate
source-second contacts or higher-order first vertices.
