# FX complete-H6 source-response lane

This isolated no-lab gravity lane closes the diagonal-source omission left
explicit in FV13.  It enumerates every diagonal closed word and fold through
H6 on the FO 180-state translation component, then composes the result with
the independently audited FW response.

Run the exact certificate from the repository root:

```bash
python3 LANE_GRA_FX_F3_Q4_COULOMB_COMPLETE_H6_SOURCE_RESPONSE_V001/derive_complete_h6_source.py
```

The run uses exact `Fraction` arithmetic for the full history enumeration
and can take several minutes.  The central result is that the generated
diagonal source renormalizes the existing pair `E` channel by `f_E(x)` and
adds identities; it creates no new nonidentity response direction.

See `THEOREM.md` for the proof and scientific ceiling, `RESULT.md` for the
short outcome, and `SELF_AUDIT.md` for overclaim controls.
