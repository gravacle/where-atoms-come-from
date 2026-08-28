# GRA-FP q4/F3 RGRL-B constraint-origin screen

This lane tests whether the fixed degree-two diamond-ice parent already owns
the rank-two constraint algebra required by RGRL-B.  It adds no fields or
interactions.

Run:

```bash
python3 LANE_GRA_FP_F3_Q4_RGRLB_CONSTRAINT_ORIGIN_SCREEN_V001/verify_rgrlb_constraint_origin.py
```

The replay uses only Python's standard library.  It verifies dependency
custody, the six-state ice/pair ranks and identities, `S4` characters, exact
periodic incidence ranks, closed-ring Gauss preservation, continuum constraint
symbols, TT dimension, and the helicity-one/helicity-two character mismatch.

Read `THEOREM.md` for the precise proof boundary and `SELF_AUDIT.md` for the
hostile ceiling review.  The lane does not edit `MODEL.md`, registers, or any
canonical integration surface.

