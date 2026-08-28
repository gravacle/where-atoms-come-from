# Preserved prior GF V001 rejection

This file preserves the earlier independent hostile verdict rather than
overwriting it with the V002 re-audit.

The rejected V001 author target had the following key hashes:

- theorem: `eda0d3785fc04c184f18e8c8ef21918c7684e09f7773b86d0537e84362a4f454`
- result: `69847f4ee869fb0cc2411bbab90f6a62e2e7c37bc8285306d6a8d2c807e62b9b`
- observable contract: `a0ab7259f8ca11d7a6a9b861501b601cf94004888f4de660dd32b5b556a5faca`
- author verifier: `a4b38592cb2983efd7d50bf38348d140a78a1014cf7269e288134be8937a2799`
- manifest: `674970ab0ac3a69c2185c448a98f59b419555435b7cc475692473f02e0d0c737`
- seal file: `02fb266d18484d7e992ea829c188e2bce44d96f95656467bd8b0630a4c3a4a0c`

The V001 hostile audit rejected that target for three material omissions:

1. spatial TT rank two was allowed to stand in for a Wigner/Poincare
   helicity `+2,-2` representation, so a covariant two-scalar doublet was not
   excluded;
2. the all-character registry applied `Pi_TT(k)` at `k=0`, where its defining
   formula is singular; and
3. load-bearing correlators selected `|0>` without typing a degenerate ground
   space basis-invariantly.

V002 is a new target with new bytes. Passing or rejecting V002 does not erase
the V001 rejection.
