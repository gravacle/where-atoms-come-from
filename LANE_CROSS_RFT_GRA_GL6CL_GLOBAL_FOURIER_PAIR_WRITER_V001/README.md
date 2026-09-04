# GL6CL Global Fourier Pair Writer

This packet derives the exact parent/child Fourier symbol of the audited
`GL6CH` complete `T2` pair-source-to-ring writer on the infinite `Q4` parent.  It proves
that the six same-parent pair directions remain accessible for sufficiently
slowly varying common fields when the locked read is included, and supplies
an explicit analytic left inverse.  The locked read supplies `A1+E`; the
complete writer supplies `T2`.  The larger unprojected canonical-direct row
is retained only as bookkeeping because arbitrary-profile off-diagonal
`A1/E` completion is not classified.

Run from this directory:

```text
python3 derive_global_fourier_pair_writer.py
python3 verify_packet.py
```

The calculation also exposes genuine limits: independent parent/child pair
fields are underdetermined, the tensor normal has an exact cubic `T2` block
at quadratic momentum order, and exact rank loss occurs at finite momentum.
Because `T2` is only part of the full `E2+T2` shear representation, that block
alone does not decide physical rotational anisotropy.  This advances
smooth-field gluing but does not establish
an autonomous source, stationary response, spacetime, metric, Ricci law,
gravity, or `G`.
