# Independent hostile audit — GL6CM

This packet independently audits the finite-component stationary spectral
response of the global order-six writer.  It verifies the Perron--Frobenius
scope, exact spectral Gram factorization and kernel, common-rescaling null,
`W^T K W` pullback, isolated-ring zero, shared-star strict response, physical
coefficient, units, and all declared claim ceilings.

Run:

```bash
python3 verify_packet.py
```

The audit disposition is `PASS` on the writer-only spectral surface.  It
does not promote the result to a complete contact-plus-spectral Hessian,
authenticated record dynamics, bulk locality, gravity, or `G`.
