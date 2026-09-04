# GL6CO Cycle-Response Tensor Matching

This packet classifies the most general reciprocal, inversion-even,
translation- and `S4`-invariant four-orientation cycle-response symbol through
quadratic coordinate momentum.  It composes that symbol with the complete
`GL6CL` `T2` writer and proves that extension of the resulting `T2-T2` block
to an `SO(3)`-covariant symmetric-tensor operator requires exactly one linear
condition:

```text
c + d = kappa/2
```

The cycle symbol is the bare susceptibility
`2 Re <0|T_c R T_c'|0>`; the `GL6CL` writer is applied exactly once.
The packet also keeps the `GL6CL` unnormalized common coordinate distinct
from the orthonormal common parent/child coordinate used by the `GL6BV`
contact; the latter converts the cycle Hessian prefactor from `mu^2` to
`mu^2/2`.

The condition is not enforced by tetrahedral symmetry, but it is compatible
with positivity.  A separately typed reconstruction of the `GL6BV` order-
`h^2` contact shows the exact same-state equation that would apply if one
completed functional owned both blocks.  A full six-component solder test
also proves that later Ricci comparison cannot be closed from the `T2-T2`
block alone.

Run from this directory:

```text
python3 derive_cycle_response_tensor_matching.py
python3 verify_packet.py
```

This is a symmetry and matching theorem.  It does not calculate the response
coefficients, select a phase, complete `E2+T2`, invert the response to a 1PI
kernel, prove masslessness/background stationarity, or prove a metric, Ricci
dynamics, gravity, or `G`.
