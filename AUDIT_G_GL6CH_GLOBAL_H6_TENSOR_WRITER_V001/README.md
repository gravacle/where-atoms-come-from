# GL6CH independent hostile audit

This packet independently audits
`LANE_CROSS_RFT_GRA_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001` without importing,
executing, or modifying its author derivation.

The replay separately reconstructs:

- both sets of `720` direct alternating-hexagon histories;
- the canonical full source gradient and its pure-`T2` projection;
- the lower `h0`, `h2`, and `h4` tensor-source exclusion by an exact
  differentiated finite-star Rayleigh calculation on all `486` locked
  radius-one neighborhoods;
- every simple six-cycle of `Q4`, the no-four-cycle condition, owner counts,
  and the infinite-parent step-difference test;
- the four-orientation tensor Gram and rank; and
- folding, sign, source-dimension, remainder, generic-graph, and
  interpretation boundaries.

Run from the repository root:

```text
python3 -B AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/verify_gl6ch_independent.py
python3 -B AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/verify_packet.py
```

