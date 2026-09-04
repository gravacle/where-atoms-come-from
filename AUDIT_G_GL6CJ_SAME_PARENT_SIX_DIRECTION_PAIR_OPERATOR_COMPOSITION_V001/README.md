# GL6CJ independent hostile audit

This packet independently audits
`LANE_CROSS_RFT_GRA_GL6CJ_SAME_PARENT_SIX_DIRECTION_PAIR_OPERATOR_COMPOSITION_V001`.
It does not import or execute the author derivation and does not modify the
target packet.

The replay reconstructs the six-pair projectors, all six locked-word
diagonal rows, every simple `Q4` hexagon, the writer incidence at all `128`
nodes, exact kernels and generalized inverses, and the combined rank-six
map.  It also checks that the two maps are selected effective-operator
derivatives of one pre-Feshbach source while remaining only an operator-jet
closure.

Run from the repository root:

```text
python3 -B AUDIT_G_GL6CJ_SAME_PARENT_SIX_DIRECTION_PAIR_OPERATOR_COMPOSITION_V001/verify_gl6cj_independent.py
python3 -B AUDIT_G_GL6CJ_SAME_PARENT_SIX_DIRECTION_PAIR_OPERATOR_COMPOSITION_V001/verify_packet.py
```

