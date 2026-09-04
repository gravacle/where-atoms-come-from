# GL6CH independent hostile-audit report

**Verdict:** `PASS`.  The independent replay passes `41466/41466` exact
checks against the frozen result.  No material mathematical, custody,
dimensional, graph-scope, or interpretation defect was found on the claimed
surface.

**Target:**
`LANE_CROSS_RFT_GRA_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001`

The target manifest is
`a895ecdb1ab5340634808c0d6d379e96a0b161f8756d44ba2460ac5c404a34e5`;
the target seal-file hash is
`e61abdab6ea225bb1c52c3d2f4e2050dc3dae4f70c0be1ba7712e369cc7ff61a`.
The audit pins all twelve target bytes and rechecks all twelve upstream
dependency hashes.

## Direct histories, sign, and source normalization

For each of the two alternating ring orientations, the audit independently
enumerates all `6!=720` orders.  Every proper prefix stays outside the locked
sector.  The nine exact denominator profiles occur with multiplicities

```text
(2,2,2,2,2):96   (2,2,2,4,2):48   (2,2,4,2,2):48
(2,2,4,4,2):96   (2,4,2,2,2):48   (2,4,2,4,2):24
(2,4,4,2,2):96   (2,4,4,4,2):192  (2,4,6,4,2):72.
```

With `W=-sum_e X_e`, six hop signs and five negative reduced resolvents give
the source-free sign.  Exact summation returns

\[
 \sum_{\pi\in S_6}w_\pi=-{63\over8}.
\]

The audit then constructs the endpoint-midpoint pair-memory score without
using the author formula.  Both alternating directions, all six vertices,
all six local port pairs, both cycle-edge assignments, and both exterior
assignments give `288` contexts.  Every context returns

\[
 g_v={105\over8}e_{ab},\qquad
 P_Tg_v={105\over16}\Theta_{v,c},\qquad
 \Theta_{v,c}=e_{ab}-e_{\overline{ab}}.
\]

Thus `||Theta||^2=2`, and the literal source direction `j=Theta` has
derivative `+105/8`.  The sign and the distinction between the full canonical
gradient and its tensor projection are correct.

## Independent lower-order and folding screen

The lower-order audit uses a distinct exact construction.  Around one
four-port node it includes its four incident edges, the twelve other edges
of the neighboring constraints, and distinct exterior endpoints implied by
girth six.  It differentiates

\[
 E_4=A B-P
\]

directly.  Here `P` is re-enumerated from all `480` four-flip identity words
whose three proper prefixes remain in `Q`.  Across all six central locked
words and all `3^4` compatible neighbor completions—`486` cases—the result is

\[
 V_v^{(2)}=-M_v,
 \qquad
 V_v^{(4)}=-{4\over9}{\bf1}_6-{37\over12}M_v.
\]

Every locked `M_v` has zero `T2` projection, so the first tensor-source
vertices at `h0`, `h2`, and `h4` vanish.  Girth six excludes a distinct
locked configuration change below six flips.  Consequently an order-six
off-diagonal fold has neither a lower off-diagonal operator nor a lower
pure-`T2` source vertex to multiply.  The direct result is therefore the
complete first off-diagonal tensor-source term at this order.

## Global geometry, owners, and tensor access

The audit constructs the `128` nodes and `256` links of `Q4`, proves
connectedness, and checks every pair of parent nodes for a four-cycle.  A
separate graph DFS obtains exactly `256` simple six-cycles and agrees
edge-for-edge with the canonical translation/three-port construction.  Each
orientation has `64` cycles, and every physical link occurs in exactly six
cycles.  Because a transition's symmetric difference fixes its cycle, the
operator sum has one owner per undirected elementary cycle.

For the infinite incidence graph, all twelve ordered nonzero differences of
the four step vectors are distinct over `Z^3`; hence two parent vertices
cannot share two children and there is no four-cycle there either.  This is
the needed infinite-parent check, not an extrapolation from the finite
quotient alone.

The independently reconstructed orientation vectors obey

\[
 \sum_du_d=0,\qquad \|u_d\|^2=24,\qquad
 u_d^Tu_{d'}=-8\ (d\ne d'),\qquad
 \sum_du_du_d^T=32P_T,
\]

and have rank three.  This supports full `T2` source access.  It does not by
itself produce a propagating tensor phase.

## Dimensions, graph scope, and physical wording

The stable target's remainder is dimensionally explicit.  With
`[h]=[U_d]=[j]=E`, all of

\[
 {h^6\over U_d^5},\quad {h^6j\over U_d^6},\quad
 {h^6j^2\over U_d^7},\quad {h^8\over U_d^7},\quad
 {h^8j\over U_d^8}
\]

have energy dimension one.  The displayed leading omitted source and hop
orders are consistent.

The generic-graph challenge also passes.  The theorem does not infer its
lower-order exclusion for every simple four-regular bipartite graph; it
expressly restricts the operator claim to the declared girth-six `Q4` and
infinite diamond-incidence parent.  `Arbitrary locked state` means every
locked basis state of that parent, not every graph.

Finally, the phrase `candidate-field-dependent future writer` is supported
in the following narrow sense: an externally charted microscopic pair field
changes the amplitude of a later physical ring transition.  The packet
explicitly does not prove that retained records generate that field, select
a phase, create a pole or common cone, produce a metric/Ricci response,
establish gravity, calculate `G`, or introduce a graviton.  Calling the
diagonal pair-memory image a `record` remains inherited program language,
not a new record-authentication theorem in GL6CH.

`PASS__GL6CH_INDEPENDENT_720_HISTORY_AND_288_CONTEXT_REPLAY__MINUS_63_OVER_8__CANONICAL_105_OVER_8_EAB__T2_105_OVER_16_THETA__DIFFERENTIATED_FINITE_STAR_H0_H2_H4_EXCLUSION__Q4_AND_INFINITE_GIRTH6__256_OWNER_ONCE_HEXAGONS__ORIENTATION_GRAM_32PT_RANK3__DIMENSIONS_AND_REMAINDERS_TYPED__CANDIDATE_WRITER_ONLY__NO_PHASE_RECORD_AUTH_METRIC_RICCI_GRAVITY_G_OR_GRAVITON`
