# GL6CQ — stationary-response moment sum rules

This packet turns every undetermined coefficient in the GL6CO analytic
cycle-response symbol into a zeroth- or second-moment projection of one bare
connected cycle susceptibility.  It then combines those observables with
the GL6BV same-state contact probability and rewrites the two GL6CO matching
conditions as exact real-space tests.

The central formulas are

```text
kappa=Z_T/3
b=-M_perp/12
c=-M_parallel/6+M_perp/12
d=-M_cross/6
```

and

```text
(mu^2/2)[-2 Z_T+M_perp-2M_parallel-2M_cross]
    +4 g_ct(2p-1)=0

-(mu^2/2)[Z_T+M_perp]+2 g_ct(1-4p)=0.
```

The cycle kernel is bare, `mu^2` is applied once, the additional half converts
to the orthonormal common parent/child source used by the contact, and `p` is
evaluated in the same stationary state.  These are quadratic-gradient
matching tests, not claims that a phase satisfies them and not a Ricci,
gravity, or `G` result.

Run from this directory:

```text
python3 derive_stationary_response_moment_sum_rules.py
python3 verify_packet.py
```
