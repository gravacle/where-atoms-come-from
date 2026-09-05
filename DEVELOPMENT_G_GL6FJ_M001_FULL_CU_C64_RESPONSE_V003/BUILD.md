# Build and use

Author verification (nonproduction):

```sh
python3 -B DEVELOPMENT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003/verify_gl6fj_v003.py
```

Preflight (never launches):

```sh
python3 -B DEVELOPMENT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003/run_gl6fj_v003_pilot.py preflight
```

The launch token is intentionally documented only as custody inside the
controller and measurement plan.  Possessing it is insufficient: the sealed
target, every dependency, and the exact distinct hostile audit must all pass
before the controller acquires an output lock.  The upstream chain is current
in V003, but this target's distinct hostile audit is deliberately absent at
author freeze, so preflight remains unauthorized.
