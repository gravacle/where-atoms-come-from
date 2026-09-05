# GL6FJ V003 independent hostile prelaunch audit

This packet independently reviews the frozen
`DEVELOPMENT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003` launch target.  It is
strictly a prelaunch audit: it authenticates the frozen target, replays its
nonproduction verification in ordinary and optimized interpreters, rebuilds
the C++ sampler independently in three modes, and checks the unaliased
`m001` character, 300-ray catalog, and six-owner checkpoint balance on
sampled rays.

It does not launch the 300-ray flight, invert a response, form a 1PI kernel,
assert a physical Ward identity, or claim gravity.  Passing it only permits
the target controller to decide whether its separately declared production
flight may start.

Run:

```text
python3 -B verify_packet.py
python3 -O -B verify_packet.py
```
