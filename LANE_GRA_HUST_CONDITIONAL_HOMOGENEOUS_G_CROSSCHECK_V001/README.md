# HUST conditional homogeneous-kernel G cross-check

This lane combines the audited HUST response forward with the independently
reconstructed nominal source kernels. It computes exactly what follows under
an explicit zero full-kernel-remainder premise, retains the open ToS correction
as an affine parameter, propagates only available local uncertainty fields, and
shows why this independently reconstructed conditional family identifies no
physical \(G\) point or deterministic compact interval from the public packet.
The paper's own processed-model \(G\) summaries remain acknowledged as
post-calculation comparators.

Run:

```bash
python3 -B LANE_GRA_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001/calculate_conditional_homogeneous_g.py
python3 -B LANE_GRA_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001/verify_conditional_homogeneous_g.py
```

The processed HUST source coefficients and already-derived \(G\) rows are
parsed as part of pinned parent payloads, but their keys are selected only in
the post-calculation comparator phase. No accepted value of \(G\) is used.
