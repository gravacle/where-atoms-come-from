# GL6CS strict-lock scale separation

Run:

```text
python3 derive_strict_lock_scale_separation.py
```

The exact replay enumerates all 24 eligible locked-node ring incidences,
verifies that every change is nonzero and pure `E2`, and checks the exact
coupling prefactors and powers used in the same-state scale-separation
theorem.  This is an author calculation pending independent hostile audit.
