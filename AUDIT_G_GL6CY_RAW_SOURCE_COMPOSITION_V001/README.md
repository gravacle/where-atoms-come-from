# Independent audit — GL6CY raw source composition

**Target:** `GL6CY_DEVELOPMENT_CU_V002_RAW_SOURCE_COMPOSITION`  
**Disposition:** `PASS` only when `independent_reconstruction.py` and
`verify_audit.py` both pass.  
**Scope:** an exact, declared-path audit of the raw Hamiltonian source
composition

\[
 {\cal R}^{\rm raw}_{sr}=H_{,AB}\eta_{A,s}\eta_{B,r}
                         +H_{,A}\eta_{A,sr}.
\]

The audit uses the independently audited GL6CU V002 source-Jet engine as a
pinned input, but does not import the target composition program.  It
independently rebuilds all diagonal and 64 alternating-cycle source owners,
constructs the declared six tensor and three literal-centre source directions,
and compares all 64 Q4 character values with a fresh target replay.

The result is deliberately bounded.  It does **not** form a stationary
connected response, a CTP or 1PI kernel, a physical owner-once quotient, a
Ward null, GL6CR completion, Einstein gravity, \(C_R\), or \(G\).

Run from the repository root:

```text
python3 -B AUDIT_G_GL6CY_RAW_SOURCE_COMPOSITION_V001/independent_reconstruction.py
python3 -B AUDIT_G_GL6CY_RAW_SOURCE_COMPOSITION_V001/verify_audit.py
```

The reconstruction intentionally executes a fresh full target replay; it is
the cost of comparing every one of the 64 character values rather than
certifying only a \(k=0\) excerpt.
