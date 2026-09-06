# GL6CY complete physical source/reembedding residual target

This is a **bounded calculation target**, not a theorem packet and not a
gravity result. It specifies the shortest honest composition after GL6CX
V002 and the repaired GL6CU source Jet:

\[
 {\cal R}_{rs}^{\rm src}
 =H_{,AB}\eta_{A,r}\eta_{B,s}+H_{,A}\eta_{A,rs}.
\]

The calculation must then add every non-source physical owner exactly once,
perform the complete connected-to-1PI/Schur/Legendre/quotient construction,
and only then score a Ward residual. An omitted owner is reported as
`UNDEFINED`; it is never set to zero.

The target separates:

1. the exact barycentric source embedding and its first and second Jets;
2. the raw GL6CU Hamiltonian source-Jet composition; and
3. the complete physical quotient residual and its
   bulk/exact/coexact/boundary/refinement classification.

Item 1 is closed by the included exact executable.  The raw source-Jet
prerequisite in item 2 is independently audited by GL6CU V002; only its
non-source physical completion remains open. Item 3 requires the still-unowned physical position/coframe, state, recoil,
apparatus, retained-field, and quotient data listed in `OWNER_SCHEMA.json`.

No nonzero one-cell residual is interpreted as a failed Ward route or as a
reason to alter A3 kinematics. It must first be completed, accumulated,
decomposed into bulk/exact/coexact/boundary pieces, and tested in the
normalized joint refinement.
