# Preliminary raw GL6CU mixed-source diagnostic

**Status:** exact development replay of the independently supported GL6CU
V001 arithmetic core; not load bearing until repaired GL6CU V002 receives a
distinct hostile `PASS`. This is only the first term of (CY02), not a complete
physical response or Ward residual.

## Frozen diagnostic convention

Use pair-row order

\[
 (01),(02),(03),(12),(13),(23)
\]

uniformly on all 128 constraint nodes. Use the three unit staggered center
fields

\[
\begin{aligned}
 t_x&=(1,0,0,0,0,-1),\\
 t_y&=(0,1,0,0,-1,0),\\
 t_z&=(0,0,1,-1,0,0),
\end{aligned}
\]

with `+t_i` at every parent node and `-t_i` at every child node. For a
dimensionless center displacement `b=e_i`, the literal physical first source
(CY15) is `-2` times this convention. The matrices below therefore isolate
the operator contraction without absorbing that physical factor.

## Exact replay

The selected-row diagonal order-two mixed Hessian vanishes:

\[
 {\cal H}^{(2)}_{\rm diag}=0_{6\times3}.
\]

At fourth order,

\[
 {\cal H}^{(4)}_{\rm diag}
 ={8032\over27}
 \begin{pmatrix}
 0&0&0\\
 0&1&0\\
 0&0&-1\\
 0&0&1\\
 0&-1&0\\
 0&0&0
 \end{pmatrix}.
\]

At sixth order the diagonal owners give

\[
 {\cal H}^{(6)}_{\rm diag}
 ={5010712\over6075}
 \begin{pmatrix}
 0&0&0\\
 0&1&0\\
 0&0&-1\\
 0&0&1\\
 0&-1&0\\
 0&0&0
 \end{pmatrix},
\]

while all 64 active alternating-cycle writers give

\[
 {\cal H}^{(6)}_{\rm cyc}
 =-528
 \begin{pmatrix}
 0&0&0\\
 0&1&0\\
 0&0&-1\\
 0&0&1\\
 0&-1&0\\
 0&0&0
 \end{pmatrix}.
\]

Thus the raw source-Hessian sum is nonzero:

\[
 \boxed{
 {\cal H}^{(6)}_{\rm diag+cyc}
 ={1803112\over6075}
 \begin{pmatrix}
 0&0&0\\
 0&1&0\\
 0&0&-1\\
 0&0&1\\
 0&-1&0\\
 0&0&0
 \end{pmatrix}.}
\]

This proves only that first-source `T2` darkness and closed-cycle
first-source cancellation do not also force the mixed Hamiltonian Hessian to
vanish. It does **not** establish a physical obstruction. The second chain
term

\[
 H_{,A}\eta_{A,sr}
\]

and the support, GD recoil, state, spectral/contact, constraint, boundary,
moving matching/writer, retained-field, and quotient owners are absent from
this matrix. The completed quotient remainder must also undergo the exact
bulk/exact/coexact/boundary and normalized-refinement tests in TARGET.md.

