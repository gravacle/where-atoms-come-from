# Result: minimal FU09b encoded-current lift

The FU09b charge-conserving flip has an exact minimal realization.  Pair each
link qubit with one reservoir qubit and encode the original `Z=-1,+1` states
as `|- ,+>_LR` and `|+ ,->_LR`.  On this fixed-total-charge subspace,

\[
 \widetilde X V=VX,
 \qquad (q_*Z+Q_R)V=0.
\]

The dressed dynamics is therefore unitarily equivalent to the inherited link
dynamics, and the construction tensors to the full q4 node and any finite link
set.  It produces a real internal U(1) exchange current satisfying
`qdot_L+I_(L->R)=0` and `Qdot_R-I_(L->R)=0`.

Reservoir placement is not determined by that algebra.  Co-location with the
FY link midpoint gives an exact nonzero-momentum link/port continuity witness;
moving the reservoir elsewhere requires an additional connector current.
Only global total charge is allocation independent.  An explicit exterior
port completes the three-factor charge ledger, but a reservoir-only active
port leaks out of the encoded hold subspace.  No spatial bond current or
vertex divergence has been derived.

On a closed response hold where every added reservoir/support/port term is one
common scalar across the full encoded `P+Q` Hilbert space and has only a common
identity spatial derivative, the complete nonidentity Coulomb pair and
flip-numerator sources intertwine exactly.  A scalar `c[j]I` is retained as a
reference term: it cancels from virtual gaps and, after shifting the Feshbach
energy reference by `c[j]`, descends as `c[j]I_P`.  Thus nonidentity H6
coefficients, folds, commutators, ranks, matrix elements, and connected
responses remain unchanged; the uncentered Hamiltonian/source are literally
equal only when `c[j]=0`.

Source independence by itself is not enough: a constant reservoir bias
descends to `Z` and changes the source-off gaps.  Shared charging descends to
cross-link `ZZ`, an active port changes the source-off parent, and
strain-dependent nonidentity reservoir or port terms require a fresh source
audit.

This advances U(1) charge-current ancestry only.  The scalar exchange current
is not `T^{0j}`, a metric Ward identity, or gravity.
