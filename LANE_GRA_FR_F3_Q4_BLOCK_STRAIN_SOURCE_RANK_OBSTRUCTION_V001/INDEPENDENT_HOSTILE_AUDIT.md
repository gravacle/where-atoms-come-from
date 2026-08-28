# Independent hostile audit

**Lane:** `GRA-FR-F3-Q4-BSSRO-V001`

**Audit date:** 2026-08-27

**Verdict:**
`PASS_AFTER_CONTACT_AND_COMPLETE_SOURCE_SCOPE_REPAIRS__ADDITIVE_EDGE_LINEAR_SOURCE_ONLY__NONEDGE_ROOT_AND_ROTATED_COFRAME_QUERIES_OPEN`

## 1. Material defects found and repaired

### 1.1 General contact derivatives

The submitted theorem originally promoted invariance of the four additive
linear contractions to invariance of the **complete** source-deformed
Hamiltonian and then stated that every CTP derivative with an `E` leg
vanishes. That did not follow from the pinned BS20a parent. BS20a explicitly
retains a general Hermitian `H_contact=O(j^2)`. Such a contact has zero first
derivative at `j=0`, but it can have a nonzero Hessian with one or two `E`
legs and can therefore break finite-source `E`-shift invariance.

The theorem, result, README, and self-audit were corrected before promotion.
They now distinguish two exact statements:

1. every contact that factors only through the four invariants `x_a=j:D_a`
   preserves the `E` null at every source order; and
2. a general prospectively frozen `O(j^2)` contact can change higher source
   derivatives, but it cannot change the linear source operator at `j=0`.

This repair does not weaken the scored obstruction. FQ's required retarded
packet and commutator moments are built from the source-off linear operator
`Q_eff=-2 partial_j H_eff|_0`. Both `E` contractions of that operator remain
identically zero. A seagull Hessian is not a missing linear conjugate and does
not supply a `Q_E` spectral residue or commutator pole. A contact genuinely
derived and frozen from the same parent is not automatically a new query;
only a post-result insertion would violate the preregistration.

### 1.2 The additive edge sector is not the complete BS20 source

The submitted packet also promoted the rank of FQ17a's one-/multi-edge
additive weights to the rank of the complete Q4-BLOCK source. The pinned BS20
and FQ texts do not permit that step. They require onsite, node, port,
boundary, and controller linear weights to be frozen separately, while FQ17a
specifies only the affine dyad rule for edge-supported terms. Their tensor
span has not yet been derived, and the packet may not silently set them to
zero or assume that they lie in the four edge-dyad span.

The corrected theorem now closes only the additive edge-supported linear
source subclass. It explicitly leaves the complete source rank open and
nominates `Q4-COMPLETE-SOURCE-RANK-AUDIT`: freeze the complete same-parent
term list, derive every remaining linear tensor weight, and rank their union
before performing the expensive CTP spectrum. This preserves the exact early
obstruction: any successful complete source must obtain both missing `E`
directions outside the additive edge sector.

## 2. Independent rank and representation replay

Using the four unnormalized tetrahedral sign vectors, the audit reconstructed
the dyads without importing the lane verifier. Their Gram matrix has diagonal
`9`, off-diagonal `1`, uniform eigenvalue `12`, three contrast eigenvalues
`8`, and rank four. The separate four-by-six source contraction map also has
rank four. Its complete kernel is

`span{diag(1,-1,0), diag(1,1,-2)}`,

the two-dimensional diagonal-traceless sector. The character of the
four-label permutation module after removing its uniform line is
`(3,1,-1,0,-1)` on the frozen class ordering, independently confirming the
`T2` label. Thus the exact local decomposition is `A1+T2`, with `E` absent.

## 3. Additive closure and blocking

The audit tested several independent signed coefficient families, a stronger
algebraic class than nonnegative occurrence counts. Every multi-edge weight
has the form `N C`; hence its rank is at most four and both `E` null vectors
survive. A three-block direct sum has rank `12=3x4`, and arbitrary
source-independent output mixing cannot increase that rank. An invertible
change of the six source coordinates also leaves the rank four. Therefore a
different **basis** or ordinary blocking map does not evade the frozen
FQ17a conclusion.

That statement is intentionally not a gluing no-go. A rational exact
orthogonal rotation generated from quaternion `(1,2,3,4)` was used to rotate
one tetrahedral coframe. Each coframe separately remains rank four, while the
union of the original and independently rotated spans has rank six. Thus a
physically derived, independently rotating block coframe could reopen the
question. It is not a basis change inside the frozen uniform q4 affine
coframe and is correctly left outside this theorem.

## 4. Root-edge query boundary

The six sibling-root dyads were reconstructed independently and have rank
six. The corresponding six additive two-edge weights have rank four. Their
difference is exactly the cross-dyad part of
`(n_b-n_a)(n_b-n_a)^T`. Therefore a prospectively derived blocked-root or
cross-dyad source is a mathematically real successor. Its existence limits
the theorem's scope; it does not alter the rank of the already frozen
additive source and cannot be substituted after Feshbach reduction.

## 5. Feshbach, CTP, and commutator replay

For a separately constructed exact rational block Hamiltonian, the audit
differentiated the fixed Schur/Feshbach map analytically. Both microscopic
`E` directional derivatives vanish, both effective derivatives vanish, and
the full six-coordinate family of effective first derivatives has rank at
most four. A quadratic `E` seagull was then added as a hostile countercase:
it changed the finite-source Hamiltonian and effective Schur complement but
remained even in the `E` amplitude, so its first derivative at source off was
exactly zero. A mixed `E-A1` contact supplied an explicit nonzero Hessian with
an `E` leg, validating the necessity of the wording repair.

A generic pulled-back noncontact response `C^T K C` has rank four and both
exact `E` nulls. Separately constructed operator-valued linear sources vanish
when contracted with either `E` vector before any state is chosen. Their
nested commutator operators and every tested commutator moment vanish through
order five. These finite representative calculations illustrate the general
operator chain-rule proof; the decisive fact is `D_E H|_0=0`, not the chosen
matrix entries.

## 6. Promotion ceiling

The audited conclusion is exactly:

`frozen uniform-coframe FQ17a additive edge-supported linear source -> rank 4 A1+T2 -> E null 2 -> fixed Feshbach preserves the source-off linear null -> the additive-edge-only six-direction/conjugacy prerequisite fails`.

It is not a theorem against:

- all collective or nonlinear variables of the H6/H8 parent;
- the complete BS20 source before its non-edge linear weights are derived;
- a prospectively derived root-edge, cross-dyad, loop, or surface source;
- an independently rotating physical coframe or new gluing architecture;
- thermodynamic tensor poles or emergent constraint algebras;
- adopted RGRL-B, gravity, or `G`.

The independent executable completed `105/105` checks, including all eleven
pinned dependency digests. The audit therefore passes the corrected packet
at its narrow, scientifically useful scope.
