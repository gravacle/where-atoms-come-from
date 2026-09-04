# Independent hostile re-audit — repaired GL6CG global order-`h4` source operator

**Target:** `LANE_CROSS_RFT_GRA_GL6CG_DENSE_PARENT_GLOBAL_H4_SOURCE_OPERATOR_V001`  
**Disposition:** **PASS**  
**Prior defect:** repaired and retained as a regression test.

The repair is scientifically sufficient.  The target now distinguishes two
claims that its first sealed version conflated:

1. the **diagonal** coefficient identity is valid on every locked state of
   every simple degree-four bipartite incidence graph with consistent local
   port labels; and
2. that diagonal result is the **complete** order-`h4` source operator only
   on the inherited degree-four, girth-at-least-six parent.

No material defect remains in the repaired claim.

## 1. Independent reconstruction

The re-audit does not import or execute the author derivation.  It combines:

- a literal Rayleigh--Schrodinger recurrence over all four irreducible
  `i,i,j,j` diagonal return paths;
- all `6*3^4=486` radius-one locked neighborhoods;
- a new full-graph evaluation of (CG.4) on all 90 labelled locked states of
  `K4,4`, all eight nodes per state, all 16 links, and all 120 unordered link
  pairs;
- an exhaustive `Q4` incidence check and an analytic infinite-lift step
  difference check; and
- the retained 24-path alternating-square counterexample plus the independent
  64-transition `Q4` commutator replay.

The exact reconstruction is frozen in `verify_gl6cg_independent.py` and
`INDEPENDENT_RESULT.json`.

## 2. Broad diagonal identity — PASS

For distinct links `i,j`, direct differentiation of the four permitted
diagonal return paths gives

\[
 {q_i+q_j\over2p_{ij}}+{q_{ij}\over p_{ij}^2}.
\]

Combining this with the normalized fourth-order subtraction gives the full
all-pair expression

\[
 V^{(4)}(s)=-{3N\over16}\sum_iq_i+
 \sum_{i<j}\left[{q_i+q_j\over2p_{ij}}+
 {q_{ij}\over p_{ij}^2}\right].
\]

Every link in a simple degree-four incidence graph has exactly six adjacent
partners.  For every disjoint pair, `p_ij=4` and `q_ij=q_i+q_j`; eliminating
those terms leaves only the central node and the other three links at each
neighboring endpoint.  The 486 possibilities therefore exhaust the local
data relevant to one node even when longer cycles correlate those choices.
Every case gives

\[
 V^{(4)}_{v,ab}=-{4\over9}-{37\over12}z_{va}z_{vb}.
\]

The raw census has three coefficient vectors, each 162 times.  Their `A1`
contraction is always `7/2`; all three `T2` contractions vanish pointwise.

The new full-graph `K4,4` replay is the hostile check against the locality
argument.  It retains every four-cycle correlation rather than treating the
four neighboring endpoint words as independent.  All `90*8=720` node
coefficients reproduce the same identity exactly, including pointwise zero
diagonal `T2`.  Thus the earlier square counterexample does not damage the
broad diagonal theorem.

## 3. Complete-operator scope — PASS

A nontrivial difference between two locked degree-two states is a disjoint
union of alternating even cycles.  An off-diagonal locked-to-locked process
using at most four microscopic flips would therefore require an alternating
four-cycle.  On a simple bipartite graph of girth at least six, none exists.
Consequently the diagonal coefficient is the complete order-`h4` operator on
that restricted domain.

The audit independently checked the inherited parents:

- `Q4` has 128 nodes, 256 links, and no pair of distinct parent vertices with
  more than one common child.  Hence it has zero four-cycles.  A literal
  six-cycle exists, so its girth is exactly six.
- In the infinite lift, the twelve nonzero ordered differences
  `D[a]-D[b]`, `a!=b`, of the four simplex steps are all distinct.  A
  nontrivial square would require one such difference twice, which is
  impossible.

The repaired theorem, result, README, ledger, operator payload, and verifier
consistently restrict complete/no-`T2` language to this inherited
degree-four, girth-at-least-six domain.  The legacy packet identifier does not
override that explicit scope.

## 4. `K4,4` counterexample regression — PASS

The earlier kill test remains valid and is now used correctly.  In `K4,4`,
take the locked occupation formed by the two perfect matchings

```text
{L0-R1,L1-R0,L2-R3,L3-R2}
{L0-R2,L1-R3,L2-R0,L3-R1}.
```

The square `L0-R0-L1-R1-L0` is alternating.  Summing all `4!=24` orders and
Hermitian-averaging the two endpoint source conventions gives

\[
 V^{(4,H)}_{ts}=(15/4,0,0,0,0,0)
\]

in pair coordinates and

\[
 (A,E_a,E_b,T_1,T_2,T_3)=
 (5/8,5/16,15/16,15/8,0,0).
\]

Thus a girth-four graph has a nonzero off-diagonal tensor-irrep source.  The
repaired target explicitly says that (CG.1)--(CG.2) remain valid diagonally
while such square operators lie outside its complete no-`T2` theorem.  This
is exactly the necessary resolution: `K4,4` is retained inside the broad
diagonal domain and excluded only from the complete-operator domain.

## 5. Selected `Q4` row and boundaries — PASS

The independent replay still matches all 64 frozen outgoing transition
profiles exactly.  It reproduces:

```text
k0 row rank / E2 rank          2 / 2
k0 kernel                     A1 plus T2
fully k0-dark transitions     16
dark transitions open at k1  16/16, in E2
full first-moment row rank    6
```

Every transition profile is nodewise zero in `A1` and `T2`; this is not a
moment cancellation.  The target continues to distinguish the diagonal pair
source from the native link `T2`, calls the selected state nonstationary, and
does not promote the commutator row to a pole or response.  It claims no
record authentication, physical coordinates, metric, Ricci law, gravity, or
value of `G`.

## 6. Custody — PASS

The target author verifier passes `484/484`, including all 24 frozen
dependency rows and a fresh `69470/69470` derivation replay.  The target
manifest authenticates every payload once, and its one-line seal authenticates
that manifest.  The audit pins the repaired target bytes, the inherited
GL6AO/GL6CC scope evidence, its own independent result, and its own replay.

## 7. Verdict

**PASS.**  The exact diagonal identity is broad; the exact complete no-`T2`
operator conclusion is narrow.  The repaired GL6CG now states and enforces
that distinction, preserves the `K4,4` falsifier as the boundary case, and is
fit for downstream use within the inherited girth-six gravity lane.
