# Composed finite pair-memory custody theorem

**Lane:** `CROSS-RFT-MGFT-FPMH-V001`  
**Extension:** serial chain and one branch/merge cell  
**Claim class:** exact finite composition of the unchanged one-pair writer,
active/quarantine route, gated link, and imported FHBC/record theorems

**Disposition:**
`FINITE_SERIAL_PAIR_MEMORY_PATH_FORMED__LAWFUL_TWO_CARRIER_BRANCH_MERGE_CYCLE_FORMED__SINGLE_UNKNOWN_CARRIER_BRANCH_AND_SINGLE_SLOT_MERGE_OBSTRUCTED__LINEAR_FACTOR_RESOURCE_SCALING__EXTENSIVE_REVERSIBLE_BREAK__REC_FHBC_DCL_COVERAGE_MODEL_WITNESSES__ROBUST_NETWORK_HIGHER_DIMENSION_AND_GRAVITY_OPEN`

## 1. Composition rule and non-negotiable content constraint

Import the one-pair content-blind handoff unitary

\[
U_{i\to j}:|c,B,0\rangle_{V_iV_jL_{ij}}
\longleftrightarrow |B,c,1\rangle_{V_iV_jL_{ij}},
\qquad c\in\{0,1\},                                  \tag{C01}
\]

with the complement fixed. Every use has a fresh relation target. It transfers
the complete two-dimensional content subspace and does not clone it.

This restriction matters at a branch. There is no isometry that takes one
arbitrary unknown content state and two blanks to two perfect copies. For two
nonorthogonal inputs with overlap `s`, such a map would require

\[
s=s^2\langle g_\phi|g_\psi\rangle,
\]

whose right side has magnitude at most `|s|^2<|s|`. A branch made by cloning one
unknown carrier is therefore rejected.

Likewise, two arbitrary two-state carriers have four orthogonal basis states.
A reversible map that blanks both inputs and stores their complete joint content
in one qutrit with no garbage would embed a four-dimensional subspace into a
three-dimensional one. That single-slot merge is impossible. A lawful arbitrary-
content branch/merge cell must use two actual source carrier slots and two
recipient slots, or declare additional garbage/content capacity. The construction
below uses two slots and no hidden garbage.

These are physics conclusions, not reasons to insert a copying or merge rule.

## 2. Finite serial custody chain

Fix any finite `M>=1`. Use carrier qutrits

\[
V_0,V_1,\ldots,V_M
\]

and, for each edge `e_i={V_(i-1),V_i}`, fresh writer, active, quarantine, and
link registers `(L_i,K_i,G_i,a_i)`. All relation and link registers begin zero.

The formation source places an arbitrary content state `rho_C` at `V_0`; every
other carrier is blank. The matched prepositioned sham puts the same `rho_C` at
`V_M` and leaves every earlier carrier blank. Apply the one common ordered
writer

\[
U_{\rm serial}=U_{M-1\to M}\cdots U_{1\to2}U_{0\to1}.
\tag{C02}
\]

Induction on (C01) gives, for arbitrary mixed or coherent `rho_C`,

\[
\begin{aligned}
F:&\quad
\rho_C^{V_0}\otimes|B\cdots B\rangle
\longmapsto
|B\cdots B\rangle\otimes\rho_C^{V_M}
\otimes|1\rangle\!\langle1|_{L}^{\otimes M},\\
S:&\quad
|B\cdots B\rangle\otimes\rho_C^{V_M}
\longmapsto
|B\cdots B\rangle\otimes\rho_C^{V_M}
\otimes|0\rangle\!\langle0|_{L}^{\otimes M}.
                                                               \tag{C03}
\end{aligned}
\]

Thus formation and sham have exactly the same final carrier state, while the
actual sequence of `M` custody handoffs forms exactly the serial pair relation

\[
G_\ell=P_{M+1}.                                      \tag{C04}
\]

The path is not supplied as a coordinate lattice: its order is the physical
predecessor/successor order of the declared custody mission. This is still a
constructed mission law, not a claim that generic lineage selects a path.

## 3. Extensive reversible KEEP/BREAK and exact serial support

For every edge choose a prospectively fixed route bit `z_i`. Route in parallel:

\[
z_i=K:\ \operatorname{SWAP}_{L_iK_i},
\qquad
z_i=B:\ \operatorname{SWAP}_{L_iG_i}.                \tag{C05}
\]

All swaps act on disjoint triples and commute. Full KEEP takes `z_i=K` for all
`i`; full BREAK takes `z_i=B` for all `i`. On the formation arm,

\[
\sum_{i=1}^{M}(L_i+K_i+G_i)=M                         \tag{C06}
\]

in every route pattern. Full BREAK therefore removes all active support while
reversibly preserving every relation excitation in explicit quarantine.

Use the common finite link Hamiltonian

\[
H_{\rm link}=-h\sum_{i=1}^{M}
 |1\rangle\!\langle1|_{K_i}\otimes X_{a_i}.           \tag{C07}
\]

The terms commute. After the same `pi hbar/(2h)` pulse, each blank link obeys

\[
a_i=K_i.                                              \tag{C08}
\]

Hence the state-sector transition support is exactly the physical active
relation pattern. Formation/full KEEP activates all `M` path edges;
formation/full BREAK and either sham activate none. Any partial route pattern
deletes exactly its quarantined edge subset.

The complete registered active query has the finite alphabet
`(K_1,...,K_M,a_1,...,a_M) in {0,1}^{2M}`. Every outcome remains present even
when its ideal probability is zero.

## 4. Exact factor-resource scaling

The serial family uses

\[
M+1\ \text{carrier qutrits},\qquad
4M\ \text{relation/link bits},                        \tag{C09}
\]

plus a constant number of source/reference factors and `O(M)` explicitly
declared pulse/history factors. Thus the number of physical tensor factors,
writer operations, route operations, and active query coordinates is `O(M)`.
The Hilbert-space dimension is of course exponential in the number of factors;
no polynomial state-vector simulation claim is made.

Physical custody makes the serial writer depth `M`. The disjoint route and link
pulses have depth one in the ideal schedule. No `N^2` potential-edge array is
allocated.

## 5. Symmetric two-stage branch/merge cell

Use four physical node identities `A,B,C,D`, but give the branch node two source
slots `(A_1,A_2)` and the merge node two recipient slots `(D_1,D_2)`. Intermediate
nodes `B,C` each have one slot. Let `rho_12` be an arbitrary joint state on the
two content subspaces; it may be entangled.

The formation source places `rho_12` at `(A_1,A_2)` and blanks `B,C,D_1,D_2`.
The matched sham puts the same `rho_12` at `(D_1,D_2)` and blanks the upstream
slots. Use two stages:

\[
U_{\rm branch}=U_{A_1\to B}U_{A_2\to C},
\qquad
U_{\rm merge}=U_{B\to D_1}U_{C\to D_2}.              \tag{C10}
\]

Within each stage the two writer unitaries have disjoint support and commute.
Linearity of (C01) gives

\[
\begin{aligned}
F:&\quad \rho_{12}^{A_1A_2}
\longmapsto \rho_{12}^{D_1D_2}
 \otimes|1111\rangle\!\langle1111|_L,\\
S:&\quad \rho_{12}^{D_1D_2}
\longmapsto \rho_{12}^{D_1D_2}
 \otimes|0000\rangle\!\langle0000|_L,                \tag{C11}
\end{aligned}
\]

with every other carrier slot blank. Arbitrary joint coherence and entanglement
are preserved exactly. Nothing is copied and the word “merge” refers to custody
arriving at the same physical node identity in two declared slots, not fusion of
two quantum states into one slot.

The four formed pair memories are

\[
E_\ell=\{AB,AC,BD,CD\}.                               \tag{C12}
\]

Their undirected support is the four-cycle `A-B-D-C-A`, with cycle rank
`|E|-|V|+1=1`. This is the minimal symmetric cell in the declared class with two
disjoint writer operations in each of two equal-depth custody stages. A lawful
asymmetric branch/rejoin can instead have triangular support, so no absolute
minimality claim is made. The constructed cell is a finite relational cycle
witness, not a spatial plaquette or a claim of two-dimensional geometry.

Applying (C05)--(C08) on all four relations gives a four-edge full KEEP and a
four-edge full BREAK with

\[
\sum_{e\in E_\ell}(L_e+K_e+G_e)=4                    \tag{C13}
\]

in both formed intervention arms. Partial route vectors produce the exact
corresponding support subgraph.

## 6. Complete ports and simultaneous content-blind custody

For each fixed serial size or branch/merge cell, the device census contains:

- the common single- or two-carrier content source and authorized formation/sham
  preparation map;
- every carrier slot and fresh writer bit;
- every active and quarantine relation destination;
- every link qubit;
- the prospectively fixed route vector and all pulse controllers/clocks/work
  references;
- one positive source-off hold after the last custody write;
- the complete active query and finite controller history; and
- an arbitrary external complement with exactly factorized future dynamics.

No writer, relation bit, quarantine, controller, failure branch, or query outcome
is removed from the device state. The registered active query acts as identity on
quarantine. Calibration missions may query quarantine, but such a port is not
silently introduced into the active scored mission.

The two branch writers and the two merge writers are simultaneously content
blind: their products commute with independent content-basis exchange on both
tokens and transport an arbitrary joint density operator unchanged. The exact
no-cloning and dimension obstructions above delimit what “simultaneous” can mean.

## 7. Existing `REC`, FHBC, DCL, and Coverage-U compose

For each fixed finite `M`, define the serial event to be formation versus
prepositioned sham, the nominated relation lineage to be
`(L_1,...,L_M)->(K_1,...,K_M)`, and the closure regime to begin after the last
writer. For the branch/merge cell use the analogous four-edge lineage.

The existing record clauses are met exactly:

1. all relation targets are initially blank;
2. only actual custody writes them, while the matched sham has the same final
   carrier state;
3. a positive source-off hold retains the relation vector;
4. full KEEP gives deterministic unit active-query contrast;
5. extensive BREAK moves the whole relation vector to quarantine and makes the
   active formation/sham laws equal;
6. the query is complete, label blind, and cannot recreate a blank active
   lineage; and
7. all causal routes and garbage are explicit.

Thus these are exact model-conditional derivative records:

\[
\operatorname{REC}(r_{\rm serial,M}),
\qquad
\operatorname{REC}(r_{\rm diamond}).                 \tag{C14}
\]

Each fixed mission also satisfies FHBC H0--H4: its device is finite, its bounded
piecewise Hamiltonian schedule is exactly closed, one joint root contains all
content correlations and blanks, only the authorized source receives `F/S`, all
later laws are arm-common for fixed route context, the terminal instrument is
complete, and write precedes query on a nonempty state-wire route. Therefore the
imported theorem gives

\[
\begin{aligned}
\operatorname{FHBC}(r_{\rm serial,M})
&\Longrightarrow DCL_{\rm phys}(r_{\rm serial,M}),\\
\operatorname{FHBC}(r_{\rm diamond})
&\Longrightarrow DCL_{\rm phys}(r_{\rm diamond}),     \tag{C15}
\end{aligned}
\]

and (C14)--(C15) give the existing per-record `COV_union` result. No new record
axiom and no universal U-DCL premise are used.

## 8. Exact boundary

The composed result earns finite physical-model formation of the serial support
used conditionally in `GRA-CH`, plus one lawful finite branch/rejoin cycle. It
also supplies an extensive, reversible whole-lineage support BREAK with linear
factor resources.

It does not prove an open thermodynamic family of physically stabilized pair
registers, resilience to missed or false relations, autonomous source/controller
ports, generic graph selection, carrier transport through the later staggered
link phase, finite-dimensional locality above one dimension, a common continuum
cone, tensor gravity, universal stress coupling, nonlinear closure, or `G`.

The next physics question is whether one same-parent local formation law can
stabilize repeated branch/rejoin cells against relation errors without specifying
a target mesh. The no-cloning and single-slot-merge obstructions must remain in
that search.
