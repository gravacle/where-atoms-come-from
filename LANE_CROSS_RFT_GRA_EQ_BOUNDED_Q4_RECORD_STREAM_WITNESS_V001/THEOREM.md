# Bounded q=4 record-stream witness theorem

**Lane ID:** `CROSS-RFT-GRA-EQ-BQ4RSW-V001`

**Official short name:** `BQ4RSW`

**Date:** 2026-08-27

**Claim class:** exact finite reversible model witness for one complete-port
four-operation compatible merger; exact reference-stable KEEP quotient; exact
same-logical-hardware selector BREAK; exact finite-cap and `S_4` covariance
theorem

**Builder status:** `SOURCE_FROZEN_PENDING_INDEPENDENT_HOSTILE_AUDIT`

**Not claimed:** that nature supplies these four streams; a derivation or
selection of four operation types or three spatial dimensions; an unbounded
total-live sector; a measured apparatus; a scalable low-resource
implementation; a physical causal-simplex binding; coherent lateral
propagation; a connection, continuum, curvature, gravity, or `G`

## 1. Exact inherited boundary

This lane uses three pinned results only. EG proves that a supplied compatible
`q`-operation complete-front merger gives a count quotient, but current
physics had supplied only its `q=2` case. EK supplies the exact
reference-stable residual-descent and complete-port conditions. EO proves what
follows from a supplied q=4 merger and exact complete-port `S_4` covariance,
while leaving that antecedent open.

This packet does not alter EO. It constructs one finite, unitary-dilatable
record-stream model satisfying the missing merger antecedent through a
declared mission cap. Whether the actual world instantiates the four supplied
streams remains experimental.

## 2. Four supplied record-counter streams

Fix a total mission cap \(R\in\mathbb N\) with \(R\ge4\), alphabet
\(\Sigma=\{1,2,3,4\}\), and depth \(0\le d\le R\). Supply four identical
record-counter live-code streams

\[
 \mathcal S_a^{\rm live}
 =\operatorname{span}\{|r\rangle_a:0\le r\le R\},
 \qquad a\in\Sigma.                                  \tag{EQ01}
\]

The state \(|r\rangle_a\) means that stream \(a\) retains \(r\)
authenticated, distinguishable, identical-format record marks. Every mark
must independently pass the inherited record predicate; this theorem does
not make a counter increment into a record by definition. The four local
stream factors, mark formats, and append circuits are identical up to label.
Concretely, \(|r\rangle_a\) may denote the canonical prefix
\(|1\rangle^{\otimes r}|0\rangle^{\otimes(R-r)}\) of \(R\) separately
addressed retained cells; the displayed \((R+1)\)-dimensional space is the
prefix-code subspace, not erasure of the individual cell lineages.
The padded unitary ambient is
\(\widetilde{\mathcal S}_a
=\operatorname{span}\{|r\rangle_a:0\le r\le R+1\}\).
The \(R+1\) level is an explicit guard/off-code level. The cap controller
prevents a live shift at total depth \(R\); cyclic action on the padded
complement is part of the dilation and never a live counter wrap.

At depth \(d\), define the active count front and complete order register

\[
 \mathcal Q_d=\operatorname{span}\{|m\rangle:
 m\in\mathbb N_0^4,\ |m|_1=d\},                    \tag{EQ02}
\]

\[
 \mathcal Z_d=\operatorname{span}\{|w\rangle:w\in\Sigma^d\}. \tag{EQ03}
\]

Here
\(|m\rangle=|m_1\rangle_1\otimes\cdots\otimes|m_4\rangle_4\).
The basis label \(|w\rangle\) is the complete word/provenance code: it includes
every varying route, target-stream, controller, work, environment, and garbage
label carrying more information than the count. Every remaining factor
\(\mathcal P_d\) is explicitly label-invariant blank/allocation scaffolding:
\[
 \mathcal P_d=\mathcal P_0\otimes
 \bigotimes_{j=1}^{d}\mathcal B_j,\qquad
 J_d|p\rangle=|p\rangle|0\rangle_{\mathcal B_{d+1}}, \tag{EQ03a}
\]
with every \(|0\rangle_{\mathcal B_j}\) fixed by \(S_4\). Invariant matched
spectators remain with the untouched reference. No varying port is omitted.

The admitted mission does not use the whole tensor product
\(\mathcal Q_d\otimes\mathcal Z_d\). Its reachable code is
\[
 \mathcal H_d^{\rm code}
 =\operatorname{span}\{
 |m(w)\rangle_Q|w\rangle_Z|p_d\rangle_P:
 w\in\Sigma^d\},                                   \tag{EQ03b}
\]
where \(m(w)\) is the word-count vector and \(|p_d\rangle\) is the invariant
allocated-blank state. It has dimension \(4^d\), starts from the unique empty
codeword, and is preserved by every admitted \(V_{a,d}\). Inconsistent
\(|m\rangle_Q|w\rangle_Z\) pairs and other dilation-complement states are
retained off-code states, not admitted physical histories and not counted as
fronts below.

For \(d<R\), define the isometries

\[
 C_{a,d}|m\rangle=|m+e_a\rangle,\qquad
 W_{a,d}|w\rangle=|wa\rangle.                       \tag{EQ04}
\]

Let \(J_d\) be the common isometry carrying the blind scaffolding and
allocating its fresh slots. The write isometry on nonreference factors is

\[
 V_{a,d}=C_{a,d}\otimes W_{a,d}\otimes J_d:
 \mathcal Q_d\otimes\mathcal Z_d\otimes\mathcal P_d
 \longrightarrow
 \mathcal Q_{d+1}\otimes\mathcal Z_{d+1}\otimes\mathcal P_{d+1}. \tag{EQ05}
\]

Each \(C_{a,d}\) writes one mark into stream \(a\); each \(W_{a,d}\) writes
the same operation label into the next fresh order slot. Because every factor
is injective, \(V_{a,d}\) is an isometry. Its unitary completion is not chosen
arbitrarily. On
\[
 \mathcal H_{a,d}^{\rm dil}
 =(\mathcal Q_d\otimes\mathcal Z_d\otimes\mathcal P_d)
  \mathbin{\dot\oplus}
  (\mathcal Q_{d+1}\otimes\mathcal Z_{d+1}\otimes\mathcal P_{d+1}),
\]
use the canonical self-adjoint swap dilation
\[
 \mathfrak D(V_{a,d})=
 \begin{pmatrix}
 0&V_{a,d}^{\dagger}\\
 V_{a,d}&I-V_{a,d}V_{a,d}^{\dagger}
 \end{pmatrix}.                                    \tag{EQ05a}
\]
It swaps every input with its image and fixes the orthogonal target
complement. Direct multiplication using \(V_{a,d}^{\dagger}V_{a,d}=I\) gives
\(\mathfrak D(V_{a,d})^2=I\). Thus invalid and off-code target states are
retained rather than erased.

Equivalently, a gate-level live implementation uses cyclic shifts on
\((R+2)\)-level counters and one fresh five-level blank/label slot per stage.
On the live code subspace it has no wraparound or erasure.

The complete map and dilation are respectively
\(V_{a,d}\otimes I_{\mathcal R}\) and
\(\mathfrak D(V_{a,d})\otimes I_{\mathcal R}\): every invariant matched
spectator and untouched reference is retained, not included in the blind trace
or omitted.

The same stream-local circuit and mark format recur with one fresh finite slot
at each stage. This is bounded reuse, not a post-hoc identification of four
unrelated boxes.

The lower bound \(R\ge4\) makes all pair, three-letter, and four-letter overlap
contexts nonvacuous. The formulas below remain algebraically valid at smaller
caps, but those caps do not by themselves witness the full overlap audit.

## 3. One common KEEP/BREAK activation block

Supply at each depth a fresh routing-front register

\[
 \mathcal F_d=\operatorname{span}
 \bigl(\{|\bot\rangle\}\cup
       \{|\widehat w\rangle:w\in\Sigma^d\}\bigr),    \tag{EQ06}
\]

where \(|\bot\rangle\) is not the empty-word code, and a prospective selector
\(s\in\{K,B\}\). For each \(w\), let \(X_w\) swap
\(|\bot\rangle\) with \(|\widehat w\rangle\) and fix all other basis states.
Define one common activation unitary

\[
 U_d^{\rm act}
 =|K\rangle\!\langle K|\otimes I_{Z_dF_d}
 +|B\rangle\!\langle B|\otimes
   \sum_{w\in\Sigma^d}|w\rangle\!\langle w|_{Z_d}\otimes X_w. \tag{EQ07}
\]

It obeys

\[
 \begin{aligned}
 U_d^{\rm act}|K,w,\bot\rangle&=|K,w,\bot\rangle,\\
 U_d^{\rm act}|B,w,\bot\rangle&=|B,w,\widehat w\rangle.
 \end{aligned}                                       \tag{EQ08}
\]

Equation (EQ07) is a direct sum of transpositions and is exactly unitary.
It copies only the declared orthogonal record basis; on superpositions it
creates the corresponding entanglement and makes no universal no-cloning
claim.
KEEP and BREAK use the same stream writes, word writes, logical registers,
activation hardware, gate schedule, and cap. Only the prospectively prepared
selector input differs. The ideal logical factors are energy-degenerate, so
the model adds no branch-specific logical energy change or storage count. A
laboratory claim must still match and census physical heat, recoil, supply,
clock, route, and control ports; that equality is not inferred from the
logical circuit.

The selector is registered and fixed for one mission. It is an intervention
label, not a fifth operation or front dimension.

The chronological stage is explicit: first \(V_{a,d}\) maps the word from
length \(d\) to \(d+1\); then \(U_{d+1}^{\rm act}\) acts on that new
\(\mathcal Z_{d+1}\) and a fresh \(\mathcal F_{d+1}\). Thus BREAK copies
\(wa\), not the pre-append word \(w\), into the new current routing front.

## 4. Theorem BQ4RSW-1 -- exact reference-stable KEEP merger

On KEEP, define the active disposition and complete blind residual

\[
 \mathcal A_d^K=\mathcal Q_d\otimes
 \operatorname{span}\{|\bot\rangle_{F_d}\}\otimes
 \operatorname{span}\{|{\rm LIVE}\rangle\},         \tag{EQ09}
\]

\[
 \mathcal B_d^K=\mathcal Z_d\otimes\mathcal P_d
 \otimes\mathcal F_{<d}.                            \tag{EQ10}
\]

The factor \(\mathcal F_{<d}\) contains all archived routing registers.
Matched common spectators remain with the untouched reference; they are not
silently traced. Define the KEEP active/archive isometry
\[
 G_{a,d}^K|m,\bot_d,{\rm LIVE}\rangle
 =|m+e_a,\bot_{d+1},{\rm LIVE}\rangle_{A_{d+1}}
  |\bot_d\rangle_{F_d^{\rm arch}},                  \tag{EQ11}
\]
and the old-blind isometry
\[
 L_{a,d}=W_{a,d}\otimes J_d\otimes I_{F_{<d}}.      \tag{EQ11a}
\]
After the canonical output-factor reordering,
\(\Phi_{a,d}^K=G_{a,d}^K\otimes L_{a,d}\).
Thus the old current routing register is transferred to the blind archive,
not erased, while the new current blank is allocated. Let
\(\bar\Phi_{a,d}^K\) be the descended active isometry obtained from (EQ11)
after tracing its fixed one-dimensional archive output.

For every state \(\rho\) on
\(\mathcal A_d^K\otimes\mathcal B_d^K\otimes\mathcal R\) and every untouched
reference \(\mathcal R\),

\[
 \boxed{
 \operatorname{Tr}_{\mathcal B_{d+1}^K}
 \bigl[(\Phi_{a,d}^K\otimes I_{\mathcal R})\rho
       (\Phi_{a,d}^{K\dagger}\otimes I_{\mathcal R})\bigr]
 =(\bar\Phi_{a,d}^K\otimes I_{\mathcal R})
   \operatorname{Tr}_{\mathcal B_d^K}(\rho)
   (\bar\Phi_{a,d}^{K\dagger}\otimes I_{\mathcal R}).} \tag{EQ12}
\]

Before outcomes, freeze the admitted KEEP future family
\(\mathfrak F_K^{(R)}\) to the live generators
\(\{\Phi_{a,d}^K:d<R\}\), the four cap-terminal generators
\(\{\Psi_{a,R}^K\}\), their finite chronological compositions, and identities
on untouched references and matched common spectators. No other channel which
reads \(\mathcal B_d^K\) is admitted. Each cap generator has product form
\[
 \Psi_{a,R}^K=\bar\Psi_{a,R}^K\otimes I_{\mathcal B_R^K}, \tag{EQ12a}
\]
where \(\bar\Psi_{a,R}^K\) changes LIVE to the explicit `CAP,a` disposition
without changing the count or current blank. Hence every generator satisfies
(EQ12), and induction under composition gives
\[
 \operatorname{Tr}_{B_{\rm out}}
 [(\Phi\otimes I_{\mathcal R})\rho(\Phi^\dagger\otimes I_{\mathcal R})]
 =(\bar\Phi\otimes I_{\mathcal R})
   \operatorname{Tr}_{B_{\rm in}}(\rho)
  (\bar\Phi^\dagger\otimes I_{\mathcal R})
 \quad\text{for every }\Phi\in\mathfrak F_K^{(R)}.  \tag{EQ12b}
\]

The cap-attempt channel is terminal, not an append and not a diamond crossing
the cap. Wherever both append orders remain admitted, namely \(d\le R-2\),

\[
 C_{b,d+1}C_{a,d}=C_{a,d+1}C_{b,d}\qquad(a\ne b),   \tag{EQ13}
\]

so one four-stream architecture realizes all six pair relations and all
overlapping live contexts. At \(d=R-1\), one final append may reach the
boundary, but a second append is not admitted; at \(d=R\), an attempted
operation produces only the explicit terminal response in section 6. Its
exact depth-\(d\) front is

\[
 \boxed{
 \{m\in\mathbb N_0^4:|m|_1=d\},\qquad
 F_d^K={d+3\choose3}.}                              \tag{EQ14}
\]

### Proof

Every blind factor in (EQ11) is carried by an isometry. Partial trace is
invariant under an isometry on the traced subsystem, proving (EQ12) on matrix
units and then by linearity, including arbitrary reference entanglement.
In selector state \(K\), every new routing register is the common blank and
can be archived blind. The count maps add standard basis vectors, proving
(EQ13). Adjacent swaps generate every word permutation. Counts are the only
descended labels, while orthogonal count states prevent extra identification.
Stars and bars gives (EQ14). QED.

Words with equal counts remain orthogonal complete states in
\(\mathcal B_d^K\). KEEP is equality of the descended future front, not
equality or erasure of complete histories. This is an explicit bounded model
witness for `Q4-MERGE`, not an actual-universe assertion.

## 5. Theorem BQ4RSW-2 -- same-logical-hardware BREAK restores the free word tree

On BREAK, the same \(V_{a,d}\) writes the same stream marks, word slots, and
provenance, after which \(U_{d+1}^{\rm act}\) copies the new word \(wa\) into
the fresh active routing front. The complete word/provenance register is now
future-feeding rather than blind. Distinct words have orthogonal future-active
outputs:

\[
 w\ne w'\Longrightarrow
 \langle\widehat w|\widehat {w'}\rangle=0.          \tag{EQ15}
\]

Refreshing the current routing front after every append gives

\[
 \widehat w\longrightarrow\widehat{wa},
 \qquad a\in\Sigma.                                 \tag{EQ16}
\]

BREAK therefore removes every interchange identification with no replacement.
On the reachable code (EQ03b), at depth \(d\) and through cap \(R\),

\[
 \boxed{
 F_d^B=4^d,\qquad
 |B_B(o,R)|={4^{R+1}-1\over3}.}                     \tag{EQ17}
\]

The KEEP rooted census is

\[
 \boxed{|B_K(o,R)|={R+4\choose4}.}                  \tag{EQ18}
\]

### Proof

The word code is a recursively future-active output, not a passive read. The
map \(w\mapsto\widehat w\) is injective, and the next stage produces
\(\widehat{wa}\), giving the free four-ary tree. A geometric sum proves
(EQ17); the hockey-stick identity proves (EQ18). QED.

BREAK neither destroys nor reverses a record. The exponential complete memory
already exists under KEEP; BREAK changes its causal role.

## 6. Theorem BQ4RSW-3 -- exact cap and complete disposition

For \(d<R\), every \(A_a\) has one deterministic LIVE append. At \(d=R\), no
counter is incremented and no word slot is overwritten. An attempted \(A_a\)
is a terminal cap query, not an admitted append, and is mapped isometrically to

\[
 |m,w,f_s(w),s,p,{\rm LIVE},0\rangle
 \longmapsto
 |m,w,f_s(w),s,p,{\rm CAP},a\rangle,                \tag{EQ19}
\]

where \(f_K(w)=\bot\), \(f_B(w)=\widehat w\), and \(p\) denotes all complete
provenance/scaffolding labels. The map retains the attempted label, selector,
word, counters, routing front, controller, and all common ports, and has no
outgoing active edge. Calling this cap isometry \(K_{a,R}^s\), its off-code
completion is the same canonical dilation \(\mathfrak D(K_{a,R}^s)\) as in
(EQ05a). Since `CAP,a` and every label-bearing input are permuted together,
this dilation is also exactly `S_4` covariant. The complete census is:

1. the four record streams and count are active under both arms;
2. word/provenance is blind under KEEP and explicitly feeding under BREAK;
3. the current routing front is common blank under KEEP and active under BREAK;
4. the selector is a registered fixed intervention input;
5. LIVE and every `CAP,a` are explicit orthogonal dispositions;
6. route, controller, clock, work, supply, heat, recoil, environment, garbage,
   validity, blank slots, and reference ports are either in complete
   provenance or are matched common spectators; and
7. all invalid and off-code completion states are retained in the explicit
   dilation complement and are either invariant or placed in complete
   label-permutation orbits.

The ideal permutation circuit has no hidden failure, unregistered overflow,
wraparound, or overwritten slot. A nonideal realization must add every failure
and invalid branch; postselecting LIVE would not instantiate the theorem. The
finite active domains are

\[
 \mathcal D_K^{(R)}=\{m\in\mathbb N_0^4:|m|_1\le R\},
 \qquad
 \mathcal D_B^{(R)}=\bigcup_{d=0}^{R}\Sigma^d.      \tag{EQ20}
\]

No finite \(R\) proves an infinite total-live sector.

## 7. Theorem BQ4RSW-4 -- exact complete-port `S_4` covariance

For \(\pi\in S_4\), define

\[
 U_\pi^Q|m_1,m_2,m_3,m_4\rangle
 =|m_{\pi^{-1}(1)},m_{\pi^{-1}(2)},
   m_{\pi^{-1}(3)},m_{\pi^{-1}(4)}\rangle,          \tag{EQ21}
\]

\[
 U_\pi^Z|a_1\cdots a_d\rangle
 =|\pi(a_1)\cdots\pi(a_d)\rangle,\qquad
 U_\pi^F|\widehat w\rangle=|\widehat{\pi w}\rangle. \tag{EQ22}
\]

This fixes \(|\bot\rangle\) and the selector, and permutes streams, operations,
routes, provenance copies, and `CAP,a`. Then

\[
 U_\pi^Q C_{a,d}=C_{\pi(a),d}U_\pi^Q,\qquad
 U_\pi^Z W_{a,d}=W_{\pi(a),d}U_\pi^Z,              \tag{EQ23}
\]

Define \(U_\pi^{P_d}=I_{P_d}\). Then
\[
 U_\pi^{P_{d+1}}J_d=J_dU_\pi^{P_d}.                \tag{EQ23a}
\]
Consequently \(V_{a,d}\) intertwines the full \(QZP\) actions. The direct-sum
action on (EQ05a) obeys
\[
 U_\pi^{\rm dil}\mathfrak D(V_{a,d})
 =\mathfrak D(V_{\pi(a),d})U_\pi^{\rm dil}.         \tag{EQ23b}
\]
The activation unitary commutes with the joint permutation action. Every
`INVALID,a` or `CAP,a` label is mapped to the corresponding
\(\pi(a)\)-label, while invariant blank and generic invalid labels are fixed.

### Proof

The relations follow on basis states. Relabelling \(w\) maps \(X_w\) to
\(X_{\pi w}\); the sum in (EQ07) is invariant. All label-bearing terminal and
provenance ports are permuted rather than discarded. Equation (EQ23a) is
immediate from the invariant blank allocation. Intertwining of \(V_{a,d}\)
then carries both its range projector and orthogonal complement covariantly,
which proves (EQ23b). Thus the off-code completion cannot break the symmetry.
QED.

This covariance makes the four supplied streams identical under relabelling.
It does not explain why there are four. The construction generalizes to other
supplied finite alphabets.

## 8. Relation to EO and ceilings

For each fixed \(R\), the model supplies

\[
 \text{four authenticated record streams}
 \longrightarrow\text{one compatible q=4 descended count front}
 \xrightarrow{\rm BREAK}\text{one active free-word front}.     \tag{EQ24}
\]

It closes EO's first q=4 merger antecedent **inside EO's explicit finite-cap
restriction for this model**. Physical sibling-diamond adjacency exists
through layer \(R-1\), not on the top layer \(R\). Joining EO yields its
finite-layer \(A_3\) contrast and tetrahedral Gram consequences only with
EO's separate calibration and physical-cell binding.

It does not supply actual-world preparation, select \(q=4\), prove an
unbounded sector, protect blindness autonomously, provide coherent sibling
transport or a connection, obtain the full Lorentz cone, or derive curvature,
stress coupling, gravity, or \(G\). The word and routing factors are finite
but have exponentially growing Hilbert dimension with \(R\); this is an exact
finite witness, not a scalable resource theorem.

## 9. Sharp next physical antecedent

\[
 \boxed{
 \begin{gathered}
 \text{Instantiate four identical authenticated record streams in one}\\
 \text{same-parent apparatus; verify reference-stable KEEP, recursive BREAK,}\\
 \text{the complete cap ledger, and all 24 complete-port permutations.}
 \end{gathered}}                                      \tag{EQ25}
\]

A passing experiment instantiates EO's q=4 premise on that bounded domain. A
failure of residual blindness, complete-port symmetry, or recursive BREAK
separation refutes the instantiation without changing the conditional theorem.

## 10. Disposition

**Disposition:**

`BOUNDED_FOUR_IDENTICAL_RECORD_STREAMS_SUPPLY_ONE_COMPLETE_Q4_MERGER_EXACTLY__KEEP_REFERENCE_STABLY_DESCENDS_TO_COUNT_FRONTS__SAME_LOGICAL_HARDWARE_SELECTOR_BREAK_RESTORES_FOUR_ARY_FREE_WORD_FRONTS__FINITE_CAP_OVERFLOW_AND_COMPLETE_PORTS_EXPLICIT__S4_COVARIANCE_EXACT__MODEL_CLOSES_EO_Q4_ANTECEDENT_ONLY_CONDITIONALLY__FOUR_STREAMS_AND_ACTUAL_WORLD_INSTANTIATION_NOT_DERIVED__NO_GRAVITY`
