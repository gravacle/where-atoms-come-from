# Faithful admission theorem for bona-fide finite-mission records

**Claim class:** exact domain reduction and physical-cover-to-admission theorem

**Not claimed:** derivation of a complete physical cover from response data,
unbounded all-future coverage, empirical authentication of inaccessible records,
physical reset accessibility, or an unconditioned quantum actualization law

## 1. Independently specified finite-mission domain

First freeze a candidate record claim

\[
 r=(\mathcal O,\mathcal D,X^\circ,Z,\mathcal I,\mathcal J,Q,
       \mathcal R,\mathcal W,\mathcal K,\varepsilon).       \tag{1}
\]

The entries identify a physical observer or patch, the event-to-query mission,
the declared event or innovation alternatives, matched contexts, allowed source
interventions, registered instruments and one label-blind query, the nominated
physical carrier lineage, the acquisition/write region, the declared closure
regime, and the resolution/error contract.  The specification is admitted to
the candidate set \(\mathfrak C_{\rm FM}\) only under the following clauses.

**D0 -- physical referents fixed independently.**  The source, target, writer,
query, clocks, raw outputs, and claimed physical roles are fixed before the
record response is scored.  Renaming a fitted latent coordinate as a carrier is
not a physical specification.

**D1 -- finite operational mission and resolution.**  The event-to-query region
and time interval are bounded.  At the declared resolution, \(X^\circ\), the
context set, the registered protocol family, and each complete registered
outcome history alphabet are finite.  Continuous quantities enter through a
prospectively fixed finite partition.  No one universal bound on size is imposed.

This is operational finiteness.  It does not say that the underlying state space
is finite-dimensional, that the full physical ancestor process has a finite
exact module cover, or that every external port has been discovered.

**D2 -- independently stated formation question.**  The claim states which
counterfactual distinction is alleged to be newly acquired at the nominated
physical coordinate, or which inherited distinction is alleged to enter a new
derivative coordinate.  It states the no-new-write reference, the acquisition
region, and whether closure means survival after source exclusion or only a
declared source-supported live interval.  An extant record whose formation lies
outside \(\mathcal D\) is not silently promoted to a formation episode.

**D3 -- complete registered response law.**  Every failure, no-click,
saturation, censoring, destructive-read result, and adaptive history in a frozen
protocol remains in the registered outcome alphabet.  No positive contrast is
created by conditioning away an outcome.

**D4 -- representation independence.**  Candidate specification and the
physical record predicate below may use source custody, controls, ancestry,
interventions, and physical theory.  They may not use membership in one of the
four admission classes, a successful factorization, quotient dimension, rank,
decoder success, or the existence of an encoding constructed later.

For a candidate \(r\), let \(P_r(h\mid x,z,j)\) denote its population law for
the complete registered history \(h\).  The notation `arm(x)` replaces
`do(x)` when \(x\) is an endogenous alternative rather than an authorized
intervention.  A necessary finite resolved read condition is

\[
 \Delta_r:=\max_{z,j,x\ne x'}
 D_{\rm TV}\!\left(P_r(\cdot\mid x,z,j),
                   P_r(\cdot\mid x',z,j)\right)>0.        \tag{2}
\]

Define \({\rm RECORD}_{\rm FM}(r)\) by the RFT physical invariant, not by (2)
alone: an interaction in \(\mathcal W\) creates the declared new
counterfactual predictive distinction in \(\mathcal R\); that distinction
survives under the declared closure regime to a registered label-blind query;
and neither the writer nor reader recreates the event label.  Novelty,
lineage, closure, and noncreation are physical clauses.  Equation (2) is their
registered positive consequence, not a substitute for them.

The independently specified bona-fide finite-mission domain is

\[
 \mathfrak R^{\rm bf}_{\rm FM}
 =\{r\in\mathfrak C_{\rm FM}:{\rm RECORD}_{\rm FM}(r)\}. \tag{3}
\]

The definition contains no preferred material, lifecycle catalogue, algebra, or
constructor.  Whether an empirical packet establishes membership is a separate
evidential question.

## 2. Faithful physical encoding and independent authentication

An encoding \(\Xi\) of \(r\) is **physically faithful at the declared
operational resolution**, written
\({\rm PhysEnc}_{\rm FM}(r,\Xi)\), when E0--E3 hold.

**E0 -- referent preservation.**  Encoded source, acquisition region, carrier
boundary, writer exclusion, and query refer to the physical systems in (1).
No new wire carrying \(x\) is invented merely to reproduce the score table.

**E1 -- exact registered-law preservation.**  For every frozen source setting,
context, registered protocol, and complete history,

\[
 P_\Xi(h\mid x,z,j)=P_r(h\mid x,z,j).                    \tag{4}
\]

The equality is at the population law of the declared finite resolution; the
empirical error contract concerns estimation of that law.

**E2 -- causal and lineage preservation.**  Every physical event-bearing route
from an authorized source to a registered query is represented.  Every crossing
memory, controller, reference, edge/central degree, classical history, and joint
correlation in the registered ancestors is retained in the boundary object or
inside an exact bounded module.  No physical bypass is deleted and no
response-defined carrier is added.

**E3 -- intervention and provenance preservation.**  The encoding has the same
authorized source operations and future settings.  After the declared cut, the
future has no event-label argument except through the complete physical boundary.
All registered outcomes remain complete.

**E4 -- outcome-independent physical authentication.**  Physical identities,
module support, ports, maps, cut completeness, source custody, and label
provenance are authenticated independently of the scored record contrast and
frozen before it.  Internal mathematical syntax is not its own external physical
certificate.

Define

\[
 {\rm AuthEnc}_{\rm FM}(r,\Xi)
 \quad\Longleftrightarrow\quad
 {\rm PhysEnc}_{\rm FM}(r,\Xi)\land E4.                  \tag{4a}
\]

Physical faithfulness is deliberately stronger than equality (4).
Authentication is stronger again: a physical representation can be correct in
fact even when no surviving evidence can independently certify that fact.

## 3. Theorem FM1 -- every finite law has a weak score-table simulator

For every normalized law in D1--D3, there is a finite classical stochastic DAG
that reproduces (4).  Let a source node output \(X=x\), let arm-common setting
nodes output \((Z,J)=(z,j)\), and let the final node draw the complete history
from the kernel

\[
 K(h\mid x,z,j):=P_r(h\mid x,z,j).                       \tag{5}
\]

Normalization of \(P_r\) makes \(K\) a normalized finite stochastic kernel.
The constructed DAG has finitely many nodes and exact response equality.

### Proof

The finite alphabets in D1 make all nodes and the history kernel finite.  For
fixed \((x,z,j)\), summing (5) over the complete history alphabet gives one by
normalization.  Therefore the joint DAG law conditioned on its three input nodes
is exactly \(P_r\).  QED.

FM1 is only a weak response representation.  Its apparent cut carries \(X\)
itself.  Unless that wire and every path are independently identified in the
physical episode, E0, E2, and E4 do not follow.  Calling FM1 universal record
coverage would make the theorem circular.

## 4. Ontic and independently authenticated complete causal covers

The remaining physical bridge can be stated without using an admission verdict.
An **ontic complete causal cover** \({\cal C}\) for \(r\), written
\({\rm OCC}_{\cup}(r,{\cal C})\), satisfies C0--C4 as physical facts, whether or
not those facts remain empirically certifiable.

**C0 -- same episode.**  The cover maps its source, modules, boundary ports,
carrier lineage, and registered query to the physical referents of (1).  Its
maps reproduce (4) because they are the actual physical maps for that
episode, not because a latent model was selected to fit (4).

**C1 -- finite exact causal cover.**  Finitely many prospectively bounded exact
physical modules cover every source-to-registered-query ancestor in the mission.
Their causal support, inputs, outputs, memory, and external ports are explicitly
contained in the cover.
No unresolved infinite ancestor chain, unbounded feedback process, or unknown
label-bearing port is renamed one module.

**C2 -- complete definite event-to-query boundary.**  A macro-partition \(G\)
places every authorized event occurrence on its past side and every registered
query on its future side.  The crossing object \(C_G\) contains every target,
memory, controller, reference, center/edge degree, outcome history, and joint
correlation crossing the partition in the registered ancestors.  Internal
indefinite order is allowed only inside a bounded higher-order module wholly on
one macro-side.  Every remaining future input has one arm-common joint rule.

**C3 -- source and label provenance.**  The actual source-generating mechanism
and every allowed label recipient are specified as physical facts.  No post-cut
component, instrument, analyst, or external port receives \(x\) outside \(C_G\);
complete outcomes prohibit hidden label-conditioned postselection.

**C4 -- one existing composability branch.**  The physical modules and
interfaces obey at least one of:

1. **K:** a finite classical standard-Borel stochastic DAG, after finite feedback
   unrolling, or a finite-dimensional finite-depth quantum circuit/causally
   ordered comb with normalized kernels or complete CP instruments;
2. **W:** a finite acyclic slot cover by unital von Neumann boundary algebras,
   normal states, normal UCP maps, and complete normal CP instruments;
3. **A:** a finite acyclic slot cover by unital C*-boundary algebras, arbitrary
   algebraic states, UCP maps, and complete CP instruments; or
4. **I:** finite-dimensional deterministic completely CP-preserving
   higher-order modules, complete CP instruments, and the definite external
   macro-cut in C2.

Classical registers and constraint/center data may be included in the relevant
algebraic boundary.  C4 is a physical composability statement, not an inference
from endpoint positivity.

**C5 -- independent authentication.**  C0--C4, including completeness of the
joint port/correlation census and correctness of the module maps, are frozen from
physical architecture, independently calibrated dynamics, custody, and causal
support evidence before scored record outcomes.  Positive contrast, a fitted
quotient, or the desired universal theorem may not establish the cover.

Define the stronger evidential predicate

\[
 {\rm IACC}_{\cup}(r,{\cal C})
 \quad\Longleftrightarrow\quad
 {\rm OCC}_{\cup}(r,{\cal C})\land C5.                   \tag{5a}
\]

An optional strengthened flag \({\rm GRAPH}_{\rm auth}({\cal C})\) supplies a
finite physical incidence graph containing every source-to-query morphism and
identifies all crossing incidences at \(G\).  It is not included silently in
base `IACC` because the algebraic and higher-order admission lanes distinguish
complete separator authentication from A5's separately exposed graph route.

## 5. Theorem FM2 -- faithful transcription into the four-class union

Define

\[
\begin{aligned}
 {\rm ADM}_{\cup}(\Xi):={}&{\rm FCN\_ADMIT}(\Xi)\lor
 {\rm ACN\_CP\_ADMIT}(\Xi)\\
 &\lor{\rm CSTAR\_CP\_ADMIT}(\Xi)\lor
 {\rm FICO\mbox{-}DM\mbox{-}ADMIT}(\Xi).                \tag{6}
\end{aligned}
\]

Let \({\rm ARCH}_{\cup}(\Xi)\) denote the same four physical mathematical
architecture branches with their structural completeness clauses true in fact,
but without asserting that an independent empirical certificate survives.

For every independently specified finite-mission candidate, including a
zero-contrast candidate,

\[
 \boxed{\begin{aligned}
 {\rm OCC}_{\cup}(r,{\cal C})
 &\Longrightarrow
 \exists\Xi\,[{\rm PhysEnc}_{\rm FM}(r,\Xi)
               \land{\rm ARCH}_{\cup}(\Xi)],              \tag{7a}\\
 {\rm IACC}_{\cup}(r,{\cal C})
 &\Longrightarrow
 \exists\Xi\,[{\rm AuthEnc}_{\rm FM}(r,\Xi)
               \land{\rm ADM}_{\cup}(\Xi)].               \tag{7b}
\end{aligned}}
\]

Under (7a), A1--A4 and the Coverage-U factorization hold as objective physical
statements.  Under (7b), they additionally have the existing independently
authenticated admission packet.  If \({\rm RECORD}_{\rm FM}(r)\) holds, its
registered event contrast is positive.  If
\({\rm GRAPH}_{\rm auth}({\cal C})\) also holds, the complete full cut earns
A5's authenticated-graph route.

### Proof by physical branch

**K.**  Copy the physical nodes or circuit slots, normalized kernels or CP maps,
complete outcomes, ports, source provenance, and \(G,C_G\) from \({\cal C}\)
into the FCN packet.  C1 gives finite architecture; C4(K) gives valid local
sewing; C2 gives the complete joint boundary; C3 gives no post-cut label
injection.  C0--C4 give the ontic K architecture.  C5 additionally gives
outcome-independent physical authentication and hence `FCN_ADMIT`.

**W.**  Copy the bounded slots, von Neumann boundary algebras, normal source
states, normal UCP maps and instruments, complete separator, ports, and
provenance into the ACN packet.  C0--C4 are exactly the physical content required
for the structural normal-W* branch; C5 upgrades the packet to
`ACN_CP_ADMIT`, without selecting an algebra from response contrast.

**A.**  Copy the bounded slots, common C*-observable system, arbitrary algebraic
states, UCP maps and instruments, complete separator, ports, and provenance into
the C*-packet.  C0--C4 give the ontic C*-UCP architecture; disjoint states need
no common selected physical folium.  C5 supplies the independent authentication
required by `CSTAR_CP_ADMIT`.

**I.**  Copy the bounded completely CP-preserving higher-order modules, complete
instruments, open ports, macro-order, complete boundary, and provenance into the
FICO packet.  Internal indefinite order is wholly within a macro-side by C2, so
the true definite macro-cut gives the ontic higher-order architecture.  C5
upgrades it to the authenticated `FICO-DM-ADMIT` packet.

In every branch C0 gives E0 and exact law equality E1; C2 supplies E2; and C3
supplies E3.  Thus C0--C4 prove the ontic physical transcription (7a).  When C5
also holds, it supplies E4 and the outcome-independent authentication clause of
the corresponding existing admission predicate, proving (7b).  Neither
transcription is the synthetic construction of FM1.

For (7b), each cited admission theorem already proves A1--A4, so the existing CU
factorization follows with an authenticated packet.  Under `GRAPH_auth`, deleting
all crossing incidences removes every finite source-to-query path by the
first-crossing argument, earning A5's graph route.  No discard/reprepare operation
follows.  QED.

## 6. Theorem FM3 -- ontic universal-quantifier reduction

The ontic physical target, with no claim of universal certifiability, is

\[
\begin{aligned}
 {\rm OUFA}_{\rm FM}:\quad
 \forall r\in\mathfrak R^{\rm bf}_{\rm FM}\ \exists\Xi\,[
 {\rm PhysEnc}_{\rm FM}(r,\Xi)\land{\rm ARCH}_{\cup}(\Xi)].
                                                               \tag{8}
\end{aligned}
\]

Let the corresponding ontic physical-cover statement be

\[
 {\rm OPCC}_{\rm FM}:\quad
 \forall r\in\mathfrak R^{\rm bf}_{\rm FM}\ \exists{\cal C}\,
 {\rm OCC}_{\cup}(r,{\cal C}).                             \tag{9}
\]

Then, for the definitions in this lane,

\[
 \boxed{{\rm OUFA}_{\rm FM}\quad\Longleftrightarrow\quad
        {\rm OPCC}_{\rm FM}.}                             \tag{10}
\]

### Proof

FM2(7a) applied record by record proves (9) implies (8).  Conversely, take the
faithful physical encoding supplied by (8).  Forget only its mathematical packet
packaging while retaining the physical referent map required by E0, exact bounded
modules and maps, complete boundary and event-bearing paths required by E2, and
provenance required by E3.  The surviving physical objects satisfy C0--C4 in the
corresponding K, W, A, or I branch, yielding the cover in (9).  QED.

Equation (10) is a quantifier reduction, not a proof of (9).  The remaining
ontic physical proposition is precisely that every bona-fide finite-mission
record has one true complete cover in the existing standard composability union.
That proposition remains open.  It does not license a URM `PASS` when the physical
facts cannot be independently certified.

The superficially stronger authenticated all-record statement

\[
 \forall r\in\mathfrak R^{\rm bf}_{\rm FM}\ \exists\Xi\,
 [\,{\rm AuthEnc}_{\rm FM}(r,\Xi)\land{\rm ADM}_{\cup}(\Xi)\,]
                                                               \tag{10a}
\]

is not the ontic target.  It is false on the domain (3) when that domain includes
unique inaccessible natural records, as FM5 proves below.

## 7. Theorem FM4 -- response-only faithful authentication is impossible

No rule whose inputs are only the finite registered response law and the visible
target labels can be both:

1. **sound:** every certified visible boundary is the complete physical
   event-to-query boundary; and
2. **complete:** it certifies such a boundary for every bona-fide finite response
   record.

### Exact observational twins

Let \(X\in\{0,1\}\).  Both worlds expose the visible target \(C=X\) and the
registered query \(Q=X\):

\[
\begin{array}{c|ccc}
 &C&H&Q\\ \hline
 \text{honest world }H_0&X&0&C\\
 \text{bypass world }H_1&X&X&H.
\end{array}                                               \tag{11}
\]

The complete registered law of \((C,Q)\) is identical in both worlds.  Under a
common replacement \(C\leftarrow0\), however,

\[
 D_{\rm TV}(Q_1,Q_0)=0\quad\text{in }H_0,
 \qquad
 D_{\rm TV}(Q_1,Q_0)=1\quad\text{in }H_1.                \tag{12}
\]

Thus the visible boundary \(\{C\}\) is complete in \(H_0\) and bypassed in
\(H_1\).  A rule receiving the common registered law must return the same answer
in both worlds.  If it certifies \(\{C\}\), it is unsound in \(H_1\); if it
withholds the certificate, it is incomplete for \(H_0\).  Naming an unobserved
\(H\) is not physical authentication.  An authenticated graph/port census or a
valid spanning intervention family distinguishes the worlds; the response law
alone does not.  QED.

Endpoint response laws also do not certify C4.  On \(M_2\), identity and
transposition agree on every diagonal state and diagonal effect.  The transpose
map is nevertheless not completely positive: its Choi matrix is the swap
operator \(F\), and for \(v=|01\rangle-|10\rangle\),

\[
 v^*Fv=-2<0.                                               \tag{13}
\]

Therefore a finite endpoint test family can fail to distinguish a composable CP
module from a merely positive non-CP rule.  Independent module physics or a
complete composability certificate is indispensable.

FM4 does not say that physical authentication is impossible.  It proves that it
cannot be manufactured from the same finite record-response law whose causal
carrier is at issue.

## 8. Theorem FM5 -- all-ontic authenticated universality is false

Consider a unique finite bounded natural episode.  In the record world \(M_R\),

\[
 X:=U_X,\qquad R:=X,\qquad Q:=R.                          \tag{14}
\]

The transient event physically writes the surviving carrier \(R\), and the later
query reads it without creating the distinction.  Thus
\({\rm RECORD}_{\rm FM}(r)\) is ontically true.  By construction, after the
episode the event, writer, every independent provenance trace, every intervention
opportunity, transported calibration, netlist/custody record, and negative-control
family is destroyed or inaccessible.  Only the unique archive and its read
remain.

Now define the latent-common-cause world \(M_C\),

\[
 U\sim P_X,\qquad X:=U,\qquad R:=U,\qquad Q:=R,            \tag{15}
\]

with no \(X\to R\) write edge.  The accessible laws agree exactly:

\[
 P_{M_R}(X,R,Q)=P_{M_C}(X,R,Q).                           \tag{16}
\]

Any outcome-independent certification rule using the permitted surviving
evidence returns the same result in both worlds.  Accepting both falsely
authenticates the \(X\to R\) lineage in \(M_C\); rejecting both leaves the
genuine record in \(M_R\) without E4.  Since the construction removes every
independent certificate by premise, \(M_R\) has no
\({\rm AuthEnc}_{\rm FM}\) even though it may have a correct
\({\rm PhysEnc}_{\rm FM}\).

Therefore

\[
 \boxed{\neg\!\left[
 \forall r\in\mathfrak R^{\rm bf}_{\rm FM}\ \exists\Xi\,
 ({\rm AuthEnc}_{\rm FM}(r,\Xi)\land{\rm ADM}_{\cup}(\Xi))
 \right]}
                                                               \tag{17}
\]

when the domain includes unique inaccessible natural records.  This is a
counterexample to universal independent **certifiability**, not to ontic physical
separator or representation existence.  A submitted packet lacking C5 must
return `INDETERMINATE` or `OUTSIDE_IACC_DOMAIN`, never `NO_RECORD`.  QED.

## 9. Exact scientific status

Proved:

1. a representation-independent finite-mission record domain;
2. universal weak finite score-table representability;
3. ontic faithful transcription from a true complete standard physical cover;
4. faithful authenticated admission into the existing four-class union from an
   independently authenticated complete physical cover;
5. the exact reduction of ontic universal faithful coverage to ontic complete
   physical-cover existence;
6. response-only impossibility for sound-and-complete cut authentication, plus an
   exact non-CP composability witness; and
7. refutation of universal independently authenticated admission over the
   all-ontic domain containing unique inaccessible records.

Open:

1. `OPCC_FM`, the ontic complete-cover proposition, for every bona-fide
   finite-mission physical record;
2. physical records with no definite event-to-query macro-cut or
   no exact bounded module cover in an existing branch;
3. exact gravitational/background-dynamical records not known to have a common
   fixed observable system or compact separator;
4. A5 physical replacement and finite-sample confirmation where claimed; and
5. the unconditioned law that produces one objective endogenous quantum outcome.

The refuted authenticated universal statement must not be relabeled `open`.
The surviving open target is ontic physical coverage.  The certification theorem
is record-by-record: C5 plus C0--C4 yields an authenticated packet, while missing
evidence yields an indeterminate/outside-domain verdict.

The open ontic item is not a request for another material taxonomy.  A new
admission theorem is warranted only if a concrete bona-fide finite-mission record
has a physically justified complete cover that fails all four existing
composability branches.
