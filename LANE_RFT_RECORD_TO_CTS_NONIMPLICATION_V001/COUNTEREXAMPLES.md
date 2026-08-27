# Exact witnesses and boundaries

## 1. Compact stochastic CTC record

The countermodel is defined in a classical Deutsch-style **distributional
fixed-point CTC theory** and uses two physical interfaces between the same
writer and query regions:

\[
 W\overset{B}{\longrightarrow}Q,
 \qquad
 Q\overset{C}{\longrightarrow}W
       \overset{C}{\longrightarrow}Q.
\]

The first channel is the ordinary record carrier. The second is an independent
closed causal controller. Its backward transition is

\[
K(0,\cdot)=(7/8,1/8),\qquad
K(1,\cdot)=(3/8,5/8).
\]

For a loop law with \(p=P(C=1)\),

\[
 p'=(1-p)\frac18+p\frac58=\frac18+\frac12p.
\]

The distributional consistency equation \(p'=p\) has the unique solution
\(p=1/4\). The contraction coefficient \(1/2\) makes uniqueness explicit. All
interface probabilities are nonnegative and normalized. This is not an ordinary
acyclic sample-path claim that one realized bit obeys both legs simultaneously;
the fixed object is the loop marginal law.

A common blank \(\bot\) is written to \(B=X\), and \(B\) remains unchanged.
The query returns \(H=(B\oplus C,C)\). The event is recovered as the parity of
the two complete outputs. There is a real write, stable hold, positive read,
complete history, and no direct label input at the query.

The reverse controller edge prevents every writer-before-query physical
macro-cut. This is stronger than an unresolved infinite process: the model is
finite-state, compact, normalized, and has a unique consistency law.

### Why the nuisance loop is not removable

The loop bit enters the registered query map and is reported in the complete
history. It is therefore a physical query ancestor, not an unobserved decorative
variable. More importantly, the controller channel physically joins the named
query region back to the named writer region. Deleting it or replacing the
query with a later external log changes the physical episode fixed by D0.

The loop carries no event label, so it does not bootstrap the record. The new
distinction is created only by \(\bot\mapsto X\) in the forward carrier. This
separates record novelty from chronology of every nuisance ancestor.

## 2. What additional axiom excludes the witness

Any one of the following standard physical assumptions excludes the particular
countermodel:

1. stable causality of the complete ancestor domain;
2. global hyperbolicity of the relevant spacetime mission;
3. an authenticated acyclic external incidence frontier between all writes and
   all registered queries; or
4. chronology protection strong enough to forbid physical closed causal
   channels.

RECORD_FM contains none of them. Adding one is legitimate physics, but it is
not a consequence of the word “record.”

## 3. Why a mere hidden bypass does not refute existential CTS

A genuine carrier \(R=X\) can coexist with a direct source broadcast \(B=X\).
A proposed cut containing only \(R\) is incomplete, but a larger cut containing
\((R,B)\) may satisfy T4 and T5. This refutes response-only authentication of
the small cut, not existence of every complete cut.

For the same reason, mutually compatible component marginals do not prove T4,
but an actual standard classical or quantum global state often supplies the
needed joint separator. No standard globally hyperbolic T4 counterexample is
claimed here.

## 4. Positive transposition on an isolated qubit

The transpose map is positive and trace preserving on \(M_2\). It preserves
the basis projectors and therefore preserves a computational-basis record.
However, its Choi operator is the swap, which is negative on the antisymmetric
subspace. It fails complete positivity under an entangled ancillary extension.

Three different conclusions must not be conflated:

1. **Basis-only mission:** a classical diagonal K cover exists; no OPCC
   counterexample.
2. **Tomographically complete isolated mission:** the actual transpose is
   distinguished from identity and lies outside standard CP sewing, but a claim
   that no faithful classical cover exists needs additional contextuality or
   locality constraints.
3. **Entangled-ancilla mission:** partial transpose is nonpositive, so the map is
   not a coherent deterministic operation of standard quantum mechanics.

The second case is a possible countermodel only in a nonstandard physical theory
whose composites are restricted enough to permit positive non-CP dynamics. The
current all-physical record definition does not explicitly exclude such a
theory, but the basis construction alone cannot establish failure of OPCC.

## 5. Why tomography alone is not a no-classical theorem

A finite tomographically complete set can identify the affine reflection
\(s_y\mapsto-s_y\). It rules out calling the hold an identity bit channel. A
finite classical score-table simulator nevertheless always exists, and a
contextual hidden-state packaging can retain preparation or measurement context.

To rule that out physically rather than verbally, the mission must freeze which
operationally equivalent preparations and effects are the same physical
referent, require an encoding to preserve convex mixing and those equivalences,
and supply a violated noncontextuality inequality. Alternatively, a Bell witness
plus authenticated spacelike separation rules out a local classical causal
cover. Neither addition is part of the basis-only transpose record.
