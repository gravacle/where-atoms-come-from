# Counterexamples and excluded cases

These examples refute tempting weakenings of CTS.  Except where explicitly
stated, they are not counterexamples to `OPCC_FM` itself.

## 1. A finite score table is not a Cauchy cover

For any finite registered law \(P(h\mid x,z,j)\), construct a DAG with a source
node carrying \(x\) and a final node that samples exactly that table.  The model
is finite and normalized, but its event-label wire need not be a physical carrier
or ancestor of the actual episode.

Thus

\[
 \text{finite mission plus exact endpoint fit}
 \not\Longrightarrow \operatorname{CTS}(r).
\]

It fails T0 and does not establish T1--T4.  This is why CTS is a physical premise
fixed independently of response success.

## 2. A convenient carrier marginal can omit an exact bypass

Let \(x\in\{0,1\}\).  The visible proposed boundary variable is \(C=x\), while an
unobserved controller also stores \(H=x\).  The registered future output is
\(Q=H\).  In ordinary operation \(C=Q=x\), so the visible marginal looks like a
perfect carrier.  Under common replacement \(C\leftarrow0\), however, \(Q=x\)
still holds through \(H\).

The proposed \(C\) is not the T4 crossing object.  The example refutes
"positive contrast plus a plausible boundary implies CTS."  It does not refute
OPCC, because the enlarged boundary \((C,H)\) may be a valid complete cover.

## 3. Bounded laboratory time does not imply a well-posed propagator

A deterministic evolution equation can have nonunique continuation from the
same Cauchy datum, blow up before the query, or require an omitted boundary
condition.  The writer and query may occur in a finite clock interval even though
there is no exact single-valued physical flow or normalized stochastic selection
law across the interval.

Such an episode fails T2 or T3.  A numerical solver that selects one branch does
not repair the physical premise.  If the actual theory supplies an additional
normalized selection kernel and all its inputs, the augmented episode may reenter
the theorem.

## 4. Non-Markov visible dynamics does not by itself defeat CTS

A colored-noise process or active memory can fail to be Markovian in a displayed
target coordinate.  That is not a counterexample when the physical bath/history
state is a standard-Borel Cauchy datum and is retained at \(G\).  It becomes an
excluded case only when the required memory or external input is omitted or no
exact causal transition law exists.

This prevents a false inference in both directions: the visible target alone is
not complete, but infinite-dimensional memory does not by itself require a new
admission class.

## 5. Type-III closures, disjoint states, and gauge centers are not automatic
counterexamples

A represented local AQFT von Neumann algebra may be type III rather than a
type-I tensor factor, and gauge constraints can produce centers, flux sectors,
edge data, or topological observables.  The W* and
C* branches do not require a tensor factor or finite-dimensional Hilbert space.
Accordingly, a W* algebra in one physical folium may be type III.  Separately,
the A branch can carry disjoint algebraic states only on one prospectively common
physical C*-observable algebra.  Those facts, nonfactorization, and centers are
not by themselves failures of CTS.

A **naive local tensor cut** that omits Gauss-law flux, a central/topological
degree, a Wilson-line dressing, or its correlations does fail T4.  Universal and
reduced Maxwell constructions can also differ in locality properties.  The lane
therefore requires the actual gauge-invariant time-slice algebra and complete
constraint data rather than asserting that every local gauge net qualifies.

Even a perfectly valid field time-slice algebra need not include the apparatus.
For example, a trivial field net with \(\mathcal A_{\rm fld}(N)=\mathbb C\)
satisfies its field time-slice condition, while a memory bit \(H=x\) crossing the
same cut and later read as \(Q=H\) requires at least a joint commutative separator
such as \(\mathbb C^2\).  Equations for the registered read cannot factor through
the field algebra \(\mathbb C\) alone.  This is why T4 requires
\(\widehat{\mathcal A}_G\) and a compatible joint state in addition to (4).

## 6. Formal perturbative algebra is not automatically a physical C* branch

A formal power-series star algebra with a formal time-slice property is not, by
that fact alone, a unital physical C*-algebra with positive normalized states and
UCP instruments.  Calling it one would erase T3's physical positivity premise.

Scoped perturbative AQFT results therefore motivate a possible construction but
do not establish `CTS` for an actual quantum-gravitational record.

## 7. Dynamical gravity is the open obstruction, not a proved counterexample

In gravity, a diffeomorphism-invariant observable can require dressing extending
far outside a nominal local laboratory.  The metric and causal relation can also
differ between event branches.  Consequently a proposed bounded \(D_r\), common
\(N_G\), or local separator algebra may fail T1 or T4; exact positive physical
instruments may also be unavailable, failing T3.

Merely naming an asymptotic port does not make a noncompact dressing fit T1's
bounded ancestor mission.  The actual physical support must meet T1.

No independently established bona-fide gravitational record with a proof that
**every** K/W*/C*/I cover fails is supplied here.  Dynamical gravity is therefore
outside the theorem and remains the sharp candidate obstruction, not an ontic
counterexample to `OPCC_FM`.

## 8. Internal indefinite order requires an external macro-cut

A finite quantum switch wholly before or wholly after \(G\) can be contracted as
an I-branch module.  By contrast, a proposed source/query pair embedded in a
process with no definite event-before-query macro-partition fails T1.  Abstract
process-matrix normalization does not manufacture the missing external order or
a physical spacetime realization.

Again this is an excluded case.  It becomes an OPCC counterexample only if a
bona-fide record is independently established and every possible complete
external macro-cut and standard cover is physically ruled out.

## 9. Actualization is upstream of this theorem

For an endogenous quantum event, CTS begins after an event alternative is
specified as part of the finite mission and addresses the record it leaves.
Branch-conditioned instruments and their recorded histories do not determine why
exactly one objective outcome occurs.  The lane neither solves nor assumes away
the strong outcome-as-record actualization problem.
