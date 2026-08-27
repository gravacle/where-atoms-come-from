# Record-front causal-diamond theorem

**Lane ID:** `CROSS-RFT-GRA-EI-RFCD-V001`

**Official short name:** `RFCD`

**Date:** 2026-08-27

**Claim class:** exact conditional finite-depth causal-order theorem; exact
dimensionless interval reconstruction; exact higher-rank product-cone
obstruction

**Not claimed:** that current F3 autonomously selects or protects the required
front quotient; that the order fiber disappears; a physical clock or cell
scale; a continuum limit; our `3+1` dimensional spacetime; curvature, a
spin-two mode, universal stress coupling, gravity, or Newton's constant

## 1. Physical antecedent

Let `A` and `B` be two reusable, content-symmetric physical operations on one
event-bearing front.  At finite depth `R`, retain every complete history and
all of its ports.  Assume the following already-isolated physical clauses.

1. **Operational front merger.**  Histories that differ only by adjacent
   interchange of `A` and `B` have one common future event-bearing front.  The
   complete order distinction is retained in explicitly named complementary
   history factors, so no complete information is erased; those factors are
   called records only if independently `REC`-qualified.  In the notation of
   the reversible operational-front quotient theorem, every admitted future
   front channel descends through `pi=Tr_Z`.
2. **No extra front identifications or shortcuts.**  Front classes are
   identified exactly by the adjacent interchange relation, and their only
   elementary future transitions append one authenticated `A` or one
   authenticated `B`.
3. **Physical direction.**  Appending an operation is a later event-bearing
   transition.  Its inverse is not silently admitted as a formation step.
4. **Complete branch custody.**  Success, failure, quarantine, boundary, and
   terminal outcomes remain explicit.  The theorem applies on every
   nonterminal branch on which clauses 1--3 hold.

Fix, in addition, a prospectively declared active-success mission `D_R`.  For
every `|w|<R`, `D_R` contains one admitted nonterminal active-front outcome and
both authenticated append outcomes `wA,wB`.  For `|w|=R`, `D_R` contains one
declared cap outcome, which may be terminal.  Early terminal, failure, and
quarantine outcomes lie outside `Q_R`; the declared depth-`R` cap states lie
inside `Q_R` but have no outgoing append in the mission.  Thus `Q_R` below is
prospective mission support, not the sequence of states visited along one
realized history.

These clauses are stronger than the existence of a comparison record or an
undirected edge between `AB` and `BA`.  They are exactly the physical role
isolated by CM and made reversibly possible by CN.  URFT governs whether the
retained histories and order/history fiber qualify as records; it neither qualifies
them without the `REC` premises nor supplies the front-merger law.

## 2. The front quotient

For a word `w` define

\[
 u(w)=\#A(w),\qquad v(w)=\#B(w).                 \tag{EI01}
\]

Adjacent interchanges connect all and only words with the same count pair, so
the event-bearing front classes are

\[
 Q_R=\{(u,v)\in\mathbb N_0^2:u+v\le R\}.         \tag{EI02}
\]

The directed elementary transitions are

\[
 (u,v)\longrightarrow(u+1,v),\qquad
 (u,v)\longrightarrow(u,v+1).                    \tag{EI03}
\]

Write `p <= q` when a directed path leads from `p` to `q`.  Then

\[
 (u,v)\preceq(u',v')
 \quad\Longleftrightarrow\quad
 u'\ge u\ \hbox{and}\ v'\ge v.                  \tag{EI04}
\]

This order is derived from future composition.  It is not obtained by first
drawing an undirected grid and then calling one coordinate time.

## 3. Theorem RFCD-1 -- exact causal-diamond order

Define the dimensionless coordinates

\[
 t=u+v,\qquad x=u-v.                              \tag{EI05}
\]

For any two front classes `p,q`, with differences `Delta t, Delta x`,

\[
 \boxed{
 p\preceq q
 \quad\Longleftrightarrow\quad
 \Delta t\ge |\Delta x|,
 \quad \Delta t\equiv\Delta x\pmod 2 .}
                                                               \tag{EI06}
\]

Under the representative quadratic form

\[
 ds^2=-dt^2+dx^2=-4\,du\,dv,                    \tag{EI07}
\]

every elementary `A` or `B` transition is future null.  The directed front
therefore has exactly the causal order of the `1+1`-dimensional Minkowski
diamond lattice in the truncated future wedge

\[
 Q_R\cong\{(t,x)\in\mathbb Z^2:
 0\le t\le R,\ |x|\le t,\ t\equiv x\pmod2\}.       \tag{EI07a}
\]

Its boundaries `x=+/-t` are null and `t=R` is the terminal mission cap.  It is
not the full translation-invariant Minkowski lattice; physical inverse/replay
steps have not been admitted.

### Proof

For `q-p=(a,b)`, equation (EI04) says `a,b>=0`.  Equation (EI05) gives

\[
 a={\Delta t+\Delta x\over2},\qquad
 b={\Delta t-\Delta x\over2}.                    \tag{EI08}
\]

Both are nonnegative integers exactly when (EI06) holds.  One `A` step has
`(Delta t,Delta x)=(1,1)` and one `B` step has `(1,-1)`, so both have zero
quadratic interval under (EI07) and point toward increasing `t`. QED.

The exact-depth front `t=n` contains `n+1` sites with
`x=-n,-n+2,...,n`.  The truncated causal future of the root through depth `R`
contains

\[
 |Q_R|={(R+1)(R+2)\over2}.                        \tag{EI09}
\]

Thus `|Q_R|` is the counting measure of the truncated operational causal-order
transition poset.  It becomes a physical event/spacetime-volume census only
under a separately proved one-class/one-event lift and an independently
calibrated uniform event-volume identification.  Neither follows from CM/CN.
The `t=n` antichain has one independent contrast label `x` and `n+1` elements;
there is no within-slice edge or earned spatial topology/metric.  The quadratic
count must not be described as two earned spatial dimensions.

## 4. Theorem RFCD-2 -- interval volume reconstructs the discrete interval

Let `p <= q`, put `a=u'-u`, `b=v'-v`, and let

\[
 I(p,q)=\{z:p\preceq z\preceq q\}.               \tag{EI10}
\]

Let `ell(p,q)` be the maximum number of covering edges in a chain from `p` to
`q`; here every maximal chain has that same edge count.  The chain-edge count
and interval census are

\[
 \ell(p,q)=a+b=\Delta t,
 \qquad |I(p,q)|=(a+1)(b+1).                      \tag{EI11}
\]

Consequently the chosen flat representative has the exact dimensionless
combinatorial interval

\[
 \boxed{
 \tau_{\rm comb}^2(p,q):=\Delta t^2-\Delta x^2
 =4ab=4\bigl(|I(p,q)|-\ell(p,q)-1\bigr).}        \tag{EI12}
\]

### Proof

The interval is the integer rectangle
`[u,u'] x [v,v']`, proving its census.  Every directed path from `p` to `q`
contains exactly `a` `A` steps and `b` `B` steps, proving the chain-edge
count.  Expanding `(a+1)(b+1)` and using
`Delta t=a+b`, `Delta x=a-b` proves (EI12). QED.

Equation (EI12) defines an exact combinatorial interval `tau_comb^2` for the
chosen flat representative (EI07).  Calling it physical proper time requires a
one-class/one-event lift, independent identification of a uniform physical
event-volume density, clock/length calibration, and common-probe use of this
order.  Absent those, the physical conformal factor--not only its global
normalization--remains open.

The directed formation order is also distinct from transport generated by an
undirected adjacency Hamiltonian.  Forgetting the arrows adds edges toward
smaller `t`, and `H=epsilon I-t_hop A` evolves in a separate external parameter;
neither fact makes its response retarded in the RFCD order.  A physical
common-probe cone therefore needs a separate propagation/descent theorem.

## 5. Constitutive intervention

Under front-merger KEEP, the exact-depth census is `n+1` and RFCD-1 holds.
Call the matched intervention that removes every adjacent-interchange
identification, keeps both operations admitted, and supplies no replacement
relation **FULL-INTERCHANGE BREAK**.  Then all binary words remain distinct:

\[
 F_n^{\rm BREAK}=2^n.                             \tag{EI13}
\]

No injective identification of these depth-`n` fronts with the RFCD diamond
slice exists for `n>=2`, because that slice has only `n+1` sites.  A query that
reads the retained order fiber while leaving all future front channels
fiber-blind does not change RFCD.  A generic **DESCENT BREAK** is any admitted
fiber-to-front coupling that invalidates RFCD autonomy.  It may produce a
partial or refined quotient and need not yield (EI13) or the free-word tree.

This is the precise memory statement: the complete record of order can remain
real, but redundant order information must not multiply independent future
influence fronts if the RFCD operational-poset organization is to survive.

## 6. Theorem RFCD-3 -- why more commuting types do not give `3+1`

Under the analogous total mission domain and no-extra-identification clauses,
`q` completely commuting reusable operation types have count front `N_0^q`
with future tangent cone

\[
 C_q=\mathbb R_+^q.                               \tag{EI14}
\]

For `q=2`, (EI05) maps this cone exactly to the `1+1` Lorentz cone.  For every
`q>=3`, no invertible linear transformation maps `C_q` onto the standard
`1+(q-1)` Lorentz cone.

### Proof

`C_q` is polyhedral and has exactly `q` extreme rays.  In dimension at least
three, the standard Lorentz cone has a continuum of null extreme rays, one for
every direction on `S^(q-2)`.  Invertible linear maps preserve extreme rays
and polyhedrality.  The cones therefore cannot be linearly equivalent. QED.

The naive rule "one commuting record type per spacetime coordinate" is thus
rejected.  A route to `3+1` must earn collective isotropization, a sufficiently
rich family of causal directions, and one common smooth quadratic cone.  A
network of authenticated causal-diamond cells remains a candidate; its gluing,
dimension, transport, and scaling law must be derived rather than inserted.
The theorem excludes only direct invertible-linear identification of the raw
orthant; it does not exclude a separately derived nonlinear or coarse-grained
isotropizing limit.

## 7. Relation to history-wise formation and gravity

If `Phi_*` is prospectively defined to include EI clauses 1--4 and the total
`D_R` domain, and `W_*` satisfies EE04/EE04a on every branch, EE conditionally
preserves that typed RFCD phase history-wise.  EF applies to a separately
authenticated endpoint projection only when UCAIC, ACL, and SAF all hold;
owned cells and bounded capacity alone are insufficient.  The stronger
orthant lift/frontier composition may be imported only after EG is
source-frozen and pinned.  Present FPMH does not derive the merger/descent
phase or a universal cell gluing law.

RFCD closes one minimal origin implication:

\[
 \boxed{
 \begin{gathered}
 \text{independently REC-qualified retained physical histories}
 +\text{future-front coalescence}
 \Longrightarrow\\
 \text{an exact discrete operational causal-order transition poset}\\
 \text{with a dimensionless interval representative in }1+1
 \end{gathered}}
                                                               \tag{EI15}
\]

It is a proof of principle for records becoming constitutive of causal
organization, not a gravity proof.  Gravity still requires an actual
higher-dimensional common Lorentzian phase, physical scale/volume, a common
probe metric, a protected tensor response, universal complete-stress coupling,
and back-reaction from the same parent.

**Disposition:**

`CONDITIONAL_OPERATIONAL_RECORD_FRONT_MERGER_GIVES_EXACT_TRUNCATED_1PLUS1_CAUSAL_DIAMOND_POSET__ORDER_PLUS_INTERVAL_COUNT_RECONSTRUCTS_DIMENSIONLESS_COMBINATORIAL_REPRESENTATIVE__FULL_INTERCHANGE_BREAK_RESTORES_EXPONENTIAL_FRONTS_WHILE_GENERIC_DESCENT_BREAK_NEED_NOT__ORDER_FIBER_MAY_REMAIN_RECORDED_IF_FUTURE_BLIND__EVENT_LIFT_PHYSICAL_VOLUME_COMMON_PROBES_AND_CONFORMAL_FACTOR_OPEN__FINITE_HIGHER_RANK_COMMUTATION_CONE_NOT_LINEARLY_LORENTZIAN__3PLUS1_METRIC_AND_GRAVITY_OPEN`
