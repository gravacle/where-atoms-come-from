# Finite custody-handoff pair-memory mission

**Lane ID:** `CROSS-RFT-MGFT-FPMH-V001`

**Claim class:** exact finite model witness joining an actual custody handoff,
an independently qualified derivative record, the existing FHBC-to-DCL theorem,
a physical binary pair-memory gate, and a reversible matched support BREAK

**Disposition:**
`EXACT_CONTENT_BLIND_HANDOFF_RELATION_FORMED__KEEP_BREAK_ROUTES_MATCHED_RELATION_EXCITATION_TO_ACTIVE_OR_QUARANTINE__EXACT_LINK_SUPPORT_CONTRAST__REC_AND_FHBC_DCL_CONDITIONAL_MODEL_WITNESS__PHYSICAL_INSTANTIATION_ROBUST_LOCALITY_AND_GRAVITY_OPEN`

**Not claimed:** a new record axiom, universal FHBC/FCGHTS membership, SEALED
irreversibility, outcome selection, a laboratory implementation, robust many-pair
support, dimension above one, a metric, a tensor carrier, gravity, or `G`

## 1. Question and inherited theorems

`GRA-CH` showed conditionally that a physical pair memory `ell_uv` can gate
transition support. Its open physical premise was the formation, retention, and
matched physical BREAK of that memory. This packet asks only for the smallest
finite witness of that missing object.

No new coverage machinery is introduced. The existing finite Hamiltonian
boundary-closure theorem supplies

\[
\operatorname{FHBC}(r)
\Longrightarrow\operatorname{FCLPD}_W(r)
\Longrightarrow DCL_{\rm phys}(r),                    \tag{1}
\]

and

\[
\operatorname{REC}(r)\land\operatorname{FHBC}(r)
\Longrightarrow\operatorname{COV}_{\cup}(r).          \tag{2}
\]

The `FCGHTS` branch is not needed for this finite witness; it remains the other
branch of the already proved sufficient class
`OSP=FHBC or FCGHTS`.

## 2. Finite physical registers

Let carriers `A` and `B` be qutrits with basis

\[
\{|B\rangle,|0\rangle,|1\rangle\}.
\]

`|B>` is a unique physical blank and `|0>,|1>` span the carried content. Define
the occupation projector

\[
q=|0\rangle\!\langle0|+|1\rangle\!\langle1|.
\]

Let `L`, `K`, and `G` be binary relation registers:

- `L` is the fresh handoff writer target;
- `K` is the active retained pair-memory register coupled to the link gate; and
- `G` is an isomorphic quarantined garbage/reservoir register with no later
  coupling to the active query or link.

Let `a` be one link qubit, initially `|0>_a`. The prospectively bounded FHBC
device also contains the source transducer, route selector `z`, pulse controller,
clock, work/reference registers, and terminal outcome history. They are finite
and are retained in the device census even when their ideal state is unchanged.

The pair relation is undirected for support purposes:

\[
\ell_{AB}=\ell_{BA}:=K.
\]

Its writer is oriented by the physical custody event `A -> B`; the stored support
relation does not encode the transported content or its value.

### 2.1 Complete ideal-model port and outcome census

The finite witness has no implicit success filter. Its prospectively frozen
ports and registered roles are:

| Object | Input/role | Fate in the mission |
|---|---|---|
| content source `C` | common arbitrary `rho_C` | authorized source map places it at `A` (`F`) or `B` (`S`); no later source port |
| carrier pair `A,B` | one occupied carrier and one blank | common `U_H`; retained inside `D` |
| relation registers `L,K,G` | all blank | writer target, active route, and explicit quarantine; none is traced out of the device state |
| route context `z` | prospectively fixed `KEEP/BREAK` | selects one of the two permutation-related route pulses; retained in controller history |
| link `a` | blank | common gated pulse; terminal active-query coordinate |
| clock/controller/work/reference | fixed arm-common ready state | exact deterministic pulse history; retained inside `D` |
| registered active query | no event-label input | complete four-outcome instrument `(K,a) in {0,1}^2`; zero-probability outcomes remain in the alphabet |
| quarantine `G` | no query or link output port after routing | identity future and identity query effect; available only in separately declared calibration missions |
| external complement `E` | arbitrary normal root correlation allowed | exactly factorized dynamics and no registered operation during the mission |

Ideal deterministic unitaries have one complete controller outcome. If a
physical implementation admits loss, no-click, pulse failure, or saturation,
those outcomes must be added to the finite history and terminal instruments;
they cannot be discarded while retaining the theorem claim.

## 3. Exact content-blind handoff writer

For each `c in {0,1}`, let `U_H` exchange the two orthogonal basis states

\[
|c,B,0\rangle_{ABL}
\longleftrightarrow
|B,c,1\rangle_{ABL},                                  \tag{3}
\]

and fix every basis state outside these two two-dimensional exchange blocks.
It is therefore a Hermitian unitary, `U_H^2=I`.

Let `C` fix `|B>` and exchange `|0>` with `|1>`. Since (3) is identical for
both contents,

\[
[U_H,C_A\otimes C_B\otimes I_L]=0.                    \tag{4}
\]

Thus the writer is exactly content blind.

Take any density operator `rho_C` on the two-dimensional content subspace. The
authorized source alternatives prepare

\[
\begin{aligned}
\rho_F^{\rm in}
 &=\rho_C^A\otimes|B\rangle\!\langle B|_B
   \otimes|0\rangle\!\langle0|_L,\\
\rho_S^{\rm in}
 &=|B\rangle\!\langle B|_A\otimes\rho_C^B
   \otimes|0\rangle\!\langle0|_L .                  \tag{5}
\end{aligned}
\]

`F` is the formation arm: custody actually begins at `A` and is handed to `B`.
`S` is the matched prepositioned sham: the same content and total carrier
occupation are already at `B`. Applying the one common writer gives

\[
\begin{aligned}
U_H\rho_F^{\rm in}U_H^\dagger
 &=|B\rangle\!\langle B|_A\otimes\rho_C^B
   \otimes|1\rangle\!\langle1|_L,\\
U_H\rho_S^{\rm in}U_H^\dagger
 &=|B\rangle\!\langle B|_A\otimes\rho_C^B
   \otimes|0\rangle\!\langle0|_L.                  \tag{6}
\end{aligned}
\]

The post-write carrier state is exactly the same in both arms. In `F`, physical
occupation changes `(q_A,q_B):(1,0)->(0,1)` while fresh `L` changes `0->1`.
In `S`, `(0,1)` and `L=0` remain unchanged. Equation (6) holds for coherent
`rho_C`: the content is transferred, not cloned, and no value is learned by `L`.

After (6), the authorized source operation is over and all source couplings are
off. A positive identity hold of duration `tau_H>0`, or any declared common
Hamiltonian commuting with `L`, retains the relation before routing.

## 4. Matched reversible KEEP and BREAK

Both `K` and `G` begin blank. Freeze two context values, `z=KEEP` and
`z=BREAK`, before the source arm. Their route unitaries are

\[
U_{\rm KEEP}=\operatorname{SWAP}_{L,K},
\qquad
U_{\rm BREAK}=\operatorname{SWAP}_{L,G}.              \tag{7}
\]

They are related by the physical permutation `K <-> G`, have the same duration
and spectrum, and move one relation excitation into an initially blank,
degenerate destination. No information is deleted and both operations are
exactly reversible.

For a formed relation,

\[
\begin{array}{c|ccc}
 &L&K&G\\ \hline
\mathrm{KEEP}&0&1&0\\
\mathrm{BREAK}&0&0&1
\end{array}                                           \tag{8}
\]

whereas the sham remains `(0,0,0)` in either context. The total relation
occupation `L+K+G` is one in both formed intervention arms. In `BREAK`, `G` is
explicitly retained in the FHBC device but quarantined: every later active
Hamiltonian and registered active query acts as the identity on `G`.

This is a matched support intervention, not thermodynamic erasure. `KEEP`
transfers the relation into the active lineage; `BREAK` transfers the same
physical relation into an inaccessible-to-support reservoir and thereby blanks
the active coordinate.

## 5. Exact pair-memory-to-link support change

Use the existing pair-gated link form on this one pair,

\[
H_{\rm gate}=-h\,|1\rangle\!\langle1|_K\otimes X_a,
\qquad h>0.                                           \tag{9}
\]

It commutes with `K`, so each retained `K` eigenpath is exact. The state-sector
transition matrix element is

\[
\langle1_a|H_{\rm gate}|0_a\rangle
=\begin{cases}
-h,&K=1,\\
0,&K=0.
\end{cases}                                           \tag{10}
\]

For the frozen duration

\[
\tau_g={\pi\hbar\over2h},                             \tag{11}
\]

the active link flips, up to a global phase, exactly when `K=1`. Starting from
`a=0`, the complete active query `(K,a)` therefore has the deterministic table

\[
\begin{array}{c|cc|cc}
&\multicolumn{2}{c|}{\mathrm{KEEP}}
&\multicolumn{2}{c}{\mathrm{BREAK}}\\
&K&a&K&a\\ \hline
F&1&1&0&0\\
S&0&0&0&0
\end{array}.                                          \tag{12}
\]

Consequently,

\[
D_{\rm TV}[Q_{F,\rm KEEP},Q_{S,\rm KEEP}]=1,
\qquad
D_{\rm TV}[Q_{F,\rm BREAK},Q_{S,\rm BREAK}]=0,       \tag{13}
\]

and, conditional on formation,

\[
D_{\rm TV}[Q_{F,\rm KEEP},Q_{F,\rm BREAK}]=1.        \tag{14}
\]

Equations (10)--(14) are an exact lineage-localization and support-change test.
The Hamiltonian is common; the physical retained pair coordinate selects
whether the link transition exists on the actual eigenpath. The result is one
edge only. It proves no extended propagation or locality phase.

## 6. Existing `REC` predicate is satisfied in the finite witness

Define `r_AB` prospectively in the `KEEP` context with event alternatives
`F/S`, derivative-event provenance, fresh target `L`, active carrier lineage
`L -> K`, positive holds, and the complete label-blind terminal query of `(K,a)`.
Then the existing finite-mission record clauses are satisfied without changing
their definition:

1. **Prior absence:** `L=K=0` in both source arms before the common writer.
2. **Causal formation:** equation (6) writes `L=1` only on the actual
   `A -> B` custody handoff; the prepositioned sham has the same post-write
   carrier state but leaves `L=0`.
3. **Physical lineage and localization:** `KEEP` transfers `L` to active `K`;
   `BREAK` transfers it instead to quarantined `G`. Equations (13)--(14) make
   the active read and link response follow the physical relation destination,
   not a label.
4. **Source-off closure and positive hold:** after `U_H` no source port remains
   coupled, and a strictly positive common hold preserves the distinction.
5. **Relevance and distinguishability:** the complete query has unit total-
   variation formation/sham contrast in `KEEP`.
6. **Noncreating read:** the terminal projective instrument is common,
   label-blind, complete over all four `(K,a)` outcomes, and acts as the identity
   on `G` and on the event source. If the active lineage is blank, it cannot
   recreate the contrast.
7. **No replay or hidden bypass:** all source, writer, route, link, controller,
   clock, reservoir, and query degrees are in the finite device; `U_H` is used
   once and there is no post-write source or label input.

Therefore, for this exact model episode,

\[
\boxed{\operatorname{REC}(r_{AB}).}                   \tag{15}
\]

This is a model-conditional derivative record theorem. It is not an empirical
BFRE certificate or a claim that every physical pair memory realizes these
ideal premises.

## 7. FHBC constructs `DCL_phys`

The same packet satisfies the already defined FHBC premises:

- **H0:** `A,B,L,K,G,a`, source, controllers, clock, work/reference, and
  outcome-history registers form one prospectively fixed finite tensor-product
  device.
- **H1:** take the external complement dynamically decoupled during the mission.
  The displayed finite unitaries have bounded Hamiltonian generators; the
  piecewise source, write, hold, route, gate, and read schedule is exact and
  D-local.
- **H2:** one common joint root contains `rho_C` and all blanks. Only the
  authorized source preparation receives `F/S`. For each fixed context `z`,
  every later unitary and query is arm-common; `z` is an explicit prospectively
  chosen route context, not an undeclared second event label.
- **H3:** the terminal `(K,a)` pointer instrument contains all four outcomes
  and has a CPTP total map. Controller failures or additional physical outcomes,
  if admitted, are retained in the finite history register rather than
  postselected.
- **H4:** the mission has finite ordered source, write, positive hold, route,
  link evolution, and terminal query modules, with a nonempty state-wire path
  from the handoff writer to the query.

The imported theorem, rather than a new DCL construction, now gives

\[
\boxed{
\operatorname{FHBC}(r_{AB})
\Longrightarrow DCL_{\rm phys}(r_{AB}).}              \tag{16}
\]

Combining (15), (16), and the existing per-record closure gives

\[
\boxed{
\operatorname{REC}(r_{AB})\land\operatorname{FHBC}(r_{AB})
\Longrightarrow\operatorname{COV}_{\cup}(r_{AB}).}   \tag{17}
\]

No U-DCL postulate is needed for this finite witness.

## 8. Exact advance and remaining boundary

This packet closes one bounded cross-lane existence question:

\[
\boxed{
\text{actual content-blind custody handoff}
\longrightarrow\text{retained physical pair relation}
\longrightarrow\text{exact active link-support change}.}
\]

It also demonstrates why memory and ALLOW can act together without treating
`ALLOW` as energy: the relation bit does not add a force; it selects whether an
already declared transition matrix element is physically active. The reversible
BREAK moves that relation to explicit quarantine while conserving the relation
excitation.

The unresolved work is downstream and physical: instantiate the writer and
quarantine with complete ports; test robustness against missed/false pair
relations; form an extended composable pair-memory network without inserting a
target graph; show carrier custody through its collective phase; derive a stable
higher-dimensional local cone; and only then test tensor character, universal
stress coupling, reciprocal feedback, nonlinear closure, and gravity identity.
