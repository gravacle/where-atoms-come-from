# GL6CG — Global owner-once diagonal order-`h4` pair-source coefficient

## 1. Result

Insert the six local port-pair sources before the F3 Feshbach reduction and
sum the microscopic histories on the complete incidence parent, counting
each link and each node-owned adjacent-link pair once.  On every locked
degree-two state of any **simple 4-regular bipartite incidence graph** with
consistent four-port labels, the exact diagonal order-`h4` first-source
coefficient at node `v` is

\[
 V^{(4)}_{v,ab}=-{4\over9}-{37\over12}z_{va}z_{vb},
 \qquad 0\leq a<b\leq3.                                  \tag{CG.1}
\]

Here `z_va=1-2n_va` and exactly two of the four `n_va` equal one.  The
coefficient (CG.1) multiplies
`r^4=(h/U_d)^4` in `dH_eff/dj` for source convention `H+j.M`.

In pair order `(01,02,03,12,13,23)` use

```text
A  = (1, 1, 1, 1, 1, 1)
Ea = (1, 1,-2,-2, 1, 1)
Eb = (1,-1, 0, 0,-1, 1)
T1 = (1, 0, 0, 0, 0,-1)
T2 = (0, 1, 0, 0,-1, 0)
T3 = (0, 0, 1,-1, 0, 0).
```

The locked identity
`z0 z1 z2 z3=1` makes complementary products equal.  Therefore

\[
 P_{T_2}V^{(4)}_v=0                                      \tag{CG.2}
\]

**pointwise, for every node and every locked state**.  The `A1` contraction
is the state-independent value `7/2`; only `E2` is state dependent.  Thus the
global owner-once **diagonal** `h4` pair-source coefficient is an `A1`
identity plus an `E2` operator.  On the inherited GL6AO domain—degree four
and girth at least six—there is no four-link locked-to-locked process, so
this diagonal coefficient is also the complete order-`h4` source operator
and contains no `T2` operator at any momentum.  The completeness conclusion
does not extend to girth-four graphs: an alternating square can carry an
off-diagonal `T2` source even though (CG.1)--(CG.2) remain true.

This does not make `V4` the first or complete `E2` response.  The diagonal
first-source hierarchy through this order is

\[
 P M_vP=M_v,\qquad V^{(2)}_v=-M_v,\qquad
 V^{(4)}_v=-{4\over9}{\bf1}_6-{37\over12}M_v.            \tag{CG.2a}
\]

Consequently the bare locked-space `E2` source is already state dependent;
`V2` and `V4` renormalize it.  No quotient removing that bare vertex is
declared here.  The tensor distinction is sharper:

\[
 P M_{T_2}P=V^{(2)}_{T_2}=V^{(4)}_{T_2}=0.               \tag{CG.2b}
\]

This is the decisive result on the actual girth-six parent.  There, the
nonzero `h4` `T2` vertices of frozen collars are partial
boundary/incomplete-owner pieces; they cancel when the full parent histories
are owned once, and no square transition remains.  It does not follow that
the complete tensor response vanishes: the same-state order-`h2` contact
remains, and an off-diagonal sourced `h6` or higher-order ring operator can
open a tensor channel.  On a different girth-four graph, sourced square
transitions must instead be included explicitly.

## 2. Complete F3 history reduction

Let `s` be any locked state, let `M(s)` collect all local pair observables,
and define

\[
 q_i=M(s\mathbin\triangle i)-M(s),\qquad
 q_{ij}=M(s\mathbin\triangle\{i,j\})-M(s).              \tag{CG.3}
\]

Every one-link intermediate has dimensionless defect energy `2`.  Let
`p_ij` be the two-link defect energy.  Differentiating the complete
intermediate-normalized fourth-order diagonal Feshbach expression gives

\[
 V^{(4)}(s)=-{3N\over16}\sum_iq_i+
 \sum_{i<j}\left[{q_i+q_j\over2p_{ij}}+
                        {q_{ij}\over p_{ij}^2}\right],  \tag{CG.4}
\]

where `N` is the number of microscopic links.  For disjoint links,
`p_ij=4` and `q_ij=q_i+q_j`.  Since every link is adjacent to six other
links in a simple 4-regular incidence graph, all disjoint terms can be
eliminated exactly.  The translation-local form is

\[
 V^{(4)}(s)=-{3\over16}\sum_iq_i+
 \sum_{\substack{i<j\\i,j\ {\rm adjacent}}}
 \left[\left({1\over2p_{ij}}-{3\over16}\right)(q_i+q_j)
             +{q_{ij}\over p_{ij}^2}\right],           \tag{CG.5}
\]

with `p_ij` exactly `2` or `6`.  Equation (CG.5) is the machine-ready
operator stencil.  It is not a sum of overlapping collar responses.

For `Q_L`, (CG.5) counts `4L^3` links, `2L^3` nodes, and `12L^3`
node-owned adjacent-link pairs.  At `L=4` these counts are `256`, `128`, and
`768`.  The exact script also evaluates (CG.4) over all `32640` unordered
link pairs of the selected GL6CC state and obtains (CG.5) byte for byte.

## 3. Universal 486-neighborhood proof

The coefficient at one node depends only on its four incident links and the
other three links at each of their neighboring endpoints.  There are six
central degree-two words.  Once the shared occupation of an incident link
is fixed, its degree-two neighbor has exactly three compatible external
three-bit words.  Hence

\[
 6\,3^4=486                                                \tag{CG.6}
\]

exhausts the radius-one local possibilities.  This enumeration is a
superset when longer graph cycles correlate neighbor choices, so it proves
the result for every simple graph in the stated domain.

For each incident link, the three external neighbor links contain one link
with the same occupation and two with the opposite occupation.  The
neighbor-owned contribution is consequently independent of which compatible
word was chosen.  Exact reduction of the central and neighbor-owned terms in
(CG.5) gives (CG.1).  The enumeration provides a second proof and has the
histogram

```text
(-127, 95, 95, 95, 95,-127)/36   count 162
(  95,-127,95, 95,-127, 95)/36   count 162
(  95, 95,-127,-127,95, 95)/36   count 162.
```

Those are precisely the three complementary-pair patterns of a locked node.
Their projected values are

```text
A1: 7/2 in all 486 cases
E2: (-37/3,-37/3), (-37/3,37/3), or (74/3,0)
T2: (0,0,0) in all 486 cases.
```

## 4. Exact GL6CC background commutator row

The repaired, audited period-four GL6CC background is locked and has `64`
outgoing elementary alternating hexagons.  With

\[
 H_6=-J A_{\cal C},\qquad J={63\over8}{h^6\over U_d^5}, \tag{CG.7}
\]

the exact matrix element for a diagonal source component is

\[
 \langle s\mathbin\triangle c|[V_X^{(4)}(k),H_6]|s\rangle
 =-J\,[V_X^{(4)}(s\mathbin\triangle c;k)-V_X^{(4)}(s;k)]. \tag{CG.8}
\]

The packet computes the complete 64-edge row from the selected basis state,
using both the infinite-lift local stencil and direct Q4 reduction.  Every
transition profile vanishes nodewise in `A1` and `T2`.  At `k=0` the
six-component row has rank two and kernel

\[
 A_1\oplus T_2.                                           \tag{CG.9}
\]

The `E2` span has rank two.  Sixteen of the 64 moves are fully dark at
`k=0`; every one of those sixteen has a nonzero first spatial moment in
`E2`.  The other 48 are `E2`-bright at `k=0`.  The complete
orientation/background-phase classification is in `EXACT_LEDGER.json`.

For the summed transition norm, the parity-odd `k1` coefficient cancels
exactly.  The raw first-moment Gram has rank `6/18`: it samples all of
`E2 tensor R^3` and has the exact 12-dimensional `A1 tensor R^3` plus
`T2 tensor R^3` kernel.  There is no `E/T` mixing because the `T2`
amplitude is identically zero before any momentum expansion.

This 64-edge row is a diagnostic, not a stationary residue.  The selected
basis state has 64 outgoing moves and is not an eigenstate of (CG.7).

## 5. Native link `T2` is a different operator

For every one of the 64 cycles, the inherited tetrahedral/Z3 coordinates
reproduce the audited GL6AS native link symbol

\[
 C_a=z_b-z_c\quad\hbox{and cyclic partners}.             \tag{CG.10}
\]

The native link density is complement odd, is separately `k=0`-dark in
each port, and opens at first spatial order; its selected-background `T2`
gradient Gram has rank `3/9`.  In contrast, (CG.1) is complement even and
its pair-memory `T2` projection is the zero operator at every `k`.

Thus the common tetrahedral representation label `T2` does not identify
pair memory with the conserved native link density.  No link-continuity or
gauge conclusion transfers to the pair source.

## 6. Finite-component obstruction and stationary continuation

An exact H6 breadth-first expansion from the GL6CC state gives new-shell
counts

```text
depth 0..4: 1, 64, 1,952, 37,968, 532,080
total through depth 4: 572,065.
```

The induced depth-zero/one sector is exactly a 65-state star with 64 edges,
but it leaks to depth two and is not invariant.  A stationary diagonalization
or a branch gap extracted from this star would therefore be spurious.  This
packet forms neither one.

For any finite connected flip component `C`, audited GL6AR gives the unique
positive Perron--Frobenius ground state `|0>` of (CG.7), energy `E0`, and a
positive component gap.  Put `Q=1-|0><0|`.  The correctly typed static
same-state continuation is

\[
 W_{AB}(k)=\langle0|C^{(2)}_{AB}(k)|0\rangle
 -2\,\mathrm{Re}\,\langle0|\delta V_A(-k)
 Q(H_6-E_0)^{-1}Q\delta V_B(k)|0\rangle.                \tag{CG.11}
\]

Here `C^(2)` is the general-state owner-once source-before-Feshbach contact,
and `V=r^4 V^(4)` is the physical first-source vertex.  Both expectations
must use the same `|0>` or the same stationary density.  Factoring
`h^2/U_d^3` from (CG.11) gives the exact pole coefficient `16/63` with
dimensionless resolvent `(rho_C-A_C)^(-1)`.

Equation (CG.2) makes every spectral term with a `T2` leg vanish at this
order.  Therefore the leading `T2` block in (CG.11) is the same-state
order-`h2` contact alone.  The formal scale
`V4^2/J=(8/63)h^2/U_d^3` is correct, but its tensor coefficient is exactly
zero.  This does not evaluate the contact expectation and does not exclude
an off-diagonal sourced `H6` term or a higher-order tensor vertex.

The resolvent is equivalently the convergent imaginary-time integral of
`Q exp[-tau(H6-E0)]Q` on the finite gapped component.  Since `-J A_C` is
stoquastic and its PF vector is positive, diagonal `V4` and contact
estimators admit nonnegative worldline sampling weights.  Complex Fourier
profiles can be split into sine and cosine observables; signs of the measured
observable are not a sampling-weight sign problem.

No real-time kernel is formed.  A later retarded Kubo calculation must keep
its explicit `1/hbar` and excitation frequencies `(En-E0)/hbar`; the static
resolvent (CG.11) itself contains no `hbar`.

## 7. Claim boundary

This packet proves a broad owner-once diagonal `h4` coefficient theorem and,
on the inherited degree-four girth-at-least-six parent, a complete `h4`
no-`T2` operator theorem.  It also gives one complete nonstationary H6
commutator row plus a matrix-free stationary prescription.  It does **not**
compute a stationary contact expectation,
stationary pole, full CTP functional, connected/Moebius weight, or `Phi_C`.
It supplies no orbit average, material law, metric, record authentication,
Ricci tensor, gravity law, threshold, or value of `G`.  The inherited
tetrahedral/Z3 coordinates are a calculation chart, not a derived physical
spacetime.
