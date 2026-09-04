# Independent hostile audit — GL6CK two-overlap `H6` tensor-writer response

**Target:** `LANE_CROSS_RFT_GRA_GL6CK_TWO_OVERLAP_H6_TENSOR_WRITER_RESPONSE_V001`  
**Verdict:** **PASS** on the stated finite, frozen-exterior,
spectral-response surface.

## 1. Independent audit surface

No target code was imported or executed by the independent science replay.
The audit rebuilt:

1. the source-free ring coefficient and Hermitian first-source gradient by
   summing literal six-flip histories;
2. the two actual Q4 cycles, a fresh 256-link degree-two completion, and the
   complete frozen-exterior locked-state census of their 29-link collar;
3. the isolated `K2` and three-state-star spectra and source matrix elements
   in exact `Q(sqrt(2))` arithmetic; and
4. every physical `h/U_d` power, coefficient, sign, and literal-versus-unit
   source normalization.

The executable is `verify_gl6ck_independent.py`; its frozen output is
`INDEPENDENT_RESULT.json`.

## 2. Global `H6` writer normalization — PASS

For each alternating orientation, every one of the `6!=720` temporal orders
was summed.  If `S_r` is the first `r` flipped cycle edges, the source-free
weight is

\[
 -\prod_{r=1}^{5}{1\over E(S_r)},
\]

and the exact total is

\[
 -{63\over8}.
\]

For a pair source at a cycle vertex, the audit independently differentiated
the five intermediate denominators relative to the Hermitian midpoint of the
two locked endpoints:

\[
 g_v=\sum_{\pi}\left(\prod_{r=1}^{5}{1\over E(S_r)}\right)
 \sum_{r=1}^{5}{M_v(S_r)-\tfrac12(M_v^0+M_v^1)\over E(S_r)}.
\]

The replay exhausts both alternating phases, six cycle vertices, all six
local port pairs, both cycle-port assignments, and both complement-port
assignments: 288 exact contexts, broader than the target's fixed-pair
48-context calibration.  Every context gives

\[
 g_v={105\over8}e_{ab},\qquad
 P_Tg_v={105\over16}\Theta_{v,c},\qquad
 \|\Theta_{v,c}\|^2=2.
\]

This verifies the sign and the `H+j.M` source convention.  A literal source
`j=s Theta` changes the aligned ring amplitude by
`(105/8)(h^6/U_d^6)s`.

## 3. Concrete two-overlap Q4 geometry — PASS

The cycles reconstructed directly from Q4 incidence are

```text
c0 = ((0,0,0),0) ((1,3,0),1) ((1,3,0),2)
     ((0,3,1),0) ((0,3,1),1) ((0,0,0),2)

c1 = ((0,1,3),0) ((1,0,3),1) ((1,0,3),2)
     ((0,0,0),0) ((0,0,0),1) ((0,1,3),2).
```

They share exactly the physical link `((0,0,0),0)` and hence exactly its
two endpoints `P(0,0,0)` and `C(1,0,0)`.  At the parent endpoint their port
pairs are respectively `{0,2}` and `{0,1}`; at the child endpoint those
pairs are interchanged.  Therefore at either endpoint

\[
 \Theta_0^2=\Theta_1^2=2,\qquad \Theta_0\cdot\Theta_1=0.
\]

This was not accepted as an abstract star.  The audit independently solved a
full Q4 degree-two completion containing both alternating cycles, froze its
exterior, and exhaustively solved the 29 active-link constraints on 28
affected nodes.  Exactly three states exist:

\[
 B,\qquad A=B\mathbin\triangle c_0,\qquad
 C=B\mathbin\triangle c_1.
\]

The double toggle is not locked, and the leaves differ on ten links.  Thus
the order-six configuration graph is the claimed physical three-state star,
not an arbitrary algebraic choice.

For the actual local source `dj/ds=Theta_0`, direct contraction gives

\[
 w_0={105\over16}{h^6\over U_d^6}\Theta_0^2
     ={105\over8}{h^6\over U_d^6},\qquad
 w_1={105\over16}{h^6\over U_d^6}
       \Theta_0\cdot\Theta_1=0.
\]

The concrete geometry therefore really realizes `w0 != w1` at either shared
endpoint.

## 4. Isolated-ring and star response — PASS

For one ring,

\[
 H_{K2}(s)=(-J+sw)\sigma_x.
\]

The source vertex and Hamiltonian commute.  The off-diagonal eigenbasis
matrix element is zero, the local ground branch is linear, and

\[
 E_{g,K2}''(0)=0.
\]

For the ordered star `(B,A,C)`, exact diagonalization gives the normalized
states

\[
 |g\rangle={\sqrt2|B\rangle+|A\rangle+|C\rangle\over2},\quad
 |m\rangle={|A\rangle-|C\rangle\over\sqrt2},\quad
 |u\rangle={\sqrt2|B\rangle-|A\rangle-|C\rangle\over2}
\]

with energies `-sqrt(2)J,0,+sqrt(2)J`.  For the source matrix with arm
weights `w0,w1`, the independently verified matrix elements are

\[
 \langle g|B|g\rangle=\sqrt2\,\bar w,\qquad
 \langle m|B|g\rangle=\delta,\qquad
 \langle u|B|g\rangle=0,
\]

where `wbar=(w0+w1)/2` and `delta=(w0-w1)/2`.  Both the spectral sum and
twice differentiating the exact branch give

\[
 \boxed{E_g''(0)=-{\sqrt2\over4J}(w_0-w_1)^2}.
\]

The sign is negative semidefinite, the factor of two in the spectral Hessian
is correct, and the response vanishes exactly when the two arm weights are
equal.

## 5. Physical scale and source normalization — PASS

Using

\[
 J={63\over8}{h^6\over U_d^5},\qquad
 w_0={105\over8}{h^6\over U_d^6},\qquad w_1=0,
\]

gives

\[
 E_g''(0)=-{175\sqrt2\over32}{h^6\over U_d^7}
\]

for the literal direction `dj/ds=Theta_0`.  Because `Theta_0^2=2`, the unit
direction `Theta_0/sqrt(2)` gives exactly half:

\[
 E_{g,unit}''(0)=-{175\sqrt2\over64}{h^6\over U_d^7}.
\]

The full local rank-one spectral tensor is equivalently

\[
 -{175\sqrt2\over128}{h^6\over U_d^7}
 (\Theta_0-\Theta_1)(\Theta_0-\Theta_1)^T.
\]

Here `s` has energy units, so `d^2E/ds^2` has inverse-energy units.  The
power count is exact: `w^2/J` scales as `h^6/U_d^7`.

## 6. Hostile boundary attacks

### Order-`h4` contamination

**Pass.**  The GL6BY finite-collar `h4` tensor is explicitly excluded.  The
audit independently counted zero four-cycles in the concrete Q4 incidence,
so no order-four locked-to-locked square process exists.  Every two-of-four
locked word has equal complementary pair products, making the bare, `h2`,
and diagonal global `h4` pair vertices pointwise `T2`-dark.  The order-six
writer is therefore the first relevant off-diagonal tensor source on this
declared parent.

### Missing contact

**Pass as a ceiling.**  The negative quantity above is the writer-induced
spectral contribution, not the full physical Hessian.  A source-second-
derivative contact and other pole terms may add to or cancel it.  The target
states this explicitly and makes no full-response or zero-crossing claim.

### Stationarity

**Pass as typed.**  The ground state is stationary in the declared
frozen-exterior three-state Hamiltonian.  The target does not claim that this
star is an invariant sector of the unfrozen full Q4 parent or a selected
thermodynamic state.

### Accumulation language

**Pass as typed.**  “Accumulated” means the exact composition of two
overlapping ring writers.  It does not mean a linked-cluster sum,
thermodynamic limit, common cone, or macroscopic phase.

### Physical promotion

**Pass.**  The packet does not authenticate either ring as a record and does
not claim a continuum solder, metric, Ricci/Einstein response, gravity, or a
value of `G`.

The upstream GL6CH author packet has no separately sealed hostile audit at
this checkpoint.  That custody fact is already disclosed by GL6CK.  It does
not weaken this audit's finite result because the exact writer coefficients
needed here were independently rederived across 288 local contexts.

## 7. Verdict

**PASS.**  GL6CK establishes exactly what it claims: the global order-six
candidate-field writer is spectrally inert on one isolated ring but produces
a nonzero, exactly normalized stationary spectral curvature on the smallest
concrete shared-edge locked star whenever its two physical ring channels are
weighted differently.  It is a real finite accumulation mechanism and a
valid building block for the next bulk/contact calculation; it is not yet a
gravity theorem.
