# Independent hostile audit: GL6BV V001

**Target:** `LANE_CROSS_RFT_GRA_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001`  
**Audit:** `AUDIT_G_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001`  
**Date:** 2026-09-02  
**Verdict:** `PASS_AFTER_AUTHOR_LITERAL_BYTE_REPAIR`  
**Target edits by auditor:** none

## 1. Verdict

The repaired `GL6BV` source-before-Feshbach theorem is correct within its
declared pure local-pair `T2` source chart and active-edge boundary.  The
order-`h^2` static functional, its denominator domain, all first and second
source derivatives, the microscopic gapped retarded normalization, the
translation/`S4` parametrization, the bounded `Q4` witness, and the fixed-
projection scope of BV35 all pass independent exact reconstruction.

The result is a nonzero **high-gap `T2` contact/response ingredient**.  It is
not a low-energy `T2` pole and cannot be combined with the selected
six-core isolated-hexagon `E2` response as though both came from one active
microscopic collar or one stationary state.  The complete eighteen-link
source-first `h^2/h^4/h^6` calculation remains open.  Consequently the
common/Ricci diagnostic has not been tested in one completed causal
coefficient, and no gravity or `G` conclusion follows.

## 2. Author repair and custody

The first frozen snapshot contained two literal theorem-byte defects:
`{cal P}` lacked the calligraphic backslash in BV01, and BV27 contained an
embedded carriage return in the intended `\rho` subscript.  The auditor
reported both without editing the target.  The author repaired them,
removed every `0x0d` byte, rebuilt the manifest and seal, and reproduced the
same `1394/1394` result.  This audit pins and replays only that corrected
snapshot.

No scientific formula changed in the repair.

## 3. Defect-frame and first-vertex reconstruction

In pair order `(01,02,03,12,13,23)`, the four degree-one/three defect vectors
are

```text
tau_0=(-1,-1,-1,+1,+1,+1)
tau_1=(-1,+1,+1,-1,-1,+1)
tau_2=(+1,-1,+1,-1,+1,-1)
tau_3=(+1,+1,-1,+1,-1,-1).
```

The audit rebuilt the `A1/E2/T2` projectors and found

\[
 P_T\tau_a=\tau_a,\qquad
 \tau_a^T\tau_b=8\delta_{ab}-2,\qquad
 \sum_a\tau_a=0,
 \qquad\sum_a\tau_a\tau_a^T=8P_T.                       \tag{A-BV01}
\]

For every one of the six locked local states, the four incident link flips
produce the four `tau_a` in a state-dependent permutation.  The independent
replay additionally exhausted every nonempty proper subset: no subset of
one, two, or three vectors has zero sum.  Thus the source-off first vertex
cancels precisely when all four incident equal-`h` flips are active; a
six-core hexagon, with only two active incident flips at a cycle node, does
not inherit that cancellation.

## 4. Static source-before-Feshbach functional

For

\[
 H(j)=U_d\sum_vq_v^2-h\sum_eX_e+\sum_vj_v^TM_v,
 \qquad P_Tj_v=j_v,                                      \tag{A-BV02}
\]

the locked pair `T2` value is zero.  One link flip creates one unit defect
at each endpoint, so its exact source-dressed cost at order `h^2` is

\[
 \Delta_e(j;\eta)=2U_d+j_v^Tt_{v,e}(\eta)
                         +j_w^Tt_{w,e}(\eta).             \tag{A-BV03}
\]

The inherited incidence is simple: two flips return to the lock only when
the same edge is flipped twice.  One physical edge therefore owns exactly
one denominator, giving

\[
 F^{(2)}(j)=-h^2\sum_{\eta,e=\{v,w\}}
 { |\eta\rangle\langle\eta|\over\Delta_e(j;\eta)}.       \tag{A-BV04}
\]

The exact domain is the connected component containing `j=0` on which every
displayed denominator is nonzero.  Since `||tau_a||^2=6`, two applications
of Cauchy--Schwarz prove the stated sufficient ball

\[
 \max_v\|j_v\|<{U_d\over\sqrt6}.                         \tag{A-BV05}
\]

Within that ball the two endpoint shifts have total magnitude strictly less
than `2U_d`, so no denominator can cross zero.

Direct differentiation fixes both sign and normalization:

\[
 {\partial F^{(2)}\over\partial j_v}\bigg|_0
 ={h^2\over4U_d^2}\sum_{e\ni v}t_{v,e},                  \tag{A-BV06}
\]

and, with all four incident flips active,

\[
 \boxed{F^{(2)}_{j_vj_v}(0)=-{2h^2\over U_d^3}P_T}.      \tag{A-BV07}
\]

For the endpoints of one edge,

\[
 \boxed{F^{(2)}_{j_vj_w}(0;\eta)
 =-{h^2\over4U_d^3}t_{v,e}(\eta)t_{w,e}(\eta)^T}.        \tag{A-BV08}
\]

Every allowed edge block is nonzero and rank one; the audit independently
found all nine locally allowed blocks per port.  The reverse block is its
transpose.  Assembling each physical edge once gives

\[
 F^{(2)\prime\prime}(0;\eta)
 =-{h^2\over4U_d^3}{\cal A}_\eta^*{\cal A}_\eta,          \tag{A-BV09}
\]

which is negative semidefinite.  The independent `Q4` sparse assembly has
128 onsite `8P_T` blocks and 256 owned rank-one edge blocks, and two unrelated
exact source fields reproduce the quadratic norm `||A j||^2`.

At `j=0` the locked branches are degenerate.  The operator/branch Hessian is
well defined, but the minimum of the branch energies can be nonsmooth.
`GL6BV` correctly does not promote BV19 to a unique multisite ground-energy
Hessian.

## 5. Retarded normalization

The dressed locked state has one-link amplitude `h/(2U_d)`.  Therefore the
transition Gram is

\[
 {h^2\over4U_d^2}{\cal A}^*{\cal A}
\]

at gap `2U_d`.  With the inherited
`(i/2) theta(t)<[M(t),M(0)]>` convention, the onsite sine coefficient is

\[
 \boxed{K^R_{vv,T}(t)
 ={2h^2\over U_d^2}\theta(t)\sin(2U_dt)P_T+\cdots .}     \tag{A-BV10}
\]

The zero-frequency sine integral contributes `1/(2U_d)`, so

\[
 \boxed{F^{(2)\prime\prime}(0)=-2K_T^R(0)}.              \tag{A-BV11}
\]

This factor and sign pass independently.  The eliminated contact and the
uneliminated high-gap Kubo history describe the same virtual defects and
must not be added as separate effects.

## 6. Translation/S4 parametrization and Q4 witness

For an edge of port `a`, each endpoint partner label lies among the other
three ports.  The port stabilizer has exactly two orbits on the ordered
endpoint-label pair: equal and unequal.  Thus a translation- and
`S4`-invariant locked density is described by one alignment probability
`p`.  Direct averaging gives

\[
 C_a=\beta P_T+\left({2\over3}-\beta\right)Q_a,
 \quad \beta=4p-{4\over3},\quad Q_a={\tau_a\tau_a^T\over6}. \tag{A-BV12}
\]

The audit reconstructed the defining period-four flow independently.  It is
locked at all 128 nodes, contains the selected alternating hexagon, and has

\[
 p={218\over256}={109\over128}.                           \tag{A-BV13}
\]

An explicit translation/port-orbit average agrees with (A-BV12), yielding

\[
 \beta={199\over96},\qquad
 \delta={2\over3}-\beta=-{45\over32},\qquad
 \gamma={77\over12},                                    \tag{A-BV14}
\]

and fixed common/relative zero-character values
`173P_T/12` and `19P_T/12`.

The parity-even quadratic coefficients in BV35 also replay exactly:

```text
fixed common:   -199/192 I2 P_T + 45/64 sum_a theta_a^2 Q_a
fixed relative: +199/192 I2 P_T - 45/64 sum_a theta_a^2 Q_a.
```

For `theta=(1,-1,0,0)`, the independent linear-character matrix is nonzero.
In the fixed common/relative basis it is off-diagonal.  BV35 is therefore a
fixed diagonal projection only, not an eigenbranch dispersion and not a
Schur-reduced `k^2` coefficient.  The target states this limitation
correctly.

## 7. Active-support boundary

The audit compared the two `Q4` locked configurations related by the target
hexagon toggle.  In the separate full `Q4` embedding every cycle node owns
four incident one-link histories and the full `8P_T` frame.  That algebraic
fact does not attach the coefficient to the isolated `K2` state.

In the six-core isolated control, each of the six cycle nodes owns exactly
two incident core flips.  Before and after the toggle those two defect
vectors form the same unordered pair, but their sum is nonzero and their
outer-product frame has rank two.  Thus the partial contact is common to the
two basis configurations but anisotropic; neither the complete first-vertex
cancellation nor the universal onsite coefficient follows.

The eighteen-link set in `GL6BW` is a **projected support**, not an already
declared active microscopic source-first collar.  Activating it changes the
finite problem: its locked subspace, source-free order-`h^4/h^6` effective
Hamiltonian, and stationary state must be recomputed.  `GL6BV` leaves the
full eighteen-link source-first `h^2/h^4/h^6` calculation explicitly open
and does not identify the selected `rho_+` with a state of that undeclared
parent.

## 8. Final scope verdict

The repaired theorem establishes a lawful finite-coupling route to nonzero
local-pair `T2` response, but only at the microscopic gap `2U_d` or as its
adiabatically eliminated contact.  The selected `E2` control lives at the
different scale `2J=O(h^6/U_d^5)` and uses a separately selected active
support/state.  Comparing those two numbers would mix Hamiltonians,
states, temporal structures, and coefficient orders.

Accordingly:

```text
pure-T source-before-Feshbach contact: PASS
common active parent and stationary E/T state: OPEN
eighteen-link h2/h4/h6 source-first calculation: OPEN
Ricci diagnostic in one completed causal coefficient: NOT TESTED
gravity or G: NOT DERIVED
```

## 9. Replay

The corrected target passes:

```text
GL6BV exact replay PASS (1394/1394)
PASS GL6BV packet checks 130/130
```

The independent reconstruction passes:

```text
PASS GL6BV independent hostile replay 2565/2565
```

`AUDITED_TARGETS.sha256` binds this verdict to every final target byte.  The
audit has its own manifest and seal.  The auditor made no target edit.

