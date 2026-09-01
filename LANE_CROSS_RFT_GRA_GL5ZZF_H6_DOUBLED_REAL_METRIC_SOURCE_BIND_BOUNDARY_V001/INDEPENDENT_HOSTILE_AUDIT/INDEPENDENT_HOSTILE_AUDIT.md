# Independent hostile audit — GL5ZZF doubled-real metric-source boundary

**Date:** 2026-08-30  
**Frozen theorem:** `2f9fbdde026765b5a4dc335d6b87777e1042efcace1bd4d6482ce0f8ac235b22`  
**Author replay:** `PASS__H6DRMSB__81/81`  
**Independent replay:** `PASS__HOSTILE_H6DRMSB__99/99`

## Verdict

**CLEAN WITHIN THE MEMBER-IDENTIFIABILITY BOUNDARY.**

The theorem proves that the same finite source-free H6 operator at fixed
`J_6` does not, by itself, select a unique metric-source derivative.  It also
constructs exact Hermitian real quadratures from the conjugate FY momentum
modes and one exact complement-even doubled candidate.  It does not prove
that two distinct microscopic F3 parents are the same theory, that the copied
source is the metric derivative of one common parent, or that any physical
metric--memory crossed kernel has been calculated.

## Exact two-member comparison

The hostile verifier independently recomputes

\[
 a(x)=\frac{8}{63x^6}
 \left(1-x^2-\frac{37}{12}x^4-\frac{16247}{900}x^6\right)
\]

at the two declared FY samples and obtains

\[
 a(2/5)=\frac{2415673}{113400},\qquad
 a(1/2)=\frac{31706}{14175},\qquad
 a(2/5)-a(1/2)=\frac{3203}{168}>0.
\]

For either `x`, setting

\[
 U_d(x)=\frac{8J_6}{63x^6},\qquad h(x)=xU_d(x)
\]

gives positive but different microscopic parameter pairs with the same
coefficient in

\[
 H_6=-J_6\sum_C B_C.
\]

This is equality of the finite source-free effective H6 operator.  It is not
equality of the microscopic parents or their higher-order terms, boundary
kernels, preparation laws, or metric couplings.

## Scope repairs and the normalized FY arrays

The pre-audit draft dropped the FY truncation superscript in the difference
equation and used “complete metric source” without qualification.  The frozen
theorem now correctly retains

\[
 \frac{Q_{\rm eff}^{(\le6)}(1;x)}{J_6}=a(x)D_1+R_1
\]

and states that the source is complete only through H6 under `FV-PURE`.

The original draft also stated only that `D_1` and `R_1` are linearly
independent.  That is not enough to cancel `R_1` between two values of `x`.
The frozen theorem now pins the stronger FY fact: on the same frozen FO
basis,

\[
 D_1=Q_{\rm pair}/U_d,qquad R_1=Q_{\rm ring}/J_6
\]

are fixed dimensionless arrays constructed before sampling `x`.  `D_1` is
diagonal, while `R_1` has a certified nonzero off-diagonal entry, so they are
also linearly independent.  Therefore

\[
 Q_{\rm eff}^{(\le6)}(1;2/5)
 -Q_{\rm eff}^{(\le6)}(1;1/2)
 =J_6\frac{3203}{168}D_1\ne0
\]

is exactly typed within the sealed FY scope.

## Real source and complement copy

From `Q_29=Q_1^dagger`, the theorem defines

\[
 Q_c=\frac{Q_1+Q_{29}}{\sqrt2},\qquad
 Q_s=\frac{Q_1-Q_{29}}{i\sqrt2},\qquad
 Q_\phi=\frac{e^{-i\phi}Q_1+e^{i\phi}Q_{29}}{\sqrt2}.
\]

The independent replay evaluates generic non-Hermitian test matrices and
confirms that all fixed real-profile combinations are Hermitian.  The
`1/sqrt(2)` real-quadrature factor and FY's inherited `1/sqrt(60)` Fourier
normalization are distinct and both load bearing.

If `F:H_+ -> H_-` is the fixed occupation-complement isometry, then

\[
 Q_\phi^{\rm dbl}=Q_\phi^+\oplus FQ_\phi^+F^\dagger
\]

is exactly Hermitian and complement even.  The ZZA law has 17 target rows and
total probability `1/2` in each component, so the doubled candidate has the
right finite carrier scope.  This is nevertheless only a prospectively
selected source extension.  The theorem correctly leaves open whether it is
the derivative of a common microscopic F3 metric coupling and retains all
projector/isometry derivatives if they are source dependent.

## Source duality and CTP sign

The FY convention gives

\[
 V=\frac{\partial H}{\partial g}=-\frac12Q.
\]

The hostile replay checks that a Hamiltonian kick `H -> H+gV` produces the
later observable derivative `+i[V,Q]`, so the sign of the displayed bulk
response and its decomposition is correct.  The closed-time-path insertion

\[
 S_g=\sum_{\sigma=\pm}\sigma\int dt\,g^\sigma V^\sigma
\]

has the corresponding branch sign.

The flattened symmetric-tensor source convention is also exact.  The
off-diagonal source coordinates carry factors of two because the Frobenius
sum includes both `(i,j)` and `(j,i)`.  A bare six-by-six representation
identity is therefore not a physical dual-source map.

## What is and is not identified

The bulk commutator can be organized into two fixed response arrays,

\[
 \widehat L^{\rm bulk}
 =-\frac{J_6}{2}\left[a(x)L^D+L^R\right].
\]

This is not the complete CTP/Schur kernel.  The theorem explicitly retains
the derivatives of the prepared state and conditional normalization,
selected/complement projectors, complement isometry, profile basis,
operators, router/work/failure/boundary coordinates, and every nuisance
direction outside the six chart coordinates.  Those terms belong to one
fixed parent and are not adjustable contacts.

Accordingly, the exact result is an identifiability stop:

- source-free H6 plus the writer and memory response do not select a unique
  metric source;
- one independent parent law must select `x`, `phi`, and the doubled source,
  or `a(x)` must remain symbolic;
- no `GL5ZF` row, post-Schur factorization, Ricci form, gravity theorem, or
  value of `G` follows from this packet.

No material defect remains within that claim boundary.
