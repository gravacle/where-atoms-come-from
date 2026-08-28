# Independent hostile audit -- q4/F3 ice `T2` Fisher-solder boundary

**Lane:** `GRA-FN-F3-Q4-ITFSB-V001`

**Audit date:** 2026-08-27

## Verdict

`ACCEPT_AFTER_SPECTRAL_DEPENDENCY_AND_FISHER_MAP_TYPING_REPAIR__EXACT_SECOND_JET_COMPLEMENT_NO_GO_GENERIC_RANK_AND_SYMMETRIZATION_RESULTS_SURVIVE__VECTOR_BACKGROUND_NOT_TENSOR_MODE__PHYSICAL_METRIC_SOLDER_AND_GRAVITY_OPEN`

The six-state calculation survives exact independent replay.  The first
Fisher derivative vanishes, the polarized second derivative is an invertible
`S4` intertwiner with the stated `A1/E/T2` eigenvalues, every
complement-preserving local diagonal background is first-order `T2` silent,
and the augmented covariance map has rank six exactly at generic nonzero
three-component mean.  The `T2` response is nevertheless the linearization
of a vector dyad, not an independently propagating tensor.

The hostile audit found two scoped custody/typing defects and repaired both
without changing the theorem.  First, the conditional pole interpretation
used the final `FL` result without pinning its theorem and audit bytes.
Second, “rank-six Fisher” could be misread as a six-parameter Fisher metric.
The repaired text defines the background-based exponential tilt explicitly
and identifies the rank as that of the map from background controls into the
six covariance components.  It does not count six source parameters or six
propagating modes.

## 1. Frozen dependency replay

The verifier now pins four load-bearing files:

```text
FK theorem                    cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98
FK independent audit          c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4
FL theorem                    98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452
FL independent audit          327bf6a4476c4c6382757dc156a96c6032233d34c25c1f7935e2582acf6c607a
```

Only `FK` is needed for the finite six-state theorem through Section 6.
`FL` is needed only for Section 7's conditional statement about photon poles
under `MAXWELL-IR`.  This prevents a circular promotion: neither the local
Fisher theorem nor the spectral screen assumes a graviton or physical metric
solder.

## 2. Six-state geometry and exact Fisher derivatives

In the orthonormal basis printed in the theorem, direct enumeration maps the
six two-in/two-out sign words to

\[
 \{\pm2e_1,\pm2e_2,\pm2e_3\}.
\]

For the odd exponential source, independent summation gives

\[
 Z(t)=2\sum_i\cosh(2t_i),\qquad
 m_i(t)={2\sinh(2t_i)\over\sum_j\cosh(2t_j)},
\]

and therefore

\[
 F(t)=4\operatorname{diag}\!\left(
 {\cosh(2t_i)\over\sum_j\cosh(2t_j)}\right)-m(t)m(t)^{\mathsf T}.
\]

At the origin, the exact moments are

\[
 \mathbb E[z_i z_j]={4\over3}\delta_{ij},\qquad
 \mathbb E[z_i z_j z_k]=0,\qquad
 \mathbb E[z_i z_j z_k z_l]
 ={16\over3}\,\mathbf1_{i=j=k=l}.
\]

The vanishing third cumulant proves all 27 components of `DF(0)=0`.
Subtracting the three pair contractions from the fourth moment reproduces
all 81 components of

\[
 D^2F_0[x,y]
 ={16\over3}\operatorname{diag}(x_i y_i)
 -{16\over9}(x\cdot y)I
 -{16\over9}(xy^{\mathsf T}+yx^{\mathsf T}).
\]

Its quadratic coefficient has trace
`-16 ||t||^2/9`, so it is nonzero for every nonzero one-direction source.
Thus the second order is genuinely the first nonzero order, not merely one
convenient nonzero component.

## 3. `S4` decomposition and failure of `O(3)` equivariance

On the exact decomposition

\[
 \operatorname{Sym}^2(V)=A_1\oplus E\oplus T_2,
\]

the Hessian operator has eigenvalues

\[
 \lambda_{A_1}=-{32\over9},\qquad
 \lambda_E={16\over9},\qquad
 \lambda_{T_2}=-{32\over9}.
\]

All are nonzero, so polarization gives rank six.  The verifier additionally
constructs all 24 coordinate-permutation actions on the ice-axis basis,
checks that they are distinct exact orthogonal signed-permutation matrices,
and verifies the intertwining identity on a six-element symmetric-matrix
basis.  Hence the positive classification is a genuine `S4` isomorphism,
not a dimension count.

The `O(3)` no-go is also exact.  A 45-degree rotation in the `e_1,e_2`
plane conjugates `diag(1,-1,0)` in `E` to the symmetric `(12)` off-diagonal
matrix in `T2`.  An `O(3)` intertwiner would give both the same `ell=2`
eigenvalue.  The values `16/9` and `-32/9` differ, so the Hessian exposes the
tetrahedral axes and is not `O(3)` equivariant.

The six-dimensional object is the polarized bilinear second jet.  A single
source `t` still supplies only the three-dimensional quadratic image of
`t t^T`; polarization does not manufacture a six-dimensional linear state
tangent.

## 4. Complement boundary and uniqueness of the local vector query

For arbitrary interior probabilities `p_i^+,p_i^-`, put

\[
 w_i=p_i^++p_i^-,\qquad m_i=2(p_i^+-p_i^-).
\]

Direct summation gives

\[
 F(p)=4\operatorname{diag}(w)-mm^{\mathsf T}.
\]

The repaired theorem makes the typing explicit: this is the covariance of
the fixed one-link statistic, equivalently the Fisher/commuting-SLD-QFI at
zero source of the tilt `p(z)e^{theta dot z}/Z_p(theta)`.  It is not the
intrinsic Fisher metric on the five-dimensional probability simplex.

Complement preservation is exactly `m=0`.  Then every normalized tangent
obeys

\[
 dF=4\operatorname{diag}(dw),\qquad \sum_i dw_i=0,
\]

whose image is exactly the two-dimensional `E` sector.  An independent
scalar adds only `A1`, raising the rank to three and leaving `T2` absent.
The six-state permutation character `(6,2,2,0,0)` contains the vector
character `(3,1,-1,0,-1)` once, so there is no second equivariant local
diagonal vector statistic that evades the same covariance boundary.

This no-go does not cover nonlocal or noncommuting queries, and the theorem
does not say otherwise.

## 5. Generic broken-background rank

At a general interior background,

\[
 dF=4\operatorname{diag}(dw)
 -(m\,dm^{\mathsf T}+dm\,m^{\mathsf T}).
\]

The map from `dm` to the three off-diagonal entries has matrix

\[
 -\begin{pmatrix}
 m_2&m_1&0\\
 m_3&0&m_1\\
 0&m_3&m_2
 \end{pmatrix},
\qquad \det=2m_1m_2m_3.
\]

Consequently its `T2` projection is invertible exactly when all three mean
components are nonzero.  The three mean directions plus two normalized
weight directions then give rank five.  Adding an independently controlled
`rho I` gives rank six because the weight directions and `rho` span all
diagonal matrices.  If any `m_i` is zero the augmented determinant vanishes;
at `m=0` the ranks are two without the scalar and three with it.

This full rank is algebraic but becomes ill-conditioned continuously as
`m_1m_2m_3` approaches zero.  More importantly, its off-diagonal response is

\[
 \delta F_{T_2}
 =-\operatorname{offdiag}(m\,\delta m^{\mathsf T}
                          +\delta m\,m^{\mathsf T}),
\]

so it is a fluctuation of a vector mean in a symmetry-broken background.
The F3 parent has not derived preparation or stabilization of that
background, and the scalar remains independently supplied.

## 6. Complement-mixture cancellation and pole interpretation

Complement sends `m` to `-m` while leaving the conditional covariance
unchanged.  For the unlabelled mixture `bar p=(p+Cp)/2`, the law of total
covariance gives exactly

\[
 F(\bar p)=F(p)+mm^{\mathsf T}=4\operatorname{diag}(w).
\]

Thus retaining the sign label retains a conditional dyad only by retaining
an external classical record.  Discarding the label restores complement
symmetry and cancels every off-diagonal `T2` entry.  Randomization does not
provide a hidden symmetric-state solder.

The final spectral interpretation agrees with the separately audited `FL`
screen.  At zero mean the nonzero second jet is an even two-vector
susceptibility, which at the Gaussian Maxwell fixed point belongs to a
two-photon channel rather than an isolated one-particle tensor pole.  At
nonzero mean its linear variation inherits the pole of the same vector
fluctuation, which under `MAXWELL-IR` is spin one.  Neither construction
earns helicity two, a tensor Ward identity, or a graviton.

## 7. Promotion ceiling

The strongest safe statement is:

> The exact ice one-link Fisher metric has a full tetrahedral `A1+E+T2`
> polarized second jet, but its `T2` component is the dyad of a
> complement-breaking vector mean.  Complement-preserving local diagonal
> states are first-order `T2` silent.  A generic broken vector background
> plus an independent scalar gives a rank-six covariance response, not an
> independent tensor mode or graviton.

Still open are physical preparation and stabilization, same-parent scalar
ownership and calibration, an `O(3)`-covariant metric response, a protected
tensor pole with TT residue and Ward identities, refinement/gluing,
universal stress coupling, RGRL-B, Einstein dynamics, gravity, and numerical
`G`.  No fitted interaction, unowned scalar, programmed background, or
conditional label can be promoted as closure.

## 8. Reproduction and canonical hashes

Run:

```text
python3 LANE_GRA_FN_F3_Q4_ICE_T2_FISHER_SOLDER_BOUNDARY_V001/verify_ice_t2_fisher_solder_boundary.py
```

Expected result: `SUMMARY 58/58 PASS`.

Canonical hashes before adding this audit and the regenerated transcript to
the manifest:

```text
THEOREM.md                                be69f15d611827db9841bd932042604deb4f82a777ff9da28b80e4493cef7596
verify_ice_t2_fisher_solder_boundary.py   1a3007f06b9f34aaf707b906c11973c17daf0c81569c3082c7aa8760c92692ae
SELF_AUDIT.md                             6e9243f4621045749778bd6251fcab18434998126c8d99d119a142290b821a4e
```

The manifest pins the final theorem, verifier, preserved builder self-audit,
this independent audit, and the regenerated verification transcript.
