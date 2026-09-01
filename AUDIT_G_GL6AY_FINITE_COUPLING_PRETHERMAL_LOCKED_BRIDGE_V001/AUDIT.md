# Independent hostile audit — GL6AY finite-coupling prethermal locked bridge

**Frozen author theorem SHA-256:**
`84c8b221974a8e463e381c40ab4827fc805ee9277f39002e25d23e4d92318cca`

**Frozen author manifest SHA-256:**
`d51c6aea006c1b5cdc7a75023dfd59e0cdb363549565c92aab0a5e1ee5083710`

**Frozen author seal-file SHA-256:**
`307b57e94c2054bc19c81a44b64f6ae10d7d68bcb3937bfc7005bd8caec3492d`

**Disposition:**
`FAIL__REPAIR_REQUIRED__GLOBAL_LOCK_PROJECTOR_AND_DRESSED_SUBSPACE_SCOPE`

No author byte was modified.

## 1. Custody and independence

All twelve author files are frozen in `AUDITED_TARGETS.sha256`.  The author
manifest resolves to all ten manifested content files, and its seal resolves.
Both author verifiers pass in normal and optimized Python modes.  The
independent replay imports no author module.

The three external theorem inputs were checked against the exact versioned
arXiv PDFs named by the author:

- `1509.05386v3`, Theorems 3.1--3.3 and equations (3.1)--(3.2);
- `1704.08703v2`, Appendix A; and
- `1105.0675v1`, Section 4.

The source URLs are version-pinned but their PDF bytes are not vendored in
the packet.  This audit records theorem/version custody, not an offline hash
of the external PDFs.

## 2. Primary theorem and F3 mapping — pass

The exact F3 coarse tensor factorization is legitimate: each cell owns four
link qubits, so its local dimension is sixteen.  The parent constraint is
one-cell supported and each child constraint is supported on four cells.
The commuting local terms

```text
q_v^2=(k_v-2)^2
```

have spectrum `{0,1,4}`.  Hence their sum `N_def` has integer spectrum and
`exp(2 pi i N_def)=I`.

For a flip `(x,a)`, the set

```text
S_(x,a)={x+d_a-d_b:b=0,1,2,3}
```

contains the full parent support and full child support of the only two
constraint terms that fail to commute with the flip.  It is connected in
the bounded-range `A3` cell adjacency, has four cells, and a fixed cell lies
in exactly sixteen translated labeled flip supports.  The crude bounds

```text
d_0<=16 exp(4kappa_0)h,
v_0<=32 exp(4kappa_0)h
```

are therefore valid potential-norm bounds, not global operator bounds.

The Else--Fendley--Kemp--Nayak strong-support extension applies because the
constraint summands commute and have uniformly finite support.  It restores
support preservation under `exp(i theta N_def)` and preserves the ADHH
commutator calculus.

Direct comparison with ADHH v3 verifies the author constants and hypotheses:

```text
nu_0=(54pi/kappa_0^2)(d_0+2v_0),
n_*=floor{(U_d/nu_0)/(1+ln(U_d/nu_0))^3}-2,
U_d>=9pi v_0/kappa_0,
kappa_*=kappa_0/[1+log(n_*+1)].
```

The exact transformed identity, `[D_hat,N_def]=0`, remainder factor
`(2/3)^n_* v_0`, and potential conjugation estimate match Theorem 3.1.

The charged-resonance correction is also necessary and correct.  A flip has

```text
Delta N_def=2delta(q_u+q_v)+2.
```

Both zero-frequency resonance classes occur outside `P`; only `P D_0 P=0`
is true.

## 3. Material defect one — global `P` is not an infinite interaction

Equations (AY.24)--(AY.26), plus the README and RESULT, promote

```text
Phi(S)=P D_hat(S)P,
P=chi(N_def=0),
```

to an infinite-volume interaction.  This is not well typed.  In infinite
volume `N_def` is extensive and no global spectral projector `P` is an
element of the quasi-local algebra.  A sequence of finite-volume projectors
does not converge in operator norm to a quasi-local projection.

The affine endpoint argument is sound, but it must be attached to a local
representative rather than to a nonexistent global `P`.

### Sufficient local repair

For each strongly supported term `D_hat(S)`, define

```text
N_S=sum_(v:supp(q_v^2) subset S) q_v^2,
P_S^0=chi(N_S=0)=product_(v:supp(q_v^2) subset S)chi(q_v^2=0),
Phi_S=P_S^0 D_hat(S) P_S^0.                              (AUD.1)
```

Strong support gives

```text
[D_hat(S),q_v^2]=0 whenever supp(q_v^2) is not contained in S.
```

Termwise pinching gives `[D_hat(S),N_def]=0`; subtracting the outside
constraint sum yields `[D_hat(S),N_S]=0`.  Hence a globally locked input is
kept locked by `Phi_S`: contained constraints remain in the nonnegative
zero eigenspace, and every noncontained constraint is individually
unchanged.  `Phi_S` is finite-support and needs no global projector.

On globally locked endpoint configurations, the GL6AX affine theorem then
gives exact conservation of all four port totals.  This is conservation of
the restricted locked interaction/local `U(1)` action; no full-Hilbert port
symmetry is implied.

The repair preserves the twist estimate.  `P_S^0` is diagonal in link
occupation and commutes with the local twist generator `A_S`, so

```text
[A_S,[A_S,Phi_S]]
 =P_S^0[A_S,[A_S,D_hat(S)]]P_S^0.                       (AUD.2)
```

Projection does not increase norm.  The author fourth-moment estimate and
wrapping-tail split therefore survive with `Phi_S` substituted for
`P D_hat(S)P`.

## 4. Finite second moment and GL6AX interface — pass after local typing

For connected `S` with `m=|S|`, bounded coordination gives coordinate
diameter `O(m)` and at most `O(m)` port-zero links.  Thus

```text
||A_S||<=C_geo m^2,
||[A_S,[A_S,Phi_S]]||<=4C_geo^2m^4||D_hat(S)||.
```

The exponential strong-support potential norm dominates this polynomial.
The author's use of `sup_m m^4 exp(-kappa_*m)` is loose but valid.  A more
economical count would use a cubic power after the per-cell incidence count;
the displayed quartic bound remains an upper bound.

The torus wrapping tail has minimum connected support proportional to the
wrapped period.  Fixed positive `kappa_*` therefore gives exponential tail
suppression up to polynomial volume.  The effective GL6AX anisotropic
dichotomy follows for the locally typed locked interaction, subject to the
same centered-sector and finite/GNS ceilings.

## 5. Order-six interface — pass

The normal-form unitary is perturbatively near identity.  Since all
locked-sector terms below order six are scalar, an allowed near-identity
within-`P` gauge cannot change the first nonscalar order-six matrix element.
The coefficient

```text
-(63/8)(h^6/U_d^5) sum_c T_c
```

therefore matches sealed GL6AO as a formal coefficient.  The author correctly
does not turn this match into convergence of the full series.

## 6. Local-observable horizon — pass

ADHH Theorem 3.3 states, in dimension `d`,

```text
0<r_1<ln(3/2)/(d+1),
||tau_t^H(O)-tau_t^(nu N+D_hat)(O)||<=K_3(O)/nu,
t<=exp(r_1 n_*).
```

With `d=3` and `nu=U_d`, AY.33--AY.35 match the source.  `K_3(O)` is
volume-independent and independent of the large frequency with the other
model parameters fixed.  The author correctly leaves `U_d/h`, clock
calibration, and the record horizon open.

## 7. Material defect two — no global `Y^*P` closeness follows

Section 8 says that the deviation of the dressed space `Y^*P` from the bare
space is controlled by AY.20.  It is not.  AY.20 is a strong potential-norm
bound, with a corresponding operator-norm statement for a fixed local
observable.  It supplies no volume-uniform estimate on

```text
||Y^* P Y-P||
```

for the extensive finite-volume spectral projector, and there is no global
`P` to conjugate in the infinite quasi-local algebra.

The failure is not merely formal.  A product of identical small local
rotations can be arbitrarily close to identity on each fixed local observable
while the rank-one product projection rotates to norm distance

```text
sqrt(1-cos(epsilon)^(2V)) ->1.
```

The repaired theorem may claim local observable/potential dressing and the
finite-volume algebraic identity.  It must delete every global projector or
subspace-closeness inference.

## 8. Local leakage versus winding return — pass

`P V_hat P=0` is an exact finite-volume pinching identity, not an invariant
space theorem.  A local remainder term can take locked data to a charged
state without winding.  If a finite contractible process returns to locked
endpoints, affine port conservation applies.  On a torus, a port-changing
locked-to-locked return requires a noncontractible symmetric difference of
at least `2L_min` links.  The author's distinction is exact once finite and
infinite projector typing is kept separate.

## 9. First winding coefficient — pass

For the simple alternating `j/3` cycle of length `r=2L_j`, every proper
nonempty subset has positive even boundary `b(S)`.  Its excitation energy is
`U_db(S)`, with `2<=b(S)<=r`.  At the first endpoint-changing order, each
changed link appears exactly once and no intermediate subset is locked.
Hence no folded or lower-endpoint term contributes, all resolvent words have
the same negative sign, and

```text
A_C^(r)
 =-h^r sum_pi product_(s=1)^(r-1)[1/(U_d b(S_(pi,s)))],
r!h^r/(rU_d)^(r-1)<=|A_C^(r)|<=r!h^r/(2U_d)^(r-1).
```

Independent dynamic programming and explicit permutation replay agree.
This proves a nonzero leading coefficient on each fixed torus, not a
volume-uniform analytic radius.

## 10. Whole-band boundary — pass with stated ceiling

The locked-to-charged gap is `2U_d`; `N_def=1` is forbidden by equality of
the total parent and child degree deviations and `N_def=2` is attained by one
link flip.  Since all `X_e` commute,

```text
||W_L||=h|E_L|.
```

A raw global Neumann/direct-rotation sufficient condition therefore scales
as `h|E_L|<cU_d` and collapses with volume.  This is exactly the proof
boundary described in Bravyi--DiVincenzo--Loss Section 4: fixed-order linked
coefficients remain meaningful while the standard global Taylor/direct
rotation construction need not be uniform.

The author correctly stops short of an impossibility theorem for every
special algebraic model.  The global-norm result is only a failure of the
standard whole-band proof route.

## 11. Required author repair

Before promotion, the author packet must:

1. replace every infinite-volume `P D_hat(S)P` interaction with the local
   construction (AUD.1), or an exactly equivalent finite-collar definition;
2. restate port conservation as a property of the resulting interaction
   restricted to globally locked configurations/local locked algebra;
3. propagate this local definition through the wrapping-tail and `D_2`
   formulas, using (AUD.2);
4. delete the AY.20 inference about global `Y^*P` or global subspace norm;
5. limit dressing claims to source-authorized local observables/potentials
   and finite-volume algebraic identities; and
6. update README, RESULT, SELF_AUDIT, verifier tokens, manifest, and seal,
   followed by a fresh independent hostile audit.

The core finite-coupling prethermal bridge is repairable without a new
theory decision.  The frozen V001 packet is not promotable as written.
