# Distinct hostile audit — GL6AZ record-authenticated prethermal mission

**Target:** `LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/`  
**Frozen theorem SHA-256:** `bcca352d9e58deba63a068a29f94a0e19fb88ae9994391ab1b5053120587ba44`  
**Frozen author-manifest SHA-256:** `a932e083b5dd629e41f0014d52bf1c65f982612b55314e1e03e082b0feb8ebc7`  
**Frozen author-seal-file SHA-256:** `620eb8bcb254b0296bd4d2f7e81a0b3e78c0d0b9ce80ad257e97a4db394f392d`  
**Disposition:** `PASS__REPAIRED_HIGH_BRANCH_PRETHERMAL_APPLICATION_AND_RECORD_MARGINAL_IDENTIFIABILITY__SELECTED_MEMBER_AND_MISSION_DATA_OPEN`

## 1. Custody and defect history

All thirteen repaired author files are byte-pinned in
`AUDITED_TARGETS.sha256`.  The independent replay imports no author module.
The author and audit verifiers run identically in normal and optimized Python
modes.  All thirty-three declared GL6AY, GL6AM, GL6V, GL6AN, GL5ZZF, and U-DCL
dependency objects resolve at their frozen hashes.

The first hostile pass found a material source-domain omission before this
audit packet was created.  The compact ADHH theorem statement lists the first
smallness inequality and `n_*>=1`, but its proof in section 4.3 also starts
from `nu_0<=nu`.  Without that condition the formula has a spurious low-ratio
branch near `nu/nu_0=e^{-1}`.  The author repaired the theorem, ledgers,
verifiers, and continuation before resealing.  This audit starts from and
passes only those repaired bytes.  No author byte was altered by the audit.

## 2. Dimensionless source-theorem application

Writing

```text
R=U_d/h,  Hbar=H/h,  s=h t_phys/hbar
```

gives the exact unitary identity

```text
exp(i H t_phys/hbar) O exp(-i H t_phys/hbar)
 = exp(i Hbar s) O exp(-i Hbar s).
```

Pinching and the strong potential norm are homogeneous, so the normalized
ADHH frequency is `R`, the normalized interaction is fixed, and the local
error is `Kbar_3(O)/R`.  The complete proof-licensed hypotheses are

```text
R >= 9 pi vbar_0/kappa_0,
R >= nubar_0,
n_*(R) >= 1,
0 < r_1 < ln(3/2)/4.
```

The restored middle inequality selects `x=R/nubar_0>=1`.  An explicit hostile
counterexample to the unrepaired wording uses `kappa_0=1/4` and
`x=e^{-0.9}`: it passes the first inequality and gives `n_*=404`, yet violates
`R>=nubar_0`.  The repaired high branch removes this false license.

## 3. Three distinct sufficient-domain floors

The one-flip matrix element below forces
`vbar_0>=exp(4 kappa_0)`.  Therefore the first smallness inequality alone has
the exact minimum

```text
R >= 36 pi e = 307.4304320162484       at kappa_0=1/4.
```

The restored scale-separation inequality separately implies

```text
nubar_0 >= 108 pi exp(4 kappa_0)/kappa_0^2
          >= 432 pi e^2 = 10028.190682380982
```

at `kappa_0=1/2`.  On the high branch, `n_*>=1` is equivalent to

```text
x/(1+ln x)^3 >= 3.
```

The unique root above `e^2` is
`x_*=1861.32559690908`.  Thus the complete source-proof domain has the
universal lower floor

```text
R >= x_* 432 pi e^2 = 18,665,728.0078.
```

Exact member norms can only strengthen this sufficient requirement.  The
sample ratios `R=2` and `R=5/2` fail even the first layer.  None of these
inequalities is necessary for the physical phase, and their failure proves no
instability, phase transition, or loss of recordhood.

## 4. Native ratio and non-identifiability

For a locked computational configuration and a chosen link flip, both
endpoint charges change from zero to `+/-1`.  Hence

```text
N_def(n^e)-N_def(n)=2,
|<n^e|H|n>|=h,
<n^e|H|n^e>-<n|H|n>=2 U_d,
R=Delta_def/(2 A_X).
```

`Delta_def` is a computational-basis diagonal separation, not a spectral gap
of the finite-`h` Hamiltonian.  The independent replay exhausts every local
locked completion.  It also reconstructs every positive `R` with the licensed
integer witness `d_*=3`, `Delta=4U_d`; no fractional `d_*` is used.  Thus the
bare lock ray admits every positive `R` and does not select one.

The exact GL5ZZF pair is also reproduced.  At fixed `J_6`, `x=2/5` and
`x=1/2` give the same source-free H6 coefficient but different microscopic
parameter pairs, and

```text
a(2/5)-a(1/2)=3203/168.
```

This is an exact descendant non-identifiability witness, not a statement that
the full microscopic members are identical.

## 5. Grouped strong support cannot cancel

For link `e=(x,a)`, the child support

```text
S_e={x+d_a-d_b : b=0,1,2,3}
```

has four distinct cells on the infinite native lattice and on injectively
lifted noncollapsed quotients.  Four incoming links to a common child can
share the same potential support, so the lower bound must survive grouping.
It does: between `|n>` and its specific one-link flip `|n^e>`, every other
link flip reaches an orthogonal configuration, while the `N_def` pinch removes
the charged matrix element.  The grouped coefficient therefore retains

```text
|<n^e|V_0(S_e)|n>|=h,  ||V_0(S_e)||>=h.
```

This proves the stated strong-norm lower bound.  The audit does not extend
“exactly four” to arbitrarily collapsed finite quotients, where GL6AY says
only “at most four.”

## 6. Authenticated read and clock typing

The complete four-link read has sixteen outcomes.  Coarsening by
`m_ab=q_a q_b` gives an exact eight-plus/eight-minus binary PVM.  For any
common system state on the selected ready/`MATCH` factor,

```text
D_TV = |Delta <M_beta>|/2 <= Kbar_3(M_beta)/(2R).
```

This is the binary pair marginal of a complete flag-retaining instrument.  It
does not bound total variation of all status, failure, and authentication
outputs.  No flag is discarded and no `MATCH` postselection is performed;
the selected factor is a stated premise.  The observable `M_beta` is not
renamed a record—the retaining terminal instrument owns the record.

The hold endpoint `t_Q` is typed as the system sampling time immediately
before the terminal read dilation.  The common read is then performed with
the response off or refocused, so its pulse duration is not falsely charged
to evolution under the F3 generator.  A same-parent clock bind still requires
the physical formation-completion and sampling timestamps plus proof that the
declared Hamiltonian governs the interval.  U-DCL supplies no numerical
`h`-scaled duration or equation tying `R` to that mission, although an
individual witness may retain clock data.

## 7. Exact continuation fork and ceilings

Calibration does not automatically close the theorem.  It decides one of two
branches:

1. **Inside-domain branch.**  If the selected physical member and mission
   satisfy every repaired source inequality and the horizon, compute the
   exact normalized norms and `Kbar_3`, then apply the binary marginal bound.
2. **Outside-domain branch.**  If any source inequality fails, GL6AY cannot
   certify that mission.  A sharper finite-mission local-observable theorem or
   a direct controlled calculation is required; clock binding and
   extrapolation of `n_*` cannot repair an out-of-domain proof.

This packet does not select a numerical `R`, a mission duration, a state, an
all-time phase, gravity, or `G`.  It introduces no graviton, Ricci template,
Einstein endpoint, or gravitational calibration.

Within those exact ceilings, every post-repair hostile attack survived.

**Hostile verdict: PASS.**
