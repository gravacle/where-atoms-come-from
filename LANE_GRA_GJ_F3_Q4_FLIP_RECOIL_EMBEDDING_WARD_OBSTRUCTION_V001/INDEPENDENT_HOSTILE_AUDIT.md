# Independent hostile audit: FRSEWO

**Date:** 2026-08-29  
**Disposition:** `PASS_WITH_STRICT_CONDITIONAL_CEILING`

No theorem-invalidating defect was found in the narrow claim.

## Independent reconstruction

The audit independently checked the exact operator typing rather than relying
on the new theorem's prose.

1. GD gives
   `W^dagger P_L W=pI-(hbar kappa/2)Z` and
   `W^dagger P_R W=pI+(hbar kappa/2)Z`.  Any scalar local, support, or Fourier
   weighting of those bare encoded factors remains configuration diagonal.
   Therefore `diag([H,D])=0` for every finite `H` and configuration-diagonal
   `D`.
2. FZ05 is one exact nonzero configuration-diagonal matrix-entry witness under
   its supplied embedding contraction.  The witness is nonzero in
   `Q(zeta_240)`, and at `x=2/5` its complete FY coefficient is
   `rho f_E=2415673/113400 != 0`.  FY's ring source is strictly off diagonal
   and cannot cancel that entry.
3. GD's off-diagonal recoil operator is a factor-edge flux, not a derived
   `T^{0j}`.  Every single-link flip exits the local ice fiber, so its direct
   ice projection vanishes.  A source-derived or Feshbach-dressed effective
   `T^{0j}` remains open.
4. `O(j^2)` and mixed `h_0i j_kl` contacts have zero first insertion when all
   sources are off.  An active boundary or another nonzero first-source term
   changes the parent.

## Replay

The final lane verifier was replayed independently:

```text
python3 LANE_GRA_GJ_F3_Q4_FLIP_RECOIL_EMBEDDING_WARD_OBSTRUCTION_V001/verify_flip_recoil_embedding_ward_obstruction.py
SUMMARY 42/42 flip-recoil embedding-Ward checks passed
CEILING conditional supplied-embedding obstruction only; native physical divergence and complete source-before-Feshbach construction remain open
```

An earlier independent working count was `65/65` because it counted the
twenty-four single-link projection cases separately.  The retained verifier
consolidates those same predicates into one all-cases check, reducing the count
by twenty-three without changing the test.

FZ's exact-fast verifier was also replayed:

```text
python3 LANE_GRA_FZ_F3_Q4_M1_CONTINUITY_CONTACT_WARD_BOUNDARY_V001/verify_continuity_contact_ward_boundary.py
SUMMARY 49/49 continuity/contact/Ward/boundary checks passed
```

## Required promotion ceiling

The admissible claim is:

> The bare, directly projected, scalar-weighted GD `P_L/P_R` density cannot
> close the frozen FY source-off identity under FZ's supplied embedding
> contraction, or under a native divergence later proved equivalent on that
> source.

It must not be promoted to an unconditional statement that "flip recoil alone
cannot close" or to a proof that pair-field/support momentum is logically
necessary.  Dynamical position-weighted localization, interaction
contributions to `T^{0j}`, source-before-Feshbach dressing, native `Delta_m`,
active boundaries, and modified spatial sources are untested.

## Upstream custody warning

GD's primary verifier currently stops on a stale hash of
`GRAVITY_NO_LAB_PROOF_TASK_PLAN_V001.md`: the expected hash begins `875797`,
while the current hash begins `2f43cb`.  The present lane directly pins and
verifies GD's unchanged `THEOREM.md` and `RESULT.md`, so this is not a
mathematical defect in FRSEWO.  It is an upstream packet-integrity warning and
must not be silently represented as a passing replay of GD's complete packet.

## Verdict

Promote FRSEWO only with the stated supplied-embedding, frozen-source, and
bare-direct-density ceiling.  It closes RF3a's bounded falsification test.  It
does not close RF3, gravity formation, Newtonian gravity, or `G`.
