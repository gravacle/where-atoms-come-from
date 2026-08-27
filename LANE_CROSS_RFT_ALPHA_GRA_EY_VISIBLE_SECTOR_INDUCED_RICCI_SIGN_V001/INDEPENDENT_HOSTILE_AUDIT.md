# Independent hostile audit -- VSIRS

**Date:** 2026-08-27

**Auditor verdict:**
`PASS__SOURCE_FREEZE_AUTHORIZED__EXACT_CONDITIONAL_ONE_LOOP_PROPER_TIME_CONTRIBUTION_ONLY__NOT_GRAVITY`

**Frozen scientific bytes audited:**

- `THEOREM.md`:
  `21c86b24025b9393008c2975c6f421146d5688da0342c19c16cf21bbab4a35b4`
- `verify_visible_sector_induced_ricci_sign.py`:
  `e14ce25b533b1becbe993d4137f6fe2c57b9a7111abfba31ac8d646f1d51677a`

Independent replay on those source bytes returned `PASS 21/21`.

## Hostile scientific findings

1. **The Visser-to-RIEHB sign dictionary is explicit and sound.**  Visser's
   convention writes the Einstein--Hilbert term as
   `-R_V/(16 pi G)` and gives
   `1/G=-(str k1) Delta/(2 pi)`.  The declared conversion
   `R_R=-R_V` therefore gives the RIEHB coefficient
   `C_R=-(str k1) Delta/(32 pi^2)` multiplying `+R_R`.  The theorem forbids
   importing an unmapped nonminimal-coupling convention.
2. **The statistics sign is now inserted exactly once.**  Visser's table
   value `k1(Weyl)=-1/6` is the traced Laplace-operator coefficient before
   the supertrace weight.  Equation EY01a supplies the single additional
   `(-1)^F`; the vector coefficient already owns its gauge/ghost census and
   is not reweighted as a fermion.  This convention is also consistent with
   Visser's net chiral-supermultiplet table entry.
3. **The visible field census is correct.**  Per generation the Standard
   Model contains `6+3+3+2+1=15` two-component Weyl fields without a
   right-handed neutrino, hence 45 over three generations.  The unbroken
   gauge basis has `8+3+1=12` massless vectors, and one complex Higgs doublet
   supplies four real scalar components.  The Weyl and vector table entries
   already sum their physical spin states, so neither multiplicity is
   doubled.
4. **All corrected arithmetic and normalization factors rederive.**  The
   exact declared-shell supertrace is
   `str k1=4(1/6-xi_H)-45(-1/6)+12(-2/3)=1/6-4 xi_H`.
   Hence positivity holds exactly when `xi_H>1/24`.  Minimal coupling gives
   `C_R=-Delta/(192 pi^2)`, while table-conformal coupling gives
   `str k1=-1/2` and therefore
   `C_R=+Delta/(64 pi^2)`.  The latter factor is `1/64`, not `1/32`.
5. **The additional-scalar control is correct.**  For `N_p` admitted minimal
   real scalars, `str k1=(1+N_p)/6-4 xi_H` and positivity requires
   `24 xi_H>1+N_p`.  Six such modes at minimal Higgs coupling give
   `str k1=7/6` and `C_R=-7 Delta/(192 pi^2)`.  The theorem does not promote
   PMMDC coordinates to propagating scalar determinants without a separate
   physical proof.
6. **The chiral and electroweak typing is adequate.**  The Weyl determinant
   census is attached to the complete anomaly-free Standard-Model
   representation set, including the mixed gauge--gravitational and global
   `SU(2)` guards.  The calculation is explicitly a zero-temperature
   ultraviolet symmetric-gauge-basis proper-time window with its infrared
   end above the resolved electroweak thresholds, not a thermal restoration
   theorem or an exact momentum-shell projector.
7. **The alpha ancestry statement is corrected.**  Hypercharge `B` and
   neutral weak `W3` are jointly identified as the photon ancestors, so the
   lane does not incorrectly attribute low-energy alpha to hypercharge alone.
   Alpha is not counted as an extra species or stress source.
8. **The scientific ceiling is intact.**  The result is regulator- and
   matching-prescription conditional, one loop, and a visible-field
   contribution on an owned `xi_H` domain.  A positive total Newton
   coefficient additionally requires the complete ultraviolet spectrum,
   thresholds, measures, remainders, and the distinguished zero-bare-term
   matching rule.  Record-to-metric soldering, refinement, constraints,
   complete stress, the cosmological term, same-parent ancestry, and
   real-world gravity remain open.

## Custody disposition

No remaining material algebraic, coefficient, statistics, field-census,
convention, anomaly, ancestry, or scope defect was found.  This audit
authorizes source freeze of exactly the theorem and verifier hashes above as
an exact conditional visible-sector one-loop proper-time sign theorem.  It
does not authorize promotion to a complete spectrum result, an absolute
prediction of Newton's constant, or a gravity theorem.
