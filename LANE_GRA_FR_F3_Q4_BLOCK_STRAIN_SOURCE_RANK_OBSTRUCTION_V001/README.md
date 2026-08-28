# GRA-FR: q4 block-strain source rank obstruction

This isolated lane performs the exact early eligibility test for the
additive edge-supported component of the `Q4-BLOCK-STRAIN-CTP` successor
frozen by `FQ`.

Run:

```bash
python3 LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/verify_block_strain_source_rank_obstruction.py
python3 LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/independent_hostile_audit.py
cd LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001
shasum -a 256 -c MANIFEST.sha256
shasum -a 256 -c SEAL.sha256
```

The verifier uses standard-library rational arithmetic, replays every pinned
dependency hash, checks the tetrahedral dyad Gram matrix and `S4` character,
proves the two `E` nulls, exhausts all additive label multiplicities through
eight occurrences, compares additive pair weights with the rank-six `A3`
root dyads, checks an exact Feshbach Schur complement, separates induced from
general `O(j^2)` seagulls, and verifies the source-off linear
effective-operator and commutator-moment nulls.

The result closes only the additive edge-supported BS20/FQ17a source
subclass. It does not rank the still-unfrozen onsite/node/port/boundary/
controller linear weights required for the complete source. It neither
changes the source-off F3/q4 Hamiltonian nor rules out a separately derived
cross-dyad, blocked-root, loop/surface, gluing, or thermodynamic collective
source. The immediate successor is `Q4-COMPLETE-SOURCE-RANK-AUDIT`.
