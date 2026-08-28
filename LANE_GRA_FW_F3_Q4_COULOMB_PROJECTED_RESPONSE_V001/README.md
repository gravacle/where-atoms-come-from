# FW Coulomb projected response

This lane composes the two-piece `FV-WITNESS` source -- direct pair plus
irreducible differentiated H6 ring -- with FO's exact 180-state translation-
closed `H6` sector.  It does not include FV's generated diagonal source
derivatives.  Run:

```bash
python3 LANE_GRA_FW_F3_Q4_COULOMB_PROJECTED_RESPONSE_V001/verify_projected_response.py
```

The calculation keeps operator, commutator, retarded, and first-moment ranks
separate.  Its scope is one homogeneous finite-sector ground-state response.
