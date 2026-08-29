# GRA-GJ: flip-recoil supplied-embedding Ward obstruction

This lane performs the bounded `RF3a` test of whether GD's already derived
direct encoded recoil momenta can fill FY/FZ's missing temporal momentum slot.
They cannot under FZ's supplied embedding contraction: the direct density is
configuration diagonal, its commutator has zero diagonal, and the complete FY
spatial-source contraction has a nonzero exact diagonal entry.

The verifier is standard-library only and checks dependency custody, GD's
encoded momenta and factor-edge current, the direct ice projection of a
single-link current, exact matrix typing, `Phi_240`, the FZ longitudinal
witness and complete coefficient, contacts, boundaries, and the theorem's
native-divergence/Feshbach ceilings.

Run:

```text
python3 LANE_GRA_GJ_F3_Q4_FLIP_RECOIL_EMBEDDING_WARD_OBSTRUCTION_V001/verify_flip_recoil_embedding_ward_obstruction.py
```

This is a conditional obstruction to one insufficient completion.  It is not
a no-go for a properly source-completed parent and is not gravity closure.
