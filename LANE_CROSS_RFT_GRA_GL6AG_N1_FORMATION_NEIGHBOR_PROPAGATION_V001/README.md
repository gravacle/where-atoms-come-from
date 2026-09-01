# GL6AG — exact matched formation-to-neighbor propagation

This mutable author packet tests one narrow finite statement in the sealed
GL6AA/GL6AB `N=1` sixteen-link parent.  Cells `1,2,3` remain formed/KEEP.
Only cell `0`'s physical `K`-support word varies, and every neighbor mean is
subtracted against the matched `0000` branch.

Fast structural replay:

```sh
python3 verify_structure_and_ledger.py
```

Independent full-basis replay:

```sh
c++ -O3 -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib \
  verify_n1_matched_formation_propagation.cpp -lgmpxx -lgmp \
  -o /tmp/verify_gl6ag_exact
/tmp/verify_gl6ag_exact
```

The C++ replay constructs every branch Hamiltonian directly on all `65536`
active computational words.  It checks the baseline, four singles, all six
pairs, the complete sixteen-pattern order-twelve receiver census, pair
Möbius terms through order sixteen, and the bridge-off term-ablation
diagnostic.  It imports
no response coefficient from an earlier packet.

Status: author frozen and sealed after clean independent hostile pre-freeze
review.  A distinct independent post-freeze custody audit is required.
