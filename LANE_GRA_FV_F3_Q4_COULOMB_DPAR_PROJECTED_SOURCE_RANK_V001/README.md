# FV projected-source rank lane

This isolated lane performs the next physics calculation after FU.  Run:

```bash
python3 verify_projected_source_rank.py
```

The executable uses exact rational arithmetic, rebuilds `G_5`, enumerates all
720 H6 paths, differentiates every numerator and denominator, constructs
global ice witnesses, and checks operator rank.  No canonical file or git
state is modified by this lane.

The theorem is conditional on FU's `S1`--`S9` physical solder and the stronger
`S10 / FV-PURE` premise identifying the complete nonidentity first source
derivative with the ideal-Coulomb pair source plus the unchanged FS one-edge
source.  FU alone permits residual kernels and does not imply this premise.
Under that explicit completion the lane closes projected off-shell source rank
through leading H6 and formally through H8; it deliberately does not promote
that rank to a CTP tensor pole or gravity response.
