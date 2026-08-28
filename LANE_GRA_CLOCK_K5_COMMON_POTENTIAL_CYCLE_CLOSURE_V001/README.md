# Clock K5 common-potential cycle-closure lane

This isolated no-laboratory lane reacquires the four official CSV files from
Zenodo DOI `10.5281/zenodo.8184043` and scores only `Fig4.csv`. It performs an
exact K5 cut/cycle decomposition and a covariance-honest marginal-box
compatibility calculation on the ten processed pairwise inferred clock heights.
The fitted node object is a scalar height vector; interpreting its common scalar
multiple as gravitational potential uses the article's already-assumed GR and
common local $g$ mapping.

Run from the repository root:

```bash
python3 LANE_GRA_CLOCK_K5_COMMON_POTENTIAL_CYCLE_CLOSURE_V001/analyze_clock_k5_cycle_closure.py
python3 LANE_GRA_CLOCK_K5_COMMON_POTENTIAL_CYCLE_CLOSURE_V001/verify_clock_k5_cycle_closure.py
```

The lane tests a necessary processed-output common-node-scalar condition only.
It does not use conventional chi-square, does not assume shared-clock
independence, and does not test GR, gravity, a common metric, record lineage, or
gravity emergence.
