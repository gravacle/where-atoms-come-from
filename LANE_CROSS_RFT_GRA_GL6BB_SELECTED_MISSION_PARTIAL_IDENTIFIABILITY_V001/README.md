# GL6BB — selected-mission partial identifiability

This packet takes the shortest physics step after sealed `GL6BA`.  It audits
the actual F3/FPSS/GL6AN/GL5ZZF/GL6AZ/GL6BA custody for the three inputs of a
direct finite-mission calculation, proves the sharp state-free obstruction,
and derives the strongest robust interval that can be reported without
inventing mission data.

It also performs the exact conditional calculation already licensed by the
prepared-blank branch: the `L=0` four-link collar reduces to a five-state
Dicke Hamiltonian for both admitted `R=2` and `R=5/2` members.  The included
calculator uses exact rational arithmetic plus analytic remainder bounds.

Run:

```text
python3 verify_selected_mission_partial_identifiability.py
python3 -O verify_selected_mission_partial_identifiability.py
python3 calculate_prepared_blank_collar0.py --ratio 2 --sigma <same-clock-value>
```

The calculator requires the caller to supply `sigma`; it has no physical
default.  The packet does not use a conventional graviton, Ricci/Einstein
ansatz, gravity identification, or `G`.
