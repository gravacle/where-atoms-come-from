# GL6AO — Complete Sixth-Order Locked Hamiltonian

This isolated packet computes the full native locked-sector effective
Hamiltonian through order `h^6/U_d^5` on GL6AN's sealed period-four quotient.

The result is exact:

```text
H_eff = scalar
  -(63/8)(h^6/U_d^5) sum_c T_c
  + O(h^8/U_d^7),
```

where `T_c` toggles one alternating six-cycle.  The scalar coefficients per
link are

```text
order 2: -1/2,
order 4: -7/24,
order 6: -893/1080.
```

The alternating six-cycle is the only configuration-changing term at order
six.  The complete diagonal order-six term is common to every locked
configuration; there is no flippable-cycle-count potential at this order.

The packet also defines the corresponding finite-range, bounded formal
linked interaction on the infinite incidence.  It does not claim an
all-orders limit, phase, pole, physical cone, photon, graviton, gravity, or
`G`.

Run:

```bash
python3 LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/verify_sixth_order_hamiltonian.py
python3 LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/verify_packet.py
```

Expected results are recorded in `VERIFICATION.txt`.
