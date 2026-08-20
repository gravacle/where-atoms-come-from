#!/usr/bin/env python3
"""GR4_DERIVE -- the physical numbers behind the clause audit.

Every number printed here is either (a) a standard constant, (b) a cited
material parameter (source named inline), or (c) an Arrhenius/Landauer
combination of the two. Nothing is fitted.
"""
import math

kB = 1.380649e-23          # J/K, exact (SI)
T300 = 300.0
kT = kB * T300             # J
eV = 1.602176634e-19       # J
YR = 3.156e7               # s

print("=== FLOORS ===")
print(f"kT(300K)            = {kT:.3e} J = {kT/eV*1000:.1f} meV")
print(f"Landauer kT ln2     = {kT*math.log(2):.3e} J")

def barrier_for_retention(t_m, nu=1e9, delta=1.0):
    """Delta_F/kT >= ln(nu * t_m / delta)."""
    return math.log(nu * t_m / delta)

print("\n=== THE RETENTION INEQUALITY  Delta_F >= kT ln(nu t_m / delta) ===")
for years, nu, delta, label in [
    (10, 1e9, 1.0, "10 yr, nu=1e9/s (Neel attempt), delta=1"),
    (10, 1e9, 1e-2, "10 yr, delta=1e-2"),
    (10, 1e12, 1e-4, "10 yr, nu=1e12/s, delta=1e-4"),
    (1e9, 1e9, 1.0, "1 Gyr geological"),
]:
    print(f"  {label:45s} -> {barrier_for_retention(years*YR, nu, delta):5.1f} kT")
print("  (industry KuV/kT >= 60 'thermal stability rule' sits inside this band)")

print("\n=== RECORD 1: HAMR HDD grain, L1_0 FePt ===")
Ku = 6.6e6        # J/m^3, L1_0 FePt bulk anisotropy (Weller et al., IEEE Trans. Magn. 36 (2000) 10)
d = 5e-9          # m, HAMR-era grain size
V = d**3
Ms = 1.14e6       # A/m, FePt saturation magnetization (same source)
KuV = Ku * V
print(f"  KuV = {KuV:.2e} J = {KuV/kT:.0f} kT  (barrier)")
tau0 = 1e-9       # s, Neel attempt time (1e-9..1e-11 conventional)
tau = tau0 * math.exp(min(KuV/kT, 700))
print(f"  Neel-Arrhenius tau = tau0 exp(KuV/kT) ~ 1e-9 * e^{KuV/kT:.0f} s  (astronomically > age of universe)")
m_grain = Ms * V
dE_earth = 2 * m_grain * 5e-5   # Zeeman split in Earth's field 50 uT
print(f"  clause-(iii) check: up/down splitting in Earth's field = 2 m B = {dE_earth:.2e} J = {dE_earth/kT:.1e} kT")
print(f"    -> NOT exactly degenerate; harmless because {dE_earth/KuV:.1e} of the barrier")
# write
H_head = 2.0      # T, max realistic head field; FePt room-T coercivity ~ 4-6 T -> needs heat assist
E_write_min = 2 * m_grain * H_head
print(f"  write: Zeeman work scale 2 m B_head = {E_write_min:.2e} J = {E_write_min/kT:.0f} kT; plus laser heat pulse to ~700 K")
drive_power, drive_rate = 6.0, 2.0e9   # W, bits/s (order: 250 MB/s drive)
print(f"  device-level write energy ~ {drive_power/drive_rate:.1e} J/bit = {drive_power/drive_rate/(kT*math.log(2)):.1e} x Landauer")
print(f"  timescale ratio tau/t_write ~ 1e17s / 1e-9s = 1e26")

print("\n=== RECORD 2: DNA base (covalent record) ===")
k_dep = 3e-11     # /s per purine, depurination at 37C pH7.4 (Lindahl & Nyberg, Biochemistry 11 (1972) 3610)
Ea = 130e3        # J/mol (~31 kcal/mol, same literature)
T310 = 310.0
print(f"  depurination k = 3e-11 /s -> tau = {1/k_dep/YR:.0f} yr per site (unrepaired)")
print(f"  Ea = 130 kJ/mol = {Ea/(8.314*T310):.0f} kT_310  (Arrhenius barrier)")
print(f"  fossil-DNA half-life 521 yr at 13C (Allentoft et al., Proc R Soc B 279 (2012) 4724)")
atp = 30.5e3 / 6.022e23    # J per ATP hydrolysis (standard biochem, ~30.5 kJ/mol)
print(f"  write: polymerase, ~2 phosphodiester-bond equivalents ~ {2*atp:.1e} J = {2*atp/kT:.0f} kT/base")
print(f"    ({2*atp/(kT*math.log(2)):.0f} x Landauer -- biology is within ~2 orders of the floor)")

print("\n=== RECORD 3: NAND flash floating gate ===")
phi = 3.1 * eV    # Si/SiO2 electron barrier
print(f"  barrier = 3.1 eV = {phi/kT:.0f} kT; thermal escape alone would give tau ~ 1e-13*e^120 s ~ 1e39 s")
print(f"  actual retention spec 10 yr (JEDEC, 55C): limited by oxide-defect tunneling (SILC), NOT thermal --")
print(f"  a reminder that real durability is set by the worst channel, not the named one")
n_e, Vpp = 1e4, 20.0
print(f"  write: ~{n_e:.0e} electrons through ~{Vpp:.0f} V ~ {n_e*eV/1.0*Vpp:.1e} J = {n_e*eV*Vpp/(kT*math.log(2)):.1e} x Landauer")

print("\n=== RECORD 4: K-Ar geochronometer (feldspar/mica) ===")
print(f"  bit(s): trapped 40Ar amount; written by 40K decay, t_half = 1.25 Gyr (weak interaction)")
Ea_Ar = 180e3     # J/mol, order of Ar diffusion activation energy in K-feldspar (40-50 kcal/mol range)
print(f"  durability: Ar diffusion Ea ~ 180 kJ/mol = {Ea_Ar/(8.314*T300):.0f} kT_300; closure T ~ 150-350 C")
print(f"  read: mass spectrometry -- DESTRUCTIVE readout; QND reading is not universal among real records")

print("\n=== RECORD 5: baryon number / B-L (the record that atoms exist) ===")
E_sph = 9.0e12 * eV   # sphaleron barrier ~9 TeV (Klinkhamer-Manton 1984; lattice ~9.1 TeV)
print(f"  B: violated only by sphalerons, barrier {E_sph/eV/1e12:.0f} TeV = {E_sph/kT:.1e} kT_300 -> rate ~ e^-3e14")
print(f"  B-L: anomaly-free in the SM -> [H_SM, B-L] = 0 EXACTLY -- the one known real record")
print(f"       satisfying clause (ii) exactly; and it is NOT WRITABLE by any SM process.")
print(f"       Exact durability and writability do not coexist under one fixed dynamics.")

print("\n=== RECORD 6: photograph (AgBr latent image) ===")
print(f"  bit: a >=4-atom Ag cluster on a grain; development amplifies x ~1e9 (standard photographic chemistry)")
print(f"  write: ~1 absorbed photon (2-3 eV = ~100 kT) per cluster atom; redundancy manufactured chemically")

print("\n=== THE TWO-SIDED WORK INEQUALITY (the real (iv)+(v)) ===")
print("  bath writes spontaneously at rate ~ nu exp(-DeltaF/kT); an agent writes by paying W >= DeltaF")
print("  PROTECTED: nu t_m exp(-DeltaF/kT) <= delta  <=>  DeltaF >= kT ln(nu t_m/delta)")
print("  WRITABLE:  a driven protocol supplies DeltaF transiently and dissipates >= kT ln 2 (Landauer)")
print("  -> durability and writability are reconciled by WORK and NONSTATIONARITY, not by symmetry.")
print("  Crooks: P(bath does a work-W write)/P(driven reverse) = exp(-W/kT) -- the protection IS the")
print("  fluctuation theorem, an improbability, never an impossibility.")

print("\n=== REDUNDANCY (the missing clause) ===")
print("  illuminated 1-um dust grain: pointer info redundancy ~1e8 independent photon fragments")
print("  after ~1 us of sunlight (Riedel & Zurek, PRL 105 (2010) 020404; order-of-magnitude)")
flux = 1e21  # photons/m^2/s, bright sunlight, order
area = 1e-6  # m^2, a 1 mm x 1 mm ink mark
print(f"  1 mm^2 ink mark in sunlight scatters ~{flux*area:.0e} photons/s -- copies/s of the bit")
