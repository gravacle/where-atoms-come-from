# LANE_GR2_PROTECTION — arithmetic check for the protection census.
# Every number below is standard-literature physics recomputed from constants;
# provenance for each input is stated in GR2_PROTECTION_CENSUS.md.
import math

kB = 1.380649e-23      # J/K
eV = 1.602176634e-19   # J
T300 = 300.0
kT300 = kB * T300      # J
print(f"kT at 300 K = {kT300:.3e} J = {kT300/eV:.4f} eV")
print(f"Landauer kT ln2 at 300 K = {kT300*math.log(2):.3e} J")

def arr(tau0, Delta):
    """Neel-Arrhenius lifetime tau = tau0 * exp(Delta)."""
    return tau0 * math.exp(Delta)

yr = 3.156e7  # s

print("\n--- 1. HDD grain (CoCrPt-oxide perpendicular media) ---")
# grain: cylinder, diameter 8 nm, height 8 nm  (Weller-Moser era design; modern PMR)
r, h = 4e-9, 8e-9
V = math.pi * r*r * h
Ku = 3e5  # J/m^3, CoCrPt alloy media (range 2-6e5)
Eb = Ku * V
print(f"V = {V:.2e} m^3, Ku = {Ku:.1e} J/m^3, Eb = {Eb:.2e} J = {Eb/eV:.2f} eV, Delta = {Eb/kT300:.1f} kT")
for D in (40, 60):
    t = arr(1e-9, D)
    print(f"  Delta={D}: tau = 1e-9*e^{D} = {t:.2e} s = {t/yr:.2e} yr")
# write: head field ~0.8-1.0 T over ~50nm; Zeeman energy to reverse grain moment
Ms = 5e5  # A/m CoCrPt
mu = Ms * V  # A m^2
B = 1.0
print(f"grain moment mu = {mu:.2e} A m^2; Zeeman 2*mu*B at 1 T = {2*mu*B:.2e} J = {2*mu*B/eV:.1f} eV")
print(f"HDD write dissipation ~10 pJ/bit => {1e-11/(kT300*math.log(2)):.1e} x Landauer")

print("\n--- 2. NAND flash floating gate ---")
print(f"Si/SiO2 electron barrier = 3.1 eV = {3.1*eV/kT300:.0f} kT300")
print(f"retention Ea ~ 1.1 eV = {1.1*eV/kT300:.1f} kT300; tau(300K) = {arr(1e-13, 1.1*eV/kT300):.2e} s = {arr(1e-13,1.1*eV/kT300)/yr:.1e} yr")
# NB attempt prefactor for charge-loss is uncertain; 1e-13 s is a phonon-time guess. Bound below.
print(f"  same Ea with tau0=1e-9: {arr(1e-9,1.1*eV/kT300)/yr:.1e} yr")

print("\n--- 3. DNA base (cytosine deamination) ---")
Ea = 1.2 * eV  # ~28 kcal/mol
print(f"Ea ~ 1.2 eV = {Ea/kT300:.1f} kT300")
t_half_ss = 200*yr   # Lindahl: ssDNA ~200 yr at 37C pH 7.4
print(f"ssDNA t1/2 ~ 200 yr (Lindahl 1993); dsDNA ~140x slower ~ {200*140:.0f} yr per base")
print(f"UV photon 254 nm = {6.626e-34*3e8/254e-9/eV:.2f} eV")

print("\n--- 5. CD-RW GST mark ---")
print(f"GST crystallization Ea ~ 2.3 eV = {2.3*eV/kT300:.0f} kT300; tau(300K,tau0=1e-13) = {arr(1e-13, 2.3*eV/kT300)/yr:.1e} yr")
print(f"write pulse 10 mW x 50 ns = {10e-3*50e-9:.1e} J")

print("\n--- 6. STT-MRAM MTJ ---")
print(f"Delta=60 => Eb = {60*kT300/eV:.2f} eV; write 100 uA x 0.5 V x 10 ns = {100e-6*0.5*10e-9:.1e} J = {100e-6*0.5*10e-9/(kT300*math.log(2)):.1e} x Landauer")

print("\n--- 7. Superconducting flux quantum (Al nanowire ring) ---")
# phase-slip barrier ~ condensation energy in a coherence-length segment: (Bc^2/2mu0)*A*xi
mu0 = 4e-7*math.pi
Bc = 0.0105  # T, Al thermodynamic critical field
A = 100e-9*100e-9  # wire cross-section 100x100 nm
xi = 100e-9        # Al coherence length ~100 nm (dirty-limit shorter; take 100nm clean-ish)
dF = (Bc**2/(2*mu0))*A*xi
print(f"condensation energy density = {Bc**2/(2*mu0):.1f} J/m^3; dF = {dF:.2e} J = {dF/eV:.2f} eV")
for T in (1.0, 4.2):
    print(f"  Delta at {T} K = {dF/(kB*T):.0f}  => tau astronomically large (e^Delta)")
print("File & Mills 1963: persistent-current decay bound, tau > ~1e5 yr (measured NMR bound)")

print("\n--- 8. sphaleron/baryon number ---")
alphaW = 1/30
print(f"T=0 B-violation suppression ~ exp(-4pi/alpha_W) = exp(-{4*math.pi/alphaW:.0f}) = 10^{-4*math.pi/alphaW/math.log(10):.0f}")

print("\n--- 9. lunar crater ---")
# 1 km crater: impactor ~50 m diameter rock at 17 km/s (pi-scaling, Melosh order-of-magnitude)
rho, d_imp, v = 3000.0, 50.0, 17e3
m = rho*(4/3)*math.pi*(d_imp/2)**3
print(f"impactor mass = {m:.2e} kg, KE = {0.5*m*v*v:.2e} J")

print("\n--- toric code at finite T (no self-correction) ---")
# anyon pair creation cost 2*Delta_gap, then free diffusion: lifetime ~ e^{2Delta/kT} INDEPENDENT of L
print("tau_toric ~ tau0 * exp(2*Dgap/kT), no L dependence (Alicki-Fannes-Horodecki 2009; Bravyi-Terhal 2009)")
print("vs HDD grain: Eb = Ku*V grows with VOLUME -> nature's barrier is extensive, code distance is not")
