# GR3_WRITING: configuration-energy and write-cost numbers for real records.
# Every input is labelled with its source class:
#   [S] standard/textbook value   [E] engineering-typical value (order of magnitude)   [B] bound, not a value
import math
kB = 1.380649e-23      # J/K  [S]
T  = 300.0             # K storage temperature
kT = kB*T              # J
eV = 1.602176634e-19   # J
print(f"kT(300K) = {kT:.3e} J = {kT/eV*1000:.1f} meV ; kT ln2 = {kT*math.log(2):.3e} J")

print("\n--- 1. HDD grain (CoCrPt-oxide perpendicular media) ---")
d  = 8e-9              # grain diameter [E, ~7-9 nm modern PMR]
t  = 10e-9             # thickness [E]
V  = math.pi/4*d*d*t
Ku = 3e5               # J/m^3 [E, CoCrPt 2-4e5]
Ms = 6e5               # A/m  [E, CoCrPt ~0.5-0.8 T/mu0]
m  = Ms*V              # A m^2
Eb = Ku*V
print(f"V = {V:.2e} m^3, moment m = {m:.2e} J/T, barrier KuV = {Eb:.2e} J = {Eb/kT:.0f} kT")
for name,B in [("Earth field 50 uT",5e-5),("neighbour-grain stray 20 mT [E]",0.02),("write field 1 T [E head field 1-2.4 T]",1.0)]:
    print(f"  Zeeman splitting 2mB @ {name}: {2*m*B:.2e} J = {2*m*B/kT:.3g} kT")
print(f"  hysteresis write loss ~ 2KuV per grain = {2*Eb/kT:.0f} kT; ~10 grains/bit -> ~{20*Eb/kT:.0f} kT/bit  [E]")

print("\n--- 2. Flash floating gate ---")
N_e = 100              # electrons stored [E; older nodes ~1000, modern 3D NAND ~50-200 per level]
dV  = 3.0              # threshold shift, V [E]
E_store = 0.5*N_e*eV*dV   # charging energy Q^2/2C = (1/2) Q dV
print(f"stored electrostatic energy ~ (1/2)N e dV = {E_store:.2e} J = {E_store/kT:.0f} kT  (splitting between record values)")
print(f"single electron across 3.1 eV SiO2 barrier: {3.1*eV/kT:.0f} kT barrier per carrier [S]")
Vpp=18.0
E_prog = N_e*eV*Vpp    # electrons pumped through ~18 V during program [E]
print(f"program energy, intrinsic charge transport ~ N e Vpp = {E_prog:.2e} J = {E_prog/kT:.0f} kT; device-level ~0.1-100 pJ/bit = {0.1e-12/kT:.1e}-{100e-12/kT:.1e} kT [E]")

print("\n--- 3. Zircon U-Pb ---")
Q = 51.7e6*eV          # total decay energy 238U -> 206Pb chain [S]
print(f"238U->206Pb chain Q = 51.7 MeV = {Q:.2e} J = {Q/kT:.2e} kT per decayed atom")
print(f"alpha-decay Coulomb barrier ~25-30 MeV [S]; t_half = 4.47 Gyr [S]")

print("\n--- 4. Photographic grain (AgBr latent image) ---")
n_ph=4; E_ph=2.5*eV
print(f"minimum latent image ~{n_ph} Ag atoms, ~{n_ph}-10 absorbed photons at 2.5 eV: write energy >= {n_ph*E_ph/kT:.0f} kT [S: Gurney-Mott]")
print(f"stored chemical energy of Ag_n cluster vs lattice: order 1 eV scale = ~40 kT/atom [B: 0.1-1 eV/atom]")

print("\n--- 5. DNA base identity ---")
print(f"sequence isomers differ by ~1-3 kcal/mol duplex stability = {1*4184/6.022e23/kT:.1f}-{3*4184/6.022e23/kT:.1f} kT [S: SantaLucia NN params]")
print(f"depurination/hydrolysis activation ~1.1-1.3 eV = {1.2*eV/kT:.0f} kT [S]; replication write cost ~ 2 ATP+2 PPi ~ {50*4184/6.022e23/kT:.0f} kT/base [S: ~30-50 kT]")

print("\n--- 6. Molecular chirality ---")
pved_J_per_mol = 1e-11 # J/mol scale [S: Quack et al., 1e-11..1e-13 J/mol]
print(f"parity-violating energy difference ~1e-11 J/mol = {pved_J_per_mol/6.022e23/kT:.1e} kT per molecule [B: <=1e-13 kT]")
print(f"racemization barrier (typical stereocentre) 100-200 kJ/mol = {100e3/6.022e23/kT:.0f}-{200e3/6.022e23/kT:.0f} kT [S]")

print("\n--- Balance-defect tanh(Delta/2kT) on the JOINT (system+bath) eigenspace ---")
for name,DkT in [("isolated grain, Earth field",0.007),("grain among neighbours",3.0),("flash",7500.0),("zircon (per atom)",2e9),("chirality",1e-13)]:
    x=DkT/2.0
    b=math.tanh(x) if x<20 else 1.0
    print(f"  {name:32s} Delta = {DkT:.3g} kT -> |Tr(P_E R)|/dim ~ {b:.3g}")

print("\n--- Arrhenius lifetime check: does splitting threaten retention? ---")
f0=1e9  # attempt frequency 1/s [S: 1e9-1e12]
for name,EbkT,DkT in [("HDD grain",45,3.0),("tape particle",1000,5.0)]:
    ex = EbkT - DkT/2
    if ex > 700:
        print(f"  {name}: barrier {EbkT} kT, splitting {DkT} kT -> higher-well lifetime > 1e290 s (astronomical)")
    else:
        tau_hi = math.exp(ex)/f0   # shallower well
        print(f"  {name}: barrier {EbkT} kT, splitting {DkT} kT -> higher-well lifetime {tau_hi:.2e} s ({tau_hi/3.15e7:.1e} yr)")

print("\n--- tape particle Earth-field splitting ---")
Ms=3.5e5; V=0.3e-6*0.05e-6*0.05e-6; m=Ms*V
print(f"gamma-Fe2O3 particle 300x50x50 nm: m={m:.2e} J/T, 2mB(Earth) = {2*m*5e-5/kT:.2f} kT")
