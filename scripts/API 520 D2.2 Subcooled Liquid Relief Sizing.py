# API 520 (Annex D.2.2) — Subcooled liquid at inlet (two-phase flashing) sizing
# MODIFIED for:
#  1) PRV only
#  2) Burst disk ONLY
#  3) Burst disk + PRV combination (via Kc factor)

# IMPORTANT:
# - Disk-only sizing (API 520 3.11.1.2) uses the SAME fluid equations, but with Kd = 0.62.
# - Disk+PRV is NOT "0.62 * 0.9".
#   For disk+PRV you size the PRV and apply Kc (often 0.90 if no published value), i.e.
#     capacity(combination) = Kc * capacity(PRV)
#   In the area equation this appears as dividing by Kc (or multiplying required area by 1/Kc).
#   Kd remains the PRV discharge coefficient (e.g., 0.65 preliminary for subcooled liquids).

import math
from CoolProp.CoolProp import PropsSI

# Inputs
fluid = "NitrousOxide"

T_C = 20.0                     # inlet liquid temperature [°C]
set_pressure_bar_g = 80.0      # set pressure [bar(g)]
overpressure_frac = 0.10       # allowable overpressure fraction (10% typical)
Pa_bar_a = 1.01325             # downstream total backpressure [bar(a)]

mdot_required = 0.56           # required relieving mass flow [kg/s]

# Choose device configuration:
#   "prv_only"      -> PRV only (Kd ~ 0.65 for subcooled liquid preliminary)
#   "disk_only"     -> Burst disk only (Kd = 0.62 per API 520 3.11.1.2 coefficient method)
#   "disk_plus_prv" -> Burst disk upstream of PRV (apply Kc ~ 0.90 if no published value)
device_mode = "disk_plus_prv"

# Coefficients / factors
Kd_prv_subcooled = 0.65        # API preliminary discharge coefficient for PRV with subcooled liquid (D.2.2)
Kd_disk = 0.62                 # API 520 3.11.1.2 coefficient method for rupture disk device
Kb = 1.0                       # backpressure correction (liquid); usually 1.0 if conventional & low built-up BP
Kc = 0.90                      # combination factor for disk+PRV if no published value (API 520 3.11.2)
# If you have published Kc for your specific disk/valve combo, replace 0.90.

# Unit conversions
T_K = T_C + 273.15

bar_to_Pa = 1e5
psi_to_Pa = 6894.757293168
Pa_to_psia = 1.0 / psi_to_Pa

kgm3_to_lbft3 = 0.0624279606
Jkg_to_Btu_lb = 0.000429922614
JkgK_to_Btu_lbR = 0.000238845897
m3kg_to_ft3lb = 16.01846337396

# Pressures (API definitions)
# Po = relieving pressure at inlet: set pressure plus allowed gauge overpressure, then convert to absolute
Pset_bar_a = set_pressure_bar_g + 1.01325
Po_bar_a = set_pressure_bar_g * (1.0 + overpressure_frac) + 1.01325
Po_Pa = Po_bar_a * bar_to_Pa
Pa_Pa = Pa_bar_a * bar_to_Pa

# Saturation pressure at To
Ps_Pa = PropsSI("P", "T", T_K, "Q", 0, fluid)

# CoolProp properties for omega_s (Eq D.8)
rho_lo = PropsSI("D", "T", T_K, "P", Po_Pa, fluid)           # kg/m^3
cp = PropsSI("Cpmass", "T", T_K, "P", Po_Pa, fluid)          # J/kg-K

rho_l_sat = PropsSI("D", "T", T_K, "Q", 0, fluid)
rho_v_sat = PropsSI("D", "T", T_K, "Q", 1, fluid)
v_l = 1.0 / rho_l_sat
v_v = 1.0 / rho_v_sat

h_l = PropsSI("Hmass", "T", T_K, "Q", 0, fluid)
h_v = PropsSI("Hmass", "T", T_K, "Q", 1, fluid)
hvls = h_v - h_l

# Eq D.8 in US units
# omega_s = [0.185 * rho_lo * Cp * To * Ps * (vvls/hvls)]^2
rho_lo_lbft3 = rho_lo * kgm3_to_lbft3
cp_Btu_lbR = cp * JkgK_to_Btu_lbR
To_R = (T_C * 9.0/5.0) + 491.67
Ps_psia = Ps_Pa * Pa_to_psia

vvls_ft3lb = (v_v - v_l) * m3kg_to_ft3lb
hvls_Btu_lb = hvls * Jkg_to_Btu_lb

omega_s = (0.185 * rho_lo_lbft3 * cp_Btu_lbR * To_R * Ps_psia * (vvls_ft3lb / hvls_Btu_lb)) ** 2

# Step 2: subcooling region
eta_st = (2.0 * omega_s) / (1.0 + 2.0 * omega_s)
Po_psia = Po_Pa * Pa_to_psia
region = "low" if (Ps_psia > eta_st * Po_psia) else "high"

# Step 3: critical vs subcritical (HIGH subcooling region)
Pa_psia = Pa_Pa * Pa_to_psia
if region == "high":
    criticality = "critical" if (Ps_Pa > Pa_Pa) else "subcritical_all_liquid"
else:
    criticality = "needs_Fig_D3"

# Step 4: mass flux (use SI equivalent of Eq D.11)
# G = sqrt(2 * rho_lo * (Po - P))  [kg/s/m^2]
# For critical: P = Ps ; for subcritical(all-liquid): P = Pa
if region != "high":
    raise SystemExit("Low-subcooling region detected. Implement Figure D.3 (ηc) path or re-check inputs.")

P_use = Ps_Pa if criticality == "critical" else Pa_Pa
dP = Po_Pa - P_use
if dP <= 0:
    raise SystemExit("Po <= P_use (Ps or Pa). Check pressures/temperature/scenario.")

G = math.sqrt(2.0 * rho_lo * dP)

# Step 5: required effective discharge area
# Baseline (API D.12): A = mdot / (Kd * Kb * Kc * G)
# For disk-only: Kd = 0.62 (API 520 3.11.1.2)
# For PRV-only:  Kd = 0.65 (prelim for subcooled liquid)
# For disk+PRV:  Kd = PRV's (0.65 prelim) AND apply Kc (0.90 default if no published value)

if device_mode == "prv_only":
    Kd_used = Kd_prv_subcooled
    Kc_used = 1.0
elif device_mode == "disk_only":
    Kd_used = Kd_disk
    Kc_used = 1.0
elif device_mode == "disk_plus_prv":
    Kd_used = Kd_prv_subcooled
    Kc_used = Kc
else:
    raise SystemExit("device_mode must be: 'prv_only', 'disk_only', or 'disk_plus_prv'.")

A_m2 = mdot_required / (Kd_used * Kb * Kc_used * G)
A_in2 = A_m2 * 1550.0031000062

# API 526 orifice areas (in^2) — for PRV selection only.
api526 = {
    "D": 0.110, "E": 0.196, "F": 0.307, "G": 0.503, "H": 0.785,
    "J": 1.287, "K": 1.838, "L": 2.853, "M": 3.600, "N": 4.340,
    "P": 6.380, "Q": 11.050, "R": 16.000, "T": 26.000
}
pick = None
for letter, area in sorted(api526.items(), key=lambda x: x[1]):
    if area >= A_in2:
        pick = letter
        break

# Results
print("=== API 520 D.2.2 Subcooled Liquid Relief Sizing (CoolProp) ===")
print(f"Fluid: {fluid}")
print(f"Mode: {device_mode}")
print()
print(f"T_inlet = {T_C:.2f} °C  ({To_R:.2f} R)")
print(f"Pset = {set_pressure_bar_g:.2f} bar(g)  | Overpressure = {overpressure_frac*100:.1f}%")
print(f"Po (relieving inlet) = {Po_bar_a:.3f} bar(a)  = {Po_psia:.2f} psia")
print(f"Pa (backpressure)    = {Pa_bar_a:.3f} bar(a)  = {Pa_psia:.2f} psia")
print(f"Psat(To)             = {(Ps_Pa/bar_to_Pa):.3f} bar(a)  = {Ps_psia:.2f} psia")
print()
print(f"rho_lo(Po,To) = {rho_lo:.2f} kg/m^3  ({rho_lo_lbft3:.2f} lb/ft^3)")
print(f"Cp(Po,To)     = {cp:.1f} J/kg-K  ({cp_Btu_lbR:.4f} Btu/lb-R)")
print(f"hvls(To)      = {hvls/1000:.2f} kJ/kg  ({hvls_Btu_lb:.2f} Btu/lb)")
print(f"vvls(To)      = {(v_v-v_l):.6e} m^3/kg  ({vvls_ft3lb:.6e} ft^3/lb)")
print()
print(f"omega_s = {omega_s:.4f}")
print(f"eta_st  = {eta_st:.4f}")
print(f"Region  = {region} subcooling")
print(f"Flow    = {criticality}")
print(f"Mass flux G = {G:.2f} kg/s/m^2")
print()
print(f"mdot_required = {mdot_required:.6f} kg/s")
print(f"Kd_used = {Kd_used:.3f} | Kb = {Kb:.3f} | Kc_used = {Kc_used:.3f}")
print(f"A_required = {A_m2:.6e} m^2  = {A_in2:.6f} in^2")
print()

if device_mode in ("prv_only", "disk_plus_prv"):
    print(f"API 526 PRV orifice pick: {pick if pick else 'Larger than T / custom'}")
else:
    print("Disk-only: select a rupture disk device whose MINIMUM NET FLOW AREA >= A_required.")
    print("           (Use manufacturer 'net flow area', not API 526 orifice letters.)")