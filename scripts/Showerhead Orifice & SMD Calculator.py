import math
from CoolProp.CoolProp import PropsSI

# Inputs
m_ox_dot_target = 0.56      # Target oxidizer mass flow rate [kg/s]
Cd = 0.31                   # Discharge coefficient [-]
N_inj = 114                 # Number of injector orifices [-]
Delta_P = 15e5              # Pressure drop across injector [Pa]

T_fluid_K = 303.15          # Fluid temperature [K]
fluid = "NitrousOxide"      # CoolProp fluid name
property_mode = "coolprop"  # "coolprop" or "paper_table_3"

# Choose how to get density for sizing
use_saturated_liquid_density_for_sizing = True
P_bulk_Pa = 6.0e6           # Only used if above is False

# Viscosity for SMD correlation. CoolProp does not provide N2O viscosity here.
mu_L = 325e-6               # [Pa*s]

# Fluid Properties
if property_mode == "coolprop":
    rho_ox_L = PropsSI("D", "T", T_fluid_K, "Q", 0, fluid)   # Saturated liquid density [kg/m^3]
    rho_ox_G = PropsSI("D", "T", T_fluid_K, "Q", 1, fluid)   # Saturated vapor density [kg/m^3]
    sigma = PropsSI("I", "T", T_fluid_K, "Q", 0, fluid)      # Surface tension [N/m]
    smd_label = "CoolProp-state estimate"
elif property_mode == "paper_table_3":
    rho_ox_L = 772.25      # [kg/m^3]
    rho_ox_G = 1.83        # [kg/m^3], gas density at ambient pressure and temperature per paper
    sigma = 0.24           # [N/m]
    mu_L = 325e-6          # [Pa*s]
    smd_label = "paper Table 3 estimate"
else:
    raise ValueError('property_mode must be "coolprop" or "paper_table_3"')

if use_saturated_liquid_density_for_sizing:
    rho_ox = rho_ox_L
else:
    rho_ox = PropsSI("D", "T", T_fluid_K, "P", P_bulk_Pa, fluid)

# Required orifice/diameter for target flow
# m_dot = Cd * A_total * sqrt(2 * rho * Delta_P)
A_inj_total_required = m_ox_dot_target / (Cd * math.sqrt(2.0 * rho_ox * Delta_P))   # [m^2]
A_single_inj_required = A_inj_total_required / N_inj                                  # [m^2]
d_inj_required = math.sqrt(4.0 * A_single_inj_required / math.pi)                     # [m]

# Velocity through injector at target condition
u_ox = m_ox_dot_target / (rho_ox * A_inj_total_required)                              # [m/s]

# SMD estimate:
SMD = (
    47.0
    * (d_inj_required / u_ox)
    * ((sigma / rho_ox_G) ** 0.25)
    * (1.0 + 331.0 * (mu_L / math.sqrt(rho_ox_L * sigma * d_inj_required)))
)

SMD_micrometers = SMD * 1e6

# Outputs
print("INPUTS:")
print(f"Fluid: {fluid}")
print(f"Fluid temperature: {T_fluid_K:.2f} K")
print(f"Property mode: {property_mode}")
print(f"Target oxidizer mass flow rate: {m_ox_dot_target:.6f} kg/s")
print(f"Discharge coefficient (Cd): {Cd:.6f}")
print(f"Number of orifices: {N_inj}")
print(f"Injector pressure drop: {Delta_P:.3e} Pa ({Delta_P/1e5:.3f} bar)")

if use_saturated_liquid_density_for_sizing:
    print("Sizing density mode: saturated liquid density at T")
else:
    print(f"Sizing density mode: density at T and P = {P_bulk_Pa:.3e} Pa")

print("\nFLUID PROPERTIES USED:")
print(f"rho_ox (sizing density): {rho_ox:.5f} kg/m^3")
print(f"rho_ox_L (liquid density): {rho_ox_L:.5f} kg/m^3")
print(f"rho_ox_G (gas density): {rho_ox_G:.5f} kg/m^3")
print(f"mu_L (liquid viscosity): {mu_L:.8e} Pa*s")
print(f"sigma (surface tension): {sigma:.8f} N/m")

print("\nREQUIRED INJECTOR GEOMETRY FOR TARGET MDOT:")
print(f"Required total injector area: {A_inj_total_required:.8e} m^2")
print(f"Required single-orifice area: {A_single_inj_required:.8e} m^2")
print(f"Required single-orifice diameter: {d_inj_required * 1e3:.5f} mm")
print(f"Oxidizer axial velocity: {u_ox:.5f} m/s")
print(f"Sauter Mean Diameter (SMD, {smd_label}): {SMD_micrometers:.2f} micrometers")
