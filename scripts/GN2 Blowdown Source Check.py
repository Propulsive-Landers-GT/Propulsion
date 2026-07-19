import CoolProp.CoolProp as CP

fluid = "Nitrogen"

# Source GN2 tank
P_source_i_bar = 80.0     # initial source tank pressure [bar absolute]
V_source_L = 50          # source tank volume [L]

# Running tank
P_run_i_bar = 1.01325      # initial running tank pressure [bar absolute]
P_run_f_bar = 65.0         # required constant running tank pressure [bar absolute]
V_run_L = 3.78            # running tank gas volume / ullage volume [L]

# Temperature assumptions
T_source_K = 293.15        # source tank temperature [K]
T_run_K = 293.15           # running tank temperature [K]

# Unit Conversions
bar_to_Pa = 1e5
L_to_m3 = 1e-3

P_source_i = P_source_i_bar * bar_to_Pa
P_run_i = P_run_i_bar * bar_to_Pa
P_run_f = P_run_f_bar * bar_to_Pa

V_source = V_source_L * L_to_m3
V_run = V_run_L * L_to_m3

# Calcs w/ Coolprop
rho_source_i = CP.PropsSI("D", "P", P_source_i, "T", T_source_K, fluid)
rho_run_i = CP.PropsSI("D", "P", P_run_i, "T", T_run_K, fluid)
rho_run_f = CP.PropsSI("D", "P", P_run_f, "T", T_run_K, fluid)

m_source_i = rho_source_i * V_source
m_run_i = rho_run_i * V_run
m_run_f = rho_run_f * V_run

m_added_to_run = m_run_f - m_run_i
m_source_f = m_source_i - m_added_to_run

rho_source_f = m_source_f / V_source

# Final source pressure from CoolProp
P_source_f = CP.PropsSI("P", "D", rho_source_f, "T", T_source_K, fluid)

# Results
print("GN2 Pressurization Result")
print(f"Initial source pressure: {P_source_i_bar:.3f} bar abs")
print(f"Final source pressure:   {P_source_f / bar_to_Pa:.3f} bar abs")
print()
print(f"Initial source mass:     {m_source_i:.4f} kg")
print(f"Final source mass:       {m_source_f:.4f} kg")
print(f"GN2 used:                {m_added_to_run:.4f} kg")
print()
print(f"Initial running mass:    {m_run_i:.4f} kg")
print(f"Final running mass:      {m_run_f:.4f} kg")

if m_source_f <= 0:
    print()
    print("WARNING: Source tank does not contain enough GN2.")
elif P_source_f < P_run_f:
    print()
    print("WARNING: Final source tank pressure is below required running tank pressure.")
    print("The source tank cannot maintain the running tank at the target pressure.")
else:
    print()
    print("Source tank has enough pressure remaining.")