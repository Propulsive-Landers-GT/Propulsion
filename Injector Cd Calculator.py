from CoolProp.CoolProp import PropsSI
import math

# Oxidizer mass flow rate [kg/s]
m_dot = 0.42

# Pressure drop across injector [bar]
delta_p_bar = 3.2

# Total injector flow area [m^2]
flow_area = 8.43859931e-05

# Nitrous temperature [C]
T_C = 24.0

# FUNCTIONS
def get_saturated_nitrous_density(T_C):
    """
    Returns saturated liquid nitrous density [kg/m^3]
    at the given temperature in Celsius.
    """
    fluid = "NitrousOxide"
    T_K = T_C + 273.15
    rho = PropsSI("D", "T", T_K, "Q", 0, fluid)
    return rho


# Calculates discharge coefficient Cd from: m_dot = Cd * A * sqrt(2 * rho * delta_p)
def calculate_cd(m_dot, delta_p, flow_area, rho):
    if delta_p <= 0:
        raise ValueError("Pressure drop must be greater than 0.")
    if flow_area <= 0:
        raise ValueError("Flow area must be greater than 0.")
    if rho <= 0:
        raise ValueError("Density must be greater than 0.")
    if m_dot <= 0:
        raise ValueError("Mass flow rate must be greater than 0.")

    cd = m_dot / (flow_area * math.sqrt(2.0 * rho * delta_p))
    return cd

# MAIN
def main():
    print("Injector Cd Calculator")

    delta_p = delta_p_bar * 1e5  # bar -> Pa
    rho = get_saturated_nitrous_density(T_C)
    cd = calculate_cd(m_dot, delta_p, flow_area, rho)

    print("\nResults:")
    print(f"Mass flow rate               : {m_dot:.6f} kg/s")
    print(f"Pressure drop               : {delta_p:.3f} Pa ({delta_p_bar:.3f} bar)")
    print(f"Nitrous temperature         : {T_C:.3f} C")
    print(f"Saturated liquid density    : {rho:.3f} kg/m^3")
    print(f"Total injector flow area    : {flow_area:.9e} m^2")
    print(f"Calculated discharge coeff. : {cd:.6f}")


if __name__ == "__main__":
    main()