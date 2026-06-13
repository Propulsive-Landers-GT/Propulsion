import CoolProp.CoolProp as CP

def calculate_nitrous_fill_mass(tank_volume_L: float = 2.8, ullage_fraction: float = 0.15, temperature_C: float = 26.0) -> dict:

    #Unit conv.
    tank_volume_m3 = tank_volume_L * 1e-3          # L  → m^3
    temperature_K  = temperature_C + 273.15         # C → K

    #Liquid volume after accounting for ullage
    liquid_fraction   = 1.0 - ullage_fraction
    liquid_volume_m3  = liquid_fraction * tank_volume_m3
    liquid_volume_L   = liquid_volume_m3 * 1e3

    #CoolProp: saturated liquid density of N2O at fill temperature
    # Density in kg/m^3
    density_kg_m3 = CP.PropsSI("D", "T", temperature_K, "Q", 0, "N2O")

    #Sat pressure
    sat_pressure_Pa  = CP.PropsSI("P", "T", temperature_K, "Q", 0, "N2O")
    sat_pressure_bar = sat_pressure_Pa / 1e5

    #Mass calc
    fill_mass_kg = density_kg_m3 * liquid_volume_m3
    fill_mass_g  = fill_mass_kg * 1e3

    results = {
        "tank_volume_L":       tank_volume_L,
        "ullage_fraction_%":   ullage_fraction * 100,
        "liquid_fraction_%":   liquid_fraction * 100,
        "liquid_volume_L":     liquid_volume_L,
        "fill_temperature_C":  temperature_C,
        "fill_temperature_K":  temperature_K,
        "sat_pressure_bar":    sat_pressure_bar,
        "liquid_density_kg_m3": density_kg_m3,
        "fill_mass_kg":        fill_mass_kg,
        "fill_mass_g":         fill_mass_g,
    }
    return results

def print_results(r: dict) -> None:
    print("       N2O Tank Fill Mass Calculator")
    print(f"  Tank volume          : {r['tank_volume_L']:.3f} L")
    print(f"  Ullage requirement   : {r['ullage_fraction_%']:.1f} %")
    print(f"  Liquid fill fraction : {r['liquid_fraction_%']:.1f} %")
    print(f"  Liquid volume        : {r['liquid_volume_L']:.4f} L")
    print("-" * 5)
    print(f"  Fill temperature     : {r['fill_temperature_C']:.1f} °C  "
          f"({r['fill_temperature_K']:.2f} K)")
    print(f"  Sat. pressure (N2O)  : {r['sat_pressure_bar']:.3f} bar")
    print(f"  Liquid density       : {r['liquid_density_kg_m3']:.2f} kg/m³")
    print("-" * 5)
    print(f"  FILL MASS         : {r['fill_mass_kg']:.4f} kg"
          f"  ({r['fill_mass_g']:.1f} g)")
    print("-" * 5)

#poop
if __name__ == "__main__":
    TANK_VOLUME_L   = 2.8    # L
    ULLAGE_FRACTION = 0.1   # 15 %
    TEMPERATURE_C   = 29.5   # °C  ------(Ambient temp/day temp)

    results = calculate_nitrous_fill_mass(
        tank_volume_L   = TANK_VOLUME_L,
        ullage_fraction = ULLAGE_FRACTION,
        temperature_C   = TEMPERATURE_C,
    )
    print_results(results)