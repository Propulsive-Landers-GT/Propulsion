from pathlib import Path
from datetime import datetime
import subprocess
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# Inputs
L = 0.2                          # m
m_dot_o = 0.45                    # kg/s

# Set this if you know the actual nozzle throat diameter.
# Leave as None to size the throat area so the average chamber pressure matches p_cham.
THROAT_DIAMETER_MM = None

# Fixed parameters
D_pi = 0.05                       # m
a = 0.0000722
n = 0.67
rho_f = 900.0                     # kg/m^3
p_cham = 40 * 100000.0          # Pa
epsilon = 5.4830
t = np.linspace(0.0, 40.0, 1000)

# Nitrous properties at desired temperature
T_n2o_C = 20.0
T_n2o_K = T_n2o_C + 273.15
rho_ox = PropsSI('D', 'T', T_n2o_K, 'Q', 0, 'NitrousOxide')
p_sat_n2o = PropsSI('P', 'T', T_n2o_K, 'Q', 0, 'NitrousOxide')

# Physics/Numerical limits defined by YOU
L_MIN = 0.1
L_MAX = 0.8

M_DOT_O_MIN = 0.2
M_DOT_O_MAX = 0.8

OF_MIN_ALLOWED = 0.2
OF_MAX_ALLOWED = 16

M_DOT_F_MIN_ALLOWED = 1e-8

# CEA Settings
CEA_EXE = Path(r"C:\Users\chald\Downloads\cea-main\cea-main\build\source\cea.exe")
CEA_DATA_DIR = Path(r"C:\Users\chald\Downloads\cea-main\cea-main\data")
BASE_DIR = Path(r"C:\Users\chald\FeedSystemDesign\pythonProject1\CEA Run")

PRESSURE_BAR = p_cham / 1e5
PE_BAR = 1.01325  # 1 atm in bar

# Start with a useful wide range
CEA_OF_MIN = 0.30
CEA_OF_MAX = 20.00
CEA_OF_STEP = 0.20

G0 = 9.80665
DEBUG = True

# Same efficiency treatment you were using
ETA_CSTAR = 0.90
ETA_CF = 0.90

# Debugging purposes

def debug_print(msg: str) -> None:
    if DEBUG:
        print(msg)

# CEA Parsing helpers
def safe_float(x, default=np.nan):
    try:
        return float(str(x).replace("D", "E"))
    except Exception:
        return default

def compute_cstar_from_isp_cf(isp_s: float, cf: float) -> float:
    if np.isnan(isp_s) or np.isnan(cf) or cf <= 0:
        return np.nan
    return isp_s * G0 / cf

def build_cea_input_text(pc_bar: float, of_ratio: float, pressure_ratio: float) -> str:
    return f"""problem rocket equilibrium
  p(bar)={pc_bar}
  o/f={of_ratio}
  pi/p={pressure_ratio}
reactants
  fuel=Paraffin wt%=100 t(k)=298.15 h,cal/mol=-444600 C 73 H 124
  oxid=N2O wt%=100 t(k)=298.15
output siunits
end
"""

def run_cea_case(case_name: str, pc_bar: float, of_ratio: float, pressure_ratio: float, results_dir: Path) -> str:
    base = results_dir / case_name
    inp_path = base.with_suffix(".inp")
    out_path = base.with_suffix(".out")

    inp_path.write_text(build_cea_input_text(pc_bar, of_ratio, pressure_ratio), encoding="utf-8")

    env = os.environ.copy()
    env["CEA_DATA_DIR"] = str(CEA_DATA_DIR)

    subprocess.run(
        [str(CEA_EXE), str(base)],
        check=True,
        cwd=str(results_dir),
        env=env
    )

    return out_path.read_text(encoding="utf-8", errors="ignore")

def parse_property_row(out_text: str, row_name: str):
    float_pattern = r"[-+]?(?:\d*\.?\d+(?:[EeDd][-+]?\d+)?|NAN)"
    row_name_up = row_name.upper()

    for line in out_text.splitlines():
        line_clean = line.strip().upper().replace("D", "E")
        if line_clean.startswith(row_name_up):
            nums = re.findall(float_pattern, line_clean)
            return [safe_float(x) for x in nums]

    return []

def parse_cea_output(out_text: str) -> dict:
    data = {
        "cf": np.nan,
        "isp_raw": np.nan,
        "ivac_raw": np.nan,
        "tc": np.nan,
        "te": np.nan,
        "gamma_c": np.nan,
        "gamma_e": np.nan,
        "mw_c": np.nan,
        "mw_e": np.nan,
        "mach_e": np.nan,
        "sonic_e": np.nan,
        "cstar": np.nan,
    }

    vals = parse_property_row(out_text, "CF")
    if len(vals) >= 2:
        data["cf"] = vals[-1]

    vals = parse_property_row(out_text, "ISP")
    if len(vals) >= 2:
        data["isp_raw"] = vals[-1]

    vals = parse_property_row(out_text, "IVAC")
    if len(vals) >= 2:
        data["ivac_raw"] = vals[-1]

    vals = parse_property_row(out_text, "T, K")
    if len(vals) >= 1:
        data["tc"] = vals[0]
    if len(vals) >= 2:
        data["te"] = vals[-1]

    vals = parse_property_row(out_text, "GAMM")
    if len(vals) >= 1:
        data["gamma_c"] = vals[0]
    if len(vals) >= 2:
        data["gamma_e"] = vals[-1]

    vals = parse_property_row(out_text, "M, (1/N)")
    if len(vals) >= 1:
        data["mw_c"] = vals[0]
    if len(vals) >= 2:
        data["mw_e"] = vals[-1]

    if np.isnan(data["mw_c"]):
        vals = parse_property_row(out_text, "MOL WT")
        if len(vals) >= 1:
            data["mw_c"] = vals[0]
        if len(vals) >= 2:
            data["mw_e"] = vals[-1]

    vals = parse_property_row(out_text, "MACH")
    if len(vals) >= 1:
        data["mach_e"] = vals[-1]

    vals = parse_property_row(out_text, "SON VEL")
    if len(vals) >= 1:
        data["sonic_e"] = vals[-1]

    isp_s = data["isp_raw"] / G0 if not np.isnan(data["isp_raw"]) else np.nan
    data["cstar"] = compute_cstar_from_isp_cf(isp_s, data["cf"])

    return data

# CEA Profile
def generate_cea_profile(pc_bar: float, pe_bar: float, of_min: float, of_max: float, of_step: float):
    run_name = datetime.now().strftime("cea_profile_%Y-%m-%d_%H-%M-%S")
    results_dir = BASE_DIR / "CEA_Results" / run_name
    results_dir.mkdir(parents=True, exist_ok=True)

    pressure_ratio = pc_bar / pe_bar
    of_values = np.arange(of_min, of_max + 0.5 * of_step, of_step)

    results = {
        "of": [],
        "cf": [],
        "isp_s": [],
        "ivac_s": [],
        "cstar": [],
        "tc": [],
    }

    print(f"Generating CEA profile in: {results_dir}")

    for of_ratio in of_values:
        try:
            case_name = f"cea_of_{of_ratio:.4f}".replace(".", "p")
            out_text = run_cea_case(case_name, pc_bar, of_ratio, pressure_ratio, results_dir)
            parsed = parse_cea_output(out_text)

            isp_s = parsed["isp_raw"] / G0 if not np.isnan(parsed["isp_raw"]) else np.nan
            ivac_s = parsed["ivac_raw"] / G0 if not np.isnan(parsed["ivac_raw"]) else np.nan

            results["of"].append(of_ratio)
            results["cf"].append(parsed["cf"])
            results["isp_s"].append(isp_s)
            results["ivac_s"].append(ivac_s)
            results["cstar"].append(parsed["cstar"])
            results["tc"].append(parsed["tc"])

            debug_print(
                f"CEA O/F={of_ratio:.3f} | "
                f"Cf={parsed['cf']:.6f}, "
                f"Isp={isp_s:.3f} s, "
                f"C*={parsed['cstar']:.3f} m/s, "
                f"Tc={parsed['tc']:.3f} K"
            )

        except Exception as e:
            debug_print(f"CEA FAILED at O/F = {of_ratio:.3f}: {e}")

    cea_profile = {key: np.array(val, dtype=float) for key, val in results.items()}

    valid = (
        np.isfinite(cea_profile["of"]) &
        np.isfinite(cea_profile["cf"]) &
        np.isfinite(cea_profile["isp_s"]) &
        np.isfinite(cea_profile["cstar"]) &
        np.isfinite(cea_profile["tc"])
    )

    for key in cea_profile:
        cea_profile[key] = cea_profile[key][valid]

    if len(cea_profile["of"]) < 4:
        raise RuntimeError("CEA sweep did not produce enough valid points to build interpolation profiles.")

    return cea_profile


def ensure_cea_range(cea_profile, of_ratio_profile):
    current_min = float(np.min(cea_profile["of"]))
    current_max = float(np.max(cea_profile["of"]))
    needed_min = float(np.min(of_ratio_profile))
    needed_max = float(np.max(of_ratio_profile))

    if needed_min >= current_min and needed_max <= current_max:
        return cea_profile, False

    new_min = min(current_min, max(0.2, 0.9 * needed_min))
    new_max = max(current_max, 1.1 * needed_max)

    debug_print(
        f"Expanding CEA O/F range from [{current_min:.3f}, {current_max:.3f}] "
        f"to [{new_min:.3f}, {new_max:.3f}]"
    )

    new_profile = generate_cea_profile(PRESSURE_BAR, PE_BAR, new_min, new_max, CEA_OF_STEP)
    return new_profile, True

# Engine performance calculator
def evaluate_current_design(L_current, m_dot_o_current, cea_profile, apply_efficiency=True, allow_expand=True):
    L_current = float(L_current)
    m_dot_o_current = float(m_dot_o_current)

    if not np.isfinite(L_current) or not np.isfinite(m_dot_o_current):
        return {"valid": False, "reason": "Fuel grain length and oxidizer flow must be finite numbers."}

    if L_current < L_MIN or L_current > L_MAX:
        return {
            "valid": False,
            "reason": f"Fuel grain length is outside the allowed range: L = {L_current:.6f} m"
        }

    if m_dot_o_current < M_DOT_O_MIN or m_dot_o_current > M_DOT_O_MAX:
        return {
            "valid": False,
            "reason": f"Oxidizer mass flow is outside the allowed range: m_dot_o = {m_dot_o_current:.6f} kg/s"
        }

    Dp_t_local = (
        (D_pi ** (2 * n + 1)) +
        ((2 * n + 1) * (2 ** (2 * n + 1)) * a * (m_dot_o_current ** n) * t) / (np.pi ** n)
    ) ** (1.0 / (2 * n + 1))

    Gox_local = (4.0 * m_dot_o_current) / (np.pi * Dp_t_local * Dp_t_local)
    r_dot_t_local = a * (Gox_local ** n)
    r_dot_t_mm_local = r_dot_t_local * 1000.0
    m_dot_f_local = rho_f * np.pi * Dp_t_local * L_current * r_dot_t_local

    if not np.all(np.isfinite(m_dot_f_local)):
        return {"valid": False, "reason": "Non-finite fuel mass flow encountered."}

    if np.min(m_dot_f_local) <= 0 or np.mean(m_dot_f_local) < M_DOT_F_MIN_ALLOWED:
        return {"valid": False, "reason": "Fuel mass flow collapsed or became non-positive."}

    of_ratio_local = m_dot_o_current / m_dot_f_local

    if not np.all(np.isfinite(of_ratio_local)):
        return {"valid": False, "reason": "Non-finite O/F encountered."}

    if np.min(of_ratio_local) < OF_MIN_ALLOWED or np.max(of_ratio_local) > OF_MAX_ALLOWED:
        return {
            "valid": False,
            "reason": f"Nonphysical O/F encountered: min = {np.min(of_ratio_local):.3f}, max = {np.max(of_ratio_local):.3f}"
        }

    if allow_expand:
        cea_profile, expanded = ensure_cea_range(cea_profile, of_ratio_local)
    else:
        expanded = False

    cea_of = cea_profile["of"]
    cea_cstar = cea_profile["cstar"]
    cea_cf = cea_profile["cf"]
    cea_isp = cea_profile["isp_s"]
    cea_ivac = cea_profile["ivac_s"]
    cea_tc = cea_profile["tc"]

    if np.min(of_ratio_local) < np.min(cea_of) or np.max(of_ratio_local) > np.max(cea_of):
        return {
            "valid": False,
            "reason": (
                f"O/F outside CEA interpolation range after attempted handling: "
                f"burn range = [{np.min(of_ratio_local):.3f}, {np.max(of_ratio_local):.3f}], "
                f"CEA range = [{np.min(cea_of):.3f}, {np.max(cea_of):.3f}]"
            )
        }

    cstar_local = np.interp(of_ratio_local, cea_of, cea_cstar)
    cf_local = np.interp(of_ratio_local, cea_of, cea_cf)
    isp_local = np.interp(of_ratio_local, cea_of, cea_isp)
    ivac_local = np.interp(of_ratio_local, cea_of, cea_ivac)
    tc_local = np.interp(of_ratio_local, cea_of, cea_tc)

    if not np.all(np.isfinite(cstar_local)) or not np.all(np.isfinite(cf_local)):
        return {"valid": False, "reason": "Interpolated CEA properties became non-finite."}

    cstar_ideal_local = cstar_local.copy()
    cf_ideal_local = cf_local.copy()
    isp_ideal_local = isp_local.copy()
    ivac_ideal_local = ivac_local.copy()

    if apply_efficiency:
        cstar_local = cstar_local * ETA_CSTAR
        cf_local = cf_local * ETA_CF
        isp_local = isp_local * ETA_CSTAR * ETA_CF
        ivac_local = ivac_local * ETA_CSTAR * ETA_CF

    thrust_local = (m_dot_o_current + m_dot_f_local) * cstar_local * cf_local

    if not np.all(np.isfinite(thrust_local)):
        return {"valid": False, "reason": "Non-finite thrust encountered."}

    avg_thrust_local = float(np.mean(thrust_local))
    avg_of_ratio_local = float(np.mean(of_ratio_local))

    return {
        "valid": True,
        "reason": "OK",
        "cea_profile": cea_profile,
        "expanded_cea": expanded,
        "L": L_current,
        "m_dot_o": m_dot_o_current,
        "Dp_t": Dp_t_local,
        "Gox": Gox_local,
        "r_dot_t": r_dot_t_local,
        "r_dot_t_mm": r_dot_t_mm_local,
        "m_dot_f": m_dot_f_local,
        "of_ratio": of_ratio_local,
        "cstar_profile": cstar_local,
        "cf_profile": cf_local,
        "isp_profile": isp_local,
        "ivac_profile": ivac_local,
        "cstar_real_profile": cstar_local,
        "cf_real_profile": cf_local,
        "isp_real_profile": isp_local,
        "ivac_real_profile": ivac_local,
        "cstar_ideal_profile": cstar_ideal_local,
        "cf_ideal_profile": cf_ideal_local,
        "isp_ideal_profile": isp_ideal_local,
        "ivac_ideal_profile": ivac_ideal_local,
        "tc_profile": tc_local,
        "thrust": thrust_local,
        "avg_thrust": avg_thrust_local,
        "avg_of_ratio": avg_of_ratio_local,
    }

# Build Initial CEA Profile
cea_profile = generate_cea_profile(PRESSURE_BAR, PE_BAR, CEA_OF_MIN, CEA_OF_MAX, CEA_OF_STEP)

print("\nCEA profile successfully generated.")
print(f"CEA O/F range: {cea_profile['of'].min():.3f} to {cea_profile['of'].max():.3f}")

print("\nEvaluating fixed engine configuration")
print(f"Fuel grain length input: {L:.6f} m")
print(f"Oxidizer mass flow input: {m_dot_o:.6f} kg/s")

final_eval = evaluate_current_design(
    L,
    m_dot_o,
    cea_profile,
    apply_efficiency=True,
    allow_expand=True
)

if not final_eval["valid"]:
    raise RuntimeError(f"Fixed configuration is invalid: {final_eval['reason']}")

cea_profile = final_eval["cea_profile"]
# Final solution
L = final_eval["L"]
m_dot_o = final_eval["m_dot_o"]
Dp_t = final_eval["Dp_t"]
Gox = final_eval["Gox"]
r_dot_t = final_eval["r_dot_t"]
r_dot_t_mm = final_eval["r_dot_t_mm"]
m_dot_f = final_eval["m_dot_f"]
of_ratio = final_eval["of_ratio"]
cstar_profile = final_eval["cstar_profile"]
cf_profile = final_eval["cf_profile"]
isp_profile = final_eval["isp_profile"]
ivac_profile = final_eval["ivac_profile"]
cstar_real_profile = final_eval["cstar_real_profile"]
cf_real_profile = final_eval["cf_real_profile"]
isp_real_profile = final_eval["isp_real_profile"]
ivac_real_profile = final_eval["ivac_real_profile"]
cstar_ideal_profile = final_eval["cstar_ideal_profile"]
cf_ideal_profile = final_eval["cf_ideal_profile"]
isp_ideal_profile = final_eval["isp_ideal_profile"]
ivac_ideal_profile = final_eval["ivac_ideal_profile"]
tc_profile = final_eval["tc_profile"]
thrust = final_eval["thrust"]
avg_thrust = final_eval["avg_thrust"]
avg_of_ratio = final_eval["avg_of_ratio"]

m_dot_total = m_dot_o + m_dot_f
fuel_mass_burned = np.trapezoid(m_dot_f, t)
ox_mass_burned = m_dot_o * (t[-1] - t[0])
total_prop_mass_burned = fuel_mass_burned + ox_mass_burned
impulse = np.trapezoid(thrust, t)

# Derived nozzle values. If no throat diameter is supplied, size the throat so
# average chamber pressure matches p_cham for this fixed configuration.
if THROAT_DIAMETER_MM is None:
    at = float(np.mean(m_dot_total * cstar_profile) / p_cham)
    throat_source = "auto-sized from average Pc"
else:
    dt_input = THROAT_DIAMETER_MM / 1000.0
    at = np.pi * (dt_input / 2.0) ** 2
    throat_source = "user input"

ae = epsilon * at
dt = 2.0 * np.sqrt(at / np.pi)
de = 2.0 * np.sqrt(ae / np.pi)
dt_mm = dt * 1000.0
de_mm = de * 1000.0

divergent_half_angle_deg = 15.0
divergent_half_angle_rad = np.deg2rad(divergent_half_angle_deg)
nozzle_length = (de - dt) / (2.0 * np.tan(divergent_half_angle_rad))
pc_profile_bar = thrust / (cf_profile * at) / 1e5

print("\nFinal Results:")
print(f'Average Thrust: {avg_thrust:.2f} N')
print(f'Average O/F Ratio: {avg_of_ratio:.4f}')
print(f'Fuel Grain Length: {L:.6f} m')
print(f'Fixed Oxidizer Mass Flow Rate: {m_dot_o:.6f} kg/s')
print(f'Nitrous Saturated Liquid Density at {T_n2o_C:.2f} C: {rho_ox:.3f} kg/m^3')
print(f'Nitrous Saturation Pressure at {T_n2o_C:.2f} C: {p_sat_n2o / 1e5:.3f} bar')
print(f'Nominal CEA Chamber Pressure: {p_cham / 1e5:.3f} bar')
print(f'Throat Area Source: {throat_source}')
print(f'Initial Chamber Pressure: {pc_profile_bar[0]:.3f} bar')
print(f'Final Chamber Pressure: {pc_profile_bar[-1]:.3f} bar')
print(f'Average Chamber Pressure: {np.mean(pc_profile_bar):.3f} bar')
print(f'Initial Port Diameter: {D_pi:.6f} m')
print(f'Final Port Diameter at burn end: {Dp_t[-1]:.6f} m')
print(f'Initial Fuel Mass Flow Rate: {m_dot_f[0]:.6f} kg/s')
print(f'Final Fuel Mass Flow Rate: {m_dot_f[-1]:.6f} kg/s')
print(f'Average Fuel Mass Flow Rate: {np.mean(m_dot_f):.6f} kg/s')
print(f'Initial O/F Ratio: {of_ratio[0]:.4f}')
print(f'Final O/F Ratio: {of_ratio[-1]:.4f}')
print(f'Initial Thrust: {thrust[0]:.2f} N')
print(f'Final Thrust: {thrust[-1]:.2f} N')
print(f'Peak Thrust: {np.max(thrust):.2f} N')
print(f'Minimum Thrust: {np.min(thrust):.2f} N')
print(f'Average C* Real/Effective: {np.mean(cstar_real_profile):.3f} m/s')
print(f'Average C* Ideal CEA: {np.mean(cstar_ideal_profile):.3f} m/s')
print(f'Average Cf Real/Effective: {np.mean(cf_real_profile):.5f}')
print(f'Average Cf Ideal CEA: {np.mean(cf_ideal_profile):.5f}')
print(f'Average Isp Real/Effective: {np.mean(isp_real_profile):.3f} s')
print(f'Average Isp Ideal CEA: {np.mean(isp_ideal_profile):.3f} s')
print(f'Average Ivac Real/Effective: {np.mean(ivac_real_profile):.3f} s')
print(f'Average Ivac Ideal CEA: {np.mean(ivac_ideal_profile):.3f} s')
print(f'Average Chamber Temperature: {np.mean(tc_profile):.2f} K')
print(f'Total Fuel Burned: {fuel_mass_burned:.4f} kg')
print(f'Total Oxidizer Burned: {ox_mass_burned:.4f} kg')
print(f'Total Propellant Burned: {total_prop_mass_burned:.4f} kg')
print(f'Total Impulse: {impulse:.2f} N*s')
print(f'Throat Area: {at:.6e} m^2')
print(f'Exit Area: {ae:.6e} m^2')
print(f'Throat Diameter: {dt_mm:.4f} mm')
print(f'Exit Diameter: {de_mm:.4f} mm')
print(f'Length between throat and exit: {nozzle_length:.6f} m')
print(f'Final CEA O/F interpolation range: {np.min(cea_profile["of"]):.3f} to {np.max(cea_profile["of"]):.3f}')

fig_engine, axes = plt.subplots(3, 3, figsize=(16, 11), constrained_layout=True)
fig_engine.suptitle('Engine Burn Summary', fontsize=16, fontweight='bold')

ax = axes[0, 0]
ax.plot(t, thrust, color='tab:blue', linewidth=1.8)
ax.axhline(avg_thrust, color='0.35', linestyle=':', linewidth=1.2, label='Actual avg')
ax.set_title('Thrust')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Thrust (N)')
ax.grid(True, alpha=0.3)
ax.legend(loc='best')

ax = axes[0, 1]
ax.plot(t, pc_profile_bar, color='tab:red', linewidth=1.8)
ax.set_title('Chamber Pressure')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Pc (bar)')
ax.grid(True, alpha=0.3)

ax = axes[0, 2]
ax.plot(t, m_dot_f, color='tab:brown', linewidth=1.8, label='Fuel')
ax.axhline(m_dot_o, color='tab:cyan', linestyle='--', linewidth=1.5, label='Oxidizer')
ax.plot(t, m_dot_total, color='tab:gray', linewidth=1.6, label='Total')
ax.set_title('Mass Flow Rates')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Mass flow (kg/s)')
ax.grid(True, alpha=0.3)
ax.legend(loc='best')

ax = axes[1, 0]
ax.plot(t, of_ratio, color='tab:orange', linewidth=1.8)
ax.axhline(avg_of_ratio, color='0.35', linestyle=':', linewidth=1.2, label='Actual avg')
ax.set_title('O/F Ratio')
ax.set_xlabel('Time (s)')
ax.set_ylabel('O/F')
ax.grid(True, alpha=0.3)
ax.legend(loc='best')

ax = axes[1, 1]
ax.plot(t, Dp_t * 1000.0, color='tab:green', linewidth=1.8)
ax.set_title('Port Diameter')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Port diameter (mm)')
ax.grid(True, alpha=0.3)

ax = axes[1, 2]
line_gox, = ax.plot(t, Gox, color='tab:purple', linewidth=1.8, label='Gox')
ax.set_title('Internal Ballistics')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Gox (kg/m^2/s)', color='tab:purple')
ax.tick_params(axis='y', labelcolor='tab:purple')
ax.grid(True, alpha=0.3)
ax2 = ax.twinx()
line_rdot, = ax2.plot(t, r_dot_t_mm, color='tab:pink', linewidth=1.8, label='Regression')
ax2.set_ylabel('Regression rate (mm/s)', color='tab:pink')
ax2.tick_params(axis='y', labelcolor='tab:pink')
ax.legend([line_gox, line_rdot], ['Gox', 'Regression'], loc='best')

ax = axes[2, 0]
line_cstar_real, = ax.plot(t, cstar_real_profile, color='tab:blue', linewidth=1.8, label='C* real')
line_cstar_ideal, = ax.plot(t, cstar_ideal_profile, color='tab:blue', linestyle='--', linewidth=1.4, label='C* ideal')
ax.set_title('Characteristic Velocity and Cf')
ax.set_xlabel('Time (s)')
ax.set_ylabel('C* (m/s)', color='tab:blue')
ax.tick_params(axis='y', labelcolor='tab:blue')
ax.grid(True, alpha=0.3)
ax2 = ax.twinx()
line_cf, = ax2.plot(t, cf_real_profile, color='tab:green', linewidth=1.6, label='Cf real')
ax2.set_ylabel('Cf', color='tab:green')
ax2.tick_params(axis='y', labelcolor='tab:green')
ax.legend([line_cstar_real, line_cstar_ideal, line_cf], ['C* Predict', 'C* ideal', 'Cf real'], loc='best')

ax = axes[2, 1]
ax.plot(t, isp_real_profile, color='tab:blue', linewidth=1.8, label='Isp Predict')
ax.plot(t, isp_ideal_profile, color='tab:blue', linestyle='--', linewidth=1.4, label='Isp ideal')
ax.plot(t, ivac_real_profile, color='tab:cyan', linewidth=1.6, label='Ivac Predict')
ax.plot(t, ivac_ideal_profile, color='tab:cyan', linestyle='--', linewidth=1.3, label='Ivac ideal')
ax.set_title('Specific Impulse')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Specific impulse (s)')
ax.grid(True, alpha=0.3)
ax.legend(loc='best')

ax = axes[2, 2]
ax.plot(t, tc_profile, color='tab:red', linewidth=1.8)
ax.set_title('Chamber Temperature')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Tc (K)')
ax.grid(True, alpha=0.3)

plt.show()
