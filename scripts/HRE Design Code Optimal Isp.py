from pathlib import Path
from datetime import datetime
import subprocess
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# Inputs
target_avg_thrust = 1000    # N
target_avg_of_ratio = None      # Computed from CEA max-Isp O/F after the CEA sweep
OF_OPTIMIZATION_PROPERTY = "isp_s"  # Use "isp_s" for ambient Isp, "ivac_s" for vacuum Isp

# Initial guesses
L = 0.30                          # m
m_dot_o = 0.04                    # kg/s

# Convergence tolerance
tolerance = 0.001

t_tolerance_trigger = 0.0005
o_tolerance_trigger = None
tcount = 0
ocount = 0

# Iteration parameters
max_iterations = 9000
iteration = 0

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
SCRIPT_DIR = Path(__file__).resolve().parent
CEA_EXE = Path(os.environ.get("CEA_EXE", SCRIPT_DIR / "cea.exe"))
CEA_DATA_DIR = Path(os.environ.get("CEA_DATA_DIR", SCRIPT_DIR / "data"))
BASE_DIR = SCRIPT_DIR

PRESSURE_BAR = p_cham / 1e5
PE_BAR = 1.01325  # 1 atm in bar

# Start with a useful wide range
CEA_OF_MIN = 0.30
CEA_OF_MAX = 20.00
CEA_OF_STEP = 0.20

G0 = 9.80665
DEBUG = True

# Efficiency factors to apply to CEA ideal performance for more realistic estimates
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
    L_current = float(np.clip(L_current, L_MIN, L_MAX))
    m_dot_o_current = float(np.clip(m_dot_o_current, M_DOT_O_MIN, M_DOT_O_MAX))

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
        "tc_profile": tc_local,
        "thrust": thrust_local,
        "avg_thrust": avg_thrust_local,
        "avg_of_ratio": avg_of_ratio_local,
    }

def find_optimal_of_ratio(cea_profile, metric_key="isp_s"):
    """Return the O/F ratio that maximizes a CEA performance metric."""
    if metric_key not in cea_profile:
        raise KeyError(f"CEA profile does not contain metric '{metric_key}'.")

    of_values = np.asarray(cea_profile["of"], dtype=float)
    metric_values = np.asarray(cea_profile[metric_key], dtype=float)

    valid = (
        np.isfinite(of_values) &
        np.isfinite(metric_values) &
        (of_values >= OF_MIN_ALLOWED) &
        (of_values <= OF_MAX_ALLOWED)
    )

    if not np.any(valid):
        raise RuntimeError(f"No valid CEA points were available to optimize {metric_key}.")

    valid_indices = np.flatnonzero(valid)
    best_idx = int(valid_indices[np.argmax(metric_values[valid])])
    optimal_of = float(of_values[best_idx])
    optimal_metric = float(metric_values[best_idx])

    if 0 < best_idx < len(of_values) - 1:
        local_indices = np.array([best_idx - 1, best_idx, best_idx + 1])
        if np.all(valid[local_indices]):
            x = of_values[local_indices]
            y = metric_values[local_indices]
            a_fit, b_fit, c_fit = np.polyfit(x, y, 2)
            if a_fit < 0:
                peak_of = -b_fit / (2.0 * a_fit)
                if x[0] <= peak_of <= x[-1]:
                    optimal_of = float(peak_of)
                    optimal_metric = float(a_fit * peak_of * peak_of + b_fit * peak_of + c_fit)

    return optimal_of, optimal_metric

# Build Initial CEA Profile
cea_profile = generate_cea_profile(PRESSURE_BAR, PE_BAR, CEA_OF_MIN, CEA_OF_MAX, CEA_OF_STEP)

print("\nCEA profile successfully generated.")
print(f"CEA O/F range: {cea_profile['of'].min():.3f} to {cea_profile['of'].max():.3f}")

target_avg_of_ratio, optimal_isp_s = find_optimal_of_ratio(cea_profile, OF_OPTIMIZATION_PROPERTY)
o_tolerance_trigger = 0.2 / target_avg_of_ratio
optimization_label = "ambient Isp" if OF_OPTIMIZATION_PROPERTY == "isp_s" else "vacuum Isp"
print(f"Optimal CEA O/F for max {optimization_label}: {target_avg_of_ratio:.4f}")
print(f"Peak ideal CEA {optimization_label} at optimal O/F: {optimal_isp_s:.3f} s")


# Sanity checks
print("\nSanity checks")
for mdot_test in [M_DOT_O_MIN, 0.40, M_DOT_O_MAX]:
    for L_test in [0.10, 0.30, 0.60]:
        check = evaluate_current_design(L_test, mdot_test, cea_profile, apply_efficiency=True, allow_expand=True)
        if check["valid"]:
            print(
                f"m_dot_o={check['m_dot_o']:.4f} kg/s, L={check['L']:.4f} m -> "
                f"avg thrust={check['avg_thrust']:.2f} N, "
                f"avg O/F={check['avg_of_ratio']:.4f}, "
                f"O/F range=[{np.min(check['of_ratio']):.3f}, {np.max(check['of_ratio']):.3f}]"
            )
            cea_profile = check["cea_profile"]
        else:
            print(
                f"m_dot_o={mdot_test:.4f} kg/s, L={L_test:.4f} m -> INVALID: {check['reason']}"
            )


def length_for_target_avg_of(m_dot_o_current, target_of):
    """For fixed oxidizer flow, solve grain length directly from average O/F."""
    m_dot_o_current = float(np.clip(m_dot_o_current, M_DOT_O_MIN, M_DOT_O_MAX))

    Dp_t_local = (
        (D_pi ** (2 * n + 1)) +
        ((2 * n + 1) * (2 ** (2 * n + 1)) * a * (m_dot_o_current ** n) * t) / (np.pi ** n)
    ) ** (1.0 / (2 * n + 1))

    Gox_local = (4.0 * m_dot_o_current) / (np.pi * Dp_t_local * Dp_t_local)
    r_dot_t_local = a * (Gox_local ** n)
    fuel_flow_per_meter = rho_f * np.pi * Dp_t_local * r_dot_t_local

    if not np.all(np.isfinite(fuel_flow_per_meter)) or np.min(fuel_flow_per_meter) <= 0:
        return np.nan

    L_required = m_dot_o_current * float(np.mean(1.0 / fuel_flow_per_meter)) / target_of
    return float(L_required)


def evaluate_target_of_design(m_dot_o_current, cea_profile, apply_efficiency=True):
    L_required = length_for_target_avg_of(m_dot_o_current, target_avg_of_ratio)

    if not np.isfinite(L_required):
        return {"valid": False, "reason": "Could not compute a finite length for target average O/F."}

    if L_required < L_MIN or L_required > L_MAX:
        return {
            "valid": False,
            "reason": f"Length required for target average O/F is out of bounds: L = {L_required:.6f} m"
        }

    return evaluate_current_design(
        L_required,
        m_dot_o_current,
        cea_profile,
        apply_efficiency=apply_efficiency,
        allow_expand=True
    )


def solve_design_for_targets(cea_profile, apply_efficiency=True):
    """Bisection on oxidizer flow, with L solved analytically for target average O/F."""
    lo = M_DOT_O_MIN
    hi = M_DOT_O_MAX

    lo_eval = evaluate_target_of_design(lo, cea_profile, apply_efficiency=apply_efficiency)
    if lo_eval["valid"]:
        cea_profile = lo_eval["cea_profile"]

    hi_eval = evaluate_target_of_design(hi, cea_profile, apply_efficiency=apply_efficiency)
    if hi_eval["valid"]:
        cea_profile = hi_eval["cea_profile"]

    if not lo_eval["valid"] or not hi_eval["valid"]:
        details = f"low end: {lo_eval['reason']} | high end: {hi_eval['reason']}"
        raise RuntimeError(f"Could not bracket a valid design inside current limits. {details}")

    lo_error = lo_eval["avg_thrust"] - target_avg_thrust
    hi_error = hi_eval["avg_thrust"] - target_avg_thrust

    if lo_error == 0:
        return lo_eval, cea_profile, True, 1, [lo_eval]
    if hi_error == 0:
        return hi_eval, cea_profile, True, 1, [hi_eval]

    if lo_error * hi_error > 0:
        best = min([lo_eval, hi_eval], key=lambda item: abs(item["avg_thrust"] - target_avg_thrust))
        raise RuntimeError(
            "Target thrust is not bracketed inside current m_dot_o bounds. "
            f"At {lo:.3f} kg/s: {lo_eval['avg_thrust']:.2f} N; "
            f"at {hi:.3f} kg/s: {hi_eval['avg_thrust']:.2f} N. "
            f"Closest valid point is {best['avg_thrust']:.2f} N at "
            f"m_dot_o={best['m_dot_o']:.6f} kg/s, L={best['L']:.6f} m."
        )

    history = []
    best_eval = None
    converged_local = False

    for i in range(1, max_iterations + 1):
        mid = 0.5 * (lo + hi)
        mid_eval = evaluate_target_of_design(mid, cea_profile, apply_efficiency=apply_efficiency)

        if not mid_eval["valid"]:
            raise RuntimeError(f"Solver entered an invalid design at iteration {i}: {mid_eval['reason']}")

        cea_profile = mid_eval["cea_profile"]
        history.append(mid_eval)
        best_eval = mid_eval

        thrust_error = mid_eval["avg_thrust"] - target_avg_thrust
        of_error = mid_eval["avg_of_ratio"] - target_avg_of_ratio

        if abs(thrust_error) < tolerance and abs(of_error) < tolerance:
            converged_local = True
            break

        if lo_error * thrust_error <= 0:
            hi = mid
            hi_error = thrust_error
        else:
            lo = mid
            lo_error = thrust_error

    return best_eval, cea_profile, converged_local, len(history), history


print("\nSolving design at CEA-optimal average O/F")
final_eval, cea_profile, converged, iteration, solve_history = solve_design_for_targets(
    cea_profile,
    apply_efficiency=True
)

L_history = np.array([item["L"] for item in solve_history], dtype=float)
m_dot_o_history = np.array([item["m_dot_o"] for item in solve_history], dtype=float)
avg_thrust_history = np.array([item["avg_thrust"] for item in solve_history], dtype=float)
avg_of_ratio_history = np.array([item["avg_of_ratio"] for item in solve_history], dtype=float)

if converged:
    print(f"Target values achieved after {iteration} iterations.")
else:
    print(
        f"Stopped after {iteration} iterations with thrust error "
        f"{avg_thrust_history[-1] - target_avg_thrust:.6f} N and O/F error "
        f"{avg_of_ratio_history[-1] - target_avg_of_ratio:.6f}."
    )
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
tc_profile = final_eval["tc_profile"]
thrust = final_eval["thrust"]
avg_thrust = final_eval["avg_thrust"]
avg_of_ratio = final_eval["avg_of_ratio"]

# Derived nozzle values
avg_cf_final = float(np.mean(cf_profile))
at = target_avg_thrust / (p_cham * avg_cf_final)
ae = epsilon * at
dt = 2.0 * np.sqrt(at / np.pi)
de = 2.0 * np.sqrt(ae / np.pi)
dt_mm = dt * 1000.0
de_mm = de * 1000.0

divergent_half_angle_deg = 15.0
divergent_half_angle_rad = np.deg2rad(divergent_half_angle_deg)
nozzle_length = (de - dt) / (2.0 * np.tan(divergent_half_angle_rad))
pc_profile_bar = thrust / (cf_profile * at) / 1e5


# Final outputs:

m_dot_total = m_dot_o + m_dot_f
fuel_mass_burned = np.trapezoid(m_dot_f, t)
ox_mass_burned = m_dot_o * (t[-1] - t[0])
total_prop_mass_burned = fuel_mass_burned + ox_mass_burned
impulse = np.trapezoid(thrust, t)

print("\nFinal Results:")
print(f'Converged: {converged}')
print(f'Final Average Thrust: {avg_thrust:.2f} N')
print(f'Final Average O/F Ratio: {avg_of_ratio:.4f}')
print(f'CEA-Optimal Target Average O/F Ratio: {target_avg_of_ratio:.4f}')
print(f'Ideal CEA {optimization_label} at Optimal O/F: {optimal_isp_s:.3f} s')
print(f'Required Length of Fuel Grain: {L:.6f} m')
print(f'Required Oxidizer Mass Flow Rate: {m_dot_o:.6f} kg/s')
print(f'Nitrous Saturated Liquid Density at {T_n2o_C:.2f} C: {rho_ox:.3f} kg/m^3')
print(f'Nitrous Saturation Pressure at {T_n2o_C:.2f} C: {p_sat_n2o / 1e5:.3f} bar')
print(f'Nominal CEA Chamber Pressure: {p_cham / 1e5:.3f} bar')
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
print(f'Average C* (effective): {np.mean(cstar_profile):.3f} m/s')
print(f'Average Cf (effective): {np.mean(cf_profile):.5f}')
print(f'Average Isp (effective): {np.mean(isp_profile):.3f} s')
print(f'Average Ivac (effective): {np.mean(ivac_profile):.3f} s')
print(f'Average Chamber Temperature: {np.mean(tc_profile):.2f} K')
print(f'Total Fuel Burned: {fuel_mass_burned:.4f} kg')
print(f'Total Oxidizer Burned: {ox_mass_burned:.4f} kg')
print(f'Total Propellant Burned: {total_prop_mass_burned:.4f} kg')
print(f'Total Impulse: {impulse:.2f} N·s')
print(f'Throat Area: {at:.6e} m^2')
print(f'Exit Area: {ae:.6e} m^2')
print(f'Throat Diameter: {dt_mm:.4f} mm')
print(f'Exit Diameter: {de_mm:.4f} mm')
print(f'Length between throat and exit: {nozzle_length:.6f} m')
print(f'Final CEA O/F interpolation range: {np.min(cea_profile["of"]):.3f} to {np.max(cea_profile["of"]):.3f}')

# Plots
iterations = np.arange(1, len(avg_thrust_history) + 1)

fig_engine, axes = plt.subplots(3, 3, figsize=(16, 11), constrained_layout=True)
fig_engine.suptitle('Engine Burn Summary', fontsize=16, fontweight='bold')

ax = axes[0, 0]
ax.plot(t, thrust, color='tab:blue', linewidth=1.8)
ax.axhline(target_avg_thrust, color='tab:red', linestyle='--', linewidth=1.2, label='Target avg')
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
ax.axhline(target_avg_of_ratio, color='tab:red', linestyle='--', linewidth=1.2, label='Target avg')
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
line_cstar, = ax.plot(t, cstar_profile, color='tab:blue', linewidth=1.8, label='C*')
ax.set_title('Characteristic Velocity and Cf')
ax.set_xlabel('Time (s)')
ax.set_ylabel('C* (m/s)', color='tab:blue')
ax.tick_params(axis='y', labelcolor='tab:blue')
ax.grid(True, alpha=0.3)
ax2 = ax.twinx()
line_cf, = ax2.plot(t, cf_profile, color='tab:green', linewidth=1.8, label='Cf')
ax2.set_ylabel('Cf', color='tab:green')
ax2.tick_params(axis='y', labelcolor='tab:green')
ax.legend([line_cstar, line_cf], ['C*', 'Cf'], loc='best')

ax = axes[2, 1]
ax.plot(t, isp_profile, color='tab:blue', linewidth=1.8, label='Isp')
ax.plot(t, ivac_profile, color='tab:cyan', linestyle='--', linewidth=1.8, label='Ivac')
ax.set_title('Effective Specific Impulse')
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

fig_conv, axes_conv = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
fig_conv.suptitle('Solver Convergence', fontsize=14, fontweight='bold')

ax = axes_conv[0]
ax.plot(iterations, avg_thrust_history, color='tab:blue', linewidth=1.8)
ax.axhline(target_avg_thrust, color='tab:red', linestyle='--', linewidth=1.2, label='Target')
ax.set_title('Average Thrust')
ax.set_xlabel('Iteration')
ax.set_ylabel('Average thrust (N)')
ax.grid(True, alpha=0.3)
ax.legend(loc='best')

ax = axes_conv[1]
line_l, = ax.plot(iterations, L_history, color='tab:green', linewidth=1.8, label='L')
ax.set_title('Solved Variables')
ax.set_xlabel('Iteration')
ax.set_ylabel('Fuel grain length (m)', color='tab:green')
ax.tick_params(axis='y', labelcolor='tab:green')
ax.grid(True, alpha=0.3)
ax2 = ax.twinx()
line_mdot, = ax2.plot(iterations, m_dot_o_history, color='tab:purple', linewidth=1.8, label='m_dot_o')
ax2.set_ylabel('Oxidizer flow (kg/s)', color='tab:purple')
ax2.tick_params(axis='y', labelcolor='tab:purple')
ax.legend([line_l, line_mdot], ['L', 'm_dot_o'], loc='best')

plt.show()
