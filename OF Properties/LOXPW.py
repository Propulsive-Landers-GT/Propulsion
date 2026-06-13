from pathlib import Path
from datetime import datetime
import subprocess
import os
import re
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
CEA_EXE = Path(os.environ.get("CEA_EXE", SCRIPT_DIR / "cea.exe"))
CEA_DATA_DIR = Path(os.environ.get("CEA_DATA_DIR", SCRIPT_DIR / "data"))
BASE_DIR = SCRIPT_DIR

PRESSURE_BAR = 40.0
PE_BAR = 1.01325  # 1 atm in bar, exit pressure
OF_VALUES = np.arange(0.35, 10.01, 0.35) # range of O/F to iterate

OX_NAME = "LOX"
FUEL_NAME = "Paraffin"

# Reactant temperatures
LOX_TEMP_K = 90.17
PARAFFIN_TEMP_K = 298.15

G0 = 9.80665
DEBUG = True


def debug_print(msg: str) -> None:
    if DEBUG:
        print(msg)


def safe_float(x, default=np.nan):
    try:
        return float(str(x).replace("D", "E"))
    except Exception:
        return default


def compute_cstar_from_isp_cf(isp_s: float, cf: float) -> float:
    if np.isnan(isp_s) or np.isnan(cf) or cf <= 0:
        return np.nan
    return isp_s * G0 / cf


# CEA input and output
def build_input_text(pc_bar: float, of_ratio: float, pressure_ratio: float) -> str:
    return f"""problem rocket equilibrium
  p(bar)={pc_bar}
  o/f={of_ratio}
  pi/p={pressure_ratio}
reactants
  fuel=Paraffin wt%=100 t(k)={PARAFFIN_TEMP_K} h,cal/mol=-444600 C 73 H 124
  oxid=O2(L) wt%=100 t(k)={LOX_TEMP_K}
output siunits
end
"""


def run_cea_case(case_name: str, pc_bar: float, of_ratio: float, pressure_ratio: float, results_dir: Path) -> str:
    base = results_dir / case_name
    inp_path = base.with_suffix(".inp")
    out_path = base.with_suffix(".out")

    inp_text = build_input_text(pc_bar, of_ratio, pressure_ratio)
    inp_path.write_text(inp_text, encoding="utf-8")

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
    """
    Parse a table row like:
      CF                ...
      Isp, M/SEC        ...
      Ivac, M/SEC       ...
      T, K              ...
      GAMMAs            ...
      M, (1/n)          ...
      MACH NUMBER       ...
      SON VEL,M/SEC     ...
    """
    float_pattern = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"

    for line in out_text.splitlines():
        line_clean = line.strip().upper().replace("D", "E")
        if line_clean.startswith(row_name.upper()):
            nums = re.findall(float_pattern, line_clean)
            return [safe_float(x) for x in nums]

    return []


def parse_cea_output(out_text: str) -> dict:
    data = {
        "cf": np.nan,
        "isp_raw": np.nan,   # CEA often gives this in m/s in siunits
        "ivac_raw": np.nan,  # same here
        "tc": np.nan,
        "te": np.nan,
        "gamma_c": np.nan,
        "gamma_e": np.nan,
        "mw_c": np.nan,
        "mw_e": np.nan,
        "mach_e": np.nan,
        "sonic_e": np.nan,
    }

    # CF
    vals = parse_property_row(out_text, "CF")
    if len(vals) >= 1:
        data["cf"] = vals[-1]

    # Isp and Ivac in SI output are often in m/s
    vals = parse_property_row(out_text, "ISP")
    if len(vals) >= 1:
        data["isp_raw"] = vals[-1]

    vals = parse_property_row(out_text, "IVAC")
    if len(vals) >= 1:
        data["ivac_raw"] = vals[-1]

    # Temperatures
    vals = parse_property_row(out_text, "T, K")
    if len(vals) >= 1:
        data["tc"] = vals[0]
    if len(vals) >= 2:
        data["te"] = vals[-1]

    # Gamma
    vals = parse_property_row(out_text, "GAMM")
    if len(vals) >= 1:
        data["gamma_c"] = vals[0]
    if len(vals) >= 2:
        data["gamma_e"] = vals[-1]

    # Molecular weight / 1/n row in CEA output
    vals = parse_property_row(out_text, "M, (1/N)")
    if len(vals) >= 1:
        data["mw_c"] = vals[0]
    if len(vals) >= 2:
        data["mw_e"] = vals[-1]

    # Alternate MW row possibility
    if np.isnan(data["mw_c"]):
        vals = parse_property_row(out_text, "MOL WT")
        if len(vals) >= 1:
            data["mw_c"] = vals[0]
        if len(vals) >= 2:
            data["mw_e"] = vals[-1]

    # Mach
    vals = parse_property_row(out_text, "MACH")
    if len(vals) >= 1:
        data["mach_e"] = vals[-1]

    # Sonic velocity
    vals = parse_property_row(out_text, "SON VEL")
    if len(vals) >= 1:
        data["sonic_e"] = vals[-1]

    return data


# MAIN
if __name__ == "__main__":
    RUN_NAME = datetime.now().strftime("run_%Y-%m-%d_%H-%M-%S")
    RESULTS_DIR = BASE_DIR / "CEA_Results" / RUN_NAME
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Saving CEA files to: {RESULTS_DIR}")

    results = {
        "of": [],
        "cf": [],
        "isp_s": [],
        "ivac_s": [],
        "cstar": [],
        "tc": [],
        "te": [],
        "gamma_c": [],
        "gamma_e": [],
        "mw_c": [],
        "mw_e": [],
        "mach_e": [],
        "sonic_e": [],
    }

    pressure_ratio = PRESSURE_BAR / PE_BAR  # keep exit pressure at 1 atm

    for of_ratio in OF_VALUES:
        try:
            case_name = f"of_props_pc_{PRESSURE_BAR:.2f}_of_{of_ratio:.3f}".replace(".", "p")
            out_text = run_cea_case(case_name, PRESSURE_BAR, of_ratio, pressure_ratio, RESULTS_DIR)
            parsed = parse_cea_output(out_text)

            # CEA SI output often gives Isp/Ivac in m/s, so divide by g0 to convert to seconds
            isp_s = parsed["isp_raw"] / G0 if not np.isnan(parsed["isp_raw"]) else np.nan
            ivac_s = parsed["ivac_raw"] / G0 if not np.isnan(parsed["ivac_raw"]) else np.nan
            cstar = compute_cstar_from_isp_cf(isp_s, parsed["cf"])

            results["of"].append(of_ratio)
            results["cf"].append(parsed["cf"])
            results["isp_s"].append(isp_s)
            results["ivac_s"].append(ivac_s)
            results["cstar"].append(cstar)
            results["tc"].append(parsed["tc"])
            results["te"].append(parsed["te"])
            results["gamma_c"].append(parsed["gamma_c"])
            results["gamma_e"].append(parsed["gamma_e"])
            results["mw_c"].append(parsed["mw_c"])
            results["mw_e"].append(parsed["mw_e"])
            results["mach_e"].append(parsed["mach_e"])
            results["sonic_e"].append(parsed["sonic_e"])

            debug_print(
                f"O/F={of_ratio:.3f} | "
                f"Cf={parsed['cf']:.6f}, "
                f"Isp={isp_s:.3f} s, "
                f"Ivac={ivac_s:.3f} s, "
                f"C*={cstar:.3f} m/s, "
                f"Tc={parsed['tc']:.3f} K, "
                f"Te={parsed['te']:.3f} K, "
                f"gamma_c={parsed['gamma_c']:.6f}, "
                f"MW_c={parsed['mw_c']:.6f}"
            )

        except Exception as e:
            debug_print(f"FAILED at O/F = {of_ratio:.3f}: {e}")

            results["of"].append(of_ratio)
            results["cf"].append(np.nan)
            results["isp_s"].append(np.nan)
            results["ivac_s"].append(np.nan)
            results["cstar"].append(np.nan)
            results["tc"].append(np.nan)
            results["te"].append(np.nan)
            results["gamma_c"].append(np.nan)
            results["gamma_e"].append(np.nan)
            results["mw_c"].append(np.nan)
            results["mw_e"].append(np.nan)
            results["mach_e"].append(np.nan)
            results["sonic_e"].append(np.nan)

    # Convert to arrays
    of = np.array(results["of"])
    cf = np.array(results["cf"])
    isp_s = np.array(results["isp_s"])
    ivac_s = np.array(results["ivac_s"])
    cstar = np.array(results["cstar"])
    tc = np.array(results["tc"])
    te = np.array(results["te"])
    gamma_c = np.array(results["gamma_c"])
    gamma_e = np.array(results["gamma_e"])
    mw_c = np.array(results["mw_c"])
    mw_e = np.array(results["mw_e"])
    mach_e = np.array(results["mach_e"])
    sonic_e = np.array(results["sonic_e"])

    # Plotting
    fig, axs = plt.subplots(3, 2, figsize=(16, 9))
    fig.suptitle(
        f"{OX_NAME} / {FUEL_NAME} O/F Sweep at Pc = {PRESSURE_BAR:.1f} bar, Pe = {PE_BAR:.5f} bar",
        fontsize=16
    )

    ax = axs[0, 0]
    ax.plot(of, cf, linewidth=2)
    ax.set_title("Cf vs O/F")
    ax.set_xlabel("O/F Ratio")
    ax.set_ylabel("Exit Thrust Coefficient, Cf")
    ax.grid(True)

    ax = axs[0, 1]
    ax.plot(of, isp_s, linewidth=2, label="Isp")
    ax.plot(of, ivac_s, "--", linewidth=2, label="Ivac")
    ax.set_title("Specific Impulse vs O/F")
    ax.set_xlabel("O/F Ratio")
    ax.set_ylabel("Specific Impulse (s)")
    ax.legend()
    ax.grid(True)

    ax = axs[1, 0]
    ax.plot(of, cstar, linewidth=2)
    ax.set_title("C* vs O/F")
    ax.set_xlabel("O/F Ratio")
    ax.set_ylabel("Characteristic Velocity, C* (m/s)")
    ax.grid(True)

    ax = axs[1, 1]
    ax.plot(of, tc, linewidth=2, label="Tc")
    ax.plot(of, te, linewidth=2, label="Te")
    ax.set_title("Temperature vs O/F")
    ax.set_xlabel("O/F Ratio")
    ax.set_ylabel("Temperature (K)")
    ax.legend()
    ax.grid(True)

    ax = axs[2, 0]
    ax.plot(of, gamma_c, linewidth=2, label="Gamma chamber")
    ax.plot(of, gamma_e, linewidth=2, label="Gamma exit")
    ax.set_title("Gamma and Molecular Weight vs O/F")
    ax.set_xlabel("O/F Ratio")
    ax.set_ylabel("Gamma")
    ax.grid(True)

    ax2 = ax.twinx()
    ax2.plot(of, mw_c, ":", linewidth=2, label="MW chamber")
    ax2.plot(of, mw_e, "--", linewidth=2, label="MW exit")
    ax2.set_ylabel("Molecular Weight")

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="best")

    ax = axs[2, 1]
    ax.plot(of, mach_e, linewidth=2, label="Mach exit")
    ax.set_title("Exit Mach and Sonic Velocity vs O/F")
    ax.set_xlabel("O/F Ratio")
    ax.set_ylabel("Exit Mach Number")
    ax.grid(True)

    ax2 = ax.twinx()
    ax2.plot(of, sonic_e, "--", linewidth=2, label="Sonic velocity exit")
    ax2.set_ylabel("Exit Sonic Velocity (m/s)")

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="best")

    plt.tight_layout()
    plt.show()

    valid_cstar = ~np.isnan(cstar)
    if np.any(valid_cstar):
        idx = np.nanargmax(cstar)
        print("\nSummary")
        print(f"Peak C*      : {cstar[idx]:.2f} m/s at O/F = {of[idx]:.2f}")
        print(f"Cf at peak   : {cf[idx]:.4f}")
        print(f"Isp at peak  : {isp_s[idx]:.2f} s")
        print(f"Ivac at peak : {ivac_s[idx]:.2f} s")
        print(f"Tc at peak   : {tc[idx]:.2f} K")
        print(f"gamma_c peak : {gamma_c[idx]:.4f}")
        print(f"MW_c peak    : {mw_c[idx]:.4f}")
    else:
        print("\nNo valid C* points were computed.")