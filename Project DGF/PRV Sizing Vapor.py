"""Modular API 520 Part I preliminary sizing for a gas/vapor PRV.

This module calculates the required effective discharge area for one
user-supplied relieving condition. It does not determine relief scenarios or
the required relieving mass flow.

API 520 SI equations implemented
---------------------------------
Critical flow:
    A = 13160 W sqrt(T Z / M) / (C Kd P1 Kb)

Subcritical flow:
    A = 17.9 W sqrt[Z T / (M P1 (P1 - P2))] / (F2 Kd)

A is in mm^2, W is in kg/hr, P1 and P2 are absolute kPa, T is K, and
M is kg/kmol. Kd = 0.975 is the API preliminary PRV value. Select and
verify the final valve using certified manufacturer capacity data.
"""

from __future__ import annotations

import math
from typing import TypedDict


# =============================================================================
# USER INPUTS - edit values only in this block
# =============================================================================
# Enter one documented relieving condition. All pressures are ABSOLUTE.
# The script gets k, Z, and molecular weight from CoolProp using fluid, P1, T.
USER_INPUTS = {
    "mass_flow_kg_s": 0.1642,                 # Required relieving flow [kg/s]
    "fluid": "Nitrogen",                      # CoolProp fluid name
    "relieving_pressure_pa_abs": 3_401_325.0, # P1 at PRV inlet [Pa abs]
    "relieving_temperature_k": 300.0,         # T at PRV inlet [K]
    "back_pressure_pa_abs": 101_325.0,        # P2 at PRV outlet [Pa abs]
    "kd": 0.975,                              # API preliminary PRV coefficient
    "kb": 1.0,                                # Back-pressure correction factor
}


class SizingResult(TypedDict):
    area_mm2: float
    equivalent_diameter_mm: float
    flow_regime: str
    critical_pressure_pa_abs: float
    k: float
    Z: float
    molecular_weight_kg_kmol: float


def api_c(k: float) -> float:
    """Return API 520 coefficient C for the heat-capacity ratio k."""
    _validate_k(k)
    return 520.0 * math.sqrt(
        k * (2.0 / (k + 1.0)) ** ((k + 1.0) / (k - 1.0))
    )


def api_f2(k: float, pressure_ratio: float) -> float:
    """Return API 520 subcritical-flow coefficient F2."""
    _validate_k(k)
    if not 0.0 < pressure_ratio < 1.0:
        raise ValueError("P2/P1 must be between zero and one.")
    return math.sqrt(
        (k / (k - 1.0))
        * pressure_ratio ** (2.0 / k)
        * (1.0 - pressure_ratio ** ((k - 1.0) / k))
        / (1.0 - pressure_ratio)
    )


def critical_pressure_ratio(k: float) -> float:
    """Return the ideal-gas critical pressure ratio Pcf/P1."""
    _validate_k(k)
    return (2.0 / (k + 1.0)) ** (k / (k - 1.0))


def properties_from_coolprop(
    fluid: str,
    pressure_pa_abs: float,
    temperature_k: float,
) -> tuple[float, float, float]:
    """Return (k, Z, molecular weight in kg/kmol) from CoolProp."""
    from CoolProp.CoolProp import PropsSI

    cp = PropsSI("CPMASS", "P", pressure_pa_abs, "T", temperature_k, fluid)
    cv = PropsSI("CVMASS", "P", pressure_pa_abs, "T", temperature_k, fluid)
    k = cp / cv
    z = PropsSI("Z", "P", pressure_pa_abs, "T", temperature_k, fluid)
    molecular_weight = PropsSI("M", fluid) * 1000.0
    return k, z, molecular_weight


def size_prv_vapor(
    mass_flow_kg_s: float,
    relieving_pressure_pa_abs: float,
    relieving_temperature_k: float,
    back_pressure_pa_abs: float,
    k: float,
    compressibility: float,
    molecular_weight_kg_kmol: float,
    *,
    kd: float = 0.975,
    kb: float = 1.0,
) -> SizingResult:
    """Calculate the preliminary API 520 effective area for one PRV.

    Supply properties at the PRV inlet relieving condition. ``kb`` is used in
    the critical-flow equation; obtain it from the applicable valve/back-
    pressure basis. For a conventional or pilot-operated valve where API 520
    permits no correction, use kb=1.0.
    """
    _validate_inputs(
        mass_flow_kg_s,
        relieving_pressure_pa_abs,
        relieving_temperature_k,
        back_pressure_pa_abs,
        k,
        compressibility,
        molecular_weight_kg_kmol,
        kd,
        kb,
    )

    mass_flow_kg_hr = mass_flow_kg_s * 3600.0
    p1_kpa = relieving_pressure_pa_abs / 1000.0
    p2_kpa = back_pressure_pa_abs / 1000.0
    pcf = critical_pressure_ratio(k) * relieving_pressure_pa_abs

    if back_pressure_pa_abs <= pcf:
        regime = "critical"
        area_mm2 = (
            13_160.0
            * mass_flow_kg_hr
            / (api_c(k) * kd * p1_kpa * kb)
            * math.sqrt(
                relieving_temperature_k
                * compressibility
                / molecular_weight_kg_kmol
            )
        )
    else:
        regime = "subcritical"
        f2 = api_f2(k, back_pressure_pa_abs / relieving_pressure_pa_abs)
        area_mm2 = (
            17.9
            * mass_flow_kg_hr
            / (f2 * kd)
            * math.sqrt(
                compressibility
                * relieving_temperature_k
                / (
                    molecular_weight_kg_kmol
                    * p1_kpa
                    * (p1_kpa - p2_kpa)
                )
            )
        )

    return {
        "area_mm2": area_mm2,
        "equivalent_diameter_mm": math.sqrt(4.0 * area_mm2 / math.pi),
        "flow_regime": regime,
        "critical_pressure_pa_abs": pcf,
        "k": k,
        "Z": compressibility,
        "molecular_weight_kg_kmol": molecular_weight_kg_kmol,
    }


def size_prv_vapor_with_coolprop(
    mass_flow_kg_s: float,
    fluid: str,
    relieving_pressure_pa_abs: float,
    relieving_temperature_k: float,
    back_pressure_pa_abs: float,
    *,
    kd: float = 0.975,
    kb: float = 1.0,
) -> SizingResult:
    """Convenience wrapper that obtains k, Z, and M from CoolProp."""
    k, z, molecular_weight = properties_from_coolprop(
        fluid, relieving_pressure_pa_abs, relieving_temperature_k
    )
    return size_prv_vapor(
        mass_flow_kg_s,
        relieving_pressure_pa_abs,
        relieving_temperature_k,
        back_pressure_pa_abs,
        k,
        z,
        molecular_weight,
        kd=kd,
        kb=kb,
    )


def _validate_k(k: float) -> None:
    if k <= 1.0:
        raise ValueError("k must be greater than 1.0 for these equations.")


def _validate_inputs(
    mass_flow_kg_s: float,
    p1: float,
    temperature: float,
    p2: float,
    k: float,
    z: float,
    molecular_weight: float,
    kd: float,
    kb: float,
) -> None:
    if mass_flow_kg_s <= 0.0:
        raise ValueError("mass_flow_kg_s must be greater than zero.")
    if p1 <= p2 or p2 < 0.0:
        raise ValueError("Pressures must satisfy P1 > P2 >= 0, in absolute Pa.")
    if temperature <= 0.0:
        raise ValueError("relieving_temperature_k must be greater than zero.")
    _validate_k(k)
    if z <= 0.0 or molecular_weight <= 0.0:
        raise ValueError("Z and molecular weight must be greater than zero.")
    if not 0.0 < kd <= 1.0 or not 0.0 < kb <= 1.0:
        raise ValueError("Kd and Kb must be in the interval (0, 1].")


if __name__ == "__main__":
    # Run the calculation using the values in USER_INPUTS above.
    result = size_prv_vapor_with_coolprop(
        **USER_INPUTS,
    )
    print("PRV preliminary API 520 sizing")
    for key, value in result.items():
        print(f"{key}: {value}")
