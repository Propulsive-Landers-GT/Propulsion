from __future__ import annotations

import csv
import json
import sys
from dataclasses import asdict, dataclass
from math import cos, log10, pi, radians, sin
from typing import Dict, Iterable, List, Optional, Union

try:
    from CoolProp.CoolProp import PropsSI
except ImportError as exc:
    raise SystemExit(
        "CoolProp is required for this script. Install it with: py -m pip install CoolProp"
    ) from exc


BAR_TO_PA = 100_000.0
IN_TO_M = 0.0254
FT_TO_M = 0.3048
LBM_TO_KG = 0.45359237
PA_TO_BAR = 1.0 / BAR_TO_PA
PA_TO_PSI = 0.0001450377377
PSI_TO_PA = 6894.757293
M3_S_TO_GPM = 15850.323141
WATER_DENSITY_KG_M3 = 999.0


# Crane 410, A-27 Appendix, "Pipe Friction Data for Schedule 40 Clean Commercial
# Steel Pipe with Flow in Zone of Complete Turbulence." Fitting K values such as
# K = 30 f_T use this table, not the actual Darcy friction factor calculated for
# pipe runs.
CRANE_FT_BY_NOMINAL_SIZE_IN: Dict[float, float] = {
    0.5: 0.026,
    0.75: 0.024,
    1.0: 0.022,
    1.25: 0.021,
    1.5: 0.020,
    2.0: 0.019,
    2.5: 0.018,
    3.0: 0.017,
    4.0: 0.016,
    6.0: 0.015,
    8.0: 0.014,
    12.0: 0.013,
    18.0: 0.012,
    30.0: 0.011,
}


CRANE_SIZE_LABELS_BY_NOMINAL_SIZE_IN: Dict[float, str] = {
    0.5: "1/2 in",
    0.75: "3/4 in",
    1.0: "1 in",
    1.25: "1 1/4 in",
    1.5: "1 1/2 in",
    2.0: "2 in",
    2.5: "2 1/2 in",
    3.0: "3 in",
    4.0: "4 in",
    6.0: "5-6 in",
    8.0: "8 in",
    12.0: "10-14 in",
    18.0: "16-22 in",
    30.0: "24-36 in",
}


CRANE_SIZE_VALUE_BY_LABEL: Dict[str, float] = {
    f"{label}  f_T={CRANE_FT_BY_NOMINAL_SIZE_IN[size]:.3f}": size
    for size, label in CRANE_SIZE_LABELS_BY_NOMINAL_SIZE_IN.items()
}


# Roughness reference: Engineering Toolbox, "Surface Roughness Coefficients of Pipes and Tubes."
# URL: https://www.engineeringtoolbox.com/surface-roughness-pipes-d_1127.html
# Values are absolute roughness k. The source table lists k in mm; these GUI
# values are converted to micrometers. For source ranges, the midpoint is used.
ROUGHNESS_UM_BY_LABEL: Dict[str, float] = {
    "Stainless steel, electro-polished - 0.45 um": 0.45,
    "Drawn copper/brass/aluminum, new - 1.5 um": 1.5,
    "Stainless steel, turned - 3.2 um": 3.2,
    "Stainless steel, bead blasted - 3.5 um": 3.5,
    "PVC/PE smooth plastic pipe - 4.25 um": 4.25,
    "Stretched steel - 15 um": 15.0,
    "Welded steel - 45 um": 45.0,
    "Commercial steel or wrought iron - 67.5 um": 67.5,
    "Galvanized steel - 150 um": 150.0,
}


# Dropdown-friendly catalog. Values with k_multiplier use K = multiplier * f_T.
# Values with fixed_k use K directly from Crane's representative entrance/exit
# data. Each Crane-backed dropdown entry has a source comment immediately above
# it so the manual page can be checked.
COMPONENT_CATALOG: Dict[str, Dict[str, object]] = {
    # Internal app component, not a Crane fitting table entry. Major loss uses Darcy-Weisbach.
    "pipe": {
        "label": "Straight pipe",
        "kind": "pipe",
    },
    # User/vendor input, not a Crane table entry. Use when a supplier gives K directly.
    "custom_k": {
        "label": "Custom K",
        "kind": "custom_k",
    },
    # User/vendor input, not a Crane table entry. Uses standard liquid Cv relation with US gpm and psi.
    "custom_cv": {
        "label": "Custom Cv (liquid)",
        "kind": "custom_cv",
    },
    # User/vendor/design input, not a Crane table entry. Uses incompressible dP from CdA.
    "custom_cda": {
        "label": "Custom CdA / orifice",
        "kind": "custom_cda",
    },
    # Crane 410, Chapter 2, pp. 2-14 through 2-16, "Hydraulic Resistance of Tees and Wyes";
    # A-30 Appendix points standard tees/wyes to these pages.
    "tee_90": {
        "label": "Crane tee - 90 deg side branch",
        "kind": "tee_wye",
        "fixed_angle_deg": 90.0,
    },
    # Crane 410, Chapter 2, pp. 2-14 through 2-16, "Hydraulic Resistance of Tees and Wyes";
    # A-30 Appendix points standard tees/wyes to these pages.
    "wye_angled": {
        "label": "Crane wye - angled branch",
        "kind": "tee_wye",
        "angle_options": (30.0, 45.0, 60.0),
    },
    # Crane 410, A-27 Appendix, "Sudden and Gradual Contraction"; uses Formula 1 or 2.
    "contraction": {
        "label": "Contraction / reducer",
        "kind": "transition",
        "transition_type": "contraction",
    },
    # Crane 410, A-27 Appendix, "Sudden and Gradual Enlargement"; uses Formula 3 or 4.
    "enlargement": {
        "label": "Enlargement / expansion",
        "kind": "transition",
        "transition_type": "enlargement",
    },
    # Crane 410, A-27 Appendix, "Formulas For Calculating K Factors for Valves and Fittings with Reduced Port."
    "reduced_port": {
        "label": "Reduced-port correction",
        "kind": "reduced_port",
    },
    # Crane 410, A-30 Appendix, "Standard Elbows": 90 deg, K = 30 f_T.
    "standard_elbow_90": {
        "label": "Standard 90 deg elbow",
        "kind": "crane_multiplier",
        "k_multiplier": 30.0,
    },
    # Crane 410, A-30 Appendix, "Standard Elbows": 45 deg, K = 16 f_T.
    "standard_elbow_45": {
        "label": "Standard 45 deg elbow",
        "kind": "crane_multiplier",
        "k_multiplier": 16.0,
    },
    # Crane 410, A-30 Appendix, "Close Pattern Return Bends": K = 50 f_T.
    "close_return_bend": {
        "label": "Close-pattern return bend",
        "kind": "crane_multiplier",
        "k_multiplier": 50.0,
    },
    # Crane 410, A-29 Appendix, "Ball Valves": full port/open, beta = 1, theta = 0, K = 3 f_T.
    "ball_valve_full_port": {
        "label": "Ball valve, full port, open",
        "kind": "crane_multiplier",
        "k_multiplier": 3.0,
    },
    # Crane 410, A-28 Appendix, "Gate Valves": wedge disc/double disc/plug type, beta = 1, theta = 0, K = 8 f_T.
    "gate_valve_wedge_open": {
        "label": "Gate valve, wedge disc, open",
        "kind": "crane_multiplier",
        "k_multiplier": 8.0,
    },
    # Crane 410, A-28 Appendix, "Globe and Angle Valves": globe valve, beta = 1, K = 340 f_T.
    "globe_valve": {
        "label": "Globe valve, open",
        "kind": "crane_multiplier",
        "k_multiplier": 340.0,
    },
    # Crane 410, A-28 Appendix, "Globe and Angle Valves": angle valve, beta = 1, K = 55 f_T.
    "angle_valve": {
        "label": "Angle valve, open",
        "kind": "crane_multiplier",
        "k_multiplier": 55.0,
    },
    # Crane 410, A-28 Appendix, "Swing Check Valves": full disc lift type, K = 100 f_T.
    "swing_check_valve_full_lift": {
        "label": "Swing check valve, full disc lift",
        "kind": "crane_multiplier",
        "k_multiplier": 100.0,
    },
    # Crane 410, A-28 Appendix, "Swing Check Valves": clearway/low-loss style shown as K = 50 f_T.
    "swing_check_valve_clearway": {
        "label": "Swing check valve, clearway type",
        "kind": "crane_multiplier",
        "k_multiplier": 50.0,
    },
    # Crane 410, A-28 Appendix, "Lift Check Valves": globe type, beta = 1, K = 600 f_T.
    "lift_check_valve_globe": {
        "label": "Lift check valve, globe type",
        "kind": "crane_multiplier",
        "k_multiplier": 600.0,
    },
    # Crane 410, A-28 Appendix, "Lift Check Valves": angle type, beta = 1, K = 55 f_T.
    "lift_check_valve_angle": {
        "label": "Lift check valve, angle type",
        "kind": "crane_multiplier",
        "k_multiplier": 55.0,
    },
    # Crane 410, A-29 Appendix, "Foot Valves with Strainer": poppet disc, K = 420 f_T.
    "foot_valve_poppet_strainer": {
        "label": "Foot valve with strainer, poppet disc",
        "kind": "crane_multiplier",
        "k_multiplier": 420.0,
    },
    # Crane 410, A-29 Appendix, "Foot Valves with Strainer": hinged disc, K = 75 f_T.
    "foot_valve_hinged_strainer": {
        "label": "Foot valve with strainer, hinged disc",
        "kind": "crane_multiplier",
        "k_multiplier": 75.0,
    },
    # Crane 410, A-30 Appendix, "Plug Valves and Cocks": straight way, beta = 1, K = 18 f_T.
    "plug_valve_straight": {
        "label": "Plug valve/cock, straight way",
        "kind": "crane_multiplier",
        "k_multiplier": 18.0,
    },
    # Crane 410, A-30 Appendix, "Plug Valves and Cocks": 3-way through flow, beta = 1, K = 30 f_T.
    "plug_valve_3way_through": {
        "label": "Plug valve/cock, 3-way through flow",
        "kind": "crane_multiplier",
        "k_multiplier": 30.0,
    },
    # Crane 410, A-30 Appendix, "Plug Valves and Cocks": 3-way branch flow, beta = 1, K = 90 f_T.
    "plug_valve_3way_branch": {
        "label": "Plug valve/cock, 3-way branch flow",
        "kind": "crane_multiplier",
        "k_multiplier": 90.0,
    },
    # Crane 410, A-29 Appendix, "Diaphragm Valves": weir type, beta = 1, K = 149 f_T.
    "diaphragm_valve_weir": {
        "label": "Diaphragm valve, weir",
        "kind": "crane_multiplier",
        "k_multiplier": 149.0,
    },
    # Crane 410, A-29 Appendix, "Diaphragm Valves": straight-through type, beta = 1, K = 39 f_T.
    "diaphragm_valve_straight": {
        "label": "Diaphragm valve, straight through",
        "kind": "crane_multiplier",
        "k_multiplier": 39.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": centric, size range 2-8 in, K = 45 f_T.
    "butterfly_valve_centric_small": {
        "label": "Butterfly valve, centric, 2-8 in",
        "kind": "crane_multiplier",
        "k_multiplier": 45.0,
        "min_size_in": 2.0,
        "max_size_in": 8.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": double offset, size range 2-8 in, K = 74 f_T.
    "butterfly_valve_double_offset_small": {
        "label": "Butterfly valve, double offset, 2-8 in",
        "kind": "crane_multiplier",
        "k_multiplier": 74.0,
        "min_size_in": 2.0,
        "max_size_in": 8.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": triple offset, size range 2-8 in, K = 218 f_T.
    "butterfly_valve_triple_offset_small": {
        "label": "Butterfly valve, triple offset, 2-8 in",
        "kind": "crane_multiplier",
        "k_multiplier": 218.0,
        "min_size_in": 2.0,
        "max_size_in": 8.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": centric, size range 10-14 in, K = 35 f_T.
    "butterfly_valve_centric_medium": {
        "label": "Butterfly valve, centric, 10-14 in",
        "kind": "crane_multiplier",
        "k_multiplier": 35.0,
        "min_size_in": 12.0,
        "max_size_in": 12.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": double offset, size range 10-14 in, K = 62 f_T.
    "butterfly_valve_double_offset_medium": {
        "label": "Butterfly valve, double offset, 10-14 in",
        "kind": "crane_multiplier",
        "k_multiplier": 62.0,
        "min_size_in": 12.0,
        "max_size_in": 12.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": triple offset, size range 10-14 in, K = 96 f_T.
    "butterfly_valve_triple_offset_medium": {
        "label": "Butterfly valve, triple offset, 10-14 in",
        "kind": "crane_multiplier",
        "k_multiplier": 96.0,
        "min_size_in": 12.0,
        "max_size_in": 12.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": centric, size range 16-24 in, K = 25 f_T.
    "butterfly_valve_centric_large": {
        "label": "Butterfly valve, centric, 16-24 in",
        "kind": "crane_multiplier",
        "k_multiplier": 25.0,
        "min_size_in": 18.0,
        "max_size_in": 30.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": double offset, size range 16-24 in, K = 43 f_T.
    "butterfly_valve_double_offset_large": {
        "label": "Butterfly valve, double offset, 16-24 in",
        "kind": "crane_multiplier",
        "k_multiplier": 43.0,
        "min_size_in": 18.0,
        "max_size_in": 30.0,
    },
    # Crane 410, A-29 Appendix, "Butterfly Valves": triple offset, size range 16-24 in, K = 55 f_T.
    "butterfly_valve_triple_offset_large": {
        "label": "Butterfly valve, triple offset, 16-24 in",
        "kind": "crane_multiplier",
        "k_multiplier": 55.0,
        "min_size_in": 18.0,
        "max_size_in": 30.0,
    },
    # Crane 410, A-30 Appendix, "Mitre Bends": 45 deg, K = 15 f_T.
    "miter_bend_45": {
        "label": "Miter bend, 45 deg",
        "kind": "crane_multiplier",
        "k_multiplier": 15.0,
    },
    # Crane 410, A-30 Appendix, "Mitre Bends": 15 deg, K = 4 f_T.
    "miter_bend_15": {
        "label": "Miter bend, 15 deg",
        "kind": "crane_multiplier",
        "k_multiplier": 4.0,
    },
    # Crane 410, A-30 Appendix, "Mitre Bends": 30 deg, K = 8 f_T.
    "miter_bend_30": {
        "label": "Miter bend, 30 deg",
        "kind": "crane_multiplier",
        "k_multiplier": 8.0,
    },
    # Crane 410, A-30 Appendix, "Mitre Bends": 60 deg, K = 25 f_T.
    "miter_bend_60": {
        "label": "Miter bend, 60 deg",
        "kind": "crane_multiplier",
        "k_multiplier": 25.0,
    },
    # Crane 410, A-30 Appendix, "Mitre Bends": 75 deg, K = 40 f_T.
    "miter_bend_75": {
        "label": "Miter bend, 75 deg",
        "kind": "crane_multiplier",
        "k_multiplier": 40.0,
    },
    # Crane 410, A-30 Appendix, "Mitre Bends": 90 deg, K = 60 f_T.
    "miter_bend_90": {
        "label": "Miter bend, 90 deg",
        "kind": "crane_multiplier",
        "k_multiplier": 60.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=1, K = 20 f_T.
    "pipe_bend_90_rd_1": {
        "label": "90 deg pipe bend, r/d=1",
        "kind": "crane_multiplier",
        "k_multiplier": 20.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=1.5, K = 14 f_T.
    "pipe_bend_90_rd_1_5": {
        "label": "90 deg pipe bend, r/d=1.5",
        "kind": "crane_multiplier",
        "k_multiplier": 14.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=2, K = 12 f_T.
    "pipe_bend_90_rd_2": {
        "label": "90 deg pipe bend, r/d=2",
        "kind": "crane_multiplier",
        "k_multiplier": 12.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=3, K = 12 f_T.
    "pipe_bend_90_rd_3": {
        "label": "90 deg pipe bend, r/d=3",
        "kind": "crane_multiplier",
        "k_multiplier": 12.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=4, K = 14 f_T.
    "pipe_bend_90_rd_4": {
        "label": "90 deg pipe bend, r/d=4",
        "kind": "crane_multiplier",
        "k_multiplier": 14.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=6, K = 17 f_T.
    "pipe_bend_90_rd_6": {
        "label": "90 deg pipe bend, r/d=6",
        "kind": "crane_multiplier",
        "k_multiplier": 17.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=8, K = 24 f_T.
    "pipe_bend_90_rd_8": {
        "label": "90 deg pipe bend, r/d=8",
        "kind": "crane_multiplier",
        "k_multiplier": 24.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=10, K = 30 f_T.
    "pipe_bend_90_rd_10": {
        "label": "90 deg pipe bend, r/d=10",
        "kind": "crane_multiplier",
        "k_multiplier": 30.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=12, K = 34 f_T.
    "pipe_bend_90_rd_12": {
        "label": "90 deg pipe bend, r/d=12",
        "kind": "crane_multiplier",
        "k_multiplier": 34.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=14, K = 38 f_T.
    "pipe_bend_90_rd_14": {
        "label": "90 deg pipe bend, r/d=14",
        "kind": "crane_multiplier",
        "k_multiplier": 38.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=16, K = 42 f_T.
    "pipe_bend_90_rd_16": {
        "label": "90 deg pipe bend, r/d=16",
        "kind": "crane_multiplier",
        "k_multiplier": 42.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows": r/d=20, K = 50 f_T.
    "pipe_bend_90_rd_20": {
        "label": "90 deg pipe bend, r/d=20",
        "kind": "crane_multiplier",
        "k_multiplier": 50.0,
    },
    # Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding 90 deg Elbows";
    # uses K_B = (n - 1)(0.25*pi*f_T*r/d + 0.5K) + K.
    "multiple_90_bends": {
        "label": "Multiple 90 deg pipe bends",
        "kind": "multiple_90_bends",
    },
    # Crane 410, A-30 Appendix, "Pipe Entrance": sharp-edged entrance, r/d = 0.00, K = 0.5.
    "pipe_entrance_sharp": {
        "label": "Pipe entrance, sharp edge",
        "kind": "fixed_k",
        "fixed_k": 0.5,
    },
    # Crane 410, A-30 Appendix, "Pipe Entrance": inward projecting entrance, K = 0.78.
    "pipe_entrance_projecting": {
        "label": "Pipe entrance, inward projecting",
        "kind": "fixed_k",
        "fixed_k": 0.78,
    },
    # Crane 410, A-30 Appendix, "Pipe Entrance": rounded/flush entrance, r/d = 0.02, K = 0.28.
    "pipe_entrance_rounded_r_d_0_02": {
        "label": "Pipe entrance, rounded r/d=0.02",
        "kind": "fixed_k",
        "fixed_k": 0.28,
    },
    # Crane 410, A-30 Appendix, "Pipe Entrance": rounded/flush entrance, r/d = 0.04, K = 0.24.
    "pipe_entrance_rounded_r_d_0_04": {
        "label": "Pipe entrance, rounded r/d=0.04",
        "kind": "fixed_k",
        "fixed_k": 0.24,
    },
    # Crane 410, A-30 Appendix, "Pipe Entrance": rounded/flush entrance, r/d = 0.06, K = 0.15.
    "pipe_entrance_rounded_r_d_0_06": {
        "label": "Pipe entrance, rounded r/d=0.06",
        "kind": "fixed_k",
        "fixed_k": 0.15,
    },
    # Crane 410, A-30 Appendix, "Pipe Entrance": rounded/flush entrance, r/d = 0.10, K = 0.09.
    "pipe_entrance_rounded_r_d_0_10": {
        "label": "Pipe entrance, rounded r/d=0.10",
        "kind": "fixed_k",
        "fixed_k": 0.09,
    },
    # Crane 410, A-30 Appendix, "Pipe Entrance": rounded/flush entrance, r/d = 0.15 and up, K = 0.04.
    "pipe_entrance_rounded_r_d_0_15": {
        "label": "Pipe entrance, rounded r/d=0.15 and up",
        "kind": "fixed_k",
        "fixed_k": 0.04,
    },
    # Crane 410, A-30 Appendix, "Pipe Exit": projecting/sharp-edged/rounded shown as K = 1.0.
    "pipe_exit": {
        "label": "Pipe exit",
        "kind": "fixed_k",
        "fixed_k": 1.0,
    },
}


FluidCatalogEntry = Dict[str, Union[str, float]]


# CoolProp is used where it has the fluid. Fixed-property entries are practical
# placeholders for propellants CoolProp does not normally include.
FLUID_CATALOG: Dict[str, FluidCatalogEntry] = {
    "LOX": {
        "label": "LOX - liquid oxygen",
        "provider": "coolprop",
        "coolprop_name": "Oxygen",
        "phase": "liquid",
        "default_temperature_k": 90.0,
        "default_pressure_bar": 60.0,
    },
    "GN2": {
        "label": "N2 - gaseous nitrogen",
        "provider": "coolprop",
        "coolprop_name": "Nitrogen",
        "phase": "gas",
        "default_temperature_k": 293.15,
        "default_pressure_bar": 60.0,
    },
    "LN2": {
        "label": "LN2 - liquid nitrogen",
        "provider": "coolprop",
        "coolprop_name": "Nitrogen",
        "phase": "liquid",
        "default_temperature_k": 77.0,
        "default_pressure_bar": 5.0,
    },
    "N2O_LIQ": {
        "label": "N2O - liquid nitrous oxide",
        "provider": "coolprop",
        "coolprop_name": "NitrousOxide",
        "phase": "liquid",
        "fallback_viscosity_pa_s": 0.00007,
        "default_temperature_k": 293.15,
        "default_pressure_bar": 60.0,
    },
    "H2O2_90": {
        "label": "H2O2 - 90 percent peroxide",
        "provider": "fixed",
        "density_kg_m3": 1390.0,
        "viscosity_pa_s": 0.00125,
        "default_temperature_k": 293.15,
        "default_pressure_bar": 30.0,
    },
    "WATER": {
        "label": "Water",
        "provider": "coolprop",
        "coolprop_name": "Water",
        "phase": "liquid",
        "default_temperature_k": 293.15,
        "default_pressure_bar": 10.0,
    },
    "ETHANOL": {
        "label": "Ethanol",
        "provider": "coolprop",
        "coolprop_name": "Ethanol",
        "phase": "liquid",
        "default_temperature_k": 293.15,
        "default_pressure_bar": 10.0,
    },
    "LCH4": {
        "label": "LCH4 - liquid methane",
        "provider": "coolprop",
        "coolprop_name": "Methane",
        "phase": "liquid",
        "default_temperature_k": 111.0,
        "default_pressure_bar": 20.0,
    },
    "GHE": {
        "label": "GHe - gaseous helium",
        "provider": "coolprop",
        "coolprop_name": "Helium",
        "phase": "gas",
        "default_temperature_k": 293.15,
        "default_pressure_bar": 200.0,
    },
    "GH2": {
        "label": "GH2 - gaseous hydrogen",
        "provider": "coolprop",
        "coolprop_name": "Hydrogen",
        "phase": "gas",
        "default_temperature_k": 293.15,
        "default_pressure_bar": 100.0,
    },
    "RP1": {
        "label": "RP-1 / kerosene approximate",
        "provider": "fixed",
        "density_kg_m3": 810.0,
        "viscosity_pa_s": 0.0018,
        "default_temperature_k": 293.15,
        "default_pressure_bar": 30.0,
    },
}


@dataclass(frozen=True)
class FlowCase:
    fluid: str = "LOX"
    temperature_k: float = 90.0
    inlet_pressure_pa: float = 60.0 * BAR_TO_PA
    mass_flow_kg_s: float = 2.0
    pipe_inner_diameter_m: float = 2.0 * IN_TO_M
    nominal_size_in: float = 2.0
    roughness_m: float = 15e-6


@dataclass(frozen=True)
class ComponentSpec:
    """Serializable component definition suitable for future GUI use."""

    name: str
    catalog_key: str
    quantity: int = 1
    length_m: Optional[float] = None
    inner_diameter_m: Optional[float] = None
    nominal_size_in: Optional[float] = None
    roughness_m: Optional[float] = None
    tee_flow_mode: Optional[str] = None
    tee_loss_path: Optional[str] = None
    tee_angle_deg: Optional[float] = None
    tee_q_branch_over_q_comb: Optional[float] = None
    tee_beta_branch: Optional[float] = None
    transition_angle_deg: Optional[float] = None
    transition_beta: Optional[float] = None
    custom_k: Optional[float] = None
    custom_cv: Optional[float] = None
    custom_cda_mm2: Optional[float] = None
    bend_count: Optional[int] = None
    bend_r_over_d: Optional[float] = None
    reduced_port_formula: Optional[int] = None
    reduced_port_k1: Optional[float] = None
    reduced_port_beta: Optional[float] = None
    reduced_port_angle_deg: Optional[float] = None


@dataclass(frozen=True)
class FluidProperties:
    density_kg_m3: float
    viscosity_pa_s: float


@dataclass(frozen=True)
class LadderRow:
    component: str
    quantity: int
    pressure_in_pa: float
    delta_p_pa: float
    pressure_out_pa: float
    density_kg_m3: float
    velocity_m_s: float
    reynolds: float
    darcy_f: Optional[float]
    k_total: Optional[float]


def fluid_properties(case: FlowCase, pressure_pa: float) -> FluidProperties:
    entry = FLUID_CATALOG[case.fluid]
    provider = str(entry["provider"])

    if provider == "fixed":
        return FluidProperties(
            density_kg_m3=float(entry["density_kg_m3"]),
            viscosity_pa_s=float(entry["viscosity_pa_s"]),
        )

    coolprop_name = str(entry["coolprop_name"])
    density = PropsSI("D", "T", case.temperature_k, "P", pressure_pa, coolprop_name)
    try:
        viscosity = PropsSI("V", "T", case.temperature_k, "P", pressure_pa, coolprop_name)
    except ValueError:
        if "fallback_viscosity_pa_s" not in entry:
            raise
        viscosity = float(entry["fallback_viscosity_pa_s"])
    return FluidProperties(density_kg_m3=density, viscosity_pa_s=viscosity)


def liquid_saturation_pressure_pa(case: FlowCase) -> Optional[float]:
    entry = FLUID_CATALOG[case.fluid]
    if entry.get("provider") != "coolprop" or entry.get("phase") != "liquid":
        return None

    coolprop_name = str(entry["coolprop_name"])
    try:
        critical_temperature_k = PropsSI("Tcrit", coolprop_name)
    except ValueError:
        return None
    if case.temperature_k >= critical_temperature_k:
        return None

    try:
        return PropsSI("P", "T", case.temperature_k, "Q", 0, coolprop_name)
    except ValueError:
        return None


def assert_liquid_pressure_valid(case: FlowCase, pressure_pa: float, location: str) -> None:
    saturation_pressure_pa = liquid_saturation_pressure_pa(case)
    if saturation_pressure_pa is None:
        return
    if pressure_pa <= saturation_pressure_pa:
        fluid_label = FLUID_CATALOG[case.fluid]["label"]
        raise ValueError(
            f"{fluid_label} at {case.temperature_k:g} K is below saturation pressure at {location}. "
            f"P = {pressure_pa * PA_TO_BAR:.3f} bar, Psat = {saturation_pressure_pa * PA_TO_BAR:.3f} bar. "
            "The liquid pressure-ladder model is no longer valid; increase pressure, lower temperature, "
            "increase pipe size, reduce flow, or use a two-phase/compressible model."
        )


def area_from_diameter(diameter_m: float) -> float:
    return pi * diameter_m**2 / 4.0


def velocity_from_mass_flow(mass_flow_kg_s: float, density_kg_m3: float, diameter_m: float) -> float:
    return mass_flow_kg_s / (density_kg_m3 * area_from_diameter(diameter_m))


def reynolds_number(density_kg_m3: float, velocity_m_s: float, diameter_m: float, viscosity_pa_s: float) -> float:
    return density_kg_m3 * velocity_m_s * diameter_m / viscosity_pa_s


def darcy_friction_factor(reynolds: float, diameter_m: float, roughness_m: float) -> float:
    if reynolds <= 0:
        raise ValueError("Reynolds number must be positive.")
    if reynolds < 2300.0:
        return 64.0 / reynolds

    relative_roughness = roughness_m / diameter_m
    # Haaland explicit approximation to Colebrook. Good enough for sizing and
    # avoids iterative solver weirdness in GUI use.
    return (-1.8 * log10((relative_roughness / 3.7) ** 1.11 + 6.9 / reynolds)) ** -2


def crane_ft(nominal_size_in: float) -> float:
    if nominal_size_in in CRANE_FT_BY_NOMINAL_SIZE_IN:
        return CRANE_FT_BY_NOMINAL_SIZE_IN[nominal_size_in]

    allowed_sizes = ", ".join(CRANE_SIZE_LABELS_BY_NOMINAL_SIZE_IN[size] for size in CRANE_FT_BY_NOMINAL_SIZE_IN)
    raise ValueError(
        f"{nominal_size_in:g} in is not in the Crane f_T table. "
        f"Select one of: {allowed_sizes}."
    )


# Crane 410, Chapter 2, p. 2-15, Eq. 2-35 and Table 2-1,
# "Constants for Equation 2-35" for converging flow through tees and wyes.
CONVERGING_TEE_CONSTANTS: Dict[float, Dict[str, Optional[tuple[Optional[float], float, float, float]]]] = {
    30.0: {
        "branch": (None, 1.0, 2.0, 1.74),
        "run": (1.0, 0.0, 1.0, 1.74),
    },
    45.0: {
        "branch": (None, 1.0, 2.0, 1.41),
        "run": (1.0, 0.0, 1.0, 1.41),
    },
    60.0: {
        "branch": (2.0, 1.0, 2.0, 1.0),
        "run": (1.0, 0.0, 1.0, 1.0),
    },
    90.0: {
        "branch": (1.0, 1.0, 2.0, 0.0),
        "run": None,
    },
}


# Crane 410, Chapter 2, p. 2-15, Table 2-2, "Values of C for Equation 2-35."
def tee_c_converging(q_ratio: float, beta_branch: float) -> float:
    beta_sq = beta_branch**2
    if beta_sq <= 0.35:
        return 1.0
    if q_ratio <= 0.35:
        return 0.9 * (1.0 - q_ratio)
    return 0.55


# Crane 410, Chapter 2, p. 2-15, Table 2-4, "Values of G for Equation 2-37."
def tee_g_diverging(q_ratio: float, beta_branch: float) -> float:
    beta_sq = beta_branch**2
    if beta_sq <= 0.35:
        if q_ratio <= 0.6:
            return 1.1 - 0.7 * q_ratio
        return 0.85
    if q_ratio <= 0.4:
        return 1.0 - 0.6 * q_ratio
    return 0.6


# Crane 410, Chapter 2, p. 2-15, Table 2-5, "Values of M for Equation 2-38."
def tee_m_diverging(q_ratio: float, beta_branch: float) -> float:
    beta_sq = beta_branch**2
    if beta_sq <= 0.4:
        return 0.4
    if q_ratio <= 0.5:
        return 2.0 * (2.0 * q_ratio - 1.0)
    return 0.3 * (2.0 * q_ratio - 1.0)


def crane_tee_wye_k(
    flow_mode: str,
    loss_path: str,
    angle_deg: float,
    q_branch_over_q_comb: float,
    beta_branch: float,
) -> float:
    # Crane 410, Chapter 2, p. 2-15:
    # converging flow uses Eq. 2-35 and Eq. 2-36 for 90 deg run;
    # diverging flow uses Eq. 2-37 for branch and Eq. 2-38 for run.
    if flow_mode not in {"converging", "diverging"}:
        raise ValueError("Tee/wye flow mode must be converging or diverging.")
    if loss_path not in {"branch", "run"}:
        raise ValueError("Tee/wye loss path must be branch or run.")
    if angle_deg not in {30.0, 45.0, 60.0, 90.0}:
        raise ValueError("Tee/wye angle must be 30, 45, 60, or 90 degrees.")
    if not 0.0 < q_branch_over_q_comb < 1.0:
        raise ValueError("Tee/wye Q_branch/Q_comb must be between 0 and 1.")
    if not 0.0 < beta_branch <= 1.0:
        raise ValueError("Tee/wye beta_branch must be greater than 0 and no more than 1.")

    q_ratio = q_branch_over_q_comb
    beta_sq = beta_branch**2
    velocity_ratio = q_ratio / beta_sq

    if flow_mode == "converging":
        if loss_path == "run" and angle_deg == 90.0:
            return 1.55 * q_ratio - q_ratio**2

        constants = CONVERGING_TEE_CONSTANTS[angle_deg][loss_path]
        if constants is None:
            raise ValueError("No converging tee/wye constants are available for this selection.")
        c_raw, d_const, e_const, f_const = constants
        c_const = tee_c_converging(q_ratio, beta_branch) if c_raw is None else c_raw
        return c_const * (
            1.0
            + d_const * velocity_ratio**2
            - e_const * (1.0 - q_ratio) ** 2
            - f_const * velocity_ratio**2
        )

    if loss_path == "run":
        return tee_m_diverging(q_ratio, beta_branch) * q_ratio**2

    if angle_deg == 90.0 and beta_sq > 2.0 / 3.0:
        g_const = 1.0 + 0.3 * q_ratio**2
        h_const = 0.3
        j_const = 0.0
    else:
        g_const = tee_g_diverging(q_ratio, beta_branch)
        h_const = 1.0
        j_const = 2.0

    return g_const * (
        1.0
        + h_const * velocity_ratio**2
        - j_const * velocity_ratio * cos(radians(angle_deg))
    )


# Crane 410, A-27 Appendix, "Formulas For Calculating K Factors for Valves and
# Fittings with Reduced Port": contraction uses Formula 1 or 2; enlargement uses
# Formula 3 or 4. The A-27 contraction/enlargement diagrams select the formula
# based on included angle.
def crane_transition_k(transition_type: str, angle_deg: float, beta_small_over_large: float) -> float:
    if transition_type not in {"contraction", "enlargement"}:
        raise ValueError("Transition type must be contraction or enlargement.")
    if not 0.0 < beta_small_over_large < 1.0:
        raise ValueError("Transition diameter ratio must be between 0 and 1.")
    if not 0.0 < angle_deg <= 180.0:
        raise ValueError("Transition angle must be greater than 0 and no more than 180 degrees.")

    beta = beta_small_over_large
    beta_sq = beta**2
    beta_fourth = beta**4
    sin_half_angle = sin(radians(angle_deg) / 2.0)

    if transition_type == "contraction":
        if angle_deg <= 45.0:
            return 0.8 * sin_half_angle * (1.0 - beta_sq) / beta_fourth
        return 0.5 * (1.0 - beta_sq) * sin_half_angle**0.5 / beta_fourth

    if angle_deg <= 45.0:
        return 2.6 * sin_half_angle * (1.0 - beta_sq) ** 2 / beta_fourth
    return (1.0 - beta_sq) ** 2 / beta_fourth


PIPE_BEND_90_MULTIPLIER_BY_R_OVER_D: Dict[float, float] = {
    1.0: 20.0,
    1.5: 14.0,
    2.0: 12.0,
    3.0: 12.0,
    4.0: 14.0,
    6.0: 17.0,
    8.0: 24.0,
    10.0: 30.0,
    12.0: 34.0,
    14.0: 38.0,
    16.0: 42.0,
    20.0: 50.0,
}


# Crane 410, A-30 Appendix, "90 deg Pipe Bends and Flanged or Butt-Welding
# 90 deg Elbows." Applies the printed multiple-bend equation.
def crane_multiple_90_bends_k(bend_count: int, r_over_d: float, nominal_size_in: float) -> float:
    if bend_count < 1:
        raise ValueError("Number of 90 deg bends must be at least 1.")
    if r_over_d not in PIPE_BEND_90_MULTIPLIER_BY_R_OVER_D:
        allowed = ", ".join(f"{value:g}" for value in PIPE_BEND_90_MULTIPLIER_BY_R_OVER_D)
        raise ValueError(f"r/d must be one of the Crane table values: {allowed}.")

    ft = crane_ft(nominal_size_in)
    single_bend_k = PIPE_BEND_90_MULTIPLIER_BY_R_OVER_D[r_over_d] * ft
    return (bend_count - 1) * (0.25 * pi * ft * r_over_d + 0.5 * single_bend_k) + single_bend_k


# Crane 410, A-27 Appendix, "Formulas For Calculating K Factors for Valves and
# Fittings with Reduced Port." Formula numbers match the A-27 figure.
def crane_reduced_port_k(formula: int, k1: float, beta: float, angle_deg: float) -> float:
    if formula not in {1, 2, 3, 4, 5, 6, 7}:
        raise ValueError("Reduced-port formula must be 1 through 7.")
    if k1 < 0.0:
        raise ValueError("Reduced-port K1 must be non-negative.")
    if not 0.0 < beta <= 1.0:
        raise ValueError("Reduced-port beta must be greater than 0 and no more than 1.")
    if not 0.0 <= angle_deg <= 180.0:
        raise ValueError("Reduced-port angle must be between 0 and 180 degrees.")

    beta_sq = beta**2
    beta_fourth = beta**4
    sin_half_angle = sin(radians(angle_deg) / 2.0)
    one_minus_beta_sq = 1.0 - beta_sq

    additive_1 = 0.8 * sin_half_angle * one_minus_beta_sq / beta_fourth
    additive_2 = 0.5 * one_minus_beta_sq * sin_half_angle**0.5 / beta_fourth
    additive_3 = 2.6 * sin_half_angle * one_minus_beta_sq**2 / beta_fourth
    additive_4 = one_minus_beta_sq**2 / beta_fourth
    k1_scaled = k1 / beta_fourth

    if formula == 1:
        return k1_scaled + additive_1
    if formula == 2:
        return k1_scaled + additive_2
    if formula == 3:
        return k1_scaled + additive_3
    if formula == 4:
        return k1_scaled + additive_4
    if formula == 5:
        return k1_scaled + additive_1 + additive_3
    if formula == 6:
        return k1_scaled + additive_2 + additive_4

    return (k1 + beta * (0.5 * one_minus_beta_sq + one_minus_beta_sq**2)) / beta_fourth


def component_k(spec: ComponentSpec, case: FlowCase) -> float:
    catalog_entry = COMPONENT_CATALOG[spec.catalog_key]
    kind = catalog_entry["kind"]

    if kind == "fixed_k":
        return float(catalog_entry["fixed_k"]) * spec.quantity

    if kind == "custom_k":
        if spec.custom_k is None:
            raise ValueError(f"Custom K component {spec.name!r} is missing K.")
        if spec.custom_k < 0.0:
            raise ValueError("Custom K must be non-negative.")
        return spec.custom_k * spec.quantity

    if kind == "crane_multiplier":
        nominal_size_in = spec.nominal_size_in or case.nominal_size_in
        min_size_in = catalog_entry.get("min_size_in")
        max_size_in = catalog_entry.get("max_size_in")
        if min_size_in is not None and nominal_size_in < float(min_size_in):
            raise ValueError(
                f"{catalog_entry['label']} is only valid for the listed Crane size range, "
                f"not {nominal_size_in:g} in."
            )
        if max_size_in is not None and nominal_size_in > float(max_size_in):
            raise ValueError(
                f"{catalog_entry['label']} is only valid for the listed Crane size range, "
                f"not {nominal_size_in:g} in."
            )
        return float(catalog_entry["k_multiplier"]) * crane_ft(nominal_size_in) * spec.quantity

    if kind == "tee_wye":
        if (
            spec.tee_flow_mode is None
            or spec.tee_loss_path is None
            or spec.tee_angle_deg is None
            or spec.tee_q_branch_over_q_comb is None
            or spec.tee_beta_branch is None
        ):
            raise ValueError(f"Tee/wye component {spec.name!r} is missing required geometry fields.")
        return (
            crane_tee_wye_k(
                spec.tee_flow_mode,
                spec.tee_loss_path,
                spec.tee_angle_deg,
                spec.tee_q_branch_over_q_comb,
                spec.tee_beta_branch,
            )
            * spec.quantity
        )

    if kind == "transition":
        if spec.transition_angle_deg is None or spec.transition_beta is None:
            raise ValueError(f"Transition component {spec.name!r} is missing angle or diameter ratio.")
        return (
            crane_transition_k(
                str(catalog_entry["transition_type"]),
                spec.transition_angle_deg,
                spec.transition_beta,
            )
            * spec.quantity
        )

    if kind == "multiple_90_bends":
        if spec.bend_count is None or spec.bend_r_over_d is None:
            raise ValueError(f"Multiple-bend component {spec.name!r} is missing bend count or r/d.")
        return crane_multiple_90_bends_k(
            spec.bend_count,
            spec.bend_r_over_d,
            spec.nominal_size_in or case.nominal_size_in,
        ) * spec.quantity

    if kind == "reduced_port":
        if (
            spec.reduced_port_formula is None
            or spec.reduced_port_k1 is None
            or spec.reduced_port_beta is None
            or spec.reduced_port_angle_deg is None
        ):
            raise ValueError(f"Reduced-port component {spec.name!r} is missing formula, K1, beta, or angle.")
        return crane_reduced_port_k(
            spec.reduced_port_formula,
            spec.reduced_port_k1,
            spec.reduced_port_beta,
            spec.reduced_port_angle_deg,
        ) * spec.quantity

    raise ValueError(f"Component {spec.catalog_key!r} is not a minor-loss component.")


def pressure_drop_custom_cv(spec: ComponentSpec, case: FlowCase, props: FluidProperties) -> float:
    if spec.custom_cv is None or spec.custom_cv <= 0.0:
        raise ValueError(f"Custom Cv component {spec.name!r} needs Cv > 0.")
    volumetric_flow_m3_s = case.mass_flow_kg_s / props.density_kg_m3
    flow_gpm = volumetric_flow_m3_s * M3_S_TO_GPM
    specific_gravity = props.density_kg_m3 / WATER_DENSITY_KG_M3
    delta_p_psi = specific_gravity * (flow_gpm / spec.custom_cv) ** 2
    return delta_p_psi * PSI_TO_PA * spec.quantity


def pressure_drop_custom_cda(spec: ComponentSpec, case: FlowCase, props: FluidProperties) -> float:
    if spec.custom_cda_mm2 is None or spec.custom_cda_mm2 <= 0.0:
        raise ValueError(f"Custom CdA component {spec.name!r} needs CdA > 0.")
    cda_m2 = spec.custom_cda_mm2 * 1e-6
    effective_velocity = case.mass_flow_kg_s / (props.density_kg_m3 * cda_m2)
    return props.density_kg_m3 * effective_velocity**2 / 2.0 * spec.quantity


def pressure_drop_pipe(spec: ComponentSpec, case: FlowCase, props: FluidProperties) -> tuple[float, float, float, float]:
    diameter_m = spec.inner_diameter_m or case.pipe_inner_diameter_m
    roughness_m = spec.roughness_m or case.roughness_m
    length_m = spec.length_m
    if length_m is None:
        raise ValueError(f"Pipe component {spec.name!r} needs length_m.")

    velocity = velocity_from_mass_flow(case.mass_flow_kg_s, props.density_kg_m3, diameter_m)
    reynolds = reynolds_number(props.density_kg_m3, velocity, diameter_m, props.viscosity_pa_s)
    darcy_f = darcy_friction_factor(reynolds, diameter_m, roughness_m)
    dynamic_pressure = props.density_kg_m3 * velocity**2 / 2.0
    delta_p = darcy_f * (length_m / diameter_m) * dynamic_pressure * spec.quantity
    return delta_p, velocity, reynolds, darcy_f


def pressure_drop_minor(spec: ComponentSpec, case: FlowCase, props: FluidProperties) -> tuple[float, float, float, float]:
    diameter_m = spec.inner_diameter_m or case.pipe_inner_diameter_m
    velocity = velocity_from_mass_flow(case.mass_flow_kg_s, props.density_kg_m3, diameter_m)
    reynolds = reynolds_number(props.density_kg_m3, velocity, diameter_m, props.viscosity_pa_s)
    dynamic_pressure = props.density_kg_m3 * velocity**2 / 2.0
    kind = COMPONENT_CATALOG[spec.catalog_key]["kind"]
    if kind == "custom_cv":
        delta_p = pressure_drop_custom_cv(spec, case, props)
        k_total = delta_p / dynamic_pressure
    elif kind == "custom_cda":
        delta_p = pressure_drop_custom_cda(spec, case, props)
        k_total = delta_p / dynamic_pressure
    else:
        k_total = component_k(spec, case)
        delta_p = k_total * dynamic_pressure
    return delta_p, velocity, reynolds, k_total


def run_pressure_ladder(case: FlowCase, components: Iterable[ComponentSpec]) -> List[LadderRow]:
    rows: List[LadderRow] = []
    pressure_pa = case.inlet_pressure_pa

    for spec in components:
        catalog_entry = COMPONENT_CATALOG[spec.catalog_key]
        assert_liquid_pressure_valid(case, pressure_pa, f"upstream of {spec.name}")
        props = fluid_properties(case, pressure_pa)

        if catalog_entry["kind"] == "pipe":
            delta_p, velocity, reynolds, darcy_f = pressure_drop_pipe(spec, case, props)
            k_total = None
        else:
            delta_p, velocity, reynolds, k_total = pressure_drop_minor(spec, case, props)
            darcy_f = None

        pressure_out_pa = pressure_pa - delta_p
        assert_liquid_pressure_valid(case, pressure_out_pa, f"downstream of {spec.name}")
        rows.append(
            LadderRow(
                component=spec.name,
                quantity=spec.quantity,
                pressure_in_pa=pressure_pa,
                delta_p_pa=delta_p,
                pressure_out_pa=pressure_out_pa,
                density_kg_m3=props.density_kg_m3,
                velocity_m_s=velocity,
                reynolds=reynolds,
                darcy_f=darcy_f,
                k_total=k_total,
            )
        )
        pressure_pa = pressure_out_pa

    return rows


def format_ladder(rows: Iterable[LadderRow]) -> str:
    lines = [
        "Component                         Qty   P_in(bar)   dP(bar)   P_out(bar)   dP(psi)   v(m/s)      Re        f/K",
        "-" * 111,
    ]

    total_dp_pa = 0.0
    final_pressure_pa = None

    for row in rows:
        total_dp_pa += row.delta_p_pa
        final_pressure_pa = row.pressure_out_pa
        fk_text = f"f={row.darcy_f:.5f}" if row.darcy_f is not None else f"K={row.k_total:.4f}"
        lines.append(
            f"{row.component:<33} {row.quantity:>3d} "
            f"{row.pressure_in_pa * PA_TO_BAR:>10.3f} "
            f"{row.delta_p_pa * PA_TO_BAR:>9.4f} "
            f"{row.pressure_out_pa * PA_TO_BAR:>11.3f} "
            f"{row.delta_p_pa * PA_TO_PSI:>9.3f} "
            f"{row.velocity_m_s:>8.3f} "
            f"{row.reynolds:>9.2e} "
            f"{fk_text:>10}"
        )

    if final_pressure_pa is not None:
        lines.append("-" * 111)
        lines.append(
            f"Total dP: {total_dp_pa * PA_TO_BAR:.4f} bar "
            f"({total_dp_pa * PA_TO_PSI:.3f} psi)"
        )
        lines.append(f"Final pressure: {final_pressure_pa * PA_TO_BAR:.3f} bar")

    return "\n".join(lines)


def default_lox_case() -> tuple[FlowCase, List[ComponentSpec]]:
    case = FlowCase(
        fluid="LOX",
        temperature_k=90.0,
        inlet_pressure_pa=60.0 * BAR_TO_PA,
        mass_flow_kg_s=2.0,
        pipe_inner_diameter_m=2.0 * IN_TO_M,
        nominal_size_in=2.0,
        roughness_m=15e-6,
    )

    components = [
        ComponentSpec("Sharp tank outlet entrance", "pipe_entrance_sharp"),
        ComponentSpec("Main 2 inch LOX run", "pipe", length_m=10.0),
        ComponentSpec("Standard 90 deg elbow", "standard_elbow_90", quantity=2),
        ComponentSpec("Full-port ball valve", "ball_valve_full_port"),
        ComponentSpec("Pipe exit / injector interface", "pipe_exit"),
    ]
    return case, components


def launch_gui() -> None:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    flow_case, component_specs = default_lox_case()
    fluid_label_to_key = {
        str(entry["label"]): key
        for key, entry in FLUID_CATALOG.items()
    }
    fluid_labels = list(fluid_label_to_key)
    label_to_key = {
        str(entry["label"]): key
        for key, entry in COMPONENT_CATALOG.items()
    }
    catalog_labels = list(label_to_key)
    crane_size_labels = list(CRANE_SIZE_VALUE_BY_LABEL)
    roughness_labels = list(ROUGHNESS_UM_BY_LABEL)

    root = tk.Tk()
    root.title("PFS Pressure Ladder")
    root.geometry("1180x760")
    root.minsize(1040, 680)

    style = ttk.Style(root)
    if "vista" in style.theme_names():
        style.theme_use("vista")

    root.columnconfigure(0, weight=1)
    root.rowconfigure(2, weight=1)
    root.rowconfigure(4, weight=2)

    def number_var(value: float) -> tk.StringVar:
        return tk.StringVar(value=f"{value:g}")

    pressure_unit = tk.StringVar(value="bar")
    mass_flow_unit = tk.StringVar(value="kg/s")
    temperature_unit = tk.StringVar(value="K")
    length_unit = tk.StringVar(value="m")

    def pressure_to_pa(value: float) -> float:
        return value * BAR_TO_PA if pressure_unit.get() == "bar" else value * PSI_TO_PA

    def pressure_from_pa(value_pa: float) -> float:
        return value_pa * PA_TO_BAR if pressure_unit.get() == "bar" else value_pa * PA_TO_PSI

    def mass_flow_to_kg_s(value: float) -> float:
        return value if mass_flow_unit.get() == "kg/s" else value * LBM_TO_KG

    def mass_flow_from_kg_s(value_kg_s: float) -> float:
        return value_kg_s if mass_flow_unit.get() == "kg/s" else value_kg_s / LBM_TO_KG

    def temperature_to_k(value: float) -> float:
        return value if temperature_unit.get() == "K" else value + 273.15

    def temperature_from_k(value_k: float) -> float:
        return value_k if temperature_unit.get() == "K" else value_k - 273.15

    def length_to_m(value: float) -> float:
        unit = length_unit.get()
        if unit == "m":
            return value
        if unit == "ft":
            return value * FT_TO_M
        return value * IN_TO_M

    def length_from_m(value_m: float) -> float:
        unit = length_unit.get()
        if unit == "m":
            return value_m
        if unit == "ft":
            return value_m / FT_TO_M
        return value_m / IN_TO_M

    selected_fluid_label = tk.StringVar(value=str(FLUID_CATALOG[flow_case.fluid]["label"]))
    inlet_pressure = number_var(pressure_from_pa(flow_case.inlet_pressure_pa))
    mass_flow = number_var(mass_flow_from_kg_s(flow_case.mass_flow_kg_s))
    temperature = number_var(temperature_from_k(flow_case.temperature_k))
    default_size_label = next(
        label
        for label, size in CRANE_SIZE_VALUE_BY_LABEL.items()
        if size == flow_case.nominal_size_in
    )
    selected_pipe_size_label = tk.StringVar(value=default_size_label)
    default_roughness_um = flow_case.roughness_m * 1e6
    default_roughness_label = min(
        roughness_labels,
        key=lambda label: abs(ROUGHNESS_UM_BY_LABEL[label] - default_roughness_um),
    )
    selected_roughness_label = tk.StringVar(value=default_roughness_label)

    selected_component_label = tk.StringVar(value=COMPONENT_CATALOG["pipe"]["label"])
    component_name = tk.StringVar(value=COMPONENT_CATALOG["pipe"]["label"])
    component_quantity = tk.StringVar(value="1")
    component_length = tk.StringVar(value="10")
    length_label_widget: Optional[ttk.Label] = None
    length_entry_widget: Optional[ttk.Entry] = None
    tee_flow_pattern = tk.StringVar(value="joining: branch + straight -> one outlet")
    tee_loss_path_display = tk.StringVar(value="side branch path")
    tee_angle_deg = tk.StringVar(value="90")
    tee_q_ratio = tk.StringVar(value="0.5")
    tee_beta_branch = tk.StringVar(value="1.0")
    transition_angle_deg = tk.StringVar(value="45")
    transition_beta = tk.StringVar(value="0.5")
    custom_k_value = tk.StringVar(value="1.0")
    custom_cv_value = tk.StringVar(value="1.0")
    custom_cda_mm2 = tk.StringVar(value="10.0")
    bend_count = tk.StringVar(value="2")
    bend_r_over_d = tk.StringVar(value="1.5")
    reduced_formula = tk.StringVar(value="6")
    reduced_k1 = tk.StringVar(value="1.0")
    reduced_beta = tk.StringVar(value="0.8")
    reduced_angle_deg = tk.StringVar(value="90")
    summary_text = tk.StringVar(value="Ready.")
    diagnostics_text = tk.StringVar(value="Run a calculation to see diagnostics.")

    tee_flow_pattern_to_mode = {
        "joining: branch + straight -> one outlet": "converging",
        "splitting: one inlet -> straight + branch": "diverging",
    }
    tee_loss_path_to_key = {
        "side branch path": "branch",
        "straight-through run path": "run",
    }

    def parse_float(var: tk.StringVar, label: str, positive: bool = True) -> float:
        try:
            value = float(var.get())
        except ValueError as exc:
            raise ValueError(f"{label} must be a number.") from exc
        if positive and value <= 0:
            raise ValueError(f"{label} must be greater than zero.")
        return value

    def parse_int(var: tk.StringVar, label: str) -> int:
        try:
            value = int(var.get())
        except ValueError as exc:
            raise ValueError(f"{label} must be an integer.") from exc
        if value <= 0:
            raise ValueError(f"{label} must be greater than zero.")
        return value

    def make_case_from_inputs() -> FlowCase:
        selected_size_in = CRANE_SIZE_VALUE_BY_LABEL[selected_pipe_size_label.get()]
        return FlowCase(
            fluid=fluid_label_to_key[selected_fluid_label.get()],
            temperature_k=temperature_to_k(parse_float(temperature, "Temperature", positive=False)),
            inlet_pressure_pa=pressure_to_pa(parse_float(inlet_pressure, "Inlet pressure")),
            mass_flow_kg_s=mass_flow_to_kg_s(parse_float(mass_flow, "Mass flow")),
            pipe_inner_diameter_m=selected_size_in * IN_TO_M,
            nominal_size_in=selected_size_in,
            roughness_m=ROUGHNESS_UM_BY_LABEL[selected_roughness_label.get()] * 1e-6,
        )

    def refresh_component_tree() -> None:
        component_tree.delete(*component_tree.get_children())
        component_tree.heading("length", text=f"Length {length_unit.get()}")
        for index, spec in enumerate(component_specs, start=1):
            label = COMPONENT_CATALOG[spec.catalog_key]["label"]
            length_text = "" if spec.length_m is None else f"{length_from_m(spec.length_m):g}"
            details_text = ""
            if COMPONENT_CATALOG[spec.catalog_key]["kind"] == "tee_wye":
                flow_text = "joining" if spec.tee_flow_mode == "converging" else "splitting"
                path_text = "side branch" if spec.tee_loss_path == "branch" else "straight run"
                details_text = (
                    f"{flow_text}, {path_text}, "
                    f"a={spec.tee_angle_deg:g}, q={spec.tee_q_branch_over_q_comb:g}, "
                    f"beta={spec.tee_beta_branch:g}"
                )
            elif COMPONENT_CATALOG[spec.catalog_key]["kind"] == "transition":
                details_text = f"angle={spec.transition_angle_deg:g}, beta={spec.transition_beta:g}"
            elif COMPONENT_CATALOG[spec.catalog_key]["kind"] == "custom_k":
                details_text = f"K={spec.custom_k:g}"
            elif COMPONENT_CATALOG[spec.catalog_key]["kind"] == "custom_cv":
                details_text = f"Cv={spec.custom_cv:g}"
            elif COMPONENT_CATALOG[spec.catalog_key]["kind"] == "custom_cda":
                details_text = f"CdA={spec.custom_cda_mm2:g} mm2"
            elif COMPONENT_CATALOG[spec.catalog_key]["kind"] == "multiple_90_bends":
                details_text = f"n={spec.bend_count}, r/d={spec.bend_r_over_d:g}"
            elif COMPONENT_CATALOG[spec.catalog_key]["kind"] == "reduced_port":
                details_text = (
                    f"Formula {spec.reduced_port_formula}, K1={spec.reduced_port_k1:g}, "
                    f"beta={spec.reduced_port_beta:g}, angle={spec.reduced_port_angle_deg:g}"
                )
            component_tree.insert(
                "",
                "end",
                iid=str(index - 1),
                values=(index, spec.name, label, spec.quantity, length_text, details_text),
            )

    def build_diagnostics(rows: List[LadderRow], case: FlowCase) -> str:
        if not rows:
            return "Status: No results\n\nAdd components and calculate."

        total_dp_pa = sum(row.delta_p_pa for row in rows)
        final_pressure_pa = rows[-1].pressure_out_pa
        inlet_pressure_pa = rows[0].pressure_in_pa
        largest_loss = max(rows, key=lambda row: row.delta_p_pa)
        max_velocity = max(rows, key=lambda row: row.velocity_m_s)
        warnings: List[str] = []
        entry = FLUID_CATALOG[case.fluid]

        total_dp_fraction = total_dp_pa / inlet_pressure_pa if inlet_pressure_pa > 0 else 0.0
        if total_dp_fraction > 0.20:
            warnings.append("High total dP fraction")

        phase = str(entry.get("phase", "unknown"))
        if phase == "liquid":
            if max_velocity.velocity_m_s > 20.0:
                warnings.append("Very high liquid velocity")
            elif max_velocity.velocity_m_s > 10.0:
                warnings.append("High liquid velocity")

        saturation_pressure_pa = liquid_saturation_pressure_pa(case)
        psat_lines: List[str] = []
        if saturation_pressure_pa is not None:
            margins = [row.pressure_out_pa - saturation_pressure_pa for row in rows]
            min_margin_pa = min(margins)
            psat_lines.append(f"Psat: {pressure_from_pa(saturation_pressure_pa):.3f} {pressure_unit.get()}")
            psat_lines.append(f"Min Psat margin: {pressure_from_pa(min_margin_pa):.3f} {pressure_unit.get()}")
            if min_margin_pa < 0.0:
                warnings.append("Below saturation")
            elif min_margin_pa < 5.0 * BAR_TO_PA:
                warnings.append("Near saturation")

        if phase == "gas":
            warnings.append("Gas dP is incompressible approximation")

        status = "OK" if not warnings else "Warning"
        lines = [
            f"Status: {status}",
            f"Fluid: {entry['label']}",
            f"Total dP: {pressure_from_pa(total_dp_pa):.4f} {pressure_unit.get()}",
            f"Final P: {pressure_from_pa(final_pressure_pa):.3f} {pressure_unit.get()}",
            f"dP / inlet P: {100.0 * total_dp_fraction:.1f}%",
            "",
            f"Largest loss: {largest_loss.component}",
            f"  {pressure_from_pa(largest_loss.delta_p_pa):.4f} {pressure_unit.get()}",
            f"Max velocity: {max_velocity.velocity_m_s:.2f} m/s",
            f"  at {max_velocity.component}",
        ]
        if psat_lines:
            lines.extend(["", *psat_lines])
        if warnings:
            lines.extend(["", "Warnings:", *[f"- {warning}" for warning in warnings]])
        else:
            lines.extend(["", "Warnings: none"])
        return "\n".join(lines)

    def refresh_results_tree(rows: Iterable[LadderRow], case: FlowCase) -> None:
        result_tree.delete(*result_tree.get_children())
        headings["p_in_bar"] = f"P in {pressure_unit.get()}"
        headings["dp_bar"] = f"dP {pressure_unit.get()}"
        headings["p_out_bar"] = f"P out {pressure_unit.get()}"
        headings["velocity"] = "v m/s"
        for column, heading in headings.items():
            result_tree.heading(column, text=heading)
        row_list = list(rows)
        total_dp_pa = 0.0
        final_pressure_pa = None

        for row in row_list:
            total_dp_pa += row.delta_p_pa
            final_pressure_pa = row.pressure_out_pa
            fk_text = f"f={row.darcy_f:.5f}" if row.darcy_f is not None else f"K={row.k_total:.4f}"
            result_tree.insert(
                "",
                "end",
                values=(
                    row.component,
                    row.quantity,
                    f"{pressure_from_pa(row.pressure_in_pa):.3f}",
                    f"{pressure_from_pa(row.delta_p_pa):.5f}",
                    f"{pressure_from_pa(row.pressure_out_pa):.3f}",
                    f"{row.delta_p_pa * PA_TO_PSI:.4f}",
                    f"{row.velocity_m_s:.3f}",
                    f"{row.reynolds:.3e}",
                    fk_text,
                ),
            )

        if final_pressure_pa is None:
            summary_text.set("No components in ladder.")
        else:
            summary_text.set(
                f"Total dP = {pressure_from_pa(total_dp_pa):.5f} {pressure_unit.get()} "
                f"({total_dp_pa * PA_TO_PSI:.4f} psi)    "
                f"Final pressure = {pressure_from_pa(final_pressure_pa):.3f} {pressure_unit.get()}"
            )
        diagnostics_text.set(build_diagnostics(row_list, case))

    def copy_result_rows(selected_only: bool = True) -> None:
        item_ids = result_tree.selection() if selected_only else result_tree.get_children()
        if selected_only and not item_ids:
            item_ids = result_tree.get_children()
        if not item_ids:
            return

        header = [headings[column] for column in result_columns]
        lines = ["\t".join(header)]
        for item_id in item_ids:
            values = result_tree.item(item_id, "values")
            lines.append("\t".join(str(value) for value in values))
        root.clipboard_clear()
        root.clipboard_append("\n".join(lines))
        summary_text.set(f"Copied {len(item_ids)} result row(s) to clipboard.")

    def component_from_form() -> ComponentSpec:
        label = selected_component_label.get()
        key = label_to_key[label]
        kind = COMPONENT_CATALOG[key]["kind"]
        quantity = parse_int(component_quantity, "Quantity")
        name = component_name.get().strip() or label
        length_m = None
        tee_kwargs: Dict[str, object] = {}
        transition_kwargs: Dict[str, object] = {}
        custom_kwargs: Dict[str, object] = {}

        if kind == "pipe":
            length_m = length_to_m(parse_float(component_length, "Pipe length"))
        elif kind == "tee_wye":
            catalog_entry = COMPONENT_CATALOG[key]
            angle = (
                float(catalog_entry["fixed_angle_deg"])
                if "fixed_angle_deg" in catalog_entry
                else parse_float(tee_angle_deg, "Tee/wye angle")
            )
            tee_kwargs = {
                "tee_flow_mode": tee_flow_pattern_to_mode[tee_flow_pattern.get()],
                "tee_loss_path": tee_loss_path_to_key[tee_loss_path_display.get()],
                "tee_angle_deg": angle,
                "tee_q_branch_over_q_comb": parse_float(tee_q_ratio, "Branch flow fraction"),
                "tee_beta_branch": parse_float(tee_beta_branch, "Branch diameter ratio"),
            }
        elif kind == "transition":
            transition_kwargs = {
                "transition_angle_deg": parse_float(transition_angle_deg, "Transition angle"),
                "transition_beta": parse_float(transition_beta, "Transition diameter ratio"),
            }
        elif kind == "custom_k":
            custom_kwargs = {"custom_k": parse_float(custom_k_value, "Custom K", positive=False)}
        elif kind == "custom_cv":
            custom_kwargs = {"custom_cv": parse_float(custom_cv_value, "Custom Cv")}
        elif kind == "custom_cda":
            custom_kwargs = {"custom_cda_mm2": parse_float(custom_cda_mm2, "Custom CdA")}
        elif kind == "multiple_90_bends":
            custom_kwargs = {
                "bend_count": parse_int(bend_count, "Bend count"),
                "bend_r_over_d": parse_float(bend_r_over_d, "Bend r/d"),
            }
        elif kind == "reduced_port":
            custom_kwargs = {
                "reduced_port_formula": parse_int(reduced_formula, "Reduced-port formula"),
                "reduced_port_k1": parse_float(reduced_k1, "Reduced-port K1", positive=False),
                "reduced_port_beta": parse_float(reduced_beta, "Reduced-port beta"),
                "reduced_port_angle_deg": parse_float(reduced_angle_deg, "Reduced-port angle", positive=False),
            }

        return ComponentSpec(
            name=name,
            catalog_key=key,
            quantity=quantity,
            length_m=length_m,
            **tee_kwargs,
            **transition_kwargs,
            **custom_kwargs,
        )

    def load_component_into_form(spec: ComponentSpec) -> None:
        catalog_entry = COMPONENT_CATALOG[spec.catalog_key]
        selected_component_label.set(str(catalog_entry["label"]))
        component_name.set(spec.name)
        component_quantity.set(str(spec.quantity))
        component_length.set("" if spec.length_m is None else f"{length_from_m(spec.length_m):g}")
        if spec.tee_flow_mode:
            for display, mode in tee_flow_pattern_to_mode.items():
                if mode == spec.tee_flow_mode:
                    tee_flow_pattern.set(display)
        if spec.tee_loss_path:
            for display, path in tee_loss_path_to_key.items():
                if path == spec.tee_loss_path:
                    tee_loss_path_display.set(display)
        if spec.tee_angle_deg is not None:
            tee_angle_deg.set(f"{spec.tee_angle_deg:g}")
        if spec.tee_q_branch_over_q_comb is not None:
            tee_q_ratio.set(f"{spec.tee_q_branch_over_q_comb:g}")
        if spec.tee_beta_branch is not None:
            tee_beta_branch.set(f"{spec.tee_beta_branch:g}")
        if spec.transition_angle_deg is not None:
            transition_angle_deg.set(f"{spec.transition_angle_deg:g}")
        if spec.transition_beta is not None:
            transition_beta.set(f"{spec.transition_beta:g}")
        if spec.custom_k is not None:
            custom_k_value.set(f"{spec.custom_k:g}")
        if spec.custom_cv is not None:
            custom_cv_value.set(f"{spec.custom_cv:g}")
        if spec.custom_cda_mm2 is not None:
            custom_cda_mm2.set(f"{spec.custom_cda_mm2:g}")
        if spec.bend_count is not None:
            bend_count.set(str(spec.bend_count))
        if spec.bend_r_over_d is not None:
            bend_r_over_d.set(f"{spec.bend_r_over_d:g}")
        if spec.reduced_port_formula is not None:
            reduced_formula.set(str(spec.reduced_port_formula))
        if spec.reduced_port_k1 is not None:
            reduced_k1.set(f"{spec.reduced_port_k1:g}")
        if spec.reduced_port_beta is not None:
            reduced_beta.set(f"{spec.reduced_port_beta:g}")
        if spec.reduced_port_angle_deg is not None:
            reduced_angle_deg.set(f"{spec.reduced_port_angle_deg:g}")
        update_component_option_visibility()

    def selected_component_index() -> Optional[int]:
        selection = component_tree.selection()
        if not selection:
            return None
        return int(selection[0])

    drag_start_index: Optional[int] = None

    def on_component_drag_start(event: tk.Event) -> None:
        nonlocal drag_start_index
        row_id = component_tree.identify_row(event.y)
        if not row_id:
            drag_start_index = None
            return
        drag_start_index = int(row_id)
        component_tree.selection_set(row_id)

    def on_component_drag_release(event: tk.Event) -> None:
        nonlocal drag_start_index
        if drag_start_index is None:
            return
        target_id = component_tree.identify_row(event.y)
        if not target_id:
            drag_start_index = None
            return

        target_index = int(target_id)
        if target_index == drag_start_index:
            drag_start_index = None
            return

        moved_spec = component_specs.pop(drag_start_index)
        if target_index > drag_start_index:
            target_index -= 1
        component_specs.insert(target_index, moved_spec)
        refresh_component_tree()
        component_tree.selection_set(str(target_index))
        drag_start_index = None

    def apply_fluid_defaults(*_: object) -> None:
        entry = FLUID_CATALOG[fluid_label_to_key[selected_fluid_label.get()]]
        temperature.set(f"{temperature_from_k(float(entry['default_temperature_k'])):g}")
        inlet_pressure.set(f"{pressure_from_pa(float(entry['default_pressure_bar']) * BAR_TO_PA):g}")

    def update_component_option_visibility() -> None:
        key = label_to_key[selected_component_label.get()]
        catalog_entry = COMPONENT_CATALOG[key]
        kind = catalog_entry["kind"]
        option_frames = (
            tee_options,
            transition_options,
            custom_k_options,
            custom_cv_options,
            custom_cda_options,
            multiple_bend_options,
            reduced_port_options,
        )
        for frame in option_frames:
            frame.grid_remove()
        if length_label_widget is not None and length_entry_widget is not None:
            if kind == "pipe":
                length_label_widget.configure(text=f"Pipe length ({length_unit.get()})")
                length_entry_widget.configure(state="normal")
            else:
                length_label_widget.configure(text=f"Pipe length ({length_unit.get()})")
                length_entry_widget.configure(state="disabled")
                component_length.set("")
        if kind == "tee_wye":
            tee_options.grid()
            if "fixed_angle_deg" in catalog_entry:
                fixed_angle = float(catalog_entry["fixed_angle_deg"])
                tee_angle_deg.set(f"{fixed_angle:g}")
                tee_angle_menu.configure(values=(f"{fixed_angle:g}",), state="disabled")
            else:
                angle_values = tuple(f"{float(angle):g}" for angle in catalog_entry.get("angle_options", (30.0, 45.0, 60.0)))
                tee_angle_menu.configure(values=angle_values, state="readonly")
                if tee_angle_deg.get() not in angle_values:
                    tee_angle_deg.set(angle_values[0])
        elif kind == "transition":
            transition_options.grid()
        elif kind == "custom_k":
            custom_k_options.grid()
        elif kind == "custom_cv":
            custom_cv_options.grid()
        elif kind == "custom_cda":
            custom_cda_options.grid()
        elif kind == "multiple_90_bends":
            multiple_bend_options.grid()
        elif kind == "reduced_port":
            reduced_port_options.grid()

    def on_catalog_selected(*_: object) -> None:
        label = selected_component_label.get()
        component_name.set(label)
        key = label_to_key[label]
        if COMPONENT_CATALOG[key]["kind"] == "pipe":
            component_length.set("10")
        else:
            component_length.set("")
        update_component_option_visibility()

    def add_component() -> None:
        try:
            component_specs.append(component_from_form())
            refresh_component_tree()
        except ValueError as exc:
            messagebox.showerror("Input error", str(exc))

    def edit_selected_component() -> None:
        index = selected_component_index()
        if index is None:
            return
        load_component_into_form(component_specs[index])
        summary_text.set("Loaded selected component into editor.")

    def update_selected_component() -> None:
        index = selected_component_index()
        if index is None:
            return
        try:
            component_specs[index] = component_from_form()
            refresh_component_tree()
            component_tree.selection_set(str(index))
            summary_text.set("Updated selected component.")
        except ValueError as exc:
            messagebox.showerror("Input error", str(exc))

    def remove_selected_component() -> None:
        index = selected_component_index()
        if index is None:
            return
        component_specs.pop(index)
        refresh_component_tree()

    def move_selected(delta: int) -> None:
        index = selected_component_index()
        if index is None:
            return
        new_index = index + delta
        if not 0 <= new_index < len(component_specs):
            return
        component_specs[index], component_specs[new_index] = component_specs[new_index], component_specs[index]
        refresh_component_tree()
        component_tree.selection_set(str(new_index))

    def load_default_components() -> None:
        _, defaults = default_lox_case()
        component_specs.clear()
        component_specs.extend(defaults)
        refresh_component_tree()
        result_tree.delete(*result_tree.get_children())
        summary_text.set("Default LOX ladder loaded.")
        diagnostics_text.set("Default ladder loaded. Run a calculation to refresh diagnostics.")

    def clear_components() -> None:
        component_specs.clear()
        refresh_component_tree()
        result_tree.delete(*result_tree.get_children())
        summary_text.set("Component ladder cleared.")
        diagnostics_text.set("No results.")

    def save_ladder() -> None:
        try:
            case = make_case_from_inputs()
        except ValueError as exc:
            messagebox.showerror("Input error", str(exc))
            return
        path = filedialog.asksaveasfilename(
            title="Save pressure ladder",
            defaultextension=".json",
            filetypes=(("JSON files", "*.json"), ("All files", "*.*")),
        )
        if not path:
            return
        data = {
            "version": 1,
            "flow_case": asdict(case),
            "components": [asdict(spec) for spec in component_specs],
            "ui_units": {
                "pressure": pressure_unit.get(),
                "mass_flow": mass_flow_unit.get(),
                "temperature": temperature_unit.get(),
                "length": length_unit.get(),
                "roughness_label": selected_roughness_label.get(),
                "pipe_size_label": selected_pipe_size_label.get(),
            },
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2)
        summary_text.set(f"Saved ladder to {path}")

    def load_ladder() -> None:
        path = filedialog.askopenfilename(
            title="Load pressure ladder",
            filetypes=(("JSON files", "*.json"), ("All files", "*.*")),
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as handle:
                data = json.load(handle)
            loaded_case = FlowCase(**data["flow_case"])
            loaded_components = [ComponentSpec(**item) for item in data["components"]]
        except Exception as exc:
            messagebox.showerror("Load error", str(exc))
            return

        ui_units = data.get("ui_units", {})
        pressure_unit.set(ui_units.get("pressure", pressure_unit.get()))
        mass_flow_unit.set(ui_units.get("mass_flow", mass_flow_unit.get()))
        temperature_unit.set(ui_units.get("temperature", temperature_unit.get()))
        length_unit.set(ui_units.get("length", length_unit.get()))
        selected_fluid_label.set(str(FLUID_CATALOG[loaded_case.fluid]["label"]))
        inlet_pressure.set(f"{pressure_from_pa(loaded_case.inlet_pressure_pa):g}")
        mass_flow.set(f"{mass_flow_from_kg_s(loaded_case.mass_flow_kg_s):g}")
        temperature.set(f"{temperature_from_k(loaded_case.temperature_k):g}")
        selected_pipe_size_label.set(
            next(label for label, size in CRANE_SIZE_VALUE_BY_LABEL.items() if size == loaded_case.nominal_size_in)
        )
        roughness_um = loaded_case.roughness_m * 1e6
        selected_roughness_label.set(ui_units.get(
            "roughness_label",
            min(roughness_labels, key=lambda label: abs(ROUGHNESS_UM_BY_LABEL[label] - roughness_um)),
        ))
        component_specs.clear()
        component_specs.extend(loaded_components)
        refresh_component_tree()
        result_tree.delete(*result_tree.get_children())
        diagnostics_text.set("Loaded ladder. Run a calculation to refresh diagnostics.")
        summary_text.set(f"Loaded ladder from {path}")

    def export_results_csv() -> None:
        item_ids = result_tree.get_children()
        if not item_ids:
            messagebox.showinfo("Export CSV", "No result rows to export.")
            return
        path = filedialog.asksaveasfilename(
            title="Export pressure ladder results",
            defaultextension=".csv",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
        )
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow([headings[column] for column in result_columns])
            for item_id in item_ids:
                writer.writerow(result_tree.item(item_id, "values"))
        summary_text.set(f"Exported results to {path}")

    def calculate() -> None:
        try:
            case = make_case_from_inputs()
            rows = run_pressure_ladder(case, component_specs)
            refresh_results_tree(rows, case)
        except Exception as exc:  # CoolProp errors are intentionally surfaced here.
            diagnostics_text.set(f"Status: Invalid\n\n{exc}")
            summary_text.set("Calculation failed.")
            messagebox.showerror("Calculation error", str(exc))

    def on_pressure_unit_changed(*_: object) -> None:
        summary_text.set("Pressure unit changed. Existing typed value will be interpreted in the selected unit.")

    def on_mass_flow_unit_changed(*_: object) -> None:
        summary_text.set("Mass flow unit changed. Existing typed value will be interpreted in the selected unit.")

    def on_temperature_unit_changed(*_: object) -> None:
        summary_text.set("Temperature unit changed. Existing typed value will be interpreted in the selected unit.")

    def on_length_unit_changed(*_: object) -> None:
        new = length_unit.get()
        if length_label_widget is not None:
            length_label_widget.configure(text=f"Pipe length ({new})")
        refresh_component_tree()
        summary_text.set("Length unit changed. Existing typed value will be interpreted in the selected unit.")

    inputs = ttk.LabelFrame(root, text="Flow case")
    inputs.grid(row=0, column=0, sticky="ew", padx=12, pady=(12, 6))
    for column in range(14):
        inputs.columnconfigure(column, weight=1)

    input_fields = [
        ("Fluid", selected_fluid_label),
        ("P inlet", inlet_pressure),
        ("mdot", mass_flow),
        ("T", temperature),
        ("Pipe size", selected_pipe_size_label),
        ("Surface roughness", selected_roughness_label),
    ]
    for column, (label, var) in enumerate(input_fields):
        ttk.Label(inputs, text=label).grid(row=0, column=column * 2, sticky="w", padx=(8, 4), pady=8)
        if label == "Fluid":
            fluid_menu = ttk.Combobox(
                inputs,
                textvariable=selected_fluid_label,
                values=fluid_labels,
                state="readonly",
                width=24,
            )
            fluid_menu.grid(row=0, column=column * 2 + 1, sticky="ew", padx=(0, 8), pady=8)
            fluid_menu.bind("<<ComboboxSelected>>", apply_fluid_defaults)
        elif label == "Pipe size":
            size_menu = ttk.Combobox(
                inputs,
                textvariable=selected_pipe_size_label,
                values=crane_size_labels,
                state="readonly",
                width=18,
            )
            size_menu.grid(row=0, column=column * 2 + 1, sticky="ew", padx=(0, 8), pady=8)
        elif label == "Surface roughness":
            roughness_menu = ttk.Combobox(
                inputs,
                textvariable=selected_roughness_label,
                values=roughness_labels,
                state="readonly",
                width=34,
            )
            roughness_menu.grid(row=0, column=column * 2 + 1, sticky="ew", padx=(0, 8), pady=8)
        else:
            ttk.Entry(inputs, textvariable=var, width=11).grid(row=0, column=column * 2 + 1, sticky="ew", padx=(0, 8), pady=8)

    unit_fields = [
        ("Pressure unit", pressure_unit, ("bar", "psi")),
        ("Mass flow unit", mass_flow_unit, ("kg/s", "lbm/s")),
        ("Temperature unit", temperature_unit, ("K", "C")),
        ("Length unit", length_unit, ("m", "ft", "in")),
    ]
    for column, (label, var, values) in enumerate(unit_fields):
        base_column = (column + 1) * 2
        ttk.Label(inputs, text=label).grid(row=1, column=base_column, sticky="w", padx=(8, 4), pady=(0, 8))
        unit_menu = ttk.Combobox(inputs, textvariable=var, values=values, state="readonly", width=9)
        unit_menu.grid(
            row=1, column=base_column + 1, sticky="ew", padx=(0, 8), pady=(0, 8)
        )
        if label == "Pressure unit":
            unit_menu.bind("<<ComboboxSelected>>", on_pressure_unit_changed)
        elif label == "Mass flow unit":
            unit_menu.bind("<<ComboboxSelected>>", on_mass_flow_unit_changed)
        elif label == "Temperature unit":
            unit_menu.bind("<<ComboboxSelected>>", on_temperature_unit_changed)
        elif label == "Length unit":
            unit_menu.bind("<<ComboboxSelected>>", on_length_unit_changed)

    builder = ttk.LabelFrame(root, text="Add component")
    builder.grid(row=1, column=0, sticky="ew", padx=12, pady=6)
    for column in range(10):
        builder.columnconfigure(column, weight=1)

    ttk.Label(builder, text="Type").grid(row=0, column=0, sticky="w", padx=(8, 4), pady=8)
    component_menu = ttk.Combobox(
        builder,
        textvariable=selected_component_label,
        values=catalog_labels,
        state="readonly",
        width=28,
    )
    component_menu.grid(row=0, column=1, sticky="ew", padx=(0, 8), pady=8)
    component_menu.bind("<<ComboboxSelected>>", on_catalog_selected)

    ttk.Label(builder, text="Name").grid(row=0, column=2, sticky="w", padx=(8, 4), pady=8)
    ttk.Entry(builder, textvariable=component_name, width=28).grid(row=0, column=3, sticky="ew", padx=(0, 8), pady=8)
    ttk.Label(builder, text="Qty").grid(row=0, column=4, sticky="w", padx=(8, 4), pady=8)
    ttk.Entry(builder, textvariable=component_quantity, width=7).grid(row=0, column=5, sticky="ew", padx=(0, 8), pady=8)
    length_label_widget = ttk.Label(builder, text=f"Pipe length ({length_unit.get()})")
    length_label_widget.grid(row=0, column=6, sticky="w", padx=(8, 4), pady=8)
    length_entry_widget = ttk.Entry(builder, textvariable=component_length, width=9)
    length_entry_widget.grid(row=0, column=7, sticky="ew", padx=(0, 8), pady=8)
    ttk.Button(builder, text="Add", command=add_component).grid(row=0, column=8, sticky="ew", padx=(8, 4), pady=8)
    ttk.Button(builder, text="Calculate", command=calculate).grid(row=0, column=9, sticky="ew", padx=(4, 8), pady=8)

    tee_options = ttk.Frame(builder)
    tee_options.grid(row=1, column=0, columnspan=10, sticky="ew", padx=8, pady=(0, 8))
    for column in range(10):
        tee_options.columnconfigure(column, weight=1)

    ttk.Label(tee_options, text="Flow pattern").grid(row=0, column=0, sticky="w", padx=(0, 4))
    ttk.Combobox(
        tee_options,
        textvariable=tee_flow_pattern,
        values=tuple(tee_flow_pattern_to_mode),
        state="readonly",
        width=34,
    ).grid(row=0, column=1, columnspan=2, sticky="ew", padx=(0, 8))
    ttk.Label(tee_options, text="Loss path").grid(row=0, column=3, sticky="w", padx=(0, 4))
    ttk.Combobox(
        tee_options,
        textvariable=tee_loss_path_display,
        values=tuple(tee_loss_path_to_key),
        state="readonly",
        width=22,
    ).grid(row=0, column=4, sticky="ew", padx=(0, 8))
    ttk.Label(tee_options, text="Branch angle deg").grid(row=0, column=5, sticky="w", padx=(0, 4))
    tee_angle_menu = ttk.Combobox(
        tee_options,
        textvariable=tee_angle_deg,
        values=("30", "45", "60", "90"),
        state="readonly",
        width=8,
    )
    tee_angle_menu.grid(row=0, column=6, sticky="ew", padx=(0, 8))
    ttk.Label(tee_options, text="Branch flow fraction").grid(row=1, column=0, sticky="w", padx=(0, 4), pady=(6, 0))
    ttk.Entry(tee_options, textvariable=tee_q_ratio, width=8).grid(row=1, column=1, sticky="ew", padx=(0, 8), pady=(6, 0))
    ttk.Label(tee_options, text="Branch diameter ratio").grid(row=1, column=2, sticky="w", padx=(0, 4), pady=(6, 0))
    ttk.Entry(tee_options, textvariable=tee_beta_branch, width=8).grid(row=1, column=3, sticky="ew", padx=(0, 8), pady=(6, 0))
    ttk.Label(tee_options, text="flow fraction = branch flow / total flow").grid(
        row=1, column=4, columnspan=3, sticky="w", padx=(0, 8), pady=(6, 0)
    )
    ttk.Label(tee_options, text="diameter ratio = branch diameter / main diameter").grid(
        row=1, column=7, columnspan=3, sticky="w", pady=(6, 0)
    )

    transition_options = ttk.Frame(builder)
    transition_options.grid(row=2, column=0, columnspan=10, sticky="ew", padx=8, pady=(0, 8))
    for column in range(10):
        transition_options.columnconfigure(column, weight=1)
    ttk.Label(transition_options, text="Included angle deg").grid(row=0, column=0, sticky="w", padx=(0, 4))
    ttk.Entry(transition_options, textvariable=transition_angle_deg, width=8).grid(row=0, column=1, sticky="ew", padx=(0, 8))
    ttk.Label(transition_options, text="Diameter ratio").grid(row=0, column=2, sticky="w", padx=(0, 4))
    ttk.Entry(transition_options, textvariable=transition_beta, width=8).grid(row=0, column=3, sticky="ew", padx=(0, 8))
    ttk.Label(transition_options, text="diameter ratio = smaller diameter / larger diameter").grid(
        row=0, column=4, columnspan=4, sticky="w", padx=(0, 8)
    )
    ttk.Label(transition_options, text="K is based on velocity in larger pipe").grid(
        row=0, column=8, columnspan=2, sticky="w"
    )

    custom_k_options = ttk.Frame(builder)
    custom_k_options.grid(row=3, column=0, columnspan=10, sticky="ew", padx=8, pady=(0, 8))
    for column in range(10):
        custom_k_options.columnconfigure(column, weight=1)
    ttk.Label(custom_k_options, text="K").grid(row=0, column=0, sticky="w", padx=(0, 4))
    ttk.Entry(custom_k_options, textvariable=custom_k_value, width=10).grid(row=0, column=1, sticky="ew", padx=(0, 8))
    ttk.Label(custom_k_options, text="Use vendor/test/custom K directly.").grid(
        row=0, column=2, columnspan=4, sticky="w"
    )

    custom_cv_options = ttk.Frame(builder)
    custom_cv_options.grid(row=4, column=0, columnspan=10, sticky="ew", padx=8, pady=(0, 8))
    for column in range(10):
        custom_cv_options.columnconfigure(column, weight=1)
    ttk.Label(custom_cv_options, text="Cv").grid(row=0, column=0, sticky="w", padx=(0, 4))
    ttk.Entry(custom_cv_options, textvariable=custom_cv_value, width=10).grid(row=0, column=1, sticky="ew", padx=(0, 8))
    ttk.Label(custom_cv_options, text="Liquid Cv mode: dP(psi) = SG * (Q_gpm / Cv)^2.").grid(
        row=0, column=2, columnspan=6, sticky="w"
    )

    custom_cda_options = ttk.Frame(builder)
    custom_cda_options.grid(row=5, column=0, columnspan=10, sticky="ew", padx=8, pady=(0, 8))
    for column in range(10):
        custom_cda_options.columnconfigure(column, weight=1)
    ttk.Label(custom_cda_options, text="CdA mm2").grid(row=0, column=0, sticky="w", padx=(0, 4))
    ttk.Entry(custom_cda_options, textvariable=custom_cda_mm2, width=10).grid(row=0, column=1, sticky="ew", padx=(0, 8))
    ttk.Label(custom_cda_options, text="Use for orifices/injector elements when CdA is known.").grid(
        row=0, column=2, columnspan=6, sticky="w"
    )

    multiple_bend_options = ttk.Frame(builder)
    multiple_bend_options.grid(row=6, column=0, columnspan=10, sticky="ew", padx=8, pady=(0, 8))
    for column in range(10):
        multiple_bend_options.columnconfigure(column, weight=1)
    ttk.Label(multiple_bend_options, text="Bend count").grid(row=0, column=0, sticky="w", padx=(0, 4))
    ttk.Entry(multiple_bend_options, textvariable=bend_count, width=10).grid(row=0, column=1, sticky="ew", padx=(0, 8))
    ttk.Label(multiple_bend_options, text="r/d").grid(row=0, column=2, sticky="w", padx=(0, 4))
    ttk.Combobox(
        multiple_bend_options,
        textvariable=bend_r_over_d,
        values=tuple(f"{value:g}" for value in PIPE_BEND_90_MULTIPLIER_BY_R_OVER_D),
        state="readonly",
        width=10,
    ).grid(row=0, column=3, sticky="ew", padx=(0, 8))
    ttk.Label(multiple_bend_options, text="Crane A-30 multiple close-spaced 90 deg bend equation.").grid(
        row=0, column=4, columnspan=5, sticky="w"
    )

    reduced_port_options = ttk.Frame(builder)
    reduced_port_options.grid(row=7, column=0, columnspan=10, sticky="ew", padx=8, pady=(0, 8))
    for column in range(10):
        reduced_port_options.columnconfigure(column, weight=1)
    ttk.Label(reduced_port_options, text="Formula").grid(row=0, column=0, sticky="w", padx=(0, 4))
    ttk.Combobox(
        reduced_port_options,
        textvariable=reduced_formula,
        values=("1", "2", "3", "4", "5", "6", "7"),
        state="readonly",
        width=8,
    ).grid(row=0, column=1, sticky="ew", padx=(0, 8))
    ttk.Label(reduced_port_options, text="K1").grid(row=0, column=2, sticky="w", padx=(0, 4))
    ttk.Entry(reduced_port_options, textvariable=reduced_k1, width=8).grid(row=0, column=3, sticky="ew", padx=(0, 8))
    ttk.Label(reduced_port_options, text="beta").grid(row=0, column=4, sticky="w", padx=(0, 4))
    ttk.Entry(reduced_port_options, textvariable=reduced_beta, width=8).grid(row=0, column=5, sticky="ew", padx=(0, 8))
    ttk.Label(reduced_port_options, text="angle deg").grid(row=0, column=6, sticky="w", padx=(0, 4))
    ttk.Entry(reduced_port_options, textvariable=reduced_angle_deg, width=8).grid(row=0, column=7, sticky="ew", padx=(0, 8))
    ttk.Label(reduced_port_options, text="Crane A-27 formula number; beta=d_small/d_large.").grid(
        row=0, column=8, columnspan=2, sticky="w"
    )
    update_component_option_visibility()

    components_frame = ttk.LabelFrame(root, text="Component ladder")
    components_frame.grid(row=2, column=0, sticky="nsew", padx=12, pady=6)
    components_frame.columnconfigure(0, weight=1)
    components_frame.rowconfigure(0, weight=1)

    component_columns = ("#", "name", "type", "qty", "length", "details")
    component_tree = ttk.Treeview(components_frame, columns=component_columns, show="headings", height=8)
    component_tree.heading("#", text="#")
    component_tree.heading("name", text="Name")
    component_tree.heading("type", text="Type")
    component_tree.heading("qty", text="Qty")
    component_tree.heading("length", text=f"Length {length_unit.get()}")
    component_tree.heading("details", text="Details")
    component_tree.column("#", width=44, anchor="center", stretch=False)
    component_tree.column("name", width=220)
    component_tree.column("type", width=230)
    component_tree.column("qty", width=70, anchor="center", stretch=False)
    component_tree.column("length", width=90, anchor="center", stretch=False)
    component_tree.column("details", width=260)
    component_tree.grid(row=0, column=0, sticky="nsew")
    component_tree.bind("<ButtonPress-1>", on_component_drag_start)
    component_tree.bind("<ButtonRelease-1>", on_component_drag_release)
    component_tree.bind("<Double-1>", lambda _event: edit_selected_component())

    component_scroll = ttk.Scrollbar(components_frame, orient="vertical", command=component_tree.yview)
    component_tree.configure(yscrollcommand=component_scroll.set)
    component_scroll.grid(row=0, column=1, sticky="ns")

    component_buttons = ttk.Frame(root)
    component_buttons.grid(row=3, column=0, sticky="ew", padx=12, pady=(0, 6))
    ttk.Button(component_buttons, text="Remove Selected", command=remove_selected_component).pack(side="left", padx=(0, 6))
    ttk.Button(component_buttons, text="Edit Selected", command=edit_selected_component).pack(side="left", padx=6)
    ttk.Button(component_buttons, text="Update Selected", command=update_selected_component).pack(side="left", padx=6)
    ttk.Button(component_buttons, text="Move Up", command=lambda: move_selected(-1)).pack(side="left", padx=6)
    ttk.Button(component_buttons, text="Move Down", command=lambda: move_selected(1)).pack(side="left", padx=6)
    ttk.Button(component_buttons, text="Load Default", command=load_default_components).pack(side="left", padx=6)
    ttk.Button(component_buttons, text="Save JSON", command=save_ladder).pack(side="left", padx=6)
    ttk.Button(component_buttons, text="Load JSON", command=load_ladder).pack(side="left", padx=6)
    ttk.Button(component_buttons, text="Export CSV", command=export_results_csv).pack(side="left", padx=6)
    ttk.Button(component_buttons, text="Clear", command=clear_components).pack(side="left", padx=6)

    results_frame = ttk.LabelFrame(root, text="Pressure ladder results")
    results_frame.grid(row=4, column=0, sticky="nsew", padx=12, pady=6)
    results_frame.columnconfigure(0, weight=1)
    results_frame.columnconfigure(2, weight=0)
    results_frame.rowconfigure(0, weight=1)

    result_columns = (
        "component",
        "qty",
        "p_in_bar",
        "dp_bar",
        "p_out_bar",
        "dp_psi",
        "velocity",
        "reynolds",
        "fk",
    )
    result_tree = ttk.Treeview(results_frame, columns=result_columns, show="headings", height=10, selectmode="extended")
    headings = {
        "component": "Component",
        "qty": "Qty",
        "p_in_bar": f"P in {pressure_unit.get()}",
        "dp_bar": f"dP {pressure_unit.get()}",
        "p_out_bar": f"P out {pressure_unit.get()}",
        "dp_psi": "dP psi",
        "velocity": "v m/s",
        "reynolds": "Re",
        "fk": "f / K",
    }
    for column, heading in headings.items():
        result_tree.heading(column, text=heading)
    result_tree.column("component", width=250)
    result_tree.column("qty", width=60, anchor="center", stretch=False)
    result_tree.column("p_in_bar", width=90, anchor="e")
    result_tree.column("dp_bar", width=90, anchor="e")
    result_tree.column("p_out_bar", width=90, anchor="e")
    result_tree.column("dp_psi", width=90, anchor="e")
    result_tree.column("velocity", width=80, anchor="e")
    result_tree.column("reynolds", width=95, anchor="e")
    result_tree.column("fk", width=90, anchor="e")
    result_tree.grid(row=0, column=0, sticky="nsew")
    result_tree.bind("<Control-c>", lambda _event: copy_result_rows(selected_only=True))
    result_tree.bind("<Control-C>", lambda _event: copy_result_rows(selected_only=True))

    result_scroll = ttk.Scrollbar(results_frame, orient="vertical", command=result_tree.yview)
    result_tree.configure(yscrollcommand=result_scroll.set)
    result_scroll.grid(row=0, column=1, sticky="ns")

    diagnostics_frame = ttk.LabelFrame(results_frame, text="Diagnostics")
    diagnostics_frame.grid(row=0, column=2, sticky="nsew", padx=(10, 0))
    diagnostics_frame.columnconfigure(0, weight=1)
    diagnostics_frame.rowconfigure(0, weight=1)
    ttk.Label(
        diagnostics_frame,
        textvariable=diagnostics_text,
        anchor="nw",
        justify="left",
        width=34,
    ).grid(row=0, column=0, sticky="nsew", padx=8, pady=8)
    ttk.Button(
        diagnostics_frame,
        text="Copy Selected Rows",
        command=lambda: copy_result_rows(selected_only=True),
    ).grid(row=1, column=0, sticky="ew", padx=8, pady=(0, 6))
    ttk.Button(
        diagnostics_frame,
        text="Copy All Rows",
        command=lambda: copy_result_rows(selected_only=False),
    ).grid(row=2, column=0, sticky="ew", padx=8, pady=(0, 8))

    ttk.Label(root, textvariable=summary_text, anchor="w").grid(row=5, column=0, sticky="ew", padx=12, pady=(2, 10))

    refresh_component_tree()
    calculate()
    root.mainloop()


if __name__ == "__main__":
    flow_case, component_list = default_lox_case()
    if "--cli" in sys.argv:
        result_rows = run_pressure_ladder(flow_case, component_list)
        print(format_ladder(result_rows))
    else:
        launch_gui()
