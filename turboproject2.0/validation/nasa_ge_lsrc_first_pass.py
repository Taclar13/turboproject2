"""Reproduce the first NASA/GE LSRC comparison using the unchanged P2 solver.

This is a transparent work/thermodynamic consistency comparison, not an
independent efficiency prediction. See README.md in this directory for data
provenance, conversions, and limitations.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compressor import Assumptions, CompressorDesign, CP, GAMMA, solve_compressor


def main() -> None:
    # NASA-CR-159498 design point. The report plots the four-stage coefficients
    # divided by four; use the published per-stage coefficients below.
    stage_count = 4
    tip_speed = 63.82  # m/s
    pressure_coefficient_per_stage = 0.555
    nasa_flow_coefficient = 0.407
    measured_overall_efficiency_used_as_model_stage_eta = 0.900

    # Reference inlet state used for this first conversion. The report's
    # performance coefficients normalize to reference density/velocity; this
    # standard-air state is an explicit approximation, not a claimed measured
    # inlet condition.
    inlet_t0 = 288.15  # K
    inlet_p0 = 101325.0  # Pa

    # Constant-annulus geometry from the NASA/GE LSRC report.
    tip_radius = 1.524 / 2.0
    hub_tip_ratio = 0.85
    hub_radius = hub_tip_ratio * tip_radius
    mean_radius = (tip_radius + hub_radius) / 2.0
    annulus_area = math.pi * (tip_radius**2 - hub_radius**2)

    omega = tip_speed / tip_radius
    rpm = omega * 60.0 / (2.0 * math.pi)
    mean_blade_speed = omega * mean_radius
    phi_p2 = nasa_flow_coefficient * tip_speed / mean_blade_speed

    # The report's pressure coefficient is isentropic enthalpy rise normalized
    # by 1/2 U_tip^2. Multiply its per-stage value by four for total rise.
    isentropic_work_per_stage = pressure_coefficient_per_stage * 0.5 * tip_speed**2
    isentropic_work_total = stage_count * isentropic_work_per_stage
    pressure_ratio_target = (
        1.0 + isentropic_work_total / (CP * inlet_t0)
    ) ** (GAMMA / (GAMMA - 1.0))

    # Convert NASA's tip-speed flow coefficient into the model's mean-radius
    # definition. Use reference density and nominal annulus area to estimate
    # mass flow for P2's continuity-based geometry calculation.
    gas_r = CP * (GAMMA - 1.0) / GAMMA
    reference_density = inlet_p0 / (gas_r * inlet_t0)
    axial_velocity_from_nasa_phi = nasa_flow_coefficient * tip_speed
    estimated_mass_flow = (
        reference_density * annulus_area * axial_velocity_from_nasa_phi
    )

    design = CompressorDesign(
        inlet_t0=inlet_t0,
        inlet_p0=inlet_p0,
        mass_flow=estimated_mass_flow,
        pressure_ratio_target=pressure_ratio_target,
        stage_count=stage_count,
        rpm=rpm,
        phi=phi_p2,
        radius_in=mean_radius,
        radius_out=mean_radius,
        radius_profile="constant",
        loading_profile="constant",
    )

    # Only broaden physical-screen ranges needed to calculate this research
    # rig's scale and low speed. No screening report is interpreted as a
    # feasibility pass/fail for this benchmark.
    assumptions = Assumptions(
        stage_eta=measured_overall_efficiency_used_as_model_stage_eta,
        max_mean_radius=1.0,
        min_mean_radius=0.01,
        max_tip_radius=2.0,
        min_rpm=0.01,
        max_rpm=25000.0,
        min_stage_loading=0.0,
    )
    result = solve_compressor(design, assumptions)

    # Measured third-stage, 50%-immersion vector data from NASA Table 13.
    # Tangential velocities and axial velocities are read from the published
    # table; U at this radius is inferred from its reported tip speed and the
    # rig's hub/tip geometry.
    nasa_stage3_midspan = {
        "immersion_percent": 50.0,
        "vtheta_in_m_s": 10.8,
        "vtheta_out_m_s": 32.6,
        "vx_in_m_s": 29.4,
        "vx_out_m_s": 28.3,
        "u_midspan_m_s": mean_blade_speed,
    }
    nasa_stage3_euler_work = nasa_stage3_midspan["u_midspan_m_s"] * (
        nasa_stage3_midspan["vtheta_out_m_s"]
        - nasa_stage3_midspan["vtheta_in_m_s"]
    )
    p2_stage3 = result["stages"][2]
    p2_stage3_euler_work = p2_stage3["delta_h0"]
    nasa_reference_specific_work = (
        isentropic_work_total / measured_overall_efficiency_used_as_model_stage_eta
    )

    report = {
        "status": "first-pass thermodynamic/work consistency comparison",
        "model_source_modified": False,
        "nasa_reference": {
            "report": "NASA-CR-159498, Core Compressor Exit Stage Study, Volume II",
            "design_point": "four-stage Rotor A/Stator A, design pressure-rise point",
            "stage_count": stage_count,
            "tip_speed_m_s": tip_speed,
            "flow_coefficient_tip_based": nasa_flow_coefficient,
            "pressure_coefficient_per_stage": pressure_coefficient_per_stage,
            "overall_efficiency_used_as_each_p2_stage_efficiency": measured_overall_efficiency_used_as_model_stage_eta,
            "reference_inlet_t0_K_assumed": inlet_t0,
            "reference_inlet_p0_Pa_assumed": inlet_p0,
            "tip_radius_m": tip_radius,
            "hub_radius_m": hub_radius,
            "mean_radius_m": mean_radius,
            "rpm_derived_from_tip_speed": rpm,
            "mass_flow_kg_s_estimated_from_reference_density_phi_and_annulus": estimated_mass_flow,
            "specific_work_j_kg_from_pressure_coefficient_and_efficiency": nasa_reference_specific_work,
            "stage3_midspan_vector_data": nasa_stage3_midspan,
            "stage3_midspan_euler_work_j_kg": nasa_stage3_euler_work,
        },
        "p2": {
            "design": result["design"],
            "assumptions": result["assumptions"],
            "achieved_pressure_ratio": result["achieved_pressure_ratio"],
            "specific_work_j_kg": result["specific_work_j_kg"],
            "shaft_power_w_using_estimated_mass_flow": result["shaft_power_w"],
            "outlet_t0_K": result["outlet_t0"],
            "stage_results": [
                {
                    "stage": s["stage"],
                    "stage_pressure_ratio": s["stage_pressure_ratio"],
                    "specific_work_j_kg": s["delta_h0"],
                    "loading_coefficient": s["psi"],
                    "u_mean_m_s": s["u_mean"],
                    "vx_m_s": s["vx"],
                    "vtheta_in_m_s": s["rotor_inlet"]["v_theta"],
                    "vtheta_out_m_s": s["rotor_exit"]["v_theta"],
                    "rotor_exit_alpha_deg": s["rotor_exit"]["alpha_deg"],
                    "rotor_exit_height_m": s["rotor_exit"]["blade_height"],
                    "rotor_exit_hub_radius_m": s["rotor_exit"]["r_hub"],
                    "rotor_exit_tip_radius_m": s["rotor_exit"]["r_tip"],
                }
                for s in result["stages"]
            ],
            "stage3_euler_work_j_kg": p2_stage3_euler_work,
        },
        "comparisons": {
            "overall_work_difference_percent_vs_coefficient_eta_reference": 100.0
            * (result["specific_work_j_kg"] / nasa_reference_specific_work - 1.0),
            "stage3_midspan_euler_work_difference_percent": 100.0
            * (p2_stage3_euler_work / nasa_stage3_euler_work - 1.0),
            "stage3_p2_vtheta_in_m_s": p2_stage3["rotor_inlet"]["v_theta"],
            "stage3_p2_vtheta_out_m_s": p2_stage3["rotor_exit"]["v_theta"],
            "stage3_p2_vtheta_change_m_s": p2_stage3["rotor_exit"]["v_theta"]
            - p2_stage3["rotor_inlet"]["v_theta"],
            "stage3_nasa_vtheta_change_m_s": nasa_stage3_midspan["vtheta_out_m_s"]
            - nasa_stage3_midspan["vtheta_in_m_s"],
            "stage1_p2_blade_height_m": result["stages"][0]["rotor_inlet"]["blade_height"],
            "stage4_p2_blade_height_m": result["stages"][-1]["rotor_inlet"]["blade_height"],
            "nasa_constant_blade_height_m": tip_radius - hub_radius,
        },
        "interpretation": [
            "The overall work agreement is not an independent prediction: measured efficiency and pressure coefficient are used to construct the comparison target.",
            "The midspan third-stage Euler comparison is more informative. It compares measured swirl change to the model's Euler work while exposing the zero-inlet-swirl difference.",
            "The P2 mean-line flow coefficient is converted from NASA's tip-speed definition. Mass flow is estimated from standard-air reference density and nominal annulus area because an exact matched mass-flow datum was not transcribed for this run.",
            "P2 retains constant cp/gamma, constant mean radius, equal stage loading, zero rotor-inlet swirl, and stator swirl removal. NASA's measured stage flow does not share all of these assumptions.",
        ],
    }

    output = Path(__file__).parent / "nasa_ge_lsrc_first_pass_results.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
