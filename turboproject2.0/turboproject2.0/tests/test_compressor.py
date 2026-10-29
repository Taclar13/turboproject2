import math
import unittest
from dataclasses import replace

from compressor import (
    Assumptions,
    CompressorDesign,
    CP,
    GAMMA,
    GAS_R,
    ModelError,
    design_profiles,
    evaluate_design,
    solve_compressor,
)


def example_design() -> CompressorDesign:
    return CompressorDesign(
        inlet_t0=300.0,
        inlet_p0=101325.0,
        mass_flow=16.0,
        pressure_ratio_target=5.0,
        stage_count=10,
        rpm=9000.0,
        phi=0.55,
        radius_in=0.23,
        radius_out=0.21,
        radius_profile="linear",
        loading_profile="constant",
        loading_shape=0.0,
    )


class CompressorVerificationTests(unittest.TestCase):
    def test_one_stage_against_independent_constant_property_equation(self):
        pr = 1.265
        eta = Assumptions().stage_eta
        design = CompressorDesign(
            inlet_t0=300.0,
            inlet_p0=101325.0,
            mass_flow=2.0,
            pressure_ratio_target=pr,
            stage_count=1,
            rpm=9000.0,
            phi=0.55,
            radius_in=0.20,
            radius_out=0.20,
            radius_profile="constant",
        )
        result = solve_compressor(design)
        stage = result["stages"][0]
        expected_dh = CP * design.inlet_t0 * (
            pr ** ((GAMMA - 1.0) / GAMMA) - 1.0
        ) / eta
        self.assertAlmostEqual(stage["delta_h0"], expected_dh, delta=1e-8)
        self.assertAlmostEqual(stage["stage_pressure_ratio"], pr, places=10)

    def test_overall_pressure_ratio_closes_at_requested_value(self):
        result = solve_compressor(example_design())
        product = math.prod(s["stage_pressure_ratio"] for s in result["stages"])
        self.assertAlmostEqual(product, 5.0, places=10)
        self.assertAlmostEqual(result["achieved_pressure_ratio"], 5.0, places=10)
        self.assertAlmostEqual(result["outlet_p0"] / 101325.0, 5.0, places=10)

    def test_euler_work_and_total_temperature_rise_per_stage(self):
        design = example_design()
        result = solve_compressor(design)
        omega = 2.0 * math.pi * design.rpm / 60.0
        for stage in result["stages"]:
            u = omega * stage["radius_mean"]
            euler_work = u * (
                stage["rotor_exit"]["v_theta"]
                - stage["rotor_inlet"]["v_theta"]
            )
            thermal_work = CP * (
                stage["t0_out"] - stage["t0_in"]
            )
            self.assertAlmostEqual(euler_work, stage["delta_h0"], delta=1e-8)
            self.assertAlmostEqual(thermal_work, stage["delta_h0"], delta=1e-8)

    def test_station_continuity_and_exact_arithmetic_mean_annulus_area(self):
        design = example_design()
        result = solve_compressor(design)
        for station in result["stations"]:
            flow = station["rho"] * station["area"] * station["vx"]
            geometric_area = math.pi * (
                station["r_tip"] ** 2 - station["r_hub"] ** 2
            )
            self.assertAlmostEqual(flow, design.mass_flow, delta=1e-9)
            self.assertAlmostEqual(geometric_area, station["area"], delta=1e-12)

    def test_rotor_is_lossless_and_stator_accounts_for_stage_entropy_generation(self):
        result = solve_compressor(example_design())
        for stage in result["stages"]:
            s_total = stage["entropy_generation_j_kgk"]
            s_rotor = (
                CP * math.log(stage["t0_rotor_exit"] / stage["t0_in"])
                - GAS_R * math.log(stage["p0_rotor_exit"] / stage["p0_in"])
            )
            s_stator = -GAS_R * math.log(
                stage["p0_out"] / stage["p0_rotor_exit"]
            )
            self.assertGreaterEqual(s_total, 0.0)
            self.assertAlmostEqual(s_rotor, 0.0, delta=1e-9)
            self.assertAlmostEqual(s_stator, s_total, delta=1e-9)

    def test_named_example_passes_selected_screens_and_default_reports_failures(self):
        example = solve_compressor(example_design())
        example_report = evaluate_design(example)
        self.assertTrue(example_report["passes_selected_screening"])
        self.assertEqual(example_report["number_failed"], 0)
        self.assertEqual(example_report["number_of_checks"], 326)

        default = solve_compressor(CompressorDesign())
        default_report = evaluate_design(default)
        self.assertAlmostEqual(default["achieved_pressure_ratio"], 14.0, places=9)
        self.assertFalse(default_report["passes_selected_screening"])
        self.assertGreater(default_report["number_failed"], 0)

    def test_loading_profiles_are_normalized_to_unit_mean(self):
        for profile, shape in (("constant", 0.0), ("linear", 0.25), ("exponential", -0.6)):
            design = replace(example_design(), loading_profile=profile, loading_shape=shape)
            _, _, weights = design_profiles(design)
            self.assertTrue(all(w > 0.0 for w in weights))
            self.assertAlmostEqual(float(weights.mean()), 1.0, places=12)

    def test_invalid_pressure_ratio_is_rejected(self):
        with self.assertRaises(ModelError):
            solve_compressor(replace(example_design(), pressure_ratio_target=1.0))


if __name__ == "__main__":
    unittest.main()
