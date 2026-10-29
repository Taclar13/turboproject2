"""Project 2: manually controlled, design-point, 1-D mean-line axial-compressor study.

One call solves ONE user-specified design. No mass minimization, hidden redesign,
automatic mass-flow sweep, off-design map, actual blade-design or validated efficiency.
All inputs and outputs use SI units except rpm, degrees and dimensionless ratios.

Model pedigree: refactored from the user's earlier Project 2 exploratory solver.
Its constant-property air, fixed stage efficiency, aerodynamic screening thresholds,
geometry approximations and mass estimates are illustrative, NOT engineering validation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq

CP = 1005.0  # J/(kg K), constant-air approximation
GAMMA = 1.4
GAS_R = CP * (GAMMA - 1.0) / GAMMA


class ModelError(ValueError):
    """Invalid input, unrepresentable geometric state, or failed numerical solve."""


@dataclass(frozen=True)
class CompressorDesign:
    """Everything the user chooses for ONE compressor; no parameter optimization."""

    inlet_t0: float = 300.0              # K
    inlet_p0: float = 101325.0           # Pa
    mass_flow: float = 16.0              # kg/s
    pressure_ratio_target: float = 14.0  # stagnation pressure ratio
    stage_count: int = 16
    rpm: float = 10500.0
    phi: float = 0.55                    # same Vx/Umean at all modeled stages
    radius_in: float = 0.25             # first-stage arithmetic annulus mean radius, m
    radius_out: float = 0.21            # last-stage arithmetic annulus mean radius, m
    radius_profile: str = "linear"       # constant, linear, power
    radius_exponent: float = 1.0        # power-profile p > 0; ignored for others
    loading_profile: str = "constant"    # constant, linear, exponential
    loading_shape: float = 0.0          # linear a: w=1+a(1-2x); exponential k: w=exp(k(x-.5))


@dataclass(frozen=True)
class Assumptions:
    """Fixed *illustrative* physics, geometry, mass and screening inputs."""

    stage_eta: float = 0.85             # ASSUMED stage stagnation isentropic efficiency
    blade_density: float = 2700.0
    casing_density: float = 2700.0
    drum_density: float = 2700.0
    shaft_density: float = 2700.0
    casing_thickness: float = 0.004
    drum_thickness: float = 0.006
    shaft_radius: float = 0.045
    rotor_disk_thickness: float = 0.006
    stator_mount_radial_thickness: float = 0.008
    stator_mount_axial_width: float = 0.010
    blade_shape_factor: float = 0.65
    rotor_aspect_ratio: float = 1.6
    stator_aspect_ratio: float = 1.5
    rotor_solidity_target: float = 1.25
    stator_solidity_target: float = 1.30
    rotor_thickness_ratio: float = 0.10
    stator_thickness_ratio: float = 0.09
    rotor_stagger_deg: float = 50.0
    stator_stagger_deg: float = 40.0
    rotor_gap_ratio: float = 0.30
    stator_gap_ratio: float = 0.35
    min_hub_radius: float = 0.075
    min_blade_height: float = 0.018
    max_tip_radius: float = 0.44
    max_wall_slope: float = 0.65  # illustrative, NOT a universal blade-design limit
    max_rotor_relative_mach: float = 0.98
    max_rotor_tip_mach_estimate: float = 1.03
    max_stator_absolute_mach: float = 0.85
    max_diffusion: float = 0.53
    min_reaction: float = 0.45
    max_reaction: float = 0.92
    min_stage_loading: float = 0.18
    max_stage_loading: float = 0.60
    min_rpm: float = 6500.0
    max_rpm: float = 18000.0
    min_mean_radius: float = 0.14
    max_mean_radius: float = 0.32
    min_phi: float = 0.38
    max_phi: float = 0.78


def _positive(name: str, value: float) -> None:
    if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ModelError(f"{name} must be finite and positive; got {value!r}")


def validate_inputs(d: CompressorDesign, a: Assumptions) -> None:
    for key in ("inlet_t0", "inlet_p0", "mass_flow", "pressure_ratio_target", "rpm",
                "phi", "radius_in", "radius_out", "radius_exponent"):
        _positive(key, getattr(d, key))
    if d.pressure_ratio_target <= 1:
        raise ModelError("pressure_ratio_target must exceed 1 for this compressor model")
    if type(d.stage_count) is not int or d.stage_count < 1:
        raise ModelError("stage_count must be a positive integer")
    if d.stage_count == 1 and not math.isclose(d.radius_in, d.radius_out, rel_tol=0, abs_tol=1e-12):
        raise ModelError("For one stage, radius_in and radius_out must be equal")
    if d.radius_profile not in ("constant", "linear", "power"):
        raise ModelError("radius_profile must be: constant, linear or power")
    if d.radius_profile == "constant" and not math.isclose(d.radius_in, d.radius_out, abs_tol=1e-12):
        raise ModelError("Constant radius requires radius_in == radius_out; inputs are never silently changed")
    if d.loading_profile not in ("constant", "linear", "exponential"):
        raise ModelError("loading_profile must be: constant, linear or exponential")
    if not math.isfinite(d.loading_shape):
        raise ModelError("loading_shape must be finite")
    if d.loading_profile == "linear" and abs(d.loading_shape) >= 1:
        raise ModelError("For linear loading, abs(loading_shape) must be less than 1")
    if d.loading_profile == "constant" and d.loading_shape != 0:
        raise ModelError("For constant loading, set loading_shape=0 (no ignored nonzero input)")
    if d.radius_profile != "power" and d.radius_exponent != 1.0:
        raise ModelError("radius_exponent is used only by radius_profile='power'; set it to 1 otherwise")
    if not 0 < a.stage_eta <= 1:
        raise ModelError("stage_eta must lie in (0,1]")
    for key in ("blade_density", "casing_density", "drum_density", "shaft_density",
                "casing_thickness", "drum_thickness", "shaft_radius", "rotor_disk_thickness",
                "stator_mount_radial_thickness", "stator_mount_axial_width", "blade_shape_factor",
                "rotor_aspect_ratio", "stator_aspect_ratio", "rotor_solidity_target",
                "stator_solidity_target", "rotor_thickness_ratio", "stator_thickness_ratio",
                "min_hub_radius", "min_blade_height", "max_tip_radius", "max_wall_slope"):
        _positive(key, getattr(a, key))
    for key in ("rotor_stagger_deg", "stator_stagger_deg"):
        if not math.isfinite(getattr(a, key)) or not 0 <= getattr(a, key) < 90:
            raise ModelError(f"{key} must be in [0,90) degrees")
    for key in ("rotor_gap_ratio", "stator_gap_ratio"):
        if not math.isfinite(getattr(a, key)) or getattr(a, key) <= 0:
            raise ModelError(f"{key} must be positive to separate axial stations")


def design_profiles(d: CompressorDesign) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return normalized stage position, stage mean radii, mean-one positive psi weights."""
    x = np.linspace(0.0, 1.0, d.stage_count) if d.stage_count > 1 else np.array([0.5])
    if d.radius_profile == "constant":
        radius = np.full(d.stage_count, d.radius_in)
    else:
        frac = x if d.radius_profile == "linear" else x ** d.radius_exponent
        radius = d.radius_in + (d.radius_out - d.radius_in) * frac
    if d.loading_profile == "constant":
        raw = np.ones(d.stage_count)
    elif d.loading_profile == "linear":
        raw = 1.0 + d.loading_shape * (1.0 - 2.0 * x)
    else:
        # Shift before exponentiation to avoid overflow for legitimate moderate k;
        # extreme k values are rejected as numerically unrepresentable.
        z = d.loading_shape * (x - 0.5)
        if np.max(z) - np.min(z) > 600:
            raise ModelError("exponential loading_shape is too extreme to represent reliably")
        raw = np.exp(z - np.max(z))
    if not np.all(np.isfinite(raw)) or np.any(raw <= 0):
        raise ModelError("loading distribution produced nonpositive/nonfinite stage weights")
    return x, radius, raw / float(np.mean(raw))


def find_loading_scale(d: CompressorDesign, a: Assumptions,
                       radii: np.ndarray, loading_weights: np.ndarray) -> float:
    """Scalar root: find psi_i = scale * weight_i so final modeled PR hits target."""
    omega = 2 * math.pi * d.rpm / 60.0
    work_per_unit_scale = loading_weights * (omega * radii) ** 2
    target_log_pr = math.log(d.pressure_ratio_target)

    def residual(scale: float) -> float:
        t0 = d.inlet_t0
        accumulated_log_pr = 0.0
        for dh in work_per_unit_scale * scale:
            dt = float(dh) / CP
            relative_isentropic_rise = a.stage_eta * dt / t0
            if not math.isfinite(relative_isentropic_rise):
                raise ModelError("Pressure-ratio root trial overflowed")
            accumulated_log_pr += GAMMA / (GAMMA - 1) * math.log1p(relative_isentropic_rise)
            t0 += dt
            if not math.isfinite(t0):
                raise ModelError("Pressure-ratio root trial produced invalid temperature")
        return accumulated_log_pr - target_log_pr

    upper = 0.1
    while residual(upper) < 0:
        upper *= 2.0
        if upper > 1.0e5:
            raise ModelError("Could not bracket the requested pressure ratio")
    return float(brentq(residual, 0, upper, xtol=1e-12, rtol=1e-12))


def static_state(t0: float, p0: float, v_abs: float) -> dict[str, float]:
    if not all(math.isfinite(v) for v in (t0, p0, v_abs)) or min(t0, p0) <= 0 or v_abs < 0:
        raise ModelError("Invalid stagnation properties or absolute velocity")
    t = t0 - v_abs * v_abs / (2 * CP)
    if t <= 0:
        raise ModelError("Flow speed implies nonpositive static temperature")
    p = p0 * (t / t0) ** (GAMMA / (GAMMA - 1))
    if p <= 0:
        raise ModelError("Static pressure is numerically unrepresentable")
    return {"t": t, "p": p, "rho": p / (GAS_R * t), "a": math.sqrt(GAMMA * GAS_R * t)}


def station(name: str, stage: int, t0: float, p0: float, d: CompressorDesign,
            rm: float, omega: float, vtheta: float) -> dict[str, Any]:
    u = omega * rm
    vx = d.phi * u
    v = math.hypot(vx, vtheta)
    state = static_state(t0, p0, v)
    area = d.mass_flow / (state["rho"] * vx)
    # Exact annulus area for ARITHMETIC mean radius: A = 2*pi*rm*(rt-rh).
    height = area / (2 * math.pi * rm)
    rh, rt = rm - height / 2, rm + height / 2
    if rh <= 0:
        raise ModelError(f"Stage {stage} {name}: annulus hub crosses centerline; cannot construct geometry")
    wtheta = vtheta - u
    w = math.hypot(vx, wtheta)
    # Only an indicative tip estimate: Vx and Vtheta assumed uniform over blade span.
    w_tip = math.hypot(vx, vtheta - omega * rt)
    return {"stage": stage, "name": name, "x": 0.0,
            "t0": t0, "p0": p0, "t": state["t"], "p": state["p"],
            "rho": state["rho"], "a": state["a"], "r_mean": rm, "r_hub": rh,
            "r_tip": rt, "blade_height": height, "area": area, "u_mean": u,
            "vx": vx, "v_theta": vtheta, "v_abs": v, "w_rel": w,
            "w_theta": wtheta, "alpha_deg": math.degrees(math.atan2(vtheta, vx)),
            "beta_deg": math.degrees(math.atan2(wtheta, vx)),
            "mach_abs": v / state["a"], "mach_rel": w / state["a"],
            "mach_tip_est": w_tip / state["a"]}


def row_geometry(rm: float, span: float, aspect: float, target_solidity: float,
                 thickness_ratio: float, stagger: float, gap_ratio: float,
                 blade_density: float, shape_factor: float) -> dict[str, Any]:
    chord = span / aspect
    requested_pitch = chord / target_solidity
    count = max(1, round(2 * math.pi * rm / requested_pitch))
    pitch = 2 * math.pi * rm / count
    axial_chord = chord * math.cos(math.radians(stagger))
    if axial_chord <= 0:
        raise ModelError("Nonpositive axial blade chord")
    # Explicit volume proxy, not a physical three-dimensional blade profile.
    mass = count * blade_density * shape_factor * chord * (thickness_ratio * chord) * span
    return {"count": count, "span": span, "chord": chord, "pitch": pitch,
            "solidity": chord / pitch, "thickness_max": thickness_ratio * chord,
            "axial_chord": axial_chord, "gap": gap_ratio * axial_chord, "mass": mass}


def diffusion(speed_in: float, speed_out: float,
              theta_in: float, theta_out: float, solidity: float) -> float:
    if min(speed_in, speed_out, solidity) <= 0:
        raise ModelError("Invalid diffusion-factor inputs")
    # Simple Lieblein-style cascade screening indicator, no 3-D loss prediction.
    return 1 - speed_out / speed_in + abs(theta_out - theta_in) / (2 * solidity * speed_in)


def shell_mass(profile: list[dict[str, Any]], radius_key: str,
               thickness: float, density: float) -> float:
    x = np.array([p["x"] for p in profile])
    radii = np.array([p[radius_key] for p in profile])
    if len(x) < 2 or np.any(np.diff(x) <= 0) or np.any(radii <= 0):
        raise ModelError("Nonmonotone axial stations or invalid shell radius")
    segment_slant = np.hypot(np.diff(x), np.diff(radii))
    return float(np.sum(math.pi * (radii[:-1] + radii[1:]) * segment_slant * thickness * density))


def solve_compressor(d: CompressorDesign, a: Assumptions = Assumptions()) -> dict[str, Any]:
    """Calculate ONE prescribed design. Do not optimize or adjust the inputs for feasibility."""
    validate_inputs(d, a)
    xnorm, radii, weights = design_profiles(d)
    scale = find_loading_scale(d, a, radii, weights)
    omega = 2 * math.pi * d.rpm / 60
    t0, p0 = d.inlet_t0, d.inlet_p0
    x_axial = 0.0
    stages: list[dict[str, Any]] = []
    profile: list[dict[str, Any]] = []
    mass_blades_rotor = mass_blades_stator = disk_mass = stator_mount_mass = 0.0

    for i, (xn, rm, weight) in enumerate(zip(xnorm, radii, weights), start=1):
        u = omega * float(rm)
        psi = scale * float(weight)
        dh = psi * u*u
        t0_out = t0 + dh / CP
        stage_pr = math.exp(GAMMA / (GAMMA - 1) * math.log1p(a.stage_eta * dh / (CP * t0)))
        p0_out = p0 * stage_pr
        # Lossless rotor: all modeled stage entropy generation is assigned to the stator.
        s_generation = CP * math.log(t0_out/t0) - GAS_R * math.log(stage_pr)
        if s_generation < -1e-8:
            raise ModelError("Calculated negative stage entropy generation")
        p0_rotor = p0 * (t0_out / t0) ** (CP / GAS_R)
        vtheta2 = dh / u  # Euler: dh=U(Vtheta2-Vtheta1), Vtheta1=0

        st1 = station("rotor_inlet", i, t0, p0, d, float(rm), omega, 0.0)
        st2 = station("rotor_exit", i, t0_out, p0_rotor, d, float(rm), omega, vtheta2)
        st3 = station("stator_exit", i, t0_out, p0_out, d, float(rm), omega, 0.0)
        rotor_span = (st1["blade_height"] + st2["blade_height"]) / 2
        stator_span = (st2["blade_height"] + st3["blade_height"]) / 2
        rotor = row_geometry(float(rm), rotor_span, a.rotor_aspect_ratio,
                             a.rotor_solidity_target, a.rotor_thickness_ratio,
                             a.rotor_stagger_deg, a.rotor_gap_ratio,
                             a.blade_density, a.blade_shape_factor)
        stator = row_geometry(float(rm), stator_span, a.stator_aspect_ratio,
                              a.stator_solidity_target, a.stator_thickness_ratio,
                              a.stator_stagger_deg, a.stator_gap_ratio,
                              a.blade_density, a.blade_shape_factor)
        mass_blades_rotor += rotor["mass"]
        mass_blades_stator += stator["mass"]
        # Simplified, NON-structurally-sized rotor disk and stator mount.
        disk_radius = (st1["r_hub"] + st2["r_hub"]) / 2 - a.drum_thickness
        if disk_radius <= a.shaft_radius:
            raise ModelError(f"Stage {i}: disk radius <= assumed shaft radius")
        disk_mass += math.pi * (disk_radius**2 - a.shaft_radius**2) * (
            a.rotor_disk_thickness * a.drum_density)
        mount_radius = (st2["r_tip"] + st3["r_tip"]) / 2 + a.casing_thickness
        stator_mount_mass += math.pi * (
            (mount_radius + a.stator_mount_radial_thickness)**2 - mount_radius**2
        ) * a.stator_mount_axial_width * a.casing_density

        st1["x"] = x_axial
        x_axial += rotor["axial_chord"]
        st2["x"] = x_axial
        x_axial += rotor["gap"]
        st2_gap = dict(st2, name="stator_inlet", x=x_axial)
        x_axial += stator["axial_chord"]
        st3["x"] = x_axial
        profile.extend((st1, st2, st2_gap, st3))
        if i < d.stage_count:
            x_axial += stator["gap"]

        rotor_d = diffusion(st1["w_rel"], st2["w_rel"],
                            st1["w_theta"], st2["w_theta"], rotor["solidity"])
        stator_d = diffusion(st2["v_abs"], st3["v_abs"],
                             st2["v_theta"], st3["v_theta"], stator["solidity"])
        actual_static_stage_rise = CP * (st3["t"] - st1["t"])
        if actual_static_stage_rise <= 0:
            raise ModelError(f"Stage {i}: nonpositive static enthalpy rise")
        reaction = CP * (st2["t"] - st1["t"]) / actual_static_stage_rise
        stages.append({"stage": i, "normalized_position": float(xn),
                       "radius_mean": float(rm), "u_mean": u, "vx": d.phi*u,
                       "loading_weight": float(weight), "psi": psi,
                       "delta_h0": dh, "stage_pressure_ratio": stage_pr,
                       "t0_in": t0, "p0_in": p0,
                       "t0_rotor_exit": t0_out, "p0_rotor_exit": p0_rotor,
                       "t0_out": t0_out, "p0_out": p0_out,
                       "entropy_generation_j_kgk": s_generation,
                       "reaction": reaction, "rotor_diffusion": rotor_d,
                       "stator_diffusion": stator_d,
                       "rotor": rotor, "stator": stator,
                       "rotor_inlet": st1, "rotor_exit": st2, "stator_exit": st3})
        t0, p0 = t0_out, p0_out

    mass_parts = {"rotor_blades": mass_blades_rotor,
                  "stator_blades": mass_blades_stator,
                  "casing": shell_mass(profile, "r_tip", a.casing_thickness, a.casing_density),
                  "rotor_drum": shell_mass(profile, "r_hub", a.drum_thickness, a.drum_density),
                  "shaft_core": math.pi * a.shaft_radius**2 * x_axial * a.shaft_density,
                  "rotor_disks": disk_mass, "stator_mounts": stator_mount_mass}
    achieved_pr = p0 / d.inlet_p0
    if not math.isclose(achieved_pr, d.pressure_ratio_target, rel_tol=1e-9, abs_tol=1e-10):
        raise ModelError("Pressure-ratio root solution did not satisfy target")
    return {"design": asdict(d), "assumptions": asdict(a),
            "loading_scale": scale, "achieved_pressure_ratio": achieved_pr,
            "outlet_t0": t0, "outlet_p0": p0,
            "specific_work_j_kg": CP * (t0-d.inlet_t0),
            "shaft_power_w": d.mass_flow * CP * (t0-d.inlet_t0),
            "length_m": x_axial, "diameter_m": 2*max(p["r_tip"] for p in profile),
            "mass_parts_kg": mass_parts, "estimated_mass_kg": sum(mass_parts.values()),
            "stages": stages, "stations": profile,
            "disclaimer": "Unvalidated mean-line screening and rough hardware mass proxy; "
                          "stage efficiency, 3D aerodynamics, structure and manufacturing are NOT predicted."}


def evaluate_design(result: dict[str, Any], a: Assumptions | None = None) -> dict[str, Any]:
    """Report EACH modeled screening condition; no hidden penalty/aggregate feasibility score.

    Default to the exact assumptions used in the supplied solve; never silently
    evaluate one physics model with a different set of screening inputs.
    """
    original_assumptions = Assumptions(**result["assumptions"])
    if a is None:
        a = original_assumptions
    elif a != original_assumptions:
        raise ModelError("Evaluation assumptions differ from those used by the solver")
    checks: list[dict[str, Any]] = []

    def check(name: str, actual: float, *, stage: int | None = None, location: str = "",
              minimum: float | None = None, maximum: float | None = None,
              units: str = "") -> None:
        for bound_name, bound, raw_margin in (("minimum", minimum, actual-(minimum or 0)),
                                               ("maximum", maximum, (maximum or 0)-actual)):
            if bound is None:
                continue
            checks.append({"name": name, "stage": stage, "location": location,
                           "actual": float(actual), "bound_type": bound_name,
                           "limit": float(bound), "margin": float(raw_margin),
                           "units": units, "pass": raw_margin >= -1e-9})

    for stage in result["stages"]:
        i = stage["stage"]
        st1, st2, st3 = stage["rotor_inlet"], stage["rotor_exit"], stage["stator_exit"]
        check("stage_loading", stage["psi"], stage=i, minimum=a.min_stage_loading,
              maximum=a.max_stage_loading)
        check("degree_of_reaction", stage["reaction"], stage=i,
              minimum=a.min_reaction, maximum=a.max_reaction)
        for label, st in (("rotor inlet", st1), ("rotor exit", st2)):
            check("rotor_relative_mach", st["mach_rel"], stage=i,
                  location=label, maximum=a.max_rotor_relative_mach)
            check("rotor_tip_mach_estimate", st["mach_tip_est"], stage=i,
                  location=label, maximum=a.max_rotor_tip_mach_estimate)
        for label, st in (("stator inlet", st2), ("stator exit", st3)):
            check("stator_absolute_mach", st["mach_abs"], stage=i,
                  location=label, maximum=a.max_stator_absolute_mach)
        check("rotor_diffusion_factor", stage["rotor_diffusion"], stage=i,
              maximum=a.max_diffusion)
        check("stator_diffusion_factor", stage["stator_diffusion"], stage=i,
              maximum=a.max_diffusion)

    stations = result["stations"]
    for st in stations:
        check("hub_radius", st["r_hub"], stage=st["stage"], location=st["name"],
              minimum=max(a.min_hub_radius, a.shaft_radius + a.drum_thickness), units="m")
        check("blade_height", st["blade_height"], stage=st["stage"],
              location=st["name"], minimum=a.min_blade_height, units="m")
        check("tip_radius", st["r_tip"], stage=st["stage"],
              location=st["name"], maximum=a.max_tip_radius, units="m")
    for left, right in zip(stations, stations[1:]):
        dz = right["x"] - left["x"]
        if dz <= 0:
            raise ModelError("Nonmonotone station x coordinates during feasibility check")
        for key, name in (("r_hub", "hub_wall_slope"), ("r_tip", "tip_wall_slope")):
            slope = abs((right[key]-left[key]) / dz)
            check(name, slope, stage=right["stage"],
                  location=f"{left['stage']}:{left['name']} -> {right['stage']}:{right['name']}",
                  maximum=a.max_wall_slope)
    d = result["design"]
    check("shaft_rpm_screen", d["rpm"], minimum=a.min_rpm, maximum=a.max_rpm, units="rpm")
    check("flow_coefficient_screen", d["phi"], minimum=a.min_phi, maximum=a.max_phi)
    for label in ("radius_in", "radius_out"):
        check(f"{label}_screen", d[label], minimum=a.min_mean_radius,
              maximum=a.max_mean_radius, units="m")
    failed = [c for c in checks if not c["pass"]]
    return {"passes_selected_screening": len(failed) == 0,
            "number_of_checks": len(checks), "number_failed": len(failed),
            "checks": checks, "failed": failed,
            "note": "Pass means only the chosen preliminary checks pass; not validation, "
                    "manufacturability, structural adequacy, or predicted real efficiency."}


def display_results(result: dict[str, Any], report: dict[str, Any]) -> None:
    d = result["design"]
    print("\nPROJECT 2 | ONE MANUALLY SPECIFIED COMPRESSOR")
    print(f"Inputs: {d['stage_count']} stages, {d['rpm']:.0f} rpm, {d['mass_flow']:.3f} kg/s, "
          f"PR target {d['pressure_ratio_target']:.5f}, phi {d['phi']:.3f}")
    print(f"Profiles: radius={d['radius_profile']} (power={d['radius_exponent']:g}); "
          f"loading={d['loading_profile']} (shape={d['loading_shape']:g})")
    print(f"Achieved PR {result['achieved_pressure_ratio']:.7f} | "
          f"T0 exit {result['outlet_t0']:.2f} K | "
          f"specific work {result['specific_work_j_kg']/1000:.2f} kJ/kg | "
          f"shaft power {result['shaft_power_w']/1e6:.3f} MW")
    print(f"Length {result['length_m']:.3f} m | maximum diameter {result['diameter_m']:.3f} m | "
          f"rough modeled mass {result['estimated_mass_kg']:.2f} kg")
    print("Preliminary screening:", "PASS" if report["passes_selected_screening"] else "FAIL",
          f"({report['number_failed']} failures among {report['number_of_checks']} individual checks)")
    for c in report["failed"][:15]:
        prefix = f"stage {c['stage']} " if c["stage"] is not None else ""
        print(f"  FAIL {prefix}{c['name']} {c['location']}: "
              f"actual={c['actual']:.5g}, {c['bound_type']}={c['limit']:.5g}, "
              f"margin={c['margin']:+.4g} {c['units']}")
    if report["number_failed"] > 15:
        print(f"  ... {report['number_failed']-15} more failures in constraints.csv")
    print("Stage | radius m | psi   | PR    | T0 out K | rotor D | stator D | reaction")
    for s in result["stages"]:
        print(f"{s['stage']:5d} | {s['radius_mean']:8.4f} | {s['psi']:5.3f} | "
              f"{s['stage_pressure_ratio']:5.3f} | {s['t0_out']:8.2f} | "
              f"{s['rotor_diffusion']:7.3f} | {s['stator_diffusion']:8.3f} | {s['reaction']:8.3f}")
    print("Caution: preliminary unvalidated physics, feasibility criteria, and mass estimates.\n")


def save_results(result: dict[str, Any], report: dict[str, Any], output: Path,
                 *, plots: bool = True) -> None:
    """Save exactly one design and its individual data; NEVER launch a sweep."""
    output.mkdir(parents=True, exist_ok=True)
    (output / "design_results.json").write_text(
        json.dumps({"result": result, "screening": report}, indent=2, allow_nan=False))

    def export_csv(name: str, rows: list[dict[str, Any]]) -> None:
        if not rows:
            return
        with (output / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    export_csv("stage_results.csv", [{k:v for k,v in s.items()
                if isinstance(v, (str, bool, int, float))} for s in result["stages"]])
    export_csv("station_results.csv", result["stations"])
    export_csv("constraints.csv", report["checks"])
    export_csv("mass_breakdown.csv", [{"component": k, "estimated_mass_kg": v}
                                       for k,v in result["mass_parts_kg"].items()])
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        stations = result["stations"]
        fig, ax = plt.subplots(figsize=(9, 4.2))
        ax.plot([p["x"] for p in stations], [p["r_hub"] for p in stations], ".-", label="Hub")
        ax.plot([p["x"] for p in stations], [p["r_tip"] for p in stations], ".-", label="Tip")
        ax.set(xlabel="Axial location (m)", ylabel="Radius (m)", title="Preliminary annulus flowpath")
        ax.grid(alpha=.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(output / "flowpath.png", dpi=170)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot([s["stage"] for s in result["stages"]],
                [s["psi"] for s in result["stages"]], "o-")
        ax.set(xlabel="Stage", ylabel="Stage loading coefficient (psi)",
               title="User-selected loading profile, scaled to target pressure ratio")
        ax.grid(alpha=.3)
        fig.tight_layout()
        fig.savefig(output / "loading_profile.png", dpi=170)
        plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--inlet-t0", type=float, default=300.0, help="Inlet stagnation T, K")
    p.add_argument("--inlet-p0", type=float, default=101325.0, help="Inlet stagnation p, Pa")
    p.add_argument("--mass-flow", type=float, default=16.0, help="Design mass flow, kg/s")
    p.add_argument("--pressure-ratio", type=float, default=14.0, help="Target overall stagnation PR")
    p.add_argument("--stage-count", type=int, default=16)
    p.add_argument("--rpm", type=float, default=10500.0)
    p.add_argument("--phi", type=float, default=0.55)
    p.add_argument("--radius-in", type=float, default=0.25, help="First-stage mean radius, m")
    p.add_argument("--radius-out", type=float, default=0.21, help="Last-stage mean radius, m")
    p.add_argument("--radius-profile", choices=("constant", "linear", "power"), default="linear")
    p.add_argument("--radius-exponent", type=float, default=1.0, help="Used for power profile")
    p.add_argument("--loading-profile", choices=("constant", "linear", "exponential"), default="constant")
    p.add_argument("--loading-shape", type=float, default=0.0,
                   help="linear a: positive means front-loaded; exponential k: negative front-loaded")
    p.add_argument("--output", type=Path, default=Path("results"))
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args()
    d = CompressorDesign(inlet_t0=args.inlet_t0, inlet_p0=args.inlet_p0,
                         mass_flow=args.mass_flow, pressure_ratio_target=args.pressure_ratio,
                         stage_count=args.stage_count, rpm=args.rpm, phi=args.phi,
                         radius_in=args.radius_in, radius_out=args.radius_out,
                         radius_profile=args.radius_profile, radius_exponent=args.radius_exponent,
                         loading_profile=args.loading_profile, loading_shape=args.loading_shape)
    try:
        result = solve_compressor(d)
        report = evaluate_design(result)
    except ModelError as exc:
        p.exit(2, f"Model could not represent this input design: {exc}\n")
    display_results(result, report)
    save_results(result, report, args.output, plots=not args.no_plots)
    print(f"Saved the single prescribed design and detailed reports to {args.output.resolve()}")


if __name__ == "__main__":
    main()
