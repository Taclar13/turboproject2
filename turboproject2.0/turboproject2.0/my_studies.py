"""Optional user-directed study scratchpad: manually change one input, rerun ONE design.

This script does not contain an optimizer, automatic sweep or mass-flow loop.
Edit `design = ...` and rerun to ask your next engineering question.
"""
from pathlib import Path
from compressor import CompressorDesign, evaluate_design, display_results, save_results, solve_compressor

# Start with a simple baseline that passes the current *illustrative* screens.
design = CompressorDesign(
    inlet_t0=300.0, inlet_p0=101325.0, mass_flow=16.0,
    pressure_ratio_target=5.0, stage_count=10, rpm=9000.0,
    phi=0.55, radius_in=0.23, radius_out=0.21,
    radius_profile="linear", radius_exponent=1.0,
    loading_profile="constant", loading_shape=0.0,
)

# For your NEXT run, manually edit just one value above, or replace the design:
# from dataclasses import replace
# design = replace(design, loading_profile="linear", loading_shape=0.25)
# design = replace(design, radius_profile="power", radius_exponent=1.5)
# design = replace(design, pressure_ratio_target=6.0)
# design = replace(design, rpm=9500.0)

result = solve_compressor(design)
report = evaluate_design(result)
display_results(result, report)
save_results(result, report, Path("my_study_results"))
