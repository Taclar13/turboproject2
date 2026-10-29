# NASA/GE LSRC first-pass comparison

## Scope

This first pass configures the unchanged P2 equations to the NASA/GE Low-Speed
Research Compressor (LSRC), four-stage Rotor A/Stator A design pressure-rise
point. P2 receives four stages, constant mean radius, constant loading, NASA's
design-point overall efficiency (0.90) applied as its common stage-efficiency
input, and the report's design-point pressure and flow coefficients. P2 retains
its zero-swirl rotor-inlet / stator-exit assumption and constant `cp`/`gamma`.

This is a thermodynamic and velocity-triangle consistency check. It is **not**
an independent prediction of compressor efficiency, and the close overall-work
match must not be described as independent validation: the pressure coefficient
and measured efficiency are used to construct the comparison work target.

## Source data and conversions

Primary source: NASA-CR-159498, *Core Compressor Exit Stage Study, Volume II:
Data and Performance Report for the Baseline Configuration*. The report PDF is
not redistributed with this project. Access it through the public NASA NTRS
record: <https://ntrs.nasa.gov/citations/19810007536>.

At the four-stage design pressure-rise point the report gives flow coefficient
0.407, pressure coefficient 0.555, and overall torque efficiency 0.900. The
report's performance plot presents the four-stage pressure/work coefficients
normalized per stage (`/4`). The run therefore uses 0.555 per stage and four
stages. NASA's tip speed is 63.82 m/s. Using the reported casing diameter of
1.524 m and hub-to-casing radius ratio 0.85 gives tip radius 0.762 m, hub radius
0.6477 m, mean radius 0.70485 m, and about 799.8 rpm.

The report's flow coefficient is tip-speed based; P2 defines `phi` against
mean-radius speed. The conversion is `phi_P2 = 0.407 * U_tip / U_mean = 0.440`.
The report's table gives the measured annulus geometry. For the P2 continuity
input, this run estimates mass flow from the annulus area, NASA flow coefficient,
and standard-air reference density at 288.15 K and 101325 Pa. Those inlet
conditions and the resulting 16.10 kg/s are explicit first-pass approximations;
they are not claimed as the report's matched measured ambient condition and
mass-flow reading.

NASA Table 13 gives third-stage design-point velocity triangles. At 50% span,
the recorded values used here are `Vtheta,in = 10.8 m/s`, `Vtheta,out = 32.6
m/s`, `Vx,in = 29.4 m/s`, and `Vx,out = 28.3 m/s`. Mean blade speed is inferred
from tip speed and the rig's annulus radii.

## Results

| Comparison | NASA / reference | P2 | Difference |
|---|---:|---:|---:|
| Overall specific work | 5.023 kJ/kg, derived from pressure coefficient and 0.90 efficiency | 5.027 kJ/kg | +0.065% (not independent) |
| Stage 3 Euler work at 50% span | 1.287 kJ/kg, from measured triangle | 1.257 kJ/kg | -2.35% |
| Stage 3 tangential-velocity change | 21.8 m/s | 21.287 m/s | -2.35% |
| Stage 3 inlet swirl | 10.8 m/s | 0 m/s | Assumption difference |
| Stage 3 exit swirl | 32.6 m/s | 21.287 m/s | -34.7% |
| Rotor-inlet axial velocity at stage 3 | 29.4 m/s at 50% span | 25.97 m/s mean-line | -11.7% |
| Annulus blade height | 114.3 mm, fixed rig geometry | 114.6 mm at stage 1 inlet; 111.5 mm at stage 4 inlet | +0.3% / -2.5% |

NASA's measured third-stage inlet swirl is substantial, but Euler work depends on
the change in swirl: `U * (Vtheta,out - Vtheta,in)`. The midspan swirl change is
21.8 m/s; P2 obtains a similar 21.287 m/s change by imposing zero inlet swirl
and a lower exit swirl. This supports the user's hypothesis for this test
point: the zero-swirl assumption creates a large velocity-angle difference,
while its effect on this stage's Euler work is much smaller. This is a
midspan-versus-mean-line comparison, not a span-integrated stage validation.

The modeled overall PR is 1.055716, constructed from the four-stage pressure
coefficient and the assumed reference inlet state. P2 solves to that target by
definition. Its four stage PRs are approximately 1.01374, 1.01368, 1.01362,
and 1.01356. Do not count overall PR closure as validation.

## Limitations and next data improvement

- Replace the approximate reference inlet state and mass-flow estimate with
  the exact matched NASA design-point ambient conditions and corrected airflow
  from the report's tabulated/graph data, if available at sufficient precision.
- The efficiency input is the reported **overall** value applied uniformly to
  each modeled stage. NASA's measured third-stage efficiency is not directly
  tabulated as a simple scalar in the material processed for this run.
- Stage 3's measured velocity triangle is spanwise; P2 is a uniform mean-line
  model. P2's equal-stage loading and constant mean radius are simplifications.
- The LSRC uses an IGV and nonzero interstage swirl; P2 does not. Rotor-exit
  swirl and rotor-exit flow angle therefore differ substantially.
- The NASA coefficient-derived work comparison is conditioned on its efficiency
  and pressure coefficient. Use the Euler-triangle and stagewise comparisons as
  the more informative results of this first pass.

## Reproduction

From this directory's parent project folder:

```bash
python validation/nasa_ge_lsrc_first_pass.py
```

The script writes `nasa_ge_lsrc_first_pass_results.json` in this directory. No
core solver equations were changed. The illustrative screening limits were
broadened only in the run assumptions to permit the LSRC's scale and speed;
they were not evaluated as proof of physical feasibility.
