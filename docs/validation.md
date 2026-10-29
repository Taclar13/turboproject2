# Project 2 verification record

## Verification scope

This record checks the implemented one-dimensional mean-line equations and the
program's own selected screening rules. It does not validate aerodynamic
performance, structural adequacy, manufacturing feasibility, or mass against a
physical compressor. Stage efficiency, gas properties, thresholds, and
hardware-mass assumptions remain prescribed or illustrative. The rotor is
assumed lossless, and all modeled stage entropy generation is assigned to the
stator; there is no adjustable rotor/stator loss split.

## Independent equation checks

The test suite independently recomputes the following from the documented
constant-property assumptions:

- **Single-stage thermodynamics:** for a compressor with prescribed stage
  efficiency, the pressure-ratio relation is
  `PR = [1 + eta_stage * delta_h0 / (cp * T0_in)]^(gamma/(gamma-1))`.
  A one-stage PR 1.265 case is checked against the algebraic inverse for
  `delta_h0`.
- **Euler work:** `delta_h0 = U * (Vtheta_out - Vtheta_in)` at each modeled
  rotor, and `cp * (T0_out - T0_in) = delta_h0`.
- **Pressure-ratio closure:** the product of stage pressure ratios equals the
  requested overall stagnation-pressure ratio.
- **Continuity and annulus geometry:** `mdot = rho * A * Vx` and
  `A = pi * (rt^2 - rh^2)` at every exported station.
- **Entropy accounting:** rotor entropy generation is zero under the lossless
  rotor assumption, and stator entropy generation equals the modeled stage
  generation.
- **Input and screening behavior:** normalized user-selected loading profiles,
  invalid input rejection, a documented feasible example, and an explicitly
  reported failing default case.

## Reproduced examples

On 2026-10-27, the documented 10-stage example (16 kg/s, 9,000 rpm, PR 5.0,
mean radius 0.23 to 0.21 m, `phi = 0.55`) achieved PR 5.000000, an outlet total
temperature of 514.16 K, 215.23 kJ/kg modeled specific work, and 3.444 MW
modeled shaft power. It passed all 326 implemented checks. Its 0.595 m length,
0.544 m maximum diameter, and 62.85 kg mass are outputs of the illustrative
geometry and component-mass proxies, not hardware predictions.

The untouched default (16 stages, 16 kg/s, 10,500 rpm, PR 14.0) achieved the
requested pressure ratio but failed 35 of 518 implemented checks, mainly
minimum blade-height checks. The solver reports these failures without
silently changing the selected design inputs.

The accompanying test suite verifies the current equations and code behavior;
passing tests are not equivalent to experimental or CFD validation.

## NASA/GE LSRC first-pass benchmark

The completed P2 benchmark is documented in [`validation/README.md`](../validation/README.md)
with its reproducible setup and output JSON. P2's overall work is within 0.065%
of the NASA coefficient-derived reference, but that is a thermodynamic
consistency result rather than independent validation because the NASA
efficiency and pressure coefficient define the comparison target.

The more informative comparison is the third-stage, 50%-span Euler-work check:
P2 is 2.35% below the work implied by NASA's measured velocity triangle. P2
uses zero inlet swirl and differs substantially in absolute exit swirl, while
the tangential-velocity change is close at this point. The test therefore
supports the implemented mean-line/Euler calculation for this selected case
under documented assumptions. It does not establish general aerodynamic
accuracy or validate the model's efficiency prediction.

This limited comparison is sufficient for P2's learning scope. Additional
aerodynamic model development and more rigorous compressor validation are
deferred to Project 2.5.
