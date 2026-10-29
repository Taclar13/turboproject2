# Historical Project 2 optimizer study

> **Archived result from an earlier model version.** This report describes the optimizer-based code that preceded the current manually configurable, one-design-at-a-time solver. The current solver does not reproduce these optimizer or mass-minimization results. Treat the results below as a historical record, not as verification of the current code.

## Study setup

For **each required design mass flow** (8, 16, and 24 kg/s), independently design a new multistage axial compressor at 300 K and 101,325 Pa inlet stagnation conditions to produce an overall **stagnation pressure ratio of 2.2**. The compressor mass is minimized **within the explicitly simplified component mass model** and screened among integer stage counts 3–8. Each stage receives the **assumed** stage efficiency 0.85. All other engineering constants and limits are listed in `results/assumptions.json`.

## Results from the shipped seeded example run

| Required flow (kg/s) | Best feasible stage count found | Estimated modeled mass (kg) | Flowpath length (m) | Maximum outer-flowpath diameter (m) | Overall PR |
|---:|---:|---:|---:|---:|---:|
| 8 | 3 | 8.06 | 0.117 | 0.328 | 2.20000 |
| 16 | 3 | 15.38 | 0.190 | 0.429 | 2.20000 |
| 24 | 3 | 23.15 | 0.226 | 0.534 | 2.20000 |

**Diameter** is twice the maximum local *flowpath tip* radius; it excludes casing-wall and external stator-mount thickness. **Mass** includes the modeled casing, rotor drum, shaft core, rotor discs, stator mounting bands, and rotor/stator blades, but excludes many real compressor assembly components (see README).

Within the searched stage counts and illustrative assumptions, the lowest modeled mass in these three cases was found at the **minimum allowed count of three**. This does *not* establish that three stages minimizes actual compressor mass: omitted aero/structural constraints and assumed fixed efficiency materially affect this tradeoff. Designs at 4–8 stages were also evaluated; their lowest feasible modeled masses and optimizer status are recorded in `results/attempts.json`.

**Validation and search status:** the seeded run's differential-evolution searches reported numerical convergence for all 18 mass-flow/stage-count searches. Independently scripted checks of saved designs pass continuity, stagnation/static energy consistency, ideal-gas density, Euler work, stage efficiency/pressure-ratio consistency, component-mass sum, station ordering, and recorded constraint margins. The test suite reports 11 passing tests. Numerical convergence does not certify a global optimum, and passing these model checks is not experimental or CFD validation.

## Essential model limitations

The blade-row aerodynamics are *screened* using mean-line and simplified tip Mach numbers, diffusion factor and reaction, rather than predicted from actual blade profiles. There is no blade stress/rotordynamics, stall/surge/choke margin, spanwise throughflow, geometry-derived loss correlation, or calibrated structural mass sizing. Material choice and thicknesses, blade aspect/solidity/stagger/shape factors, and loss allocation are illustrative and should not be treated as manufacturing-ready design data. These numerical results are demonstrations of the **implemented preliminary model**, not performance or weight claims for real machines.

## Reproducibility

```bash
python -m pip install -r requirements.txt
python compressor.py --flows 8 16 24 --pressure-ratio 2.2 --stages 3 4 5 6 7 8 --maxiter 65 --popsize 7 --seed 21 --output results
pytest -q test_compressor.py
python verify_results.py results/designs.json
```
