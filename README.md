# IEID vs ADRC: Nonlinear PWR Power Control Comparison

This project provides a comparative simulation framework between **Improved Equivalent Input Disturbance (IEID)** and **Active Disturbance Rejection Control (ADRC)** on a nonlinear model of a Pressurized Water Reactor (PWR).

## Project Structure

The codebase is organized into modular folders:

- **`plant/`**: Contains the nonlinear reactor dynamics and reactivity equations.
  - `pwr_nonlinear_dynamics.m`: Implementation of the 4-state nonlinear model.
  - `total_reactivity.m`: Temperature feedback and rod reactivity sum.
- **`ieid/`**: Implementation of the IEID controller.
  - `run_ieid_closed_loop.m`: ODE wrapper for IEID.
  - `ieid_estimator_fcn.m`: Disturbance estimation logic.
- **`adrc/`**: Implementation of the ADRC controller.
  - `run_adrc_closed_loop.m`: ODE wrapper for ADRC.
  - `eso_dynamics.m`: Extended State Observer (ESO) dynamics.
- **`shared/`**: Common utilities and parameters.
  - `get_nonlinear_params.m`: Centralized parameters for the new plant and controllers.
  - `disturbance_fcn.m`: Shared disturbance signal.
  - `plotting.m`: Comparison visualization.
  - `metrics.m`: Performance metric calculation (Rise Time, Overshoot, etc.).
- **`main_compare_ieid_adrc.m`**: The main entry point to run the comparison study.

## How to Run

1. Open MATLAB.
2. Navigate to the `e:\BTP\END SEM\` directory.
3. Execute the comparison script:
   ```matlab
   main_compare_ieid_adrc
   ```

## Comparison Methodology

Both controllers are evaluated under identical conditions:
- **Plant**: Nonlinear 4-state kinetics + Rod reactivity dynamics (5 states total).
- **Reference**: Step change from 0.5 to 0.6 power at $t=20$s.
- **Disturbance**: Complex signal including sinusoidal oscillations and transitions.
- **Solver**: `ode15s` with identical tolerances.

## Metrics Evaluated
The simulation automatically calculates and displays:
- RMS Tracking Error
- Rise Time (10% to 90%)
- Settling Time (2% band)
- Percentage Overshoot
- Steady-State Error
- Control Energy
