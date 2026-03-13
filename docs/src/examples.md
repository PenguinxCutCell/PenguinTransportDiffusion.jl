**Examples**

The repository contains runnable example scripts under `examples/`.

Run from repository root:

```bash
julia --project=. examples/gaussian_nobody_convergence.jl
julia --project=. examples/gaussian_outside_circle.jl
julia --project=. examples/moving_mono_nobody_mms_convergence.jl
julia --project=. examples/moving_diph_planar_interface_mms_convergence.jl
```

Available examples

- `examples/gaussian_nobody_convergence.jl`
  - 2D no-body (full domain), periodic boundaries.
  - Advects and diffuses a Gaussian with known analytical solution.
  - Uses `solve_unsteady!(...; scheme=:CN)` and reports observed convergence order.

- `examples/gaussian_outside_circle.jl`
  - 2D outside a circular embedded body with rotational velocity field.
  - Periodic box + interface Robin condition `Robin(0, 1, 0)`.
  - Advances with backward Euler and tracks total scalar mass drift.
  - Writes `gaussian_outside_circle_mass.csv`.

- `examples/embedded_wall_hot_cylinder.jl`
  - Embedded fixed wall with advection-diffusion coupling.
  - Demonstrates mixed model BC handling and interface diffusion law.

- `examples/diph_planar_interface_constant_jump.jl`
  - Fixed diphasic planar interface.
  - Uses `InterfaceConditions` for scalar/flux coupling.

- `examples/diph_planar_interface_mms.jl`
  - Fixed diphasic MMS benchmark.
  - Reports refinement errors and observed order.

- `examples/moving_mono_nobody_mms_convergence.jl`
  - Moving-model path without embedded interface (`body=-1`).
  - Verifies moving mono assembly convergence in periodic box.

- `examples/moving_diph_planar_interface_mms_convergence.jl`
  - Moving diphasic planar interface MMS.
  - Verifies robust convergence with translating interface.

- `examples/moving_static_recovery_compare.jl`
  - Compares fixed vs moving (with static interface, `wγ=0`) final states.
  - Prints `L2`, `L∞`, and max-abs differences.

- `examples/coupled_darcy_advdiff_oneway_1d.jl`
  - First real one-way coupling: Darcy velocity drives advection-diffusion.
  - Uses `OneWayCoupling` with `CouplingMap(:darcy, :transport, :velocity)`.
  - Compares coupled result to direct velocity injection reference.

- `examples/coupled_darcy_advdiff_twoway_1d.jl`
  - First real two-way coupling: concentration feeds back to Darcy mobility.
  - Uses `TwoWayCoupling` with outer Picard iterations and relaxation.
  - Includes map `CouplingMap(:transport, :darcy, :concentration)`.

First coupled workflow notes

- One-way: flow drives scalar transport.
- Two-way: scalar updates mobility `λ(c)`, which updates flow, which updates scalar.
- SolverCore handles orchestration only; wrappers expose coupling fields.
- Current Darcy unsteady hook is quasi-steady per transport time step.

Tests

Run the package test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Current tests cover:

- no-body Gaussian convergence under periodic flow (`:CN`),
- outside-circle rotating-flow case with mass-drift bound,
- moving static-geometry recovery vs fixed (`mono` and `diph`),
- moving zero-advection reduction to `PenguinDiffusion` (`mono` and `diph`),
- moving MMS checks (`mono` and `diph`).

MMS matrix sweep:

```bash
julia --project=. test/test_mms_matrix_report.jl
```

Validated orders (2026-03-13):

| kind | scheme | embedded | moving | p12 | p23 |
|---|---|---|---|---:|---:|
| mono | CN | false | false | 1.96 | 1.99 |
| mono | CN | false | true | 0.99 | 2.99 |
| mono | CN | true | false | 1.96 | 1.99 |
| mono | CN | true | true | 0.99 | 2.99 |
| mono | BE | false | false | 1.96 | 1.99 |
| mono | BE | false | true | 1.99 | 1.99 |
| mono | BE | true | false | 1.96 | 1.99 |
| mono | BE | true | true | 1.99 | 1.99 |
| diph | CN | false | false | 1.96 | 1.99 |
| diph | CN | false | true | 0.99 | 2.99 |
| diph | CN | true | false | 1.96 | 1.99 |
| diph | CN | true | true | 0.99 | 2.99 |
| diph | BE | false | false | 1.96 | 1.99 |
| diph | BE | false | true | 1.99 | 1.99 |
| diph | BE | true | false | 1.96 | 1.99 |
| diph | BE | true | true | 1.99 | 1.99 |
