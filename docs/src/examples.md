**Examples**

The repository contains runnable example scripts under `examples/`.

Run from repository root:

```bash
julia --project=. examples/gaussian_nobody_convergence.jl
julia --project=. examples/gaussian_outside_circle.jl
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
- outside-circle rotating-flow case with mass-drift bound.
