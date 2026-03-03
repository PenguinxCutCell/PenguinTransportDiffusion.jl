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

Tests

Run the package test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Current tests cover:

- no-body Gaussian convergence under periodic flow (`:CN`),
- outside-circle rotating-flow case with mass-drift bound.
