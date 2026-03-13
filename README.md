# PenguinTransportDiffusion.jl

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PenguinxCutCell.github.io/PenguinTransportDiffusion.jl/dev)
![CI](https://github.com/PenguinxCutCell/PenguinTransportDiffusion.jl/actions/workflows/ci.yml/badge.svg)
![Coverage](https://codecov.io/gh/PenguinxCutCell/PenguinTransportDiffusion.jl/branch/main/graph/badge.svg)


`PenguinTransportDiffusion.jl` is a cut-cell advection-diffusion package built on top of `PenguinDiffusion.jl` and `PenguinTransport.jl`.
It also provides a SolverCore coupling adapter (`AdvDiffCoupledModelMono`) for Darcy-to-transport workflows.

It solves

```math
\partial_t T + \nabla \cdot (uT) = \nabla \cdot (D\nabla T) + f
```

for mono-phase and diphasic problems on Cartesian cut-cell grids using coupled `(ω, γ)` formulations.

## Documentation

- [Home](https://PenguinxCutCell.github.io/PenguinTransportDiffusion.jl/dev/)
- [Advection-Diffusion model](https://PenguinxCutCell.github.io/PenguinTransportDiffusion.jl/dev/transport_diffusion/)
- [Moving geometry](https://PenguinxCutCell.github.io/PenguinTransportDiffusion.jl/dev/moving/)
- [Algorithms](https://PenguinxCutCell.github.io/PenguinTransportDiffusion.jl/dev/algorithms/)
- [API reference](https://PenguinxCutCell.github.io/PenguinTransportDiffusion.jl/dev/api/)
- [Examples guide](https://PenguinxCutCell.github.io/PenguinTransportDiffusion.jl/dev/examples/)

## Install

```julia
using Pkg
Pkg.add(url="https://github.com/PenguinxCutCell/PenguinTransportDiffusion.jl")
```

## Quick Start

```julia
using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinTransportDiffusion

grid = (range(0.0, 1.0; length=65), range(0.0, 1.0; length=65))
moms = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)
cap = assembled_capacity(moms; bc=0.0)

bc = BorderConditions(; left=Periodic(), right=Periodic(), bottom=Periodic(), top=Periodic())
model = AdvDiffModelMono(cap, 1e-3, (1.0, 0.0), (1.0, 0.0); source=0.0, bc=bc)

u0 = zeros(cap.ntotal)
sol = solve_unsteady!(model, u0, (0.0, 0.01); dt=1e-4, scheme=:BE, method=:direct)
```

## Time Schemes

- Fixed geometry unsteady solves: `scheme=:BE`, `scheme=:CN`, or numeric `θ` with `0 ≤ θ ≤ 1`.
- Moving geometry unsteady solves: same accepted scheme values via `solve_unsteady_moving!`.
- `:BE` means `θ=1`; `:CN` means `θ=1/2`.

## Dependency / Environment Notes

- Root runtime dependencies are in `Project.toml`.
- Docs environment is in `docs/Project.toml`.
- Example scripts run from root project (`julia --project=.`).

## Examples

- `examples/gaussian_nobody_convergence.jl`
- `examples/gaussian_outside_circle.jl`
- `examples/embedded_wall_hot_cylinder.jl`
- `examples/diph_planar_interface_constant_jump.jl`
- `examples/diph_planar_interface_mms.jl`
- `examples/moving_mono_nobody_mms_convergence.jl`
- `examples/moving_diph_planar_interface_mms_convergence.jl`
- `examples/moving_static_recovery_compare.jl`
- `examples/coupled_darcy_advdiff_oneway_1d.jl`
- `examples/coupled_darcy_advdiff_twoway_1d.jl`

## First Real Coupled Workflow: Darcy + Advection-Diffusion

This first real multiphysics workflow uses:

- `DarcyCoupledModelMono` from `PenguinDarcy.jl` as producer of `:velocity`
- `AdvDiffCoupledModelMono` from `PenguinTransportDiffusion.jl` as consumer of `:velocity` and producer of `:concentration`
- Solver orchestration from `PenguinSolverCore.jl`

One-way (`Darcy -> transport`):

```julia
using PenguinDarcy, PenguinTransportDiffusion, PenguinSolverCore

darcy_block = CoupledBlock(:darcy, DarcyCoupledModelMono(darcy_model); init=nothing, cache=Dict{Symbol,Any}())
transport_block = CoupledBlock(:transport, AdvDiffCoupledModelMono(advdiff_model); init=(concentration=c0,), cache=Dict{Symbol,Any}())

problem = CoupledProblem(
    [darcy_block, transport_block],
    OneWayCoupling([:darcy, :transport]);
    maps=[CouplingMap(:darcy, :transport, :velocity)],
)

step_coupled!(problem, t, dt; method=:direct)
```

Two-way (Picard feedback with mobility law `λ(c)`):

```julia
mode = TwoWayCoupling(
    [:darcy, :transport];
    maxiter=30,
    atol=1e-10,
    rtol=1e-8,
    relaxation=Dict(:darcy => 0.7, :transport => 1.0),
    norm_fields=[:velocity, :concentration],
)

problem = CoupledProblem(
    [darcy_block, transport_block],
    mode;
    maps=[
        CouplingMap(:darcy, :transport, :velocity),
        CouplingMap(:transport, :darcy, :concentration),
    ],
)

step_coupled!(problem, t, dt; method=:direct)
```

Current Darcy coupling in this workflow is **quasi-steady per transport time step**.

## Feature Release Matrix

| Area | Item | Status | Notes |
|---|---|---|---|
| Models | Mono-phase steady advection-diffusion | Implemented | `AdvDiffModelMono` + `assemble_steady_mono!` |
| Models | Mono-phase unsteady advection-diffusion | Implemented | `assemble_unsteady_mono!` + `solve_unsteady!` (`:BE`, `:CN`, numeric `θ`) |
| Models | Diphasic advection-diffusion | Implemented | `AdvDiffModelDiph` + `assemble_steady_diph!` / `assemble_unsteady_diph!` |
| Models | Moving mono-phase advection-diffusion | Implemented | `MovingAdvDiffModelMono` + `assemble_unsteady_mono_moving!` + `solve_unsteady_moving!` |
| Models | Moving diphasic advection-diffusion | Implemented | `MovingAdvDiffModelDiph` + `assemble_unsteady_diph_moving!` + `solve_unsteady_moving!` |
| Coupling | Diffusion assembly reuse | Implemented | Delegates to `PenguinDiffusion.assemble_steady_mono!` |
| Coupling | Advection operator reuse | Implemented | Delegates advection ops/BCs to `PenguinTransport` |
| Coefficients | Constant diffusion coefficient | Implemented | Scalar `D` supported through diffusion model |
| Coefficients | Space/time variable diffusion coefficient | Implemented | Callback coefficients supported through diffusion model |
| Velocity input | Constant velocity components | Implemented | Scalar per component supported |
| Velocity input | Nodal vector velocity components | Implemented | Per-component vector of length `ntotal` supported |
| Velocity input | Callback velocity components | Implemented | `(x...)` and `(x..., t)` supported |
| Interface | Diffusion interface Robin condition | Implemented | `bc_interface_diff::Union{Nothing,Robin}` |
| Interface | Pure Neumann interface gauge fixing helper | Partial | Nullspace exists; user pins one active `γ` row in examples/tests |
| Outer BCs | Periodic | Implemented | Supported on diffusion and advection sides |
| Outer BCs | Dirichlet/Neumann in mixed model | Implemented | Diffusion uses given BC; advection side mapped to `Outflow` |
| Outer BCs | Inflow/Outflow in mixed model | Implemented | Advection uses given BC; diffusion side mapped to homogeneous Neumann |
| Time integration | Backward Euler (`:BE`) | Implemented | Via theta path (`θ=1`) |
| Time integration | Crank-Nicolson (`:CN`) | Implemented | Via theta path (`θ=1/2`) |
| Time integration | Generic theta (`θ`) | Implemented | Numeric `scheme` accepted |
| Solver path | Direct / iterative backend use | Implemented | Uses `PenguinSolverCore.solve!` with `method` keyword |
| Solver path | Constant-operator/factorization reuse in march | Implemented | `solve_unsteady!` reuses matrix/factorization when coefficients and BC data are time-independent |
| Utilities | Geometry rebuild | Implemented | `rebuild!(model, moments; ...)` |
| Utilities | Advection-op refresh API | Implemented | `update_advection_ops!` |
| Validation | Regression tests | Implemented | Convergence and mass-drift tests in `test/runtests.jl` |
| Validation | Static-geometry recovery (moving vs fixed) | Implemented | `test/test_moving_*_static_recovery.jl` |
| Validation | Zero-advection reduction to moving diffusion | Implemented | `test/test_moving_*_zero_advection_matches_diffusion.jl` |
| Validation | Moving MMS convergence examples | Implemented | `examples/moving_*_mms_convergence.jl` |
| Validation | Standalone examples | Implemented | Ten runnable scripts in `examples/` |

## MMS Matrix (Validated)

Validated on **2026-03-13** with:

```bash
julia --project=. test/test_mms_matrix_report.jl
```

`p12` and `p23` are successive refinement orders.

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

## Local Docs Build

Local build:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

## Development

Run tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
