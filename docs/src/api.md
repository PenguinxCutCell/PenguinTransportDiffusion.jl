# API & Types

Types

- `AdvDiffModelMono{N,T,...}`
- `AdvDiffModelDiph{N,T,...}`
- `MovingAdvDiffModelMono{N,T,...}`
- `MovingAdvDiffModelDiph{N,T,...}`
- `AdvDiffCoupledModelMono`
  - Wrapper/adaptor for `PenguinSolverCore` block coupling.
  - Exposes `:concentration` and accepts incoming `:velocity`.
  - Uses `advance_steady!` / `advance_unsteady!` hooks through SolverCore.

Constructors (fixed geometry)

- `AdvDiffModelMono(cap, D, uω, uγ; source=..., bc=..., bc_interface_diff=..., layout=..., coeff_mode=:harmonic, scheme=Centered())`
- `AdvDiffModelDiph(cap1, D1, u1ω, u1γ, cap2, D2, u2ω, u2γ; source=..., bc=..., ic=..., bc_interface=..., layout=..., coeff_mode=:harmonic, scheme=Centered())`
- `AdvDiffModelDiph(cap, D1, D2, u1ω, u1γ, u2ω, u2γ; ...)` (shared-capacity shorthand)

Constructors (moving geometry)

- `MovingAdvDiffModelMono(grid, body, D, uω, uγ; wγ=..., source=..., bc=..., bc_interface_diff=..., layout=..., coeff_mode=:harmonic, scheme=Centered(), geom_method=:vofijul)`
- `MovingAdvDiffModelDiph(grid, body1, D1, u1ω, u1γ, D2, u2ω, u2γ; wγ=..., source=..., body2=..., bc=..., ic=..., bc_interface=..., layout=..., coeff_mode=:harmonic, scheme=Centered(), geom_method=:vofijul)`
- Composition overloads from pre-built moving diffusion models are also available.

Public functions

- `assemble_steady_mono!(sys::LinearSystem, model::AdvDiffModelMono, t)`
- `assemble_steady_diph!(sys::LinearSystem, model::AdvDiffModelDiph, t)`
- `assemble_unsteady_mono!(sys::LinearSystem, model::AdvDiffModelMono, uⁿ, t, dt, scheme_or_theta)`
- `assemble_unsteady_diph!(sys::LinearSystem, model::AdvDiffModelDiph, uⁿ, t, dt, scheme_or_theta)`
- `assemble_unsteady_mono_moving!(sys::LinearSystem, model::MovingAdvDiffModelMono, uⁿ, t, dt, scheme_or_theta)`
- `assemble_unsteady_diph_moving!(sys::LinearSystem, model::MovingAdvDiffModelDiph, uⁿ, t, dt, scheme_or_theta)`
- `solve_steady!(model::AdvDiffModelMono; t=zero(T), method=:direct, kwargs...)`
- `solve_steady!(model::AdvDiffModelDiph; t=zero(T), method=:direct, kwargs...)`
- `solve_unsteady!(model::AdvDiffModelMono, u0, tspan; dt, scheme=:BE|:CN|θ, method=:direct, save_history=true, kwargs...)`
- `solve_unsteady!(model::AdvDiffModelDiph, u0, tspan; dt, scheme=:BE|:CN|θ, method=:direct, save_history=true, kwargs...)`
- `solve_unsteady_moving!(model::MovingAdvDiffModelMono, u0, tspan; dt, scheme=:BE|:CN|θ, method=:direct, save_history=true, kwargs...)`
- `solve_unsteady_moving!(model::MovingAdvDiffModelDiph, u0, tspan; dt, scheme=:BE|:CN|θ, method=:direct, save_history=true, kwargs...)`
- `rebuild!(model::AdvDiffModelMono, moments; bc=zero(T), t=zero(T))`
- `rebuild!(model::AdvDiffModelDiph, moments; bc=zero(T), t=zero(T))`
- `update_advection_ops!(model::AdvDiffModelMono; t=zero(T))`
- `update_advection_ops!(model::AdvDiffModelDiph; t=zero(T))`

Inputs and callbacks

- `D` can be scalar or callback (delegated to `PenguinDiffusion`).
- Velocity inputs (`uω`, `uγ`) accept tuple/vector with `N` components.
- Per-component velocity can be:
  - scalar constant,
  - vector of nodal values (`length == cap.ntotal`),
  - callback `(x...)` or `(x..., t)`,
  - `Ref` to any of the above.

Boundary semantics

- `bc` accepts diffusion and transport boundary types in one object.
- Internally split by side:
  - diffusion side: `Dirichlet`, `Neumann`, `Periodic`.
  - advection side: `Inflow`, `Outflow`, `Periodic`.
- `Dirichlet`/`Neumann` on `bc` imply advection `Outflow` at that side.
- `Inflow`/`Outflow` on `bc` imply diffusion homogeneous `Neumann` at that side.

Notes

- Unknown layout follows `layout_mono(cap.ntotal)` unless overridden.
- For fixed geometry, `solve_unsteady!` can reuse a constant assembled operator when data are time-independent.
- For moving geometry, `solve_unsteady_moving!` reassembles each step.
