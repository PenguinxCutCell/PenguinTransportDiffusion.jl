**API & Types**

Types

- `AdvDiffModelMono{N,T,...}`
  - Fields: `diff`, `uω`, `uγ`, `bc`, `scheme`, `periodic`.
  - Construct with:
    - `AdvDiffModelMono(cap, D, uω, uγ; source=..., bc=..., bc_interface_diff=..., layout=..., coeff_mode=:harmonic, scheme=Centered())`
- `AdvDiffCoupledModelMono`
  - Wrapper/adaptor for `PenguinSolverCore` block coupling.
  - Exposes `:concentration` and accepts incoming `:velocity`.
  - Uses `advance_steady!` / `advance_unsteady!` hooks through SolverCore.

Public functions

- `assemble_steady_mono!(sys::LinearSystem, model::AdvDiffModelMono, t)`
- `assemble_unsteady_mono!(sys::LinearSystem, model::AdvDiffModelMono, uⁿ, t, dt, scheme_or_theta)`
- `solve_steady!(model::AdvDiffModelMono; t=zero(T), method=:direct, kwargs...)`
- `solve_unsteady!(model::AdvDiffModelMono, u0, tspan; dt, scheme=:BE|:CN|θ, method=:direct, save_history=true, kwargs...)`
- `rebuild!(model::AdvDiffModelMono, moments; bc=zero(T), t=zero(T))`
- `update_advection_ops!(model::AdvDiffModelMono; t=zero(T))`

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
- `solve_unsteady!` currently does per-step reassembly (`reused_constant_operator=false`).
