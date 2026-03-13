**Advection-Diffusion Model**

`PenguinTransportDiffusion.jl` provides fixed and moving advection-diffusion models:

- `AdvDiffModelMono`, `AdvDiffModelDiph`
- `MovingAdvDiffModelMono`, `MovingAdvDiffModelDiph`

Continuous PDE

```math
\partial_t T + \nabla \cdot (u T) = \nabla \cdot (D \nabla T) + f.
```

Discretization split

- Diffusion + mass terms are delegated to `PenguinDiffusion`.
- Advection operators and advection-side box BC injection are delegated to `PenguinTransport`.
- This package composes both into one linear system over `(\omega, \gamma)` (mono) or `(\omega_1,\gamma_1,\omega_2,\gamma_2)` (diph).

Boundary handling

A single `BorderConditions` input is accepted by all models, then internally split:

- `Dirichlet` / `Neumann` / `Periodic` are applied on the diffusion side.
- `Inflow` / `Outflow` / `Periodic` are applied on the advection side.
- Mixed unsupported boundary kinds raise an `ArgumentError`.

Interface condition

- Mono diffusion interface coupling is optional via `bc_interface_diff::Union{Nothing, PenguinBCs.Robin}`.
- Diph interface coupling is passed through `ic::Union{Nothing,InterfaceConditions}`.
- A pure Neumann interface can introduce a `\gamma`-gauge nullspace; one `\gamma` row may need pinning in custom solves.

Time integration

- Fixed geometry: `assemble_unsteady_*` / `solve_unsteady!` support `:BE`, `:CN`, numeric `\theta`.
- Moving geometry: `assemble_unsteady_*_moving!` / `solve_unsteady_moving!` support the same values.
- For `\theta \neq 1`, the explicit correction `(1-\theta) A u^n` is moved to the RHS.
- Fixed mass contribution uses cell volumes `V / \Delta t` on `\omega`.
- Moving assembly uses slab-reduced capacities from `PenguinDiffusion` and advection with relative interface speed `(u_\gamma - w_\gamma) \cdot n_\gamma`.
