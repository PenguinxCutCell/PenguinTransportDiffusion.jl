**Advection-Diffusion Model**

`PenguinTransportDiffusion.jl` provides a mono-phase cut-cell advection-diffusion model through `AdvDiffModelMono`.

Continuous PDE

```math
\partial_t T + \nabla \cdot (u T) = \nabla \cdot (D \nabla T) + f.
```

Discretization split

- Diffusion assembly is delegated to `PenguinDiffusion` (`assemble_steady_mono!` on the diffusion model).
- Advection operators and advection-side box BC injection are delegated to `PenguinTransport`.
- The package combines both contributions into one linear system over the `(\omega, \gamma)` unknown layout.

Boundary handling

A single `BorderConditions` input is accepted by `AdvDiffModelMono`, then internally split:

- `Dirichlet` / `Neumann` / `Periodic` are applied on the diffusion side.
- `Inflow` / `Outflow` / `Periodic` are applied on the advection side.
- Mixed unsupported boundary kinds raise an `ArgumentError`.

Interface condition

- Diffusion interface coupling is optional via `bc_interface_diff::Union{Nothing, PenguinBCs.Robin}`.
- A pure Neumann interface can introduce a `\gamma`-gauge nullspace; one `\gamma` row may need pinning in custom solves.

Time integration

- Supports `:BE`, `:CN`, or numeric `\theta` in `assemble_unsteady_mono!` / `solve_unsteady!`.
- For `\theta \neq 1`, the explicit correction `(1-\theta) A u^n` is moved to the RHS.
- Mass contribution uses cell volumes `V / \Delta t` on the `\omega` block.
