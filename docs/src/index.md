# PenguinTransportDiffusion.jl

`PenguinTransportDiffusion.jl` solves cut-cell advection-diffusion on Cartesian grids by composing:

- `PenguinDiffusion.jl` for diffusion assembly (fixed and moving geometry),
- `PenguinTransport.jl` for advection operators and transport-side BC closures.

Supported model families:

- mono and diphasic,
- fixed and moving embedded interfaces,
- steady and unsteady (`θ`) formulations.

Continuous PDE:

```math
\partial_t T + \nabla \cdot (u T) = \nabla \cdot (D \nabla T) + f.
```

## Start Here

- [Advection-Diffusion Model](transport_diffusion.md)
- [Moving Geometry](moving.md)
- [Numerical Algorithms](algorithms.md)
- [API & Types](api.md)
- [Examples & Tests](examples.md)

## Feature Matrix (Summary)

| Area | Support |
|---|---|
| Fixed mono steady/unsteady | Yes |
| Fixed diph steady/unsteady | Yes |
| Moving mono unsteady | Yes |
| Moving diph unsteady | Yes |
| Time schemes | `:BE`, `:CN`, numeric `θ ∈ [0,1]` |
| Outer BCs (combined API) | Diffusion: Dirichlet/Neumann/Robin/Periodic, Advection: Inflow/Outflow/Periodic |
| Interface laws | Diffusion side (`Robin` mono, `InterfaceConditions` diph) |

## Build notes

This repository uses Documenter.jl for site generation. To preview locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```
