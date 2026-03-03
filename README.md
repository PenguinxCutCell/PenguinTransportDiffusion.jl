# PenguinTransportDiffusion.jl

Thin advection-diffusion package for cut-cell Cartesian grids:

```
∂t T + ∇·(uT) = ∇·(D∇T) + f
```

Design:
- diffusion assembly delegated to `PenguinDiffusion`
- advection operators and inflow/outflow boundary injection delegated to `CartesianOperators` / `PenguinTransport`
- one mono model (`ω,γ` layout) with unsteady support (`:BE`, `:CN`, or numeric `θ`).

Boundary condition note:
- one `BorderConditions` is accepted by `AdvDiffModelMono`.