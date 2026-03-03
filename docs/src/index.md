# PenguinTransportDiffusion.jl

Welcome to the PenguinTransportDiffusion documentation. This package couples cut-cell diffusion (`PenguinDiffusion.jl`) and cut-cell transport (`PenguinTransport.jl`) into a mono-phase advection-diffusion workflow.

The continuous model is

```math
\partial_t T + \nabla \cdot (u T) = \nabla \cdot (D \nabla T) + f.
```

Pages

- [Advection-Diffusion Model](transport_diffusion.md)
- [Numerical Algorithms](algorithms.md)
- [API & Types](api.md)
- [Examples & Tests](examples.md)

Build notes

This repository uses Documenter.jl for site generation. To preview locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```
