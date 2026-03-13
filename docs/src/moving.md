# Moving Geometry (Space-Time Slabs)

Moving advection-diffusion models:

- `MovingAdvDiffModelMono`
- `MovingAdvDiffModelDiph`

compose moving diffusion from `PenguinDiffusion` with moving advection from `PenguinTransport`.

## Assembly Composition

Per step on `[t^n, t^{n+1}]`:

1. `PenguinDiffusion.assemble_unsteady_*_moving!` builds slab geometry terms, moving mass/GCL part, and diffusion/interface contributions.
2. `PenguinTransportDiffusion` evaluates advection at `t^n + \theta \Delta t`.
3. Moving advection blocks are appended:
   - `conv_bulk = sum(opsA.C)`
   - `conv_iface = 0.5 * sum(opsA.K)`
4. Outer transport BCs are applied.
5. Active-row identity constraints are re-applied.

## Relative Interface Velocity

For moving interfaces, advection uses:

```math
u_{\gamma,\mathrm{rel}} = u_\gamma - w_\gamma
```

where `wγ` is the interface velocity sampled at `γ` points.

## API

```julia
sol = solve_unsteady_moving!(model, u0, (t0, tf); dt=Δt, scheme=:BE, method=:direct)
```

Accepted `scheme` values:

- `:BE` (`θ=1`)
- `:CN` (`θ=1/2`)
- numeric `θ` with `0 ≤ θ ≤ 1`

## Regression Coverage

- Static-geometry recovery (`wγ=0`, time-independent body) against fixed models.
- Zero-advection reduction against moving diffusion.
- MMS convergence for moving mono (no interface) and moving diphasic (translating planar interface).
