# Numerical Algorithms & Routines

This page summarizes the high-level routines implemented by the package.

1) Fixed steady assembly (`assemble_steady_mono!`, `assemble_steady_diph!`)

- Build the diffusion system with `PenguinDiffusion.assemble_steady_*!`.
- Build advection operators at time `t`.
- Add transport blocks (`ω`/`γ` per phase).
- Apply transport box BC rows on `ω` rows.
- Re-apply row-identity constraints for inactive/halo rows.

2) Fixed unsteady assembly (`assemble_unsteady_mono!`, `assemble_unsteady_diph!`)

- Assemble at `t + θΔt`.
- For `θ != 1`, scale `ω` rows by `θ` and add `-(1-θ)A*uⁿ` correction to RHS.
- Add mass matrix on each `ω` block: `M = V / Δt` and `bω += M .* uωⁿ`.
- Re-apply activity constraints and clear solver cache.

3) Moving unsteady assembly (`assemble_unsteady_mono_moving!`, `assemble_unsteady_diph_moving!`)

- Delegate slab geometry + moving diffusion/mass terms to `PenguinDiffusion.assemble_unsteady_*_moving!`.
- Evaluate advection at `τ = t + θΔt`.
- Build moving advection operators on slab-reduced capacities.
- Use relative interface velocity on `γ`: `uγ_rel = uγ - wγ`.
- Append advection blocks, apply transport BCs, enforce row identities.

4) Time marching (`solve_unsteady!`, `solve_unsteady_moving!`)

- Accepts `u0` as either phase `ω` blocks or full system vector.
- Advances with variable last step (`min(dt, tend-t)`).
- Solves each step with `PenguinSolverCore.solve!`.
- Fixed path can reuse constant operators when data are time-independent.
- Moving path reassembles every step.

5) Geometry updates (`rebuild!`)

- Rebuilds capacities (`CartesianOperators.rebuild!`) from new moments.
- Reconstructs diffusion operators/model and periodic flags.
- Refreshes advection operators for requested time.
