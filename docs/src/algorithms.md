**Numerical Algorithms & Routines**

This page summarizes the high-level routines implemented by the package.

1) Steady assembly (`assemble_steady_mono!`)

- Build the diffusion system with `PenguinDiffusion.assemble_steady_mono!`.
- Build advection operators at time `t` from `uω`, `uγ`.
- Add transport blocks:
  - `Aωω += Σ(C_d) + 0.5 Σ(K_d)`
  - `Aωγ += 0.5 Σ(K_d)`
- Apply transport box BC rows using `PenguinTransport.apply_box_bc_transport_mono!`.
- Re-apply row-identity constraints for inactive/halo rows.

2) Unsteady assembly (`assemble_unsteady_mono!`)

- Assemble at `t + θΔt`.
- For `θ != 1`, scale `ω` rows by `θ` and add `-(1-θ)A*uⁿ` correction to RHS.
- Add mass matrix on `ω`: `M = V / Δt` and `bω += M .* uωⁿ`.
- Re-apply activity constraints and clear solver cache.

3) Time marching (`solve_unsteady!`)

- Accepts `u0` as either `ω`-only vector (`ntotal`) or full system vector.
- Advances with variable last step (`min(dt, tend-t)`).
- Solves each step with `PenguinSolverCore.solve!`.
- Returns `(times, states, system, reused_constant_operator=false)`.

4) Geometry updates (`rebuild!`)

- Rebuilds capacities (`CartesianOperators.rebuild!`) from new moments.
- Reconstructs diffusion operators/model and periodic flags.
- Refreshes advection operators for requested time.
