using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinTransportDiffusion
import PenguinDiffusion

function report_diff(label, d)
    l2 = norm(d)
    linf = norm(d, Inf)
    maxabs = maximum(abs.(d))
    println("  $label: L2=$l2, Linf=$linf, maxabs=$maxabs")
end

function compare_mono()
    xyz = (range(0.0, 1.0; length=65),)
    grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))
    body(x, t) = -1.0

    moms = geometric_moments(x -> -1.0, xyz, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)

    D = 0.06
    u = (0.3,)
    source(x, t) = 0.1 * sin(2pi * x) * exp(0.2 * t)
    bc = BorderConditions(; left=Periodic(), right=Periodic())

    fixed = AdvDiffModelMono(cap, D, u, u; source=source, bc=bc, scheme=Centered())
    moving = MovingAdvDiffModelMono(grid, body, D, u, u; wγ=(0.0,), source=source, bc=bc, scheme=Centered())

    u0 = [sin(2pi * cap.C_ω[i][1]) for i in 1:cap.ntotal]
    tspan = (0.0, 0.06)
    dt = 0.003
    sol_fixed = solve_unsteady!(fixed, u0, tspan; dt=dt, scheme=:CN, method=:direct, save_history=false)
    sol_moving = solve_unsteady_moving!(moving, u0, tspan; dt=dt, scheme=:CN, method=:direct, save_history=false)

    lay = fixed.diff.layout.offsets
    dx_full = sol_moving.system.x .- sol_fixed.system.x
    dx_omega = dx_full[lay.ω]
    println("Mono static-recovery difference (moving - fixed):")
    report_diff("omega", dx_omega)
    report_diff("full", dx_full)
end

function compare_diph()
    xyz = (range(0.0, 1.0; length=81),)
    grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))
    xγ = 0.46
    body1(x, t) = x - xγ

    moms1 = geometric_moments(x -> x - xγ, xyz, Float64, nan; method=:vofijul)
    moms2 = geometric_moments(x -> -(x - xγ), xyz, Float64, nan; method=:vofijul)
    cap1 = assembled_capacity(moms1; bc=0.0)
    cap2 = assembled_capacity(moms2; bc=0.0)

    D1 = 0.09
    D2 = 0.05
    u1 = (0.2,)
    u2 = (-0.1,)
    source = (
        (x, t) -> 0.12 * exp(-0.1 * t) * cos(2pi * x),
        (x, t) -> -0.08 * exp(-0.2 * t) * sin(2pi * x),
    )
    bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
    ic = InterfaceConditions(; scalar=ScalarJump(1.0, 1.0, 0.0), flux=FluxJump(1.0, 1.0, 0.0))

    fixed = AdvDiffModelDiph(
        cap1,
        D1,
        u1,
        u1,
        cap2,
        D2,
        u2,
        u2;
        source=source,
        bc=bc,
        ic=ic,
        scheme=Centered(),
    )
    moving = MovingAdvDiffModelDiph(
        grid,
        body1,
        D1,
        u1,
        u1,
        D2,
        u2,
        u2;
        wγ=(0.0,),
        source=source,
        bc=bc,
        ic=ic,
        scheme=Centered(),
    )

    nt = cap1.ntotal
    u01 = [0.8 + 0.2 * sin(2pi * cap1.C_ω[i][1]) for i in 1:nt]
    u02 = [0.3 + 0.1 * cos(2pi * cap2.C_ω[i][1]) for i in 1:nt]
    u0 = vcat(u01, u02)

    tspan = (0.0, 0.05)
    dt = 0.0025
    sol_fixed = solve_unsteady!(fixed, u0, tspan; dt=dt, scheme=:CN, method=:direct, save_history=false)
    sol_moving = solve_unsteady_moving!(moving, u0, tspan; dt=dt, scheme=:CN, method=:direct, save_history=false)

    lay = fixed.diff.layout.offsets
    dx_full = sol_moving.system.x .- sol_fixed.system.x
    dx_omega = vcat(dx_full[lay.ω1], dx_full[lay.ω2])
    println("Diphasic static-recovery difference (moving - fixed):")
    report_diff("omega1+omega2", dx_omega)
    report_diff("full", dx_full)
end

compare_mono()
compare_diph()
