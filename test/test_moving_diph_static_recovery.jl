@testset "Moving diphasic static recovery vs fixed" begin
    xyz = (range(0.0, 1.0; length=81),)
    grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))
    xγ = 0.46
    body1(x, t) = x - xγ

    moms1 = geometric_moments(x -> x - xγ, xyz, Float64, nan; method=:vofijul)
    moms2 = geometric_moments(x -> -(x - xγ), xyz, Float64, nan; method=:vofijul)
    cap1 = assembled_capacity(moms1; bc=0.0)
    cap2 = assembled_capacity(moms2; bc=0.0)

    D1 = 0.11
    D2 = 0.07
    u1 = (0.25,)
    u2 = (-0.15,)
    C = 0.0
    source = (0.0, 0.0)
    bc = BorderConditions(; left=Periodic(), right=Periodic())
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
    u01 = fill(C, nt)
    u02 = fill(C, nt)
    u0 = vcat(u01, u02)

    lay = fixed.diff.layout.offsets

    tspan = (0.0, 0.06)
    dtm = 0.003
    sol_fixed = solve_unsteady!(fixed, u0, tspan; dt=dtm, scheme=:CN, method=:direct, save_history=false)
    sol_moving = solve_unsteady_moving!(moving, u0, tspan; dt=dtm, scheme=:CN, method=:direct, save_history=false)

    xf = sol_fixed.system.x
    xm = sol_moving.system.x
    @test maximum(abs.(xm[lay.ω1] .- xf[lay.ω1])) < 2e-2
    @test maximum(abs.(xm[lay.ω2] .- xf[lay.ω2])) < 2e-2
    @test maximum(abs.(xm .- xf)) < 3e-2
end
