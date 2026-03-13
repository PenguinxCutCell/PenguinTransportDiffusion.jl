@testset "Moving mono static recovery vs fixed" begin
    xyz = (range(0.0, 1.0; length=65),)
    grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))
    body(x, t) = -1.0

    moms = geometric_moments(x -> -1.0, xyz, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)

    D = 0.08
    u = (0.35,)
    source(x, t) = 0.1 * sin(2pi * x) * exp(0.3 * t)
    bc = BorderConditions(; left=Periodic(), right=Periodic())

    fixed = AdvDiffModelMono(
        cap,
        D,
        u,
        u;
        source=source,
        bc=bc,
        scheme=Centered(),
    )
    moving = MovingAdvDiffModelMono(
        grid,
        body,
        D,
        u,
        u;
        wγ=(0.0,),
        source=source,
        bc=bc,
        scheme=Centered(),
    )

    nt = cap.ntotal
    u0 = [sin(2pi * cap.C_ω[i][1]) + 0.2 * cos(4pi * cap.C_ω[i][1]) for i in 1:nt]

    lay = fixed.diff.layout.offsets

    tspan = (0.0, 0.08)
    dtm = 0.004
    sol_fixed = solve_unsteady!(fixed, u0, tspan; dt=dtm, scheme=:CN, method=:direct, save_history=false)
    sol_moving = solve_unsteady_moving!(moving, u0, tspan; dt=dtm, scheme=:CN, method=:direct, save_history=false)

    xf = sol_fixed.system.x
    xm = sol_moving.system.x
    @test maximum(abs.(xm[lay.ω] .- xf[lay.ω])) < 3e-3
    @test maximum(abs.(xm .- xf)) < 4e-3
end
