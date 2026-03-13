@testset "Moving diphasic zero-advection reduces to moving diffusion" begin
    xyz = (range(0.0, 1.0; length=97),)
    grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))
    xγ0 = 0.48
    sγ = 0.07
    body1(x, t) = x - (xγ0 + sγ * t)

    D1 = 0.04
    D2 = 0.09
    source = (
        (x, t) -> 0.05 * exp(-0.5 * t) * sin(2pi * x),
        (x, t) -> -0.08 * exp(-0.2 * t) * cos(2pi * x),
    )
    bc = BorderConditions(; left=Periodic(), right=Periodic())
    ic = InterfaceConditions(; scalar=ScalarJump(1.0, 1.0, 0.0), flux=FluxJump(1.0, 1.0, 0.0))

    model_diff = PenguinDiffusion.MovingDiffusionModelDiph(
        grid,
        body1,
        D1,
        D2;
        source=source,
        bc_border=bc,
        ic=ic,
    )
    model_advdiff = MovingAdvDiffModelDiph(
        grid,
        body1,
        D1,
        (0.0,),
        (0.0,),
        D2,
        (0.0,),
        (0.0,);
        wγ=(0.0,),
        source=source,
        bc=bc,
        ic=ic,
        scheme=Centered(),
    )

    nt = prod(grid.n)
    u01 = [0.7 + 0.2 * sin(2pi * (i - 1) / nt) for i in 1:nt]
    u02 = [0.4 + 0.15 * cos(2pi * (i - 1) / nt) for i in 1:nt]
    u0 = vcat(u01, u02)

    lay = model_advdiff.diff.layout.offsets
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    t = 0.015
    dt = 0.003
    θ = 0.65

    sys_diff = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    sys_adv = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    PenguinDiffusion.assemble_unsteady_diph_moving!(sys_diff, model_diff, u0, t, dt; scheme=θ)
    PenguinTransportDiffusion.assemble_unsteady_diph_moving!(sys_adv, model_advdiff, u0, t, dt, θ)

    @test isapprox(norm(sys_adv.A - sys_diff.A), 0.0; atol=3e-10, rtol=3e-10)
    @test isapprox(norm(sys_adv.b - sys_diff.b), 0.0; atol=3e-10, rtol=3e-10)

    tspan = (0.0, 0.04)
    dtm = 0.002
    sol_diff = PenguinDiffusion.solve_unsteady_moving!(model_diff, u0, tspan; dt=dtm, scheme=:CN, method=:direct, save_history=false)
    sol_adv = solve_unsteady_moving!(model_advdiff, u0, tspan; dt=dtm, scheme=:CN, method=:direct, save_history=false)

    xd = sol_diff.system.x
    xa = sol_adv.system.x
    @test maximum(abs.(xa[lay.ω1] .- xd[lay.ω1])) < 3e-7
    @test maximum(abs.(xa[lay.ω2] .- xd[lay.ω2])) < 3e-7
    @test maximum(abs.(xa .- xd)) < 9e-7
end
