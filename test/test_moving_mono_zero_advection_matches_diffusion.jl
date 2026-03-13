@testset "Moving mono zero-advection reduces to moving diffusion" begin
    xyz = (range(0.0, 1.0; length=73),)
    grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))
    body(x, t) = x - (0.42 + 0.06 * t)

    D = 0.05
    source(x, t) = 0.1 * exp(-t) * cos(2pi * x)
    bc = BorderConditions(; left=Periodic(), right=Periodic())

    model_diff = PenguinDiffusion.MovingDiffusionModelMono(
        grid,
        body,
        D;
        source=source,
        bc_border=bc,
        bc_interface=nothing,
    )
    model_advdiff = MovingAdvDiffModelMono(
        grid,
        body,
        D,
        (0.0,),
        (0.0,);
        wγ=(0.0,),
        source=source,
        bc=bc,
        bc_interface_diff=nothing,
        scheme=Centered(),
    )

    nt = prod(grid.n)
    xnodes = collect(xyz[1])
    u0 = [sin(2pi * xnodes[i]) for i in 1:nt]

    lay = model_advdiff.diff.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))
    t = 0.01
    dt = 0.004
    θ = 0.7

    sys_diff = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    sys_adv = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    PenguinDiffusion.assemble_unsteady_mono_moving!(sys_diff, model_diff, u0, t, dt; scheme=θ)
    PenguinTransportDiffusion.assemble_unsteady_mono_moving!(sys_adv, model_advdiff, u0, t, dt, θ)

    @test isapprox(norm(sys_adv.A - sys_diff.A), 0.0; atol=2e-10, rtol=2e-10)
    @test isapprox(norm(sys_adv.b - sys_diff.b), 0.0; atol=2e-10, rtol=2e-10)

    tspan = (0.0, 0.05)
    dtm = 0.003
    sol_diff = PenguinDiffusion.solve_unsteady_moving!(model_diff, u0, tspan; dt=dtm, scheme=:CN, method=:direct, save_history=false)
    sol_adv = solve_unsteady_moving!(model_advdiff, u0, tspan; dt=dtm, scheme=:CN, method=:direct, save_history=false)

    xd = sol_diff.system.x
    xa = sol_adv.system.x
    @test maximum(abs.(xa[lay.ω] .- xd[lay.ω])) < 2e-7
    @test maximum(abs.(xa .- xd)) < 6e-7
end
