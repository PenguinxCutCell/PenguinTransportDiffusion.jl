function _active_idx_phase(cap)
    idx = Int[]
    LI = LinearIndices(cap.nnodes)
    N = length(cap.nnodes)
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        if all(d -> I[d] < cap.nnodes[d], 1:N)
            v = cap.buf.V[i]
            if isfinite(v) && v > 0.0
                push!(idx, i)
            end
        end
    end
    return idx
end

@testset "Moving diphasic MMS convergence with translating interface" begin
    D = 0.03
    u = 0.35
    c = 0.12
    xγ0 = 0.5
    sγ = 0.08
    k = 2pi
    tfinal = 0.03

    exact(x, t) = exp(t) * sin(k * (x - c * t))
    source_phase(x, t) = exp(t) * ((1.0 + D * k^2) * sin(k * (x - c * t)) + k * (u - c) * cos(k * (x - c * t)))
    bc = BorderConditions(; left=Periodic(), right=Periodic())
    ic = InterfaceConditions(; scalar=ScalarJump(1.0, 1.0, 0.0), flux=FluxJump(1.0, 1.0, 0.0))
    body1(x, t) = x - (xγ0 + sγ * t)

    errs = Float64[]
    hs = Float64[]
    for n in (33, 65, 129)
        xyz = (range(0.0, 1.0; length=n),)
        grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))

        moms1 = geometric_moments(x -> x - xγ0, xyz, Float64, nan; method=:vofijul)
        moms2 = geometric_moments(x -> -(x - xγ0), xyz, Float64, nan; method=:vofijul)
        cap1 = assembled_capacity(moms1; bc=0.0)
        cap2 = assembled_capacity(moms2; bc=0.0)
        idx1 = _active_idx_phase(cap1)
        idx2 = _active_idx_phase(cap2)

        model = MovingAdvDiffModelDiph(
            grid,
            body1,
            D,
            (u,),
            (u,),
            D,
            (u,),
            (u,);
            wγ=(sγ,),
            source=(source_phase, source_phase),
            bc=bc,
            ic=ic,
            scheme=Centered(),
        )

        u01 = [exact(cap1.C_ω[i][1], 0.0) for i in 1:cap1.ntotal]
        u02 = [exact(cap2.C_ω[i][1], 0.0) for i in 1:cap2.ntotal]
        u0 = vcat(u01, u02)

        h = step(xyz[1])
        dt = 0.2 * h
        sol = solve_unsteady_moving!(model, u0, (0.0, tfinal); dt=dt, scheme=:CN, method=:direct, save_history=false)

        lay = model.diff.layout.offsets
        u1 = sol.system.x[lay.ω1]
        u2 = sol.system.x[lay.ω2]
        err1 = sum(cap1.buf.V[i] * (u1[i] - exact(cap1.C_ω[i][1], tfinal))^2 for i in idx1)
        err2 = sum(cap2.buf.V[i] * (u2[i] - exact(cap2.C_ω[i][1], tfinal))^2 for i in idx2)
        push!(errs, sqrt(err1 + err2))
        push!(hs, h)
    end

    @test all(isfinite, errs)
    @test maximum(errs) < 0.1
    @test errs[end] <= 1.1 * errs[1]
end
