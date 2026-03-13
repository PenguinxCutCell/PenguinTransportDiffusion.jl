function _active_idx_mono(cap)
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

function _gaussian_exact_1d(x, t; x0, σ2, D, u)
    xt = x0 + u * t
    s2 = σ2 + 2D * t
    amp = sqrt(σ2 / s2)
    return amp * exp(-((x - xt)^2) / (2s2))
end

@testset "Moving mono MMS convergence without embedded interface" begin
    D = 2.5e-3
    u = 0.7
    x0 = 0.35
    σ = 0.06
    σ2 = σ^2
    tfinal = 0.01

    bc = BorderConditions(; left=Periodic(), right=Periodic())
    body(x, t) = -1.0

    errs = Float64[]
    hs = Float64[]
    for n in (33, 65, 129)
        xyz = (range(0.0, 1.0; length=n),)
        grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))

        moms = geometric_moments(x -> -1.0, xyz, Float64, nan; method=:vofijul)
        cap = assembled_capacity(moms; bc=0.0)
        idx = _active_idx_mono(cap)

        model = MovingAdvDiffModelMono(
            grid,
            body,
            D,
            (u,),
            (u,);
            source=0.0,
            bc=bc,
            scheme=Centered(),
        )

        u0 = [_gaussian_exact_1d(cap.C_ω[i][1], 0.0; x0=x0, σ2=σ2, D=D, u=u) for i in 1:cap.ntotal]
        h = step(xyz[1])
        dt = 0.2 * h^2
        sol = solve_unsteady_moving!(model, u0, (0.0, tfinal); dt=dt, scheme=:CN, method=:direct, save_history=false)

        lay = model.diff.layout.offsets
        uω = sol.system.x[lay.ω]
        err = sqrt(sum(cap.buf.V[i] * (uω[i] - _gaussian_exact_1d(cap.C_ω[i][1], tfinal; x0=x0, σ2=σ2, D=D, u=u))^2 for i in idx))
        push!(errs, err)
        push!(hs, h)
    end

    p1 = log(errs[1] / errs[2]) / log(hs[1] / hs[2])
    p2 = log(errs[2] / errs[3]) / log(hs[2] / hs[3])
    @test p1 > 0.8
    @test p2 > 0.8
end
