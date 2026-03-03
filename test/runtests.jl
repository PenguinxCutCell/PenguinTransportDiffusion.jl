using Test
using LinearAlgebra
using SparseArrays

using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinTransportDiffusion
using PenguinSolverCore

full_moments(grid) = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)

function active_physical_indices(cap)
    LI = LinearIndices(cap.nnodes)
    idx = Int[]
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

function weighted_l2_error(cap, uω, uexact, idx)
    num = 0.0
    for i in idx
        w = cap.buf.V[i]
        x = cap.C_ω[i]
        d = uω[i] - uexact(x...)
        num += w * d^2
    end
    return sqrt(num)
end

function gaussian_exact(x, y, t; x0, y0, σ2, D, ux, uy)
    xt = x0 + ux * t
    yt = y0 + uy * t
    s2 = σ2 + 2D * t
    amp = σ2 / s2
    r2 = (x - xt)^2 + (y - yt)^2
    return amp * exp(-r2 / (2s2))
end

@testset "No-body Gaussian convergence (periodic, CN)" begin
    D = 2.5e-3
    ux, uy = 1.0, 0.7
    x0, y0 = 0.35, 0.40
    σ = 0.06
    σ2 = σ^2
    tfinal = 0.006

    errs = Float64[]
    hs = Float64[]

    for n in (17, 33, 65)
        grid = (range(0.0, 1.0; length=n), range(0.0, 1.0; length=n))
        moms = full_moments(grid)
        cap = assembled_capacity(moms; bc=0.0)

        bc = BorderConditions(; left=Periodic(), right=Periodic(), bottom=Periodic(), top=Periodic())
        model = AdvDiffModelMono(
            cap,
            D,
            (ux, uy),
            (ux, uy);
            source=0.0,
            bc=bc,
            bc_interface_diff=nothing,
            scheme=Centered(),
        )

        u0 = [gaussian_exact(cap.C_ω[i][1], cap.C_ω[i][2], 0.0; x0=x0, y0=y0, σ2=σ2, D=D, ux=ux, uy=uy) for i in 1:cap.ntotal]
        dt = 0.2 * step(grid[1])^2
        sol = solve_unsteady!(model, u0, (0.0, tfinal); dt=dt, scheme=:CN, method=:direct, save_history=false)

        lay = model.diff.layout.offsets
        uω = sol.system.x[lay.ω]
        idx = active_physical_indices(cap)
        err = weighted_l2_error(cap, uω, (x, y) -> gaussian_exact(x, y, tfinal; x0=x0, y0=y0, σ2=σ2, D=D, ux=ux, uy=uy), idx)
        push!(errs, err)
        push!(hs, step(grid[1]))
    end

    p1 = log(errs[1] / errs[2]) / log(hs[1] / hs[2])
    p2 = log(errs[2] / errs[3]) / log(hs[2] / hs[3])
    @test p1 > 1.5
    @test p2 > 1.5
end

@testset "Gaussian outside circle, rotating flow, mass conservation" begin
    Ω = 1.4
    D = 1e-3
    R = 0.35
    xc, yc = 0.0, 0.0

    grid = (range(-1.0, 1.0; length=65), range(-1.0, 1.0; length=65))
    body(x, y) = R - sqrt((x - xc)^2 + (y - yc)^2)

    moms = geometric_moments(body, grid, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)

    bc = BorderConditions(; left=Periodic(), right=Periodic(), bottom=Periodic(), top=Periodic())
    bc_interface = PenguinBCs.Robin(0.0, 1.0, 0.0)

    uω = (
        (x, y, t) -> -Ω * (y - yc),
        (x, y, t) -> Ω * (x - xc),
    )
    uγ = (0.0, 0.0)

    model = AdvDiffModelMono(
        cap,
        D,
        uω,
        uγ;
        source=0.0,
        bc=bc,
        bc_interface_diff=bc_interface,
        scheme=Centered(),
    )

    xg, yg = 0.65, 0.0
    σ = 0.13
    u0 = [exp(-((cap.C_ω[i][1] - xg)^2 + (cap.C_ω[i][2] - yg)^2) / (2σ^2)) for i in 1:cap.ntotal]

    lay = model.diff.layout.offsets
    idx = active_physical_indices(cap)
    nsys = maximum((last(lay.ω), last(lay.γ)))
    nt = cap.ntotal

    # Pure Neumann interface introduces a γ-gauge nullspace: pin one active γ row.
    idxγ = Int[]
    LI = LinearIndices(cap.nnodes)
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        if all(d -> I[d] < cap.nnodes[d], 1:2) && isfinite(cap.buf.Γ[i]) && cap.buf.Γ[i] > 0.0
            push!(idxγ, i)
        end
    end
    @test !isempty(idxγ)
    rowγ = lay.γ[first(idxγ)]

    u = zeros(Float64, nsys)
    u[lay.ω] .= u0
    sys = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys); x=copy(u))
    tspan = (0.0, 0.05)
    dt = 0.0015
    t = tspan[1]
    tol = sqrt(eps(Float64))
    times = Float64[t]
    states = Vector{Vector{Float64}}(undef, 1)
    states[1] = copy(u)

    while t < tspan[2] - tol
        dt_step = min(dt, tspan[2] - t)
        assemble_unsteady_mono!(sys, model, u, t, dt_step, :BE)
        for j in 1:nsys
            sys.A[rowγ, j] = 0.0
        end
        sys.A[rowγ, rowγ] = 1.0
        sys.b[rowγ] = 0.0
        solve!(sys; method=:direct, reuse_factorization=false)
        u .= sys.x
        t += dt_step
        push!(times, t)
        push!(states, copy(u))
    end

    masses = Float64[]
    for state in states
        uωs = state[lay.ω]
        push!(masses, sum(cap.buf.V[i] * uωs[i] for i in idx))
    end

    @test all(isfinite, masses)
    @test all(isfinite, sys.x)
    rel_drift = maximum(abs.(masses .- first(masses))) / max(abs(first(masses)), 1e-14)
    @test rel_drift < 2e-3
end
