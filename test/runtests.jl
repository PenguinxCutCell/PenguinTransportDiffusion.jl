using Test
using LinearAlgebra
using SparseArrays

using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
import PenguinDiffusion
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

function interface_indices(cap)
    idx = Int[]
    LI = LinearIndices(cap.nnodes)
    N = length(cap.nnodes)
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        if all(d -> I[d] < cap.nnodes[d], 1:N)
            γ = cap.buf.Γ[i]
            if isfinite(γ) && γ > 0.0
                push!(idx, i)
            end
        end
    end
    return idx
end

function gaussian_exact(x, y, t; x0, y0, σ2, D, ux, uy)
    xt = x0 + ux * t
    yt = y0 + uy * t
    s2 = σ2 + 2D * t
    amp = σ2 / s2
    r2 = (x - xt)^2 + (y - yt)^2
    return amp * exp(-r2 / (2s2))
end

if "--mms-matrix" in ARGS
    include("test_mms_matrix_report.jl")
end

@testset "Mono PTD u=0 reduction to diffusion (steady+unsteady, outer BC variants)" begin
    grid = (range(0.0, 1.0; length=33),)
    cap = assembled_capacity(full_moments(grid); bc=0.0)
    nt = cap.ntotal

    bc_cases = (
        BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0)),
        BorderConditions(; left=Neumann(0.3), right=Neumann(-0.1)),
        BorderConditions(; left=Robin(2.0, 1.0, 0.4), right=Robin(1.1, 0.7, -0.2)),
    )

    u_prev = rand(nt)
    D = 0.3
    source = 0.2
    t = 0.2
    dt = 0.01
    θ = 0.6

    for bc in bc_cases
        ops = DiffusionOps(cap; periodic=periodic_flags(bc, 1))
        model_diff = PenguinDiffusion.DiffusionModelMono(cap, ops, D; source=source, bc_border=bc, bc_interface=nothing)
        model_ptd = AdvDiffModelMono(
            cap,
            D,
            (0.0,),
            (0.0,);
            source=source,
            bc=bc,
            bc_interface_diff=nothing,
            scheme=Centered(),
        )

        lay = model_ptd.diff.layout.offsets
        nsys = maximum((last(lay.ω), last(lay.γ)))

        sdiff = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
        sptd = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
        PenguinDiffusion.assemble_steady_mono!(sdiff, model_diff, t)
        PenguinTransportDiffusion.assemble_steady_mono!(sptd, model_ptd, t)
        @test isapprox(norm(sptd.A - sdiff.A), 0.0; atol=1e-12)
        @test isapprox(norm(sptd.b - sdiff.b), 0.0; atol=1e-12)

        udiff = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
        uptd = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
        PenguinDiffusion.assemble_unsteady_mono!(udiff, model_diff, u_prev, t, dt, θ)
        PenguinTransportDiffusion.assemble_unsteady_mono!(uptd, model_ptd, u_prev, t, dt, θ)
        @test isapprox(norm(uptd.A - udiff.A), 0.0; atol=1e-12)
        @test isapprox(norm(uptd.b - udiff.b), 0.0; atol=1e-12)
    end
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

@testset "Mono embedded boundary: advection leaves gamma rows unchanged" begin
    Ω = 1.4
    D = 1e-3
    R = 0.35
    grid = (range(-1.0, 1.0; length=49), range(-1.0, 1.0; length=49))
    body(x, y) = R - sqrt(x^2 + y^2)

    moms = geometric_moments(body, grid, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)

    bc = BorderConditions(; left=Periodic(), right=Periodic(), bottom=Periodic(), top=Periodic())
    bc_interface = PenguinBCs.Robin(0.0, 1.0, 0.0)
    ops = DiffusionOps(cap; periodic=periodic_flags(bc, 2))

    model_diff = PenguinDiffusion.DiffusionModelMono(cap, ops, D; source=0.0, bc_border=bc, bc_interface=bc_interface)
    model_zero = AdvDiffModelMono(cap, D, (0.0, 0.0), (0.0, 0.0); source=0.0, bc=bc, bc_interface_diff=bc_interface, scheme=Centered())
    model_adv = AdvDiffModelMono(
        cap,
        D,
        ((x, y, t) -> -Ω * y, (x, y, t) -> Ω * x),
        (0.0, 0.0);
        source=0.0,
        bc=bc,
        bc_interface_diff=bc_interface,
        scheme=Centered(),
    )

    lay = model_zero.diff.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))
    sys_diff = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    sys_zero = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    sys_adv = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))

    PenguinDiffusion.assemble_steady_mono!(sys_diff, model_diff, 0.0)
    PenguinTransportDiffusion.assemble_steady_mono!(sys_zero, model_zero, 0.0)
    PenguinTransportDiffusion.assemble_steady_mono!(sys_adv, model_adv, 0.0)

    @test isapprox(norm(sys_zero.A - sys_diff.A), 0.0; atol=1e-12)
    @test isapprox(norm(sys_zero.b - sys_diff.b), 0.0; atol=1e-12)

    @test isapprox(norm(sys_adv.A[lay.γ, :] - sys_diff.A[lay.γ, :]), 0.0; atol=1e-12)
    @test isapprox(norm(sys_adv.b[lay.γ] - sys_diff.b[lay.γ]), 0.0; atol=1e-12)
    @test norm(sys_adv.A[lay.ω, :] - sys_diff.A[lay.ω, :]) > 1e-10
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

@testset "Diphasic PTD structure and reuse" begin
    grid = (range(0.0, 1.0; length=61),)
    xγ = 0.47
    moms1 = geometric_moments((x) -> x - xγ, grid, Float64, nan; method=:vofijul)
    moms2 = geometric_moments((x) -> -(x - xγ), grid, Float64, nan; method=:vofijul)
    cap1 = assembled_capacity(moms1; bc=0.0)
    cap2 = assembled_capacity(moms2; bc=0.0)

    bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
    ic_cont = InterfaceConditions(; scalar=ScalarJump(1.0, 1.0, 0.0), flux=FluxJump(1.0, 1.0, 0.0))
    pflags = periodic_flags(bc, 1)
    ops1 = DiffusionOps(cap1; periodic=pflags)
    ops2 = DiffusionOps(cap2; periodic=pflags)

    diff_model = PenguinDiffusion.DiffusionModelDiph(
        cap1, ops1, 1.2, 0.0,
        cap2, ops2, 0.8, 0.0;
        bc_border=bc,
        ic=ic_cont,
    )
    ptd_zero = AdvDiffModelDiph(
        cap1, 1.2, (0.0,), (0.0,),
        cap2, 0.8, (0.0,), (0.0,);
        source=(0.0, 0.0),
        bc=bc,
        ic=ic_cont,
        scheme=Centered(),
    )

    lay = ptd_zero.diff.layout.offsets
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    sys_diff = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    sys_ptd0 = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    PenguinDiffusion.assemble_steady_diph!(sys_diff, diff_model, 0.0)
    PenguinTransportDiffusion.assemble_steady_diph!(sys_ptd0, ptd_zero, 0.0)

    @test isapprox(norm(sys_ptd0.A - sys_diff.A), 0.0; atol=1e-12)
    @test isapprox(norm(sys_ptd0.b - sys_diff.b), 0.0; atol=1e-12)

    ptd_adv = AdvDiffModelDiph(
        cap1, 1.2, (0.8,), (0.8,),
        cap2, 0.8, (-0.4,), (-0.4,);
        source=(0.0, 0.0),
        bc=bc,
        ic=ic_cont,
        scheme=Centered(),
    )
    sys_adv = LinearSystem(spzeros(Float64, nsys, nsys), zeros(Float64, nsys))
    PenguinTransportDiffusion.assemble_steady_diph!(sys_adv, ptd_adv, 0.0)

    γrows = vcat(collect(lay.γ1), collect(lay.γ2))
    ωrows = vcat(collect(lay.ω1), collect(lay.ω2))
    @test isapprox(norm(sys_adv.A[γrows, :] - sys_diff.A[γrows, :]), 0.0; atol=1e-12)
    @test isapprox(norm(sys_adv.b[γrows] - sys_diff.b[γrows]), 0.0; atol=1e-12)
    @test norm(sys_adv.A[ωrows, :] - sys_diff.A[ωrows, :]) > 1e-10

    C = 1.7
    bc_const = BorderConditions(; left=Dirichlet(C), right=Dirichlet(C))
    model_const = AdvDiffModelDiph(
        cap1, 1.1, (0.9,), (0.9,),
        cap2, 0.6, (-0.3,), (-0.3,);
        source=(0.0, 0.0),
        bc=bc_const,
        ic=ic_cont,
        scheme=Centered(),
    )
    sol_const = solve_steady!(model_const; method=:direct)
    u1c = sol_const.solution[lay.ω1]
    u2c = sol_const.solution[lay.ω2]
    idx1 = active_physical_indices(cap1)
    idx2 = active_physical_indices(cap2)
    @test !isempty(idx1)
    @test !isempty(idx2)
    @test maximum(abs.(u1c[idx1] .- C)) < 7e-2
    @test maximum(abs.(u2c[idx2] .- C)) < 7e-2

    m = 2.5
    C1 = 0.8
    bc_henry = BorderConditions(; left=Dirichlet(C1), right=Dirichlet(m * C1))
    ic_henry = InterfaceConditions(; scalar=ScalarJump(m, 1.0, 0.0), flux=FluxJump(1.0, 1.0, 0.0))
    model_henry = AdvDiffModelDiph(
        cap1, 0.9, (0.6,), (0.6,),
        cap2, 1.3, (-0.2,), (-0.2,);
        source=(0.0, 0.0),
        bc=bc_henry,
        ic=ic_henry,
        scheme=Centered(),
    )
    sol_henry = solve_steady!(model_henry; method=:direct)
    u1h = sol_henry.solution[lay.ω1]
    u2h = sol_henry.solution[lay.ω2]
    @test maximum(abs.(u1h[idx1] .- C1)) < 7e-2
    @test maximum(abs.(u2h[idx2] .- (m * C1))) < 7e-2
    idxγ = interface_indices(cap1)
    @test !isempty(idxγ)
    @test maximum(abs.(sol_henry.solution[lay.γ2][idxγ] .- m .* sol_henry.solution[lay.γ1][idxγ])) < 2e-6

    cap_mono = assembled_capacity(full_moments(grid); bc=0.0)
    model_mono_const = AdvDiffModelMono(
        cap_mono,
        0.7,
        (0.4,),
        (0.4,);
        source=0.0,
        bc=BorderConditions(; left=Dirichlet(0.5), right=Dirichlet(0.5)),
        bc_interface_diff=nothing,
        scheme=Centered(),
    )
    res_mono = solve_unsteady!(model_mono_const, fill(0.5, cap_mono.ntotal), (0.0, 0.1); dt=0.02, scheme=:BE, save_history=false)
    @test res_mono.reused_constant_operator

    res_diph = solve_unsteady!(model_const, vcat(fill(C, cap1.ntotal), fill(C, cap2.ntotal)), (0.0, 0.1); dt=0.02, scheme=:BE, save_history=false)
    @test res_diph.reused_constant_operator
end

include("test_coupling_adapter.jl")
include("test_coupling_integration.jl")
include("test_moving_mono_static_recovery.jl")
include("test_moving_diph_static_recovery.jl")
include("test_moving_mono_zero_advection_matches_diffusion.jl")
include("test_moving_diph_zero_advection_matches_diffusion.jl")
include("test_moving_mono_mms_convergence.jl")
include("test_moving_diph_mms_convergence.jl")
