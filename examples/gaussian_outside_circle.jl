using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using SparseArrays
using PenguinBCs
using PenguinTransportDiffusion
using PenguinSolverCore

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

function main()
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

    # Pure Neumann interface introduces a γ-gauge nullspace: pin one active γ row.
    idxγ = Int[]
    LI = LinearIndices(cap.nnodes)
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        if all(d -> I[d] < cap.nnodes[d], 1:2) && isfinite(cap.buf.Γ[i]) && cap.buf.Γ[i] > 0.0
            push!(idxγ, i)
        end
    end
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

    println("Gaussian outside circle (rotating flow)")
    println("  grid: 65x65, tspan=$tspan")
    println("  mass(t0) = ", first(masses), ", mass(tf) = ", last(masses))
    println("  relative drift = ", maximum(abs.(masses .- first(masses))) / max(abs(first(masses)), 1e-14))

    uωf = sys.x[lay.ω]
    println("  min(T) = ", minimum(uωf[idx]), ", max(T) = ", maximum(uωf[idx]))

    open("gaussian_outside_circle_mass.csv", "w") do io
        println(io, "time,mass")
        for (tt, m) in zip(times, masses)
            println(io, string(tt, ",", m))
        end
    end
    println("  wrote gaussian_outside_circle_mass.csv")
end

main()
