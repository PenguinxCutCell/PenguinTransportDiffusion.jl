using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinTransportDiffusion

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

D = 2.5e-3
ux, uy = 1.0, 0.7
x0, y0 = 0.35, 0.40
σ = 0.06
σ2 = σ^2
tfinal = 0.006

println("Gaussian no-body adv-diff convergence")
println("  D=$D, u=($ux,$uy), t_final=$tfinal")

errs = Float64[]
hs = Float64[]
for n in (17, 33, 65)
    grid = (range(0.0, 1.0; length=n), range(0.0, 1.0; length=n))
    moms = full_moments(grid)
    cap = assembled_capacity(moms; bc=0.0)

    bc = BorderConditions(; left=Periodic(), right=Periodic(), bottom=Periodic(), top=Periodic())
    model = AdvDiffModelMono(cap, D, (ux, uy), (ux, uy); source=0.0, bc=bc, scheme=Centered())

    u0 = [gaussian_exact(cap.C_ω[i][1], cap.C_ω[i][2], 0.0; x0=x0, y0=y0, σ2=σ2, D=D, ux=ux, uy=uy) for i in 1:cap.ntotal]
    dt = 0.2 * step(grid[1])^2
    sol = solve_unsteady!(model, u0, (0.0, tfinal); dt=dt, scheme=:CN, method=:direct, save_history=false)

    lay = model.diff.layout.offsets
    uω = sol.system.x[lay.ω]
    idx = active_physical_indices(cap)
    err = weighted_l2_error(cap, uω, (x, y) -> gaussian_exact(x, y, tfinal; x0=x0, y0=y0, σ2=σ2, D=D, ux=ux, uy=uy), idx)

    push!(errs, err)
    push!(hs, step(grid[1]))
    println("  n=$n  h=$(hs[end])  L2=$err")
end

for k in 1:2
    p = log(errs[k] / errs[k + 1]) / log(hs[k] / hs[k + 1])
    println("  order[$k] = $p")
end
