using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinTransportDiffusion
import PenguinDiffusion

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

D = 0.02
a = 0.45
k = 2pi
tfinal = 0.04

exact(x, t) = exp(t) * sin(k * (x - a * t))
source(x, t) = (1.0 + D * k^2) * exact(x, t)

bc = BorderConditions(; left=Periodic(), right=Periodic())
body(x, t) = -1.0

println("Moving mono MMS convergence (no embedded interface)")
println("  D=$D, a=$a, tfinal=$tfinal")
println("  scheme=CN, advection=Centered")

errs = Float64[]
hs = Float64[]
for n in (33, 65, 129)
    xyz = (range(0.0, 1.0; length=n),)
    grid = PenguinDiffusion.CartesianGrid((first(xyz[1]),), (last(xyz[1]),), (length(xyz[1]),))

    moms = geometric_moments(x -> -1.0, xyz, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)
    idx = active_physical_indices(cap)

    model = MovingAdvDiffModelMono(
        grid,
        body,
        D,
        (a,),
        (a,);
        source=source,
        bc=bc,
        scheme=Centered(),
    )

    u0 = [exact(cap.C_ω[i][1], 0.0) for i in 1:cap.ntotal]
    h = step(xyz[1])
    dt = 0.25 * h

    sol = solve_unsteady_moving!(model, u0, (0.0, tfinal); dt=dt, scheme=:CN, method=:direct, save_history=false)
    lay = model.diff.layout.offsets
    uω = sol.system.x[lay.ω]
    err = sqrt(sum(cap.buf.V[i] * (uω[i] - exact(cap.C_ω[i][1], tfinal))^2 for i in idx))

    push!(errs, err)
    push!(hs, h)
    println("  n=$n  h=$h  L2=$err")
end

println("  convergence:")
for kidx in 1:2
    p = log(errs[kidx] / errs[kidx + 1]) / log(hs[kidx] / hs[kidx + 1])
    println("    p[$kidx] = $p")
end
