using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinTransportDiffusion

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
    D = 5e-3
    U = 0.1
    T_in = 0.0
    T_wall = 1.0

    xc, yc = 0.5, 0.5
    R = 0.18
    grid = (range(0.0, 1.0; length=81), range(0.0, 1.0; length=81))
    body(x, y) = R - sqrt((x - xc)^2 + (y - yc)^2)

    moms = geometric_moments(body, grid, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)

    bc = BorderConditions(
        ;
        left=Inflow(T_in),
        right=Outflow(),
        bottom=Outflow(),
        top=Outflow(),
    )
    bc_interface = PenguinBCs.Robin(1.0, 0.0, T_wall)

    model = AdvDiffModelMono(
        cap,
        D,
        (U, 0.0),
        (U, 0.0);
        source=0.0,
        bc=bc,
        bc_interface_diff=bc_interface,
        scheme=Upwind1(),
    )

    u0 = fill(T_in, cap.ntotal)
    sol = solve_unsteady!(model, u0, (0.0, 0.05); dt=5e-4, scheme=:BE, method=:direct, save_history=false)

    lay = model.diff.layout.offsets
    u = sol.system.x[lay.ω]
    idx = active_physical_indices(cap)
    vals = u[idx]

    println("Hot cylinder in cross-flow (mono PTD)")
    println("  grid: 81x81, R=$R, D=$D, U=$U")
    println("  inlet T=$T_in, wall T=$T_wall")
    println("  final min/max T: ", extrema(vals))
    println("  final mean T: ", sum(cap.buf.V[i] * u[i] for i in idx) / sum(cap.buf.V[i] for i in idx))
end

main()
