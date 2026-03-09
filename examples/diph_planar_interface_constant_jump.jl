using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinTransportDiffusion

function active_phase_indices(cap)
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
    L = 1.0
    x_gamma = 0.47
    m = 2.5
    C1 = 0.8
    C2 = m * C1

    grid = (range(0.0, L; length=121),)
    moms1 = geometric_moments((x) -> x - x_gamma, grid, Float64, nan; method=:vofijul)
    moms2 = geometric_moments((x) -> -(x - x_gamma), grid, Float64, nan; method=:vofijul)
    cap1 = assembled_capacity(moms1; bc=0.0)
    cap2 = assembled_capacity(moms2; bc=0.0)

    bc = BorderConditions(; left=Dirichlet(C1), right=Dirichlet(C2))
    ic = InterfaceConditions(; scalar=ScalarJump(m, 1.0, 0.0), flux=FluxJump(1.0, 1.0, 0.0))

    model = AdvDiffModelDiph(
        cap1,
        1.0,
        (0.0,),
        (0.0,),
        cap2,
        1.0,
        (0.0,),
        (0.0,);
        source=(0.0, 0.0),
        bc=bc,
        ic=ic,
        scheme=Centered(),
    )

    sol = solve_steady!(model; method=:direct)
    lay = model.diff.layout.offsets
    u1 = sol.solution[lay.ω1]
    u2 = sol.solution[lay.ω2]
    idx1 = active_phase_indices(cap1)
    idx2 = active_phase_indices(cap2)
    idx_gamma = [i for i in 1:cap1.ntotal if isfinite(cap1.buf.Γ[i]) && cap1.buf.Γ[i] > 0.0]

    e1 = maximum(abs.(u1[idx1] .- C1))
    e2 = maximum(abs.(u2[idx2] .- C2))
    ejump = maximum(abs.(sol.solution[lay.γ2][idx_gamma] .- m .* sol.solution[lay.γ1][idx_gamma]))

    println("Diph planar interface: constant Henry jump")
    println("  x_gamma=$x_gamma, m=$m")
    println("  exact: T1=$C1, T2=$C2")
    println("  max |T1-C1| = $e1")
    println("  max |T2-C2| = $e2")
    println("  max |gamma2-m*gamma1| = $ejump")
end

main()
