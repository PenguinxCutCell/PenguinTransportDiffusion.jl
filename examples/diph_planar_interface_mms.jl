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
    x_gamma = 0.43
    grid = (range(0.0, L; length=121),)

    moms1 = geometric_moments((x) -> x - x_gamma, grid, Float64, nan; method=:vofijul)
    moms2 = geometric_moments((x) -> -(x - x_gamma), grid, Float64, nan; method=:vofijul)
    cap1 = assembled_capacity(moms1; bc=0.0)
    cap2 = assembled_capacity(moms2; bc=0.0)

    D1, D2 = 0.8, 1.6
    U1, U2 = 0.5, -0.3
    q = 0.4

    # Piecewise linear exact profiles with continuity and diffusion-style flux continuity.
    B1 = q / D1
    B2 = q / D2
    A1 = 0.2
    A2 = A1 + (B1 - B2) * x_gamma

    T1(x) = A1 + B1 * x
    T2(x) = A2 + B2 * x

    s1 = U1 * B1
    s2 = U2 * B2

    bc = BorderConditions(; left=Dirichlet(T1(0.0)), right=Dirichlet(T2(L)))
    ic = InterfaceConditions(; scalar=ScalarJump(1.0, 1.0, 0.0), flux=FluxJump(1.0, 1.0, 0.0))

    model = AdvDiffModelDiph(
        cap1,
        D1,
        (U1,),
        (U1,),
        cap2,
        D2,
        (U2,),
        (U2,);
        source=(s1, s2),
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

    e1 = maximum(abs.(u1[idx1] .- T1.(getindex.(cap1.C_ω[idx1], 1))))
    e2 = maximum(abs.(u2[idx2] .- T2.(getindex.(cap2.C_ω[idx2], 1))))

    println("Diph planar interface MMS (piecewise linear)")
    println("  x_gamma=$x_gamma, D1=$D1, D2=$D2, U1=$U1, U2=$U2")
    println("  source1=$s1, source2=$s2")
    println("  max error phase 1: $e1")
    println("  max error phase 2: $e2")
end

main()
