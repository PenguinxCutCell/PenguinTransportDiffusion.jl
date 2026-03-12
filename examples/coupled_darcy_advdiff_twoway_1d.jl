using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy
using PenguinSolverCore
using PenguinTransportDiffusion

full_moments(grid) = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)

n = 65
grid = (range(0.0, 1.0; length=n),)
cap = assembled_capacity(full_moments(grid); bc=0.0)
nt = cap.ntotal

bc_darcy = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
ops = DiffusionOps(cap; periodic=periodic_flags(bc_darcy, 1))

λ0 = 1.0
α = 0.6
mobility_cb = function (base_model, concentration, t)
    lay = base_model.layout.offsets
    cω = length(concentration) >= last(lay.ω) ? concentration[lay.ω] : concentration
    csum = 0.0
    nfinite = 0
    @inbounds for value in cω
        isfinite(value) || continue
        csum += value
        nfinite += 1
    end
    cavg = nfinite == 0 ? 0.0 : csum / nfinite
    return max(1e-8, λ0 * (1 + α * cavg))
end

darcy_base = DarcyModelMono(cap, ops, λ0; source=(x, t) -> 0.0, bc_border=bc_darcy)
darcy_block = CoupledBlock(:darcy,
                           DarcyCoupledModelMono(darcy_base; mobility_update! = mobility_cb);
                           init=nothing,
                           cache=Dict{Symbol,Any}())

bc_transport = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
uzero = (zeros(nt),)
adv_base = AdvDiffModelMono(cap, 2e-3, uzero, uzero; source=0.0, bc=bc_transport, scheme=Centered())

lay = adv_base.diff.layout.offsets
nsys = maximum((last(lay.ω), last(lay.γ)))
c0 = zeros(Float64, nsys)
@inbounds for i in 1:nt
    c0[lay.ω[i]] = cap.C_ω[i][1]
end

transport_block = CoupledBlock(:transport,
                               AdvDiffCoupledModelMono(adv_base);
                               init=(concentration=c0,),
                               cache=Dict{Symbol,Any}())

mode = TwoWayCoupling(
    [:darcy, :transport];
    maxiter=30,
    atol=1e-10,
    rtol=1e-8,
    relaxation=Dict(:darcy => 0.7, :transport => 1.0),
    norm_fields=[:velocity, :concentration],
    sweep=:GaussSeidel,
    verbose=true,
)

problem = CoupledProblem(
    [darcy_block, transport_block],
    mode;
    maps=[
        CouplingMap(:darcy, :transport, :velocity),
        CouplingMap(:transport, :darcy, :concentration),
    ],
)

dt = 0.01
_, history = step_coupled!(problem, 0.0, dt; method=:direct, return_history=true)

println("Two-way Darcy <-> AdvDiff 1D")
println("  dt = $dt, alpha = $α")
println("  outer iterations = ", length(history.iterations))
println("  residual history = ", history.residuals)
println("  mean concentration = ", sum(transport_block.state.concentration[lay.ω]) / length(lay.ω))
